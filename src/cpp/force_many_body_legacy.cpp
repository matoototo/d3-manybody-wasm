#include <emscripten/bind.h>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdint>
#ifdef __EMSCRIPTEN_PTHREADS__
#include <array>
#include <atomic>
#include <emscripten/threading.h>
#include <pthread.h>
#endif

// Function to calculate the Morton code (Z-order curve) for 2D coordinates
unsigned int legacySpreadBits(unsigned int value) {
    value &= 0x0000ffffu;
    value = (value | (value << 8)) & 0x00ff00ffu;
    value = (value | (value << 4)) & 0x0f0f0f0fu;
    value = (value | (value << 2)) & 0x33333333u;
    value = (value | (value << 1)) & 0x55555555u;
    return value;
}

unsigned int legacyMortonCode(double x, double y, double x0, double y0, double s, unsigned int scaleFactor) {
    unsigned int scaledX = static_cast<unsigned int>((x - x0) / s * scaleFactor);
    unsigned int scaledY = static_cast<unsigned int>((y - y0) / s * scaleFactor);
    return legacySpreadBits(scaledX) | (legacySpreadBits(scaledY) << 1);
}

class ForceManyBodyLegacy {
private:
    emscripten::val nodes;
    std::vector<float> strengths;
    std::function<double(const emscripten::val&, int, const emscripten::val&)> strength;
    double distanceMin2 = 1;
    double distanceMax2 = std::numeric_limits<double>::infinity();
    double theta2 = 0.9;
    double alpha;

    struct BodyData {
        float x, y, vx, vy;
    };
    std::vector<BodyData> bodyData;
    std::vector<BodyData> reorderedBodyData;

    struct Acceleration {
        float x, y;
    };
    std::vector<Acceleration> cachedAcceleration;
    std::vector<Acceleration> accelerationTrend;
    bool cachedAccelerationValid = false;
    int replayAge = 0;

    struct SortEntry {
        unsigned int code;
        int index;
    };
    std::vector<SortEntry> sortedIndices;
    std::vector<SortEntry> radixScratch;

    double* nodeBuffer = nullptr;
    int nodeCount = 0;
    static constexpr int stride = 4;

    std::function<double()> random;
    bool axesEnabled = false;
    float axisX = 0;
    float axisY = 0;
    float axisStrength = 0;

    struct QuadtreeNode {
        double cx, cy;  // Center of the node
        double s;       // Half-size of the region
        double value = 0;  // Total mass
        double x = 0, y = 0;  // Center of mass
        int firstChild = -1;  // Index of first child, or -1 if leaf node
        std::uint8_t childMask = 0;

        QuadtreeNode(double cx_, double cy_, double s_) : cx(cx_), cy(cy_), s(s_) {}
    };

    std::vector<QuadtreeNode> quadtreeNodes;
    void sortByMortonCode() {
        radixScratch.resize(sortedIndices.size());
        for (unsigned int shift = 0; shift < 32; shift += 8) {
            int counts[256] = {};
            for (const SortEntry& entry : sortedIndices) ++counts[(entry.code >> shift) & 0xffu];
            int offset = 0;
            for (int bucket = 0; bucket < 256; ++bucket) {
                const int count = counts[bucket];
                counts[bucket] = offset;
                offset += count;
            }
            for (const SortEntry& entry : sortedIndices) {
                radixScratch[counts[(entry.code >> shift) & 0xffu]++] = entry;
            }
            sortedIndices.swap(radixScratch);
        }
    }

    // Find bounding box of all nodes
    void findExtent(double& x0, double& y0, double& x1, double& y1) {
        int n = nodeCount;
        x0 = y0 = std::numeric_limits<double>::infinity();
        x1 = y1 = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            double xi = bodyData[i].x;
            double yi = bodyData[i].y;
            if (xi < x0) x0 = xi;
            if (xi > x1) x1 = xi;
            if (yi < y0) y0 = yi;
            if (yi > y1) y1 = yi;
        }
        // Slightly expand the bounds to avoid precision issues
        double dx = x1 - x0;
        double dy = y1 - y0;
        if (dx == 0) dx = 1;
        if (dy == 0) dy = 1;
        x0 -= dx * 0.1;
        x1 += dx * 0.1;
        y0 -= dy * 0.1;
        y1 += dy * 0.1;
    }

    // Insert a node into the quadtree
    void insertNode(int nodeIndex, int dataIndex) {
        QuadtreeNode& node = quadtreeNodes[nodeIndex];
        double x = bodyData[dataIndex].x;
        double y = bodyData[dataIndex].y;

        // If it's a leaf node
        if (node.firstChild < 0) {
            if (node.value == 0) {
                // First insertion
                node.x = x;
                node.y = y;
                node.value = strengths[dataIndex];
                return;
            }
            // Subdivide the node
            subdivideNode(nodeIndex);
        }

        // Determine the quadrant to insert into
        int quadIndex = getQuadrant(node, x, y);
        insertNode(node.firstChild + quadIndex, dataIndex);
    }

    // Subdivide a node into four quadrants
    void subdivideNode(int nodeIndex) {
        QuadtreeNode& node = quadtreeNodes[nodeIndex];
        double halfSize = node.s / 2;
        int firstChildIndex = quadtreeNodes.size();
        node.firstChild = firstChildIndex;

        // Create four children
        quadtreeNodes.emplace_back(node.cx - halfSize, node.cy - halfSize, halfSize);  // Bottom-left
        quadtreeNodes.emplace_back(node.cx + halfSize, node.cy - halfSize, halfSize);  // Bottom-right
        quadtreeNodes.emplace_back(node.cx - halfSize, node.cy + halfSize, halfSize);  // Top-left
        quadtreeNodes.emplace_back(node.cx + halfSize, node.cy + halfSize, halfSize);  // Top-right
    }

    // Determine the quadrant of a point (x, y) relative to node's center
    int getQuadrant(const QuadtreeNode& node, double x, double y) {
        int quad = 0;
        if (x >= node.cx) quad += 1;
        if (y >= node.cy) quad += 2;
        return quad;
    }

    // Build the quadtree
    void buildQuadtree(double x0, double y0, double x1, double y1) {
        double cx = (x0 + x1) / 2;
        double cy = (y0 + y1) / 2;
        double s = std::max(x1 - x0, y1 - y0) / 2 * 1.1;

        quadtreeNodes.clear();
        quadtreeNodes.emplace_back(cx, cy, s);

        int n = nodeCount;
        for (int i = 0; i < n; ++i) {
            insertNode(0, i);
        }
    }

    // Propagate masses and centers of mass upwards through the tree
    void propagate() {
        for (int i = quadtreeNodes.size() - 1; i >= 0; --i) {
            QuadtreeNode& node = quadtreeNodes[i];
            if (node.firstChild < 0) continue;

            // Combine the masses and centers of the children
            double mass = 0, x = 0, y = 0;
            std::uint8_t childMask = 0;
            for (int j = 0; j < 4; ++j) {
                QuadtreeNode& child = quadtreeNodes[node.firstChild + j];
                mass += child.value;
                x += child.x * child.value;
                y += child.y * child.value;
                if (child.value != 0) childMask |= static_cast<std::uint8_t>(1u << j);
            }

            if (mass > 0) {
                node.x = x / mass;
                node.y = y / mass;
            }
            node.value = mass;
            node.childMask = childMask;
        }
    }

    // Apply forces from the quadtree to a node
    void apply(int nodeIndex, BodyData& body) {
        QuadtreeNode& quad = quadtreeNodes[nodeIndex];

        if (quad.value == 0) return;

        double dx = quad.x - body.x;
        double dy = quad.y - body.y;
        double w = quad.s * 2;
        double d2 = dx * dx + dy * dy;

        if (w * w / theta2 < d2) {
            if (d2 < distanceMax2) {
                if (d2 == 0) {
                    dx = (random() - 0.5) * 1e-6;
                    dy = (random() - 0.5) * 1e-6;
                    d2 = dx * dx + dy * dy;
                }
                if (d2 < distanceMin2) d2 = distanceMin2;
                double factor = quad.value * alpha / d2;
                body.vx += static_cast<float>(dx * factor);
                body.vy += static_cast<float>(dy * factor);
            }
        } else if (quad.firstChild >= 0) {
            for (int i = 0; i < 4; ++i) {
                if (quad.childMask & (1u << i)) apply(quad.firstChild + i, body);
            }
        }
    }

    void applyAxes(BodyData& body) const {
        if (!axesEnabled) return;
        const float deltaX = static_cast<float>((axisX - body.x) * axisStrength * alpha);
        const float deltaY = static_cast<float>((axisY - body.y) * axisStrength * alpha);
        body.vx = static_cast<float>(body.vx + deltaX);
        body.vy = static_cast<float>(body.vy + deltaY);
    }

    void precomputeInternal() {
        const int n = nodeCount;
        bodyData.resize(n);
        for (int index = 0; index < n; ++index) {
            const int offset = index * stride;
            bodyData[index].x = static_cast<float>(nodeBuffer[offset]);
            bodyData[index].y = static_cast<float>(nodeBuffer[offset + 1]);
            bodyData[index].vx = 0;
            bodyData[index].vy = 0;
        }

        double x0, y0, x1, y1;
        findExtent(x0, y0, x1, y1);
        const double size = std::max(x1 - x0, y1 - y0);
        constexpr unsigned int scaleFactor = 65536;
        sortedIndices.resize(n);
        for (int index = 0; index < n; ++index) {
            sortedIndices[index] = {
                legacyMortonCode(bodyData[index].x, bodyData[index].y, x0, y0, size, scaleFactor),
                index
            };
        }
        sortByMortonCode();

        reorderedBodyData.resize(n);
        for (int index = 0; index < n; ++index) {
            reorderedBodyData[index] = bodyData[sortedIndices[index].index];
        }
        bodyData.swap(reorderedBodyData);
        buildQuadtree(x0, y0, x1, y1);
        propagate();
    }

#ifdef __EMSCRIPTEN_PTHREADS__
    struct ApplyTask {
        ForceManyBodyLegacy* force;
        int begin;
        int end;
    };

    static constexpr int kMaxThreadCount = 8;
    static constexpr int kApplyChunkSize = 32;
    std::array<pthread_t, kMaxThreadCount - 1> workerThreads;
    std::array<ApplyTask, kMaxThreadCount - 1> workerTasks;
    pthread_mutex_t workerMutex;
    pthread_cond_t workAvailable;
    pthread_cond_t workComplete;
    int workGeneration = 0;
    int completedWorkers = 0;
    int activeWorkerCount = 0;
    bool workersStarted = false;
    bool workersStopping = false;
    bool precomputePending = false;
    std::atomic<int> nextApplyIndex = 0;
    int applyEnd = 0;

    static void applyRange(ApplyTask& task) {
        if (task.begin < 0) {
            if (task.begin == -1) {
                task.force->precomputeInternal();
            } else {
                while (true) {
                    const int begin = task.force->nextApplyIndex.fetch_add(kApplyChunkSize, std::memory_order_relaxed);
                    if (begin >= task.force->applyEnd) break;
                    const int end = std::min(task.force->applyEnd, begin + kApplyChunkSize);
                    for (int index = begin; index < end; ++index) {
                        task.force->apply(0, task.force->bodyData[index]);
                        task.force->applyAxes(task.force->bodyData[index]);
                    }
                }
            }
            return;
        }
        for (int index = task.begin; index < task.end; ++index) {
            task.force->apply(0, task.force->bodyData[index]);
            task.force->applyAxes(task.force->bodyData[index]);
        }
    }

    static void* workerEntry(void* pointer) {
        ApplyTask& task = *static_cast<ApplyTask*>(pointer);
        ForceManyBodyLegacy& force = *task.force;
        int observedGeneration = 0;
        pthread_mutex_lock(&force.workerMutex);
        while (true) {
            while (!force.workersStopping && observedGeneration == force.workGeneration) {
                pthread_cond_wait(&force.workAvailable, &force.workerMutex);
            }
            if (force.workersStopping) break;
            observedGeneration = force.workGeneration;
            pthread_mutex_unlock(&force.workerMutex);
            applyRange(task);
            pthread_mutex_lock(&force.workerMutex);
            if (++force.completedWorkers == force.activeWorkerCount) {
                pthread_cond_signal(&force.workComplete);
            }
        }
        pthread_mutex_unlock(&force.workerMutex);
        return nullptr;
    }
#endif

    void applyPrepared() {
        const int n = nodeCount;
        for (int index = 0; index < n; ++index) {
            const int originalIndex = sortedIndices[index].index;
            const int offset = originalIndex * stride;
            bodyData[index].vx = static_cast<float>(nodeBuffer[offset + 2]);
            bodyData[index].vy = static_cast<float>(nodeBuffer[offset + 3]);
            nodeBuffer[offset] = bodyData[index].x;
            nodeBuffer[offset + 1] = bodyData[index].y;
            nodeBuffer[offset + 2] = bodyData[index].vx;
            nodeBuffer[offset + 3] = bodyData[index].vy;
        }

#ifdef __EMSCRIPTEN_PTHREADS__
        prepareWorkers();
        if (activeWorkerCount == 0) {
            for (int index = 0; index < n; ++index) {
                apply(0, bodyData[index]);
                applyAxes(bodyData[index]);
            }
        } else {
            pthread_mutex_lock(&workerMutex);
            completedWorkers = 0;
            nextApplyIndex.store(0, std::memory_order_relaxed);
            applyEnd = n;
            for (int worker = 0; worker < activeWorkerCount; ++worker) {
                workerTasks[worker].begin = -2;
                workerTasks[worker].end = -2;
            }
            ++workGeneration;
            pthread_cond_broadcast(&workAvailable);
            pthread_mutex_unlock(&workerMutex);
            ApplyTask mainTask = {this, -2, -2};
            applyRange(mainTask);
            pthread_mutex_lock(&workerMutex);
            while (completedWorkers < activeWorkerCount) {
                pthread_cond_wait(&workComplete, &workerMutex);
            }
            pthread_mutex_unlock(&workerMutex);
        }
#else
        for (int index = 0; index < n; ++index) {
            apply(0, bodyData[index]);
            applyAxes(bodyData[index]);
        }
#endif

        const bool hadCachedAcceleration = cachedAccelerationValid;
        const float trendDivisor = static_cast<float>(replayAge + 1);
        for (int index = 0; index < n; ++index) {
            const int originalIndex = sortedIndices[index].index;
            const int offset = originalIndex * stride;
            if (alpha != 0) {
                const float accelerationX = static_cast<float>(
                    (bodyData[index].vx - nodeBuffer[offset + 2]) / alpha
                );
                const float accelerationY = static_cast<float>(
                    (bodyData[index].vy - nodeBuffer[offset + 3]) / alpha
                );
                accelerationTrend[originalIndex].x = hadCachedAcceleration
                    ? (accelerationX - cachedAcceleration[originalIndex].x) / trendDivisor
                    : 0;
                accelerationTrend[originalIndex].y = hadCachedAcceleration
                    ? (accelerationY - cachedAcceleration[originalIndex].y) / trendDivisor
                    : 0;
                cachedAcceleration[originalIndex] = {accelerationX, accelerationY};
            }
            nodeBuffer[offset + 2] = bodyData[index].vx;
            nodeBuffer[offset + 3] = bodyData[index].vy;
        }
        cachedAccelerationValid = alpha != 0;
        replayAge = 0;
        quadtreeNodes.clear();
    }

public:
    ForceManyBodyLegacy()
        : strength([](const emscripten::val&, int, const emscripten::val&) { return -30; }),
          random([]() { return std::rand() / double(RAND_MAX); })
    {
#ifdef __EMSCRIPTEN_PTHREADS__
        pthread_mutex_init(&workerMutex, nullptr);
        pthread_cond_init(&workAvailable, nullptr);
        pthread_cond_init(&workComplete, nullptr);
#endif
    }

    ~ForceManyBodyLegacy() {
#ifdef __EMSCRIPTEN_PTHREADS__
        shutdownWorkers();
        pthread_cond_destroy(&workComplete);
        pthread_cond_destroy(&workAvailable);
        pthread_mutex_destroy(&workerMutex);
#endif
    }

    void prepareWorkers() {
#ifdef __EMSCRIPTEN_PTHREADS__
        if (workersStarted) return;
        workersStopping = false;
        completedWorkers = 0;
        workGeneration = 0;
        const int logicalCoreCount = std::max(1, emscripten_num_logical_cores());
        const int usefulThreadCount = std::max(1, (nodeCount + 255) / 256);
        const int desiredThreadCount = std::min({kMaxThreadCount, logicalCoreCount, usefulThreadCount});
        activeWorkerCount = 0;
        for (int worker = 0; worker < desiredThreadCount - 1; ++worker) {
            workerTasks[worker] = {this, 0, 0};
            if (pthread_create(&workerThreads[worker], nullptr, workerEntry, &workerTasks[worker]) != 0) break;
            ++activeWorkerCount;
        }
        workersStarted = true;
#endif
    }

    void shutdownWorkers() {
#ifdef __EMSCRIPTEN_PTHREADS__
        if (!workersStarted) return;
        pthread_mutex_lock(&workerMutex);
        workersStopping = true;
        ++workGeneration;
        pthread_cond_broadcast(&workAvailable);
        pthread_mutex_unlock(&workerMutex);
        for (int worker = 0; worker < activeWorkerCount; ++worker) pthread_join(workerThreads[worker], nullptr);
        activeWorkerCount = 0;
        workersStarted = false;
#endif
    }

    void beginPrecompute(double alpha_) {
        alpha = alpha_;
        if (nodeBuffer == nullptr || nodeCount == 0) return;
#ifdef __EMSCRIPTEN_PTHREADS__
        prepareWorkers();
        if (activeWorkerCount == 0) {
            precomputeInternal();
            return;
        }
        pthread_mutex_lock(&workerMutex);
        completedWorkers = 0;
        workerTasks[0].begin = -1;
        workerTasks[0].end = -1;
        for (int worker = 1; worker < activeWorkerCount; ++worker) {
            workerTasks[worker].begin = 0;
            workerTasks[worker].end = 0;
        }
        ++workGeneration;
        pthread_cond_broadcast(&workAvailable);
        pthread_mutex_unlock(&workerMutex);
        precomputePending = true;
#else
        precomputeInternal();
#endif
    }

    void force(double alpha_) {
        alpha = alpha_;
        if (nodeBuffer == nullptr || nodeCount == 0) return;
#ifdef __EMSCRIPTEN_PTHREADS__
        if (precomputePending) {
            pthread_mutex_lock(&workerMutex);
            while (completedWorkers < activeWorkerCount) {
                pthread_cond_wait(&workComplete, &workerMutex);
            }
            pthread_mutex_unlock(&workerMutex);
            precomputePending = false;
        } else {
            precomputeInternal();
        }
#else
        precomputeInternal();
#endif
        applyPrepared();
    }

    void replay(double alpha_) {
        if (nodeBuffer == nullptr || !cachedAccelerationValid) return;
        const float replayAlpha = static_cast<float>(alpha_);
        const int age = advanceReplayAge();
        const int predictionAge = std::min(age, 4);
        for (int index = 0; index < nodeCount; ++index) {
            const int offset = index * stride;
            nodeBuffer[offset + 2] += (
                cachedAcceleration[index].x + accelerationTrend[index].x * predictionAge
            ) * replayAlpha;
            nodeBuffer[offset + 3] += (
                cachedAcceleration[index].y + accelerationTrend[index].y * predictionAge
            ) * replayAlpha;
        }
    }

    int advanceReplayAge() { return ++replayAge; }
    uintptr_t getCachedAccelerationPointer() const {
        return reinterpret_cast<uintptr_t>(cachedAcceleration.data());
    }
    uintptr_t getAccelerationTrendPointer() const {
        return reinterpret_cast<uintptr_t>(accelerationTrend.data());
    }

    void initialize() {
        if (nodes.isUndefined()) return;
        int n = nodes["length"].as<int>();
        nodeCount = n;
        strengths.resize(n);
        cachedAcceleration.resize(n);
        accelerationTrend.resize(n);
        cachedAccelerationValid = false;
        replayAge = 0;
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            strengths[i] = static_cast<float>(strength(node, i, nodes));
        }
    }

    void setNodes(const emscripten::val& _nodes) {
        nodes = _nodes;
        int n = nodes["length"].as<int>();
        bodyData.resize(n);
        nodeCount = n;
        initialize();
    }

    void setNodeBuffer(uintptr_t ptr, int count) {
        nodeBuffer = reinterpret_cast<double*>(ptr);
        nodeCount = count;
    }

    void setStrength(const emscripten::val& _strength) {
        if (_strength.typeOf().as<std::string>() == "function") {
            strength = _strength.as<std::function<double(const emscripten::val&, int, const emscripten::val&)>>();
        } else {
            double s = _strength.as<double>();
            strength = [s](const emscripten::val&, int, const emscripten::val&) { return s; };
        }
        initialize();
    }

    void setStrengths(const emscripten::val& values) {
        const int count = values["length"].as<int>();
        strengths.resize(count);
        for (int index = 0; index < count; ++index) {
            strengths[index] = values[index].as<float>();
        }
    }

    emscripten::val getStrength() const {
        return emscripten::val(strength);
    }

    void setDistanceMin(double d) {
        distanceMin2 = static_cast<float>(d * d);
    }

    double getDistanceMin() const {
        return static_cast<double>(std::sqrt(distanceMin2));
    }

    void setDistanceMax(double d) {
        distanceMax2 = static_cast<float>(d * d);
    }

    double getDistanceMax() const {
        return static_cast<double>(std::sqrt(distanceMax2));
    }

    void setTheta(double t) {
        theta2 = static_cast<float>(t * t);
    }

    double getTheta() const {
        return static_cast<double>(std::sqrt(theta2));
    }

    void setAxes(double x, double y, double strengthValue) {
        axisX = static_cast<float>(x);
        axisY = static_cast<float>(y);
        axisStrength = static_cast<float>(strengthValue);
        axesEnabled = true;
    }

};

ForceManyBodyLegacy* createForceManyBodyLegacy() {
    return new ForceManyBodyLegacy();
}

EMSCRIPTEN_BINDINGS(force_many_body_legacy_module) {
    emscripten::class_<ForceManyBodyLegacy>("ForceManyBodyLegacy")
        .constructor<>()
        .function("force", &ForceManyBodyLegacy::force)
        .function("replay", &ForceManyBodyLegacy::replay)
        .function("advanceReplayAge", &ForceManyBodyLegacy::advanceReplayAge)
        .function("getCachedAccelerationPointer", &ForceManyBodyLegacy::getCachedAccelerationPointer)
        .function("getAccelerationTrendPointer", &ForceManyBodyLegacy::getAccelerationTrendPointer)
        .function("setNodes", &ForceManyBodyLegacy::setNodes)
        .function("setNodeBuffer", &ForceManyBodyLegacy::setNodeBuffer)
        .function("setStrength", &ForceManyBodyLegacy::setStrength)
        .function("setStrengths", &ForceManyBodyLegacy::setStrengths)
        .function("getStrength", &ForceManyBodyLegacy::getStrength)
        .function("setDistanceMin", &ForceManyBodyLegacy::setDistanceMin)
        .function("getDistanceMin", &ForceManyBodyLegacy::getDistanceMin)
        .function("setDistanceMax", &ForceManyBodyLegacy::setDistanceMax)
        .function("getDistanceMax", &ForceManyBodyLegacy::getDistanceMax)
        .function("setTheta", &ForceManyBodyLegacy::setTheta)
        .function("getTheta", &ForceManyBodyLegacy::getTheta)
        .function("setAxes", &ForceManyBodyLegacy::setAxes)
        .function("beginPrecompute", &ForceManyBodyLegacy::beginPrecompute)
        .function("prepareWorkers", &ForceManyBodyLegacy::prepareWorkers)
        .function("shutdownWorkers", &ForceManyBodyLegacy::shutdownWorkers);

    emscripten::function("createForceManyBodyLegacy", &createForceManyBodyLegacy, emscripten::allow_raw_pointers());
}
