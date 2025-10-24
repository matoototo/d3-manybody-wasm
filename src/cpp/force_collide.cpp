#include <emscripten/bind.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <functional>
#include <cstdlib>

struct CollideQuadNode {
    float x0, y0, x1, y1;
    int dataIndex;
    int children[4];
    float maxRadius;

    CollideQuadNode()
        : x0(0), y0(0), x1(0), y1(0), dataIndex(-1), children{-1, -1, -1, -1}, maxRadius(0) {}

    CollideQuadNode(float x0_, float y0_, float x1_, float y1_)
        : x0(x0_), y0(y0_), x1(x1_), y1(y1_), dataIndex(-1), children{-1, -1, -1, -1}, maxRadius(0) {}
};

class ForceCollide {
private:
    emscripten::val nodes;
    std::function<double(const emscripten::val&, int, const emscripten::val&)> radiusFunc =
        [](const emscripten::val&, int, const emscripten::val&) { return 1.0; };
    bool radiusFuncIsCustom = false;
    float radiusConstant = 1.0f;
    float strength = 1.0f;
    int iterations = 1;

    std::vector<float> radii;
    std::vector<float> pxCache;
    std::vector<float> pyCache;
    std::vector<CollideQuadNode> quadtree;
    std::vector<std::pair<unsigned int, int>> mortonIndices;
    std::vector<int> sortedIndices;

    float* nodeBuffer = nullptr;
    float* radiusBufferPtr = nullptr;
    int nodeCount = 0;
    int radiusBufferCount = 0;
    static constexpr int stride = 4;

    std::function<float()> random = []() {
        return static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
    };

    static unsigned int mortonCode(float x, float y, float x0, float y0, float span, unsigned int scaleFactor) {
        float nx = span != 0 ? (x - x0) / span : 0.0f;
        if (nx < 0.0) nx = 0.0;
        if (nx > 0.999999) nx = 0.999999f;
        float ny = span != 0 ? (y - y0) / span : 0.0f;
        if (ny < 0.0) ny = 0.0;
        if (ny > 0.999999) ny = 0.999999f;
        unsigned int scaledX = static_cast<unsigned int>(nx * scaleFactor);
        unsigned int scaledY = static_cast<unsigned int>(ny * scaleFactor);
        unsigned int morton = 0;
        for (unsigned int i = 0; i < sizeof(unsigned int) * 8 / 2; ++i) {
            morton |= ((scaledX & (1u << i)) << i) | ((scaledY & (1u << i)) << (i + 1));
        }
        return morton;
    }

    void computeRadii() {
        if (radiusBufferPtr && radiusBufferCount >= nodeCount && nodeCount > 0) {
            return;
        }
        if (nodeCount <= 0 || nodes.isUndefined()) {
            radii.clear();
            return;
        }
            radii.resize(nodeCount);
        if (radiusFuncIsCustom) {
            for (int i = 0; i < nodeCount; ++i) {
                double value = radiusFunc(nodes[i], i, nodes);
                radii[i] = static_cast<float>(std::max(0.0, value));
            }
        } else {
            float constant = std::max(0.0f, radiusConstant);
            std::fill(radii.begin(), radii.end(), constant);
        }
    }

    int childIndexFor(const CollideQuadNode& quad, int dataIndex) const {
        float midX = 0.5f * (quad.x0 + quad.x1);
        float midY = 0.5f * (quad.y0 + quad.y1);
        int child = 0;
        if (pxCache[dataIndex] >= midX) child |= 1;
        if (pyCache[dataIndex] >= midY) child |= 2;
        return child;
    }

    void subdivide(int quadIndex) {
        CollideQuadNode& quad = quadtree[quadIndex];

        float midX = 0.5f * (quad.x0 + quad.x1);
        float midY = 0.5f * (quad.y0 + quad.y1);

        int baseIndex = quadtree.size();
        quad.children[0] = baseIndex;
        quad.children[1] = baseIndex + 1;
        quad.children[2] = baseIndex + 2;
        quad.children[3] = baseIndex + 3;

        quadtree.emplace_back(quad.x0, quad.y0, midX, midY);
        quadtree.emplace_back(midX, quad.y0, quad.x1, midY);
        quadtree.emplace_back(quad.x0, midY, midX, quad.y1);
        quadtree.emplace_back(midX, midY, quad.x1, quad.y1);
    }

    void insertNode(int quadIndex, int dataIndex) {
        CollideQuadNode& quad = quadtree[quadIndex];

        if (quad.children[0] == -1) {
            if (quad.dataIndex == -1) {
                quad.dataIndex = dataIndex;
                return;
            }

            int existingIndex = quad.dataIndex;
            if (pxCache[existingIndex] == pxCache[dataIndex] &&
                pyCache[existingIndex] == pyCache[dataIndex]) {
                pxCache[dataIndex] += (random() - 0.5) * 1e-6;
                pyCache[dataIndex] += (random() - 0.5) * 1e-6;
            }

            quad.dataIndex = -1;
            if (quadtree.capacity() < quadtree.size() + 4) {
                quadtree.reserve((quadtree.size() + 4) * 2);
            }
            subdivide(quadIndex);
            int existingChild = childIndexFor(quad, existingIndex);
            insertNode(quad.children[existingChild], existingIndex);
        }

        int child = childIndexFor(quad, dataIndex);
        insertNode(quad.children[child], dataIndex);
    }

    void updateMaxRadius(int quadIndex, const float* radiiValues) {
        CollideQuadNode& quad = quadtree[quadIndex];
        if (quad.children[0] == -1) {
            quad.maxRadius = (quad.dataIndex >= 0) ? radiiValues[quad.dataIndex] : 0.0;
            return;
        }
        float maxR = 0.0f;
        for (int i = 0; i < 4; ++i) {
            int child = quad.children[i];
            if (child >= 0) {
                updateMaxRadius(child, radiiValues);
                if (quadtree[child].maxRadius > maxR) {
                    maxR = quadtree[child].maxRadius;
                }
            }
        }
        quad.maxRadius = maxR;
    }

    void resolvePair(int nodeIndex, int otherIndex, const float* radiiValues) {
        if (nodeIndex == otherIndex || otherIndex < nodeIndex) return;

        float dx = pxCache[nodeIndex] - pxCache[otherIndex];
        float dy = pyCache[nodeIndex] - pyCache[otherIndex];
        float r = radiiValues[nodeIndex] + radiiValues[otherIndex];
        if (r <= 0.0) return;

        float r2 = r * r;
        float l2 = dx * dx + dy * dy;
        if (l2 >= r2) return;

        if (l2 == 0.0f) {
            dx = (random() - 0.5f) * 1e-6f;
            dy = (random() - 0.5f) * 1e-6f;
            l2 = dx * dx + dy * dy;
            if (l2 == 0.0f) return;
        }

        float l = std::sqrt(l2);
        float adjustment = (l - r) / l * strength;
        dx *= adjustment;
        dy *= adjustment;

        int offsetNode = nodeIndex * stride;
        nodeBuffer[offsetNode + 2] -= dx;
        nodeBuffer[offsetNode + 3] -= dy;

        int offsetOther = otherIndex * stride;
        nodeBuffer[offsetOther + 2] += dx;
        nodeBuffer[offsetOther + 3] += dy;
    }

    void visit(int quadIndex, int nodeIndex, const float* radiiValues) {
        CollideQuadNode& quad = quadtree[quadIndex];
        float r = radiiValues[nodeIndex] + quad.maxRadius;
        if (r <= 0.0f) return;

        float px = pxCache[nodeIndex];
        float py = pyCache[nodeIndex];
        if (quad.x0 > px + r || quad.x1 < px - r ||
            quad.y0 > py + r || quad.y1 < py - r) {
            return;
        }

        if (quad.children[0] == -1) {
            if (quad.dataIndex >= 0) {
                resolvePair(nodeIndex, quad.dataIndex, radiiValues);
            }
            return;
        }

        for (int i = 0; i < 4; ++i) {
            int child = quad.children[i];
            if (child >= 0) {
                visit(child, nodeIndex, radiiValues);
            }
        }
    }

public:
    ForceCollide() = default;

    void force(double /*alpha*/) {
        if (nodeBuffer == nullptr || nodeCount == 0) return;

        const float* radiiValues = nullptr;
        if (radiusBufferPtr && radiusBufferCount >= nodeCount) {
            radiiValues = radiusBufferPtr;
        } else {
            if (static_cast<int>(radii.size()) < nodeCount) {
                computeRadii();
            }
            if (static_cast<int>(radii.size()) < nodeCount) {
                radii.resize(nodeCount, std::max(0.0f, radiusConstant));
            }
            radiiValues = radii.data();
        }

        pxCache.resize(nodeCount);
        pyCache.resize(nodeCount);
        mortonIndices.resize(nodeCount);
        sortedIndices.resize(nodeCount);

        for (int iter = 0; iter < iterations; ++iter) {
            float x0 = std::numeric_limits<float>::infinity();
            float y0 = std::numeric_limits<float>::infinity();
            float x1 = -std::numeric_limits<float>::infinity();
            float y1 = -std::numeric_limits<float>::infinity();

            for (int i = 0; i < nodeCount; ++i) {
                int offset = i * stride;
                float x = nodeBuffer[offset];
                float y = nodeBuffer[offset + 1];
                float vx = nodeBuffer[offset + 2];
                float vy = nodeBuffer[offset + 3];
                float px = x + vx;
                float py = y + vy;
                pxCache[i] = px;
                pyCache[i] = py;
                if (px < x0) x0 = px;
                if (px > x1) x1 = px;
                if (py < y0) y0 = py;
                if (py > y1) y1 = py;
            }

            if (!std::isfinite(x0)) {
                x0 = y0 = x1 = y1 = 0.0f;
            }

            float dx = x1 - x0;
            float dy = y1 - y0;
            if (dx == 0.0f) dx = 1.0f;
            if (dy == 0.0f) dy = 1.0f;
            x0 -= dx * 0.1;
            x1 += dx * 0.1;
            y0 -= dy * 0.1;
            y1 += dy * 0.1;

            float span = std::max(x1 - x0, y1 - y0);
            if (span == 0.0f) span = 1.0f;
            unsigned int scaleFactor = 65535u;

            for (int i = 0; i < nodeCount; ++i) {
                unsigned int code = mortonCode(pxCache[i], pyCache[i], x0, y0, span, scaleFactor);
                mortonIndices[i] = {code, i};
            }

            std::sort(mortonIndices.begin(), mortonIndices.end(),
                      [](const auto& a, const auto& b) {
                          return a.first < b.first;
                      });

            for (int i = 0; i < nodeCount; ++i) {
                sortedIndices[i] = mortonIndices[i].second;
            }

            quadtree.clear();
            quadtree.reserve(nodeCount * 2 + 1);
            quadtree.emplace_back(x0, y0, x1, y1);

            for (int index : sortedIndices) {
                insertNode(0, index);
            }

            updateMaxRadius(0, radiiValues);

            for (int index : sortedIndices) {
                visit(0, index, radiiValues);
            }
        }
    }

    void setNodes(const emscripten::val& _nodes) {
        nodes = _nodes;
        nodeCount = nodes.isUndefined() ? 0 : nodes["length"].as<int>();
        if (!(radiusBufferPtr && radiusBufferCount >= nodeCount)) {
            computeRadii();
        }
    }

    void setNodeBuffer(uintptr_t ptr, int count) {
        nodeBuffer = (ptr && count > 0) ? reinterpret_cast<float*>(ptr) : nullptr;
        nodeCount = count;
    }

    void setRadiusBuffer(uintptr_t ptr, int count) {
        radiusBufferPtr = (ptr && count > 0) ? reinterpret_cast<float*>(ptr) : nullptr;
        radiusBufferCount = count;
    }

    void setRadius(const emscripten::val& _radius) {
        if (_radius.typeOf().as<std::string>() == "function") {
            radiusFunc = _radius.as<std::function<double(const emscripten::val&, int, const emscripten::val&)>>();
            radiusFuncIsCustom = true;
        } else {
            radiusConstant = _radius.as<double>();
            radiusFuncIsCustom = false;
        }
        if (!(radiusBufferPtr && radiusBufferCount >= nodeCount)) {
            computeRadii();
        }
    }

    emscripten::val getRadius() const {
        if (radiusFuncIsCustom) {
            return emscripten::val::undefined();
        }
        return emscripten::val(radiusConstant);
    }

    void setStrength(double s) {
        strength = s;
    }

    double getStrength() const {
        return strength;
    }

    void setIterations(int iters) {
        iterations = std::max(1, iters);
    }

    int getIterations() const {
        return iterations;
    }
};

ForceCollide* createForceCollide() {
    return new ForceCollide();
}

EMSCRIPTEN_BINDINGS(force_collide_module) {
    emscripten::class_<ForceCollide>("ForceCollide")
        .constructor<>()
        .function("force", &ForceCollide::force)
        .function("setNodes", &ForceCollide::setNodes)
        .function("setNodeBuffer", &ForceCollide::setNodeBuffer)
        .function("setRadiusBuffer", &ForceCollide::setRadiusBuffer)
        .function("setRadius", &ForceCollide::setRadius)
        .function("getRadius", &ForceCollide::getRadius)
        .function("setStrength", &ForceCollide::setStrength)
        .function("getStrength", &ForceCollide::getStrength)
        .function("setIterations", &ForceCollide::setIterations)
        .function("getIterations", &ForceCollide::getIterations);

    emscripten::function("createForceCollide", &createForceCollide, emscripten::allow_raw_pointers());
}
