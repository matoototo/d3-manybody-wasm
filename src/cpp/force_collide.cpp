#include <emscripten/bind.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr int kStride = 4;

struct Quad {
    std::array<int, 4> children = {-1, -1, -1, -1};
    int data = -1;
    int next = -1;
    double radius = 0;

    bool isLeaf() const { return data >= 0; }
};

}  // namespace

class ForceCollide {
private:
    emscripten::val nodes = emscripten::val::undefined();
    std::function<double(const emscripten::val&, int, const emscripten::val&)> radiusFunction =
        [](const emscripten::val&, int, const emscripten::val&) { return 1.0; };
    bool customRadius = false;
    float radiusConstant = 1.0f;
    double strength = 1.0;
    int iterations = 1;

    std::vector<float> radii;
    std::vector<double> positionsX;
    std::vector<double> positionsY;
    std::vector<double> treeX;
    std::vector<double> treeY;
    std::vector<double> velocitiesX;
    std::vector<double> velocitiesY;
    std::vector<Quad> quads;
    double* nodeBuffer = nullptr;
    float* radiusBuffer = nullptr;
    int nodeCount = 0;
    int radiusBufferCount = 0;
    uint32_t randomState = 1;

    double jiggle() {
        randomState = 1664525u * randomState + 1013904223u;
        return (static_cast<double>(randomState) / 4294967296.0 - 0.5) * 1e-6;
    }

    void computeRadii() {
        if (radiusBuffer && radiusBufferCount >= nodeCount) return;
        radii.resize(nodeCount);
        if (customRadius && !nodes.isUndefined() && !nodes.isNull()) {
            for (int index = 0; index < nodeCount; ++index) {
                radii[index] = static_cast<float>(radiusFunction(nodes[index], index, nodes));
            }
        } else {
            std::fill(radii.begin(), radii.end(), radiusConstant);
        }
    }

    const float* radiusValues() {
        if (radiusBuffer && radiusBufferCount >= nodeCount) return radiusBuffer;
        if (static_cast<int>(radii.size()) != nodeCount) computeRadii();
        return radii.data();
    }

    int makeLeaf(int data, int next = -1) {
        quads.emplace_back();
        quads.back().data = data;
        quads.back().next = next;
        return static_cast<int>(quads.size()) - 1;
    }

    int makeInternal() {
        quads.emplace_back();
        return static_cast<int>(quads.size()) - 1;
    }

    void insert(int dataIndex, int& root, double x0, double y0, double x1, double y1) {
        const double x = treeX[dataIndex];
        const double y = treeY[dataIndex];
        if (!std::isfinite(x) || !std::isfinite(y)) return;

        const int leaf = makeLeaf(dataIndex);
        if (root < 0) {
            root = leaf;
            return;
        }

        int parent = -1;
        int parentQuadrant = -1;
        int node = root;
        while (!quads[node].isLeaf()) {
            const double xm = (x0 + x1) / 2;
            const double ym = (y0 + y1) / 2;
            const bool right = x >= xm;
            const bool bottom = y >= ym;
            const int quadrant = (bottom ? 2 : 0) | (right ? 1 : 0);
            if (right) x0 = xm; else x1 = xm;
            if (bottom) y0 = ym; else y1 = ym;
            parent = node;
            parentQuadrant = quadrant;
            const int child = quads[node].children[quadrant];
            if (child < 0) {
                quads[node].children[quadrant] = leaf;
                return;
            }
            node = child;
        }

        const int oldData = quads[node].data;
        const double oldX = treeX[oldData];
        const double oldY = treeY[oldData];
        if (x == oldX && y == oldY) {
            quads[leaf].next = node;
            if (parent < 0) root = leaf;
            else quads[parent].children[parentQuadrant] = leaf;
            return;
        }

        int branchParent = parent;
        int branchQuadrant = parentQuadrant;
        while (true) {
            const int branch = makeInternal();
            if (branchParent < 0) root = branch;
            else quads[branchParent].children[branchQuadrant] = branch;

            const double xm = (x0 + x1) / 2;
            const double ym = (y0 + y1) / 2;
            const bool right = x >= xm;
            const bool bottom = y >= ym;
            const int newQuadrant = (bottom ? 2 : 0) | (right ? 1 : 0);
            const int oldQuadrant = (oldY >= ym ? 2 : 0) | (oldX >= xm ? 1 : 0);
            if (newQuadrant != oldQuadrant) {
                quads[branch].children[oldQuadrant] = node;
                quads[branch].children[newQuadrant] = leaf;
                return;
            }

            if (right) x0 = xm; else x1 = xm;
            if (bottom) y0 = ym; else y1 = ym;
            branchParent = branch;
            branchQuadrant = newQuadrant;
        }
    }

    void prepare(int node) {
        Quad& quad = quads[node];
        if (quad.isLeaf()) {
            quad.radius = radii[quad.data];
            return;
        }
        quad.radius = 0;
        for (int quadrant = 0; quadrant < 4; ++quadrant) {
            const int child = quad.children[quadrant];
            if (child < 0) continue;
            prepare(child);
            quad.radius = std::max(quad.radius, quads[child].radius);
        }
    }

    int buildTree(double& x0, double& y0, double& x1, double& y1) {
        double minX = std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();
        for (int index = 0; index < nodeCount; ++index) {
            const double x = treeX[index];
            const double y = treeY[index];
            if (!std::isfinite(x) || !std::isfinite(y)) continue;
            minX = std::min(minX, x);
            minY = std::min(minY, y);
            maxX = std::max(maxX, x);
            maxY = std::max(maxY, y);
        }
        if (minX > maxX || minY > maxY) return -1;

        x0 = std::floor(minX);
        y0 = std::floor(minY);
        double size = 1;
        while (maxX >= x0 + size || maxY >= y0 + size) size *= 2;
        x1 = x0 + size;
        y1 = y0 + size;

        quads.clear();
        quads.reserve(static_cast<size_t>(nodeCount) * 2);
        int root = -1;
        for (int index = 0; index < nodeCount; ++index) {
            insert(index, root, x0, y0, x1, y1);
        }
        if (root >= 0) prepare(root);
        return root;
    }

    void applyNode(int left, int root, double x0, double y0, double x1, double y1, const float* radius) {
        struct Visit {
            int node;
            double x0;
            double y0;
            double x1;
            double y1;
        };
        std::vector<Visit> stack;
        stack.reserve(64);
        stack.push_back({root, x0, y0, x1, y1});

        const double ri = radius[left];
        const double ri2 = ri * ri;
        double xi = positionsX[left] + velocitiesX[left];
        double yi = positionsY[left] + velocitiesY[left];

        while (!stack.empty()) {
            const Visit visit = stack.back();
            stack.pop_back();
            const Quad& quad = quads[visit.node];
            const double rjInitial = quad.radius;
            const double combined = ri + rjInitial;

            if (quad.isLeaf()) {
                const int right = quad.data;
                if (right > left) {
                    double dx = xi - positionsX[right] - velocitiesX[right];
                    double dy = yi - positionsY[right] - velocitiesY[right];
                    double distance2 = dx * dx + dy * dy;
                    if (distance2 < combined * combined) {
                        if (dx == 0) {
                            dx = jiggle();
                            distance2 += dx * dx;
                        }
                        if (dy == 0) {
                            dy = jiggle();
                            distance2 += dy * dy;
                        }
                        double distance = std::sqrt(distance2);
                        const double overlap = (combined - distance) / distance * strength;
                        dx *= overlap;
                        dy *= overlap;
                        const double rj2 = rjInitial * rjInitial;
                        const double leftWeight = rj2 / (ri2 + rj2);
                        velocitiesX[left] += dx * leftWeight;
                        velocitiesY[left] += dy * leftWeight;
                        velocitiesX[right] -= dx * (1 - leftWeight);
                        velocitiesY[right] -= dy * (1 - leftWeight);
                        xi = positionsX[left] + velocitiesX[left];
                        yi = positionsY[left] + velocitiesY[left];
                    }
                }
                continue;
            }

            if (visit.x0 > xi + combined || visit.x1 < xi - combined ||
                visit.y0 > yi + combined || visit.y1 < yi - combined) {
                continue;
            }

            const double xm = (visit.x0 + visit.x1) / 2;
            const double ym = (visit.y0 + visit.y1) / 2;
            for (int quadrant = 3; quadrant >= 0; --quadrant) {
                const int child = quad.children[quadrant];
                if (child < 0) continue;
                const bool right = (quadrant & 1) != 0;
                const bool bottom = (quadrant & 2) != 0;
                stack.push_back({
                    child,
                    right ? xm : visit.x0,
                    bottom ? ym : visit.y0,
                    right ? visit.x1 : xm,
                    bottom ? visit.y1 : ym
                });
            }
        }
    }

public:
    void force(double /*alpha*/) {
        if (!nodeBuffer || nodeCount <= 0) return;
        const float* radius = radiusValues();
        positionsX.resize(nodeCount);
        positionsY.resize(nodeCount);
        treeX.resize(nodeCount);
        treeY.resize(nodeCount);
        velocitiesX.resize(nodeCount);
        velocitiesY.resize(nodeCount);

        for (int index = 0; index < nodeCount; ++index) {
            const int offset = index * kStride;
            positionsX[index] = nodeBuffer[offset];
            positionsY[index] = nodeBuffer[offset + 1];
            velocitiesX[index] = nodeBuffer[offset + 2];
            velocitiesY[index] = nodeBuffer[offset + 3];
        }

        for (int iteration = 0; iteration < iterations; ++iteration) {
            for (int index = 0; index < nodeCount; ++index) {
                treeX[index] = positionsX[index] + velocitiesX[index];
                treeY[index] = positionsY[index] + velocitiesY[index];
            }
            double x0, y0, x1, y1;
            const int root = buildTree(x0, y0, x1, y1);
            if (root < 0) break;
            for (int index = 0; index < nodeCount; ++index) {
                applyNode(index, root, x0, y0, x1, y1, radius);
            }
        }

        for (int index = 0; index < nodeCount; ++index) {
            const int offset = index * kStride;
            nodeBuffer[offset + 2] = static_cast<float>(velocitiesX[index]);
            nodeBuffer[offset + 3] = static_cast<float>(velocitiesY[index]);
        }
    }

    void setNodes(const emscripten::val& value) {
        nodes = value;
        nodeCount = nodes.isUndefined() || nodes.isNull() ? 0 : nodes["length"].as<int>();
        computeRadii();
    }

    void setNodeBuffer(uintptr_t pointer, int count) {
        nodeBuffer = pointer && count > 0 ? reinterpret_cast<double*>(pointer) : nullptr;
        nodeCount = count;
    }

    void setRadiusBuffer(uintptr_t pointer, int count) {
        radiusBuffer = pointer && count > 0 ? reinterpret_cast<float*>(pointer) : nullptr;
        radiusBufferCount = count;
    }

    void setRadii(const emscripten::val& values) {
        radiusBuffer = nullptr;
        radiusBufferCount = 0;
        const int count = values["length"].as<int>();
        radii.resize(count);
        for (int index = 0; index < count; ++index) radii[index] = values[index].as<float>();
    }

    void setRadius(const emscripten::val& value) {
        if (value.typeOf().as<std::string>() == "function") {
            radiusFunction = value.as<std::function<double(const emscripten::val&, int, const emscripten::val&)>>();
            customRadius = true;
        } else {
            radiusConstant = static_cast<float>(value.as<double>());
            customRadius = false;
        }
        computeRadii();
    }

    emscripten::val getRadius() const {
        return customRadius ? emscripten::val::undefined() : emscripten::val(radiusConstant);
    }

    void setStrength(double value) { strength = value; }
    double getStrength() const { return strength; }
    void setIterations(int value) { iterations = std::max(1, value); }
    int getIterations() const { return iterations; }
};

ForceCollide* createForceCollide() { return new ForceCollide(); }

EMSCRIPTEN_BINDINGS(force_collide_module) {
    emscripten::class_<ForceCollide>("ForceCollide")
        .constructor<>()
        .function("force", &ForceCollide::force)
        .function("setNodes", &ForceCollide::setNodes)
        .function("setNodeBuffer", &ForceCollide::setNodeBuffer)
        .function("setRadiusBuffer", &ForceCollide::setRadiusBuffer)
        .function("setRadii", &ForceCollide::setRadii)
        .function("setRadius", &ForceCollide::setRadius)
        .function("getRadius", &ForceCollide::getRadius)
        .function("setStrength", &ForceCollide::setStrength)
        .function("getStrength", &ForceCollide::getStrength)
        .function("setIterations", &ForceCollide::setIterations)
        .function("getIterations", &ForceCollide::getIterations);

    emscripten::function("createForceCollide", &createForceCollide, emscripten::allow_raw_pointers());
}
