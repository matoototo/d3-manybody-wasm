#include <emscripten/bind.h>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>

// Function to calculate the Morton code (Z-order curve) for 2D coordinates
unsigned int mortonCode(double x, double y, double x0, double y0, double s, unsigned int scaleFactor) {
    unsigned int scaledX = static_cast<unsigned int>((x - x0) / s * scaleFactor);
    unsigned int scaledY = static_cast<unsigned int>((y - y0) / s * scaleFactor);

    unsigned int morton = 0;
    for (unsigned int i = 0; i < sizeof(unsigned int) * 8 / 2; ++i) {
        morton |= ((scaledX & (1 << i)) << i) | ((scaledY & (1 << i)) << (i + 1));
    }

    return morton;
}

class ForceManyBody {
private:
    emscripten::val nodes;
    std::vector<double> strengths;
    std::function<double(const emscripten::val&, int, const emscripten::val&)> strength;
    double distanceMin2 = 1;
    double distanceMax2 = std::numeric_limits<double>::infinity();
    double theta2 = 0.9;
    double alpha;

    struct BodyData {
        double x, y, vx, vy;
    };
    std::vector<BodyData> bodyData;

    std::function<double()> random;

    struct QuadtreeNode {
        double cx, cy;  // Center of the node
        double s;       // Half-size of the region
        double value = 0;  // Total mass
        double x = 0, y = 0;  // Center of mass
        int firstChild = -1;  // Index of first child, or -1 if leaf node

        QuadtreeNode(double cx_, double cy_, double s_) : cx(cx_), cy(cy_), s(s_) {}
    };

    std::vector<QuadtreeNode> quadtreeNodes;

    // Find bounding box of all nodes
    void findExtent(double& x0, double& y0, double& x1, double& y1) {
        int n = nodes["length"].as<int>();
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
    void buildQuadtree() {
        double x0, y0, x1, y1;
        findExtent(x0, y0, x1, y1);
        double cx = (x0 + x1) / 2;
        double cy = (y0 + y1) / 2;
        double s = std::max(x1 - x0, y1 - y0) / 2 * 1.1;

        quadtreeNodes.clear();
        quadtreeNodes.emplace_back(cx, cy, s);

        int n = nodes["length"].as<int>();
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
            for (int j = 0; j < 4; ++j) {
                QuadtreeNode& child = quadtreeNodes[node.firstChild + j];
                mass += child.value;
                x += child.x * child.value;
                y += child.y * child.value;
            }

            if (mass > 0) {
                node.x = x / mass;
                node.y = y / mass;
            }
            node.value = mass;
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
                body.vx += dx * factor;
                body.vy += dy * factor;
            }
        } else if (quad.firstChild >= 0) {
            for (int i = 0; i < 4; ++i) {
                apply(quad.firstChild + i, body);
            }
        }
    }

public:
    ForceManyBody()
        : strength([](const emscripten::val&, int, const emscripten::val&) { return -30; }),
          random([]() { return std::rand() / double(RAND_MAX); })
    {}

    void force(double alpha_) {
        alpha = alpha_;

        // Sync data from JavaScript to C++
        int n = nodes["length"].as<int>();
        bodyData.resize(n);
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            bodyData[i].x = node["x"].as<double>();
            bodyData[i].y = node["y"].as<double>();
            bodyData[i].vx = node["vx"].as<double>();
            bodyData[i].vy = node["vy"].as<double>();
        }

        // Calculate Morton codes and sort nodes
        double x0, y0, x1, y1;
        findExtent(x0, y0, x1, y1);
        double s = std::max(x1 - x0, y1 - y0);
        unsigned int scaleFactor = 65536;

        std::vector<std::pair<unsigned int, int>> sortedIndices(n);
        for (int i = 0; i < n; ++i) {
            sortedIndices[i] = {mortonCode(bodyData[i].x, bodyData[i].y, x0, y0, s, scaleFactor), i};
        }

        std::sort(sortedIndices.begin(), sortedIndices.end());

        // Reorder bodyData based on Morton code
        std::vector<BodyData> newBodyData(n);
        for (int i = 0; i < n; ++i) {
            int oldIndex = sortedIndices[i].second;
            newBodyData[i] = bodyData[oldIndex];
        }

        bodyData = newBodyData;

        buildQuadtree();
        propagate();

        // Apply forces using C++ data
        for (int i = 0; i < n; ++i) {
            apply(0, bodyData[i]);
        }

        // Sync data back to JavaScript
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[sortedIndices[i].second];
            node.set("vx", bodyData[i].vx);
            node.set("vy", bodyData[i].vy);
        }

        quadtreeNodes.clear();
    }

    void initialize() {
        if (nodes.isUndefined()) return;
        int n = nodes["length"].as<int>();
        strengths.resize(n);
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            strengths[i] = strength(node, i, nodes);
        }
    }

    void setNodes(const emscripten::val& _nodes) {
        nodes = _nodes;
        int n = nodes["length"].as<int>();
        bodyData.resize(n);
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            bodyData[i].x = node["x"].as<double>();
            bodyData[i].y = node["y"].as<double>();
            bodyData[i].vx = node["vx"].as<double>();
            bodyData[i].vy = node["vy"].as<double>();
        }
        initialize();
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

    emscripten::val getStrength() const {
        return emscripten::val(strength);
    }

    void setDistanceMin(double d) {
        distanceMin2 = d * d;
    }

    double getDistanceMin() const {
        return std::sqrt(distanceMin2);
    }

    void setDistanceMax(double d) {
        distanceMax2 = d * d;
    }

    double getDistanceMax() const {
        return std::sqrt(distanceMax2);
    }

    void setTheta(double t) {
        theta2 = t * t;
    }

    double getTheta() const {
        return std::sqrt(theta2);
    }
};

// Factory function
ForceManyBody* createForceManyBody() {
    return new ForceManyBody();
}

// Bindings
EMSCRIPTEN_BINDINGS(force_many_body_module) {
    emscripten::class_<ForceManyBody>("ForceManyBody")
        .constructor<>()
        .function("force", &ForceManyBody::force)
        .function("setNodes", &ForceManyBody::setNodes)
        .function("setStrength", &ForceManyBody::setStrength)
        .function("getStrength", &ForceManyBody::getStrength)
        .function("setDistanceMin", &ForceManyBody::setDistanceMin)
        .function("getDistanceMin", &ForceManyBody::getDistanceMin)
        .function("setDistanceMax", &ForceManyBody::setDistanceMax)
        .function("getDistanceMax", &ForceManyBody::getDistanceMax)
        .function("setTheta", &ForceManyBody::setTheta)
        .function("getTheta", &ForceManyBody::getTheta);

    emscripten::function("createForceManyBody", &createForceManyBody, emscripten::allow_raw_pointers());
}