#include <emscripten/bind.h>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>

class ForceManyBody {
private:
    emscripten::val nodes;
    std::vector<double> strengths;
    std::function<double(const emscripten::val&, int, const emscripten::val&)> strength;
    double distanceMin2 = 1;
    double distanceMax2 = std::numeric_limits<double>::infinity();
    double theta2 = 0.81;
    double alpha;

    // Random number generator for jiggle
    std::function<double()> random;

    // Quadtree node structure
    struct QuadtreeNode {
        double x0, y0, x1, y1; // Bounds
        double value = 0;      // Total strength
        double x = 0, y = 0;   // Center of mass
        QuadtreeNode* children[4] = { nullptr, nullptr, nullptr, nullptr };
        emscripten::val data = emscripten::val::undefined(); // For leaf nodes
        QuadtreeNode* next = nullptr; // For linked list of coincident nodes

        ~QuadtreeNode() {
            for (int i = 0; i < 4; ++i) {
                delete children[i];
            }
        }
    };

    QuadtreeNode* quadtreeRoot = nullptr;

    void findExtent(double& x0, double& y0, double& x1, double& y1) {
        int n = nodes["length"].as<int>();
        x0 = y0 = std::numeric_limits<double>::infinity();
        x1 = y1 = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            double xi = node["x"].as<double>();
            double yi = node["y"].as<double>();
            if (xi < x0) x0 = xi;
            if (xi > x1) x1 = xi;
            if (yi < y0) y0 = yi;
            if (yi > y1) y1 = yi;
        }
        // Expand the bounds slightly to avoid precision errors
        double dx = x1 - x0;
        double dy = y1 - y0;
        if (dx == 0) dx = 1;
        if (dy == 0) dy = 1;
        x0 -= dx * 0.1;
        x1 += dx * 0.1;
        y0 -= dy * 0.1;
        y1 += dy * 0.1;
    }

    void insertNode(QuadtreeNode* node, emscripten::val& dataNode) {
        double x = dataNode["x"].as<double>();
        double y = dataNode["y"].as<double>();

        // Leaf node
        if (node->data.isUndefined()) {
            node->data = dataNode;
            return;
        }

        // Internal node or leaf with data, need to split
        if (node->children[0] == nullptr) {
            // Subdivide node
            subdivideNode(node);
            // Re-insert existing data node
            emscripten::val existingDataNode = node->data;
            node->data = emscripten::val::undefined();
            insertNode(node, existingDataNode);
        }

        // Determine which quadrant to insert the node into
        int i = getQuadrant(node, x, y);
        insertNode(node->children[i], dataNode);
    }

    void subdivideNode(QuadtreeNode* node) {
        double x0 = node->x0;
        double y0 = node->y0;
        double x1 = node->x1;
        double y1 = node->y1;
        double xm = (x0 + x1) / 2;
        double ym = (y0 + y1) / 2;

        node->children[0] = new QuadtreeNode(); // NW
        node->children[0]->x0 = x0;
        node->children[0]->y0 = y0;
        node->children[0]->x1 = xm;
        node->children[0]->y1 = ym;

        node->children[1] = new QuadtreeNode(); // NE
        node->children[1]->x0 = xm;
        node->children[1]->y0 = y0;
        node->children[1]->x1 = x1;
        node->children[1]->y1 = ym;

        node->children[2] = new QuadtreeNode(); // SW
        node->children[2]->x0 = x0;
        node->children[2]->y0 = ym;
        node->children[2]->x1 = xm;
        node->children[2]->y1 = y1;

        node->children[3] = new QuadtreeNode(); // SE
        node->children[3]->x0 = xm;
        node->children[3]->y0 = ym;
        node->children[3]->x1 = x1;
        node->children[3]->y1 = y1;
    }

    int getQuadrant(QuadtreeNode* node, double x, double y) {
        double xm = (node->x0 + node->x1) / 2;
        double ym = (node->y0 + node->y1) / 2;
        return (y >= ym ? 2 : 0) + (x >= xm ? 1 : 0);
    }

    void buildQuadtree() {
        double x0, y0, x1, y1;
        findExtent(x0, y0, x1, y1);
        quadtreeRoot = new QuadtreeNode();
        quadtreeRoot->x0 = x0;
        quadtreeRoot->y0 = y0;
        quadtreeRoot->x1 = x1;
        quadtreeRoot->y1 = y1;

        int n = nodes["length"].as<int>();
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            insertNode(quadtreeRoot, node);
        }

        accumulate(quadtreeRoot);
    }

    void accumulate(QuadtreeNode* node) {
        double strength = 0.0;
        double x = 0.0, y = 0.0;
        double weight = 0.0;

        if (node->children[0] != nullptr) {
            // Internal node
            for (int i = 0; i < 4; ++i) {
                QuadtreeNode* child = node->children[i];
                if (child) {
                    accumulate(child);
                    double c = std::abs(child->value);
                    if (c > 0) {
                        strength += child->value;
                        weight += c;
                        x += c * child->x;
                        y += c * child->y;
                    }
                }
            }
            node->x = x / weight;
            node->y = y / weight;
        } else if (!node->data.isUndefined()) {
            // Leaf node
            node->x = node->data["x"].as<double>();
            node->y = node->data["y"].as<double>();
            int index = node->data["index"].as<int>();
            strength += strengths[index];
            // Handle coincident nodes if any
            QuadtreeNode* q = node->next;
            while (q) {
                int idx = q->data["index"].as<int>();
                strength += strengths[idx];
                q = q->next;
            }
        }
        node->value = strength;
    }

    void apply(QuadtreeNode* quad, emscripten::val& node) {
        if (quad->value == 0) return;

        double dx = quad->x - node["x"].as<double>();
        double dy = quad->y - node["y"].as<double>();
        double w = quad->x1 - quad->x0;
        double l = dx * dx + dy * dy;

        // Apply the Barnes-Hut approximation if possible.
        if (w * w / theta2 < l) {
            if (l < distanceMax2) {
                if (dx == 0) dx = jiggle(), l += dx * dx;
                if (dy == 0) dy = jiggle(), l += dy * dy;
                if (l < distanceMin2) l = distanceMin2;
                double factor = quad->value * alpha / l;
                node.set("vx", node["vx"].as<double>() + dx * factor);
                node.set("vy", node["vy"].as<double>() + dy * factor);
            }
        } else if (quad->children[0] != nullptr) {
            // Internal node, recurse into children
            for (int i = 0; i < 4; ++i) {
                if (quad->children[i]) {
                    apply(quad->children[i], node);
                }
            }
        } else if (!quad->data.isUndefined() && l < distanceMax2) {
            // Leaf node
            if (quad->data != node || quad->next != nullptr) {
                if (dx == 0) dx = jiggle(), l += dx * dx;
                if (dy == 0) dy = jiggle(), l += dy * dy;
                if (l < distanceMin2) l = distanceMin2;
            }
            QuadtreeNode* q = quad;
            do {
                if (q->data != node) {
                    int index = q->data["index"].as<int>();
                    double w = strengths[index] * alpha / l;
                    node.set("vx", node["vx"].as<double>() + dx * w);
                    node.set("vy", node["vy"].as<double>() + dy * w);
                }
                q = q->next;
            } while (q != nullptr);
        }
    }

    double jiggle() {
        return (random() - 0.5) * 1e-6;
    }

public:
    ForceManyBody()
        : strength([](const emscripten::val&, int, const emscripten::val&) { return -30; }),
          random([]() { return std::rand() / double(RAND_MAX); })
    {}

    void force(double alpha_) {
        alpha = alpha_;
        buildQuadtree();
        int n = nodes["length"].as<int>();
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            apply(quadtreeRoot, node);
        }
        // Clean up
        delete quadtreeRoot;
        quadtreeRoot = nullptr;
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

    void setRandom(const emscripten::val& _random) {
        random = _random.as<std::function<double()>>();
    }

    emscripten::val getRandom() const {
        return emscripten::val::global("Math")["random"];
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
        .function("getTheta", &ForceManyBody::getTheta)
        .function("setRandom", &ForceManyBody::setRandom)
        .function("getRandom", &ForceManyBody::getRandom);

    emscripten::function("createForceManyBody", &createForceManyBody, emscripten::allow_raw_pointers());
}
