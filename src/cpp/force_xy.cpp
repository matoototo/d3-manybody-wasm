#include <emscripten/bind.h>
#include <vector>
#include <functional>
#include <cmath>

class ForceXY {
protected:
    std::function<double(const emscripten::val&, int, const emscripten::val&)> coordinate;
    std::function<double(const emscripten::val&, int, const emscripten::val&)> strength;
    std::vector<double> strengths;
    std::vector<double> coordz;
    emscripten::val nodes;

public:
    ForceXY(const emscripten::val& coord) : 
        coordinate([coord](const emscripten::val&, int, const emscripten::val&) { return coord.as<double>(); }),
        strength([](const emscripten::val&, int, const emscripten::val&) { return 0.1; }) {}

    virtual ~ForceXY() {} // Virtual destructor

    void force(double alpha) {
        int n = nodes["length"].as<int>();
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            updateNodeVelocity(node, i, alpha);
        }
    }

    virtual void updateNodeVelocity(emscripten::val& node, int i, double alpha) = 0;

    void initialize() {
        if (nodes.isUndefined()) return;
        int n = nodes["length"].as<int>();
        strengths.resize(n);
        coordz.resize(n);
        for (int i = 0; i < n; ++i) {
            coordz[i] = coordinate(nodes[i], i, nodes);
            strengths[i] = std::isnan(coordz[i]) ? 0 : strength(nodes[i], i, nodes);
        }
    }

    void setNodes(const emscripten::val& _nodes) {
        nodes = _nodes;
        initialize();
    }

    emscripten::val getStrength() const {
        return emscripten::val(strength);
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

    emscripten::val getCoordinate() const {
        return emscripten::val(coordinate);
    }

    void setCoordinate(const emscripten::val& _coordinate) {
        if (_coordinate.typeOf().as<std::string>() == "function") {
            coordinate = _coordinate.as<std::function<double(const emscripten::val&, int, const emscripten::val&)>>();
        } else {
            double c = _coordinate.as<double>();
            coordinate = [c](const emscripten::val&, int, const emscripten::val&) { return c; };
        }
        initialize();
    }
};

class ForceX : public ForceXY {
public:
    ForceX(const emscripten::val& x) : ForceXY(x) {}

    void updateNodeVelocity(emscripten::val& node, int i, double alpha) override {
        node.set("vx", node["vx"].as<double>() + (coordz[i] - node["x"].as<double>()) * strengths[i] * alpha);
    }
};

class ForceY : public ForceXY {
public:
    ForceY(const emscripten::val& y) : ForceXY(y) {}

    void updateNodeVelocity(emscripten::val& node, int i, double alpha) override {
        node.set("vy", node["vy"].as<double>() + (coordz[i] - node["y"].as<double>()) * strengths[i] * alpha);
    }
};

// Factory functions
ForceX* createForceX(const emscripten::val& x) {
    return new ForceX(x);
}

ForceY* createForceY(const emscripten::val& y) {
    return new ForceY(y);
}

// Bindings
EMSCRIPTEN_BINDINGS(force_xy_module) {
    emscripten::class_<ForceXY>("ForceXY")
        .function("force", &ForceXY::force)
        .function("initialize", &ForceXY::initialize)
        .function("setNodes", &ForceXY::setNodes)
        .function("getStrength", &ForceXY::getStrength)
        .function("setStrength", &ForceXY::setStrength)
        .function("getCoordinate", &ForceXY::getCoordinate)
        .function("setCoordinate", &ForceXY::setCoordinate);

    emscripten::class_<ForceX, emscripten::base<ForceXY>>("ForceX")
        .constructor<emscripten::val>();

    emscripten::class_<ForceY, emscripten::base<ForceXY>>("ForceY")
        .constructor<emscripten::val>();

    emscripten::function("createForceX", &createForceX, emscripten::allow_raw_pointers());
    emscripten::function("createForceY", &createForceY, emscripten::allow_raw_pointers());
}