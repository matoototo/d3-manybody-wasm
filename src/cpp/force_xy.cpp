#include <emscripten/bind.h>
#include <vector>
#include <functional>
#include <cmath>
#include <cstdint>
#include <algorithm>

class ForceXY {
protected:
    std::function<double(const emscripten::val&, int, const emscripten::val&)> coordinate;
    std::function<double(const emscripten::val&, int, const emscripten::val&)> strength;
    std::vector<float> strengths;
    std::vector<float> coordz;
    emscripten::val nodes;
    double* nodeBuffer = nullptr;
    int nodeCount = 0;
    static constexpr int stride = 4;

public:
    ForceXY(const emscripten::val& coord) :
        coordinate([coord](const emscripten::val&, int, const emscripten::val&) { return coord.as<double>(); }),
        strength([](const emscripten::val&, int, const emscripten::val&) { return 0.1; }) {}

    virtual ~ForceXY() {}

    void force(double alpha) {
        if (nodeBuffer != nullptr && nodeCount > 0) {
            for (int i = 0; i < nodeCount; ++i) {
                updateNodeVelocityBuffer(i, alpha);
            }
            return;
        }

        if (nodes.isUndefined()) return;
        int n = nodes["length"].as<int>();
        for (int i = 0; i < n; ++i) {
            emscripten::val node = nodes[i];
            updateNodeVelocity(node, i, alpha);
        }
    }

    virtual void updateNodeVelocity(emscripten::val& node, int i, double alpha) = 0;
    virtual void updateNodeVelocityBuffer(int index, double alpha) = 0;

    void initialize() {
        if (nodes.isUndefined()) return;
        int n = nodes["length"].as<int>();
        nodeCount = n;
        strengths.resize(n);
        coordz.resize(n);
        for (int i = 0; i < n; ++i) {
            double cval = coordinate(nodes[i], i, nodes);
            double sval = strength(nodes[i], i, nodes);
            coordz[i] = static_cast<float>(cval);
            strengths[i] = std::isnan(cval) ? 0.0f : static_cast<float>(sval);
        }
    }

    void setNodes(const emscripten::val& _nodes) {
        nodes = _nodes;
        if (!nodes.isUndefined()) {
            nodeCount = nodes["length"].as<int>();
        } else {
            nodeCount = 0;
        }
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

    void setNodeBuffer(uintptr_t ptr, int count) {
        nodeBuffer = (ptr != 0 && count > 0) ? reinterpret_cast<double*>(ptr) : nullptr;
        nodeCount = count;
    }
};

class ForceX : public ForceXY {
public:
    ForceX(const emscripten::val& x) : ForceXY(x) {}

    void updateNodeVelocity(emscripten::val& node, int i, double alpha) override {
        node.set("vx", node["vx"].as<double>() + (coordz[i] - node["x"].as<double>()) * strengths[i] * alpha);
    }

    void updateNodeVelocityBuffer(int index, double alpha) override {
        if (nodeBuffer == nullptr) return;
        float coord = coordz[index];
        float strengthValue = strengths[index];
        if (strengthValue == 0 || std::isnan(coord)) return;
        int offset = index * stride;
        float currentX = static_cast<float>(nodeBuffer[offset]);
        float currentVelocity = static_cast<float>(nodeBuffer[offset + 2]);
        float delta = static_cast<float>((coord - currentX) * strengthValue * alpha);
        nodeBuffer[offset + 2] = static_cast<float>(currentVelocity + delta);
    }
};

class ForceY : public ForceXY {
public:
    ForceY(const emscripten::val& y) : ForceXY(y) {}

    void updateNodeVelocity(emscripten::val& node, int i, double alpha) override {
        node.set("vy", node["vy"].as<double>() + (coordz[i] - node["y"].as<double>()) * strengths[i] * alpha);
    }

    void updateNodeVelocityBuffer(int index, double alpha) override {
        if (nodeBuffer == nullptr) return;
        float coord = coordz[index];
        float strengthValue = strengths[index];
        if (strengthValue == 0 || std::isnan(coord)) return;
        int offset = index * stride;
        float currentY = static_cast<float>(nodeBuffer[offset + 1]);
        float currentVelocity = static_cast<float>(nodeBuffer[offset + 3]);
        float delta = static_cast<float>((coord - currentY) * strengthValue * alpha);
        nodeBuffer[offset + 3] = static_cast<float>(currentVelocity + delta);
    }
};

ForceX* createForceX(const emscripten::val& x) {
    return new ForceX(x);
}

ForceY* createForceY(const emscripten::val& y) {
    return new ForceY(y);
}

void integrateNodeBuffer(uintptr_t pointer, int count, double velocityRetention) {
    if (pointer == 0 || count <= 0) return;
    double* nodeBuffer = reinterpret_cast<double*>(pointer);
    for (int index = 0; index < count; ++index) {
        const int offset = index * 4;
        nodeBuffer[offset + 2] *= velocityRetention;
        nodeBuffer[offset + 3] *= velocityRetention;
        nodeBuffer[offset] += nodeBuffer[offset + 2];
        nodeBuffer[offset + 1] += nodeBuffer[offset + 3];
    }
}

struct ReplayAcceleration {
    float x, y;
};

void replayTwoFieldsAndIntegrate(
    uintptr_t nodePointer,
    int count,
    double alphaValue,
    double velocityRetention,
    uintptr_t firstAccelerationPointer,
    uintptr_t firstTrendPointer,
    int firstAge,
    uintptr_t secondAccelerationPointer,
    uintptr_t secondTrendPointer,
    int secondAge
) {
    if (nodePointer == 0 || count <= 0) return;
    double* nodeBuffer = reinterpret_cast<double*>(nodePointer);
    const ReplayAcceleration* firstAcceleration = reinterpret_cast<ReplayAcceleration*>(firstAccelerationPointer);
    const ReplayAcceleration* firstTrend = reinterpret_cast<ReplayAcceleration*>(firstTrendPointer);
    const ReplayAcceleration* secondAcceleration = reinterpret_cast<ReplayAcceleration*>(secondAccelerationPointer);
    const ReplayAcceleration* secondTrend = reinterpret_cast<ReplayAcceleration*>(secondTrendPointer);
    const float alpha = static_cast<float>(alphaValue);
    const int firstPredictionAge = std::min(firstAge, 4);
    const int secondPredictionAge = std::min(secondAge, 4);
    for (int index = 0; index < count; ++index) {
        const int offset = index * 4;
        nodeBuffer[offset + 2] += (
            firstAcceleration[index].x + firstTrend[index].x * firstPredictionAge
        ) * alpha;
        nodeBuffer[offset + 3] += (
            firstAcceleration[index].y + firstTrend[index].y * firstPredictionAge
        ) * alpha;
        nodeBuffer[offset + 2] += (
            secondAcceleration[index].x + secondTrend[index].x * secondPredictionAge
        ) * alpha;
        nodeBuffer[offset + 3] += (
            secondAcceleration[index].y + secondTrend[index].y * secondPredictionAge
        ) * alpha;
        nodeBuffer[offset + 2] *= velocityRetention;
        nodeBuffer[offset + 3] *= velocityRetention;
        nodeBuffer[offset] += nodeBuffer[offset + 2];
        nodeBuffer[offset + 1] += nodeBuffer[offset + 3];
    }
}

EMSCRIPTEN_BINDINGS(force_xy_module) {
    emscripten::class_<ForceXY>("ForceXY")
        .function("force", &ForceXY::force)
        .function("initialize", &ForceXY::initialize)
        .function("setNodes", &ForceXY::setNodes)
        .function("setNodeBuffer", &ForceXY::setNodeBuffer)
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
    emscripten::function("integrateNodeBuffer", &integrateNodeBuffer);
    emscripten::function("replayTwoFieldsAndIntegrate", &replayTwoFieldsAndIntegrate);
}
