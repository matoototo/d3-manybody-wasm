#include <emscripten/bind.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

constexpr int kStride = 4;

}  // namespace

class ForceLink {
private:
    std::vector<int> sources;
    std::vector<int> targets;
    std::vector<double> biases;
    std::vector<double> strengths;
    std::vector<double> distances;

    struct Velocity {
        double x, y;
    };
    struct Acceleration {
        float x, y;
    };
    std::vector<Velocity> initialVelocity;
    std::vector<Acceleration> cachedAcceleration;
    std::vector<Acceleration> accelerationTrend;
    bool cachedAccelerationValid = false;
    int replayAge = 0;

    double* nodeBuffer = nullptr;
    int nodeCount = 0;
    int iterations = 1;
    uint32_t randomState = 1;

    double jiggle() {
        randomState = 1664525u * randomState + 1013904223u;
        return (static_cast<double>(randomState) / 4294967296.0 - 0.5) * 1e-6;
    }

public:
    void force(double alphaValue) {
        if (!nodeBuffer || nodeCount <= 0) return;
        const double alpha = alphaValue;
        for (int index = 0; index < nodeCount; ++index) {
            const int offset = index * kStride;
            initialVelocity[index] = {nodeBuffer[offset + 2], nodeBuffer[offset + 3]};
        }
        const int linkCount = static_cast<int>(sources.size());
        for (int iteration = 0; iteration < iterations; ++iteration) {
            for (int index = 0; index < linkCount; ++index) {
                const int source = sources[index];
                const int target = targets[index];
                if (source < 0 || source >= nodeCount || target < 0 || target >= nodeCount) continue;
                const int sourceOffset = source * kStride;
                const int targetOffset = target * kStride;
                double dx = nodeBuffer[targetOffset] + nodeBuffer[targetOffset + 2] -
                           nodeBuffer[sourceOffset] - nodeBuffer[sourceOffset + 2];
                double dy = nodeBuffer[targetOffset + 1] + nodeBuffer[targetOffset + 3] -
                           nodeBuffer[sourceOffset + 1] - nodeBuffer[sourceOffset + 3];
                if (dx == 0) dx = jiggle();
                if (dy == 0) dy = jiggle();
                const double length = std::sqrt(dx * dx + dy * dy);
                if (length == 0) continue;
                const double scale = (length - distances[index]) / length * alpha * strengths[index];
                dx *= scale;
                dy *= scale;
                const double targetWeight = biases[index];
                nodeBuffer[targetOffset + 2] -= dx * targetWeight;
                nodeBuffer[targetOffset + 3] -= dy * targetWeight;
                const double sourceWeight = 1.0 - targetWeight;
                nodeBuffer[sourceOffset + 2] += dx * sourceWeight;
                nodeBuffer[sourceOffset + 3] += dy * sourceWeight;
            }
        }

        const bool hadCachedAcceleration = cachedAccelerationValid;
        const float trendDivisor = static_cast<float>(replayAge + 1);
        if (alpha != 0) {
            for (int index = 0; index < nodeCount; ++index) {
                const int offset = index * kStride;
                const Acceleration acceleration = {
                    static_cast<float>((nodeBuffer[offset + 2] - initialVelocity[index].x) / alpha),
                    static_cast<float>((nodeBuffer[offset + 3] - initialVelocity[index].y) / alpha)
                };
                accelerationTrend[index] = hadCachedAcceleration
                    ? Acceleration{
                        (acceleration.x - cachedAcceleration[index].x) / trendDivisor,
                        (acceleration.y - cachedAcceleration[index].y) / trendDivisor
                    }
                    : Acceleration{0, 0};
                cachedAcceleration[index] = acceleration;
            }
        }
        cachedAccelerationValid = alpha != 0;
        replayAge = 0;
    }

    void replay(double alpha) {
        if (!nodeBuffer || !cachedAccelerationValid) return;
        const int age = advanceReplayAge();
        const int predictionAge = std::min(age, 4);
        const float replayAlpha = static_cast<float>(alpha);
        for (int index = 0; index < nodeCount; ++index) {
            const int offset = index * kStride;
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

    void setLinks(
        const emscripten::val& sourceValues,
        const emscripten::val& targetValues,
        const emscripten::val& biasValues,
        const emscripten::val& strengthValues,
        const emscripten::val& distanceValues
    ) {
        const int count = sourceValues["length"].as<int>();
        sources.resize(count);
        targets.resize(count);
        biases.resize(count);
        strengths.resize(count);
        distances.resize(count);
        for (int index = 0; index < count; ++index) {
            sources[index] = sourceValues[index].as<int>();
            targets[index] = targetValues[index].as<int>();
            biases[index] = biasValues[index].as<double>();
            strengths[index] = strengthValues[index].as<double>();
            distances[index] = distanceValues[index].as<double>();
        }
    }

    void setNodeBuffer(uintptr_t pointer, int count) {
        nodeBuffer = pointer && count > 0 ? reinterpret_cast<double*>(pointer) : nullptr;
        nodeCount = count;
        initialVelocity.resize(count);
        cachedAcceleration.resize(count);
        accelerationTrend.resize(count);
        cachedAccelerationValid = false;
        replayAge = 0;
    }

    void setIterations(int value) { iterations = std::max(1, value); }
    int getIterations() const { return iterations; }
};

ForceLink* createForceLink() { return new ForceLink(); }

EMSCRIPTEN_BINDINGS(force_link_module) {
    emscripten::class_<ForceLink>("ForceLink")
        .constructor<>()
        .function("force", &ForceLink::force)
        .function("replay", &ForceLink::replay)
        .function("advanceReplayAge", &ForceLink::advanceReplayAge)
        .function("getCachedAccelerationPointer", &ForceLink::getCachedAccelerationPointer)
        .function("getAccelerationTrendPointer", &ForceLink::getAccelerationTrendPointer)
        .function("setLinks", &ForceLink::setLinks)
        .function("setNodeBuffer", &ForceLink::setNodeBuffer)
        .function("setIterations", &ForceLink::setIterations)
        .function("getIterations", &ForceLink::getIterations);

    emscripten::function("createForceLink", &createForceLink, emscripten::allow_raw_pointers());
}
