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
    }

    void setIterations(int value) { iterations = std::max(1, value); }
    int getIterations() const { return iterations; }
};

ForceLink* createForceLink() { return new ForceLink(); }

EMSCRIPTEN_BINDINGS(force_link_module) {
    emscripten::class_<ForceLink>("ForceLink")
        .constructor<>()
        .function("force", &ForceLink::force)
        .function("setLinks", &ForceLink::setLinks)
        .function("setNodeBuffer", &ForceLink::setNodeBuffer)
        .function("setIterations", &ForceLink::setIterations)
        .function("getIterations", &ForceLink::getIterations);

    emscripten::function("createForceLink", &createForceLink, emscripten::allow_raw_pointers());
}
