#include "../examples/endpoint_exact_arc.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string_view>

namespace {
    constexpr double TOLERANCE = 1e-12;

    void require(const bool condition, const std::string_view message) {
        if (!condition) {
            std::cerr << message << '\n';
            std::exit(1);
        }
    }

    void requirePosition(const path_tempo::example::Vector3 &actual, const path_tempo::example::Vector3 &expected, const std::string_view message) {
        for (std::size_t component = 0; component < actual.size(); ++component) {
            require(std::abs(actual[component] - expected[component]) <= TOLERANCE, message);
        }
    }
}

int main() {
    const path_tempo::example::ArcGeometry source {
        .from = {1.0, 0.0, 3.0},
        .to = {0.0, 2.0, 5.0},
        .center = {0.0, 0.0, 1.0},
        .axis = {0.0, 0.0, 2.0},
    };
    const auto arc = path_tempo::example::EndpointExactArc::create(source);

    require(arc.has_value(), "an axially offset helical arc should be accepted");
    requirePosition(arc->positionAtParameter(0.0), source.from, "the arc should reproduce its exact start point");
    requirePosition(arc->positionAtParameter(1.0), source.to, "the arc should reproduce its exact end point");
    require(arc->length() > 0.0, "the arc should have positive length");

    return 0;
}
