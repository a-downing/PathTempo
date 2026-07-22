#include "../examples/endpoint_exact_arc.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
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

    constexpr double tinyRadius = 1e-10;
    const path_tempo::example::ArcGeometry tinyQuarterCircle {
        .from = {tinyRadius, 0.0, 0.0},
        .to = {0.0, tinyRadius, 0.0},
        .center = {0.0, 0.0, 0.0},
        .axis = {0.0, 0.0, 1.0},
    };
    const auto tinyQuarterArc = path_tempo::example::EndpointExactArc::create(tinyQuarterCircle);

    require(tinyQuarterArc.has_value(), "a tiny quarter-circle should be accepted");
    require(std::abs(tinyQuarterArc->length() / tinyRadius - std::numbers::pi / 2.0) <= TOLERANCE,
        "a tiny quarter-circle should not be classified as a complete revolution");

    auto tinyCompleteCircle = tinyQuarterCircle;
    tinyCompleteCircle.to = tinyCompleteCircle.from;
    const auto tinyCompleteArc = path_tempo::example::EndpointExactArc::create(tinyCompleteCircle);

    require(tinyCompleteArc.has_value(), "a tiny complete circle should be accepted");
    require(std::abs(tinyCompleteArc->length() / tinyRadius - 2.0 * std::numbers::pi) <= TOLERANCE,
        "scale-relative closure detection should preserve a tiny complete revolution");

    return 0;
}
