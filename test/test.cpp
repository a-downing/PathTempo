#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string_view>

#include "path_tempo/Sampling.h"
#include "path_tempo/Types.h"
#include "path_tempo/Version.h"

namespace {
    void require(const bool condition, const std::string_view message) {
        if (!condition) {
            std::cerr << message << '\n';
            std::exit(1);
        }
    }

    void testLineSampling() {
        const path_tempo::Line<3> line {
            .from = {0.0, 0.0, 0.0},
            .to = {3.0, 4.0, 0.0},
        };
        const auto sampled = path_tempo::sampleLine(line, 7, 2.0, {.intervals = 4});

        require(sampled.has_value(), "line sampling should succeed");
        require(sampled->id == 7, "line sampling should preserve piece identity");
        require(std::abs(sampled->length - 5.0) <= 1e-12, "line sampling should calculate length");
        require(sampled->stations.size() == 5, "four intervals should produce five samples");
        require(sampled->stations.front().distance == 0.0, "first sample should begin at zero");
        require(sampled->stations.back().distance == sampled->length, "last sample should end at piece length");
        require(std::abs(sampled->stations.front().tangent[0] - 0.6) <= 1e-12, "line tangent X should be normalized");
        require(std::abs(sampled->stations.front().tangent[1] - 0.8) <= 1e-12, "line tangent Y should be normalized");
    }

    void testLineSamplingErrors() {
        const path_tempo::Line<3> line {
            .from = {0.0, 0.0, 0.0},
            .to = {3.0, 4.0, 0.0},
        };
        auto invalidVelocity = path_tempo::sampleLine(line, 1, 0.0);

        require(!invalidVelocity.has_value(), "line sampling should reject zero programmed velocity");
        require(invalidVelocity.error() == path_tempo::SamplingError::InvalidProgrammedVelocity, "zero programmed velocity should have a specific error");

        auto nonFiniteLine = line;
        nonFiniteLine.to[0] = std::numeric_limits<double>::infinity();
        const auto nonFinite = path_tempo::sampleLine(nonFiniteLine, 1, 1.0);

        require(!nonFinite.has_value(), "line sampling should reject non-finite geometry");
        require(nonFinite.error() == path_tempo::SamplingError::NonFiniteGeometry, "non-finite geometry should have a specific error");
    }

    void testCubicTimeSegment() {
        const path_tempo::CubicTimeSegment segment {
            .piece = 4,
            .duration = 2.0,
            .c0 = 1.0,
            .c1 = 2.0,
            .c2 = 3.0,
            .c3 = 4.0,
        };

        require(segment.position(2.0) == 49.0, "cubic position evaluation should use physical time");
        require(segment.velocity(2.0) == 62.0, "cubic velocity evaluation should be exact");
        require(segment.acceleration(2.0) == 54.0, "cubic acceleration evaluation should be exact");
        require(segment.jerk() == 24.0, "cubic jerk should be constant");
    }
}

int main() {
    require(path_tempo::version() == "0.1.0", "version should match the project version");
    testLineSampling();
    testLineSamplingErrors();
    testCubicTimeSegment();
}
