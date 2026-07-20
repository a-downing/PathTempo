#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string_view>

#include "path_tempo/Sampling.h"
#include "path_tempo/Planner.h"
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

    void testRestToRestTransition() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 12,
            .length = 10.0,
            .beginning = {},
            .ending = {},
            .maximumVelocity = 4.0,
            .maximumAcceleration = 3.0,
            .maximumJerk = 8.0,
        });

        require(transition.has_value(), "rest-to-rest transition should solve");
        require(!transition->empty(), "rest-to-rest transition should contain cubic segments");
        require(transition->duration() > 0.0, "rest-to-rest transition should have positive duration");
        require(transition->front().piece == 12, "transition segments should preserve piece identity");
        require(std::abs(transition->front().position(0.0)) <= 1e-12, "transition should start at zero distance");
        require(std::abs(transition->front().velocity(0.0)) <= 1e-12, "transition should start at rest");

        const auto &last = transition->back();
        require(std::abs(last.position(last.duration) - 10.0) <= 1e-9, "transition should reach its requested distance");
        require(std::abs(last.velocity(last.duration)) <= 1e-9, "transition should end at rest");
        require(std::abs(last.acceleration(last.duration)) <= 1e-9, "transition should end at zero acceleration");

        auto previousDistance = 0.0;

        for (const auto &segment : *transition) {
            require(segment.duration > 0.0, "every transition segment should have positive duration");
            require(std::abs(segment.position(0.0) - previousDistance) <= 1e-9, "transition segments should be position-continuous");
            previousDistance = segment.position(segment.duration);
        }
    }

    void testMovingBoundaryTransition() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 18,
            .length = 2.5,
            .beginning = {.velocity = 0.5, .acceleration = 0.25},
            .ending = {.velocity = 0.75, .acceleration = -0.1},
            .maximumVelocity = 2.0,
            .maximumAcceleration = 1.0,
            .maximumJerk = 4.0,
        });

        require(transition.has_value(), "moving-boundary transition should solve");
        require(std::abs(transition->front().velocity(0.0) - 0.5) <= 1e-12, "transition should preserve initial velocity");
        require(std::abs(transition->front().acceleration(0.0) - 0.25) <= 1e-12, "transition should preserve initial acceleration");

        const auto &last = transition->back();
        require(std::abs(last.position(last.duration) - 2.5) <= 1e-9, "moving-boundary transition should reach its distance");
        require(std::abs(last.velocity(last.duration) - 0.75) <= 1e-9, "moving-boundary transition should reach its final velocity");
        require(std::abs(last.acceleration(last.duration) + 0.1) <= 1e-9, "moving-boundary transition should reach its final acceleration");
    }

    void testUnboundedAccelerationCruise() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 20,
            .length = 3.0,
            .beginning = {.velocity = 1.5},
            .ending = {.velocity = 1.5},
            .maximumVelocity = 2.0,
            .maximumAcceleration = std::numeric_limits<double>::infinity(),
            .maximumJerk = 4.0,
        });

        require(transition.has_value(), "unbounded-acceleration cruise should solve");
        require(transition->size() == 1, "unbounded-acceleration cruise should contain one segment");
        require(transition->duration() == 2.0, "unbounded-acceleration cruise should use distance divided by velocity");
    }

    void testInvalidTransition() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 21,
            .length = 0.0,
            .beginning = {},
            .ending = {},
            .maximumVelocity = 2.0,
            .maximumAcceleration = 1.0,
            .maximumJerk = 4.0,
        });

        require(!transition.has_value(), "zero-length transition should be rejected");
        require(transition.error().code == path_tempo::PlanningErrorCode::InvalidInput, "invalid transition should return a typed input error");
    }
}

int main() {
    require(path_tempo::version() == "0.1.0", "version should match the project version");
    testLineSampling();
    testLineSamplingErrors();
    testCubicTimeSegment();
    testRestToRestTransition();
    testMovingBoundaryTransition();
    testUnboundedAccelerationCruise();
    testInvalidTransition();
}
