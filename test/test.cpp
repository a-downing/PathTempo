#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <string_view>
#include <vector>

#include "path_tempo/Sampling.h"
#include "path_tempo/LinearOptimization.h"
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

    void testPersistentLinearOptimization() {
        path_tempo::PersistentLinearSolver solver;
        const auto configured = solver.configure(0.5);

        require(configured.has_value(), "persistent linear solver should configure");

        path_tempo::SparseLinearProgram first(2);
        first.columnLower(0) = 0.0;
        first.columnLower(1) = 0.0;
        first.columnCost(0) = 1.0;
        first.columnCost(1) = 2.0;
        first.addRow(3.0, path_tempo::linearProgramInfinity(), {{0, 1.0}, {1, 1.0}});
        const auto initial = solver.solve(first, {.simplexIterationLimit = 128});

        require(initial.has_value(), "initial linear program should solve");
        require(initial->status == path_tempo::LinearSolveStatus::Optimal, "initial linear program should be optimal");
        require(initial->values.size() == 2, "linear solution should contain every column");
        require(std::abs(initial->values[0] - 3.0) <= 1e-9, "linear solution should select the cheaper first column");
        require(std::abs(initial->values[1]) <= 1e-9, "linear solution should leave the expensive second column at zero");

        path_tempo::SparseLinearProgram updated(2);
        updated.columnLower(0) = 0.0;
        updated.columnLower(1) = 0.0;
        updated.columnCost(0) = 2.0;
        updated.columnCost(1) = 1.0;
        updated.addRow(3.0, path_tempo::linearProgramInfinity(), {{0, 1.0}, {1, 1.0}});
        const auto resolved = solver.solve(updated, {.simplexIterationLimit = 128});

        require(resolved.has_value(), "updated linear program should solve");
        require(resolved->status == path_tempo::LinearSolveStatus::Optimal, "updated linear program should be optimal");
        require(resolved->diagnostics.modelUpdateAttempted, "second structure-stable solve should attempt a model update");
        require(resolved->diagnostics.modelUpdateApplied, "second structure-stable solve should update the persistent model");
        require(std::abs(resolved->values[0]) <= 1e-9, "updated solution should leave the expensive first column at zero");
        require(std::abs(resolved->values[1] - 3.0) <= 1e-9, "updated solution should select the cheaper second column");
    }

    void testVelocityReachability() {
        const auto transitionDistance = path_tempo::velocityTransitionDistance(0.0, 2.0, 1.0, 2.0);

        require(std::abs(transitionDistance - 2.5) <= 1e-12, "velocity transition distance should include jerk-limited ramps");
        require(path_tempo::reachableVelocity(1.0, 0.5, 10.0, 1.0, 2.0) == 0.5, "a lower velocity cap should be immediately reachable");

        const auto reachable = path_tempo::reachableVelocity(0.0, 2.0, 1.0, 1.0, 2.0);

        require(reachable > 0.0 && reachable < 2.0, "a short piece should reduce its reachable velocity cap");
        require(path_tempo::velocityTransitionDistance(0.0, reachable, 1.0, 2.0) <= 1.0, "reachable velocity should fit within the available distance");
    }

    path_tempo::Limits<3> testPathLimits() {
        return {
            .pathAcceleration = 2.0,
            .pathJerk = 5.0,
            .axisVelocity = {4.0, 4.0, 4.0},
            .axisAcceleration = {2.0, 2.0, 2.0},
            .axisJerk = {5.0, 5.0, 5.0},
        };
    }

    void testMultiPiecePathPlanning() {
        const auto first = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {0.0, 0.0, 0.0},
            .to = {10.0, 0.0, 0.0},
        }, 31, 3.0);
        const auto second = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {10.0, 0.0, 0.0},
            .to = {20.0, 0.0, 0.0},
        }, 32, 1.5);

        require(first.has_value() && second.has_value(), "multi-piece test lines should sample");

        const std::array pieces {first->view(), second->view()};
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        });

        require(planned.has_value(), "a continuous two-piece line path should solve");
        require(!planned->timeLaw.segments.empty(), "multi-piece planning should produce cubic segments");
        require(planned->pieceBoundaries.size() == 3, "two pieces should produce three boundary states");
        require(planned->pieceBoundaries[1].velocity > 0.0, "continuous planning should retain motion at an internal boundary");
        require(planned->pieceBoundaries[1].velocity <= 1.5 + 1e-10, "an internal boundary should honor both programmed velocities");
        require(planned->diagnostics.velocitySeedDuration > 0.0, "multi-piece planning should report its materialized duration");

        auto previousDistance = 0.0;
        auto previousVelocity = 0.0;
        auto previousAcceleration = 0.0;
        auto sawFirstPiece = false;
        auto sawSecondPiece = false;

        for (const auto &segment : planned->timeLaw.segments) {
            require(std::abs(segment.position(0.0) - previousDistance) <= 1e-8, "multi-piece cubics should be globally position-continuous");
            require(std::abs(segment.velocity(0.0) - previousVelocity) <= 1e-8, "multi-piece cubics should be velocity-continuous");
            require(std::abs(segment.acceleration(0.0) - previousAcceleration) <= 1e-8, "multi-piece cubics should be acceleration-continuous");
            require(std::abs(segment.jerk()) <= 5.0 + 1e-10, "multi-piece cubics should obey the scalar and axis jerk limit");

            const auto maximumProgrammedVelocity = segment.piece == 31 ? 3.0 : 1.5;
            auto maximumSegmentVelocity = std::max(segment.velocity(0.0), segment.velocity(segment.duration));

            if (std::abs(segment.jerk()) > 1e-15) {
                const auto velocityExtremum = -segment.acceleration(0.0) / segment.jerk();

                if (velocityExtremum > 0.0 && velocityExtremum < segment.duration) {
                    maximumSegmentVelocity = std::max(maximumSegmentVelocity, segment.velocity(velocityExtremum));
                }
            }

            require(maximumSegmentVelocity <= maximumProgrammedVelocity + 1e-9, "multi-piece cubics should obey each programmed velocity");
            require(std::abs(segment.acceleration(0.0)) <= 2.0 + 1e-9 && std::abs(segment.acceleration(segment.duration)) <= 2.0 + 1e-9, "multi-piece cubics should obey the scalar and axis acceleration limit");

            previousDistance = segment.position(segment.duration);
            previousVelocity = segment.velocity(segment.duration);
            previousAcceleration = segment.acceleration(segment.duration);
            sawFirstPiece = sawFirstPiece || segment.piece == 31;
            sawSecondPiece = sawSecondPiece || segment.piece == 32;
        }

        require(sawFirstPiece && sawSecondPiece, "multi-piece segments should preserve each piece identity");
        require(std::abs(previousDistance - 20.0) <= 1e-8, "multi-piece planning should cover the complete path distance");
        require(std::abs(previousVelocity) <= 1e-8, "rest-to-rest multi-piece planning should end at rest");
        require(std::abs(previousAcceleration) <= 1e-8, "rest-to-rest multi-piece planning should end at zero acceleration");

        const auto replanned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        });

        require(replanned.has_value(), "a repeated path solve should succeed");
        require(replanned->diagnostics.linearSolverBasisReused, "a repeated path solve should reuse the HiGHS basis");
    }

    void testMultiPiecePathRejectsTangentDiscontinuity() {
        const auto first = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {0.0, 0.0, 0.0},
            .to = {2.0, 0.0, 0.0},
        }, 41, 1.0);
        const auto second = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {2.0, 0.0, 0.0},
            .to = {2.0, 2.0, 0.0},
        }, 42, 1.0);

        require(first.has_value() && second.has_value(), "discontinuous test lines should sample");

        const std::array pieces {first->view(), second->view()};
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        });

        require(!planned.has_value(), "multi-piece planning should reject a tangent discontinuity");
        require(planned.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a tangent discontinuity should be an input error");
    }

    path_tempo::SampledPathPiece<3> sampledUnitCircleInterval(const path_tempo::PathPieceId id, const double fromAngle, const double toAngle, const double programmedVelocity, const std::size_t intervals = 32) {
        const auto length = toAngle - fromAngle;
        path_tempo::SampledPathPiece<3> piece {
            .id = id,
            .length = length,
            .programmedVelocity = programmedVelocity,
            .initialLimits = {},
            .stations = {},
        };
        piece.stations.reserve(intervals + 1);

        for (std::size_t sample = 0; sample <= intervals; ++sample) {
            const auto distance = length * static_cast<double>(sample) / static_cast<double>(intervals);
            const auto angle = fromAngle + distance;
            piece.stations.push_back({
                .distance = distance,
                .tangent = {-std::sin(angle), std::cos(angle), 0.0},
                .curvature = {-std::cos(angle), -std::sin(angle), 0.0},
                .thirdDerivative = {std::sin(angle), -std::cos(angle), 0.0},
            });
        }

        return piece;
    }

    path_tempo::SampledPathPiece<3> sampledQuarterCircle(const path_tempo::PathPieceId id, const double programmedVelocity, const std::size_t intervals = 32) {
        return sampledUnitCircleInterval(id, 0.0, std::numbers::pi / 2.0, programmedVelocity, intervals);
    }

    void testCurvedPathPlanning() {
        const auto curve = sampledQuarterCircle(51, 4.0);
        const std::array pieces {curve.view()};
        path_tempo::PathPlanner planner;
        const path_tempo::Limits<3> limits {
            .pathAcceleration = 0.8,
            .pathJerk = 1.2,
            .axisVelocity = {4.0, 4.0, 4.0},
            .axisAcceleration = {0.8, 0.8, 0.8},
            .axisJerk = {1.2, 1.2, 1.2},
        };
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = {
                .linearSolveTimeLimit = 0.25,
                .simplexIterationLimit = 4096,
                .maximumCorrectionPasses = 12,
            },
        });

        require(planned.has_value(), "a sampled curved path should solve");
        require(!planned->timeLaw.segments.empty(), "curved planning should produce a scalar time law");
        require(planned->diagnostics.correctionPasses >= 1, "curved planning should report its bounded correction passes");
        require(planned->diagnostics.correctedPieces == 1, "coupled curved constraints should trigger a local correction for this path");

        for (const auto &station : curve.stations) {
            const auto tolerance = 1e-10;

            for (const auto &segment : planned->timeLaw.segments) {
                const auto localFrom = segment.position(0.0);
                const auto localTo = segment.position(segment.duration);

                if (station.distance < localFrom - tolerance || station.distance > localTo + tolerance) {
                    continue;
                }

                auto low = 0.0;
                auto high = segment.duration;

                for (unsigned iteration = 0; iteration < 56; ++iteration) {
                    const auto middle = std::midpoint(low, high);

                    if (segment.position(middle) < station.distance) {
                        low = middle;
                    } else {
                        high = middle;
                    }
                }

                const auto time = std::midpoint(low, high);
                const auto velocity = segment.velocity(time);
                const auto acceleration = segment.acceleration(time);
                const auto jerk = segment.jerk();
                auto coupledAccelerationSquared = 0.0;
                auto coupledJerkSquared = 0.0;

                for (std::size_t axis = 0; axis < 3; ++axis) {
                    const auto coupledAcceleration = station.tangent[axis] * acceleration + station.curvature[axis] * velocity * velocity;
                    const auto coupledJerk = station.tangent[axis] * jerk + 3.0 * station.curvature[axis] * velocity * acceleration + station.thirdDerivative[axis] * velocity * velocity * velocity;
                    coupledAccelerationSquared += coupledAcceleration * coupledAcceleration;
                    coupledJerkSquared += coupledJerk * coupledJerk;
                    require(std::abs(station.tangent[axis] * velocity) <= limits.axisVelocity[axis] * (1.0 + 1e-8), "curved timing should obey sampled axis velocity");
                    require(std::abs(coupledAcceleration) <= limits.axisAcceleration[axis] * (1.0 + 1e-8), "curved timing should obey sampled coupled axis acceleration");
                    require(std::abs(coupledJerk) <= limits.axisJerk[axis] * (1.0 + 1e-8), "curved timing should obey sampled coupled axis jerk");
                }

                require(std::sqrt(coupledAccelerationSquared) <= limits.pathAcceleration * (1.0 + 1e-8), "curved timing should obey sampled coupled path acceleration");
                require(std::sqrt(coupledJerkSquared) <= limits.pathJerk * (1.0 + 1e-8), "curved timing should obey sampled coupled path jerk");
            }
        }
    }

    void testCurvedBoundaryRejectsGeometricJerkViolation() {
        const auto curve = sampledQuarterCircle(61, 2.0);
        const std::array pieces {curve.view()};
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {.velocity = 0.5},
            .ending = {.velocity = 0.5},
            .limits = {
                .pathAcceleration = 1.0,
                .pathJerk = 0.1,
                .axisVelocity = {2.0, 2.0, 2.0},
                .axisAcceleration = {1.0, 1.0, 1.0},
                .axisJerk = {0.1, 0.1, 0.1},
            },
            .settings = {},
        });

        require(!planned.has_value(), "a fixed curved boundary velocity above the geometric jerk cap should be rejected");
        require(planned.error().code == path_tempo::PlanningErrorCode::InvalidInput, "an infeasible curved boundary should be an input error");
    }

    void testMultiPieceCurvedRefinement() {
        const auto first = sampledUnitCircleInterval(71, 0.0, std::numbers::pi / 4.0, 2.0, 16);
        const auto second = sampledUnitCircleInterval(72, std::numbers::pi / 4.0, std::numbers::pi / 2.0, 2.0, 16);
        const std::array pieces {first.view(), second.view()};
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = {
                .pathAcceleration = 1.0,
                .pathJerk = 2.0,
                .axisVelocity = {3.0, 3.0, 3.0},
                .axisAcceleration = {1.0, 1.0, 1.0},
                .axisJerk = {2.0, 2.0, 2.0},
            },
            .settings = {},
        });

        require(planned.has_value(), "a C2 two-piece curved chain should solve");
        require(planned->pieceBoundaries[1].velocity > 0.0, "a curved internal boundary should retain continuous motion");
        require(planned->diagnostics.sequentialSolves >= 1, "a multi-piece curved path should run coupled sequential refinement");

        auto previousDistance = 0.0;
        auto previousVelocity = 0.0;
        auto previousAcceleration = 0.0;

        for (const auto &segment : planned->timeLaw.segments) {
            require(std::abs(segment.position(0.0) - previousDistance) <= 1e-8, "refined curved cubics should be position-continuous");
            require(std::abs(segment.velocity(0.0) - previousVelocity) <= 1e-8, "refined curved cubics should be velocity-continuous");
            require(std::abs(segment.acceleration(0.0) - previousAcceleration) <= 1e-8, "refined curved cubics should be acceleration-continuous");
            previousDistance = segment.position(segment.duration);
            previousVelocity = segment.velocity(segment.duration);
            previousAcceleration = segment.acceleration(segment.duration);
        }
    }

    void testCurvedMovingBoundaryStates() {
        const auto curve = sampledQuarterCircle(81, 1.5);
        const std::array pieces {curve.view()};
        const path_tempo::BoundaryState beginning {.velocity = 0.2, .acceleration = 0.05};
        const path_tempo::BoundaryState ending {.velocity = 0.25, .acceleration = -0.04};
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = beginning,
            .ending = ending,
            .limits = {
                .pathAcceleration = 2.0,
                .pathJerk = 4.0,
                .axisVelocity = {2.0, 2.0, 2.0},
                .axisAcceleration = {2.0, 2.0, 2.0},
                .axisJerk = {4.0, 4.0, 4.0},
            },
            .settings = {},
        });

        require(planned.has_value(), "a curved path with nonzero boundary PVA should solve");
        require(std::abs(planned->timeLaw.segments.front().velocity(0.0) - beginning.velocity) <= 1e-10, "curved planning should preserve beginning velocity");
        require(std::abs(planned->timeLaw.segments.front().acceleration(0.0) - beginning.acceleration) <= 1e-10, "curved planning should preserve beginning acceleration");

        const auto &last = planned->timeLaw.segments.back();
        require(std::abs(last.velocity(last.duration) - ending.velocity) <= 1e-9, "curved planning should preserve ending velocity");
        require(std::abs(last.acceleration(last.duration) - ending.acceleration) <= 1e-9, "curved planning should preserve ending acceleration");
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
    testPersistentLinearOptimization();
    testVelocityReachability();
    testMultiPiecePathPlanning();
    testMultiPiecePathRejectsTangentDiscontinuity();
    testCurvedPathPlanning();
    testCurvedBoundaryRejectsGeometricJerkViolation();
    testMultiPieceCurvedRefinement();
    testCurvedMovingBoundaryStates();
}
