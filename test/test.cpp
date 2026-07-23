#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <string_view>
#include <vector>

#include "path_tempo/Sampling.h"
#include "path_tempo/Planner.h"
#include "path_tempo/Types.h"
#include "path_tempo/Version.h"

#include "../src/BoundaryRefinement.h"

namespace {
    // Closed-form line and transition inputs should agree to near machine
    // precision; this allowance covers only ordinary arithmetic roundoff.
    constexpr double ANALYTIC_RESULT_TOLERANCE = 1e-12;
    // Ruckig reconstructs boundary states through several polynomial phases,
    // so endpoint assertions allow modest accumulated solver roundoff.
    constexpr double SOLVER_ENDPOINT_TOLERANCE = 1e-9;
    // Multi-piece and refined paths accumulate error across many phases; this
    // remains tight enough to expose a visible continuity defect.
    constexpr double TIME_LAW_CONTINUITY_TOLERANCE = 1e-8;
    // Coupled limits are sampled after several transformations, so allow a
    // small relative margin beyond the stricter planner-side feasibility check.
    constexpr double SAMPLED_LIMIT_RELATIVE_TOLERANCE = 1e-8;
    // Boundary-limit assertions need only absorb the planner's relative
    // constraint-comparison roundoff.
    constexpr double BOUNDARY_LIMIT_TOLERANCE = 1e-10;
    // Treat jerk below this scale as zero before dividing to locate an extremum.
    constexpr double JERK_ZERO_TOLERANCE = 1e-15;
    // Station association uses the scalar transition's accepted endpoint error.
    constexpr double STATION_MATCH_ABSOLUTE_TOLERANCE = 1e-12;
    constexpr double STATION_MATCH_RELATIVE_TOLERANCE = 1e-9;

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
        const auto sampled = path_tempo::sampleLine(line, 7, 2.0);

        require(sampled.has_value(), "line sampling should succeed");
        require(sampled->id == 7, "line sampling should preserve piece identity");
        require(std::abs(sampled->length - 5.0) <= ANALYTIC_RESULT_TOLERANCE, "line sampling should calculate length");
        require(sampled->stations.size() == 2, "a line should require only its endpoint samples");
        require(sampled->stations.front().distance == 0.0, "first sample should begin at zero");
        require(sampled->stations.back().distance == sampled->length, "last sample should end at piece length");
        require(std::abs(sampled->stations.front().tangent[0] - 0.6) <= ANALYTIC_RESULT_TOLERANCE, "line tangent X should be normalized");
        require(std::abs(sampled->stations.front().tangent[1] - 0.8) <= ANALYTIC_RESULT_TOLERANCE, "line tangent Y should be normalized");
    }

    void testTinyLineSampling() {
        constexpr double scale = 1e-200;
        const path_tempo::Line<2> line {
            .from = {0.0, 0.0},
            .to = {3.0 * scale, 4.0 * scale},
        };
        const auto sampled = path_tempo::sampleLine(line, 8, 2.0);

        require(sampled.has_value(), "a representable tiny line should sample without norm underflow");
        require(std::abs(sampled->length / scale - 5.0) <= ANALYTIC_RESULT_TOLERANCE, "tiny line length should retain its geometric scale");
        require(std::abs(sampled->stations.front().tangent[0] - 0.6) <= ANALYTIC_RESULT_TOLERANCE, "tiny line tangent X should be normalized");
        require(std::abs(sampled->stations.front().tangent[1] - 0.8) <= ANALYTIC_RESULT_TOLERANCE, "tiny line tangent Y should be normalized");
    }

    void testLineSamplingErrors() {
        const path_tempo::Line<3> line {
            .from = {0.0, 0.0, 0.0},
            .to = {3.0, 4.0, 0.0},
        };
        auto invalidVelocity = path_tempo::sampleLine(line, 1, 0.0);

        require(!invalidVelocity.has_value(), "line sampling should reject a zero maximum velocity");
        require(invalidVelocity.error() == path_tempo::SamplingError::InvalidMaxVelocity, "a zero maximum velocity should have a specific error");

        auto nonFiniteLine = line;
        nonFiniteLine.to[0] = std::numeric_limits<double>::infinity();
        const auto nonFinite = path_tempo::sampleLine(nonFiniteLine, 1, 1.0);

        require(!nonFinite.has_value(), "line sampling should reject non-finite geometry");
        require(nonFinite.error() == path_tempo::SamplingError::NonFiniteGeometry, "non-finite geometry should have a specific error");
    }

    void testArcSampling() {
        const path_tempo::Arc<3> arc {
            .origin = {1.0, 0.0, 0.0},
            .radial = {1.0, 0.0, 0.0},
            .tangential = {0.0, 1.0, 0.0},
            .axial = {0.0, 0.0, 2.0},
            .sweep = std::numbers::pi / 2.0,
        };
        const auto sampled = path_tempo::sampleArc(arc, 10, 3.0, 16);

        require(sampled.has_value(), "helical arc sampling should succeed");
        require(sampled->stations.size() == 17, "arc sampling should honor its interval count");
        require(std::abs(sampled->length - std::hypot(std::numbers::pi / 2.0, 2.0)) <= SOLVER_ENDPOINT_TOLERANCE, "helical arc sampling should integrate its spatial length");
        require(sampled->stations.front().distance == 0.0 && sampled->stations.back().distance == sampled->length, "arc stations should cover the complete curve");
    }

    void testBSplineSampling() {
        const std::array<path_tempo::Vector<3>, 7> controls {{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {2.0, 0.0, 0.0},
            {3.0, 0.0, 0.0},
            {4.0, 0.0, 0.0},
            {5.0, 0.0, 0.0},
            {6.0, 0.0, 0.0},
        }};
        const std::array knots {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0};
        const auto sampled = path_tempo::sampleBSpline(path_tempo::BSpline<3> {
            .degree = 3,
            .controls = controls,
            .knots = knots,
        }, 20, 2.0, 8);

        require(sampled.has_value(), "cubic B-spline sampling should succeed");
        require(sampled->size() == 4, "a cubic B-spline should produce one piece per non-empty knot interval");
        require((*sampled)[0].id == 20 && (*sampled)[3].id == 23, "B-spline pieces should receive consecutive identities");
        require(std::ranges::all_of(*sampled, [](const auto &piece) { return piece.stations.size() == 9; }), "every B-spline interval should honor the requested sample count");
        const auto length = std::accumulate(sampled->begin(), sampled->end(), 0.0, [](const double total, const auto &piece) { return total + piece.length; });
        require(std::abs(length - 6.0) <= SOLVER_ENDPOINT_TOLERANCE, "collinear B-spline intervals should cover the complete curve length");
    }

    void testUnclampedNonUniformBSplineSampling() {
        const std::array<path_tempo::Vector<2>, 5> controls {{
            {0.0, 0.0},
            {1.0, 0.0},
            {2.0, 0.0},
            {3.0, 0.0},
            {4.0, 0.0},
        }};
        const std::array knots {0.0, 0.5, 1.0, 2.0, 4.0, 5.0, 6.0, 7.0};
        const auto sampled = path_tempo::sampleBSpline(path_tempo::BSpline<2> {
            .degree = 2,
            .controls = controls,
            .knots = knots,
        }, 24, 2.0, 4);

        require(sampled.has_value(), "unclamped non-uniform B-spline sampling should succeed");
        require(sampled->size() == 3, "an unclamped B-spline should still produce one piece per non-empty knot interval");
        require((*sampled)[0].id == 24 && (*sampled)[2].id == 26, "unclamped B-spline pieces should receive consecutive identities");
    }

    void testBSplineKnotScaleInvariance() {
        const std::array<path_tempo::Vector<2>, 2> controls {{{0.0, 0.0}, {1.0, 0.0}}};
        const std::array unitKnots {0.0, 0.0, 1.0, 1.0};
        const std::array scaledKnots {0.0, 0.0, 1e16, 1e16};
        const auto unit = path_tempo::sampleBSpline(path_tempo::BSpline<2> {
            .degree = 1,
            .controls = controls,
            .knots = unitKnots,
        }, 27, 2.0, 4);
        const auto scaled = path_tempo::sampleBSpline(path_tempo::BSpline<2> {
            .degree = 1,
            .controls = controls,
            .knots = scaledKnots,
        }, 28, 2.0, 4);

        require(unit.has_value() && scaled.has_value(), "uniform knot rescaling should preserve a valid B-spline");
        require(std::abs(unit->front().length - scaled->front().length) <= ANALYTIC_RESULT_TOLERANCE, "uniform knot rescaling should preserve B-spline length");
        require(std::abs(unit->front().stations.front().tangent[0] - scaled->front().stations.front().tangent[0]) <= ANALYTIC_RESULT_TOLERANCE, "uniform knot rescaling should preserve B-spline differential geometry");
    }

    void testNurbsSampling() {
        const std::array<path_tempo::Vector<2>, 3> controls {{
            {1.0, 0.0},
            {1.0, 1.0},
            {0.0, 1.0},
        }};
        const std::array weights {1.0, std::numbers::sqrt2 / 2.0, 1.0};
        const std::array knots {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
        const auto sampled = path_tempo::sampleNurbs(path_tempo::Nurbs<2> {
            .degree = 2,
            .controls = controls,
            .weights = weights,
            .knots = knots,
        }, 30, 2.0, 16);

        require(sampled.has_value(), "quadratic NURBS sampling should succeed");
        require(sampled->size() == 1 && sampled->front().stations.size() == 17, "a single NURBS knot interval should produce one sampled piece");
        require(std::abs(sampled->front().length - std::numbers::pi / 2.0) <= SOLVER_ENDPOINT_TOLERANCE, "rational quarter-circle sampling should recover exact arc length");
    }

    void testNurbsGeometricScaleInvariance() {
        constexpr double scale = 1e-9;
        const std::array<path_tempo::Vector<2>, 3> controls {{
            {scale, 0.0},
            {scale, scale},
            {0.0, scale},
        }};
        const std::array weights {1.0, std::numbers::sqrt2 / 2.0, 1.0};
        const std::array knots {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
        const auto sampled = path_tempo::sampleNurbs(path_tempo::Nurbs<2> {
            .degree = 2,
            .controls = controls,
            .weights = weights,
            .knots = knots,
        }, 32, 2.0, 16);

        require(sampled.has_value(), "a uniformly tiny NURBS curve should sample");
        require(std::abs(sampled->front().length / scale - std::numbers::pi / 2.0) <= SOLVER_ENDPOINT_TOLERANCE,
            "arc-length integration should retain relative accuracy below unit scale");
    }

    void testNurbsWeightScaleInvariance() {
        const std::array<path_tempo::Vector<2>, 3> controls {{{1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}}};
        const std::array weights {1e-100, std::numbers::sqrt2 / 2.0 * 1e-100, 1e-100};
        const std::array knots {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
        const auto sampled = path_tempo::sampleNurbs(path_tempo::Nurbs<2> {
            .degree = 2,
            .controls = controls,
            .weights = weights,
            .knots = knots,
        }, 31, 2.0, 16);

        require(sampled.has_value(), "uniformly tiny positive NURBS weights should preserve a valid curve");
        require(std::abs(sampled->front().length - std::numbers::pi / 2.0) <= SOLVER_ENDPOINT_TOLERANCE, "uniform NURBS weight scaling should preserve curve length");
    }

    void testParametricSamplingUsesSpeedOnlyEvaluator() {
        struct CountingEvaluator {
            mutable std::size_t stateCalls = 0;
            mutable std::size_t speedCalls = 0;

            std::expected<std::array<path_tempo::Vector<1>, 4>, path_tempo::SamplingError> operator()(const double parameter, const std::size_t) const {
                ++stateCalls;

                return std::array<path_tempo::Vector<1>, 4> {{
                    {parameter},
                    {1.0},
                    {},
                    {},
                }};
            }

            std::expected<double, path_tempo::SamplingError> parametricSpeed(const double, const std::size_t) const {
                ++speedCalls;

                return 1.0;
            }
        };

        constexpr std::size_t sampleIntervals = 4;
        const CountingEvaluator evaluator;
        const auto sampled = path_tempo::detail::sampleParametricInterval<1>(evaluator, 0.0, 1.0, 0, 33, 2.0, sampleIntervals);

        require(sampled.has_value(), "parametric sampling with a speed-only evaluator should succeed");
        require(evaluator.stateCalls == sampleIntervals + 1, "full derivative evaluation should be limited to differential stations");
        require(evaluator.speedCalls > 0, "arc-length integration should use the speed-only evaluator");
        require(std::abs(sampled->length - 1.0) <= ANALYTIC_RESULT_TOLERANCE, "speed-only integration should preserve curve length");
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
        require(std::abs(transition->front().position(0.0)) <= ANALYTIC_RESULT_TOLERANCE, "transition should start at zero distance");
        require(std::abs(transition->front().velocity(0.0)) <= ANALYTIC_RESULT_TOLERANCE, "transition should start at rest");

        const auto &last = transition->back();
        require(std::abs(last.position(last.duration) - 10.0) <= SOLVER_ENDPOINT_TOLERANCE, "transition should reach its requested distance");
        require(std::abs(last.velocity(last.duration)) <= SOLVER_ENDPOINT_TOLERANCE, "transition should end at rest");
        require(std::abs(last.acceleration(last.duration)) <= SOLVER_ENDPOINT_TOLERANCE, "transition should end at zero acceleration");

        auto previousDistance = 0.0;

        for (const auto &segment : *transition) {
            require(segment.duration > 0.0, "every transition segment should have positive duration");
            require(std::abs(segment.position(0.0) - previousDistance) <= SOLVER_ENDPOINT_TOLERANCE, "transition segments should be position-continuous");
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
        require(std::abs(transition->front().velocity(0.0) - 0.5) <= ANALYTIC_RESULT_TOLERANCE, "transition should preserve initial velocity");
        require(std::abs(transition->front().acceleration(0.0) - 0.25) <= ANALYTIC_RESULT_TOLERANCE, "transition should preserve initial acceleration");

        const auto &last = transition->back();
        require(std::abs(last.position(last.duration) - 2.5) <= SOLVER_ENDPOINT_TOLERANCE, "moving-boundary transition should reach its distance");
        require(std::abs(last.velocity(last.duration) - 0.75) <= SOLVER_ENDPOINT_TOLERANCE, "moving-boundary transition should reach its final velocity");
        require(std::abs(last.acceleration(last.duration) + 0.1) <= SOLVER_ENDPOINT_TOLERANCE, "moving-boundary transition should reach its final acceleration");
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

    void testUnboundedAccelerationCruiseValidation() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto unequalVelocity = planner.solve({
            .piece = 25,
            .length = 3.0,
            .beginning = {.velocity = 1.5},
            .ending = {.velocity = std::nextafter(1.5, 2.0)},
            .maximumVelocity = 2.0,
            .maximumAcceleration = std::numeric_limits<double>::infinity(),
            .maximumJerk = 4.0,
        });

        require(!unequalVelocity.has_value(), "an unbounded-acceleration cruise should reject unequal boundary velocities");
        require(unequalVelocity.error().code == path_tempo::PlanningErrorCode::InvalidInput, "unequal cruise velocities should return an input error");

        const auto nonzeroAcceleration = planner.solve({
            .piece = 26,
            .length = 3.0,
            .beginning = {.velocity = 1.5, .acceleration = std::numeric_limits<double>::denorm_min()},
            .ending = {.velocity = 1.5},
            .maximumVelocity = 2.0,
            .maximumAcceleration = std::numeric_limits<double>::infinity(),
            .maximumJerk = 4.0,
        });

        require(!nonzeroAcceleration.has_value(), "an unbounded-acceleration cruise should reject nonzero boundary acceleration");
        require(nonzeroAcceleration.error().code == path_tempo::PlanningErrorCode::InvalidInput, "nonzero cruise acceleration should return an input error");

        const auto infiniteDuration = planner.solve({
            .piece = 27,
            .length = std::numeric_limits<double>::max(),
            .beginning = {.velocity = std::numeric_limits<double>::denorm_min()},
            .ending = {.velocity = std::numeric_limits<double>::denorm_min()},
            .maximumVelocity = 1.0,
            .maximumAcceleration = std::numeric_limits<double>::infinity(),
            .maximumJerk = 4.0,
        });

        require(!infiniteDuration.has_value(), "an unbounded-acceleration cruise should reject a non-finite duration");
        require(infiniteDuration.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a non-finite cruise duration should return an input error");

        const auto zeroDuration = planner.solve({
            .piece = 29,
            .length = std::numeric_limits<double>::denorm_min(),
            .beginning = {.velocity = std::numeric_limits<double>::max()},
            .ending = {.velocity = std::numeric_limits<double>::max()},
            .maximumVelocity = std::numeric_limits<double>::max(),
            .maximumAcceleration = std::numeric_limits<double>::infinity(),
            .maximumJerk = 4.0,
        });

        require(!zeroDuration.has_value(), "an unbounded-acceleration cruise should reject an underflowed zero duration");
        require(zeroDuration.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a zero cruise duration should return an input error");
    }

    void testSubPicosecondCruisePhase() {
        constexpr double length = 5e-13;
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 28,
            .length = length,
            .beginning = {.velocity = 1.0},
            .ending = {.velocity = 1.0},
            .maximumVelocity = 2.0,
            .maximumAcceleration = 1.0,
            .maximumJerk = 4.0,
        });

        require(transition.has_value(), "a positive sub-picosecond cruise phase should be retained");
        require(!transition->empty(), "a positive sub-picosecond cruise should emit a segment");
        require(transition->duration() > 0.0 && transition->duration() <= 1e-12, "the retained cruise should preserve its short positive duration");
        require(std::abs(transition->back().position(transition->back().duration) - length) <= std::numeric_limits<double>::epsilon() * length, "the retained cruise should preserve its endpoint");
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

    void testUnrecoverableTransitionHasDetailedError() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 30,
            .length = 1e-3,
            .beginning = {.velocity = 1.0, .acceleration = 1.0},
            .ending = {.velocity = 1.0, .acceleration = -1.0},
            .maximumVelocity = 1.0,
            .maximumAcceleration = 1.0,
            .maximumJerk = 1.0,
        });

        require(!transition.has_value(), "an infeasible public transition should fail");
        require(!transition.error().message.empty(), "an unrecoverable public transition failure should retain a detailed message");
        require(transition.error().message.contains("Ruckig"), "an unrecoverable public transition failure should identify the underlying solver");
        require(transition.error().message.contains("length="), "an unrecoverable public transition failure should include request context");
    }

    void testTransitionLimitValidation() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto overLimitCruise = planner.solve({
            .piece = 22,
            .length = 3.0,
            .beginning = {.velocity = 1.5},
            .ending = {.velocity = 1.5},
            .maximumVelocity = 1.0,
            .maximumAcceleration = std::numeric_limits<double>::infinity(),
            .maximumJerk = 4.0,
        });

        require(!overLimitCruise.has_value(), "an unbounded-acceleration cruise above its velocity limit should be rejected");
        require(overLimitCruise.error().code == path_tempo::PlanningErrorCode::InvalidInput, "an over-limit cruise should return an input error");

        const auto infiniteVelocity = planner.solve({
            .piece = 23,
            .length = 3.0,
            .beginning = {},
            .ending = {},
            .maximumVelocity = std::numeric_limits<double>::infinity(),
            .maximumAcceleration = 1.0,
            .maximumJerk = 4.0,
        });

        require(!infiniteVelocity.has_value(), "an infinite velocity limit should be rejected before it disables transition validation");
        require(infiniteVelocity.error().code == path_tempo::PlanningErrorCode::InvalidInput, "an infinite velocity limit should return an input error");
    }

    void testShortTransitionVelocityTolerance() {
        path_tempo::ScalarTransitionPlanner planner;
        const auto transition = planner.solve({
            .piece = 24,
            .length = 0.01,
            .beginning = {.velocity = 2500.0},
            .ending = {},
            .maximumVelocity = 10000.0,
            .maximumAcceleration = 1e10,
            .maximumJerk = 1e16,
        });

        require(transition.has_value(), "a short transition should use velocity tolerance for tiny negative phase-boundary velocity roundoff");
    }

    void testVelocityReachability() {
        const auto transitionDistance = path_tempo::velocityTransitionDistance(0.0, 2.0, 1.0, 2.0);

        require(std::abs(transitionDistance - 2.5) <= ANALYTIC_RESULT_TOLERANCE, "velocity transition distance should include jerk-limited ramps");
        require(path_tempo::reachableVelocity(1.0, 0.5, 10.0, 1.0, 2.0) == 0.5, "a lower velocity cap should be immediately reachable");

        const auto reachable = path_tempo::reachableVelocity(0.0, 2.0, 1.0, 1.0, 2.0);

        require(reachable > 0.0 && reachable < 2.0, "a short piece should reduce its reachable velocity cap");
        require(path_tempo::velocityTransitionDistance(0.0, reachable, 1.0, 2.0) <= 1.0, "reachable velocity should fit within the available distance");
    }

    path_tempo::Limits<3> testPathLimits() {
        return {
            .pathAcceleration = 2.0,
            .pathJerk = 5.0,
            .coordinateVelocity = {4.0, 4.0, 4.0},
            .coordinateAcceleration = {2.0, 2.0, 2.0},
            .coordinateJerk = {5.0, 5.0, 5.0},
        };
    }

    void testPathInputValidation() {
        const auto sampled = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {0.0, 0.0, 0.0},
            .to = {3.0, 0.0, 0.0},
        }, 30, 2.0);

        require(sampled.has_value(), "path input validation should have a valid sampled line");

        auto nanLimits = *sampled;
        nanLimits.initialLimits.acceleration = std::numeric_limits<double>::quiet_NaN();
        const std::array nanPieces {nanLimits.view()};
        path_tempo::PathPlanner planner;
        const auto nanResult = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = nanPieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        });

        require(!nanResult.has_value(), "a NaN initial piece limit should be rejected");
        require(nanResult.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a NaN initial piece limit should return an input error");

        const std::array duplicatePieces {sampled->view(), sampled->view()};
        const auto duplicateResult = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = duplicatePieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        });

        require(!duplicateResult.has_value(), "duplicate path piece IDs should be rejected");
        require(duplicateResult.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a duplicate path piece ID should return an input error");

        auto tolerantStations = *sampled;
        tolerantStations.stations.front().distance = -5e-13;
        tolerantStations.stations.back().distance = tolerantStations.length + 5e-11;
        const std::array tolerantPieces {tolerantStations.view()};
        const auto tolerantResult = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = tolerantPieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {.boundaryAccelerationMode = path_tempo::BoundaryAccelerationMode::Zero},
        });

        require(tolerantResult.has_value(), "station endpoints within the documented coverage tolerance should be normalized and accepted");

        const std::array validPieces {sampled->view()};
        const auto invalidMode = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = validPieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {.boundaryAccelerationMode = static_cast<path_tempo::BoundaryAccelerationMode>(127)},
        });

        require(!invalidMode.has_value(), "an invalid boundary-acceleration mode should be rejected");
        require(invalidMode.error().code == path_tempo::PlanningErrorCode::InvalidInput, "an invalid boundary-acceleration mode should return an input error");
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
        require(planned->pieceLimits.size() == pieces.size(), "multi-piece planning should report one effective limit set per piece");
        require(planned->pieceLimits[0].velocity <= 3.0 && planned->pieceLimits[1].velocity <= 1.5, "effective piece limits should retain input maximum velocities");
        require(planned->pieceBoundaries[1].velocity > 0.0, "continuous planning should retain motion at an internal boundary");
        require(planned->pieceBoundaries[1].velocity <= 1.5 + BOUNDARY_LIMIT_TOLERANCE, "an internal boundary should honor both piece maximum velocities");
        require(planned->diagnostics.trajectoryDuration > 0.0, "multi-piece planning should report its materialized duration");

        auto previousDistance = 0.0;
        auto previousVelocity = 0.0;
        auto previousAcceleration = 0.0;
        auto sawFirstPiece = false;
        auto sawSecondPiece = false;

        for (const auto &segment : planned->timeLaw.segments) {
            require(std::abs(segment.position(0.0) - previousDistance) <= TIME_LAW_CONTINUITY_TOLERANCE, "multi-piece cubics should be globally position-continuous");
            require(std::abs(segment.velocity(0.0) - previousVelocity) <= TIME_LAW_CONTINUITY_TOLERANCE, "multi-piece cubics should be velocity-continuous");
            require(std::abs(segment.acceleration(0.0) - previousAcceleration) <= TIME_LAW_CONTINUITY_TOLERANCE, "multi-piece cubics should be acceleration-continuous");
            require(std::abs(segment.jerk()) <= 5.0 + BOUNDARY_LIMIT_TOLERANCE, "multi-piece cubics should obey the scalar and coordinate jerk limit");

            const auto pieceMaxVelocity = segment.piece == 31 ? 3.0 : 1.5;
            auto maximumSegmentVelocity = std::max(segment.velocity(0.0), segment.velocity(segment.duration));

            if (std::abs(segment.jerk()) > JERK_ZERO_TOLERANCE) {
                const auto velocityExtremum = -segment.acceleration(0.0) / segment.jerk();

                if (velocityExtremum > 0.0 && velocityExtremum < segment.duration) {
                    maximumSegmentVelocity = std::max(maximumSegmentVelocity, segment.velocity(velocityExtremum));
                }
            }

            require(maximumSegmentVelocity <= pieceMaxVelocity + SOLVER_ENDPOINT_TOLERANCE, "multi-piece cubics should obey each piece maximum velocity");
            require(std::abs(segment.acceleration(0.0)) <= 2.0 + SOLVER_ENDPOINT_TOLERANCE && std::abs(segment.acceleration(segment.duration)) <= 2.0 + SOLVER_ENDPOINT_TOLERANCE, "multi-piece cubics should obey the scalar and coordinate acceleration limit");

            previousDistance = segment.position(segment.duration);
            previousVelocity = segment.velocity(segment.duration);
            previousAcceleration = segment.acceleration(segment.duration);
            sawFirstPiece = sawFirstPiece || segment.piece == 31;
            sawSecondPiece = sawSecondPiece || segment.piece == 32;
        }

        require(sawFirstPiece && sawSecondPiece, "multi-piece segments should preserve each piece identity");
        require(std::abs(previousDistance - 20.0) <= TIME_LAW_CONTINUITY_TOLERANCE, "multi-piece planning should cover the complete path distance");
        require(std::abs(previousVelocity) <= TIME_LAW_CONTINUITY_TOLERANCE, "rest-to-rest multi-piece planning should end at rest");
        require(std::abs(previousAcceleration) <= TIME_LAW_CONTINUITY_TOLERANCE, "rest-to-rest multi-piece planning should end at zero acceleration");

        const auto replanned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        });

        require(replanned.has_value(), "a repeated path solve should succeed");
        require(!replanned->timeLaw.segments.empty(), "a repeated path solve should produce a scalar time law");
    }

    void testMultiPiecePathUsesNonzeroBoundaryAcceleration() {
        const auto first = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {0.0, 0.0, 0.0},
            .to = {1.0, 0.0, 0.0},
        }, 111, 10.0);
        const auto second = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {1.0, 0.0, 0.0},
            .to = {2.0, 0.0, 0.0},
        }, 112, 10.0);
        const auto third = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {2.0, 0.0, 0.0},
            .to = {3.0, 0.0, 0.0},
        }, 113, 10.0);
        const auto fourth = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {3.0, 0.0, 0.0},
            .to = {4.0, 0.0, 0.0},
        }, 114, 10.0);

        require(first.has_value() && second.has_value() && third.has_value() && fourth.has_value(), "nonzero-boundary test lines should sample");

        const std::array pieces {first->view(), second->view(), third->view(), fourth->view()};
        const path_tempo::Limits<3> limits {
            .pathAcceleration = 2.0,
            .pathJerk = 4.0,
            .coordinateVelocity = {10.0, 10.0, 10.0},
            .coordinateAcceleration = {2.0, 2.0, 2.0},
            .coordinateJerk = {4.0, 4.0, 4.0},
        };
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = {.applySampledCorrections = false},
        });

        require(planned.has_value(), "a collinear four-piece path should solve with continuous boundary acceleration");
        require(std::ranges::any_of(planned->pieceBoundaries.begin() + 1, planned->pieceBoundaries.end() - 1,
            [](const path_tempo::BoundaryState &state) { return std::abs(state.acceleration) > 1e-6; }),
            "a collinear multi-piece path should retain nonzero acceleration across artificial piece boundaries");

        auto boundaryIndex = std::size_t {1};

        for (std::size_t segmentIndex = 0; segmentIndex + 1 < planned->timeLaw.segments.size(); ++segmentIndex) {
            const auto &left = planned->timeLaw.segments[segmentIndex];
            const auto &right = planned->timeLaw.segments[segmentIndex + 1];

            if (left.piece == right.piece) {
                continue;
            }

            const auto &boundary = planned->pieceBoundaries[boundaryIndex++];
            require(std::abs(left.velocity(left.duration) - boundary.velocity) <= SOLVER_ENDPOINT_TOLERANCE,
                "the preceding piece should reach the proposed nonzero-acceleration boundary velocity");
            require(std::abs(left.acceleration(left.duration) - boundary.acceleration) <= SOLVER_ENDPOINT_TOLERANCE,
                "the preceding piece should reach the proposed nonzero boundary acceleration");
            require(std::abs(right.velocity(0.0) - boundary.velocity) <= SOLVER_ENDPOINT_TOLERANCE,
                "the following piece should begin at the proposed nonzero-acceleration boundary velocity");
            require(std::abs(right.acceleration(0.0) - boundary.acceleration) <= SOLVER_ENDPOINT_TOLERANCE,
                "the following piece should begin at the proposed nonzero boundary acceleration");
        }

        require(boundaryIndex == pieces.size(), "the nonzero-boundary test should materialize every piece boundary");

        const auto zeroBoundary = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = {
                .applySampledCorrections = false,
                .boundaryAccelerationMode = path_tempo::BoundaryAccelerationMode::Zero,
            },
        });

        require(zeroBoundary.has_value(), "the explicit zero-boundary mode should solve");
        require(std::ranges::all_of(zeroBoundary->pieceBoundaries.begin() + 1, zeroBoundary->pieceBoundaries.end() - 1,
            [](const path_tempo::BoundaryState &state) { return std::abs(state.acceleration) <= TIME_LAW_CONTINUITY_TOLERANCE; }),
            "zero-boundary mode should force every internal acceleration to zero");

        const auto nonzeroBoundaryDuration = std::accumulate(planned->timeLaw.segments.begin(), planned->timeLaw.segments.end(), 0.0,
            [](const double total, const path_tempo::CubicTimeSegment &segment) { return total + segment.duration; });
        const auto zeroBoundaryDuration = std::accumulate(zeroBoundary->timeLaw.segments.begin(), zeroBoundary->timeLaw.segments.end(), 0.0,
            [](const double total, const path_tempo::CubicTimeSegment &segment) { return total + segment.duration; });
        require(nonzeroBoundaryDuration + TIME_LAW_CONTINUITY_TOLERANCE < zeroBoundaryDuration,
            "nonzero acceleration across artificial boundaries should improve on the zero-acceleration profile");
    }

    void testBoundaryAccelerationCandidatePreference() {
        constexpr double durationResolution = 1e-6;
        constexpr double microDuration = 0.023501909228039557;
        constexpr double zeroDuration = 0.023501909657201489;

        require(path_tempo::detail::preferBoundaryAccelerationCandidate(
            zeroDuration, 0.0, microDuration, 1.741278423776903e-6, durationResolution),
            "duration-equivalent zero acceleration should replace a micro-acceleration state");
        require(path_tempo::detail::preferBoundaryAccelerationCandidate(
            microDuration - 2.0 * durationResolution, 0.25, microDuration, 0.0, durationResolution),
            "a meaningfully faster nonzero acceleration should remain preferable");
        require(!path_tempo::detail::preferBoundaryAccelerationCandidate(
            microDuration + 2.0 * durationResolution, 0.0, microDuration, 0.25, durationResolution),
            "a meaningfully slower zero acceleration should not replace a useful nonzero state");
        require(path_tempo::detail::resolveBoundaryRepairWeight(
            1.0309278350515464e-6, 2.0624815246445016, 95.95, durationResolution) == 0.0,
            "a repaired boundary weight that creates only a nano-phase should collapse to the zero-state baseline");
        require(path_tempo::detail::resolveBoundaryRepairWeight(
            0.1, 2.0624815246445016, 95.95, durationResolution) == 0.1,
            "a boundary weight with a resolvable acceleration ramp should be retained");
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

    void testCurvatureContinuityUsesGeometricScale() {
        const auto solveBoundary = [](const double length, const double previousCurvature, const double currentCurvature) {
            const std::array previousStations {
                path_tempo::DifferentialStation<2> {.distance = 0.0, .tangent = {1.0, 0.0}, .curvature = {0.0, previousCurvature}},
                path_tempo::DifferentialStation<2> {.distance = length, .tangent = {1.0, 0.0}, .curvature = {0.0, previousCurvature}},
            };
            const std::array currentStations {
                path_tempo::DifferentialStation<2> {.distance = 0.0, .tangent = {1.0, 0.0}, .curvature = {0.0, currentCurvature}},
                path_tempo::DifferentialStation<2> {.distance = length, .tangent = {1.0, 0.0}, .curvature = {0.0, currentCurvature}},
            };
            const std::array pieces {
                path_tempo::PathPiece<2> {.id = 43, .length = length, .maxVelocity = 1.0, .initialLimits = {}, .stations = previousStations},
                path_tempo::PathPiece<2> {.id = 44, .length = length, .maxVelocity = 1.0, .initialLimits = {}, .stations = currentStations},
            };
            path_tempo::PathPlanner planner;

            return planner.solve(path_tempo::PathPlanningRequest<2> {
                .pieces = pieces,
                .beginning = {},
                .ending = {},
                .limits = {
                    .pathAcceleration = 1e6,
                    .pathJerk = 1e9,
                    .coordinateVelocity = {1.0, 1.0},
                    .coordinateAcceleration = {1e6, 1e6},
                    .coordinateJerk = {1e9, 1e9},
                },
                .settings = {.applySampledCorrections = false},
            });
        };

        const auto baseUnits = solveBoundary(1.0, 1.0, 1.0 + 5e-9);
        const auto smallerUnits = solveBoundary(1e-3, 1e3, 1e3 * (1.0 + 5e-9));
        const auto nearZeroNoise = solveBoundary(1.0, 0.0, 2.5e-10);

        require(baseUnits.has_value() && smallerUnits.has_value(), "equivalent relative curvature noise should be accepted in either geometric scale");
        require(nearZeroNoise.has_value(), "near-zero spline curvature noise should be accepted");

        const auto smallJump = solveBoundary(1.0, 0.0, 5e-9);

        require(!smallJump.has_value(), "a small but meaningful curvature jump should not be hidden by a fixed absolute tolerance");
        require(smallJump.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a curvature discontinuity should be an input error");
    }

    path_tempo::SampledPathPiece<3> sampledUnitCircleInterval(const path_tempo::PathPieceId id, const double fromAngle, const double toAngle, const double maxVelocity, const std::size_t intervals = 32) {
        const auto length = toAngle - fromAngle;
        path_tempo::SampledPathPiece<3> piece {
            .id = id,
            .length = length,
            .maxVelocity = maxVelocity,
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

    path_tempo::SampledPathPiece<3> sampledQuarterCircle(const path_tempo::PathPieceId id, const double maxVelocity, const std::size_t intervals = 32) {
        return sampledUnitCircleInterval(id, 0.0, std::numbers::pi / 2.0, maxVelocity, intervals);
    }

    void testCurvedPathPlanning() {
        const auto curve = sampledQuarterCircle(51, 4.0);
        const std::array pieces {curve.view()};
        path_tempo::PathPlanner planner;
        const path_tempo::Limits<3> limits {
            .pathAcceleration = 0.8,
            .pathJerk = 1.2,
            .coordinateVelocity = {4.0, 4.0, 4.0},
            .coordinateAcceleration = {0.8, 0.8, 0.8},
            .coordinateJerk = {1.2, 1.2, 1.2},
        };
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = {
                .maximumCorrectionPasses = 12,
            },
        });

        require(planned.has_value(), "a sampled curved path should solve");
        require(!planned->timeLaw.segments.empty(), "curved planning should produce a scalar time law");
        require(planned->diagnostics.correctionPasses >= 1, "curved planning should report its bounded correction passes");
        require(planned->diagnostics.correctedPieces == 1, "coupled curved constraints should trigger a local correction for this path");

        for (const auto &station : curve.stations) {
            const auto tolerance = std::max(STATION_MATCH_ABSOLUTE_TOLERANCE, curve.length * STATION_MATCH_RELATIVE_TOLERANCE);
            auto visited = false;

            for (const auto &segment : planned->timeLaw.segments) {
                const auto localFrom = segment.position(0.0);
                const auto localTo = segment.position(segment.duration);

                if (station.distance < localFrom - tolerance || station.distance > localTo + tolerance) {
                    continue;
                }

                visited = true;
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
                    require(std::abs(station.tangent[axis] * velocity) <= limits.coordinateVelocity[axis] * (1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE), "curved timing should obey sampled coordinate velocity");
                    require(std::abs(coupledAcceleration) <= limits.coordinateAcceleration[axis] * (1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE), "curved timing should obey sampled coupled coordinate acceleration");
                    require(std::abs(coupledJerk) <= limits.coordinateJerk[axis] * (1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE), "curved timing should obey sampled coupled coordinate jerk");
                }

                require(std::sqrt(coupledAccelerationSquared) <= limits.pathAcceleration * (1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE), "curved timing should obey sampled coupled path acceleration");
                require(std::sqrt(coupledJerkSquared) <= limits.pathJerk * (1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE), "curved timing should obey sampled coupled path jerk");
            }

            require(visited, "every curved differential station should match at least one scalar phase");
        }
    }

    void testCoupledStationsAtTransitionBoundaries() {
        auto curve = sampledQuarterCircle(52, 4.0, 4);
        const std::array referencePieces {curve.view()};
        const auto unbounded = std::numeric_limits<double>::infinity();
        const path_tempo::Limits<3> limits {
            .pathAcceleration = 0.8,
            .pathJerk = 1.2,
            .coordinateVelocity = {unbounded, unbounded, unbounded},
            .coordinateAcceleration = {unbounded, unbounded, unbounded},
            .coordinateJerk = {unbounded, unbounded, unbounded},
        };
        const path_tempo::PathPlanningSettings settings {
            .maximumCorrectionPasses = 12,
            .boundaryAccelerationMode = path_tempo::BoundaryAccelerationMode::Zero,
        };
        path_tempo::PathPlanner referencePlanner;
        const auto reference = referencePlanner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = referencePieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = settings,
        });

        require(reference.has_value(), "the phase-boundary station reference path should solve");

        for (const auto &segment : reference->timeLaw.segments) {
            const auto distance = segment.position(segment.duration);

            if (distance <= 0.0 || distance >= curve.length) {
                continue;
            }

            const path_tempo::DifferentialStation<3> station {
                .distance = distance,
                .tangent = {-std::sin(distance), std::cos(distance), 0.0},
                .curvature = {-std::cos(distance), -std::sin(distance), 0.0},
                .thirdDerivative = {std::sin(distance), -std::cos(distance), 0.0},
            };
            curve.stations.push_back(station);
        }

        std::ranges::sort(curve.stations, {}, &path_tempo::DifferentialStation<3>::distance);
        const std::array boundaryPieces {curve.view()};
        path_tempo::PathPlanner boundaryPlanner;
        const auto boundary = boundaryPlanner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = boundaryPieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = settings,
        });

        require(boundary.has_value(), "coupled stations at scalar phase boundaries should remain covered");

        auto duplicated = curve;
        duplicated.stations.insert(duplicated.stations.end(), curve.stations.begin(), curve.stations.end());
        std::ranges::sort(duplicated.stations, {}, &path_tempo::DifferentialStation<3>::distance);
        const std::array duplicatedPieces {duplicated.view()};
        path_tempo::PathPlanner duplicatedPlanner;
        const auto enriched = duplicatedPlanner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = duplicatedPieces,
            .beginning = {},
            .ending = {},
            .limits = limits,
            .settings = settings,
        });

        require(enriched.has_value(), "duplicate coupled stations at scalar phase boundaries should remain covered");
        require(enriched->diagnostics.correctionPasses == boundary->diagnostics.correctionPasses, "duplicate phase-boundary stations should preserve correction convergence");
        require(std::abs(enriched->diagnostics.trajectoryDuration - boundary->diagnostics.trajectoryDuration) <= SOLVER_ENDPOINT_TOLERANCE, "duplicate phase-boundary stations should preserve trajectory duration");
    }

    void testAcceptedTransitionEndpointChecksCoupledStation() {
        constexpr double length = 0.019140457806488653;
        constexpr double maximumVelocity = 2.2670593417302345;
        constexpr double maximumAcceleration = 56048.752530983584;
        constexpr double maximumJerk = 13212749508.561106;
        constexpr double thirdDerivative = 1.0e9;
        const std::array stations {
            path_tempo::DifferentialStation<1> {
                .distance = 0.0,
                .tangent = {1.0},
                .thirdDerivative = {thirdDerivative},
            },
            path_tempo::DifferentialStation<1> {
                .distance = length,
                .tangent = {1.0},
                .thirdDerivative = {thirdDerivative},
            },
        };
        const std::array pieces {
            path_tempo::PathPiece<1> {
                .id = 52,
                .length = length,
                .maxVelocity = maximumVelocity,
                .initialLimits = {
                    .velocity = maximumVelocity,
                    .acceleration = maximumAcceleration,
                    .jerk = maximumJerk,
                },
                .stations = stations,
            },
        };
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(path_tempo::PathPlanningRequest<1> {
            .pieces = pieces,
            .beginning = {.velocity = 0.75667090528164205},
            .ending = {.velocity = 0.5804097813273188},
            .limits = {
                .pathAcceleration = maximumAcceleration,
                .pathJerk = maximumJerk,
                .coordinateVelocity = {maximumVelocity},
                .coordinateAcceleration = {maximumAcceleration},
                .coordinateJerk = {maximumJerk},
            },
            .settings = {.boundaryAccelerationMode = path_tempo::BoundaryAccelerationMode::Zero},
        });

        require(planned.has_value(), "an endpoint within the accepted transition tolerance should remain covered by coupled checking");
        require(planned->diagnostics.correctedPieces == 1, "a coupled violation at a slightly undershot endpoint should trigger correction");

        const auto &last = planned->timeLaw.segments.back();
        const auto velocity = last.velocity(last.duration);
        const auto coupledJerk = last.jerk() + thirdDerivative * velocity * velocity * velocity;

        require(std::abs(coupledJerk) <= maximumJerk * (1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE), "the corrected endpoint should satisfy its coupled jerk limit");
    }

    void testShortJerkPhaseUsesBracketingDifferentialStations() {
        constexpr double length = 0.009551481852077569;
        constexpr double velocity = 0.35;
        constexpr double pathJerk = 101.0;
        const std::array stations {
            path_tempo::DifferentialStation<2> {
                .distance = 0.0,
                .tangent = {1.0, 0.0},
                .thirdDerivative = {0.0, 221.07195407750322},
            },
            path_tempo::DifferentialStation<2> {
                .distance = 0.008954514236322719,
                .tangent = {1.0, 0.0},
                .thirdDerivative = {0.0, 1205.7670652463494},
            },
            path_tempo::DifferentialStation<2> {
                .distance = length,
                .tangent = {1.0, 0.0},
                .thirdDerivative = {0.0, 221.07195407750322},
            },
        };
        const std::array pieces {
            path_tempo::PathPiece<2> {
                .id = 71,
                .length = length,
                .maxVelocity = 0.5,
                .initialLimits = {
                    .velocity = 0.5,
                    .acceleration = 4.845,
                    .jerk = 95.95,
                },
                .stations = stations,
            },
        };
        auto request = path_tempo::PathPlanningRequest<2> {
            .pieces = pieces,
            .beginning = {.velocity = velocity},
            .ending = {.velocity = velocity, .acceleration = 0.03196919167034261},
            .limits = {
                .pathAcceleration = 5.1,
                .pathJerk = pathJerk,
                .coordinateVelocity = {10.0, 10.0},
                .coordinateAcceleration = {10.0, 10.0},
                .coordinateJerk = {pathJerk, pathJerk},
            },
            .settings = {
                .maximumCorrectionPasses = 32,
                .applySampledCorrections = true,
                .boundaryAccelerationMode = path_tempo::BoundaryAccelerationMode::Optimized,
            },
        };
        auto callbackCalls = std::size_t {0};
        const auto correction = [&](const path_tempo::PlannedPath &candidate)
                -> std::expected<std::vector<path_tempo::PieceCorrection>, std::string> {
            ++callbackCalls;
            auto required = 1.0;

            for (const auto &segment : candidate.timeLaw.segments) {
                for (const auto time : std::array {0.0, segment.duration}) {
                    const auto distance = std::clamp(segment.position(time), 0.0, length);
                    const auto upper = std::ranges::lower_bound(stations, distance, {},
                        &path_tempo::DifferentialStation<2>::distance);
                    const auto check = [&](const auto &station) {
                        const auto scalarVelocity = segment.velocity(time);
                        const auto geometricJerk = station.thirdDerivative[1]
                            * scalarVelocity * scalarVelocity * scalarVelocity;
                        const auto coupledJerk = std::hypot(segment.jerk(), geometricJerk);
                        required = std::max(required, std::cbrt(coupledJerk / pathJerk));
                    };

                    if (upper == stations.begin()) {
                        check(*upper);
                    } else if (upper == stations.end()) {
                        check(stations.back());
                    } else {
                        check(*std::prev(upper));
                        check(*upper);
                    }
                }
            }

            if (required > 1.0 + SAMPLED_LIMIT_RELATIVE_TOLERANCE) {
                return std::vector {path_tempo::PieceCorrection {
                    .piece = 71,
                    .requiredTimeScale = required * 1.01,
                }};
            }

            return std::vector<path_tempo::PieceCorrection> {};
        };
        request.settings.applySampledCorrections = false;
        path_tempo::PathPlanner materializationOnlyPlanner;
        const auto materializationOnly = materializationOnlyPlanner.solve(request, correction);

        require(materializationOnly.has_value(), "materialization should correct an otherwise missed short jerk phase");
        require(callbackCalls > 1, "station-only checking should leave the short jerk phase for materialization");

        request.settings.applySampledCorrections = true;
        callbackCalls = 0;
        path_tempo::PathPlanner planner;
        const auto planned = planner.solve(request, correction);

        require(planned.has_value(), "a short jerk phase between differential stations should solve");
        require(callbackCalls == 1, "sampled checking should catch a short jerk phase before materialization");
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
                .coordinateVelocity = {2.0, 2.0, 2.0},
                .coordinateAcceleration = {1.0, 1.0, 1.0},
                .coordinateJerk = {0.1, 0.1, 0.1},
            },
            .settings = {},
        });

        require(!planned.has_value(), "a fixed curved boundary velocity above the geometric jerk cap should be rejected");
        require(planned.error().code == path_tempo::PlanningErrorCode::InvalidInput, "an infeasible curved boundary should be an input error");
    }

    void testTinyNonzeroGeometryConstrainsBoundaryVelocity() {
        const std::array curvatureStations {
            path_tempo::DifferentialStation<2> {.distance = 0.0, .tangent = {1.0, 0.0}, .curvature = {0.0, 1e-16}},
            path_tempo::DifferentialStation<2> {.distance = 1.0, .tangent = {1.0, 0.0}, .curvature = {0.0, 1e-16}},
        };
        const path_tempo::PathPiece<2> curvaturePiece {
            .id = 62,
            .length = 1.0,
            .maxVelocity = 1e9,
            .initialLimits = {},
            .stations = curvatureStations,
        };
        const std::array curvaturePieces {curvaturePiece};
        path_tempo::PathPlanner planner;
        const auto curvatureLimited = planner.solve(path_tempo::PathPlanningRequest<2> {
            .pieces = curvaturePieces,
            .beginning = {.velocity = 1e9},
            .ending = {.velocity = 1e9},
            .limits = {
                .pathAcceleration = 1e300,
                .pathJerk = 1e300,
                .coordinateVelocity = {1e300, 1e300},
                .coordinateAcceleration = {1e300, 1.0},
                .coordinateJerk = {1e300, 1e300},
            },
            .settings = {},
        });

        require(!curvatureLimited.has_value(), "tiny nonzero curvature should constrain fixed boundary velocity");
        require(curvatureLimited.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a boundary above a tiny-curvature velocity cap should be an input error");

        const std::array thirdDerivativeStations {
            path_tempo::DifferentialStation<2> {.distance = 0.0, .tangent = {1.0, 0.0}, .thirdDerivative = {0.0, 1e-18}},
            path_tempo::DifferentialStation<2> {.distance = 1.0, .tangent = {1.0, 0.0}, .thirdDerivative = {0.0, 1e-18}},
        };
        const path_tempo::PathPiece<2> thirdDerivativePiece {
            .id = 63,
            .length = 1.0,
            .maxVelocity = 1e7,
            .initialLimits = {},
            .stations = thirdDerivativeStations,
        };
        const std::array thirdDerivativePieces {thirdDerivativePiece};
        const auto jerkLimited = planner.solve(path_tempo::PathPlanningRequest<2> {
            .pieces = thirdDerivativePieces,
            .beginning = {.velocity = 1e7},
            .ending = {.velocity = 1e7},
            .limits = {
                .pathAcceleration = 1e300,
                .pathJerk = 1e300,
                .coordinateVelocity = {1e300, 1e300},
                .coordinateAcceleration = {1e300, 1e300},
                .coordinateJerk = {1e300, 1.0},
            },
            .settings = {},
        });

        require(!jerkLimited.has_value(), "tiny nonzero third derivative should constrain fixed boundary velocity");
        require(jerkLimited.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a boundary above a tiny-third-derivative velocity cap should be an input error");
    }

    void testMultiPieceCurvedPlanning() {
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
                .coordinateVelocity = {3.0, 3.0, 3.0},
                .coordinateAcceleration = {1.0, 1.0, 1.0},
                .coordinateJerk = {2.0, 2.0, 2.0},
            },
            .settings = {},
        });

        require(planned.has_value(), "a C2 two-piece curved chain should solve");
        require(planned->pieceBoundaries[1].velocity > 0.5, "a curved velocity-cap boundary should not collapse to a stop");
        require(std::abs(planned->pieceBoundaries[1].acceleration) <= 1.0 + BOUNDARY_LIMIT_TOLERANCE,
            "a curved internal boundary should obey the shared scalar acceleration limit");

        auto previousDistance = 0.0;
        auto previousVelocity = 0.0;
        auto previousAcceleration = 0.0;

        for (const auto &segment : planned->timeLaw.segments) {
            require(std::abs(segment.position(0.0) - previousDistance) <= TIME_LAW_CONTINUITY_TOLERANCE, "refined curved cubics should be position-continuous");
            require(std::abs(segment.velocity(0.0) - previousVelocity) <= TIME_LAW_CONTINUITY_TOLERANCE, "refined curved cubics should be velocity-continuous");
            require(std::abs(segment.acceleration(0.0) - previousAcceleration) <= TIME_LAW_CONTINUITY_TOLERANCE, "refined curved cubics should be acceleration-continuous");
            previousDistance = segment.position(segment.duration);
            previousVelocity = segment.velocity(segment.duration);
            previousAcceleration = segment.acceleration(segment.duration);
        }

        const auto repeated = planner.solve(path_tempo::PathPlanningRequest<3> {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = {
                .pathAcceleration = 1.0,
                .pathJerk = 2.0,
                .coordinateVelocity = {3.0, 3.0, 3.0},
                .coordinateAcceleration = {1.0, 1.0, 1.0},
                .coordinateJerk = {2.0, 2.0, 2.0},
            },
            .settings = {},
        });

        require(repeated.has_value(), "a repeated curved path should solve");
        require(!repeated->timeLaw.segments.empty(), "a repeated curved path should produce a scalar time law");
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
                .coordinateVelocity = {2.0, 2.0, 2.0},
                .coordinateAcceleration = {2.0, 2.0, 2.0},
                .coordinateJerk = {4.0, 4.0, 4.0},
            },
            .settings = {},
        });

        require(planned.has_value(), "a curved path with nonzero boundary PVA should solve");
        require(std::abs(planned->timeLaw.segments.front().velocity(0.0) - beginning.velocity) <= BOUNDARY_LIMIT_TOLERANCE, "curved planning should preserve beginning velocity");
        require(std::abs(planned->timeLaw.segments.front().acceleration(0.0) - beginning.acceleration) <= BOUNDARY_LIMIT_TOLERANCE, "curved planning should preserve beginning acceleration");

        const auto &last = planned->timeLaw.segments.back();
        require(std::abs(last.velocity(last.duration) - ending.velocity) <= SOLVER_ENDPOINT_TOLERANCE, "curved planning should preserve ending velocity");
        require(std::abs(last.acceleration(last.duration) - ending.acceleration) <= SOLVER_ENDPOINT_TOLERANCE, "curved planning should preserve ending acceleration");
    }

    void testMaterializationCorrectionCallback() {
        const auto firstLine = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {0.0, 0.0, 0.0},
            .to = {10.0, 0.0, 0.0},
        }, 91, 3.0);
        const auto secondLine = path_tempo::sampleLine(path_tempo::Line<3> {
            .from = {10.0, 0.0, 0.0},
            .to = {20.0, 0.0, 0.0},
        }, 92, 3.0);

        require(firstLine.has_value() && secondLine.has_value(), "materialization-correction test lines should sample");

        const std::array pieces {firstLine->view(), secondLine->view()};
        const path_tempo::PathPlanningRequest<3> request {
            .pieces = pieces,
            .beginning = {},
            .ending = {},
            .limits = testPathLimits(),
            .settings = {},
        };
        path_tempo::PathPlanner baselinePlanner;
        const auto baseline = baselinePlanner.solve(request);

        require(baseline.has_value(), "materialization-correction baseline should solve");

        auto callbackCalls = std::size_t {0};
        const path_tempo::MaterializationCorrection correction = [&](const auto &candidate) -> std::expected<std::vector<path_tempo::PieceCorrection>, std::string> {
            ++callbackCalls;
            require(!candidate.timeLaw.segments.empty(), "materialization callback should receive a complete candidate");

            if (callbackCalls == 1) {
                return std::vector {path_tempo::PieceCorrection {.piece = 92, .requiredTimeScale = 2.0}};
            }

            return std::vector<path_tempo::PieceCorrection> {};
        };
        path_tempo::PathPlanner correctedPlanner;
        const auto corrected = correctedPlanner.solve(request, correction);

        require(corrected.has_value(), "a valid materialization correction should re-solve");
        require(callbackCalls == 2, "materialization correction should receive the corrected candidate");
        require(corrected->diagnostics.correctedPieces == 1, "materialization correction should retain corrected-piece diagnostics");
        require(corrected->pieceLimits.size() == 2, "materialization correction should report every piece limit");
        require(corrected->pieceLimits[0].velocity == baseline->pieceLimits[0].velocity, "materialization correction should not tighten the unrequested piece");
        require(corrected->pieceLimits[1].velocity <= baseline->pieceLimits[1].velocity / 2.0, "materialization correction should tighten the requested non-first piece");
        require(corrected->diagnostics.trajectoryDuration > baseline->diagnostics.trajectoryDuration, "a two-times local time scale should slow the corrected path");

        path_tempo::PathPlanner invalidPlanner;
        const auto invalid = invalidPlanner.solve(request, [](const auto &) -> std::expected<std::vector<path_tempo::PieceCorrection>, std::string> {
            return std::vector {path_tempo::PieceCorrection {.piece = 999, .requiredTimeScale = 2.0}};
        });

        require(!invalid.has_value(), "an unknown materialization correction piece should be rejected");
        require(invalid.error().code == path_tempo::PlanningErrorCode::InvalidInput, "an invalid materialization correction should return an input error");

        path_tempo::PathPlanner duplicatePlanner;
        const auto duplicate = duplicatePlanner.solve(request, [](const auto &) -> std::expected<std::vector<path_tempo::PieceCorrection>, std::string> {
            return std::vector {
                path_tempo::PieceCorrection {.piece = 92, .requiredTimeScale = 1.1},
                path_tempo::PieceCorrection {.piece = 92, .requiredTimeScale = 1.2},
            };
        });

        require(!duplicate.has_value(), "duplicate materialization corrections should be rejected");
        require(duplicate.error().code == path_tempo::PlanningErrorCode::InvalidInput, "a duplicate materialization correction should return an input error");

        auto onePassRequest = request;
        onePassRequest.settings.maximumCorrectionPasses = 1;
        path_tempo::PathPlanner exhaustedPlanner;
        const auto exhausted = exhaustedPlanner.solve(onePassRequest, [](const auto &) -> std::expected<std::vector<path_tempo::PieceCorrection>, std::string> {
            return std::vector {path_tempo::PieceCorrection {.piece = 92, .requiredTimeScale = 2.0}};
        });

        require(!exhausted.has_value(), "a materialization correction should fail when its pass budget is exhausted");
        require(exhausted.error().message.contains("materialization correction"), "correction exhaustion should identify materialization as its source");
    }
}
int main() {
    require(path_tempo::version() == PATH_TEMPO_EXPECTED_VERSION, "version should match the project version");
    testLineSampling();
    testTinyLineSampling();
    testLineSamplingErrors();
    testArcSampling();
    testBSplineSampling();
    testUnclampedNonUniformBSplineSampling();
    testBSplineKnotScaleInvariance();
    testNurbsSampling();
    testNurbsGeometricScaleInvariance();
    testNurbsWeightScaleInvariance();
    testParametricSamplingUsesSpeedOnlyEvaluator();
    testCubicTimeSegment();
    testRestToRestTransition();
    testMovingBoundaryTransition();
    testUnboundedAccelerationCruise();
    testUnboundedAccelerationCruiseValidation();
    testSubPicosecondCruisePhase();
    testInvalidTransition();
    testUnrecoverableTransitionHasDetailedError();
    testTransitionLimitValidation();
    testShortTransitionVelocityTolerance();
    testVelocityReachability();
    testPathInputValidation();
    testMultiPiecePathPlanning();
    testMultiPiecePathUsesNonzeroBoundaryAcceleration();
    testBoundaryAccelerationCandidatePreference();
    testMultiPiecePathRejectsTangentDiscontinuity();
    testCurvatureContinuityUsesGeometricScale();
    testCurvedPathPlanning();
    testCoupledStationsAtTransitionBoundaries();
    testAcceptedTransitionEndpointChecksCoupledStation();
    testShortJerkPhaseUsesBracketingDifferentialStations();
    testCurvedBoundaryRejectsGeometricJerkViolation();
    testTinyNonzeroGeometryConstrainsBoundaryVelocity();
    testMultiPieceCurvedPlanning();
    testCurvedMovingBoundaryStates();
    testMaterializationCorrectionCallback();
}
