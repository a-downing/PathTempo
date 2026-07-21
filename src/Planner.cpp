#include "path_tempo/Planner.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <numeric>
#include <utility>

#include <ruckig/ruckig.hpp>

namespace path_tempo {
    namespace {
        // Eight extra forward/backward passes let endpoint reductions propagate
        // beyond the theoretical two passes per piece before convergence testing.
        constexpr std::size_t REACHABILITY_EXTRA_PASSES = 8;
        // Fifty-two bisections span the significand precision of a double.
        constexpr unsigned VELOCITY_BISECTION_ITERATIONS = 52;
        // Four extra steps make station-time inversion insensitive to endpoint
        // rounding after the double-precision position interval has collapsed.
        constexpr unsigned STATION_TIME_BISECTION_ITERATIONS = 56;
        // Values at this scale are below useful double-precision resolution for
        // the planner's subtraction and polynomial-extremum calculations.
        constexpr double NUMERICAL_ZERO_TOLERANCE = 1e-15;
        // Preserve a small amount of piece length so a reachability estimate is
        // conservative when the resulting boundary is materialized by Ruckig.
        constexpr double REACHABILITY_ABSOLUTE_DISTANCE_MARGIN = 1e-12;
        constexpr double REACHABILITY_RELATIVE_DISTANCE_MARGIN = 1e-6;
        // Constraint comparisons allow relative roundoff from coupled polynomial
        // evaluation without accepting a materially over-limit state.
        constexpr double CONSTRAINT_RELATIVE_TOLERANCE = 1e-10;
        // Infinite-acceleration transitions still require numerically identical
        // boundary states; this absorbs only absolute representation noise.
        constexpr double BOUNDARY_STATE_ABSOLUTE_TOLERANCE = 1e-12;
        // Phases this short cannot materially affect a double-precision time law
        // and retaining them can create non-increasing phase boundaries.
        constexpr double MINIMUM_PHASE_DURATION = 1e-12;
        // Ruckig endpoint reconstruction accumulates error with distance and
        // speed, so validation uses both an absolute floor and relative term.
        constexpr double TRANSITION_ABSOLUTE_TOLERANCE = 1e-12;
        constexpr double TRANSITION_DISTANCE_RELATIVE_TOLERANCE = 1e-9;
        constexpr double TRANSITION_VELOCITY_RELATIVE_TOLERANCE = 1e-10;
        // Station lookup combines independently sampled geometry and time-law
        // evaluation, requiring a small absolute and length-relative allowance.
        constexpr double STATION_LOOKUP_ABSOLUTE_TOLERANCE = 1e-12;
        constexpr double STATION_LOOKUP_RELATIVE_TOLERANCE = 1e-10;
        // Reachability passes stop once further boundary-speed changes are below
        // meaningful solver resolution; fixed endpoints get a slightly looser test.
        constexpr double REACHABILITY_CONVERGENCE_TOLERANCE = 1e-11;
        constexpr double FIXED_ENDPOINT_REACHABILITY_TOLERANCE = 1e-10;
        // Time scales within one part per billion of unity are numerical noise.
        constexpr double CORRECTION_SIGNIFICANCE_TOLERANCE = 1e-9;
        // Add 1% when a sampled violation is real so the next solve lands inside
        // the feasible region instead of repeating at the same numerical boundary.
        constexpr double CORRECTION_SAFETY_FACTOR = 1.01;
    }

    double velocityTransitionDistance(const double fromVelocity, const double toVelocity, const double acceleration, const double jerk) {
        const auto change = std::abs(toVelocity - fromVelocity);

        if (change <= NUMERICAL_ZERO_TOLERANCE) {
            return 0.0;
        }

        const auto threshold = acceleration * acceleration / jerk;
        const auto duration = change <= threshold ? 2.0 * std::sqrt(change / jerk) : change / acceleration + acceleration / jerk;

        return std::midpoint(fromVelocity, toVelocity) * duration;
    }

    double reachableVelocity(const double fixedVelocity, const double cap, const double length, const double acceleration, const double jerk) {
        if (cap <= fixedVelocity) {
            return cap;
        }

        const auto available = std::max(0.0, length - std::max(REACHABILITY_ABSOLUTE_DISTANCE_MARGIN, length * REACHABILITY_RELATIVE_DISTANCE_MARGIN));

        if (velocityTransitionDistance(fixedVelocity, cap, acceleration, jerk) <= available) {
            return cap;
        }

        auto low = fixedVelocity;
        auto high = cap;

        for (unsigned iteration = 0; iteration < VELOCITY_BISECTION_ITERATIONS; ++iteration) {
            const auto middle = std::midpoint(low, high);

            if (velocityTransitionDistance(fixedVelocity, middle, acceleration, jerk) <= available) {
                low = middle;
            } else {
                high = middle;
            }
        }

        return low;
    }

    struct ScalarTransitionPlanner::Implementation {
        ruckig::InputParameter<1> input;
        ruckig::Ruckig<1> generator;
        ruckig::Trajectory<1> trajectory;
    };

    struct PathPlanner::Implementation {
        ScalarTransitionPlanner transitionPlanner;
    };

    std::size_t ScalarTransition::size() const noexcept {
        return m_size;
    }

    bool ScalarTransition::empty() const noexcept {
        return m_size == 0;
    }

    double ScalarTransition::duration() const noexcept {
        auto result = 0.0;

        for (const auto &segment : *this) {
            result += segment.duration;
        }

        return result;
    }

    const CubicTimeSegment &ScalarTransition::front() const {
        assert(!empty());

        return m_segments.front();
    }

    const CubicTimeSegment &ScalarTransition::back() const {
        assert(!empty());

        return m_segments[m_size - 1];
    }

    const CubicTimeSegment &ScalarTransition::operator[](const std::size_t index) const {
        assert(index < m_size);

        return m_segments[index];

    }

    bool ScalarTransition::push(const CubicTimeSegment &segment) noexcept {
        if (m_size == CAPACITY) {
            return false;
        }

        m_segments[m_size++] = segment;

        return true;
    }

    ScalarTransitionPlanner::ScalarTransitionPlanner() : m_implementation(std::make_unique<Implementation>()) {
    }

    ScalarTransitionPlanner::~ScalarTransitionPlanner() = default;

    ScalarTransitionPlanner::ScalarTransitionPlanner(ScalarTransitionPlanner &&) noexcept = default;

    ScalarTransitionPlanner &ScalarTransitionPlanner::operator=(ScalarTransitionPlanner &&) noexcept = default;

    std::expected<ScalarTransition, PlanningError> ScalarTransitionPlanner::solve(const ScalarTransitionRequest &request) {
        const auto &from = request.beginning;
        const auto &to = request.ending;

        if (!std::isfinite(request.length) || !std::isfinite(from.velocity) || !std::isfinite(from.acceleration) || !std::isfinite(to.velocity) || !std::isfinite(to.acceleration) || std::isnan(request.maximumVelocity) || std::isnan(request.maximumAcceleration) || std::isnan(request.maximumJerk) || request.length <= 0.0 || from.velocity < 0.0 || to.velocity < 0.0 || request.maximumVelocity <= 0.0 || request.maximumAcceleration <= 0.0 || request.maximumJerk <= 0.0 || std::abs(from.acceleration) > request.maximumAcceleration * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || std::abs(to.acceleration) > request.maximumAcceleration * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE)) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::InvalidInput,
                .message = "local time law received invalid distance, state, or limits",
            });
        }

        if (std::isinf(request.maximumAcceleration)) {
            if (std::abs(from.velocity - to.velocity) > BOUNDARY_STATE_ABSOLUTE_TOLERANCE || from.velocity <= 0.0 || std::abs(from.acceleration) > BOUNDARY_STATE_ABSOLUTE_TOLERANCE || std::abs(to.acceleration) > BOUNDARY_STATE_ABSOLUTE_TOLERANCE) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = "unbounded-acceleration local time law requires equal positive speeds and zero acceleration",
                });
            }

            ScalarTransition result;
            result.beginning = from;
            result.ending = to;
            result.push({
                .piece = request.piece,
                .duration = request.length / from.velocity,
                .c0 = 0.0,
                .c1 = from.velocity,
            });

            return result;
        }

        auto &input = m_implementation->input;
        input.current_position = {0.0};
        input.current_velocity = {from.velocity};
        input.current_acceleration = {from.acceleration};
        input.target_position = {request.length};
        input.target_velocity = {to.velocity};
        input.target_acceleration = {to.acceleration};
        input.max_velocity = {request.maximumVelocity};
        input.max_acceleration = {request.maximumAcceleration};
        input.max_jerk = {request.maximumJerk};

        auto &trajectory = m_implementation->trajectory;

        if (m_implementation->generator.calculate(input, trajectory) != ruckig::Result::Working) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::SolverFailure,
                .message = std::format("Ruckig failed local position timing: length={} state (v={}, a={}) -> (v={}, a={}) limits v={} a={} j={}", request.length, from.velocity, from.acceleration, to.velocity, to.acceleration, request.maximumVelocity, request.maximumAcceleration, request.maximumJerk),
            });
        }

        const auto profiles = trajectory.get_profiles();
        const auto &profile = profiles.front().front();
        ScalarTransition result;
        result.beginning = from;
        result.ending = to;
        auto distance = 0.0;
        auto velocity = from.velocity;
        auto acceleration = from.acceleration;

        const auto appendPhase = [&](const double duration, const double jerk, const bool initialStateCorrection = false) -> bool {
            if (duration <= MINIMUM_PHASE_DURATION) {
                return true;
            }

            if (!result.push({
                .piece = request.piece,
                .duration = duration,
                .c0 = distance,
                .c1 = velocity,
                .c2 = acceleration / 2.0,
                .c3 = jerk / 6.0,
                .initialStateCorrection = initialStateCorrection,
            })) {
                return false;
            }

            distance += velocity * duration + acceleration * duration * duration / 2.0 + jerk * duration * duration * duration / 6.0;
            velocity += acceleration * duration + jerk * duration * duration / 2.0;
            acceleration += jerk * duration;

            return true;
        };

        for (std::size_t phase = 0; phase < profile.brake.t.size(); ++phase) {
            if (!appendPhase(profile.brake.t[phase], profile.brake.j[phase], true)) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::CapacityExceeded,
                    .message = "local time law exceeded its fixed phase capacity",
                });
            }
        }

        for (std::size_t phase = 0; phase < profile.t.size(); ++phase) {
            if (!appendPhase(profile.t[phase], profile.j[phase])) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::CapacityExceeded,
                    .message = "local time law exceeded its fixed phase capacity",
                });
            }
        }

        if (!appendPhase(trajectory.get_duration() - result.duration(), 0.0)) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::CapacityExceeded,
                .message = "local time law exceeded its fixed phase capacity",
            });
        }

        const auto distanceTolerance = std::max(TRANSITION_ABSOLUTE_TOLERANCE, request.length * TRANSITION_DISTANCE_RELATIVE_TOLERANCE);
        const auto velocityTolerance = std::max(TRANSITION_ABSOLUTE_TOLERANCE, request.maximumVelocity * TRANSITION_VELOCITY_RELATIVE_TOLERANCE);
        auto previousTime = 0.0;
        auto previousDistance = 0.0;

        for (std::size_t phase = 0; phase < result.size(); ++phase) {
            const auto &segment = result[phase];
            const auto time = previousTime + segment.duration;
            const auto phaseDistance = phase + 1 == result.size() ? request.length : segment.position(segment.duration);
            const auto phaseVelocity = phase + 1 == result.size() ? to.velocity : segment.velocity(segment.duration);
            const auto phaseAcceleration = phase + 1 == result.size() ? to.acceleration : segment.acceleration(segment.duration);

            if (!std::isfinite(time) || !std::isfinite(phaseDistance) || !std::isfinite(phaseVelocity) || !std::isfinite(phaseAcceleration) || phaseDistance < -distanceTolerance || phaseDistance > request.length + distanceTolerance || phaseVelocity < -distanceTolerance || time <= previousTime || phaseDistance < previousDistance - distanceTolerance) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::NonMonotoneResult,
                    .message = std::format("Ruckig produced a non-monotone local path law at boundary {} of {}: time={} distance={} velocity={} acceleration={} previous time={} previous distance={} length={} tolerance={}", phase + 1, result.size() + 1, time, phaseDistance, phaseVelocity, phaseAcceleration, previousTime, previousDistance, request.length, distanceTolerance),
                });
            }

            auto minimumVelocity = std::min(segment.c1, phaseVelocity);

            if (std::abs(segment.jerk()) > NUMERICAL_ZERO_TOLERANCE) {
                const auto extremum = -segment.acceleration(0.0) / segment.jerk();

                if (extremum > 0.0 && extremum < segment.duration) {
                    minimumVelocity = std::min(minimumVelocity, segment.velocity(extremum));
                }
            }

            if (minimumVelocity < -velocityTolerance) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::DirectionReversal,
                    .message = std::format("Ruckig local position timing reverses path direction in phase {}: minimum_velocity={} tolerance={} duration={} length={} state (v={}, a={}) -> (v={}, a={}) limits v={} a={} j={}", phase, minimumVelocity, velocityTolerance, segment.duration, request.length, from.velocity, from.acceleration, to.velocity, to.acceleration, request.maximumVelocity, request.maximumAcceleration, request.maximumJerk),
                });
            }

            previousTime = time;
            previousDistance = phaseDistance;
        }

        return result;
    }

    PathPlanner::PathPlanner() : m_implementation(std::make_unique<Implementation>()) {
    }

    PathPlanner::~PathPlanner() = default;

    PathPlanner::PathPlanner(PathPlanner &&) noexcept = default;

    PathPlanner &PathPlanner::operator=(PathPlanner &&) noexcept = default;

    std::expected<PlannedPath, PlanningError> PathPlanner::solveLocal(const std::span<const LocalPiece> pieces, const BoundaryState beginning, const BoundaryState ending, const CoupledLimits &limits, const PathPlanningSettings &settings, const MaterializationCorrection &materializationCorrection) {
        auto correctedPieces = std::vector<bool>(pieces.size(), false);
        auto localPieces = std::vector<LocalPiece>(pieces.begin(), pieces.end());
        const auto stationCount = pieces.size() + 1;

        const auto coupledTimeScale = [&](const LocalPiece &piece, const ScalarTransition &transition) {
            auto required = 1.0;
            const auto consider = [&](const double measured, const double limit, const double exponent) {
                if (std::isfinite(limit) && measured > limit) {
                    required = std::max(required, std::pow(measured / limit, exponent));
                }
            };

            for (const auto &station : piece.stations) {
                const auto tolerance = std::max(STATION_LOOKUP_ABSOLUTE_TOLERANCE, piece.length * STATION_LOOKUP_RELATIVE_TOLERANCE);

                for (const auto &segment : transition) {
                    const auto from = segment.position(0.0);
                    const auto to = segment.position(segment.duration);

                    if (station.distance < from - tolerance || station.distance > to + tolerance) {
                        continue;
                    }

                    auto low = 0.0;
                    auto high = segment.duration;

                    for (unsigned iteration = 0; iteration < STATION_TIME_BISECTION_ITERATIONS; ++iteration) {
                        const auto middle = std::midpoint(low, high);

                        if (segment.position(middle) < station.distance) {
                            low = middle;
                        } else {
                            high = middle;
                        }
                    }

                    const auto time = std::midpoint(low, high);
                    const auto velocity = std::max(0.0, segment.velocity(time));
                    const auto acceleration = segment.acceleration(time);
                    const auto jerk = segment.jerk();
                    auto accelerationSquared = 0.0;
                    auto jerkSquared = 0.0;

                    for (std::size_t axis = 0; axis < station.tangent.size(); ++axis) {
                        const auto coordinateVelocity = station.tangent[axis] * velocity;
                        const auto coordinateAcceleration = station.tangent[axis] * acceleration + station.curvature[axis] * velocity * velocity;
                        const auto coordinateJerk = station.tangent[axis] * jerk + 3.0 * station.curvature[axis] * velocity * acceleration + station.thirdDerivative[axis] * velocity * velocity * velocity;
                        accelerationSquared += coordinateAcceleration * coordinateAcceleration;
                        jerkSquared += coordinateJerk * coordinateJerk;
                        consider(std::abs(coordinateVelocity), limits.coordinateVelocity[axis], 1.0);
                        consider(std::abs(coordinateAcceleration), limits.coordinateAcceleration[axis], 0.5);
                        consider(std::abs(coordinateJerk), limits.coordinateJerk[axis], 1.0 / 3.0);
                    }

                    consider(std::sqrt(accelerationSquared), limits.pathAcceleration, 0.5);
                    consider(std::sqrt(jerkSquared), limits.pathJerk, 1.0 / 3.0);
                }
            }

            return required;
        };

        for (std::size_t correctionPass = 0; correctionPass < settings.maximumCorrectionPasses; ++correctionPass) {
            std::vector<double> stationCaps(stationCount, std::numeric_limits<double>::infinity());
            stationCaps.front() = beginning.velocity;
            stationCaps.back() = ending.velocity;

            for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                stationCaps[station] = std::min(localPieces[station - 1].maximumVelocity, localPieces[station].maximumVelocity);
            }

            if (beginning.velocity > localPieces.front().maximumVelocity * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || ending.velocity > localPieces.back().maximumVelocity * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || std::abs(beginning.acceleration) > localPieces.front().maximumAcceleration * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || std::abs(ending.acceleration) > localPieces.back().maximumAcceleration * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE)) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = correctionPass == 0 ? "path boundary state exceeds a local piece limit" : "fixed path boundary state cannot satisfy the coupled curved-path limits",
                });
            }

            auto stationVelocity = stationCaps;
            const auto reachabilityPasses = 2 * pieces.size() + REACHABILITY_EXTRA_PASSES;

            for (std::size_t pass = 0; pass < reachabilityPasses; ++pass) {
                auto maximumChange = 0.0;

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    const auto &piece = localPieces[pieceIndex];
                    const auto reachable = reachableVelocity(stationVelocity[pieceIndex], std::min(stationCaps[pieceIndex + 1], piece.maximumVelocity), piece.length, piece.maximumAcceleration, piece.maximumJerk);

                    if (pieceIndex + 1 == pieces.size()) {
                        if (reachable + FIXED_ENDPOINT_REACHABILITY_TOLERANCE < ending.velocity) {
                            return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = "fixed ending velocity is not forward reachable"});
                        }
                    } else {
                        const auto reduced = std::min(stationVelocity[pieceIndex + 1], reachable);
                        maximumChange = std::max(maximumChange, stationVelocity[pieceIndex + 1] - reduced);
                        stationVelocity[pieceIndex + 1] = reduced;
                    }
                }

                for (std::size_t pieceIndex = pieces.size(); pieceIndex-- > 0;) {
                    const auto &piece = localPieces[pieceIndex];
                    const auto reachable = reachableVelocity(stationVelocity[pieceIndex + 1], std::min(stationCaps[pieceIndex], piece.maximumVelocity), piece.length, piece.maximumAcceleration, piece.maximumJerk);

                    if (pieceIndex == 0) {
                        if (reachable + FIXED_ENDPOINT_REACHABILITY_TOLERANCE < beginning.velocity) {
                            return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = "fixed beginning velocity cannot reach the remaining path"});
                        }
                    } else {
                        const auto reduced = std::min(stationVelocity[pieceIndex], reachable);
                        maximumChange = std::max(maximumChange, stationVelocity[pieceIndex] - reduced);
                        stationVelocity[pieceIndex] = reduced;
                    }
                }

                if (maximumChange <= REACHABILITY_CONVERGENCE_TOLERANCE) {
                    break;
                }
            }

            PlannedPath result;
            result.timeLaw.beginning = beginning;
            result.timeLaw.ending = ending;
            result.pieceBoundaries.resize(stationCount);
            result.pieceLimits.reserve(pieces.size());
            result.pieceBoundaries.front() = beginning;
            result.pieceBoundaries.back() = ending;
            result.diagnostics.correctionPasses = correctionPass + 1;
            std::vector<ScalarTransition> transitions;
            transitions.reserve(pieces.size());

            for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                result.pieceBoundaries[station] = {.velocity = stationVelocity[station]};
            }

            for (const auto &piece : localPieces) {
                result.pieceLimits.push_back({
                    .velocity = piece.maximumVelocity,
                    .acceleration = piece.maximumAcceleration,
                    .jerk = piece.maximumJerk,
                });
            }

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                const auto &piece = localPieces[pieceIndex];
                const auto transition = m_implementation->transitionPlanner.solve({
                    .piece = piece.id,
                    .length = piece.length,
                    .beginning = result.pieceBoundaries[pieceIndex],
                    .ending = result.pieceBoundaries[pieceIndex + 1],
                    .maximumVelocity = piece.maximumVelocity,
                    .maximumAcceleration = piece.maximumAcceleration,
                    .maximumJerk = piece.maximumJerk,
                });

                if (!transition) {
                    return std::unexpected(PlanningError {
                        .code = transition.error().code,
                        .message = std::format("could not materialize path piece {}: {}", pieceIndex, transition.error().message),
                    });
                }

                transitions.push_back(*transition);
            }

            std::vector<double> corrections(pieces.size(), 1.0);
            auto maximumCorrection = 1.0;

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                const auto sampledCorrection = settings.applySampledCorrections ? coupledTimeScale(localPieces[pieceIndex], transitions[pieceIndex]) : 1.0;
                corrections[pieceIndex] = sampledCorrection > 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE ? sampledCorrection * CORRECTION_SAFETY_FACTOR : 1.0;
                maximumCorrection = std::max(maximumCorrection, corrections[pieceIndex]);
            }

            auto refinedPathDistance = 0.0;

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                result.diagnostics.velocitySeedDuration += transitions[pieceIndex].duration();

                for (auto segment : transitions[pieceIndex]) {
                    segment.c0 += refinedPathDistance;
                    result.timeLaw.segments.push_back(segment);
                }

                refinedPathDistance += pieces[pieceIndex].length;
            }

            if (materializationCorrection) {
                const auto requested = materializationCorrection(result);

                if (!requested) {
                    return std::unexpected(PlanningError {
                        .code = PlanningErrorCode::SolverFailure,
                        .message = std::format("materialization correction failed: {}", requested.error()),
                    });
                }

                std::vector<bool> supplied(pieces.size(), false);

                for (const auto &correction : *requested) {
                    const auto found = std::ranges::find(localPieces, correction.piece, &LocalPiece::id);

                    if (found == localPieces.end() || !std::isfinite(correction.requiredTimeScale) || correction.requiredTimeScale < 1.0) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("materialization returned an invalid correction for piece {} with time scale {}", correction.piece, correction.requiredTimeScale),
                        });
                    }

                    const auto pieceIndex = static_cast<std::size_t>(found - localPieces.begin());

                    if (supplied[pieceIndex]) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("materialization returned duplicate corrections for piece {}", correction.piece),
                        });
                    }

                    supplied[pieceIndex] = true;
                    corrections[pieceIndex] = std::max(corrections[pieceIndex], correction.requiredTimeScale);
                    maximumCorrection = std::max(maximumCorrection, corrections[pieceIndex]);
                }
            }

            if (maximumCorrection <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE) {
                result.diagnostics.correctedPieces = std::ranges::count(correctedPieces, true);

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    result.diagnostics.maximumAppliedTimeScale = std::max(result.diagnostics.maximumAppliedTimeScale, pieces[pieceIndex].maximumVelocity / localPieces[pieceIndex].maximumVelocity);
                }

                return result;
            }

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                if (corrections[pieceIndex] <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE) {
                    continue;
                }

                correctedPieces[pieceIndex] = true;
                const auto factor = corrections[pieceIndex];
                localPieces[pieceIndex].maximumVelocity /= factor;
                localPieces[pieceIndex].maximumAcceleration /= factor * factor;
                localPieces[pieceIndex].maximumJerk /= factor * factor * factor;
            }
        }

        return std::unexpected(PlanningError {
            .code = PlanningErrorCode::SolverFailure,
            .message = std::format("coupled curved-path constraint correction did not converge after {} passes", settings.maximumCorrectionPasses),
        });
    }
}
