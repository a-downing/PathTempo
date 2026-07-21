#include "path_tempo/Planner.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <numeric>
#include <utility>

#include <ruckig/ruckig.hpp>

#include "path_tempo/LinearOptimization.h"

namespace path_tempo {
    double velocityTransitionDistance(const double fromVelocity, const double toVelocity, const double acceleration, const double jerk) {
        const auto change = std::abs(toVelocity - fromVelocity);

        if (change <= 1e-15) {
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

        const auto available = std::max(0.0, length - std::max(1e-12, length * 1e-6));

        if (velocityTransitionDistance(fixedVelocity, cap, acceleration, jerk) <= available) {
            return cap;
        }

        auto low = fixedVelocity;
        auto high = cap;

        for (unsigned iteration = 0; iteration < 52; ++iteration) {
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
        PersistentLinearSolver refinementSolver;
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

        if (!std::isfinite(request.length) || !std::isfinite(from.velocity) || !std::isfinite(from.acceleration) || !std::isfinite(to.velocity) || !std::isfinite(to.acceleration) || std::isnan(request.maximumVelocity) || std::isnan(request.maximumAcceleration) || std::isnan(request.maximumJerk) || request.length <= 0.0 || from.velocity < 0.0 || to.velocity < 0.0 || request.maximumVelocity <= 0.0 || request.maximumAcceleration <= 0.0 || request.maximumJerk <= 0.0 || std::abs(from.acceleration) > request.maximumAcceleration * (1.0 + 1e-10) || std::abs(to.acceleration) > request.maximumAcceleration * (1.0 + 1e-10)) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::InvalidInput,
                .message = "local time law received invalid distance, state, or limits",
            });
        }

        if (std::isinf(request.maximumAcceleration)) {
            if (std::abs(from.velocity - to.velocity) > 1e-12 || from.velocity <= 0.0 || std::abs(from.acceleration) > 1e-12 || std::abs(to.acceleration) > 1e-12) {
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
            if (duration <= 1e-12) {
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

        const auto distanceTolerance = std::max(1e-12, request.length * 1e-9);
        const auto velocityTolerance = std::max(1e-12, request.maximumVelocity * 1e-10);
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

            if (std::abs(segment.jerk()) > 1e-15) {
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
        if (auto configured = m_implementation->refinementSolver.configure(settings.linearSolveTimeLimit); !configured) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::SolverFailure,
                .message = std::format("could not configure the coupled refinement solver: {}", configured.error()),
            });
        }
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
                const auto tolerance = std::max(1e-12, piece.length * 1e-10);

                for (const auto &segment : transition) {
                    const auto from = segment.position(0.0);
                    const auto to = segment.position(segment.duration);

                    if (station.distance < from - tolerance || station.distance > to + tolerance) {
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
                    const auto velocity = std::max(0.0, segment.velocity(time));
                    const auto acceleration = segment.acceleration(time);
                    const auto jerk = segment.jerk();
                    auto accelerationSquared = 0.0;
                    auto jerkSquared = 0.0;

                    for (std::size_t axis = 0; axis < station.tangent.size(); ++axis) {
                        const auto axisVelocity = station.tangent[axis] * velocity;
                        const auto axisAcceleration = station.tangent[axis] * acceleration + station.curvature[axis] * velocity * velocity;
                        const auto axisJerk = station.tangent[axis] * jerk + 3.0 * station.curvature[axis] * velocity * acceleration + station.thirdDerivative[axis] * velocity * velocity * velocity;
                        accelerationSquared += axisAcceleration * axisAcceleration;
                        jerkSquared += axisJerk * axisJerk;
                        consider(std::abs(axisVelocity), limits.axisVelocity[axis], 1.0);
                        consider(std::abs(axisAcceleration), limits.axisAcceleration[axis], 0.5);
                        consider(std::abs(axisJerk), limits.axisJerk[axis], 1.0 / 3.0);
                    }

                    consider(std::sqrt(accelerationSquared), limits.pathAcceleration, 0.5);
                    consider(std::sqrt(jerkSquared), limits.pathJerk, 1.0 / 3.0);
                }
            }

            return required;
        };

        const auto feasibleEndpointJerk = [&](const LocalPiece &piece, const LocalPiece::Station &geometry, const double velocity, const double acceleration) -> std::optional<double> {
            std::vector<double> geometricJerk(geometry.tangent.size());
            auto lower = -piece.maximumJerk;
            auto upper = piece.maximumJerk;

            for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                geometricJerk[axis] = 3.0 * geometry.curvature[axis] * velocity * acceleration + geometry.thirdDerivative[axis] * velocity * velocity * velocity;
                const auto tangent = geometry.tangent[axis];
                const auto limit = limits.axisJerk[axis];

                if (!std::isfinite(limit)) {
                    continue;
                }

                if (std::abs(tangent) <= 1e-15) {
                    if (std::abs(geometricJerk[axis]) > limit * (1.0 + 1e-10)) {
                        return std::nullopt;
                    }

                    continue;
                }

                auto componentLower = (-limit - geometricJerk[axis]) / tangent;
                auto componentUpper = (limit - geometricJerk[axis]) / tangent;

                if (componentLower > componentUpper) {
                    std::swap(componentLower, componentUpper);
                }

                lower = std::max(lower, componentLower);
                upper = std::min(upper, componentUpper);
            }

            if (lower > upper) {
                return std::nullopt;
            }

            auto projection = 0.0;

            for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                projection += geometry.tangent[axis] * geometricJerk[axis];
            }

            const auto scalarJerk = std::clamp(-projection, lower, upper);
            auto coupledSquared = 0.0;

            for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                const auto coupled = geometry.tangent[axis] * scalarJerk + geometricJerk[axis];
                coupledSquared += coupled * coupled;
            }

            if (std::isfinite(limits.pathJerk) && std::sqrt(coupledSquared) > limits.pathJerk * (1.0 + 1e-10)) {
                return std::nullopt;
            }

            return scalarJerk;
        };

        const auto endpointFeasible = [&](const LocalPiece &piece, const LocalPiece::Station &geometry, const double velocity, const double acceleration) {
            auto accelerationSquared = 0.0;

            for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                const auto coupled = geometry.tangent[axis] * acceleration + geometry.curvature[axis] * velocity * velocity;

                if (std::isfinite(limits.axisAcceleration[axis]) && std::abs(coupled) > limits.axisAcceleration[axis] * (1.0 + 1e-10)) {
                    return false;
                }

                accelerationSquared += coupled * coupled;
            }

            if (std::isfinite(limits.pathAcceleration) && std::sqrt(accelerationSquared) > limits.pathAcceleration * (1.0 + 1e-10)) {
                return false;
            }

            return feasibleEndpointJerk(piece, geometry, velocity, acceleration).has_value();
        };

        for (std::size_t correctionPass = 0; correctionPass < settings.maximumCorrectionPasses; ++correctionPass) {
            std::vector<double> stationCaps(stationCount, std::numeric_limits<double>::infinity());
            stationCaps.front() = beginning.velocity;
            stationCaps.back() = ending.velocity;

            for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                stationCaps[station] = std::min(localPieces[station - 1].maximumVelocity, localPieces[station].maximumVelocity);
            }

            if (beginning.velocity > localPieces.front().maximumVelocity * (1.0 + 1e-10) || ending.velocity > localPieces.back().maximumVelocity * (1.0 + 1e-10) || std::abs(beginning.acceleration) > localPieces.front().maximumAcceleration * (1.0 + 1e-10) || std::abs(ending.acceleration) > localPieces.back().maximumAcceleration * (1.0 + 1e-10)) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = correctionPass == 0 ? "path boundary state exceeds a local piece limit" : "fixed path boundary state cannot satisfy the coupled curved-path limits",
                });
            }

            auto stationVelocity = stationCaps;
            const auto reachabilityPasses = 2 * pieces.size() + 8;

            for (std::size_t pass = 0; pass < reachabilityPasses; ++pass) {
                auto maximumChange = 0.0;

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    const auto &piece = localPieces[pieceIndex];
                    const auto reachable = reachableVelocity(stationVelocity[pieceIndex], std::min(stationCaps[pieceIndex + 1], piece.maximumVelocity), piece.length, piece.maximumAcceleration, piece.maximumJerk);

                    if (pieceIndex + 1 == pieces.size()) {
                        if (reachable + 1e-10 < ending.velocity) {
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
                        if (reachable + 1e-10 < beginning.velocity) {
                            return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = "fixed beginning velocity cannot reach the remaining path"});
                        }
                    } else {
                        const auto reduced = std::min(stationVelocity[pieceIndex], reachable);
                        maximumChange = std::max(maximumChange, stationVelocity[pieceIndex] - reduced);
                        stationVelocity[pieceIndex] = reduced;
                    }
                }

                if (maximumChange <= 1e-11) {
                    break;
                }
            }

            PlannedPath result;
            result.timeLaw.beginning = beginning;
            result.timeLaw.ending = ending;
            result.pieceBoundaries.resize(stationCount);
            result.pieceBoundaries.front() = beginning;
            result.pieceBoundaries.back() = ending;
            result.diagnostics.correctionPasses = correctionPass + 1;
            std::vector<ScalarTransition> transitions;
            transitions.reserve(pieces.size());
            auto pathDistance = 0.0;

            for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                result.pieceBoundaries[station] = {.velocity = stationVelocity[station]};
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

                result.diagnostics.velocitySeedDuration += transition->duration();
                transitions.push_back(*transition);

                for (auto segment : *transition) {
                    segment.c0 += pathDistance;
                    result.timeLaw.segments.push_back(segment);
                }

                pathDistance += piece.length;
            }

            for (std::size_t sequentialIteration = 0; sequentialIteration < settings.sequentialIterations && stationCount > 2; ++sequentialIteration) {
                const auto velocityColumn = [](const std::size_t station) { return station; };
                const auto accelerationColumn = [stationCount](const std::size_t station) { return stationCount + station; };
                const auto jerkColumn = [stationCount](const std::size_t pieceIndex, const bool end) { return 2 * stationCount + 2 * pieceIndex + (end ? 1U : 0U); };
                const auto deviationColumn = [stationCount, pieceCount = pieces.size()](const std::size_t station) { return 2 * stationCount + 2 * pieceCount + station; };
                const auto variableCount = 3 * stationCount + 2 * pieces.size();
                SparseLinearProgram refinement(variableCount);
                std::vector<double> accelerationTargets(stationCount);

                for (std::size_t station = 0; station < stationCount; ++station) {
                    const auto referenceVelocity = result.pieceBoundaries[station].velocity;
                    const auto referenceAcceleration = result.pieceBoundaries[station].acceleration;
                    auto accelerationLimit = std::numeric_limits<double>::infinity();

                    if (station > 0) {
                        accelerationLimit = std::min(accelerationLimit, localPieces[station - 1].maximumAcceleration);
                    }

                    if (station < pieces.size()) {
                        accelerationLimit = std::min(accelerationLimit, localPieces[station].maximumAcceleration);
                    }

                    accelerationLimit *= 0.95;
                    auto velocityLower = std::max(0.0, referenceVelocity - settings.velocityTrustFraction * std::max(1e-6, stationCaps[station]));
                    auto velocityUpper = std::min(stationCaps[station], referenceVelocity + settings.velocityTrustFraction * std::max(1e-6, stationCaps[station]));
                    auto accelerationLower = std::max(-accelerationLimit, referenceAcceleration - settings.accelerationTrustFraction * accelerationLimit);
                    auto accelerationUpper = std::min(accelerationLimit, referenceAcceleration + settings.accelerationTrustFraction * accelerationLimit);

                    if (station == 0 || station + 1 == stationCount) {
                        velocityLower = velocityUpper = referenceVelocity;
                        accelerationLower = accelerationUpper = referenceAcceleration;
                    }

                    refinement.columnLower(velocityColumn(station)) = velocityLower;
                    refinement.columnUpper(velocityColumn(station)) = velocityUpper;
                    refinement.columnLower(accelerationColumn(station)) = accelerationLower;
                    refinement.columnUpper(accelerationColumn(station)) = accelerationUpper;
                    refinement.columnLower(deviationColumn(station)) = 0.0;
                    refinement.columnUpper(deviationColumn(station)) = linearProgramInfinity();
                    refinement.columnCost(deviationColumn(station)) = 1e-6;
                    auto targetAcceleration = referenceAcceleration;

                    if (station > 0 && station < pieces.size()) {
                        const auto leftSlope = (referenceVelocity * referenceVelocity - result.pieceBoundaries[station - 1].velocity * result.pieceBoundaries[station - 1].velocity) / (2.0 * pieces[station - 1].length);
                        const auto rightSlope = (result.pieceBoundaries[station + 1].velocity * result.pieceBoundaries[station + 1].velocity - referenceVelocity * referenceVelocity) / (2.0 * pieces[station].length);
                        targetAcceleration = leftSlope * rightSlope > 0.0 ? std::copysign(std::min(std::abs(leftSlope), std::abs(rightSlope)), leftSlope) : 0.0;
                        targetAcceleration = std::clamp(targetAcceleration, accelerationLower, accelerationUpper);
                    }

                    accelerationTargets[station] = targetAcceleration;
                    refinement.addRow(-targetAcceleration, linearProgramInfinity(), {{accelerationColumn(station), -1.0}, {deviationColumn(station), 1.0}});
                    refinement.addRow(targetAcceleration, linearProgramInfinity(), {{accelerationColumn(station), 1.0}, {deviationColumn(station), 1.0}});
                }

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    refinement.columnLower(jerkColumn(pieceIndex, false)) = -localPieces[pieceIndex].maximumJerk;
                    refinement.columnUpper(jerkColumn(pieceIndex, false)) = localPieces[pieceIndex].maximumJerk;
                    refinement.columnLower(jerkColumn(pieceIndex, true)) = -localPieces[pieceIndex].maximumJerk;
                    refinement.columnUpper(jerkColumn(pieceIndex, true)) = localPieces[pieceIndex].maximumJerk;
                    const auto speedSum = std::max(1e-6, result.pieceBoundaries[pieceIndex].velocity + result.pieceBoundaries[pieceIndex + 1].velocity);
                    const auto objectiveDerivative = -2.0 * pieces[pieceIndex].length / (speedSum * speedSum);
                    refinement.columnCost(velocityColumn(pieceIndex)) += objectiveDerivative;
                    refinement.columnCost(velocityColumn(pieceIndex + 1)) += objectiveDerivative;
                }

                const auto addAccelerationConstraints = [&](const LocalPiece::Station &geometry, const std::size_t station) {
                    const auto referenceVelocity = result.pieceBoundaries[station].velocity;
                    const auto referenceAcceleration = result.pieceBoundaries[station].acceleration;
                    std::vector<double> referenceVector(geometry.tangent.size());
                    auto referenceNormSquared = 0.0;

                    for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                        const auto velocityCoefficient = 2.0 * geometry.curvature[axis] * referenceVelocity;
                        const auto constant = -geometry.curvature[axis] * referenceVelocity * referenceVelocity;

                        if (std::isfinite(limits.axisAcceleration[axis])) {
                            refinement.addRow(-limits.axisAcceleration[axis] - constant, limits.axisAcceleration[axis] - constant, {{velocityColumn(station), velocityCoefficient}, {accelerationColumn(station), geometry.tangent[axis]}});
                        }

                        referenceVector[axis] = geometry.tangent[axis] * referenceAcceleration + geometry.curvature[axis] * referenceVelocity * referenceVelocity;
                        referenceNormSquared += referenceVector[axis] * referenceVector[axis];
                    }

                    const auto referenceNorm = std::sqrt(referenceNormSquared);

                    if (std::isfinite(limits.pathAcceleration) && referenceNorm > 1e-12) {
                        auto velocityCoefficient = 0.0;
                        auto accelerationCoefficient = 0.0;
                        auto constant = 0.0;

                        for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                            const auto direction = referenceVector[axis] / referenceNorm;
                            velocityCoefficient += direction * 2.0 * geometry.curvature[axis] * referenceVelocity;
                            accelerationCoefficient += direction * geometry.tangent[axis];
                            constant -= direction * geometry.curvature[axis] * referenceVelocity * referenceVelocity;
                        }

                        refinement.addRow(-linearProgramInfinity(), limits.pathAcceleration - constant, {{velocityColumn(station), velocityCoefficient}, {accelerationColumn(station), accelerationCoefficient}});
                    }
                };

                const auto addJerkConstraints = [&](const LocalPiece::Station &geometry, const std::size_t station, const std::size_t pieceIndex, const bool end, const double referenceJerk) {
                    const auto referenceVelocity = result.pieceBoundaries[station].velocity;
                    const auto referenceAcceleration = result.pieceBoundaries[station].acceleration;
                    std::vector<double> velocityDerivative(geometry.tangent.size());
                    std::vector<double> accelerationDerivative(geometry.tangent.size());
                    std::vector<double> constant(geometry.tangent.size());
                    std::vector<double> referenceVector(geometry.tangent.size());
                    auto referenceNormSquared = 0.0;

                    for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                        velocityDerivative[axis] = 3.0 * geometry.curvature[axis] * referenceAcceleration + 3.0 * geometry.thirdDerivative[axis] * referenceVelocity * referenceVelocity;
                        accelerationDerivative[axis] = 3.0 * geometry.curvature[axis] * referenceVelocity;
                        const auto geometricReference = 3.0 * geometry.curvature[axis] * referenceVelocity * referenceAcceleration + geometry.thirdDerivative[axis] * referenceVelocity * referenceVelocity * referenceVelocity;
                        constant[axis] = geometricReference - velocityDerivative[axis] * referenceVelocity - accelerationDerivative[axis] * referenceAcceleration;
                        referenceVector[axis] = geometry.tangent[axis] * referenceJerk + geometricReference;
                        referenceNormSquared += referenceVector[axis] * referenceVector[axis];

                        if (std::isfinite(limits.axisJerk[axis])) {
                            refinement.addRow(-limits.axisJerk[axis] - constant[axis], limits.axisJerk[axis] - constant[axis], {{velocityColumn(station), velocityDerivative[axis]}, {accelerationColumn(station), accelerationDerivative[axis]}, {jerkColumn(pieceIndex, end), geometry.tangent[axis]}});
                        }
                    }

                    const auto referenceNorm = std::sqrt(referenceNormSquared);

                    if (std::isfinite(limits.pathJerk) && referenceNorm > 1e-12) {
                        auto velocityCoefficient = 0.0;
                        auto accelerationCoefficient = 0.0;
                        auto jerkCoefficient = 0.0;
                        auto projectedConstant = 0.0;

                        for (std::size_t axis = 0; axis < geometry.tangent.size(); ++axis) {
                            const auto direction = referenceVector[axis] / referenceNorm;
                            velocityCoefficient += direction * velocityDerivative[axis];
                            accelerationCoefficient += direction * accelerationDerivative[axis];
                            jerkCoefficient += direction * geometry.tangent[axis];
                            projectedConstant += direction * constant[axis];
                        }

                        refinement.addRow(-linearProgramInfinity(), limits.pathJerk - projectedConstant, {{velocityColumn(station), velocityCoefficient}, {accelerationColumn(station), accelerationCoefficient}, {jerkColumn(pieceIndex, end), jerkCoefficient}});
                    }
                };

                std::vector<double> referenceStartJerk(pieces.size());
                std::vector<double> referenceEndJerk(pieces.size());

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    const auto &piece = localPieces[pieceIndex];
                    const auto startJerk = feasibleEndpointJerk(piece, piece.stations.front(), result.pieceBoundaries[pieceIndex].velocity, result.pieceBoundaries[pieceIndex].acceleration);
                    const auto endJerk = feasibleEndpointJerk(piece, piece.stations.back(), result.pieceBoundaries[pieceIndex + 1].velocity, result.pieceBoundaries[pieceIndex + 1].acceleration);

                    if (!startJerk || !endJerk) {
                        return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = "coupled sequential reference has no feasible endpoint jerk"});
                    }

                    referenceStartJerk[pieceIndex] = *startJerk;
                    referenceEndJerk[pieceIndex] = *endJerk;
                    addAccelerationConstraints(piece.stations.front(), pieceIndex);
                    addAccelerationConstraints(piece.stations.back(), pieceIndex + 1);
                    addJerkConstraints(piece.stations.front(), pieceIndex, pieceIndex, false, *startJerk);
                    addJerkConstraints(piece.stations.back(), pieceIndex + 1, pieceIndex, true, *endJerk);
                }

                std::vector<double> referenceValues(variableCount);

                for (std::size_t station = 0; station < stationCount; ++station) {
                    referenceValues[velocityColumn(station)] = result.pieceBoundaries[station].velocity;
                    referenceValues[accelerationColumn(station)] = result.pieceBoundaries[station].acceleration;
                    referenceValues[deviationColumn(station)] = std::abs(result.pieceBoundaries[station].acceleration - accelerationTargets[station]);
                }

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    referenceValues[jerkColumn(pieceIndex, false)] = referenceStartJerk[pieceIndex];
                    referenceValues[jerkColumn(pieceIndex, true)] = referenceEndJerk[pieceIndex];
                }

                for (std::size_t row = 0; row < refinement.rowCount(); ++row) {
                    auto value = 0.0;

                    for (auto entry = refinement.rowBegin(row); entry < refinement.rowEnd(row); ++entry) {
                        value += refinement.entryValue(entry) * referenceValues[refinement.entryColumn(entry)];
                    }

                    const auto violation = std::max({0.0, refinement.rowLower(row) - value, value - refinement.rowUpper(row)});

                    if (violation > 1e-7) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::SolverFailure,
                            .message = std::format("coupled sequential linearization excludes its reference at row {}: violation={} value={} bounds=[{},{}]", row, violation, value, refinement.rowLower(row), refinement.rowUpper(row)),
                        });
                    }
                }

                const auto refined = m_implementation->refinementSolver.solve(refinement, {.simplexIterationLimit = settings.simplexIterationLimit});

                if (!refined) {
                    return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = std::format("coupled sequential refinement failed: {}", refined.error())});
                }

                ++result.diagnostics.sequentialSolves;
                result.diagnostics.linearSolverIterations += refined->diagnostics.simplexIterations;
                result.diagnostics.linearSolverBasisReused = result.diagnostics.linearSolverBasisReused || refined->diagnostics.basisReuseApplied;

                if (refined->status != LinearSolveStatus::Optimal) {
                    break;
                }

                auto acceptedAny = false;

                for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                    const auto proposedVelocity = std::clamp(refined->values[velocityColumn(station)], 0.0, stationCaps[station]);
                    const auto proposedAcceleration = refined->values[accelerationColumn(station)];
                    const auto oldDuration = transitions[station - 1].duration() + transitions[station].duration();

                    for (std::size_t lineSearch = 0; lineSearch < settings.lineSearchSteps; ++lineSearch) {
                        ++result.diagnostics.lineSearchTrials;
                        const auto fraction = std::ldexp(1.0, -static_cast<int>(lineSearch));
                        const BoundaryState trial{.velocity = std::lerp(result.pieceBoundaries[station].velocity, proposedVelocity, fraction), .acceleration = std::lerp(result.pieceBoundaries[station].acceleration, proposedAcceleration, fraction)};
                        const auto &leftPiece = localPieces[station - 1];
                        const auto &rightPiece = localPieces[station];

                        if (!endpointFeasible(leftPiece, leftPiece.stations.back(), trial.velocity, trial.acceleration) || !endpointFeasible(rightPiece, rightPiece.stations.front(), trial.velocity, trial.acceleration)) {
                            continue;
                        }

                        const auto left = m_implementation->transitionPlanner.solve({leftPiece.id, leftPiece.length, result.pieceBoundaries[station - 1], trial, leftPiece.maximumVelocity, leftPiece.maximumAcceleration, leftPiece.maximumJerk});

                        if (!left) {
                            continue;
                        }

                        const auto right = m_implementation->transitionPlanner.solve({rightPiece.id, rightPiece.length, trial, result.pieceBoundaries[station + 1], rightPiece.maximumVelocity, rightPiece.maximumAcceleration, rightPiece.maximumJerk});

                        if (!right || left->duration() + right->duration() > oldDuration * (1.0 + 1e-10)) {
                            continue;
                        }

                        result.pieceBoundaries[station] = trial;
                        transitions[station - 1] = *left;
                        transitions[station] = *right;
                        ++result.diagnostics.acceptedRefinements;
                        acceptedAny = true;
                        break;
                    }
                }

                if (!acceptedAny) {
                    break;
                }
            }

            std::vector<double> corrections(pieces.size(), 1.0);
            auto maximumCorrection = 1.0;

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                const auto sampledCorrection = coupledTimeScale(localPieces[pieceIndex], transitions[pieceIndex]);
                corrections[pieceIndex] = sampledCorrection > 1.0 + 1e-9 ? sampledCorrection * 1.01 : 1.0;
                maximumCorrection = std::max(maximumCorrection, corrections[pieceIndex]);
            }

            result.timeLaw.segments.clear();
            result.diagnostics.velocitySeedDuration = 0.0;
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

            if (maximumCorrection <= 1.0 + 1e-9) {
                result.diagnostics.correctedPieces = std::ranges::count(correctedPieces, true);

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    result.diagnostics.maximumAppliedTimeScale = std::max(result.diagnostics.maximumAppliedTimeScale, pieces[pieceIndex].maximumVelocity / localPieces[pieceIndex].maximumVelocity);
                }

                return result;
            }

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                if (corrections[pieceIndex] <= 1.0 + 1e-9) {
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
