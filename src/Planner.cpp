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
        // Safeguarded Newton inversion normally converges in a few steps. If it
        // does not, forty bisections retain the old near-machine resolution.
        constexpr unsigned STATION_TIME_NEWTON_ITERATIONS = 12;
        constexpr unsigned STATION_TIME_INVERSION_MAXIMUM_ITERATIONS = 52;
        constexpr double STATION_TIME_INVERSION_EPSILON_MULTIPLIER = 8.0;
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
        // Limit the proposed difference between a boundary acceleration and a
        // piece's average acceleration to a conservative fraction of the
        // jerk-limited change available during its crossing time.
        constexpr double SLACK_TENT_FRACTION = 1.0 / 16.0;
        // At the acceleration bound the jerk tent must stay inside the bound,
        // so bang regions get only half the nominal slack.
        constexpr double SLACK_BOUND_FRACTION = 1.0 / 32.0;
        // Seeded slopes stay this fraction inside the acceleration limit so the
        // per-piece solver always has room to structure its jerk phases.
        constexpr double ACCELERATION_SEED_RESERVE = 1e-2;
        // Alternating Gauss-Seidel sweeps converge the slope-rate proposal;
        // the cap bounds the loop while leaving headroom for long runs of tiny
        // corner pieces to propagate the projection.
        constexpr std::size_t SLOPE_SMOOTHING_MAXIMUM_SWEEPS = 64;
        // Smoothing stops once the residual slope-rate violation is below
        // meaningful acceleration resolution.
        constexpr double SLOPE_SMOOTHING_CONVERGENCE_TOLERANCE = 1e-9;
        // Failed non-zero boundary-state proposals are repeatedly moved halfway
        // toward the conservative zero-acceleration profile. Extra rounds leave
        // room for those reductions to propagate along the complete path.
        constexpr std::size_t BOUNDARY_REPAIR_EXTRA_ROUNDS = 24;
        constexpr double BOUNDARY_REPAIR_MINIMUM_WEIGHT = 1.0 / (1u << 20);
        // Retaining sixty-five percent per repair step gives useful proposal
        // resolution without requiring the many exact solves of a finer decay.
        constexpr double BOUNDARY_REPAIR_WEIGHT_FACTOR = 0.65;
        // Alternating local refinements propagate improvements through the
        // look-ahead window. Six bisections resolve each boundary weight to
        // better than two percent of its remaining interval.
        constexpr std::size_t BOUNDARY_REFINEMENT_MAXIMUM_SWEEPS = 8;
        constexpr unsigned BOUNDARY_REFINEMENT_BISECTIONS = 6;
        constexpr double BOUNDARY_REFINEMENT_DURATION_TOLERANCE = 1e-12;
        constexpr double BOUNDARY_REFINEMENT_CONVERGENCE_TOLERANCE = 1e-2;
        constexpr std::size_t BOUNDARY_STATE_REFINEMENT_MAXIMUM_SWEEPS = 12;
        constexpr double BOUNDARY_STATE_REFINEMENT_CONVERGENCE_TOLERANCE = 1e-6;

        double stationTime(const CubicTimeSegment &segment, const double requestedDistance) {
            const auto from = segment.position(0.0);
            const auto to = segment.position(segment.duration);
            const auto distance = std::clamp(requestedDistance, from, to);

            if (distance == from) {
                return 0.0;
            }
            if (distance == to) {
                return segment.duration;
            }

            auto low = 0.0;
            auto high = segment.duration;
            auto time = segment.duration * (distance - from) / (to - from);
            const auto tolerance = STATION_TIME_INVERSION_EPSILON_MULTIPLIER * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(distance));

            for (unsigned iteration = 0; iteration < STATION_TIME_INVERSION_MAXIMUM_ITERATIONS; ++iteration) {
                const auto position = segment.position(time);
                const auto error = position - distance;

                if (std::abs(error) <= tolerance) {
                    return time;
                }

                if (error < 0.0) {
                    low = time;
                } else {
                    high = time;
                }

                const auto velocity = segment.velocity(time);
                const auto newton = velocity > NUMERICAL_ZERO_TOLERANCE ? time - error / velocity : time;
                time = iteration < STATION_TIME_NEWTON_ITERATIONS && newton > low && newton < high ? newton : std::midpoint(low, high);
            }

            return std::midpoint(low, high);
        }
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

    // End-acceleration treatment used by the seeding distance models: Zero
    // pins acceleration to zero at the boundary; Free lets the model pick any
    // magnitude up to the acceleration limit.
    enum class EndAcceleration {
        Zero,
        Free,
    };

    // Minimum distance to raise velocity from fromVelocity to toVelocity
    // (toVelocity >= fromVelocity) under acceleration and jerk limits with the
    // given end-acceleration treatments.
    double acceleratingTransitionDistance(const double fromVelocity, const double toVelocity, const double acceleration, const double jerk, const EndAcceleration from, const EndAcceleration to) {
        const auto change = toVelocity - fromVelocity;

        if (change <= NUMERICAL_ZERO_TOLERANCE) {
            return 0.0;
        }

        if (from == EndAcceleration::Free && to == EndAcceleration::Free) {
            return (toVelocity * toVelocity - fromVelocity * fromVelocity) / (2.0 * acceleration);
        }

        if (from == EndAcceleration::Zero && to == EndAcceleration::Zero) {
            return velocityTransitionDistance(fromVelocity, toVelocity, acceleration, jerk);
        }

        const auto rampChange = acceleration * acceleration / (2.0 * jerk);
        const auto rampDuration = acceleration / jerk;

        if (from == EndAcceleration::Zero) {
            // Acceleration ramps up from zero, then holds at the limit.
            if (change >= rampChange) {
                const auto rampDistance = fromVelocity * rampDuration + acceleration * rampDuration * rampDuration / 6.0;
                const auto holdDuration = (change - rampChange) / acceleration;
                const auto holdEntry = fromVelocity + rampChange;

                return rampDistance + holdEntry * holdDuration + acceleration * holdDuration * holdDuration / 2.0;
            }

            const auto duration = std::sqrt(2.0 * change / jerk);

            return fromVelocity * duration + jerk * duration * duration * duration / 6.0;
        }

        // Acceleration holds at the limit, then ramps down to zero.
        if (change >= rampChange) {
            const auto holdDuration = (change - rampChange) / acceleration;
            const auto holdDistance = fromVelocity * holdDuration + acceleration * holdDuration * holdDuration / 2.0;
            const auto rampEntry = fromVelocity + acceleration * holdDuration;

            return holdDistance + rampEntry * rampDuration + acceleration * rampDuration * rampDuration / 3.0;
        }

        const auto peak = std::sqrt(2.0 * change * jerk);
        const auto duration = peak / jerk;

        return fromVelocity * duration + peak * duration * duration / 2.0 - jerk * duration * duration * duration / 6.0;
    }

    // Minimum transition distance between two velocities with the given
    // end-acceleration treatments; deceleration mirrors acceleration in time.
    double freeBoundaryTransitionDistance(const double fromVelocity, const double toVelocity, const double acceleration, const double jerk, const EndAcceleration from, const EndAcceleration to) {
        if (toVelocity >= fromVelocity) {
            return acceleratingTransitionDistance(fromVelocity, toVelocity, acceleration, jerk, from, to);
        }

        return acceleratingTransitionDistance(toVelocity, fromVelocity, acceleration, jerk, to, from);
    }

    // Highest end velocity at or above fromVelocity reachable over length with
    // the modeled end accelerations.
    double reachableEndVelocity(const double fromVelocity, const double cap, const double length, const double acceleration, const double jerk, const EndAcceleration from, const EndAcceleration to) {
        if (cap <= fromVelocity) {
            return cap;
        }

        const auto available = std::max(0.0, length - std::max(REACHABILITY_ABSOLUTE_DISTANCE_MARGIN, length * REACHABILITY_RELATIVE_DISTANCE_MARGIN));

        if (freeBoundaryTransitionDistance(fromVelocity, cap, acceleration, jerk, from, to) <= available) {
            return cap;
        }

        auto low = fromVelocity;
        auto high = cap;

        for (unsigned iteration = 0; iteration < VELOCITY_BISECTION_ITERATIONS; ++iteration) {
            const auto middle = std::midpoint(low, high);

            if (freeBoundaryTransitionDistance(fromVelocity, middle, acceleration, jerk, from, to) <= available) {
                low = middle;
            } else {
                high = middle;
            }
        }

        return low;
    }

    // Highest start velocity at or above toVelocity that can still reach the
    // given end velocity over length with the modeled end accelerations.
    double reachableStartVelocity(const double toVelocity, const double cap, const double length, const double acceleration, const double jerk, const EndAcceleration from, const EndAcceleration to) {
        if (cap <= toVelocity) {
            return cap;
        }

        const auto available = std::max(0.0, length - std::max(REACHABILITY_ABSOLUTE_DISTANCE_MARGIN, length * REACHABILITY_RELATIVE_DISTANCE_MARGIN));

        if (freeBoundaryTransitionDistance(cap, toVelocity, acceleration, jerk, from, to) <= available) {
            return cap;
        }

        auto low = toVelocity;
        auto high = cap;

        for (unsigned iteration = 0; iteration < VELOCITY_BISECTION_ITERATIONS; ++iteration) {
            const auto middle = std::midpoint(low, high);

            if (freeBoundaryTransitionDistance(middle, toVelocity, acceleration, jerk, from, to) <= available) {
                low = middle;
            } else {
                high = middle;
            }
        }

        return low;
    }

    // Piece seeding data reduced to the scalar timing inputs.
    struct SeedPiece {
        double length = 0.0;
        double maximumVelocity = 0.0;
        double maximumAcceleration = 0.0;
        double maximumJerk = 0.0;
    };

    // Acceleration bound applied while seeding so requested profiles stay
    // strictly inside the per-piece limit.
    double seedAcceleration(const SeedPiece &piece) {
        return piece.maximumAcceleration * (1.0 - ACCELERATION_SEED_RESERVE);
    }

    // Highest deceleration slope a piece may end with: the per-piece solver
    // requires the peak velocity of the final jerk ramp, v + a^2/(2j), to stay
    // under the velocity cap, so |slope| <= sqrt(2 j (vMax - v)).
    double decelerationSlopeCap(const SeedPiece &piece, const double endVelocity) {
        return std::sqrt(2.0 * piece.maximumJerk * std::max(piece.maximumVelocity - endVelocity, 0.0));
    }

    // Forward/backward boundary-velocity envelope for squared velocities under
    // the per-boundary end-acceleration plan. Reductions alternate until they
    // reach the greatest fixed point of the per-piece reachability bounds.
    std::expected<std::vector<double>, PlanningError> boundaryVelocityEnvelope(const std::span<const SeedPiece> pieces, const std::span<const double> stationCaps, const std::span<const EndAcceleration> boundaryPlan, const BoundaryState &beginning, const BoundaryState &ending) {
        const auto stationCount = pieces.size() + 1;
        std::vector<double> squared(stationCount, 0.0);

        for (std::size_t station = 0; station < stationCount; ++station) {
            squared[station] = stationCaps[station] * stationCaps[station];
        }

        squared.front() = beginning.velocity * beginning.velocity;
        squared.back() = ending.velocity * ending.velocity;
        const auto reachabilityPasses = 2 * pieces.size() + REACHABILITY_EXTRA_PASSES;

        for (std::size_t pass = 0; pass < reachabilityPasses; ++pass) {
            auto maximumChange = 0.0;

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                const auto &piece = pieces[pieceIndex];
                const auto acceleration = seedAcceleration(piece);
                double reachable = 0.0;

                if (boundaryPlan[pieceIndex] == EndAcceleration::Free && boundaryPlan[pieceIndex + 1] == EndAcceleration::Free) {
                    reachable = std::min(std::sqrt(squared[pieceIndex] + 2.0 * acceleration * piece.length), stationCaps[pieceIndex + 1]);
                } else {
                    reachable = reachableEndVelocity(std::sqrt(squared[pieceIndex]), stationCaps[pieceIndex + 1], piece.length, acceleration, piece.maximumJerk, boundaryPlan[pieceIndex], boundaryPlan[pieceIndex + 1]);
                }

                if (pieceIndex + 1 == pieces.size()) {
                    if (reachable + FIXED_ENDPOINT_REACHABILITY_TOLERANCE < ending.velocity) {
                        return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = "fixed ending velocity is not forward reachable"});
                    }
                } else {
                    const auto reduced = std::min(squared[pieceIndex + 1], reachable * reachable);
                    maximumChange = std::max(maximumChange, squared[pieceIndex + 1] - reduced);
                    squared[pieceIndex + 1] = reduced;
                }
            }

            for (std::size_t pieceIndex = pieces.size(); pieceIndex-- > 0;) {
                const auto &piece = pieces[pieceIndex];
                const auto acceleration = seedAcceleration(piece);
                double reachable = 0.0;

                if (boundaryPlan[pieceIndex] == EndAcceleration::Free && boundaryPlan[pieceIndex + 1] == EndAcceleration::Free) {
                    reachable = std::min(std::sqrt(squared[pieceIndex + 1] + 2.0 * acceleration * piece.length), stationCaps[pieceIndex]);
                } else {
                    reachable = reachableStartVelocity(std::sqrt(squared[pieceIndex + 1]), stationCaps[pieceIndex], piece.length, acceleration, piece.maximumJerk, boundaryPlan[pieceIndex], boundaryPlan[pieceIndex + 1]);
                }

                if (pieceIndex == 0) {
                    if (reachable + FIXED_ENDPOINT_REACHABILITY_TOLERANCE < beginning.velocity) {
                        return std::unexpected(PlanningError {.code = PlanningErrorCode::SolverFailure, .message = "fixed beginning velocity cannot reach the remaining path"});
                    }
                } else {
                    const auto reduced = std::min(squared[pieceIndex], reachable * reachable);
                    maximumChange = std::max(maximumChange, squared[pieceIndex] - reduced);
                    squared[pieceIndex] = reduced;
                }
            }

            if (maximumChange <= REACHABILITY_CONVERGENCE_TOLERANCE) {
                break;
            }
        }

        return squared;
    }

    // How far a proposed boundary acceleration may deviate from a piece's own
    // average acceleration. Exact feasibility is checked by Ruckig later.
    double pieceSlack(const SeedPiece &piece, const double slope, const double crossingTime) {
        const auto tent = piece.maximumJerk * crossingTime;
        const auto headroom = std::max(0.0, seedAcceleration(piece) - std::abs(slope));

        return std::min(tent * SLACK_TENT_FRACTION, headroom + tent * SLACK_BOUND_FRACTION);
    }

    // Smooths the squared-velocity envelope so adjacent piece-average
    // accelerations produce useful continuous boundary-acceleration proposals.
    // Interior boundary values move within their energy-consistent windows,
    // keeping the fixed end states and per-piece acceleration bounds intact.
    void smoothVelocityEnvelope(const std::span<const SeedPiece> pieces, const std::span<const double> stationCaps, const std::span<double> squared) {
        const auto pieceCount = pieces.size();
        const auto stationCount = pieceCount + 1;

        const auto smoothBoundary = [&](const std::size_t station) -> double {
            const auto &left = pieces[station - 1];
            const auto &right = pieces[station];
            const auto velocityLeft = std::sqrt(squared[station - 1]);
            const auto velocity = std::sqrt(squared[station]);
            const auto velocityRight = std::sqrt(squared[station + 1]);
            const auto slopeLeft = (squared[station] - squared[station - 1]) / (2.0 * left.length);
            const auto slopeRight = (squared[station + 1] - squared[station]) / (2.0 * right.length);
            const auto timeLeft = velocityLeft + velocity > NUMERICAL_ZERO_TOLERANCE ? 2.0 * left.length / (velocityLeft + velocity) : std::numeric_limits<double>::infinity();
            const auto timeRight = velocity + velocityRight > NUMERICAL_ZERO_TOLERANCE ? 2.0 * right.length / (velocity + velocityRight) : std::numeric_limits<double>::infinity();
            const auto rate = pieceSlack(left, slopeLeft, timeLeft) + pieceSlack(right, slopeRight, timeRight);
            const auto gap = slopeRight - slopeLeft;
            const auto violation = std::abs(gap) - rate;

            if (violation <= 0.0) {
                return 0.0;
            }

            const auto excess = gap - std::copysign(rate, gap);
            const auto adjusted = squared[station] + excess * (2.0 * left.length * right.length) / (left.length + right.length);

            // Keep the boundary inside the energy window shared with both
            // neighbors so per-piece acceleration bounds and caps survive.
            const auto floorLeft = std::max(0.0, squared[station - 1] - 2.0 * seedAcceleration(left) * left.length);
            const auto floorRight = std::max(0.0, squared[station + 1] - 2.0 * seedAcceleration(right) * right.length);
            const auto ceilLeft = squared[station - 1] + 2.0 * seedAcceleration(left) * left.length;
            const auto ceilRight = squared[station + 1] + 2.0 * seedAcceleration(right) * right.length;
            const auto lower = std::max(floorLeft, floorRight);
            const auto upper = std::min({stationCaps[station] * stationCaps[station], ceilLeft, ceilRight});

            if (lower > upper) {
                return violation;
            }

            const auto clamped = std::clamp(adjusted, lower, upper);
            squared[station] = clamped;

            return violation;
        };

        for (std::size_t sweep = 0; sweep < SLOPE_SMOOTHING_MAXIMUM_SWEEPS; ++sweep) {
            auto maximumViolation = 0.0;

            if (sweep % 2 == 0) {
                for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                    maximumViolation = std::max(maximumViolation, smoothBoundary(station));
                }
            } else {
                for (std::size_t station = stationCount - 1; station-- > 1;) {
                    maximumViolation = std::max(maximumViolation, smoothBoundary(station));
                }
            }

            if (maximumViolation <= SLOPE_SMOOTHING_CONVERGENCE_TOLERANCE) {
                break;
            }
        }
    }

    // Proposes scalar boundary accelerations from the smoothed velocity-profile
    // slopes, clamped by the shared acceleration limit, cap-binding signs, and
    // the solver's peak-velocity reserves. Exact adjacent-piece feasibility is
    // established later by the transition solver.
    void assignBoundaryAccelerations(const std::span<const SeedPiece> pieces, const std::span<const double> squared, const std::span<BoundaryState> boundaries) {
        const auto stationCount = boundaries.size();

        for (std::size_t station = 1; station + 1 < stationCount; ++station) {
            const auto &left = pieces[station - 1];
            const auto &right = pieces[station];
            const auto velocity = std::sqrt(squared[station]);

            if (velocity <= NUMERICAL_ZERO_TOLERANCE) {
                // A full stop may neither accelerate into the next piece nor
                // decelerate out of the previous one.
                boundaries[station].acceleration = 0.0;
                continue;
            }

            const auto velocityLeft = std::sqrt(squared[station - 1]);
            const auto velocityRight = std::sqrt(squared[station + 1]);
            const auto slopeLeft = (squared[station] - squared[station - 1]) / (2.0 * left.length);
            const auto slopeRight = (squared[station + 1] - squared[station]) / (2.0 * right.length);
            const auto timeLeft = velocityLeft + velocity > NUMERICAL_ZERO_TOLERANCE ? 2.0 * left.length / (velocityLeft + velocity) : std::numeric_limits<double>::infinity();
            const auto timeRight = velocity + velocityRight > NUMERICAL_ZERO_TOLERANCE ? 2.0 * right.length / (velocity + velocityRight) : std::numeric_limits<double>::infinity();
            const auto slackLeft = pieceSlack(left, slopeLeft, timeLeft);
            const auto slackRight = pieceSlack(right, slopeRight, timeRight);
            const auto accelerationLimit = std::min(left.maximumAcceleration, right.maximumAcceleration);
            const auto lower = std::max({slopeLeft - slackLeft, slopeRight - slackRight, -accelerationLimit});
            const auto upper = std::min({slopeLeft + slackLeft, slopeRight + slackRight, accelerationLimit});
            const auto desired = 0.5 * (slopeLeft + slopeRight);
            auto acceleration = lower <= upper ? std::clamp(desired, lower, upper) : std::clamp(desired, -accelerationLimit, accelerationLimit);

            // A boundary riding on one adjacent piece's velocity cap may only
            // accelerate away from that piece's side.
            if (velocity >= left.maximumVelocity * (1.0 - CONSTRAINT_RELATIVE_TOLERANCE)) {
                acceleration = std::max(acceleration, 0.0);
            }

            if (velocity >= right.maximumVelocity * (1.0 - CONSTRAINT_RELATIVE_TOLERANCE)) {
                acceleration = std::min(acceleration, 0.0);
            }

            // The per-piece solver validates its target state's peak velocity
            // during the final jerk ramp against the piece velocity cap.
            if (acceleration < 0.0) {
                acceleration = std::max(acceleration, -decelerationSlopeCap(left, velocity) * (1.0 - CONSTRAINT_RELATIVE_TOLERANCE));
            }

            boundaries[station].acceleration = acceleration;
        }
    }

    // Conservative zero-acceleration boundary velocities used as the feasible
    // fallback floor for local repair.
    std::expected<std::vector<double>, PlanningError> conservativeVelocityFloor(const std::span<const SeedPiece> pieces, const std::span<const double> stationCaps, const BoundaryState &beginning, const BoundaryState &ending) {
        std::vector<double> stationVelocity(stationCaps.begin(), stationCaps.end());
        const auto reachabilityPasses = 2 * pieces.size() + REACHABILITY_EXTRA_PASSES;

        for (std::size_t pass = 0; pass < reachabilityPasses; ++pass) {
            auto maximumChange = 0.0;

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                const auto &piece = pieces[pieceIndex];
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
                const auto &piece = pieces[pieceIndex];
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

        return stationVelocity;
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

        if (!std::isfinite(request.length) || !std::isfinite(from.velocity) || !std::isfinite(from.acceleration) || !std::isfinite(to.velocity) || !std::isfinite(to.acceleration) || !std::isfinite(request.maximumVelocity) || std::isnan(request.maximumAcceleration) || std::isnan(request.maximumJerk) || request.length <= 0.0 || from.velocity < 0.0 || to.velocity < 0.0 || request.maximumVelocity <= 0.0 || request.maximumAcceleration <= 0.0 || request.maximumJerk <= 0.0 || from.velocity > request.maximumVelocity * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || to.velocity > request.maximumVelocity * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || std::abs(from.acceleration) > request.maximumAcceleration * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE) || std::abs(to.acceleration) > request.maximumAcceleration * (1.0 + CONSTRAINT_RELATIVE_TOLERANCE)) {
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

        auto invalidPhase = false;
        const auto appendPhase = [&](const double duration, const double jerk, const bool initialStateCorrection = false) -> bool {
            if (!std::isfinite(duration) || !std::isfinite(jerk) || duration < -MINIMUM_PHASE_DURATION) {
                invalidPhase = true;

                return false;
            }

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
                    .code = invalidPhase ? PlanningErrorCode::SolverFailure : PlanningErrorCode::CapacityExceeded,
                    .message = invalidPhase ? "Ruckig produced an invalid local phase" : "local time law exceeded its fixed phase capacity",
                });
            }
        }

        for (std::size_t phase = 0; phase < profile.t.size(); ++phase) {
            if (!appendPhase(profile.t[phase], profile.j[phase])) {
                return std::unexpected(PlanningError {
                    .code = invalidPhase ? PlanningErrorCode::SolverFailure : PlanningErrorCode::CapacityExceeded,
                    .message = invalidPhase ? "Ruckig produced an invalid local phase" : "local time law exceeded its fixed phase capacity",
                });
            }
        }

        if (!appendPhase(trajectory.get_duration() - result.duration(), 0.0)) {
            return std::unexpected(PlanningError {
                .code = invalidPhase ? PlanningErrorCode::SolverFailure : PlanningErrorCode::CapacityExceeded,
                .message = invalidPhase ? "Ruckig produced an invalid local phase" : "local time law exceeded its fixed phase capacity",
            });
        }

        const auto distanceTolerance = std::max(TRANSITION_ABSOLUTE_TOLERANCE, request.length * TRANSITION_DISTANCE_RELATIVE_TOLERANCE);
        const auto velocityTolerance = std::max(TRANSITION_ABSOLUTE_TOLERANCE, request.maximumVelocity * TRANSITION_VELOCITY_RELATIVE_TOLERANCE);
        const auto accelerationTolerance = std::max(TRANSITION_ABSOLUTE_TOLERANCE, request.maximumAcceleration * TRANSITION_VELOCITY_RELATIVE_TOLERANCE);
        auto previousTime = 0.0;
        auto previousDistance = 0.0;

        for (std::size_t phase = 0; phase < result.size(); ++phase) {
            const auto &segment = result[phase];
            const auto time = previousTime + segment.duration;
            const auto phaseDistance = segment.position(segment.duration);
            const auto phaseVelocity = segment.velocity(segment.duration);
            const auto phaseAcceleration = segment.acceleration(segment.duration);

            if (!std::isfinite(time) || !std::isfinite(phaseDistance) || !std::isfinite(phaseVelocity) || !std::isfinite(phaseAcceleration) || phaseDistance < -distanceTolerance || phaseDistance > request.length + distanceTolerance || phaseVelocity < -velocityTolerance || time <= previousTime || phaseDistance < previousDistance - distanceTolerance) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::NonMonotoneResult,
                    .message = std::format("Ruckig produced a non-monotone local path law at boundary {} of {}: time={} distance={} velocity={} acceleration={} previous time={} previous distance={} length={} distance tolerance={} velocity tolerance={}", phase + 1, result.size() + 1, time, phaseDistance, phaseVelocity, phaseAcceleration, previousTime, previousDistance, request.length, distanceTolerance, velocityTolerance),
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

        if (result.empty() || std::abs(previousDistance - request.length) > distanceTolerance || std::abs(result.back().velocity(result.back().duration) - to.velocity) > velocityTolerance || std::abs(result.back().acceleration(result.back().duration) - to.acceleration) > accelerationTolerance) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::SolverFailure,
                .message = std::format("Ruckig local position timing missed its requested endpoint: distance={} requested={} velocity={} requested={} acceleration={} requested={}", previousDistance, request.length, result.empty() ? from.velocity : result.back().velocity(result.back().duration), to.velocity, result.empty() ? from.acceleration : result.back().acceleration(result.back().duration), to.acceleration),
            });
        }

        return result;
    }

    PathPlanner::PathPlanner() : m_implementation(std::make_unique<Implementation>()) {
    }

    PathPlanner::~PathPlanner() = default;

    PathPlanner::PathPlanner(PathPlanner &&) noexcept = default;

    PathPlanner &PathPlanner::operator=(PathPlanner &&) noexcept = default;

    std::expected<PlannedPath, PlanningError> PathPlanner::solveLocal(const std::span<const LocalPiece> pieces, const PieceIndexMap &pieceIndices, const BoundaryState beginning, const BoundaryState ending, const CoupledLimits &limits, const PathPlanningSettings &settings, const MaterializationCorrection &materializationCorrection) {
        auto correctedPieces = std::vector<bool>(pieces.size(), false);
        auto localPieces = std::vector<LocalPiece>(pieces.begin(), pieces.end());
        auto appliedSampledCorrection = false;
        auto appliedMaterializationCorrection = false;
        const auto stationCount = pieces.size() + 1;
        const auto correctionDescription = [&] {
            if (appliedSampledCorrection && appliedMaterializationCorrection) {
                return "sampled coupled-path and materialization corrections";
            }
            if (appliedMaterializationCorrection) {
                return "materialization correction";
            }

            return "sampled coupled-path constraint correction";
        };

        const auto coupledTimeScale = [&](const LocalPiece &piece, const ScalarTransition &transition) {
            if (!piece.requiresCoupledCheck) {
                return 1.0;
            }

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

                    const auto time = stationTime(segment, station.distance);
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
                    .message = correctionPass == 0 ? "path boundary state exceeds a local piece limit" : std::format("fixed path boundary state cannot satisfy the limits after {}", correctionDescription()),
                });
            }

            std::vector<SeedPiece> seedPieces;
            seedPieces.reserve(pieces.size());

            for (const auto &piece : localPieces) {
                seedPieces.push_back({piece.length, piece.maximumVelocity, piece.maximumAcceleration, piece.maximumJerk});
            }

            const auto floor = conservativeVelocityFloor(seedPieces, stationCaps, beginning, ending);

            if (!floor) {
                return std::unexpected(floor.error());
            }

            std::vector<BoundaryState> baselineBoundaryStates(stationCount);
            baselineBoundaryStates.front() = beginning;
            baselineBoundaryStates.back() = ending;

            for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                baselineBoundaryStates[station] = {.velocity = (*floor)[station]};
            }

            std::vector<ScalarTransition> baselineTransitions(pieces.size());

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                const auto &piece = localPieces[pieceIndex];
                const auto transition = m_implementation->transitionPlanner.solve({
                    .piece = piece.id,
                    .length = piece.length,
                    .beginning = baselineBoundaryStates[pieceIndex],
                    .ending = baselineBoundaryStates[pieceIndex + 1],
                    .maximumVelocity = piece.maximumVelocity,
                    .maximumAcceleration = piece.maximumAcceleration,
                    .maximumJerk = piece.maximumJerk,
                });

                if (!transition) {
                    return std::unexpected(PlanningError {
                        .code = transition.error().code,
                        .message = std::format("could not materialize conservative path piece {}: {}", pieceIndex, transition.error().message),
                    });
                }

                baselineTransitions[pieceIndex] = *transition;
            }

            std::vector<double> baselineCorrections(pieces.size(), 1.0);
            auto maximumBaselineCorrection = 1.0;

            if (settings.applySampledCorrections) {
                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    const auto sampledCorrection = coupledTimeScale(localPieces[pieceIndex], baselineTransitions[pieceIndex]);
                    baselineCorrections[pieceIndex] = sampledCorrection > 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE ? sampledCorrection * CORRECTION_SAFETY_FACTOR : 1.0;
                    maximumBaselineCorrection = std::max(maximumBaselineCorrection, baselineCorrections[pieceIndex]);
                }
            }

            if (maximumBaselineCorrection > 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE) {
                appliedSampledCorrection = true;

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    if (baselineCorrections[pieceIndex] <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE) {
                        continue;
                    }

                    correctedPieces[pieceIndex] = true;
                    const auto factor = baselineCorrections[pieceIndex];
                    localPieces[pieceIndex].maximumVelocity /= factor;
                    localPieces[pieceIndex].maximumAcceleration /= factor * factor;
                    localPieces[pieceIndex].maximumJerk /= factor * factor * factor;
                }

                continue;
            }

            PlannedPath result;
            result.timeLaw.beginning = beginning;
            result.timeLaw.ending = ending;
            result.pieceLimits.reserve(pieces.size());
            result.diagnostics.correctionPasses = correctionPass + 1;
            auto boundaryStates = baselineBoundaryStates;
            auto transitions = baselineTransitions;

            if (settings.boundaryAccelerationMode == BoundaryAccelerationMode::Optimized) {
                // Seed boundary velocities with the free-acceleration envelope,
                // rounded at the fixed path ends, then smooth the profile so its
                // acceleration changes at jerk-limited rates across boundaries.
                std::vector<EndAcceleration> boundaryPlan(stationCount, EndAcceleration::Free);

                if (std::abs(beginning.acceleration) <= NUMERICAL_ZERO_TOLERANCE) {
                    boundaryPlan.front() = EndAcceleration::Zero;
                }

                if (std::abs(ending.acceleration) <= NUMERICAL_ZERO_TOLERANCE) {
                    boundaryPlan.back() = EndAcceleration::Zero;
                }

                auto envelope = boundaryVelocityEnvelope(seedPieces, stationCaps, boundaryPlan, beginning, ending);
                std::vector<double> squared(stationCount, 0.0);

                if (envelope) {
                    squared = std::move(*envelope);
                } else {
                    for (std::size_t station = 0; station < stationCount; ++station) {
                        squared[station] = (*floor)[station] * (*floor)[station];
                    }
                }

                smoothVelocityEnvelope(seedPieces, stationCaps, squared);

                for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                    squared[station] = std::max(squared[station], (*floor)[station] * (*floor)[station]);
                }

                boundaryStates = std::vector<BoundaryState>(stationCount);
                boundaryStates.front() = beginning;
                boundaryStates.back() = ending;

                for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                    boundaryStates[station] = {.velocity = std::sqrt(squared[station])};
                }

                assignBoundaryAccelerations(seedPieces, squared, boundaryStates);
                const auto proposedBoundaryStates = boundaryStates;
                std::vector<double> boundaryWeights(stationCount, 1.0);
                boundaryWeights.front() = 0.0;
                boundaryWeights.back() = 0.0;
                const auto setBoundaryWeight = [&](const std::size_t station, const double weight) {
                    const auto floorSquared = (*floor)[station] * (*floor)[station];
                    const auto proposedSquared = proposedBoundaryStates[station].velocity * proposedBoundaryStates[station].velocity;
                    boundaryWeights[station] = weight;
                    boundaryStates[station].velocity = std::sqrt(std::lerp(floorSquared, proposedSquared, weight));
                    boundaryStates[station].acceleration = proposedBoundaryStates[station].acceleration * weight;
                };

                std::vector<std::size_t> failed;

                const auto solvePiece = [&](const std::size_t pieceIndex) -> bool {
                    const auto &piece = localPieces[pieceIndex];
                    const auto transition = m_implementation->transitionPlanner.solve({
                        .piece = piece.id,
                        .length = piece.length,
                        .beginning = boundaryStates[pieceIndex],
                        .ending = boundaryStates[pieceIndex + 1],
                        .maximumVelocity = piece.maximumVelocity,
                        .maximumAcceleration = piece.maximumAcceleration,
                        .maximumJerk = piece.maximumJerk,
                    });

                    if (!transition) {
                        return false;
                    }

                    transitions[pieceIndex] = *transition;

                    return true;
                };

                for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                    if (!solvePiece(pieceIndex)) {
                        failed.push_back(pieceIndex);
                    } else if (settings.applySampledCorrections && coupledTimeScale(localPieces[pieceIndex], transitions[pieceIndex]) > 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE) {
                        failed.push_back(pieceIndex);
                    }
                }

                // Repair infeasible proposals through a local homotopy. Each
                // affected boundary moves halfway toward the conservative profile;
                // exact transition solves decide feasibility at every step. The
                // zero weight is the original zero-acceleration construction.
                const auto maximumRepairRounds = 2 * pieces.size() + BOUNDARY_REPAIR_EXTRA_ROUNDS;

                for (std::size_t repairRound = 0; !failed.empty(); ++repairRound) {
                    if (repairRound == maximumRepairRounds) {
                        boundaryStates = baselineBoundaryStates;
                        transitions = baselineTransitions;
                        std::fill(boundaryWeights.begin() + 1, boundaryWeights.end() - 1, 0.0);
                        failed.clear();
                        break;
                    }

                    std::vector<char> reduceBoundary(stationCount, false);
                    std::vector<char> dirty(pieces.size(), false);

                    for (const auto pieceIndex : failed) {
                        dirty[pieceIndex] = true;

                        for (const auto station : {pieceIndex, pieceIndex + 1}) {
                            if (station > 0 && station < pieces.size() && boundaryWeights[station] > 0.0) {
                                reduceBoundary[station] = true;
                            }
                        }
                    }

                    auto changed = false;

                    for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                        if (!reduceBoundary[station]) {
                            continue;
                        }

                        auto &weight = boundaryWeights[station];
                        weight = weight <= BOUNDARY_REPAIR_MINIMUM_WEIGHT / BOUNDARY_REPAIR_WEIGHT_FACTOR ? 0.0 : weight * BOUNDARY_REPAIR_WEIGHT_FACTOR;
                        setBoundaryWeight(station, weight);
                        changed = true;
                        dirty[station - 1] = true;
                        dirty[station] = true;
                    }

                    if (!changed) {
                        boundaryStates = baselineBoundaryStates;
                        transitions = baselineTransitions;
                        std::fill(boundaryWeights.begin() + 1, boundaryWeights.end() - 1, 0.0);
                        failed.clear();
                        break;
                    }

                    std::vector<std::size_t> stillFailing;

                    for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                        if (!dirty[pieceIndex]) {
                            continue;
                        }

                        if (!solvePiece(pieceIndex)) {
                            stillFailing.push_back(pieceIndex);
                        } else if (settings.applySampledCorrections && coupledTimeScale(localPieces[pieceIndex], transitions[pieceIndex]) > 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE) {
                            stillFailing.push_back(pieceIndex);
                        }
                    }

                    failed = std::move(stillFailing);
                }

                const auto refineBoundary = [&](const std::size_t station) -> double {
                    if (boundaryWeights[station] >= 1.0) {
                        return 0.0;
                    }

                    const auto leftPiece = station - 1;
                    const auto rightPiece = station;
                    const auto originalWeight = boundaryWeights[station];
                    const auto originalState = boundaryStates[station];
                    const auto originalLeft = transitions[leftPiece];
                    const auto originalRight = transitions[rightPiece];
                    auto bestWeight = originalWeight;
                    auto bestState = originalState;
                    auto bestLeft = originalLeft;
                    auto bestRight = originalRight;
                    auto bestDuration = originalLeft.duration() + originalRight.duration();
                    auto low = originalWeight;
                    auto high = 1.0;

                    for (unsigned iteration = 0; iteration < BOUNDARY_REFINEMENT_BISECTIONS; ++iteration) {
                        const auto trialWeight = std::midpoint(low, high);
                        setBoundaryWeight(station, trialWeight);
                        auto feasible = solvePiece(leftPiece) && solvePiece(rightPiece);

                        if (feasible && settings.applySampledCorrections) {
                            feasible = coupledTimeScale(localPieces[leftPiece], transitions[leftPiece]) <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE && coupledTimeScale(localPieces[rightPiece], transitions[rightPiece]) <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE;
                        }

                        if (!feasible) {
                            high = trialWeight;
                            continue;
                        }

                        low = trialWeight;
                        const auto duration = transitions[leftPiece].duration() + transitions[rightPiece].duration();

                        if (duration + BOUNDARY_REFINEMENT_DURATION_TOLERANCE < bestDuration) {
                            bestWeight = trialWeight;
                            bestState = boundaryStates[station];
                            bestLeft = transitions[leftPiece];
                            bestRight = transitions[rightPiece];
                            bestDuration = duration;
                        }
                    }

                    boundaryWeights[station] = bestWeight;
                    boundaryStates[station] = bestState;
                    transitions[leftPiece] = bestLeft;
                    transitions[rightPiece] = bestRight;

                    return originalLeft.duration() + originalRight.duration() - bestDuration;
                };

                auto cycleImprovement = 0.0;

                for (std::size_t sweep = 0; sweep < BOUNDARY_REFINEMENT_MAXIMUM_SWEEPS; ++sweep) {
                    if (sweep % 2 == 0) {
                        for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                            cycleImprovement += refineBoundary(station);
                        }
                    } else {
                        for (std::size_t station = stationCount - 1; station-- > 1;) {
                            cycleImprovement += refineBoundary(station);
                        }

                        if (cycleImprovement <= BOUNDARY_REFINEMENT_CONVERGENCE_TOLERANCE) {
                            break;
                        }

                        cycleImprovement = 0.0;
                    }
                }

                const auto refineAcceleration = [&](const std::size_t station) -> double {
                    const auto leftPiece = station - 1;
                    const auto rightPiece = station;
                    const auto originalState = boundaryStates[station];
                    const auto originalLeft = transitions[leftPiece];
                    const auto originalRight = transitions[rightPiece];
                    const auto accelerationLimit = std::min(localPieces[leftPiece].maximumAcceleration, localPieces[rightPiece].maximumAcceleration);
                    const std::array candidates{
                        0.0,
                        originalState.acceleration / 2.0,
                        std::midpoint(originalState.acceleration, proposedBoundaryStates[station].acceleration),
                        proposedBoundaryStates[station].acceleration,
                        std::clamp(originalState.acceleration - accelerationLimit / 4.0, -accelerationLimit, accelerationLimit),
                        std::clamp(originalState.acceleration + accelerationLimit / 4.0, -accelerationLimit, accelerationLimit),
                    };
                    auto bestState = originalState;
                    auto bestLeft = originalLeft;
                    auto bestRight = originalRight;
                    auto bestDuration = originalLeft.duration() + originalRight.duration();

                    for (std::size_t candidateIndex = 0; candidateIndex < candidates.size(); ++candidateIndex) {
                        const auto acceleration = candidates[candidateIndex];
                        const auto preceding = std::span {candidates}.first(candidateIndex);

                        if (acceleration == originalState.acceleration || std::ranges::find(preceding, acceleration) != preceding.end()) {
                            continue;
                        }

                        boundaryStates[station] = originalState;
                        boundaryStates[station].acceleration = acceleration;
                        auto feasible = solvePiece(leftPiece) && solvePiece(rightPiece);

                        if (feasible && settings.applySampledCorrections) {
                            feasible = coupledTimeScale(localPieces[leftPiece], transitions[leftPiece]) <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE && coupledTimeScale(localPieces[rightPiece], transitions[rightPiece]) <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE;
                        }

                        if (!feasible) {
                            continue;
                        }

                        const auto duration = transitions[leftPiece].duration() + transitions[rightPiece].duration();

                        if (duration + BOUNDARY_REFINEMENT_DURATION_TOLERANCE < bestDuration) {
                            bestState = boundaryStates[station];
                            bestLeft = transitions[leftPiece];
                            bestRight = transitions[rightPiece];
                            bestDuration = duration;
                        }
                    }

                    boundaryStates[station] = bestState;
                    transitions[leftPiece] = bestLeft;
                    transitions[rightPiece] = bestRight;

                    return originalLeft.duration() + originalRight.duration() - bestDuration;
                };

                const auto refineVelocity = [&](const std::size_t station) -> double {
                    const auto leftPiece = station - 1;
                    const auto rightPiece = station;
                    const auto originalState = boundaryStates[station];
                    const auto originalLeft = transitions[leftPiece];
                    const auto originalRight = transitions[rightPiece];
                    auto bestState = originalState;
                    auto bestLeft = originalLeft;
                    auto bestRight = originalRight;
                    auto bestDuration = originalLeft.duration() + originalRight.duration();
                    auto low = originalState.velocity;
                    auto high = proposedBoundaryStates[station].velocity;

                    if (high <= low + NUMERICAL_ZERO_TOLERANCE) {
                        return 0.0;
                    }

                    for (unsigned iteration = 0; iteration < BOUNDARY_REFINEMENT_BISECTIONS; ++iteration) {
                        const auto velocity = std::midpoint(low, high);
                        boundaryStates[station] = originalState;
                        boundaryStates[station].velocity = velocity;
                        auto feasible = solvePiece(leftPiece) && solvePiece(rightPiece);

                        if (feasible && settings.applySampledCorrections) {
                            feasible = coupledTimeScale(localPieces[leftPiece], transitions[leftPiece]) <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE && coupledTimeScale(localPieces[rightPiece], transitions[rightPiece]) <= 1.0 + CORRECTION_SIGNIFICANCE_TOLERANCE;
                        }

                        if (!feasible) {
                            high = velocity;
                            continue;
                        }

                        low = velocity;
                        const auto duration = transitions[leftPiece].duration() + transitions[rightPiece].duration();

                        if (duration + BOUNDARY_REFINEMENT_DURATION_TOLERANCE < bestDuration) {
                            bestState = boundaryStates[station];
                            bestLeft = transitions[leftPiece];
                            bestRight = transitions[rightPiece];
                            bestDuration = duration;
                        }
                    }

                    boundaryStates[station] = bestState;
                    transitions[leftPiece] = bestLeft;
                    transitions[rightPiece] = bestRight;

                    return originalLeft.duration() + originalRight.duration() - bestDuration;
                };

                auto stateCycleImprovement = 0.0;

                for (std::size_t sweep = 0; sweep < BOUNDARY_STATE_REFINEMENT_MAXIMUM_SWEEPS; ++sweep) {
                    const auto refineState = [&](const std::size_t station) {
                        stateCycleImprovement += refineAcceleration(station);
                        stateCycleImprovement += refineVelocity(station);
                    };

                    if (sweep % 2 == 0) {
                        for (std::size_t station = 1; station + 1 < stationCount; ++station) {
                            refineState(station);
                        }
                    } else {
                        for (std::size_t station = stationCount - 1; station-- > 1;) {
                            refineState(station);
                        }

                        if (stateCycleImprovement <= BOUNDARY_STATE_REFINEMENT_CONVERGENCE_TOLERANCE) {
                            break;
                        }

                        stateCycleImprovement = 0.0;
                    }
                }

                const auto transitionDuration = [](const std::span<const ScalarTransition> candidate) {
                    return std::accumulate(candidate.begin(), candidate.end(), 0.0, [](const double total, const ScalarTransition &transition) { return total + transition.duration(); });
                };

                if (transitionDuration(transitions) >= transitionDuration(baselineTransitions)) {
                    boundaryStates = baselineBoundaryStates;
                    transitions = baselineTransitions;
                }
            }

            result.pieceBoundaries = boundaryStates;

            for (const auto &piece : localPieces) {
                result.pieceLimits.push_back({
                    .velocity = piece.maximumVelocity,
                    .acceleration = piece.maximumAcceleration,
                    .jerk = piece.maximumJerk,
                });
            }

            std::vector<double> corrections(pieces.size(), 1.0);
            auto maximumCorrection = 1.0;

            auto refinedPathDistance = 0.0;

            for (std::size_t pieceIndex = 0; pieceIndex < pieces.size(); ++pieceIndex) {
                result.diagnostics.trajectoryDuration += transitions[pieceIndex].duration();

                for (auto segment : transitions[pieceIndex]) {
                    segment.c0 += refinedPathDistance;
                    result.timeLaw.segments.push_back(segment);
                }

                refinedPathDistance += pieces[pieceIndex].length;
            }

            result.diagnostics.velocitySeedDuration = result.diagnostics.trajectoryDuration;

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
                    const auto found = pieceIndices.find(correction.piece);

                    if (found == pieceIndices.end() || !std::isfinite(correction.requiredTimeScale) || correction.requiredTimeScale < 1.0) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("materialization returned an invalid correction for piece {} with time scale {}", correction.piece, correction.requiredTimeScale),
                        });
                    }

                    const auto pieceIndex = found->second;

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

            appliedMaterializationCorrection = true;

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
            .message = std::format("{} did not converge after {} passes", correctionDescription(), settings.maximumCorrectionPasses),
        });
    }
}
