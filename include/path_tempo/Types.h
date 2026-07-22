#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace path_tempo {
    using PathPieceId = std::uint64_t;

    template<std::size_t DoF>
    using Vector = std::array<double, DoF>;

    template<std::size_t DoF>
    struct DifferentialStation {
        double distance = 0.0;
        Vector<DoF> tangent {};
        Vector<DoF> curvature {};
        Vector<DoF> thirdDerivative {};
    };

    struct InitialPieceLimits {
        double velocity = std::numeric_limits<double>::infinity();
        double acceleration = std::numeric_limits<double>::infinity();
        double jerk = std::numeric_limits<double>::infinity();
    };

    template<std::size_t DoF>
    struct PathPiece {
        PathPieceId id = 0;
        double length = 0.0;
        double maxVelocity = 0.0;
        InitialPieceLimits initialLimits;
        std::span<const DifferentialStation<DoF>> stations;
    };

    template<std::size_t DoF>
    struct SampledPathPiece {
        PathPieceId id = 0;
        double length = 0.0;
        double maxVelocity = 0.0;
        InitialPieceLimits initialLimits;
        std::vector<DifferentialStation<DoF>> stations;

        PathPiece<DoF> view() const {
            return {id, length, maxVelocity, initialLimits, stations};
        }
    };

    struct BoundaryState {
        double velocity = 0.0;
        double acceleration = 0.0;
    };

    template<std::size_t DoF>
    struct Limits {
        double pathAcceleration = std::numeric_limits<double>::infinity();
        double pathJerk = std::numeric_limits<double>::infinity();
        Vector<DoF> coordinateVelocity {};
        Vector<DoF> coordinateAcceleration {};
        Vector<DoF> coordinateJerk {};
    };

    struct CubicTimeSegment {
        PathPieceId piece = 0;
        double duration = 0.0;
        double c0 = 0.0;
        double c1 = 0.0;
        double c2 = 0.0;
        double c3 = 0.0;

        double position(const double time) const {
            return ((c3 * time + c2) * time + c1) * time + c0;
        }

        double velocity(const double time) const {
            return (3.0 * c3 * time + 2.0 * c2) * time + c1;
        }

        double acceleration(const double time) const {
            return 6.0 * c3 * time + 2.0 * c2;
        }

        double jerk() const {
            return 6.0 * c3;
        }
    };

    struct TimeLaw {
        std::vector<CubicTimeSegment> segments;
        BoundaryState beginning;
        BoundaryState ending;
    };

    struct PieceCorrection {
        PathPieceId piece = 0;
        double requiredTimeScale = 1.0;
    };
}
