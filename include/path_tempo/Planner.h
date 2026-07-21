#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <functional>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "path_tempo/Types.h"

namespace path_tempo {
    namespace detail {
        // Sampled violations normally settle in one or two passes; eight
        // leaves headroom for coupled constraints without permitting a loop.
        inline constexpr std::size_t DEFAULT_MAXIMUM_CORRECTION_PASSES = 8;
    }

    double velocityTransitionDistance(double fromVelocity, double toVelocity, double acceleration, double jerk);
    double reachableVelocity(double fixedVelocity, double cap, double length, double acceleration, double jerk);

    struct ScalarTransitionRequest {
        PathPieceId piece = 0;
        double length = 0.0;
        BoundaryState beginning;
        BoundaryState ending;
        double maximumVelocity = 0.0;
        double maximumAcceleration = 0.0;
        double maximumJerk = 0.0;
    };

    enum class PlanningErrorCode {
        InvalidInput,
        SolverFailure,
        NonMonotoneResult,
        DirectionReversal,
        CapacityExceeded,
    };

    struct PlanningError {
        PlanningErrorCode code = PlanningErrorCode::InvalidInput;
        std::string message;
    };

    struct PathPlanningSettings {
        std::size_t maximumCorrectionPasses = detail::DEFAULT_MAXIMUM_CORRECTION_PASSES;
        bool applySampledCorrections = true;
    };

    template<std::size_t DoF>
    struct PathPlanningRequest {
        std::span<const PathPiece<DoF>> pieces;
        BoundaryState beginning;
        BoundaryState ending;
        Limits<DoF> limits;
        PathPlanningSettings settings;
    };

    struct PathPlanningDiagnostics {
        double velocitySeedDuration = 0.0;
        std::size_t correctionPasses = 0;
        std::size_t correctedPieces = 0;
        double maximumAppliedTimeScale = 1.0;
    };

    struct PlannedPath {
        TimeLaw timeLaw;
        std::vector<BoundaryState> pieceBoundaries;
        std::vector<InitialPieceLimits> pieceLimits;
        PathPlanningDiagnostics diagnostics;
    };

    using MaterializationCorrection = std::function<std::expected<std::vector<PieceCorrection>, std::string>(const PlannedPath &candidate)>;

    class ScalarTransition {
        // Ruckig's one-DoF position interface emits at most seven phases. Ten
        // leaves room for zero-duration filtering and a final residual phase.
        static constexpr std::size_t CAPACITY = 10;
        std::array<CubicTimeSegment, CAPACITY> m_segments {};
        std::size_t m_size = 0;

    public:
        BoundaryState beginning;
        BoundaryState ending;

        [[nodiscard]] std::size_t size() const noexcept;
        [[nodiscard]] bool empty() const noexcept;
        [[nodiscard]] double duration() const noexcept;
        [[nodiscard]] const CubicTimeSegment &front() const;
        [[nodiscard]] const CubicTimeSegment &back() const;
        [[nodiscard]] const CubicTimeSegment &operator[](std::size_t index) const;
        [[nodiscard]] auto begin() const noexcept { return m_segments.begin(); }
        [[nodiscard]] auto end() const noexcept { return m_segments.begin() + m_size; }

    private:
        friend class ScalarTransitionPlanner;
        bool push(const CubicTimeSegment &segment) noexcept;
    };

    class ScalarTransitionPlanner {
        struct Implementation;
        std::unique_ptr<Implementation> m_implementation;

    public:
        ScalarTransitionPlanner();
        ~ScalarTransitionPlanner();
        ScalarTransitionPlanner(ScalarTransitionPlanner &&) noexcept;
        ScalarTransitionPlanner &operator=(ScalarTransitionPlanner &&) noexcept;
        ScalarTransitionPlanner(const ScalarTransitionPlanner &) = delete;
        ScalarTransitionPlanner &operator=(const ScalarTransitionPlanner &) = delete;

        std::expected<ScalarTransition, PlanningError> solve(const ScalarTransitionRequest &request);
    };

    class PathPlanner {
        struct Implementation;
        std::unique_ptr<Implementation> m_implementation;

        struct LocalPiece {
            PathPieceId id = 0;
            double length = 0.0;
            double maximumVelocity = 0.0;
            double maximumAcceleration = 0.0;
            double maximumJerk = 0.0;
            struct Station {
                double distance = 0.0;
                std::span<const double> tangent;
                std::span<const double> curvature;
                std::span<const double> thirdDerivative;
            };
            std::vector<Station> stations;
        };

        struct CoupledLimits {
            double pathAcceleration = 0.0;
            double pathJerk = 0.0;
            std::vector<double> coordinateVelocity;
            std::vector<double> coordinateAcceleration;
            std::vector<double> coordinateJerk;
        };

        std::expected<PlannedPath, PlanningError> solveLocal(std::span<const LocalPiece> pieces, BoundaryState beginning, BoundaryState ending, const CoupledLimits &limits, const PathPlanningSettings &settings, const MaterializationCorrection &materializationCorrection);

    public:
        PathPlanner();
        ~PathPlanner();
        PathPlanner(PathPlanner &&) noexcept;
        PathPlanner &operator=(PathPlanner &&) noexcept;
        PathPlanner(const PathPlanner &) = delete;
        PathPlanner &operator=(const PathPlanner &) = delete;

        template<std::size_t DoF>
        std::expected<PlannedPath, PlanningError> solve(const PathPlanningRequest<DoF> &request, const MaterializationCorrection &materializationCorrection = {});
    };
}

#include "path_tempo/detail/PathPlanner.inl"
