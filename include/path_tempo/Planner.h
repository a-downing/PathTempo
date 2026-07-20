#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <memory>
#include <string>

#include "path_tempo/Types.h"

namespace path_tempo {
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

    class ScalarTransition {
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
}
