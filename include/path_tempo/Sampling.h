#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>

#include "path_tempo/Types.h"

namespace path_tempo {
    enum class SamplingError {
        InvalidLength,
        InvalidMaxVelocity,
        NonFiniteGeometry,
        IrregularCurve,
    };

    template<std::size_t DoF>
    struct DifferentialState {
        Vector<DoF> tangent {};
        Vector<DoF> curvature {};
        Vector<DoF> thirdDerivative {};
    };

    template<std::size_t DoF>
    struct Line {
        Vector<DoF> from {};
        Vector<DoF> to {};
    };

    namespace detail {
        template<std::size_t DoF>
        double magnitude(const Vector<DoF> &value) {
            auto squared = 0.0;

            for (const auto component : value) {
                squared += component * component;
            }

            return std::sqrt(squared);
        }

        template<std::size_t DoF>
        bool finite(const Vector<DoF> &value) {
            return std::ranges::all_of(value, [](const double component) {
                return std::isfinite(component);
            });
        }
    }

    template<std::size_t DoF>
    std::expected<SampledPathPiece<DoF>, SamplingError> sampleLine(const Line<DoF> &line, const PathPieceId id, const double maxVelocity) {
        if (!std::isfinite(maxVelocity) || maxVelocity <= 0.0) {
            return std::unexpected(SamplingError::InvalidMaxVelocity);
        }

        if (!detail::finite(line.from) || !detail::finite(line.to)) {
            return std::unexpected(SamplingError::NonFiniteGeometry);
        }

        Vector<DoF> delta {};

        for (std::size_t axis = 0; axis < DoF; ++axis) {
            delta[axis] = line.to[axis] - line.from[axis];
        }

        const auto length = detail::magnitude(delta);

        if (!std::isfinite(length) || length <= 0.0) {
            return std::unexpected(SamplingError::InvalidLength);
        }

        Vector<DoF> tangent {};

        for (std::size_t axis = 0; axis < DoF; ++axis) {
            tangent[axis] = delta[axis] / length;
        }

        SampledPathPiece<DoF> result;
        result.id = id;
        result.length = length;
        result.maxVelocity = maxVelocity;
        result.stations = {
            {.distance = 0.0, .tangent = tangent},
            {.distance = length, .tangent = tangent},
        };

        return result;
    }
}
