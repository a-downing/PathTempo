#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>

#include "path_tempo/Types.h"

namespace path_tempo {
    struct UniformArcLengthSampling {
        std::size_t intervals = 64;
    };

    struct PerKnotIntervalSampling {
        std::size_t intervals = 16;
    };

    enum class SamplingError {
        InvalidLength,
        InvalidProgrammedVelocity,
        InvalidIntervalCount,
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
    std::expected<SampledPathPiece<DoF>, SamplingError> sampleLine(const Line<DoF> &line, const PathPieceId id, const double programmedVelocity, const UniformArcLengthSampling sampling = {.intervals = 1}) {
        if (sampling.intervals == 0) {
            return std::unexpected(SamplingError::InvalidIntervalCount);
        }

        if (!std::isfinite(programmedVelocity) || programmedVelocity <= 0.0) {
            return std::unexpected(SamplingError::InvalidProgrammedVelocity);
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
        result.programmedVelocity = programmedVelocity;
        result.stations.reserve(sampling.intervals + 1);

        for (std::size_t sample = 0; sample <= sampling.intervals; ++sample) {
            result.stations.push_back({
                .distance = length * static_cast<double>(sample) / static_cast<double>(sampling.intervals),
                .tangent = tangent,
            });
        }

        return result;
    }
}
