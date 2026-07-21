#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <limits>
#include <numeric>
#include <span>
#include <utility>
#include <vector>

#include "path_tempo/Types.h"

namespace path_tempo {
    enum class SamplingError {
        InvalidLength,
        InvalidMaxVelocity,
        InvalidSampleCount,
        InvalidKnots,
        InvalidWeights,
        InvalidPieceIds,
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

    template<std::size_t DoF>
    struct Arc {
        // At parameter zero the curve is at origin. Radial and tangential are
        // equal-length orthogonal basis vectors; axial is the total helical offset.
        Vector<DoF> origin {};
        Vector<DoF> radial {};
        Vector<DoF> tangential {};
        Vector<DoF> axial {};
        double sweep = 0.0;
    };

    template<std::size_t DoF>
    struct BSpline {
        std::size_t degree = 0;
        std::span<const Vector<DoF>> controls;
        std::span<const double> knots;
    };

    template<std::size_t DoF>
    struct Nurbs {
        std::size_t degree = 0;
        std::span<const Vector<DoF>> controls;
        std::span<const double> weights;
        std::span<const double> knots;
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

        // Derivative magnitudes below this value cannot be normalized reliably
        // in double precision and indicate a singular curve parameterization.
        inline constexpr double MINIMUM_CURVE_SPEED = 1e-15;
        // Arc frames normally come from normalized construction geometry. This
        // relative tolerance accepts its roundoff without treating an ellipse as an arc.
        inline constexpr double ARC_FRAME_RELATIVE_TOLERANCE = 1e-10;
        // Adaptive Simpson integration at this relative scale keeps station
        // distances well below the planner's geometric continuity tolerances.
        inline constexpr double ARC_LENGTH_RELATIVE_TOLERANCE = 1e-12;
        // Twenty subdivisions provide a hard recursion bound while permitting
        // over a million leaf intervals for unusually distorted parameterizations.
        inline constexpr unsigned ARC_LENGTH_MAXIMUM_DEPTH = 20;

        template<std::size_t Size>
        double dot(const std::array<double, Size> &left, const std::array<double, Size> &right) {
            return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
        }

        template<std::size_t Size>
        std::array<double, Size> scaled(const std::array<double, Size> &value, const double amount) {
            std::array<double, Size> result {};

            for (std::size_t component = 0; component < Size; ++component) {
                result[component] = value[component] * amount;
            }

            return result;
        }

        template<std::size_t Size>
        std::array<double, Size> add(const std::array<double, Size> &left, const std::array<double, Size> &right) {
            std::array<double, Size> result {};

            for (std::size_t component = 0; component < Size; ++component) {
                result[component] = left[component] + right[component];
            }

            return result;
        }

        template<std::size_t Size>
        std::array<double, Size> subtract(const std::array<double, Size> &left, const std::array<double, Size> &right) {
            std::array<double, Size> result {};

            for (std::size_t component = 0; component < Size; ++component) {
                result[component] = left[component] - right[component];
            }

            return result;
        }

        template<typename Function>
        double adaptiveSimpson(const Function &function, const double from, const double to, const double tolerance, const double fromValue, const double middleValue, const double toValue, const double estimate, const unsigned depth) {
            const auto middle = std::midpoint(from, to);
            const auto leftMiddle = std::midpoint(from, middle);
            const auto rightMiddle = std::midpoint(middle, to);
            const auto leftMiddleValue = function(leftMiddle);
            const auto rightMiddleValue = function(rightMiddle);
            const auto left = (middle - from) * (fromValue + 4.0 * leftMiddleValue + middleValue) / 6.0;
            const auto right = (to - middle) * (middleValue + 4.0 * rightMiddleValue + toValue) / 6.0;
            const auto refined = left + right;

            if (depth == ARC_LENGTH_MAXIMUM_DEPTH || std::abs(refined - estimate) <= 15.0 * tolerance) {
                return refined + (refined - estimate) / 15.0;
            }

            return adaptiveSimpson(function, from, middle, tolerance * 0.5, fromValue, leftMiddleValue, middleValue, left, depth + 1)
                + adaptiveSimpson(function, middle, to, tolerance * 0.5, middleValue, rightMiddleValue, toValue, right, depth + 1);
        }

        template<typename Function>
        double integrate(const Function &function, const double from, const double to) {
            if (to <= from) {
                return 0.0;
            }

            const auto middle = std::midpoint(from, to);
            const auto fromValue = function(from);
            const auto middleValue = function(middle);
            const auto toValue = function(to);
            const auto estimate = (to - from) * (fromValue + 4.0 * middleValue + toValue) / 6.0;
            const auto tolerance = ARC_LENGTH_RELATIVE_TOLERANCE * std::max(1.0, std::abs(estimate));

            return adaptiveSimpson(function, from, to, tolerance, fromValue, middleValue, toValue, estimate, 0);
        }

        template<std::size_t DoF>
        std::expected<DifferentialState<DoF>, SamplingError> arcLengthState(const Vector<DoF> &first, const Vector<DoF> &second, const Vector<DoF> &third) {
            if (!finite(first) || !finite(second) || !finite(third)) {
                return std::unexpected(SamplingError::NonFiniteGeometry);
            }

            const auto speed = magnitude(first);

            if (!std::isfinite(speed) || speed <= MINIMUM_CURVE_SPEED) {
                return std::unexpected(SamplingError::IrregularCurve);
            }

            const auto firstSecond = dot(first, second);
            const auto secondSquared = dot(second, second);
            const auto firstThird = dot(first, third);
            const auto inverseSpeed = 1.0 / speed;
            const auto inverseSpeed2 = inverseSpeed * inverseSpeed;
            const auto inverseSpeed4 = inverseSpeed2 * inverseSpeed2;
            const auto inverseSpeed6 = inverseSpeed4 * inverseSpeed2;
            const auto tangent = scaled(first, inverseSpeed);
            const auto curvature = subtract(scaled(second, inverseSpeed2), scaled(first, firstSecond * inverseSpeed4));
            auto parameterDerivative = scaled(third, inverseSpeed2);
            parameterDerivative = add(parameterDerivative, scaled(second, -3.0 * firstSecond * inverseSpeed4));
            parameterDerivative = add(parameterDerivative, scaled(first, -(secondSquared + firstThird) * inverseSpeed4));
            parameterDerivative = add(parameterDerivative, scaled(first, 4.0 * firstSecond * firstSecond * inverseSpeed6));

            return DifferentialState<DoF> {
                .tangent = tangent,
                .curvature = curvature,
                .thirdDerivative = scaled(parameterDerivative, inverseSpeed),
            };
        }

        template<std::size_t Size>
        std::array<double, Size> evaluateBSpline(const std::size_t degree, const std::span<const std::array<double, Size>> controls, const std::span<const double> knots, const double parameter, const std::size_t span) {
            std::vector<std::array<double, Size>> work(degree + 1);

            for (std::size_t index = 0; index <= degree; ++index) {
                work[index] = controls[span - degree + index];
            }

            for (std::size_t level = 1; level <= degree; ++level) {
                for (std::size_t index = degree; index >= level; --index) {
                    const auto knot = span - degree + index;
                    const auto denominator = knots[knot + degree - level + 1] - knots[knot];
                    const auto alpha = denominator > 0.0 ? (parameter - knots[knot]) / denominator : 0.0;
                    work[index] = add(scaled(work[index - 1], 1.0 - alpha), scaled(work[index], alpha));
                }
            }

            return work[degree];
        }

        template<std::size_t Size>
        struct DerivativeSpline {
            std::size_t degree = 0;
            std::vector<std::array<double, Size>> controls;
            std::vector<double> knots;
        };

        template<std::size_t Size>
        DerivativeSpline<Size> derivativeSpline(const DerivativeSpline<Size> &source) {
            DerivativeSpline<Size> result;

            if (source.degree == 0 || source.controls.size() < 2) {
                return result;
            }

            result.degree = source.degree - 1;
            result.controls.reserve(source.controls.size() - 1);
            result.knots.assign(source.knots.begin() + 1, source.knots.end() - 1);

            for (std::size_t index = 0; index + 1 < source.controls.size(); ++index) {
                const auto denominator = source.knots[index + source.degree + 1] - source.knots[index + 1];
                result.controls.push_back(denominator > 0.0
                    ? scaled(subtract(source.controls[index + 1], source.controls[index]), static_cast<double>(source.degree) / denominator)
                    : std::array<double, Size> {});
            }

            return result;
        }

        template<std::size_t DoF>
        bool validSpline(const std::size_t degree, const std::span<const Vector<DoF>> controls, const std::span<const double> knots) {
            if (degree == 0 || controls.size() <= degree || knots.size() != controls.size() + degree + 1) {
                return false;
            }

            if (std::ranges::any_of(controls, [](const auto &control) { return !finite(control); })
                || std::ranges::any_of(knots, [](const double knot) { return !std::isfinite(knot); })) {
                return false;
            }

            if (!std::ranges::is_sorted(knots)) {
                return false;
            }

            return knots[degree] < knots[controls.size()];
        }

        template<std::size_t DoF>
        struct RationalEvaluator {
            static constexpr std::size_t HOMOGENEOUS_SIZE = DoF + 1;
            std::array<DerivativeSpline<HOMOGENEOUS_SIZE>, 4> derivatives;

            std::expected<std::array<Vector<DoF>, 4>, SamplingError> operator()(const double parameter, const std::size_t span) const {
                std::array<std::array<double, HOMOGENEOUS_SIZE>, 4> homogeneous {};

                for (std::size_t order = 0; order < derivatives.size(); ++order) {
                    const auto &derivative = derivatives[order];

                    if (derivative.controls.empty()) {
                        continue;
                    }

                    homogeneous[order] = evaluateBSpline(derivative.degree, std::span<const std::array<double, HOMOGENEOUS_SIZE>> {derivative.controls}, derivative.knots, parameter, span - order);
                }

                const auto weight = homogeneous[0][DoF];

                if (!std::isfinite(weight) || weight <= MINIMUM_CURVE_SPEED) {
                    return std::unexpected(SamplingError::IrregularCurve);
                }

                std::array<Vector<DoF>, 4> result {};

                for (std::size_t component = 0; component < DoF; ++component) {
                    result[0][component] = homogeneous[0][component] / weight;
                    result[1][component] = (homogeneous[1][component] - homogeneous[1][DoF] * result[0][component]) / weight;
                    result[2][component] = (homogeneous[2][component] - 2.0 * homogeneous[1][DoF] * result[1][component] - homogeneous[2][DoF] * result[0][component]) / weight;
                    result[3][component] = (homogeneous[3][component] - 3.0 * homogeneous[1][DoF] * result[2][component] - 3.0 * homogeneous[2][DoF] * result[1][component] - homogeneous[3][DoF] * result[0][component]) / weight;
                }

                if (std::ranges::any_of(result, [](const auto &value) { return !finite(value); })) {
                    return std::unexpected(SamplingError::NonFiniteGeometry);
                }

                return result;
            }
        };

        template<std::size_t DoF>
        std::expected<RationalEvaluator<DoF>, SamplingError> rationalEvaluator(const std::size_t degree, const std::span<const Vector<DoF>> controls, const std::span<const double> weights, const std::span<const double> knots) {
            if (!validSpline(degree, controls, knots)) {
                return std::unexpected(SamplingError::InvalidKnots);
            }

            if (!weights.empty() && (weights.size() != controls.size() || std::ranges::any_of(weights, [](const double weight) { return !std::isfinite(weight) || weight <= 0.0; }))) {
                return std::unexpected(SamplingError::InvalidWeights);
            }

            RationalEvaluator<DoF> result;
            auto &spline = result.derivatives[0];
            spline.degree = degree;
            spline.knots.assign(knots.begin(), knots.end());
            spline.controls.resize(controls.size());

            for (std::size_t control = 0; control < controls.size(); ++control) {
                const auto weight = weights.empty() ? 1.0 : weights[control];

                for (std::size_t component = 0; component < DoF; ++component) {
                    spline.controls[control][component] = controls[control][component] * weight;
                }

                spline.controls[control][DoF] = weight;
            }

            for (std::size_t order = 1; order < result.derivatives.size(); ++order) {
                result.derivatives[order] = derivativeSpline(result.derivatives[order - 1]);
            }

            return result;
        }

        template<std::size_t DoF, typename Evaluator>
        std::expected<SampledPathPiece<DoF>, SamplingError> sampleParametricInterval(const Evaluator &evaluate, const double from, const double to, const std::size_t span, const PathPieceId id, const double maxVelocity, const std::size_t sampleIntervals) {
            if (!std::isfinite(maxVelocity) || maxVelocity <= 0.0) {
                return std::unexpected(SamplingError::InvalidMaxVelocity);
            }

            if (sampleIntervals == 0) {
                return std::unexpected(SamplingError::InvalidSampleCount);
            }

            const auto speed = [&](const double parameter) {
                const auto derivatives = evaluate(parameter, span);
                return derivatives ? magnitude((*derivatives)[1]) : std::numeric_limits<double>::quiet_NaN();
            };
            SampledPathPiece<DoF> result;
            result.id = id;
            result.maxVelocity = maxVelocity;
            result.stations.reserve(sampleIntervals + 1);
            auto distance = 0.0;

            for (std::size_t sample = 0; sample <= sampleIntervals; ++sample) {
                const auto parameter = std::lerp(from, to, static_cast<double>(sample) / static_cast<double>(sampleIntervals));
                const auto derivatives = evaluate(parameter, span);

                if (!derivatives) {
                    return std::unexpected(derivatives.error());
                }

                if (sample > 0) {
                    const auto previous = std::lerp(from, to, static_cast<double>(sample - 1) / static_cast<double>(sampleIntervals));
                    distance += integrate(speed, previous, parameter);
                }

                const auto state = arcLengthState((*derivatives)[1], (*derivatives)[2], (*derivatives)[3]);

                if (!state || !std::isfinite(distance)) {
                    return std::unexpected(state ? SamplingError::IrregularCurve : state.error());
                }

                result.stations.push_back({
                    .distance = distance,
                    .tangent = state->tangent,
                    .curvature = state->curvature,
                    .thirdDerivative = state->thirdDerivative,
                });
            }

            result.length = distance;

            if (!std::isfinite(result.length) || result.length <= 0.0) {
                return std::unexpected(SamplingError::InvalidLength);
            }

            result.stations.back().distance = result.length;
            return result;
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

    template<std::size_t DoF, typename Evaluator>
    std::expected<SampledPathPiece<DoF>, SamplingError> sampleArcLengthCurve(const double length, const PathPieceId id, const double maxVelocity, const std::size_t sampleIntervals, const Evaluator &evaluate) {
        if (!std::isfinite(length) || length <= 0.0) {
            return std::unexpected(SamplingError::InvalidLength);
        }

        if (!std::isfinite(maxVelocity) || maxVelocity <= 0.0) {
            return std::unexpected(SamplingError::InvalidMaxVelocity);
        }

        if (sampleIntervals == 0) {
            return std::unexpected(SamplingError::InvalidSampleCount);
        }

        SampledPathPiece<DoF> result;
        result.id = id;
        result.length = length;
        result.maxVelocity = maxVelocity;
        result.stations.reserve(sampleIntervals + 1);

        for (std::size_t sample = 0; sample <= sampleIntervals; ++sample) {
            const auto distance = sample == sampleIntervals
                ? length
                : length * static_cast<double>(sample) / static_cast<double>(sampleIntervals);
            const auto state = evaluate(distance);

            if (!detail::finite(state.tangent) || !detail::finite(state.curvature) || !detail::finite(state.thirdDerivative)) {
                return std::unexpected(SamplingError::NonFiniteGeometry);
            }

            result.stations.push_back({
                .distance = distance,
                .tangent = state.tangent,
                .curvature = state.curvature,
                .thirdDerivative = state.thirdDerivative,
            });
        }

        return result;
    }

    template<std::size_t DoF>
    std::expected<SampledPathPiece<DoF>, SamplingError> sampleArc(const Arc<DoF> &arc, const PathPieceId id, const double maxVelocity, const std::size_t sampleIntervals) {
        if (!detail::finite(arc.origin) || !detail::finite(arc.radial) || !detail::finite(arc.tangential) || !detail::finite(arc.axial) || !std::isfinite(arc.sweep)) {
            return std::unexpected(SamplingError::NonFiniteGeometry);
        }

        if (arc.sweep == 0.0) {
            return std::unexpected(SamplingError::InvalidLength);
        }

        const auto radialSquared = detail::dot(arc.radial, arc.radial);
        const auto tangentialSquared = detail::dot(arc.tangential, arc.tangential);
        const auto axialSquared = detail::dot(arc.axial, arc.axial);
        const auto frameScale = std::max(radialSquared, tangentialSquared);
        const auto axialFrameScale = std::sqrt(frameScale * axialSquared);

        if (frameScale <= detail::MINIMUM_CURVE_SPEED * detail::MINIMUM_CURVE_SPEED
            || std::abs(radialSquared - tangentialSquared) > detail::ARC_FRAME_RELATIVE_TOLERANCE * frameScale
            || std::abs(detail::dot(arc.radial, arc.tangential)) > detail::ARC_FRAME_RELATIVE_TOLERANCE * frameScale
            || std::abs(detail::dot(arc.radial, arc.axial)) > detail::ARC_FRAME_RELATIVE_TOLERANCE * axialFrameScale
            || std::abs(detail::dot(arc.tangential, arc.axial)) > detail::ARC_FRAME_RELATIVE_TOLERANCE * axialFrameScale) {
            return std::unexpected(SamplingError::IrregularCurve);
        }

        const auto evaluate = [&](const double parameter, const std::size_t) -> std::expected<std::array<Vector<DoF>, 4>, SamplingError> {
            const auto angle = arc.sweep * parameter;
            const auto cosine = std::cos(angle);
            const auto sine = std::sin(angle);
            std::array<Vector<DoF>, 4> result {};

            for (std::size_t component = 0; component < DoF; ++component) {
                result[0][component] = arc.origin[component] + arc.radial[component] * (cosine - 1.0) + arc.tangential[component] * sine + arc.axial[component] * parameter;
                result[1][component] = arc.sweep * (-arc.radial[component] * sine + arc.tangential[component] * cosine) + arc.axial[component];
                result[2][component] = -arc.sweep * arc.sweep * (arc.radial[component] * cosine + arc.tangential[component] * sine);
                result[3][component] = arc.sweep * arc.sweep * arc.sweep * (arc.radial[component] * sine - arc.tangential[component] * cosine);
            }

            return result;
        };

        return detail::sampleParametricInterval<DoF>(evaluate, 0.0, 1.0, 0, id, maxVelocity, sampleIntervals);
    }

    template<std::size_t DoF>
    std::expected<std::vector<SampledPathPiece<DoF>>, SamplingError> sampleNurbs(const Nurbs<DoF> &spline, const PathPieceId firstId, const double maxVelocity, const std::size_t samplesPerInterval) {
        const auto evaluator = detail::rationalEvaluator(spline.degree, spline.controls, spline.weights, spline.knots);

        if (!evaluator) {
            return std::unexpected(evaluator.error());
        }

        std::vector<SampledPathPiece<DoF>> result;

        for (std::size_t span = spline.degree; span < spline.controls.size(); ++span) {
            if (spline.knots[span + 1] <= spline.knots[span]) {
                continue;
            }

            if (result.size() > std::numeric_limits<PathPieceId>::max() - firstId) {
                return std::unexpected(SamplingError::InvalidPieceIds);
            }

            auto piece = detail::sampleParametricInterval<DoF>(*evaluator, spline.knots[span], spline.knots[span + 1], span, firstId + static_cast<PathPieceId>(result.size()), maxVelocity, samplesPerInterval);

            if (!piece) {
                return std::unexpected(piece.error());
            }

            result.push_back(std::move(*piece));
        }

        if (result.empty()) {
            return std::unexpected(SamplingError::InvalidKnots);
        }

        return result;
    }

    template<std::size_t DoF>
    std::expected<std::vector<SampledPathPiece<DoF>>, SamplingError> sampleBSpline(const BSpline<DoF> &spline, const PathPieceId firstId, const double maxVelocity, const std::size_t samplesPerInterval) {
        return sampleNurbs(Nurbs<DoF> {
            .degree = spline.degree,
            .controls = spline.controls,
            .weights = {},
            .knots = spline.knots,
        }, firstId, maxVelocity, samplesPerInterval);
    }
}
