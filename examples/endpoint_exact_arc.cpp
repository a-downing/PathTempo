#include "endpoint_exact_arc.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <numbers>

namespace path_tempo::example {
    namespace {
        constexpr double ARC_CLOSURE_ANGULAR_TOLERANCE = 1e-12;
        constexpr double ARC_CLOSURE_RADIAL_RELATIVE_TOLERANCE = 1e-12;

        double dot(const Vector3 &left, const Vector3 &right) {
            return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
        }

        Vector3 add(const Vector3 &left, const Vector3 &right) {
            return {left[0] + right[0], left[1] + right[1], left[2] + right[2]};
        }

        Vector3 subtract(const Vector3 &left, const Vector3 &right) {
            return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
        }

        Vector3 scale(const Vector3 &value, const double amount) {
            return {value[0] * amount, value[1] * amount, value[2] * amount};
        }

        Vector3 cross(const Vector3 &left, const Vector3 &right) {
            return {
                left[1] * right[2] - left[2] * right[1],
                left[2] * right[0] - left[0] * right[2],
                left[0] * right[1] - left[1] * right[0],
            };
        }

        double magnitude(const Vector3 &value) {
            return std::sqrt(dot(value, value));
        }

        Vector3 rotate(const Vector3 &value, const double angle, const Vector3 &axis) {
            return add(add(scale(value, std::cos(angle)), scale(cross(axis, value), std::sin(angle))),
                       scale(axis, dot(axis, value) * (1.0 - std::cos(angle))));
        }

        template<typename Function>
        double adaptiveSimpson(const Function &function, const double from, const double to, const double tolerance,
                               const double fromValue, const double middleValue, const double toValue,
                               const double estimate, const unsigned depth) {
            const auto middle = std::midpoint(from, to);
            const auto leftMiddle = std::midpoint(from, middle);
            const auto rightMiddle = std::midpoint(middle, to);
            const auto leftMiddleValue = function(leftMiddle);
            const auto rightMiddleValue = function(rightMiddle);
            const auto left = (middle - from) * (fromValue + 4.0 * leftMiddleValue + middleValue) / 6.0;
            const auto right = (to - middle) * (middleValue + 4.0 * rightMiddleValue + toValue) / 6.0;
            const auto refined = left + right;

            if (depth == 20 || std::abs(refined - estimate) <= 15.0 * tolerance) {
                return refined + (refined - estimate) / 15.0;
            }

            return adaptiveSimpson(function, from, middle, tolerance * 0.5, fromValue, leftMiddleValue,
                                   middleValue, left, depth + 1)
                + adaptiveSimpson(function, middle, to, tolerance * 0.5, middleValue, rightMiddleValue,
                                  toValue, right, depth + 1);
        }
    }

    std::array<Vector<3>, 4> EndpointExactArc::derivatives(const double requestedParameter) const {
        const auto parameter = std::clamp(requestedParameter, 0.0, 1.0);
        const auto start = rotate(m_startArm, m_sweep * parameter, m_axis);
        const auto end = rotate(m_endArm, -m_sweep * (1.0 - parameter), m_axis);
        const auto startFirst = scale(cross(m_axis, start), m_sweep);
        const auto endFirst = scale(cross(m_axis, end), m_sweep);
        const auto startSecond = scale(cross(m_axis, startFirst), m_sweep);
        const auto endSecond = scale(cross(m_axis, endFirst), m_sweep);
        const auto startThird = scale(cross(m_axis, startSecond), m_sweep);
        const auto endThird = scale(cross(m_axis, endSecond), m_sweep);
        std::array<Vector<3>, 4> result {};

        const auto radial = add(scale(start, 1.0 - parameter), scale(end, parameter));
        const auto xyz = add(add(add(m_source.center, m_startAxial), radial), scale(m_axialTravel, parameter));

        for (std::size_t component = 0; component < result[0].size(); ++component) {
            result[0][component] = xyz[component];
        }

        const auto radialFirst = add(add(scale(start, -1.0), scale(startFirst, 1.0 - parameter)),
                                     add(end, scale(endFirst, parameter)));
        const auto first = add(radialFirst, m_axialTravel);
        const auto second = add(add(scale(startFirst, -2.0), scale(startSecond, 1.0 - parameter)),
                                add(scale(endFirst, 2.0), scale(endSecond, parameter)));
        const auto third = add(add(scale(startSecond, -3.0), scale(startThird, 1.0 - parameter)),
                               add(scale(endSecond, 3.0), scale(endThird, parameter)));

        for (std::size_t component = 0; component < result[0].size(); ++component) {
            result[1][component] = first[component];
            result[2][component] = second[component];
            result[3][component] = third[component];
        }

        return result;
    }

    double EndpointExactArc::speed(const double parameter) const {
        return magnitude(derivatives(parameter)[1]);
    }

    double EndpointExactArc::integratedLength(const double from, const double to) const {
        if (to <= from) {
            return 0.0;
        }

        const auto function = [this](const double parameter) { return speed(parameter); };
        const auto middle = std::midpoint(from, to);
        const auto fromValue = function(from);
        const auto middleValue = function(middle);
        const auto toValue = function(to);
        const auto estimate = (to - from) * (fromValue + 4.0 * middleValue + toValue) / 6.0;

        return adaptiveSimpson(function, from, to, 1e-13 * std::max(1.0, estimate), fromValue, middleValue, toValue, estimate, 0);
    }

    double EndpointExactArc::parameterAtDistance(const double requestedDistance) const {
        const auto distance = std::clamp(requestedDistance, 0.0, m_length);

        if (distance == 0.0) {
            return 0.0;
        }
        if (distance == m_length) {
            return 1.0;
        }

        const auto upper = std::lower_bound(m_lengthNodes.begin(), m_lengthNodes.end(), distance,
            [](const LengthNode &node, const double value) { return node.distance < value; });
        const auto lower = upper - 1;
        auto from = lower->parameter;
        auto to = upper->parameter;
        auto parameter = std::lerp(from, to, (distance - lower->distance) / (upper->distance - lower->distance));
        const auto tolerance = 1e-12 * std::max(1.0, m_length);

        for (unsigned iteration = 0; iteration < 16; ++iteration) {
            const auto currentDistance = lower->distance + integratedLength(lower->parameter, parameter);
            const auto error = currentDistance - distance;

            if (std::abs(error) <= tolerance) {
                return parameter;
            }

            if (error < 0.0) {
                from = parameter;
            } else {
                to = parameter;
            }

            const auto currentSpeed = speed(parameter);
            const auto newton = currentSpeed > 1e-15 ? parameter - error / currentSpeed : parameter;
            parameter = newton > from && newton < to ? newton : std::midpoint(from, to);
        }

        return parameter;
    }

    std::expected<EndpointExactArc, std::string> EndpointExactArc::create(const ArcGeometry &source) {
        EndpointExactArc result;
        result.m_source = source;

        const auto axisLength = magnitude(source.axis);

        if (!std::isfinite(axisLength) || axisLength <= 0.0) {
            return std::unexpected("arc axis must be finite and nonzero");
        }

        result.m_axis = scale(source.axis, 1.0 / axisLength);
        const Vector3 start {source.from[0], source.from[1], source.from[2]};
        const Vector3 end {source.to[0], source.to[1], source.to[2]};
        const auto startDelta = subtract(start, source.center);
        const auto endDelta = subtract(end, source.center);
        result.m_startAxial = scale(result.m_axis, dot(startDelta, result.m_axis));
        result.m_startArm = subtract(startDelta, result.m_startAxial);
        result.m_endArm = subtract(endDelta, scale(result.m_axis, dot(endDelta, result.m_axis)));
        const auto startRadius = magnitude(result.m_startArm);
        const auto endRadius = magnitude(result.m_endArm);

        if (!std::isfinite(startRadius) || !std::isfinite(endRadius) || startRadius <= 0.0 || endRadius <= 0.0) {
            return std::unexpected("arc radius must be finite and nonzero");
        }

        result.m_axialTravel = scale(result.m_axis, dot(subtract(end, start), result.m_axis));
        const auto startUnit = scale(result.m_startArm, 1.0 / startRadius);
        const auto endUnit = scale(result.m_endArm, 1.0 / endRadius);
        result.m_sweep = std::atan2(dot(result.m_axis, cross(startUnit, endUnit)), dot(startUnit, endUnit));

        if (result.m_sweep < 0.0) {
            result.m_sweep += 2.0 * std::numbers::pi;
        }
        const auto angularSeparation = magnitude(subtract(startUnit, endUnit));
        const auto radialTolerance = ARC_CLOSURE_RADIAL_RELATIVE_TOLERANCE * std::max(startRadius, endRadius);

        if (angularSeparation <= ARC_CLOSURE_ANGULAR_TOLERANCE && std::abs(startRadius - endRadius) <= radialTolerance) {
            result.m_sweep = 2.0 * std::numbers::pi;
        }
        if (result.m_sweep <= 0.0) {
            return std::unexpected("arc sweep must be positive");
        }

        result.m_lengthNodes[0] = {};
        auto cumulative = 0.0;

        for (std::size_t index = 1; index <= LENGTH_TABLE_INTERVALS; ++index) {
            const auto parameter = static_cast<double>(index) / LENGTH_TABLE_INTERVALS;
            const auto previous = static_cast<double>(index - 1) / LENGTH_TABLE_INTERVALS;

            cumulative += result.integratedLength(previous, parameter);
            result.m_lengthNodes[index] = {parameter, cumulative};
        }

        result.m_length = cumulative;

        return result;
    }

    double EndpointExactArc::length() const {
        return m_length;
    }

    Vector3 EndpointExactArc::positionAtParameter(const double parameter) const {
        return derivatives(parameter)[0];
    }

    DifferentialState<3> EndpointExactArc::stateAtDistance(const double distance) const {
        const auto values = derivatives(parameterAtDistance(distance));
        const auto &first = values[1];
        const auto &second = values[2];
        const auto &third = values[3];
        const auto speedValue = magnitude(first);
        const auto firstSecond = dot(first, second);
        const auto secondSquared = dot(second, second);
        const auto firstThird = dot(first, third);
        const auto speed2 = speedValue * speedValue;
        const auto speed3 = speed2 * speedValue;
        const auto speed4 = speed2 * speed2;
        const auto speed5 = speed4 * speedValue;
        const auto speed7 = speed5 * speed2;
        const auto tangent = scale(first, 1.0 / speedValue);
        const auto curvature = add(scale(second, 1.0 / speed2), scale(first, -firstSecond / speed4));
        auto thirdDerivative = scale(third, 1.0 / speed3);
        thirdDerivative = add(thirdDerivative, scale(second, -3.0 * firstSecond / speed5));
        thirdDerivative = add(thirdDerivative,
            scale(first, -(secondSquared + firstThird) / speed5 + 4.0 * firstSecond * firstSecond / speed7));

        return {tangent, curvature, thirdDerivative};
    }
}
