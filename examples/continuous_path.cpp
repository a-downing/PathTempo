#include "geometry_text_loader.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdio>
#include <expected>
#include <filesystem>
#include <format>
#include <iterator>
#include <numbers>
#include <numeric>
#include <print>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "path_tempo/Planner.h"
#include "path_tempo/Sampling.h"

namespace {
    constexpr std::size_t DOF = 3;
    constexpr std::size_t CURVE_SAMPLE_INTERVALS = 16;
    constexpr std::size_t ARC_LENGTH_TABLE_INTERVALS = 64;
    using Position = path_tempo::Vector<DOF>;
    using Vector3 = path_tempo::example::Vector3;

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

    double positionDot(const Position &left, const Position &right) {
        return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
    }

    double positionMagnitude(const Position &value) {
        return std::sqrt(positionDot(value, value));
    }

    Position positionScale(const Position &value, const double amount) {
        Position result {};
        for (std::size_t component = 0; component < DOF; ++component) {
            result[component] = value[component] * amount;
        }

        return result;
    }

    Position positionAdd(const Position &left, const Position &right) {
        Position result {};
        for (std::size_t component = 0; component < DOF; ++component) {
            result[component] = left[component] + right[component];
        }

        return result;
    }

    class EndpointExactArc {
        struct LengthNode {
            double parameter = 0.0;
            double distance = 0.0;
        };

        path_tempo::example::ArcGeometry m_source;
        Vector3 m_axis {};
        Vector3 m_startArm {};
        Vector3 m_endArm {};
        Vector3 m_axial {};
        double m_sweep = 0.0;
        double m_length = 0.0;
        std::array<LengthNode, ARC_LENGTH_TABLE_INTERVALS + 1> m_lengthNodes {};

        std::array<Position, 4> derivatives(const double requestedParameter) const {
            const auto parameter = std::clamp(requestedParameter, 0.0, 1.0);
            const auto start = rotate(m_startArm, m_sweep * parameter, m_axis);
            const auto end = rotate(m_endArm, -m_sweep * (1.0 - parameter), m_axis);
            const auto startFirst = scale(cross(m_axis, start), m_sweep);
            const auto endFirst = scale(cross(m_axis, end), m_sweep);
            const auto startSecond = scale(cross(m_axis, startFirst), m_sweep);
            const auto endSecond = scale(cross(m_axis, endFirst), m_sweep);
            const auto startThird = scale(cross(m_axis, startSecond), m_sweep);
            const auto endThird = scale(cross(m_axis, endSecond), m_sweep);
            std::array<Position, 4> result {};

            const auto radial = add(scale(start, 1.0 - parameter), scale(end, parameter));
            const auto xyz = add(add(m_source.center, radial), scale(m_axial, parameter));
            for (std::size_t component = 0; component < 3; ++component) {
                result[0][component] = xyz[component];
            }
            const auto radialFirst = add(add(scale(start, -1.0), scale(startFirst, 1.0 - parameter)),
                                         add(end, scale(endFirst, parameter)));
            const auto first = add(radialFirst, m_axial);
            const auto second = add(add(scale(startFirst, -2.0), scale(startSecond, 1.0 - parameter)),
                                    add(scale(endFirst, 2.0), scale(endSecond, parameter)));
            const auto third = add(add(scale(startSecond, -3.0), scale(startThird, 1.0 - parameter)),
                                   add(scale(endSecond, 3.0), scale(endThird, parameter)));
            for (std::size_t component = 0; component < 3; ++component) {
                result[1][component] = first[component];
                result[2][component] = second[component];
                result[3][component] = third[component];
            }
            return result;
        }

        double speed(const double parameter) const {
            return positionMagnitude(derivatives(parameter)[1]);
        }

        template<typename Function>
        static double adaptiveSimpson(const Function &function, const double from, const double to, const double tolerance,
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

        double integratedLength(const double from, const double to) const {
            if (to <= from) {
                return 0.0;
            }

            const auto function = [this](const double parameter) { return speed(parameter); };
            const auto middle = std::midpoint(from, to);
            const auto fromValue = function(from);
            const auto middleValue = function(middle);
            const auto toValue = function(to);
            const auto estimate = (to - from) * (fromValue + 4.0 * middleValue + toValue) / 6.0;

            return adaptiveSimpson(function, from, to, 1e-13 * std::max(1.0, estimate),
                                   fromValue, middleValue, toValue, estimate, 0);
        }

        double parameterAtDistance(const double requestedDistance) const {
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

    public:
        static std::expected<EndpointExactArc, std::string> create(const path_tempo::example::ArcGeometry &source) {
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
            result.m_startArm = subtract(startDelta, scale(result.m_axis, dot(startDelta, result.m_axis)));
            result.m_endArm = subtract(endDelta, scale(result.m_axis, dot(endDelta, result.m_axis)));
            const auto startRadius = magnitude(result.m_startArm);
            const auto endRadius = magnitude(result.m_endArm);
            if (!std::isfinite(startRadius) || !std::isfinite(endRadius) || startRadius <= 0.0 || endRadius <= 0.0) {
                return std::unexpected("arc radius must be finite and nonzero");
            }

            result.m_axial = scale(result.m_axis, dot(subtract(end, start), result.m_axis));
            const auto startUnit = scale(result.m_startArm, 1.0 / startRadius);
            const auto endUnit = scale(result.m_endArm, 1.0 / endRadius);
            result.m_sweep = std::atan2(dot(result.m_axis, cross(startUnit, endUnit)), dot(startUnit, endUnit));
            if (result.m_sweep < 0.0) {
                result.m_sweep += 2.0 * std::numbers::pi;
            }

            if (magnitude(subtract(result.m_startArm, result.m_endArm)) < 1e-9) {
                result.m_sweep = 2.0 * std::numbers::pi;
            }

            if (result.m_sweep <= 0.0) {
                return std::unexpected("arc sweep must be positive");
            }

            result.m_lengthNodes[0] = {};
            auto cumulative = 0.0;
            for (std::size_t index = 1; index <= ARC_LENGTH_TABLE_INTERVALS; ++index) {
                const auto parameter = static_cast<double>(index) / ARC_LENGTH_TABLE_INTERVALS;
                const auto previous = static_cast<double>(index - 1) / ARC_LENGTH_TABLE_INTERVALS;

                cumulative += result.integratedLength(previous, parameter);
                result.m_lengthNodes[index] = {parameter, cumulative};
            }

            result.m_length = cumulative;

            return result;
        }

        double length() const { return m_length; }

        path_tempo::DifferentialState<DOF> stateAtDistance(const double distance) const {
            const auto values = derivatives(parameterAtDistance(distance));
            const auto &first = values[1];
            const auto &second = values[2];
            const auto &third = values[3];
            const auto speedValue = positionMagnitude(first);
            const auto firstSecond = positionDot(first, second);
            const auto secondSquared = positionDot(second, second);
            const auto firstThird = positionDot(first, third);
            const auto speed2 = speedValue * speedValue;
            const auto speed3 = speed2 * speedValue;
            const auto speed4 = speed2 * speed2;
            const auto speed5 = speed4 * speedValue;
            const auto speed7 = speed5 * speed2;
            const auto tangent = positionScale(first, 1.0 / speedValue);
            const auto curvature = positionAdd(positionScale(second, 1.0 / speed2),
                                                positionScale(first, -firstSecond / speed4));
            auto thirdDerivative = positionScale(third, 1.0 / speed3);
            thirdDerivative = positionAdd(thirdDerivative, positionScale(second, -3.0 * firstSecond / speed5));
            thirdDerivative = positionAdd(thirdDerivative,
                positionScale(first, -(secondSquared + firstThird) / speed5 + 4.0 * firstSecond * firstSecond / speed7));

            return {tangent, curvature, thirdDerivative};
        }
    };

    std::string samplingError(const path_tempo::SamplingError error) {
        return std::format("geometry sampling failed with code {}", static_cast<int>(error));
    }

    std::expected<std::vector<path_tempo::SampledPathPiece<DOF>>, std::string> sampleGeometry(const path_tempo::example::GeometryFile &geometry) {
        const auto pieceCapacity = std::accumulate(geometry.curves.begin(), geometry.curves.end(), std::size_t {0},
            [](const std::size_t total, const auto &curve) { return total + curve.feeds.size(); });
        std::vector<path_tempo::SampledPathPiece<DOF>> result;
        result.reserve(pieceCapacity);
        auto nextId = path_tempo::PathPieceId {1};

        for (const auto &curve : geometry.curves) {
            auto sampled = std::visit([&](const auto &value) -> std::expected<std::vector<path_tempo::SampledPathPiece<DOF>>, std::string> {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::same_as<T, path_tempo::example::LineGeometry>) {
                    if (curve.feeds.size() != 1) {
                        return std::unexpected("a line must have exactly one feed");
                    }

                    Position delta {};
                    for (std::size_t component = 0; component < DOF; ++component) {
                        delta[component] = value.to[component] - value.from[component];
                    }
                    const auto fullLength = positionMagnitude(delta);
                    if (curve.toDistance > fullLength + 1e-10 * std::max(1.0, fullLength)) {
                        return std::unexpected("line interval exceeds its geometry length");
                    }

                    path_tempo::Line<DOF> line;
                    for (std::size_t component = 0; component < DOF; ++component) {
                        line.from[component] = value.from[component] + delta[component] * curve.fromDistance / fullLength;
                        line.to[component] = value.from[component] + delta[component] * curve.toDistance / fullLength;
                    }

                    auto piece = path_tempo::sampleLine(line, nextId, curve.feeds.front());
                    if (!piece) {
                        return std::unexpected(samplingError(piece.error()));
                    }

                    return std::vector {std::move(*piece)};
                } else if constexpr (std::same_as<T, path_tempo::example::ArcGeometry>) {
                    if (curve.feeds.size() != 1) {
                        return std::unexpected("an arc must have exactly one feed");
                    }

                    auto arc = EndpointExactArc::create(value);
                    if (!arc) {
                        return std::unexpected(arc.error());
                    }

                    if (curve.toDistance > arc->length() + 1e-10 * std::max(1.0, arc->length())) {
                        return std::unexpected("arc interval exceeds its geometry length");
                    }

                    auto piece = path_tempo::sampleArcLengthCurve<DOF>(
                        curve.toDistance - curve.fromDistance, nextId, curve.feeds.front(),
                        CURVE_SAMPLE_INTERVALS, [&](const double distance) {
                            return arc->stateAtDistance(curve.fromDistance + distance);
                        });
                    if (!piece) {
                        return std::unexpected(samplingError(piece.error()));
                    }

                    return std::vector {std::move(*piece)};
                } else {
                    auto pieces = path_tempo::sampleBSpline(path_tempo::BSpline<DOF> {
                        .degree = value.degree,
                        .controls = value.controls,
                        .knots = value.knots,
                    }, nextId, curve.feeds.front(), CURVE_SAMPLE_INTERVALS);
                    if (!pieces) {
                        return std::unexpected(samplingError(pieces.error()));
                    }

                    if (pieces->size() != curve.feeds.size()) {
                        return std::unexpected("B-spline feed count does not match its non-empty knot intervals");
                    }

                    const auto length = std::accumulate(pieces->begin(), pieces->end(), 0.0,
                        [](const double total, const auto &piece) { return total + piece.length; });
                    const auto expectedLength = curve.toDistance - curve.fromDistance;
                    if (curve.fromDistance > 1e-10 || std::abs(length - expectedLength) > 1e-9 * std::max(1.0, expectedLength)) {
                        return std::unexpected("partial or length-mismatched B-spline geometry is not supported by this example");
                    }

                    for (std::size_t index = 0; index < pieces->size(); ++index) {
                        (*pieces)[index].maxVelocity = curve.feeds[index];
                    }

                    return std::move(*pieces);
                }
            }, curve.geometry);
            if (!sampled) {
                return std::unexpected(std::format("curve {}: {}", &curve - geometry.curves.data() + 1, sampled.error()));
            }

            nextId += sampled->size();
            result.insert(result.end(), std::make_move_iterator(sampled->begin()), std::make_move_iterator(sampled->end()));
        }

        return result;
    }

    path_tempo::Limits<DOF> exampleLimits() {
        return {
            .pathAcceleration = 25.1,
            .pathJerk = 101.0,
            .coordinateVelocity = {3.33333333333, 3.33333333333, 3.33333333333},
            .coordinateAcceleration = {5.1, 5.1, 5.1},
            .coordinateJerk = {101.0, 101.0, 101.0},
        };
    }
}

int main(const int argc, char **argv) {
    if (argc > 2) {
        std::println(stderr, "usage: path_tempo_continuous_path_example [geometry.txt]");

        return 2;
    }

    const auto path = argc == 2
        ? std::filesystem::path {argv[1]}
        : std::filesystem::path {PATH_TEMPO_EXAMPLE_DATA_DIR} / "continuous_path.txt";

    auto geometry = path_tempo::example::loadGeometryFile(path);
    if (!geometry) {
        std::println(stderr, "{}", geometry.error());

        return 1;
    }

    auto sampled = sampleGeometry(*geometry);
    if (!sampled) {
        std::println(stderr, "{}", sampled.error());

        return 1;
    }

    std::vector<path_tempo::PathPiece<DOF>> pieces;
    pieces.reserve(sampled->size());
    for (const auto &piece : *sampled) {
        pieces.push_back(piece.view());
    }

    const auto pathLength = std::accumulate(pieces.begin(), pieces.end(), 0.0,
        [](const double total, const auto &piece) { return total + piece.length; });
    const auto unit = geometry->unit == path_tempo::example::Unit::Inch ? "in" : "mm";

    std::println("loaded {} curves as {} PathPieces; path length: {:.6f} {}", geometry->curves.size(), pieces.size(), pathLength, unit);

    path_tempo::PathPlanningSettings settings;
    settings.applySampledCorrections = true;

    path_tempo::PathPlanner planner;
    const auto planningStarted = std::chrono::steady_clock::now();
    const auto planned = planner.solve(path_tempo::PathPlanningRequest<DOF> {
        .pieces = pieces,
        .beginning = {},
        .ending = {},
        .limits = exampleLimits(),
        .settings = settings,
    });
    const auto planningSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - planningStarted).count();
    if (!planned) {
        std::println(stderr, "planning failed after {:.6f} s: {}", planningSeconds, planned.error().message);

        return 1;
    }

    const auto trajectoryDuration = std::accumulate(planned->timeLaw.segments.begin(), planned->timeLaw.segments.end(), 0.0,
        [](const double total, const auto &segment) { return total + segment.duration; });

    std::println("trajectory duration: {:.6f} s", trajectoryDuration);
    std::println("cubic time segments: {}", planned->timeLaw.segments.size());
    std::println("sampled correction passes: {}", planned->diagnostics.correctionPasses);
    std::println("planning wall time: {:.6f} s", planningSeconds);

    return 0;
}
