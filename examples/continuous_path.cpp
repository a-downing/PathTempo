#include "continuous_path_options.h"
#include "endpoint_exact_arc.h"
#include "geometry_text_loader.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <expected>
#include <format>
#include <iterator>
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
    using Position = path_tempo::Vector<DOF>;

    double positionDot(const Position &left, const Position &right) {
        return std::inner_product(left.begin(), left.end(), right.begin(), 0.0);
    }

    double positionMagnitude(const Position &value) {
        return std::sqrt(positionDot(value, value));
    }

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

                    auto arc = path_tempo::example::EndpointExactArc::create(value);
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

    path_tempo::Limits<DOF> exampleLimits(const double maximumAcceleration, const double maximumJerk) {
        constexpr auto UNBOUNDED = std::numeric_limits<double>::infinity();

        return {
            .pathAcceleration = maximumAcceleration,
            .pathJerk = maximumJerk,
            .coordinateVelocity = {UNBOUNDED, UNBOUNDED, UNBOUNDED},
            .coordinateAcceleration = {maximumAcceleration, maximumAcceleration, maximumAcceleration},
            .coordinateJerk = {maximumJerk, maximumJerk, maximumJerk},
        };
    }
}

int main(const int argc, char **argv) {
    std::vector<std::string_view> arguments;
    arguments.reserve(argc > 1 ? static_cast<std::size_t>(argc - 1) : 0);
    for (auto index = 1; index < argc; ++index) {
        arguments.emplace_back(argv[index]);
    }

    const auto options = path_tempo::example::parseContinuousPathOptions(arguments);
    if (!options) {
        std::println(stderr, "{}", options.error());
        std::println(stderr, "{}", path_tempo::example::CONTINUOUS_PATH_USAGE);

        return 2;
    }

    auto geometry = path_tempo::example::loadGeometryFile(options->geometryPath);
    if (!geometry) {
        std::println(stderr, "{}", geometry.error());

        return 1;
    }

    if (geometry->unit != options->unit) {
        std::println(stderr, "geometry file uses {}, but -u specifies {}",
            path_tempo::example::unitArgument(geometry->unit), path_tempo::example::unitArgument(options->unit));

        return 2;
    }

    auto sampled = sampleGeometry(*geometry);
    if (!sampled) {
        std::println(stderr, "{}", sampled.error());

        return 1;
    }

    for (auto &piece : *sampled) {
        const auto scaled = path_tempo::example::scaleVelocityLimit(piece.maxVelocity, options->velocityScale);

        if (!scaled) {
            std::println(stderr, "path piece {}: {}", piece.id, scaled.error());

            return 1;
        }

        piece.maxVelocity = *scaled;
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
    std::println("mode: {}; limits: acceleration {:.6f} {}/s^2; jerk {:.6f} {}/s^3; PathPiece velocity scale {:.6f}",
        path_tempo::example::modeArgument(options->mode), options->maximumAcceleration, unit, options->maximumJerk, unit, options->velocityScale);

    path_tempo::PathPlanningSettings settings;
    settings.applySampledCorrections = true;
    settings.maximumCorrectionPasses = 128;
    settings.boundaryAccelerationMode = options->mode == path_tempo::example::ContinuousPathMode::Zero
        ? path_tempo::BoundaryAccelerationMode::Zero
        : path_tempo::BoundaryAccelerationMode::Optimized;

    path_tempo::PathPlanner planner;
    const auto planningStarted = std::chrono::steady_clock::now();
    const auto planned = planner.solve(path_tempo::PathPlanningRequest<DOF> {
        .pieces = pieces,
        .beginning = {},
        .ending = {},
        .limits = exampleLimits(options->maximumAcceleration, options->maximumJerk),
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
