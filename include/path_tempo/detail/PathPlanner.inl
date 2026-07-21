#pragma once

#include <algorithm>
#include <cmath>
#include <format>

namespace path_tempo {
    namespace detail {
        // Arc-length derivatives below this magnitude do not produce a stable
        // finite constraint quotient in double precision and are treated as zero.
        inline constexpr double GEOMETRY_COMPONENT_ZERO_TOLERANCE = 1e-15;
        // The squared tangent norm accumulates one rounding contribution per
        // coordinate; 1e-9 accepts sampled unit vectors without masking bad input.
        inline constexpr double UNIT_TANGENT_SQUARED_TOLERANCE = 1e-9;
        // Station endpoints originate in independent arc-length calculations,
        // so accept a small absolute roundoff floor near the distance origin.
        inline constexpr double STATION_DISTANCE_ABSOLUTE_TOLERANCE = 1e-12;
        // Scale endpoint coverage tolerance with piece length to accommodate
        // accumulated integration and inversion error on long pieces.
        inline constexpr double STATION_DISTANCE_RELATIVE_TOLERANCE = 1e-10;
        // Adjacent samples may differ by ordinary normalization noise, but a
        // larger tangent change represents a real C1 discontinuity.
        inline constexpr double TANGENT_CONTINUITY_TOLERANCE = 1e-9;
        // Curvature is a numerically derived second derivative, so its C2 check
        // is intentionally one order looser than the tangent check.
        inline constexpr double CURVATURE_CONTINUITY_TOLERANCE = 1e-8;
    }

    template<std::size_t DoF>
    std::expected<PlannedPath, PlanningError> PathPlanner::solve(const PathPlanningRequest<DoF> &request, const MaterializationCorrection &materializationCorrection) {
        if (request.pieces.empty()) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::InvalidInput,
                .message = "path planning requires at least one piece",
            });
        }

        if (!std::isfinite(request.beginning.velocity) || !std::isfinite(request.beginning.acceleration) || !std::isfinite(request.ending.velocity) || !std::isfinite(request.ending.acceleration) || request.beginning.velocity < 0.0 || request.ending.velocity < 0.0 || !std::isfinite(request.settings.linearSolveTimeLimit) || request.settings.linearSolveTimeLimit <= 0.0 || request.settings.simplexIterationLimit == 0 || request.settings.maximumCorrectionPasses == 0 || (request.settings.sequentialIterations > 0 && request.settings.lineSearchSteps == 0) || !std::isfinite(request.settings.velocityTrustFraction) || request.settings.velocityTrustFraction <= 0.0 || request.settings.velocityTrustFraction > 1.0 || !std::isfinite(request.settings.accelerationTrustFraction) || request.settings.accelerationTrustFraction <= 0.0 || request.settings.accelerationTrustFraction > 1.0) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::InvalidInput,
                .message = "path planning received invalid boundary state or solver settings",
            });
        }

        if (std::isnan(request.limits.pathAcceleration) || request.limits.pathAcceleration <= 0.0 || std::isnan(request.limits.pathJerk) || request.limits.pathJerk <= 0.0) {
            return std::unexpected(PlanningError {
                .code = PlanningErrorCode::InvalidInput,
                .message = "path planning requires positive aggregate acceleration and jerk limits",
            });
        }

        std::vector<LocalPiece> localPieces;
        localPieces.reserve(request.pieces.size());

        for (std::size_t pieceIndex = 0; pieceIndex < request.pieces.size(); ++pieceIndex) {
            const auto &piece = request.pieces[pieceIndex];

            if (std::ranges::any_of(localPieces, [&](const LocalPiece &existing) { return existing.id == piece.id; })) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = std::format("path piece {} repeats piece ID {}", pieceIndex, piece.id),
                });
            }

            if (!std::isfinite(piece.length) || piece.length <= 0.0 || !std::isfinite(piece.maxVelocity) || piece.maxVelocity <= 0.0 || piece.stations.empty()) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = std::format("path piece {} has invalid length, maximum velocity, or stations", pieceIndex),
                });
            }

            auto maximumVelocity = std::min(piece.maxVelocity, piece.initialLimits.velocity);
            auto maximumAcceleration = std::min(request.limits.pathAcceleration, piece.initialLimits.acceleration);
            auto maximumJerk = std::min(request.limits.pathJerk, piece.initialLimits.jerk);
            auto previousDistance = -1.0;
            std::vector<LocalPiece::Station> stations;
            stations.reserve(piece.stations.size());

            for (std::size_t stationIndex = 0; stationIndex < piece.stations.size(); ++stationIndex) {
                const auto &station = piece.stations[stationIndex];

                if (!std::isfinite(station.distance) || station.distance < previousDistance || station.distance < 0.0 || station.distance > piece.length) {
                    return std::unexpected(PlanningError {
                        .code = PlanningErrorCode::InvalidInput,
                        .message = std::format("path piece {} has invalid station ordering at station {}: distance={} previous={} length={}", pieceIndex, stationIndex, station.distance, previousDistance, piece.length),
                    });
                }

                auto tangentSquared = 0.0;

                for (std::size_t axis = 0; axis < DoF; ++axis) {
                    const auto tangent = station.tangent[axis];
                    const auto curvature = station.curvature[axis];
                    const auto thirdDerivative = station.thirdDerivative[axis];
                    const auto coordinateVelocity = request.limits.coordinateVelocity[axis];
                    const auto coordinateAcceleration = request.limits.coordinateAcceleration[axis];
                    const auto coordinateJerk = request.limits.coordinateJerk[axis];

                    if (!std::isfinite(tangent) || !std::isfinite(curvature) || !std::isfinite(thirdDerivative) || std::isnan(coordinateVelocity) || coordinateVelocity <= 0.0 || std::isnan(coordinateAcceleration) || coordinateAcceleration <= 0.0 || std::isnan(coordinateJerk) || coordinateJerk <= 0.0) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("path piece {} has non-finite geometry or invalid coordinate limits", pieceIndex),
                        });
                    }

                    tangentSquared += tangent * tangent;

                    if (std::abs(tangent) > detail::GEOMETRY_COMPONENT_ZERO_TOLERANCE) {
                        maximumVelocity = std::min(maximumVelocity, coordinateVelocity / std::abs(tangent));
                        maximumAcceleration = std::min(maximumAcceleration, coordinateAcceleration / std::abs(tangent));
                        maximumJerk = std::min(maximumJerk, coordinateJerk / std::abs(tangent));
                    }

                    if (std::abs(curvature) > detail::GEOMETRY_COMPONENT_ZERO_TOLERANCE) {
                        maximumVelocity = std::min(maximumVelocity, std::sqrt(coordinateAcceleration / std::abs(curvature)));
                    }

                    if (std::abs(thirdDerivative) > detail::GEOMETRY_COMPONENT_ZERO_TOLERANCE) {
                        maximumVelocity = std::min(maximumVelocity, std::cbrt(coordinateJerk / std::abs(thirdDerivative)));
                    }
                }

                if (std::abs(tangentSquared - 1.0) > detail::UNIT_TANGENT_SQUARED_TOLERANCE) {
                    return std::unexpected(PlanningError {
                        .code = PlanningErrorCode::InvalidInput,
                        .message = std::format("path piece {} has a non-unit station tangent", pieceIndex),
                    });
                }

                auto curvatureSquared = 0.0;
                auto thirdDerivativeSquared = 0.0;

                for (std::size_t axis = 0; axis < DoF; ++axis) {
                    curvatureSquared += station.curvature[axis] * station.curvature[axis];
                    thirdDerivativeSquared += station.thirdDerivative[axis] * station.thirdDerivative[axis];
                }

                const auto curvatureMagnitude = std::sqrt(curvatureSquared);
                const auto thirdDerivativeMagnitude = std::sqrt(thirdDerivativeSquared);

                if (curvatureMagnitude > detail::GEOMETRY_COMPONENT_ZERO_TOLERANCE) {
                    maximumVelocity = std::min(maximumVelocity, std::sqrt(request.limits.pathAcceleration / curvatureMagnitude));
                }

                if (thirdDerivativeMagnitude > detail::GEOMETRY_COMPONENT_ZERO_TOLERANCE) {
                    maximumVelocity = std::min(maximumVelocity, std::cbrt(request.limits.pathJerk / thirdDerivativeMagnitude));
                }

                stations.push_back({
                    .distance = station.distance,
                    .tangent = station.tangent,
                    .curvature = station.curvature,
                    .thirdDerivative = station.thirdDerivative,
                });

                previousDistance = station.distance;
            }

            if (std::abs(piece.stations.front().distance) > detail::STATION_DISTANCE_ABSOLUTE_TOLERANCE || std::abs(piece.stations.back().distance - piece.length) > std::max(detail::STATION_DISTANCE_ABSOLUTE_TOLERANCE, piece.length * detail::STATION_DISTANCE_RELATIVE_TOLERANCE)) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = std::format("path piece {} stations do not cover the complete piece", pieceIndex),
                });
            }

            if (!std::isfinite(maximumVelocity) || maximumVelocity <= 0.0 || !std::isfinite(maximumAcceleration) || maximumAcceleration <= 0.0 || !std::isfinite(maximumJerk) || maximumJerk <= 0.0) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = std::format("path piece {} has no finite positive scalar limits", pieceIndex),
                });
            }

            localPieces.push_back({piece.id, piece.length, maximumVelocity, maximumAcceleration, maximumJerk, std::move(stations)});

            if (pieceIndex > 0) {
                const auto &previous = request.pieces[pieceIndex - 1].stations.back();
                const auto &current = piece.stations.front();

                for (std::size_t axis = 0; axis < DoF; ++axis) {
                    if (std::abs(previous.tangent[axis] - current.tangent[axis]) > detail::TANGENT_CONTINUITY_TOLERANCE) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("path pieces {} and {} are not tangent-continuous", pieceIndex - 1, pieceIndex),
                        });
                    }

                    if (std::abs(previous.curvature[axis] - current.curvature[axis]) > detail::CURVATURE_CONTINUITY_TOLERANCE) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("path pieces {} and {} are not curvature-continuous", pieceIndex - 1, pieceIndex),
                        });
                    }
                }
            }
        }

        CoupledLimits coupledLimits {
            .pathAcceleration = request.limits.pathAcceleration,
            .pathJerk = request.limits.pathJerk,
            .coordinateVelocity = {request.limits.coordinateVelocity.begin(), request.limits.coordinateVelocity.end()},
            .coordinateAcceleration = {request.limits.coordinateAcceleration.begin(), request.limits.coordinateAcceleration.end()},
            .coordinateJerk = {request.limits.coordinateJerk.begin(), request.limits.coordinateJerk.end()},
        };

        return solveLocal(localPieces, request.beginning, request.ending, coupledLimits, request.settings, materializationCorrection);
    }
}
