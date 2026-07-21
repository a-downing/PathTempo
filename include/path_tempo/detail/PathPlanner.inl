#pragma once

#include <algorithm>
#include <cmath>
#include <format>

namespace path_tempo {
    template<std::size_t DoF>
    std::expected<PlannedPath, PlanningError> PathPlanner::solve(const PathPlanningRequest<DoF> &request) {
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

            if (!std::isfinite(piece.length) || piece.length <= 0.0 || !std::isfinite(piece.programmedVelocity) || piece.programmedVelocity <= 0.0 || piece.stations.empty()) {
                return std::unexpected(PlanningError {
                    .code = PlanningErrorCode::InvalidInput,
                    .message = std::format("path piece {} has invalid length, programmed velocity, or stations", pieceIndex),
                });
            }

            auto maximumVelocity = std::min(piece.programmedVelocity, piece.initialLimits.velocity);
            auto maximumAcceleration = std::min(request.limits.pathAcceleration, piece.initialLimits.acceleration);
            auto maximumJerk = std::min(request.limits.pathJerk, piece.initialLimits.jerk);
            auto previousDistance = -1.0;
            std::vector<LocalPiece::Station> stations;
            stations.reserve(piece.stations.size());

            for (const auto &station : piece.stations) {
                if (!std::isfinite(station.distance) || station.distance < previousDistance || station.distance < 0.0 || station.distance > piece.length) {
                    return std::unexpected(PlanningError {
                        .code = PlanningErrorCode::InvalidInput,
                        .message = std::format("path piece {} has invalid station ordering", pieceIndex),
                    });
                }

                auto tangentSquared = 0.0;

                for (std::size_t axis = 0; axis < DoF; ++axis) {
                    const auto tangent = station.tangent[axis];
                    const auto curvature = station.curvature[axis];
                    const auto thirdDerivative = station.thirdDerivative[axis];
                    const auto axisVelocity = request.limits.axisVelocity[axis];
                    const auto axisAcceleration = request.limits.axisAcceleration[axis];
                    const auto axisJerk = request.limits.axisJerk[axis];

                    if (!std::isfinite(tangent) || !std::isfinite(curvature) || !std::isfinite(thirdDerivative) || std::isnan(axisVelocity) || axisVelocity <= 0.0 || std::isnan(axisAcceleration) || axisAcceleration <= 0.0 || std::isnan(axisJerk) || axisJerk <= 0.0) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("path piece {} has non-finite geometry or invalid axis limits", pieceIndex),
                        });
                    }

                    tangentSquared += tangent * tangent;

                    if (std::abs(tangent) > 1e-15) {
                        maximumVelocity = std::min(maximumVelocity, axisVelocity / std::abs(tangent));
                        maximumAcceleration = std::min(maximumAcceleration, axisAcceleration / std::abs(tangent));
                        maximumJerk = std::min(maximumJerk, axisJerk / std::abs(tangent));
                    }

                    if (std::abs(curvature) > 1e-15) {
                        maximumVelocity = std::min(maximumVelocity, std::sqrt(axisAcceleration / std::abs(curvature)));
                    }

                    if (std::abs(thirdDerivative) > 1e-15) {
                        maximumVelocity = std::min(maximumVelocity, std::cbrt(axisJerk / std::abs(thirdDerivative)));
                    }
                }

                if (std::abs(tangentSquared - 1.0) > 1e-9) {
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

                if (curvatureMagnitude > 1e-15) {
                    maximumVelocity = std::min(maximumVelocity, std::sqrt(request.limits.pathAcceleration / curvatureMagnitude));
                }

                if (thirdDerivativeMagnitude > 1e-15) {
                    maximumVelocity = std::min(maximumVelocity, std::cbrt(request.limits.pathJerk / thirdDerivativeMagnitude));
                }

                stations.push_back({
                    .distance = station.distance,
                    .tangent = {station.tangent.begin(), station.tangent.end()},
                    .curvature = {station.curvature.begin(), station.curvature.end()},
                    .thirdDerivative = {station.thirdDerivative.begin(), station.thirdDerivative.end()},
                });

                previousDistance = station.distance;
            }

            if (std::abs(piece.stations.front().distance) > 1e-12 || std::abs(piece.stations.back().distance - piece.length) > std::max(1e-12, piece.length * 1e-10)) {
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
                    if (std::abs(previous.tangent[axis] - current.tangent[axis]) > 1e-9) {
                        return std::unexpected(PlanningError {
                            .code = PlanningErrorCode::InvalidInput,
                            .message = std::format("path pieces {} and {} are not tangent-continuous", pieceIndex - 1, pieceIndex),
                        });
                    }

                    if (std::abs(previous.curvature[axis] - current.curvature[axis]) > 1e-8) {
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
            .axisVelocity = {request.limits.axisVelocity.begin(), request.limits.axisVelocity.end()},
            .axisAcceleration = {request.limits.axisAcceleration.begin(), request.limits.axisAcceleration.end()},
            .axisJerk = {request.limits.axisJerk.begin(), request.limits.axisJerk.end()},
        };

        return solveLocal(localPieces, request.beginning, request.ending, coupledLimits, request.settings);
    }
}
