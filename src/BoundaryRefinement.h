#pragma once

#include <cmath>

namespace path_tempo::detail {
    inline double resolveBoundaryRepairWeight(const double weight, const double proposedAcceleration,
        const double maximumJerk, const double durationResolution) {
        if (weight >= 1.0 || !std::isfinite(maximumJerk)) {
            return weight;
        }

        const auto accelerationRampDuration = std::abs(proposedAcceleration * weight) / maximumJerk;

        return accelerationRampDuration <= durationResolution ? 0.0 : weight;
    }

    inline bool preferBoundaryAccelerationCandidate(const double candidateDuration, const double candidateAcceleration,
        const double bestDuration, const double bestAcceleration, const double durationResolution) {
        if (candidateDuration + durationResolution < bestDuration) {
            return true;
        }

        if (bestDuration + durationResolution < candidateDuration) {
            return false;
        }

        return std::abs(candidateAcceleration) < std::abs(bestAcceleration);
    }
}
