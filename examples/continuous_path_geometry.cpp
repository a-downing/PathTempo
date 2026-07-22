#include "continuous_path_geometry.h"

#include <algorithm>
#include <cmath>

namespace path_tempo::example {
    std::expected<double, std::string> canonicalCurveIntervalEnd(const double fromDistance, const double toDistance, const double geometryLength, const std::string_view geometryType) {
        if (!std::isfinite(geometryLength) || geometryLength <= 0.0) {
            return std::unexpected(std::string(geometryType) + " geometry must have finite positive length");
        }

        if (fromDistance >= geometryLength) {
            return std::unexpected(std::string(geometryType) + " interval begins at or beyond its geometry length");
        }

        const auto tolerance = 1e-10 * std::max(1.0, geometryLength);

        if (toDistance > geometryLength + tolerance) {
            return std::unexpected(std::string(geometryType) + " interval exceeds its geometry length");
        }

        return std::min(toDistance, geometryLength);
    }
}
