#pragma once

#include <expected>
#include <string>
#include <string_view>

namespace path_tempo::example {
    std::expected<double, std::string> canonicalCurveIntervalEnd(double fromDistance, double toDistance, double geometryLength, std::string_view geometryType);
}
