#pragma once

#include <expected>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>

#include "geometry_text_loader.h"

namespace path_tempo::example {
    inline constexpr std::string_view CONTINUOUS_PATH_USAGE = "usage: path_tempo_continuous_path_example -m <zero|optimized> -u <inch|mm> -a <acceleration> -j <jerk> [-v <velocity-scale>] <geometry.txt>";

    enum class ContinuousPathMode {
        Zero,
        Optimized,
    };

    struct ContinuousPathOptions {
        ContinuousPathMode mode = ContinuousPathMode::Zero;
        Unit unit = Unit::Inch;
        double maximumAcceleration = 0.0;
        double maximumJerk = 0.0;
        double velocityScale = 1.0;
        std::filesystem::path geometryPath;
    };

    std::expected<ContinuousPathOptions, std::string> parseContinuousPathOptions(std::span<const std::string_view> arguments);
    std::expected<double, std::string> scaleVelocityLimit(double velocity, double scale);
    std::string_view modeArgument(ContinuousPathMode mode);
    std::string_view unitArgument(Unit unit);
}
