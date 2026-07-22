#include "continuous_path_options.h"

#include <charconv>
#include <cmath>
#include <optional>
#include <utility>

namespace path_tempo::example {
    namespace {
        std::expected<double, std::string> positiveDouble(const std::string_view option, const std::string_view argument) {
            double value = 0.0;
            const auto parsed = std::from_chars(argument.data(), argument.data() + argument.size(), value);

            if (parsed.ec != std::errc {} || parsed.ptr != argument.data() + argument.size() || !std::isfinite(value) || value <= 0.0) {
                return std::unexpected("option " + std::string(option) + " requires a finite positive number, found '" + std::string(argument) + "'");
            }

            return value;
        }

        std::expected<std::string_view, std::string> optionValue(const std::span<const std::string_view> arguments, std::size_t &index) {
            if (index + 1 >= arguments.size()) {
                return std::unexpected("option " + std::string(arguments[index]) + " requires a value");
            }

            return arguments[++index];
        }
    }

    std::expected<ContinuousPathOptions, std::string> parseContinuousPathOptions(const std::span<const std::string_view> arguments) {
        std::optional<ContinuousPathMode> mode;
        std::optional<Unit> unit;
        std::optional<double> maximumAcceleration;
        std::optional<double> maximumJerk;
        std::optional<double> velocityScale;
        std::optional<std::filesystem::path> geometryPath;
        auto optionsEnded = false;

        for (std::size_t index = 0; index < arguments.size(); ++index) {
            const auto argument = arguments[index];

            if (!optionsEnded && argument == "--") {
                optionsEnded = true;
                continue;
            }

            if (!optionsEnded && argument == "-u") {
                if (unit) {
                    return std::unexpected("option -u may be supplied only once");
                }

                auto value = optionValue(arguments, index);
                if (!value) {
                    return std::unexpected(value.error());
                }

                if (*value == "inch") {
                    unit = Unit::Inch;
                } else if (*value == "mm") {
                    unit = Unit::Millimeter;
                } else {
                    return std::unexpected("option -u must be 'inch' or 'mm', found '" + std::string(*value) + "'");
                }

                continue;
            }

            if (!optionsEnded && argument == "-m") {
                if (mode) {
                    return std::unexpected("option -m may be supplied only once");
                }

                auto value = optionValue(arguments, index);
                if (!value) {
                    return std::unexpected(value.error());
                }

                if (*value == "zero") {
                    mode = ContinuousPathMode::Zero;
                } else if (*value == "optimized") {
                    mode = ContinuousPathMode::Optimized;
                } else {
                    return std::unexpected("option -m must be 'zero' or 'optimized', found '" + std::string(*value) + "'");
                }

                continue;
            }

            if (!optionsEnded && (argument == "-a" || argument == "-j" || argument == "-v")) {
                auto value = optionValue(arguments, index);
                if (!value) {
                    return std::unexpected(value.error());
                }

                auto parsed = positiveDouble(argument, *value);
                if (!parsed) {
                    return std::unexpected(parsed.error());
                }

                if (argument == "-a") {
                    if (maximumAcceleration) {
                        return std::unexpected("option -a may be supplied only once");
                    }
                    maximumAcceleration = *parsed;
                } else if (argument == "-j") {
                    if (maximumJerk) {
                        return std::unexpected("option -j may be supplied only once");
                    }
                    maximumJerk = *parsed;
                } else {
                    if (velocityScale) {
                        return std::unexpected("option -v may be supplied only once");
                    }
                    velocityScale = *parsed;
                }

                continue;
            }

            if (!optionsEnded && argument.starts_with('-')) {
                return std::unexpected("unknown option '" + std::string(argument) + "'");
            }

            if (geometryPath) {
                return std::unexpected("exactly one geometry file is required");
            }

            if (argument.empty()) {
                return std::unexpected("geometry file path must not be empty");
            }

            geometryPath = std::filesystem::path {argument};
        }

        if (!mode) {
            return std::unexpected("required option -m is missing");
        }
        if (!unit) {
            return std::unexpected("required option -u is missing");
        }
        if (!maximumAcceleration) {
            return std::unexpected("required option -a is missing");
        }
        if (!maximumJerk) {
            return std::unexpected("required option -j is missing");
        }
        if (!geometryPath) {
            return std::unexpected("required geometry file is missing");
        }

        return ContinuousPathOptions {
            .mode = *mode,
            .unit = *unit,
            .maximumAcceleration = *maximumAcceleration,
            .maximumJerk = *maximumJerk,
            .velocityScale = velocityScale.value_or(1.0),
            .geometryPath = std::move(*geometryPath),
        };
    }

    std::expected<double, std::string> scaleVelocityLimit(const double velocity, const double scale) {
        const auto scaled = velocity * scale;

        if (!std::isfinite(velocity) || velocity <= 0.0 || !std::isfinite(scale) || scale <= 0.0 || !std::isfinite(scaled) || scaled <= 0.0) {
            return std::unexpected("PathPiece velocity limit and scale must produce a finite positive value");
        }

        return scaled;
    }

    std::string_view modeArgument(const ContinuousPathMode mode) {
        return mode == ContinuousPathMode::Zero ? "zero" : "optimized";
    }

    std::string_view unitArgument(const Unit unit) {
        return unit == Unit::Inch ? "inch" : "mm";
    }
}
