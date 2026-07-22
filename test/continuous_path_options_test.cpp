#include "../examples/continuous_path_options.h"

#include <cstdlib>
#include <expected>
#include <initializer_list>
#include <limits>
#include <print>
#include <string_view>
#include <vector>

namespace {
    void require(const bool condition, const std::string_view message) {
        if (!condition) {
            std::println(stderr, "test failure: {}", message);
            std::exit(1);
        }
    }

    std::expected<path_tempo::example::ContinuousPathOptions, std::string> parse(const std::initializer_list<std::string_view> arguments) {
        return path_tempo::example::parseContinuousPathOptions(std::vector<std::string_view> {arguments});
    }

    void testRequiredOptionsAndDefaultVelocityScale() {
        const auto parsed = parse({"reference.txt", "-j", "12.5", "-m", "zero", "-u", "inch", "-a", "4.25"});

        require(parsed.has_value(), "required options should parse in any order");
        require(parsed->mode == path_tempo::example::ContinuousPathMode::Zero, "zero should select zero boundary acceleration");
        require(parsed->unit == path_tempo::example::Unit::Inch, "inch should select inch units");
        require(parsed->maximumAcceleration == 4.25, "-a should set the acceleration limit");
        require(parsed->maximumJerk == 12.5, "-j should set the jerk limit");
        require(parsed->velocityScale == 1.0, "the PathPiece velocity scale should default to one");
        require(parsed->geometryPath == "reference.txt", "the positional argument should set the geometry path");
    }

    void testExplicitVelocityScaleAndMillimeters() {
        const auto parsed = parse({"-u", "mm", "-m", "optimized", "-a", "250", "-v", "1.75", "geometry.txt", "-j", "1000"});

        require(parsed.has_value(), "the optional velocity scale should parse");
        require(parsed->mode == path_tempo::example::ContinuousPathMode::Optimized, "optimized should select optimized boundary acceleration");
        require(parsed->unit == path_tempo::example::Unit::Millimeter, "mm should select millimeter units");
        require(parsed->velocityScale == 1.75, "-v should set the PathPiece velocity scale");
    }

    void testRequiredArguments() {
        require(!parse({"-u", "inch", "-a", "1", "-j", "2", "geometry.txt"}), "-m should be required");
        require(!parse({"-m", "zero", "-a", "1", "-j", "2", "geometry.txt"}), "-u should be required");
        require(!parse({"-m", "zero", "-u", "inch", "-j", "2", "geometry.txt"}), "-a should be required");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "geometry.txt"}), "-j should be required");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "-j", "2"}), "the geometry path should be required");
    }

    void testInvalidArguments() {
        require(!parse({"-m", "fast", "-u", "inch", "-a", "1", "-j", "2", "geometry.txt"}), "unknown modes should fail");
        require(!parse({"-m", "zero", "-u", "meter", "-a", "1", "-j", "2", "geometry.txt"}), "unknown units should fail");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "0", "-j", "2", "geometry.txt"}), "zero acceleration should fail");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "-j", "nan", "geometry.txt"}), "non-finite jerk should fail");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "-j", "2", "-v", "-1", "geometry.txt"}), "negative velocity scale should fail");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "-a", "2", "-j", "3", "geometry.txt"}), "duplicate options should fail");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "-j", "2", "first.txt", "second.txt"}), "multiple geometry files should fail");
        require(!parse({"-m", "zero", "-u", "inch", "-a", "1", "-j", "2", "-x", "geometry.txt"}), "unknown options should fail");
    }

    void testVelocityLimitScaling() {
        const auto scaled = path_tempo::example::scaleVelocityLimit(12.0, 0.5);

        require(scaled.has_value() && *scaled == 6.0, "finite velocity scaling should succeed");
        require(!path_tempo::example::scaleVelocityLimit(std::numeric_limits<double>::max(), 2.0), "overflowing velocity scaling should fail");
    }
}

int main() {
    testRequiredOptionsAndDefaultVelocityScale();
    testExplicitVelocityScaleAndMillimeters();
    testRequiredArguments();
    testInvalidArguments();
    testVelocityLimitScaling();

    return 0;
}
