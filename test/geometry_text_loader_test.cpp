#include "../examples/geometry_text_loader.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <print>
#include <string>
#include <string_view>
#include <system_error>
#include <variant>

namespace {
    void require(const bool condition, const std::string_view message) {
        if (!condition) {
            std::println(stderr, "test failure: {}", message);
            std::exit(1);
        }
    }

    class TemporaryGeometryFile {
        std::filesystem::path m_path;

    public:
        explicit TemporaryGeometryFile(const std::string_view contents) {
            const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
            m_path = std::filesystem::temp_directory_path() / ("path_tempo_geometry_loader_" + std::to_string(suffix) + ".txt");
            std::ofstream output(m_path);
            output << contents;

            require(output.good(), "the temporary geometry file should be writable");
        }

        ~TemporaryGeometryFile() {
            std::error_code error;
            std::filesystem::remove(m_path, error);
        }

        const std::filesystem::path &path() const {
            return m_path;
        }
    };

    void testValidLineGeometry() {
        const TemporaryGeometryFile file {R"(
ngc_g64_geometry 1
units inch
curve_count 1
curve_interval 0 5
feed_count 1
feed 2.5
curve line
from 0 0 0
to 3 4 0
end_curve
end_geometry
)"};
        const auto geometry = path_tempo::example::loadGeometryFile(file.path());

        require(geometry.has_value(), "a valid line geometry file should load");
        require(geometry->unit == path_tempo::example::Unit::Inch, "the geometry unit should be preserved");
        require(geometry->curves.size() == 1, "the declared line curve should be loaded");
        require(geometry->curves.front().feeds.size() == 1 && geometry->curves.front().feeds.front() == 2.5, "the line feed should be preserved");
        require(std::holds_alternative<path_tempo::example::LineGeometry>(geometry->curves.front().geometry), "the line geometry type should be preserved");
    }

    void testTrailingTokenRejected() {
        const TemporaryGeometryFile file {R"(
ngc_g64_geometry 1
units millimeter
curve_count 0
end_geometry
trailing
)"};

        require(!path_tempo::example::loadGeometryFile(file.path()), "tokens after end_geometry should be rejected");
    }

    void testNonmonotoneSplineKnotsRejected() {
        const TemporaryGeometryFile file {R"(
ngc_g64_geometry 1
units inch
curve_count 1
curve_interval 0 1
feed_count 1
feed 1
curve bspline
degree 2
control_count 3
control 0 0 0
control 0.5 0 0
control 1 0 0
knot_count 6
knot 0
knot 0
knot 1
knot 0.5
knot 1
knot 1
end_curve
end_geometry
)"};

        require(!path_tempo::example::loadGeometryFile(file.path()), "nonmonotone B-spline knots should be rejected");
    }
}

int main() {
    testValidLineGeometry();
    testTrailingTokenRejected();
    testNonmonotoneSplineKnotsRejected();

    return 0;
}
