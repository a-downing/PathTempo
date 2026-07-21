#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <filesystem>
#include <string>
#include <variant>
#include <vector>

namespace path_tempo::example {
    using Position = std::array<double, 3>;
    using Vector3 = std::array<double, 3>;

    enum class Unit {
        Inch,
        Millimeter,
    };

    struct LineGeometry {
        Position from {};
        Position to {};
    };

    struct ArcGeometry {
        Position from {};
        Position to {};
        Vector3 center {};
        Vector3 axis {};
    };

    struct BSplineGeometry {
        std::size_t degree = 0;
        std::vector<Position> controls;
        std::vector<double> knots;
    };

    using CurveGeometry = std::variant<LineGeometry, ArcGeometry, BSplineGeometry>;

    struct Curve {
        double fromDistance = 0.0;
        double toDistance = 0.0;
        std::vector<double> feeds;
        CurveGeometry geometry;
    };

    struct GeometryFile {
        Unit unit = Unit::Inch;
        std::vector<Curve> curves;
    };

    std::expected<GeometryFile, std::string> loadGeometryFile(const std::filesystem::path &path);
}
