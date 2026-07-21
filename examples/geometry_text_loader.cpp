#include "geometry_text_loader.h"

#include <charconv>
#include <cmath>
#include <fstream>
#include <string_view>
#include <utility>

namespace path_tempo::example {
    namespace {
        class TokenReader {
            std::ifstream m_input;
            std::filesystem::path m_path;
            std::size_t m_token = 0;

        public:
            explicit TokenReader(std::filesystem::path path) : m_input(path), m_path(std::move(path)) { }

            bool open() const { return m_input.is_open(); }

            std::expected<std::string, std::string> read() {
                std::string value;
                if (!(m_input >> value)) {
                    return std::unexpected(error("unexpected end of file"));
                }

                ++m_token;

                return value;
            }

            std::expected<void, std::string> expect(const std::string_view expected) {
                auto value = read();
                if (!value) {
                    return std::unexpected(value.error());
                }

                if (*value != expected) {
                    return std::unexpected(error("expected '" + std::string(expected) + "', found '" + *value + "'"));
                }

                return {};
            }

            std::expected<double, std::string> readDouble() {
                auto token = read();
                if (!token) {
                    return std::unexpected(token.error());
                }

                double value = 0.0;
                const auto result = std::from_chars(token->data(), token->data() + token->size(), value);
                if (result.ec != std::errc {} || result.ptr != token->data() + token->size() || !std::isfinite(value)) {
                    return std::unexpected(error("expected a finite double, found '" + *token + "'"));
                }

                return value;
            }

            std::expected<std::size_t, std::string> readCount() {
                auto token = read();
                if (!token) {
                    return std::unexpected(token.error());
                }

                std::size_t value = 0;
                const auto result = std::from_chars(token->data(), token->data() + token->size(), value);
                if (result.ec != std::errc {} || result.ptr != token->data() + token->size()) {
                    return std::unexpected(error("expected a non-negative integer, found '" + *token + "'"));
                }

                return value;
            }

            std::expected<void, std::string> expectEnd() {
                std::string trailing;
                if (m_input >> trailing) {
                    return std::unexpected(error("unexpected trailing token '" + trailing + "'"));
                }

                if (!m_input.eof()) {
                    return std::unexpected(error("failed while reading the file"));
                }

                return {};
            }

            std::string error(const std::string_view message) const {
                return m_path.string() + ": token " + std::to_string(m_token + 1) + ": " + std::string(message);
            }
        };

        template<std::size_t Size>
        std::expected<std::array<double, Size>, std::string> readVector(TokenReader &reader, const std::string_view label) {
            if (auto expected = reader.expect(label); !expected) {
                return std::unexpected(expected.error());
            }

            std::array<double, Size> result {};
            for (auto &component : result) {
                auto value = reader.readDouble();
                if (!value) {
                    return std::unexpected(value.error());
                }

                component = *value;
            }

            return result;
        }

        std::expected<CurveGeometry, std::string> readLine(TokenReader &reader) {
            auto from = readVector<3>(reader, "from");
            if (!from) {
                return std::unexpected(from.error());
            }

            auto to = readVector<3>(reader, "to");
            if (!to) {
                return std::unexpected(to.error());
            }

            return LineGeometry {*from, *to};
        }

        std::expected<CurveGeometry, std::string> readArc(TokenReader &reader) {
            auto from = readVector<3>(reader, "from");
            if (!from) {
                return std::unexpected(from.error());
            }

            auto to = readVector<3>(reader, "to");
            if (!to) {
                return std::unexpected(to.error());
            }

            auto center = readVector<3>(reader, "center");
            if (!center) {
                return std::unexpected(center.error());
            }

            auto axis = readVector<3>(reader, "axis");
            if (!axis) {
                return std::unexpected(axis.error());
            }

            return ArcGeometry {*from, *to, *center, *axis};
        }

        std::expected<CurveGeometry, std::string> readBSpline(TokenReader &reader) {
            if (auto expected = reader.expect("degree"); !expected) {
                return std::unexpected(expected.error());
            }

            auto degree = reader.readCount();
            if (!degree) {
                return std::unexpected(degree.error());
            }

            if (auto expected = reader.expect("control_count"); !expected) {
                return std::unexpected(expected.error());
            }

            auto controlCount = reader.readCount();
            if (!controlCount) {
                return std::unexpected(controlCount.error());
            }

            if (*degree == 0 || *controlCount <= *degree) {
                return std::unexpected(reader.error("B-spline degree and control count are inconsistent"));
            }

            std::vector<Position> controls;
            controls.reserve(*controlCount);
            for (std::size_t index = 0; index < *controlCount; ++index) {
                auto control = readVector<3>(reader, "control");
                if (!control) {
                    return std::unexpected(control.error());
                }

                controls.push_back(*control);
            }

            if (auto expected = reader.expect("knot_count"); !expected) {
                return std::unexpected(expected.error());
            }

            auto knotCount = reader.readCount();
            if (!knotCount) {
                return std::unexpected(knotCount.error());
            }

            if (*knotCount != controls.size() + *degree + 1) {
                return std::unexpected(reader.error("B-spline knot count must equal control count plus degree plus one"));
            }

            std::vector<double> knots;
            knots.reserve(*knotCount);
            for (std::size_t index = 0; index < *knotCount; ++index) {
                if (auto expected = reader.expect("knot"); !expected) {
                    return std::unexpected(expected.error());
                }

                auto knot = reader.readDouble();
                if (!knot) {
                    return std::unexpected(knot.error());
                }

                if (!knots.empty() && *knot < knots.back()) {
                    return std::unexpected(reader.error("B-spline knots must be nondecreasing"));
                }

                knots.push_back(*knot);
            }

            return BSplineGeometry {*degree, std::move(controls), std::move(knots)};
        }

        std::expected<Curve, std::string> readCurve(TokenReader &reader) {
            if (auto expected = reader.expect("curve_interval"); !expected) {
                return std::unexpected(expected.error());
            }

            auto fromDistance = reader.readDouble();
            if (!fromDistance) {
                return std::unexpected(fromDistance.error());
            }

            auto toDistance = reader.readDouble();
            if (!toDistance) {
                return std::unexpected(toDistance.error());
            }

            if (*fromDistance < 0.0 || *toDistance <= *fromDistance) {
                return std::unexpected(reader.error("curve interval must have increasing non-negative distances"));
            }

            if (auto expected = reader.expect("feed_count"); !expected) {
                return std::unexpected(expected.error());
            }

            auto feedCount = reader.readCount();
            if (!feedCount) {
                return std::unexpected(feedCount.error());
            }

            if (*feedCount == 0) {
                return std::unexpected(reader.error("every curve requires at least one feed"));
            }

            std::vector<double> feeds;
            feeds.reserve(*feedCount);
            for (std::size_t index = 0; index < *feedCount; ++index) {
                if (auto expected = reader.expect("feed"); !expected) {
                    return std::unexpected(expected.error());
                }

                auto feed = reader.readDouble();
                if (!feed) {
                    return std::unexpected(feed.error());
                }

                if (*feed <= 0.0) {
                    return std::unexpected(reader.error("feeds must be positive"));
                }

                feeds.push_back(*feed);
            }

            if (auto expected = reader.expect("curve"); !expected) {
                return std::unexpected(expected.error());
            }

            auto type = reader.read();
            if (!type) {
                return std::unexpected(type.error());
            }

            std::expected<CurveGeometry, std::string> geometry = std::unexpected("unknown curve type");
            if (*type == "line") {
                geometry = readLine(reader);
            } else if (*type == "arc") {
                geometry = readArc(reader);
            } else if (*type == "bspline") {
                geometry = readBSpline(reader);
            } else {
                return std::unexpected(reader.error("unknown curve type '" + *type + "'"));
            }

            if (!geometry) {
                return std::unexpected(geometry.error());
            }

            if (auto expected = reader.expect("end_curve"); !expected) {
                return std::unexpected(expected.error());
            }

            return Curve {*fromDistance, *toDistance, std::move(feeds), std::move(*geometry)};
        }
    }

    std::expected<GeometryFile, std::string> loadGeometryFile(const std::filesystem::path &path) {
        TokenReader reader(path);
        if (!reader.open()) {
            return std::unexpected("failed to open geometry file " + path.string());
        }

        if (auto expected = reader.expect("ngc_g64_geometry"); !expected) {
            return std::unexpected(expected.error());
        }

        auto version = reader.readCount();
        if (!version) {
            return std::unexpected(version.error());
        }

        if (*version != 1) {
            return std::unexpected(reader.error("unsupported geometry format version " + std::to_string(*version)));
        }

        if (auto expected = reader.expect("units"); !expected) {
            return std::unexpected(expected.error());
        }

        auto unitToken = reader.read();
        if (!unitToken) {
            return std::unexpected(unitToken.error());
        }

        Unit unit;
        if (*unitToken == "inch") {
            unit = Unit::Inch;
        } else if (*unitToken == "millimeter") {
            unit = Unit::Millimeter;
        } else {
            return std::unexpected(reader.error("unknown unit '" + *unitToken + "'"));
        }

        if (auto expected = reader.expect("curve_count"); !expected) {
            return std::unexpected(expected.error());
        }

        auto curveCount = reader.readCount();
        if (!curveCount) {
            return std::unexpected(curveCount.error());
        }

        GeometryFile result {.unit = unit, .curves = {}};
        result.curves.reserve(*curveCount);
        for (std::size_t index = 0; index < *curveCount; ++index) {
            auto curve = readCurve(reader);
            if (!curve) {
                return std::unexpected(curve.error());
            }

            result.curves.push_back(std::move(*curve));
        }

        if (auto expected = reader.expect("end_geometry"); !expected) {
            return std::unexpected(expected.error());
        }

        if (auto ended = reader.expectEnd(); !ended) {
            return std::unexpected(ended.error());
        }

        return result;
    }
}
