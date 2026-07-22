#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <string>

#include "geometry_text_loader.h"

#include "path_tempo/Sampling.h"

namespace path_tempo::example {
    class EndpointExactArc {
        static constexpr std::size_t LENGTH_TABLE_INTERVALS = 64;

        struct LengthNode {
            double parameter = 0.0;
            double distance = 0.0;
        };

        ArcGeometry m_source;
        Vector3 m_axis {};
        Vector3 m_startArm {};
        Vector3 m_endArm {};
        Vector3 m_startAxial {};
        Vector3 m_axialTravel {};
        double m_sweep = 0.0;
        double m_length = 0.0;
        std::array<LengthNode, LENGTH_TABLE_INTERVALS + 1> m_lengthNodes {};

        std::array<Vector<3>, 4> derivatives(double parameter) const;
        double speed(double parameter) const;
        double integratedLength(double from, double to) const;
        double parameterAtDistance(double distance) const;

    public:
        static std::expected<EndpointExactArc, std::string> create(const ArcGeometry &source);

        [[nodiscard]] double length() const;
        [[nodiscard]] Vector3 positionAtParameter(double parameter) const;
        [[nodiscard]] DifferentialState<3> stateAtDistance(double distance) const;
    };
}
