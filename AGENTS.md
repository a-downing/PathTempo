# PathTempo Repository Guide

PathTempo is a C++23 library for jerk-limited timing of fixed geometric paths. Keep the public API independent of CNC, G-code, HiGHS, and Ruckig types. HiGHS and Ruckig are private implementation dependencies.

Follow [docs/cpp_style.md](docs/cpp_style.md) for all new and modified C++ code. Preserve unrelated changes and add focused tests for every behavioral change.

The core timing model uses scalar path position, velocity, acceleration, and jerk. A path piece supplies its length, maximum velocity, and ordered differential constraint stations containing the arc-length derivatives `q'`, `q''`, and `q'''`. Lines and arcs normally produce one path piece. For a B-spline, callers normally produce one path piece for each non-empty knot interval; clamping and uniform knot spacing are not required.

The primary solver output is a sequence of scalar time-domain cubic path-position polynomials. Coordinate-polynomial materialization and application-specific packetization remain separate layers.
