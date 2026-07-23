# PathTempo Repository Guide

PathTempo is a C++23 library for jerk-limited timing of fixed geometric paths. Keep the public API independent of CNC, G-code, and Ruckig types. Ruckig is a private implementation dependency.

Do not override CMake's configuration-specific optimization flags. Preserve the toolchain defaults for Debug, Release, RelWithDebInfo, and MinSizeRel, and do not leak warning or debug-symbol flags into dependencies.

Follow [docs/cpp_style.md](docs/cpp_style.md) for all new and modified C++ code. Preserve unrelated changes and add focused tests for every behavioral change.

The core timing model uses scalar path position, velocity, acceleration, and jerk. A path piece supplies its length, maximum velocity, and ordered differential constraint stations containing the arc-length derivatives `q'`, `q''`, and `q'''`. Geometry helpers produce one path piece for a line or helical arc and one for each non-empty B-spline or positive-weight NURBS knot interval. Any positive spline degree is accepted; clamping and uniform knot spacing are not constrained.

Sampled coupled-constraint checking evaluates every supplied differential station and both differential stations bracketing each scalar jerk-phase endpoint. Preserve the phase-endpoint checks so a short scalar phase cannot escape merely because it contains no station; these discrete checks do not establish continuous geometric proof between stations.

The primary solver output is a sequence of scalar time-domain cubic path-position polynomials. Coordinate-polynomial materialization and application-specific packetization remain separate layers.
