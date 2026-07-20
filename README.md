# PathTempo

Jerk-limited path timing for fixed geometric paths, with coupled velocity, acceleration, and jerk constraints and cubic time-law output.

PathTempo is an early-stage C++23 library. Its public API is being extracted from a CNC trajectory planner and is not yet stable.

## Design

PathTempo treats geometry as ordered timing pieces. Each piece supplies its arc length, programmed velocity, and differential constraint stations containing `q'`, `q''`, and the full `q'''`. Lines and arcs normally create one piece. Every non-empty cubic or quintic spline knot interval creates one piece.

The planner operates on scalar path velocity, acceleration, and jerk while enforcing coupled aggregate and per-coordinate limits. Its primary output is a sequence of scalar path-position cubic polynomials in physical time.

See [docs/design.md](docs/design.md) for the current architecture.

## Dependencies

PathTempo uses unmodified upstream releases:

- HiGHS v1.14.0
- Ruckig v0.19.4

They are pinned Git submodules for reproducible source builds and remain private implementation dependencies. Initialize them with:

```powershell
git submodule update --init --recursive
```

Installed CMake packages can be used instead by configuring with `-DPATH_TEMPO_USE_BUNDLED_DEPENDENCIES=OFF`.

## Build

```powershell
cmake -S . -B build -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_C_COMPILER=clang `
  -DCMAKE_CXX_COMPILER=clang++
cmake --build build
ctest --test-dir build --output-on-failure
```

## CMake consumption

From a source checkout:

```cmake
add_subdirectory(path_tempo)
target_link_libraries(application PRIVATE PathTempo::PathTempo)
```

From an installed package:

```cmake
find_package(PathTempo CONFIG REQUIRED)
target_link_libraries(application PRIVATE PathTempo::PathTempo)
```

## License

PathTempo is available under the MIT License. HiGHS and Ruckig retain their own copyright and license notices.
