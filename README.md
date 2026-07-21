# PathTempo

C++23 building blocks for jerk-limited timing of fixed geometric paths.

PathTempo is early-stage software and its public API is not yet stable.

## Current capabilities

- fixed-distance scalar transitions with specified beginning and ending velocity and acceleration;
- fixed-capacity physical-time cubic output for scalar transitions;
- continuous multi-piece timing for sampled C2 paths;
- HiGHS-optimized velocity envelopes with jerk-limited forward/backward reachability;
- HiGHS-backed sequential coupled-acceleration refinement with exact Ruckig line-search acceptance;
- bounded local correction passes for sampled coupled acceleration and full jerk constraints;
- globally positioned scalar cubic output preserving path-piece identity;
- monotonicity and direction-reversal validation;
- velocity-transition distance and reachability calculations;
- line sampling into differential constraint stations; and
- a solver-neutral sparse linear-program interface with persistent HiGHS model updates, basis reuse, resource-limit classification, and primal extraction.

## Planned capabilities

- arc, cubic-spline, and quintic-spline construction helpers;
- continuous geometric proof between differential stations; and
- materialization-driven correction from coordinate-polynomial extrema.

## Design

PathTempo's path model represents geometry as ordered timing pieces. Each piece supplies its arc length, programmed velocity, and differential constraint stations containing `q'`, `q''`, and the full `q'''`. Lines and arcs normally create one piece. Every non-empty cubic or quintic spline knot interval creates one piece.

`PathPlanner` solves ordered pieces with tangent- and curvature-continuous boundaries. It applies programmed velocity plus aggregate and per-coordinate limits, optimizes a squared-velocity envelope through HiGHS, tightens that envelope with jerk-limited reachability, and materializes every piece through Ruckig. Curved candidates are checked at every supplied differential station using complete coupled acceleration and jerk. Violating pieces receive a bounded local time-scale correction and are re-solved. A sequential HiGHS refinement proposes internal velocity and acceleration changes; station-by-station line search accepts them only when adjacent Ruckig transitions remain feasible, constraint-valid, and no slower.

These checks guarantee the supplied differential stations. They do not prove continuous geometric behavior between stations; callers must provide sufficient sampling or perform a later continuous/materialized proof and request correction.

```cpp
path_tempo::PathPlanner planner;
auto result = planner.solve(path_tempo::PathPlanningRequest<3> {
    .pieces = pieces,
    .beginning = {},
    .ending = {},
    .limits = limits,
});
```

## Scalar transitions

`ScalarTransitionPlanner` calculates a jerk-limited transition over a fixed path distance with specified beginning and ending velocity and acceleration. It owns a reusable Ruckig workspace and returns a fixed-capacity sequence of cubic path-position segments without exposing Ruckig types.

```cpp
path_tempo::ScalarTransitionPlanner planner;
auto transition = planner.solve({
    .piece = 1,
    .length = 10.0,
    .beginning = {},
    .ending = {},
    .maximumVelocity = 4.0,
    .maximumAcceleration = 3.0,
    .maximumJerk = 8.0,
});
```

The solver-neutral `SparseLinearProgram` and `PersistentLinearSolver` types similarly keep HiGHS model, basis, update, and resource-limit details behind PathTempo's implementation boundary.

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

Installed-package consumption currently requires compatible HiGHS and Ruckig CMake packages to be installed separately.

## License

PathTempo is available under the MIT License. HiGHS and Ruckig retain their own copyright and license notices.
