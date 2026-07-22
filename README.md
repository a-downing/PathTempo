# PathTempo

C++23 building blocks for jerk-limited timing of fixed geometric paths.

PathTempo is early-stage software and its public API is not yet stable.

## Current capabilities

- fixed-distance scalar transitions with specified beginning and ending velocity and acceleration;
- fixed-capacity physical-time cubic output for scalar transitions;
- continuous multi-piece timing for sampled C2 paths;
- jerk-limited forward/backward velocity reachability;
- Ruckig-backed materialization of jerk-limited piece transitions;
- bounded local correction passes for sampled coupled acceleration and full jerk constraints;
- caller-provided materialization correction with bounded re-solves;
- globally positioned scalar cubic output preserving path-piece identity;
- monotonicity and direction-reversal validation;
- velocity-transition distance and reachability calculations;
- line and helical-arc sampling into differential constraint stations;
- arbitrary-degree B-spline and positive-weight NURBS sampling, producing one path piece per non-empty knot interval;
- caller-selected station counts for arcs and spline knot intervals;
- sampling of caller-owned arc-length curves without replacing their geometry representation.

## Planned capabilities

- continuous geometric proof between differential stations.

## Design

PathTempo's path model represents geometry as ordered timing pieces. Each piece supplies its arc length, maximum velocity, and differential constraint stations containing `q'`, `q''`, and the full `q'''`. The geometry helpers create one piece for a line or helical arc and one piece for each non-empty B-spline or NURBS knot interval. Any positive spline degree is accepted; clamping and knot spacing are unrestricted. The non-rational B-spline overload is equivalent to supplying unit NURBS weights.

`PathPlanner` solves ordered pieces with tangent- and curvature-continuous boundaries. Both boundary-acceleration modes first construct and correct a conservative zero-acceleration boundary profile. `BoundaryAccelerationMode::Zero` emits that profile directly. `BoundaryAccelerationMode::Optimized`, the library default, proposes internal velocity and non-zero acceleration states with jerk-limited forward/backward reachability under the corrected limits. Every piece is materialized through Ruckig, and curved candidates are checked at every supplied differential station using complete coupled acceleration and jerk. Infeasible states are backed off locally toward the conservative profile. Alternating look-ahead sweeps optimize boundary velocity and acceleration against the combined duration of both adjacent pieces, retaining only exactly feasible improvements. The optimized result is used only when it shortens the complete trajectory. The requested beginning and ending accelerations are preserved. The result reports the effective velocity, acceleration, and jerk limits used for each piece after correction.

These checks guarantee the supplied differential stations. They do not prove continuous geometric behavior between stations; callers must provide sufficient sampling or perform a later continuous/materialized proof and request correction.

`PathPlanningSettings::applySampledCorrections` defaults to `true`. A caller may disable it only when a stronger materialization proof checks the complete emitted coordinate polynomials and re-solves with per-piece corrections.

The optional `MaterializationCorrection` callback receives each complete candidate `PlannedPath`. A caller may return per-piece required time scales discovered by coordinate-polynomial extrema or another stronger proof. PathTempo combines them with any enabled sampled corrections, tightens only the affected scalar limits, and repeats the same bounded solve.

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

See [docs/design.md](docs/design.md) for the current architecture.

## Dependencies

PathTempo uses the unmodified upstream Ruckig v0.19.4 release.

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

PathTempo does not override CMake's configuration-specific optimization flags. Debug retains the toolchain's standard unoptimized behavior, runtime checks, assertions, and symbols; Release, RelWithDebInfo, and MinSizeRel retain their toolchain defaults. Use RelWithDebInfo when profiling optimized code with debug symbols.

## Continuous-path geometry example

The continuous-path example loads exported line, endpoint-exact arc, cubic B-spline, and quintic B-spline geometry without parsing G-code. It constructs sampled path pieces with the geometry helpers and solves the complete three-coordinate path using the serialized programmed feeds. The source program is retained as `examples/data/continuous_path_reference.ngc` for reference; the example reads `examples/data/continuous_path.txt` and does not parse G-code. The complete version-1 token grammar and geometry requirements are documented in [examples/data/README.md](examples/data/README.md).

```powershell
cmake --build build --target path_tempo_continuous_path_example
.\build\path_tempo_continuous_path_example.exe -m optimized -u inch -a 5.1 -j 101 .\examples\data\continuous_path.txt
```

The geometry file, `-m <zero|optimized>`, `-u <inch|mm>`, `-a <acceleration>`, and `-j <jerk>` are required. Mode `zero` uses the fast zero-internal-acceleration profile; mode `optimized` enables non-zero boundary acceleration and look-ahead refinement. The selected unit must match the unit declared by the geometry file; acceleration and jerk are interpreted in that unit per second squared and cubed. Optional `-v <scale>` multiplies every serialized PathPiece velocity limit and defaults to `1`. The executable is part of the default build but is not run by the test suite because it is a large planning workload. The geometry text loader has a focused test target.

The example enables sampled corrections and allows up to 128 correction passes because changing its command-line limits can propagate coupled corrections across thousands of adjacent PathPieces. The library default for `maximumCorrectionPasses` remains 8. After each candidate is built, PathTempo checks the supplied differential stations and requests per-piece time scaling for sampled coupled-limit violations. Set `applySampledCorrections` to `false` to skip those checks. `maximumCorrectionPasses` bounds the candidate-and-correction loop shared by sampled corrections and the optional `MaterializationCorrection` callback and must be positive. A value of one does not disable correction: if the first candidate requests a correction, PathTempo applies the tighter limits but then reports non-convergence because no pass remains to solve them. `diagnostics.correctionPasses` counts candidate passes, so a value of one means that the initial candidate required no additional correction solve.

The optional `MaterializationCorrection` callback is independent of `applySampledCorrections`. Supplying it asks the caller to check each complete scalar timing candidate using a stronger application-specific proof and return per-piece time scales. Omitting the callback disables materialization correction.

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

Installed-package consumption requires a compatible Ruckig CMake package to be installed separately.

## License

PathTempo is available under the MIT License. Ruckig retains its own copyright and license notices.
