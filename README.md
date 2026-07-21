# PathTempo

C++23 building blocks for jerk-limited timing of fixed geometric paths.

PathTempo is early-stage software and its public API is not yet stable.

## Current capabilities

- fixed-distance scalar transitions with specified beginning and ending velocity and acceleration;
- fixed-capacity physical-time cubic output for scalar transitions;
- continuous multi-piece timing for sampled C2 paths;
- jerk-limited forward/backward velocity reachability;
- HiGHS-backed sequential coupled acceleration and full jerk refinement with exact Ruckig line-search acceptance;
- bounded local correction passes for sampled coupled acceleration and full jerk constraints;
- caller-provided materialization correction with bounded re-solves;
- globally positioned scalar cubic output preserving path-piece identity;
- monotonicity and direction-reversal validation;
- velocity-transition distance and reachability calculations;
- line and helical-arc sampling into differential constraint stations;
- arbitrary-degree B-spline and positive-weight NURBS sampling, producing one path piece per non-empty knot interval;
- caller-selected station counts for arcs and spline knot intervals;
- sampling of caller-owned arc-length curves without replacing their geometry representation; and
- a solver-neutral sparse linear-program interface with persistent HiGHS model updates, basis reuse, resource-limit classification, and primal extraction.

## Planned capabilities

- continuous geometric proof between differential stations.

## Design

PathTempo's path model represents geometry as ordered timing pieces. Each piece supplies its arc length, maximum velocity, and differential constraint stations containing `q'`, `q''`, and the full `q'''`. The geometry helpers create one piece for a line or helical arc and one piece for each non-empty B-spline or NURBS knot interval. Any positive spline degree is accepted; clamping and knot spacing are unrestricted. The non-rational B-spline overload is equivalent to supplying unit NURBS weights.

`PathPlanner` solves ordered pieces with tangent- and curvature-continuous boundaries. It applies each piece's maximum velocity plus aggregate and per-coordinate limits, constructs a jerk-limited forward/backward velocity reference, and materializes every piece through Ruckig. Its HiGHS sequential program uses station velocity, acceleration, one-sided scalar jerk, and acceleration-deviation variables with linearized coupled endpoint constraints. Station-by-station line search accepts internal state changes only when adjacent Ruckig transitions remain feasible and no slower. Curved candidates are then checked at every supplied differential station using complete coupled acceleration and jerk. Violating pieces receive a bounded local time-scale correction and are re-solved. The result reports the effective velocity, acceleration, and jerk limits used for each piece after those corrections.

These checks guarantee the supplied differential stations. They do not prove continuous geometric behavior between stations; callers must provide sufficient sampling or perform a later continuous/materialized proof and request correction.

`PathPlanningSettings::applySampledCorrections` defaults to `true`. A caller may disable it only when a stronger materialization proof checks the complete emitted coordinate polynomials and re-solves with per-piece corrections. PathTempo caches exact-key scalar line-search candidates, including failures, and rematerializes accepted cached durations through Ruckig before publishing them.

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

## Adaptive-pockets example

The adaptive-pockets example loads exported line, endpoint-exact arc, cubic B-spline, and quintic B-spline geometry without parsing G-code. It constructs sampled path pieces with the geometry helpers and solves the complete three-coordinate path using the serialized programmed feeds. The original `examples/data/adaptive_pockets.ngc` file is included only as a reference; the example reads `examples/data/adaptive_pockets.1.txt` and does not parse G-code.

```powershell
cmake --build build --target path_tempo_adaptive_pockets_example
.\build\path_tempo_adaptive_pockets_example.exe
```

The executable uses `examples/data/adaptive_pockets.1.txt` by default. A different geometry file can be passed as its only argument. The full example intentionally remains outside the default build and test suite because it is a large planning problem.

The example overrides these planning settings:

- `linearSolveTimeLimit = 0.5` overrides the 0.25-second default and sets the HiGHS wall-time limit for each individual sequential linear-program solve. It is not a limit on the complete `PathPlanner::solve` call. If HiGHS reaches this limit, PathTempo discards the non-optimal linear result, records `PlanningResourceLimit::Time` in the diagnostics, stops further sequential refinements for the current correction pass, and continues with the last Ruckig-feasible timing candidate. Sampled and caller-provided materialization checks still run, so the plan can succeed and a later correction pass can perform another independently limited linear solve. Multiple linear solves, correction passes, transition construction, and checking can therefore make total planning time substantially longer than 0.5 seconds. The solver's time limit is not a real-time deadline and may be exceeded slightly before HiGHS observes it.
- `simplexIterationLimit = 1'000'000` overrides the default of 4096 and limits simplex iterations in each individual HiGHS solve. Reaching it is handled like the time limit: the non-optimal linear result is not used, the resource-limit diagnostics are updated, and planning continues from the last feasible candidate.
- `sequentialIterations = 1` overrides the default of 4 and permits at most one HiGHS sequential-refinement solve in each correction pass. It does not control sample checking or the number of correction passes. A value of zero skips sequential linear refinement.

Sample correction and re-solving are controlled separately. `applySampledCorrections` defaults to `true`; after each candidate is built, PathTempo checks the supplied differential stations and requests per-piece time scaling for sampled coupled-limit violations. Set it to `false` to skip those checks. `maximumCorrectionPasses`, which defaults to 8, bounds the candidate-and-correction loop shared by sampled corrections and the optional `MaterializationCorrection` callback. It must be positive. A value of one does not disable correction: if the first candidate requests a correction, PathTempo applies the tighter limits but then reports non-convergence because no pass remains to solve them. `diagnostics.correctionPasses` counts candidate passes, so a value of one means that the initial candidate required no additional correction solve.

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

Installed-package consumption currently requires compatible HiGHS and Ruckig CMake packages to be installed separately.

## License

PathTempo is available under the MIT License. HiGHS and Ruckig retain their own copyright and license notices.
