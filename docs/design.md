# PathTempo Design

## Scope

PathTempo is designed to compute jerk-limited scalar timing for a fixed geometric path without exposing solver dependency types.

The target architecture includes:

- generic fixed-degree-of-freedom path timing types;
- differential constraint stations;
- common geometry sampling helpers;
- forward and backward reachability;
- HiGHS-backed sequential convex programming;
- Ruckig-backed local state transitions;
- correction-pass state and reusable solver workspaces; and
- scalar time-domain cubic path-position output.

## Path pieces

A path piece is one timing interval. Lines and arcs normally produce one path piece. Each non-empty knot interval of a cubic or quintic spline produces one path piece. Increasing sampling density inside a piece must not manufacture new semantic piece boundaries.

Each differential station contains local arc distance and the geometry-only coefficients needed to map scalar path motion into coordinate motion:

```text
coordinate velocity     = q' v
coordinate acceleration = q' a + q'' v^2
coordinate jerk         = q' j + 3 q'' v a + q''' v^3
```

The third derivative is the full arc-length derivative `q'''`, including its tangential component. It is not normal sharpness.

Fixed samples support initial limits, solver constraints, caching, and diagnostics. They are not by themselves proof of continuous behavior between stations.

## Boundary state

A planning window begins and ends with scalar velocity and acceleration. Position is implied by the beginning and end of the supplied path. Jerk is a piecewise-constant control, not a boundary state, so jerk may change across a C2 window boundary.

## Output and correction

The multi-piece output is a sequence of cubic scalar path-position segments in physical time. Every segment belongs to exactly one path piece.

Coordinate materialization may discover an exact polynomial constraint violation. The optional `MaterializationCorrection` callback receives the complete candidate and returns required time-scale corrections for owning path pieces. PathTempo validates piece identity and scale values, combines those requests with sampled corrections, tightens local velocity, acceleration, and jerk limits, and re-solves while retaining reusable HiGHS and Ruckig state. Geometry-tolerance failures that can be repaired by subdividing coordinate cubics do not require a scalar re-solve.

## Current implementation

`ScalarTransitionPlanner` calculates one fixed-distance transition between scalar velocity/acceleration boundary states, validates monotonic forward motion, and returns its constant-jerk phases as physical-time cubic path-position segments. Its fixed-capacity result and reusable private workspace avoid per-call allocation.

`PathPlanner` calculates a continuous scalar time law across ordered C2 pieces. It constructs a velocity-cap reference with jerk-limited forward and backward reachability and uses `ScalarTransitionPlanner` to materialize each piece. HiGHS then proposes internal boundary velocity and acceleration changes using a local-duration objective, acceleration-deviation penalty, trust regions, and linearized coupled endpoint acceleration and full jerk constraints. One-sided scalar jerk variables establish endpoint jerk feasibility without making jerk a boundary state. A station-local working-set line search accepts a proposal only when both adjacent Ruckig transitions remain feasible and their combined duration does not increase. The output cubics use global path distance and retain their owning piece IDs.

For curved pieces, PathTempo evaluates `q' v`, `q' a + q'' v^2`, and `q' j + 3 q'' v a + q''' v^3` at every supplied differential station and on both scalar phases when a station coincides with a phase boundary. A violation produces a per-piece time-scale correction that tightens scalar velocity, acceleration, and jerk by the corresponding first, second, and third powers before a bounded re-solve.

`PersistentLinearSolver` owns HiGHS model storage, structure-stable updates, basis reuse, resource-limit classification, and primal extraction. Callers build a solver-neutral row-wise `SparseLinearProgram`; no HiGHS type crosses the public boundary.

Velocity-transition distance and reachable-velocity helpers are implemented. Differential stations remain discrete constraints, not continuous proof between samples. PathTempo provides the correction callback but leaves coordinate-polynomial materialization and continuous proof to the caller.

## Dependencies

HiGHS and Ruckig are unmodified upstream dependencies and remain private to the implementation. Source builds use pinned Git submodules by default. Parent-provided or installed CMake package targets are also supported.
