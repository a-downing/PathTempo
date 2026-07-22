# PathTempo Design

## Scope

PathTempo is designed to compute jerk-limited scalar timing for a fixed geometric path without exposing solver dependency types.

The target architecture includes:

- generic fixed-degree-of-freedom path timing types;
- differential constraint stations;
- common geometry sampling helpers;
- forward and backward reachability;
- Ruckig-backed local state transitions;
- correction-pass state and reusable solver workspaces; and
- scalar time-domain cubic path-position output.

## Path pieces

A path piece is one timing interval. The geometry helpers produce one path piece for a line or helical arc and one path piece for each non-empty B-spline or positive-weight NURBS knot interval, regardless of degree, clamping, or uniformity. Increasing sampling density inside a piece must not manufacture new semantic piece boundaries.

Each differential station contains local arc distance and the geometry-only coefficients needed to map scalar path motion into coordinate motion:

```text
coordinate velocity     = q' v
coordinate acceleration = q' a + q'' v^2
coordinate jerk         = q' j + 3 q'' v a + q''' v^3
```

The third derivative is the full arc-length derivative `q'''`, including its tangential component. It is not normal sharpness.

Fixed samples support initial limits, solver constraints, and diagnostics. They are not by themselves proof of continuous behavior between stations.

## Boundary state

A planning window begins and ends with scalar velocity and acceleration. Position is implied by the beginning and end of the supplied path. Jerk is a piecewise-constant control, not a boundary state, so jerk may change across a C2 window boundary.

## Output and correction

The multi-piece output is a sequence of cubic scalar path-position segments in physical time. Every segment belongs to exactly one path piece.

Coordinate materialization may discover an exact polynomial constraint violation. The optional `MaterializationCorrection` callback receives the complete candidate and returns required time-scale corrections for owning path pieces. PathTempo validates piece identity and scale values, combines those requests with any enabled sampled corrections, tightens local velocity, acceleration, and jerk limits, and re-solves while retaining its reusable Ruckig state. Geometry-tolerance failures that can be repaired by subdividing coordinate cubics do not require a scalar re-solve.

## Current implementation

`ScalarTransitionPlanner` calculates one fixed-distance transition between scalar velocity/acceleration boundary states, validates monotonic forward motion, and returns its constant-jerk phases as physical-time cubic path-position segments. Its fixed-capacity result and reusable private workspace avoid per-call allocation.

`PathPlanner` calculates a continuous scalar time law across ordered C2 pieces. Both modes first solve and correct a conservative zero-acceleration boundary profile. Zero mode emits that profile without further boundary optimization. Optimized mode constructs a forward/backward envelope of squared boundary velocities, smooths the piece-average accelerations, and proposes continuous non-zero scalar acceleration at internal boundaries under the corrected limits. Each adjacent piece is materialized and checked by `ScalarTransitionPlanner`, and curved-path constraints are checked before acceptance. Infeasible proposals move locally toward the conservative profile. Alternating look-ahead sweeps revisit every internal station in both directions, first refining the proposed velocity and acceleration together and then refining each state component independently. A change is retained only when both adjacent pieces remain feasible and their combined duration decreases. The optimized proposal is used only when its complete duration is shorter than the conservative profile. The requested beginning and ending acceleration remain fixed. The output cubics use global path distance and retain their owning piece IDs.

The language-independent [path planner algorithm specification](path_planner_algorithm.md) defines the reachability formulas, proposal construction, repair process, look-ahead refinement, coupled checks, numerical constants, correction loop, and complete pseudocode required for an independent implementation.

Geometry helpers sample lines, helical arcs, arbitrary-degree B-splines, and positive-weight NURBS into arc-length differential stations. Arc and spline callers select the number of intervals sampled inside each timing piece. A separate arc-length-curve helper accepts caller-owned geometry evaluators, allowing applications to retain an authoritative curve and inversion policy while using the same station construction contract.

For curved pieces, PathTempo evaluates `q' v`, `q' a + q'' v^2`, and `q' j + 3 q'' v a + q''' v^3` at every supplied differential station and on both scalar phases when a station coincides with a phase boundary. A violation produces a per-piece time-scale correction that tightens scalar velocity, acceleration, and jerk by the corresponding first, second, and third powers before a bounded re-solve.

Applications with a stronger continuous coordinate-polynomial proof may disable sampled correction while retaining the differential stations for initial limits and sampled coupled checks performed by that stronger proof.

Velocity-transition distance and reachable-velocity helpers are implemented. Differential stations remain discrete constraints, not continuous proof between samples. PathTempo provides the correction callback but leaves coordinate-polynomial materialization and continuous proof to the caller.

## Dependencies

Ruckig is an unmodified upstream dependency and remains private to the implementation. Source builds use the pinned Git submodule by default. A parent-provided or installed CMake package target is also supported.
