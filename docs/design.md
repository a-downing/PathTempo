# PathTempo Design

## Scope

PathTempo computes jerk-limited scalar timing for a fixed geometric path. It does not interpret G-code, construct CNC blends, own presentation metadata, packetize trajectories for a motion backend, or expose solver dependency types.

The library owns:

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

The primary output is a sequence of cubic scalar path-position segments in physical time. Every segment belongs to exactly one path piece.

Coordinate materialization may discover an exact polynomial constraint violation. A caller returns a required time-scale correction for the owning path piece. PathTempo tightens its local velocity, acceleration, and jerk limits and re-solves while retaining reusable HiGHS and Ruckig state. Geometry-tolerance failures that can be repaired by subdividing coordinate cubics do not require a scalar re-solve.

## Dependencies

HiGHS and Ruckig are unmodified upstream dependencies and remain private to the implementation. Source builds use pinned Git submodules by default. Parent-provided or installed CMake package targets are also supported.
