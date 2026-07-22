# Path Planner Algorithm

## Purpose

This document specifies the multi-piece scalar path-timing algorithm in enough detail to reproduce its behavior independently of any programming language or particular trajectory-generation library.

The algorithm times a fixed geometric path. It does not change the path geometry. Its output is a continuous scalar path state consisting of position, velocity, and acceleration, with piecewise-constant scalar jerk. Internal path-piece boundaries may carry non-zero velocity and non-zero acceleration. Jerk is not required to be continuous at a boundary.

Optimized mode is deliberately conservative:

1. Construct and correct a feasible profile whose internal boundary accelerations are zero.
2. Construct a faster proposal with continuous non-zero boundary accelerations.
3. Repair the proposal toward the conservative profile wherever it is infeasible.
4. Improve it with alternating local look-ahead sweeps.
5. Use it only if the complete proposal is faster than the conservative profile.

The result is not claimed to be a globally time-optimal solution. Every accepted local change is feasible under the checks described below and strictly reduces the duration of its two adjacent pieces.

## Boundary-acceleration modes

The planner exposes two modes that share input reduction, conservative reachability, sampled correction, output construction, and optional stronger materialization correction:

- `Zero` emits the corrected conservative profile directly. Every internal boundary acceleration is zero, and proposal construction and look-ahead refinement are skipped.
- `Optimized` continues from the corrected conservative profile through the non-zero proposal, repair, and refinement stages specified below. This is the default library mode.

The selected mode does not change the caller-specified beginning or ending acceleration.

## Terminology and indexing

For `n` path pieces:

- Piece `i` has index `0 .. n-1`.
- Boundary station `i` has index `0 .. n`.
- Piece `i` connects boundary station `i` to boundary station `i+1`.
- A boundary state is `(velocity, acceleration)`.
- `L[i]`, `V[i]`, `A[i]`, and `J[i]` are the length and effective scalar velocity, acceleration, and jerk limits of piece `i`.
- `cap[i]` is the scalar velocity cap at boundary station `i`.
- `w[i]` is squared scalar velocity at boundary station `i`.
- A differential constraint station is a geometric sample inside a piece. It is distinct from a boundary station.

The path begins at a caller-specified boundary state `beginning` and ends at a caller-specified boundary state `ending`. Those two states are fixed throughout planning.

## Required local transition primitive

The algorithm requires a one-dimensional fixed-distance transition solver:

```text
SolveTransition(length, beginningState, endingState, V, A, J)
```

The solver must either fail or return a minimum-duration, forward-only transition that:

- covers exactly `length`;
- begins and ends at the requested velocity and acceleration;
- never exceeds the scalar velocity, acceleration, or jerk limits;
- never reverses scalar path direction; and
- exposes constant-jerk phases so position, velocity, and acceleration can be evaluated at any local time.

The multi-piece algorithm treats this primitive as the authoritative test of scalar feasibility. A proposal inferred from reachability formulas is never accepted merely because those formulas say it is reachable.

## Numerical constants

The following values are part of the specified behavior:

| Name | Value | Purpose |
|---|---:|---|
| Numerical zero | `1e-15` | Negligible scalar values |
| Geometry component zero | `1e-15` | Negligible tangent components in scalar-limit quotients |
| Unit-tangent squared tolerance | `1e-9` | Validates arc-length parameterization |
| Tangent continuity tolerance | `1e-9` | Validates adjacent pieces |
| Curvature continuity tolerance | `1e-8` | Validates adjacent pieces |
| Station-distance absolute tolerance | `1e-12` | Validates endpoint coverage |
| Station-distance relative tolerance | `1e-10` | Validates endpoint coverage |
| Constraint-station lookup absolute tolerance | `1e-12` | Matches samples to scalar phases |
| Constraint-station lookup relative tolerance | `1e-10` | Matches samples to scalar phases |
| Constraint relative tolerance | `1e-10` | Limit comparison allowance |
| Fixed-end reachability tolerance | `1e-10` | Endpoint feasibility allowance |
| Reachability convergence tolerance | `1e-11` | Forward/backward envelope convergence |
| Absolute reachability distance margin | `1e-12` | Conservative length reserve |
| Relative reachability distance margin | `1e-6` | Conservative length reserve |
| Reachability bisections | `52` | Reachable-velocity inversion |
| Reachability extra passes | `8` | Fixed-point propagation reserve |
| Acceleration seed reserve | `1e-2` | Keeps proposal slopes inside acceleration limits |
| Slope smoothing sweeps | `64` | Maximum alternating smoothing sweeps |
| Slope smoothing tolerance | `1e-9` | Residual slope-rate violation |
| Slack tent fraction | `1/16` | Jerk-derived acceleration slack |
| Slack bound fraction | `1/32` | Slack near an acceleration bound |
| Repair weight factor | `0.65` | Local movement toward the conservative profile |
| Minimum repair weight | `1 / 2^20` | Weight below which the conservative state is used |
| Repair extra rounds | `24` | Propagation reserve |
| Joint refinement sweeps | `8` | Maximum paired velocity/acceleration sweeps |
| Independent refinement sweeps | `12` | Maximum component-wise sweeps |
| Refinement bisections | `6` | Trial resolution |
| Duration comparison tolerance | `1e-12` | Required local improvement margin |
| Joint sweep convergence | `1e-2` | Duration improvement over a forward/backward pair |
| Independent sweep convergence | `1e-6` | Duration improvement over a forward/backward pair |
| Coupled correction significance | `1e-9` | Ignores insignificant time scales |
| Coupled correction safety factor | `1.01` | Moves a corrected baseline inside the feasible region |
| Constraint-station Newton steps | `12` | Accelerates scalar position inversion |
| Constraint-station inversion cap | `52` | Leaves forty bisection fallbacks after Newton |

The maximum number of outer correction passes is a setting. Its default is `8`, and it must be positive.

## Input reduction to effective scalar limits

Each piece supplies differential constraint stations with arc-length derivatives `q1`, `q2`, and `q3`, corresponding to the first, second, and third coordinate derivatives with respect to scalar path distance.

Initialize each piece's scalar limits with:

```text
V = min(programmedVelocityLimit, initialScalarVelocityLimit)
A = min(pathAccelerationLimit, initialScalarAccelerationLimit)
J = min(pathJerkLimit, initialScalarJerkLimit)
```

At every differential station and for every coordinate, tighten the tangent-derived scalar limits whenever `abs(q1)` is larger than the geometry-component-zero tolerance:

```text
V = min(V, coordinateVelocityLimit / abs(q1))
A = min(A, coordinateAccelerationLimit / abs(q1))
J = min(J, coordinateJerkLimit / abs(q1))
```

Every nonzero curvature or third-derivative component is constraining, regardless of magnitude:

```text
V = min(V, sqrt(coordinateAccelerationLimit / abs(q2)))
V = min(V, cbrt(coordinateJerkLimit / abs(q3)))
```

Also tighten velocity using aggregate geometric magnitudes:

```text
V = min(V, sqrt(pathAccelerationLimit / norm(q2)))
V = min(V, cbrt(pathJerkLimit / norm(q3)))
```

These limits intentionally separate the derivative terms. They are safe seeds, not a proof of the complete coupled acceleration and jerk expressions. The coupled check later evaluates all terms together.

The effective `V`, `A`, and `J` of every piece must be finite and positive. Each station tangent must have squared norm within `1e-9` of one. Differential stations must cover the complete piece in nondecreasing distance order: the first distance must be within `1e-12` of zero, and the last must be within `max(1e-12, 1e-10 * pieceLength)` of the piece length. Adjacent tangent components must agree within `1e-9`, and adjacent curvature components within `1e-8`, because a single velocity and acceleration state is shared at their boundary.

Define boundary velocity caps:

```text
cap[0] = beginning.velocity
cap[n] = ending.velocity
cap[i] = min(V[i-1], V[i])        for 0 < i < n
```

## Complete coupled-constraint check

Given a materialized transition for one piece, evaluate every supplied differential constraint station. Find every constant-jerk phase whose position interval contains the station, allowing an absolute distance tolerance of `1e-12` and a relative tolerance of `1e-10 * pieceLength`. Evaluating every containing phase intentionally evaluates both sides when a station coincides with a phase boundary.

If every station has exactly zero curvature and exactly zero third derivative, the scalar caps already imply every coordinate and aggregate constraint. Skip the coupled station-time evaluation for that piece. Any nonzero curvature or third derivative, including a value below the tangent-component tolerance, requires the coupled station-time evaluation.

Invert the monotone scalar position polynomial inside the phase with safeguarded Newton iteration. Start from linear position interpolation, retain a valid time bracket, and accept a Newton step only when it remains strictly inside that bracket; otherwise bisect. Stop when position agrees within eight double-precision epsilons scaled by `max(1, abs(stationDistance))`. Try at most `12` Newton steps, followed when necessary by up to `40` pure bisections. At the resulting time, obtain scalar velocity `v`, acceleration `a`, and jerk `j`. Clamp a tiny negative scalar velocity to zero before evaluating the coordinate constraints.

For every coordinate, calculate:

```text
coordinateVelocity     = q1 * v
coordinateAcceleration = q1 * a + q2 * v^2
coordinateJerk         = q1 * j + 3 * q2 * v * a + q3 * v^3
```

Also calculate the Euclidean norms of the complete coordinate acceleration and jerk vectors.

Convert every violation into the time scale required to remove it:

```text
velocity scale     = measuredVelocity / velocityLimit
acceleration scale = sqrt(measuredAcceleration / accelerationLimit)
jerk scale         = cbrt(measuredJerk / jerkLimit)
```

The required scale for the piece is the maximum of `1` and all coordinate and aggregate scales. The powers follow uniform time scaling: velocity changes with the first power of time, acceleration with the second, and jerk with the third.

When sampled coupled corrections are enabled, a transition is accepted by a proposal or refinement only if its required scale is no greater than `1 + 1e-9`.

## Conservative zero-acceleration profile

The conservative profile is both a safe fallback and the lower endpoint of every proposal-repair homotopy. Its internal boundary accelerations are zero. The requested beginning and ending accelerations are preserved.

### Zero-to-zero transition distance estimate

For a velocity change `delta = abs(toVelocity - fromVelocity)` under limits `A` and `J`:

```text
if delta <= 1e-15:
    transitionTime = 0
else if delta <= A^2 / J:
    transitionTime = 2 * sqrt(delta / J)
else:
    transitionTime = delta / A + A / J

transitionDistance = midpoint(fromVelocity, toVelocity) * transitionTime
```

This is the minimum distance estimate for a symmetric jerk-limited velocity transition whose endpoint accelerations are zero. The conservative reachability seed uses this estimate even when a requested path-end acceleration is non-zero; materialization with the authoritative transition solver then checks the actual requested endpoint state.

### Reachable velocity

Reserve a small part of each piece length:

```text
availableLength = max(0, L - max(1e-12, L * 1e-6))
```

To find the highest velocity between a fixed velocity and a cap, first test the cap with the transition-distance estimate. If it does not fit, bisect the velocity interval `52` times and retain the highest fitting value.

### Forward/backward fixed point

Initialize every boundary velocity from its cap and overwrite the path endpoints with their fixed velocities. Alternate:

1. A forward pass that reduces boundary `i+1` to what is reachable from boundary `i` through piece `i`.
2. A backward pass that reduces boundary `i` to what can reach boundary `i+1` through piece `i`.

Do not alter the fixed endpoints. Instead, fail if an endpoint is unreachable by more than `1e-10`.

Stop when the largest velocity reduction in a complete forward/backward pair is no greater than `1e-11`, or after `2*n + 8` pairs. The resulting velocities are `floorVelocity[i]`.

Materialize every conservative piece with `SolveTransition`, using zero acceleration at internal boundaries. Failure here is a planning failure, because all later proposal states depend on this baseline.

### Correct the baseline before optimization

If sampled coupled corrections are enabled, calculate the required coupled time scale for every conservative transition. For each piece whose scale exceeds `1 + 1e-9`, multiply that scale by `1.01` and tighten only that piece:

```text
V = V / scale
A = A / scale^2
J = J / scale^3
```

Then restart the outer correction pass from boundary caps and rebuild the conservative profile. The non-zero proposal is not constructed until the conservative profile passes the sampled coupled checks under the current effective limits.

## Non-zero boundary-state proposal

### Boundary acceleration plan

Mark every internal boundary as `Free`. At the beginning or ending boundary, use `Zero` if the absolute requested acceleration is no greater than `1e-15`; otherwise use `Free`.

`Zero` and `Free` affect only the analytical reachability seed. The fixed endpoint acceleration itself is always passed unchanged to the authoritative transition solver.

### Reachability distance with free endpoint acceleration

Set the proposal acceleration bound for each piece to:

```text
seedA = A * (1 - 0.01)
```

For an increasing velocity transition from `v0` to `v1`, let `delta = v1 - v0`.

If both endpoint accelerations are free:

```text
distance = (v1^2 - v0^2) / (2 * seedA)
```

If both endpoint accelerations are zero, use the conservative zero-to-zero transition distance with `seedA` and `J`.

For a zero-acceleration beginning and free ending, define:

```text
rampChange = seedA^2 / (2 * J)
rampTime = seedA / J
```

Then:

```text
if delta >= rampChange:
    rampDistance = v0 * rampTime + seedA * rampTime^2 / 6
    holdTime = (delta - rampChange) / seedA
    holdEntryVelocity = v0 + rampChange
    distance = rampDistance
             + holdEntryVelocity * holdTime
             + seedA * holdTime^2 / 2
else:
    time = sqrt(2 * delta / J)
    distance = v0 * time + J * time^3 / 6
```

For a free beginning and zero-acceleration ending:

```text
if delta >= rampChange:
    holdTime = (delta - rampChange) / seedA
    holdDistance = v0 * holdTime + seedA * holdTime^2 / 2
    rampEntryVelocity = v0 + seedA * holdTime
    distance = holdDistance
             + rampEntryVelocity * rampTime
             + seedA * rampTime^2 / 3
else:
    peakAcceleration = sqrt(2 * delta * J)
    time = peakAcceleration / J
    distance = v0 * time
             + peakAcceleration * time^2 / 2
             - J * time^3 / 6
```

For a decreasing velocity transition, swap the two velocities and swap the beginning/end acceleration treatments. This mirrors the accelerating transition in time.

For cases that require inversion of these distance formulas, use the conservative available-length margin and `52` bisections from the conservative reachability calculation. The direct free-to-free squared-velocity update below intentionally uses the complete piece length.

### Squared-velocity envelope

Initialize `w[i] = cap[i]^2`, then overwrite `w[0]` and `w[n]` with the squared fixed endpoint velocities.

Alternate forward and backward reductions. If both ends of a piece are free, the direct forward limit is:

```text
reachable = min(sqrt(w[i] + 2 * seedA[i] * L[i]), cap[i+1])
```

The backward expression is symmetric. Otherwise, invert the applicable free/zero transition-distance formula with bisection.

Stop when the largest reduction in squared velocity is no greater than `1e-11`, or after `2*n + 8` forward/backward pairs. If this envelope reports a fixed endpoint as unreachable, use `floorVelocity[i]^2` as the proposal envelope instead of failing the complete plan.

## Squared-velocity smoothing

The derivative of squared velocity with respect to twice the path distance is the piece-average acceleration. For piece `i`:

```text
slope[i] = (w[i+1] - w[i]) / (2 * L[i])
crossingTime[i] = 2 * L[i] / (sqrt(w[i]) + sqrt(w[i+1]))
```

If the velocity sum is no greater than `1e-15`, treat the crossing time as infinite.

For a piece with average acceleration `slope` and crossing time `time`, define:

```text
tent = J * time
headroom = max(0, seedA - abs(slope))
slack = min(tent / 16, headroom + tent / 32)
```

At internal boundary `i`, calculate:

```text
leftSlope  = slope[i-1]
rightSlope = slope[i]
allowedRate = leftSlack + rightSlack
gap = rightSlope - leftSlope
violation = abs(gap) - allowedRate
```

If `violation <= 0`, do nothing. Otherwise:

```text
excess = gap - copySign(allowedRate, gap)
adjusted = w[i] + excess * (2 * leftLength * rightLength)
                          / (leftLength + rightLength)
```

Constrain the adjusted value to the energy window shared by both neighboring pieces:

```text
lower = max(
    0,
    w[i-1] - 2 * leftSeedA * leftLength,
    w[i+1] - 2 * rightSeedA * rightLength)

upper = min(
    cap[i]^2,
    w[i-1] + 2 * leftSeedA * leftLength,
    w[i+1] + 2 * rightSeedA * rightLength)
```

If `lower <= upper`, replace `w[i]` with `clamp(adjusted, lower, upper)`. Otherwise leave it unchanged.

Perform at most `64` Gauss-Seidel sweeps, visiting internal boundaries forward on even sweeps and backward on odd sweeps. Stop when the maximum pre-adjustment violation in a sweep is no greater than `1e-9`.

After smoothing, enforce the conservative floor:

```text
w[i] = max(w[i], floorVelocity[i]^2)       for every internal boundary
```

The authoritative transition solver and repair stage handle any feasibility lost by this final floor operation.

## Proposed boundary accelerations

For each internal boundary `i`, set `velocity = sqrt(w[i])`. If velocity is no greater than `1e-15`, set acceleration to zero.

Otherwise, recompute the left and right slopes, crossing times, and slacks. Define:

```text
sharedA = min(A[i-1], A[i])
lower = max(leftSlope - leftSlack,
            rightSlope - rightSlack,
            -sharedA)
upper = min(leftSlope + leftSlack,
            rightSlope + rightSlack,
            sharedA)
desired = (leftSlope + rightSlope) / 2
```

If the intervals overlap, clamp `desired` to `[lower, upper]`. If they do not overlap, clamp `desired` only to `[-sharedA, sharedA]` and let the exact feasibility stage repair it.

Apply velocity-cap sign restrictions:

```text
if velocity >= V[i-1] * (1 - 1e-10):
    acceleration = max(acceleration, 0)

if velocity >= V[i] * (1 - 1e-10):
    acceleration = min(acceleration, 0)
```

Finally, protect the peak velocity of the final jerk ramp in the left piece. If acceleration is negative:

```text
decelerationCap = sqrt(2 * J[i-1] * max(V[i-1] - velocity, 0))
acceleration = max(acceleration, -decelerationCap * (1 - 1e-10))
```

The resulting internal states are the full proposal `proposed[i]`.

## Local homotopy repair

Associate a weight with every internal boundary. Initialize internal weights to `1`; endpoint weights remain `0` because the endpoints are fixed.

For an internal boundary and weight `x`, interpolate from the conservative state to the proposal in squared-velocity space:

```text
velocity(x) = sqrt(lerp(floorVelocity^2, proposedVelocity^2, x))
acceleration(x) = proposedAcceleration * x
```

Materialize all pieces and mark a piece failed if:

- the local transition solver fails, or
- sampled coupled corrections are enabled and the piece's required scale exceeds `1 + 1e-9`.

While failed pieces remain:

1. Mark both internal boundaries of every failed piece for reduction.
2. For every marked boundary with non-zero weight, multiply its weight by `0.65`.
3. If the old weight is no greater than `(1 / 2^20) / 0.65`, set the new weight to zero.
4. Recompute only pieces adjacent to changed boundaries.
5. Repeat for at most `2*n + 24` rounds.

If no boundary can change, or the round limit is reached, replace the entire candidate with the conservative boundary states and transitions. This guarantees a known feasible starting point for refinement.

## Alternating look-ahead refinement

Every refinement at boundary `i` considers both piece `i-1` and piece `i`. Its local objective is:

```text
localDuration = duration(piece i-1) + duration(piece i)
```

A trial is eligible only if both adjacent transitions solve and, when sampled checks are enabled, both required coupled scales are no greater than `1 + 1e-9`. An eligible trial replaces the current state only if:

```text
trialDuration + 1e-12 < bestDuration
```

Consequently, committed changes never increase total path duration: all other pieces are unchanged, and the two changed pieces become strictly faster.

### Joint velocity/acceleration refinement

If the current boundary weight is below `1`, bisect the interval from the current weight to `1` six times. Each trial uses the same squared-velocity and acceleration interpolation as the repair stage.

- On failure, move the upper bisection endpoint down to the trial weight.
- On feasibility, move the lower endpoint up to the trial weight.
- Independently of feasibility interval movement, retain the trial only if it improves the local duration.

Perform at most `8` alternating sweeps. Visit internal boundaries forward on even sweeps and backward on odd sweeps. After every backward sweep, stop if the accumulated duration improvement across that forward/backward pair is no greater than `1e-2`; otherwise reset the accumulator for the next pair.

### Independent acceleration refinement

At each internal boundary, keep velocity fixed and test these six acceleration candidates, in this order:

```text
0
currentAcceleration / 2
midpoint(currentAcceleration, proposedAcceleration)
proposedAcceleration
clamp(currentAcceleration - sharedA / 4, -sharedA, sharedA)
clamp(currentAcceleration + sharedA / 4, -sharedA, sharedA)
```

Retain the feasible candidate with the lowest adjacent-piece duration, subject to the `1e-12` improvement margin.

### Independent velocity refinement

Keep acceleration fixed. If `proposedVelocity` is greater than `currentVelocity + 1e-15`, bisect that interval six times.

- On failure, lower the upper endpoint to the trial velocity.
- On feasibility, raise the lower endpoint to the trial velocity.
- Retain a feasible velocity only when it improves adjacent-piece duration by more than `1e-12`.

### Component-wise sweep schedule

At every visited boundary, refine acceleration first and velocity second. Perform at most `12` alternating forward/backward sweeps. After every backward sweep, stop if total improvement across the forward/backward pair is no greater than `1e-6`; otherwise reset the accumulator.

## Global acceptance and output

Sum the duration of all optimized transitions and all conservative transitions. If the optimized duration is greater than or equal to the conservative duration, discard it and use the conservative profile. Equality deliberately favors the simpler conservative result.

The selected local transitions form the scalar time law. Shift each piece's local position polynomial by the accumulated lengths of all preceding pieces. Preserve the owning piece identity on every constant-jerk segment. Boundary velocity and acceleration are shared exactly between adjacent pieces.

## Optional stronger materialization correction

After the complete candidate has been assembled, an optional stronger verifier may inspect the emitted time law, for example by finding exact extrema of fully materialized coordinate polynomials between differential samples.

The verifier may return one required time scale per affected piece. Reject invalid scales below `1`, non-finite scales, unknown piece identities, or duplicate corrections. For every significant scale, tighten the affected piece using:

```text
V = V / scale
A = A / scale^2
J = J / scale^3
```

Then restart the outer correction pass. Unlike the sampled-baseline correction, this externally supplied scale is not multiplied by the `1.01` safety factor.

Return the selected path when no significant correction remains. Fail if the configured maximum number of correction passes is exhausted.

## Top-level pseudocode

```text
function PlanPath(pieces, beginning, ending, limits, settings, strongerVerifier):
    validate geometry, limits, boundary states, and C2 piece continuity
    localPieces = reduce differential constraints to scalar V, A, and J
    correctedPieces = empty set

    for correctionPass in 0 .. settings.maximumCorrectionPasses - 1:
        caps = shared boundary velocity caps(localPieces, beginning, ending)
        validate fixed endpoint states against adjacent local limits

        floorVelocity = conservative forward/backward fixed point(caps)
        baselineStates = fixed endpoints plus zero-acceleration internal states
        baselineTransitions = SolveTransition for every piece

        if sampled coupled checks are enabled:
            sampledScales = coupled scales of baselineTransitions
            if any sampledScale > 1 + 1e-9:
                tighten affected pieces by sampledScale * 1.01
                add affected pieces to correctedPieces
                continue

        candidate = baselineTransitions

        if settings.boundaryAccelerationMode is Optimized:
            proposalSquaredVelocity = free-acceleration forward/backward envelope
            if that envelope cannot satisfy a fixed endpoint:
                proposalSquaredVelocity = floorVelocity^2

            smooth proposalSquaredVelocity with alternating slope-rate projections
            raise internal proposalSquaredVelocity to floorVelocity^2
            proposedStates = assign boundary accelerations from adjacent slopes

            candidate = proposedStates
            repair failed pieces by reducing adjacent boundary weights toward zero
            jointly refine boundary weights with alternating look-ahead sweeps
            independently refine acceleration and then velocity with alternating sweeps

        if duration(candidate) >= duration(baselineTransitions):
            candidate = baselineTransitions

        result = assemble globally positioned scalar cubic segments(candidate)

        if strongerVerifier exists:
            externalScales = strongerVerifier(result)
            validate externalScales
            if any externalScale > 1 + 1e-9:
                tighten affected pieces by externalScale
                add affected pieces to correctedPieces
                continue

        return result

    fail and identify whether sampled coupled-path correction, materialization correction, or both did not converge
```

## Invariants

A conforming implementation preserves these invariants:

- Beginning and ending velocity and acceleration never change.
- Every internal boundary has one shared velocity and acceleration for both adjacent pieces.
- Zero mode keeps every internal boundary acceleration at zero and skips all proposal refinement.
- Every accepted piece is feasible according to the authoritative fixed-distance transition solver.
- When sampled coupled checks are enabled, every accepted proposal transition passes them within `1e-9` relative time-scale significance.
- Refinement commits only strict local duration reductions.
- The final non-zero profile is used only if it is strictly faster than the conservative profile.
- Corrections affect only the owning pieces and restart the complete boundary-state solve.
- Differential stations constrain only the supplied samples; they do not prove continuous geometric feasibility between samples.
