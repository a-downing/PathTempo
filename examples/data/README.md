# Geometry Text Format

## Overview

The continuous-path example reads a strict, whitespace-delimited geometry format. Whitespace separates tokens, so indentation and line breaks are for readability only. Labels and symbolic values are case-sensitive. Comments, quoted strings, and trailing tokens are not supported.

All real values must be finite decimal numbers accepted by a standard floating-point parser. Counts must be non-negative decimal integers. Positions, distances, and feeds use the unit declared in the file; knots are parameter values and do not carry a distance unit.

Only format version `1` is supported.

## File structure

```text
ngc_g64_geometry 1
units <inch|millimeter>
curve_count <number-of-curves>

<curve record 1>
<curve record 2>
...
<curve record N>

end_geometry
```

`curve_count` must exactly match the number of curve records. No tokens may follow `end_geometry`.

The file uses `inch` or `millimeter`. The continuous-path executable spells the corresponding command-line values `-u inch` and `-u mm`; the command-line unit must match the file.

## Curve record

Every curve begins with its arc-length interval and feeds, followed by one geometry definition:

```text
curve_interval <from-distance> <to-distance>
feed_count <number-of-feeds>
feed <positive-feed>
...
curve <line|arc|bspline>
<geometry fields>
end_curve
```

Requirements:

- `from-distance` must be non-negative.
- `to-distance` must be greater than `from-distance`.
- `feed_count` must be positive.
- Exactly `feed_count` positive, finite `feed` values must follow.
- Feeds are PathPiece velocity limits in the file's distance unit per second.
- The executable's optional `-v` value multiplies these feeds after geometry sampling.

`curve_interval` identifies the portion of the source geometry used by this record. Distances are measured from the source geometry's beginning.

## Line

```text
curve line
from <x> <y> <z>
to <x> <y> <z>
```

Example:

```text
curve_interval 0 10
feed_count 1
feed 2.5
curve line
from 0 0 0
to 10 0 0
end_curve
```

The source line must have positive length. The interval may select a subsegment, but `to-distance` may not materially exceed the source line length. An excess within `1e-10 * max(1, source-length)` is treated as serialized roundoff and canonicalized to the exact source length. The continuous-path example requires exactly one feed and produces one PathPiece for the selected interval.

## Arc

```text
curve arc
from <x> <y> <z>
to <x> <y> <z>
center <x> <y> <z>
axis <x> <y> <z>
```

Example of a positive quarter-circle sweep about the Z axis:

```text
curve_interval 0 1.5707963267948966
feed_count 1
feed 2.5
curve arc
from 1 0 0
to 0 1 0
center 0 0 0
axis 0 0 1
end_curve
```

The axis must be finite and non-zero. The radial components of both endpoints relative to the center and axis must be non-zero. The axis establishes the positive sweep direction by the right-hand rule. Coincident start and end radial arms represent a complete positive revolution. A displacement along the axis produces a helical arc.

The interval may select a subinterval, but `to-distance` may not materially exceed the computed arc length. An excess within `1e-10 * max(1, arc-length)` is treated as serialized roundoff and canonicalized to the exact arc length. The continuous-path example requires exactly one feed and produces one PathPiece for the selected interval.

## B-spline

```text
curve bspline
degree <positive-degree>
control_count <number-of-controls>
control <x> <y> <z>
...
knot_count <number-of-knots>
knot <value>
...
```

Example of a quadratic, clamped, straight B-spline:

```text
curve_interval 0 2
feed_count 1
feed 2.5
curve bspline
degree 2
control_count 3
control 0 0 0
control 1 0 0
control 2 0 0
knot_count 6
knot 0
knot 0
knot 0
knot 1
knot 1
knot 1
end_curve
```

Requirements:

- `degree` must be positive.
- `control_count` must be greater than `degree`.
- Exactly `control_count` control points must follow.
- `knot_count` must equal `control_count + degree + 1`.
- Exactly `knot_count` finite, nondecreasing knot values must follow.
- Clamping and uniform knot spacing are not required.
- This format describes non-rational B-splines; it has no weight field.

The geometry sampler produces one PathPiece for every non-empty knot interval, in knot order. `feed_count` must equal that number, and feed `i` becomes the velocity limit of PathPiece `i`.

The continuous-path example currently accepts only a complete B-spline record: `from-distance` must be zero within `1e-10`, and `to-distance - from-distance` must match the sampled spline length within `1e-9 * max(1, interval-length)`. Partial B-spline intervals are rejected.

## Complete minimal file

```text
ngc_g64_geometry 1
units inch
curve_count 1
curve_interval 0 1
feed_count 1
feed 1.5
curve line
from 0 0 0
to 1 0 0
end_curve
end_geometry
```

It can be planned in zero-boundary mode with:

```powershell
.\path_tempo_continuous_path_example.exe -m zero -u inch -a 5 -j 100 geometry.txt
```

Use `-m optimized` to enable non-zero internal boundary acceleration and look-ahead refinement.
