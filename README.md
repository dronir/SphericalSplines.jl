
# SphericalSplines

Spherical spline interpolation.

## Usage

Given a number of direction vectors, `directions`, which is an `(N,3)` array, and a number
of real values associated with those directions, `values`, compute the iterpolated value
in `new_direction`:

```
using SphericalSplines

S = InterpolationSpline(directions, values)

y = S(new_direction)

```

The rows of `directions` need to be normalized to be unit vectors.


## Caveats

1. This is not (yet) well tested.
2. This is not (yet) optimized for performance.
