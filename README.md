
# SphericalSplines

Spherical spline interpolation. Quick and dirty implementation, could possibly be optimized
a little.

Based on:

- [Blick, C. & Freeden, W. (2011). Spherical Spline Application to Radio Occultation Data.
  Journal of Geodetic Science, 1(4), pp. 379-395. Retrieved 20 Nov. 2018, from
  doi:10.2478/v10156-011-0017-8.][1]
  
- [Blick, C. (2012). Spherical Spline Approximation and its Application to Radio Occultation
  Data for Climate Monitoring. 2nd International Radio Occultation Workshop 2012.][2]

[1]: https://www.degruyter.com/view/j/jogs.2011.1.issue-4/v10156-011-0017-8/v10156-011-0017-8.xml
[2]: http://irowg.org/wpcms/wp-content/uploads/2013/12/blick.pdf

## Usage

### Interpolation spline

Given a number of direction vectors, `directions`, which is an `(N,3)` array, and a number
of real values associated with those directions, `values`, return an interpolating spline
solution (which passes through all the values in the given directions), and compute the
value in `new_direction`:

```
using SphericalSplines

S = InterpolationSpline(directions, values)

y = S(new_direction)

```

The rows of `directions` need to be normalized to be unit vectors.


### Smoothing spline

Given the direction vectors (as an `(N,3)` array) and values (as a `(N,)` vector) as above,
and additionally a `(N,)` shaped list of weights for each point and a smoothing
coefficient, return a smoothing spline solution and compute its value in `new_direction`:

```
using SphericalSplines

S = SmoothingSpline(directions, values, weights, smoothing_coef)

y = S(new_direction)

```

## Caveats

1. This is not (yet) well tested.
2. This is not (yet) optimized for performance.
