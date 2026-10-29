# Fit the plane of intersection of a core to its outline

Fits a plane to the peripheral outline of a core by total least squares
and finds the centre of the region the outline encloses. The plane and
the centre are the reference for the Centripetality Index (CPI) and the
Plane-Relative Angle (PRA) of
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md).

## Usage

``` r
fit_intersection_plane(outline, K = 64)
```

## Arguments

- outline:

  The outline of one core: a numeric matrix or data frame with three
  columns (X, Y, Z), one point per row, **ordered along the loop**, as
  returned by
  [`read_outline_iges()`](https://peiyuanxiao.github.io/spharmlithic/reference/read_outline_iges.md).
  The loop is closed automatically if its last point differs from its
  first. A fitted `intersection_plane` object is returned unchanged.

- K:

  Number of equally spaced points (by arc length) the plane is fitted
  to. Default 64.

## Value

An object of class `intersection_plane`: a list with

- centre:

  The centroid of the area enclosed by the outline (the centre against
  which CPI is measured).

- normal:

  Unit normal of the fitted plane.

- plane_centroid:

  The centroid of the `K` resampled points (a point on the plane).

- basis:

  3 x 3 matrix whose rows are two in-plane axes and the normal
  (right-handed).

- flatness:

  The smallest singular value of the resampled outline over the middle
  one: 0 for a perfectly flat outline; the larger it is, the more the
  single plane is an approximation.

- rms_out_of_plane:

  Root mean square distance of the outline points from the plane.

- diagonal:

  Length of the outline's bounding-box diagonal.

- loop:

  The closed outline (first point repeated as the last).

## Details

The plane is fitted to `K` points spaced equally along the outline, so
that densely sampled stretches do not dominate the fit. The centre is
the centroid of the enclosed **area**, not of the outline points: the
loop is projected onto the fitted plane, the polygon centroid is taken
there, and the result is lifted back.

The points must follow the outline in order. If the projected outline
crosses itself, a warning is given: the points may be out of order, and
the centre would then be wrong.

## See also

[`read_outline_iges()`](https://peiyuanxiao.github.io/spharmlithic/reference/read_outline_iges.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md),
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md)

## Examples

``` r
# An ellipse tilted out of the XY plane
t <- seq(0, 2 * pi, length.out = 200)
loop <- cbind(30 * cos(t), 20 * sin(t), 3 * cos(t))
fit_intersection_plane(loop)
#> Warning: collapsing to unique 'x' values
#> Warning: collapsing to unique 'x' values
#> Warning: collapsing to unique 'x' values
#> <intersection_plane>
#>   centre (area centroid): 2.306e-15  8.783e-16  2.276e-16 
#>   normal:                 0.0995  2.204e-18  -0.995 
#>   flatness:               3.31e-17 
#>   rms out of plane:       5.575e-16 
#>   outline:                200 points, diagonal 72.36 
```
