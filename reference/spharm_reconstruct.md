# Inverse spherical harmonic transform: coefficients to spherical density

**Low-level utility.** For interactive 3D exploration of reconstructed
surfaces,
[`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md)
is usually more convenient. Use this function when you need the
reconstructed density as numerical data — for example, computing
point-wise differences between specimens, extracting density peaks, or
building custom visualisations with `plotly` or `ggplot2`.

## Usage

``` r
spharm_reconstruct(coefficients, grid_size = 64)
```

## Arguments

- coefficients:

  Numeric array of shape `(2, lmax+1, lmax+1)` — the `coefficients`
  element of one specimen's `spharm_from_*()` result (real coefficients
  in pyshtools 4pi normalization).

- grid_size:

  Integer. Latitude resolution of the reconstruction grid; longitude
  resolution is `2 * grid_size`. Default 64. Larger values give smoother
  visualisations at higher cost.

## Value

A list with:

- density:

  Numeric matrix of shape `(grid_size, 2 * grid_size)` - the
  reconstructed density on the Driscoll-Healy grid. Negative values from
  finite-degree truncation are clipped to zero.

- lon:

  Numeric vector of length `2 * grid_size`. Longitude (radians, range
  `[0, 2*pi)`).

- lat:

  Numeric vector of length `grid_size`. Latitude (radians; the first row
  is the north pole, `pi/2`, and the south pole is excluded; equator =
  0).

- xyz:

  Numeric matrix of shape `(grid_size * 2 * grid_size, 3)` - unit-sphere
  Cartesian coordinates for each grid cell, in the same (column-major)
  order as `as.vector(density)`. Convenient input for
  `plotly::plot_ly(type = "surface", surfacecolor = ...)`.

## Details

Given a single specimen's spherical-harmonic coefficient array (as
returned in `spharm_from_*()$<id>$coefficients`), reconstructs the
corresponding density on a Driscoll-Healy grid on the unit sphere. The
returned object also includes the longitude / latitude vectors and an
`(n_lat * n_lon, 3)` matrix of unit-sphere Cartesian coordinates, for
direct plotting with
[`plotly::plot_ly()`](https://rdrr.io/pkg/plotly/man/plot_ly.html).

## References

Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working with
spherical harmonics. *Geochemistry, Geophysics, Geosystems*, **19**(8),
2574–2592.

## See also

[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md),
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md),
[`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- spharm_from_directions(my_aligned_data, lmax = 20)

# Reconstruct the spherical density for one specimen
rec <- spharm_reconstruct(result[[1]]$coefficients, grid_size = 64)

# Quick plotly visualisation
library(plotly)
nlat <- length(rec$lat); nlon <- length(rec$lon)
x_mat <- matrix(rec$xyz[, "x"], nlat, nlon)
y_mat <- matrix(rec$xyz[, "y"], nlat, nlon)
z_mat <- matrix(rec$xyz[, "z"], nlat, nlon)
plot_ly(x = x_mat, y = y_mat, z = z_mat,
        surfacecolor = rec$density,
        type = "surface", colorscale = "Hot",
        showscale = TRUE)

# Batch reconstruction with lapply()
all_recs <- lapply(result, function(s) spharm_reconstruct(s$coefficients))
} # }
```
