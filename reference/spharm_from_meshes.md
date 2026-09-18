# Spherical harmonic analysis from STL meshes

Pipeline (morphological analysis), per STL file:

1.  Read mesh via `open3d`. Files exceeding 3M faces are pre-decimated
    by streaming face sub-sampling to avoid memory pressure.

2.  Quadric-error decimation to `target_faces`.

3.  Laplacian smoothing (`smooth_iterations` passes).

4.  Volume-centroid normalization to a unit sphere.

5.  Area-weighted PCA alignment with energy-based sign convention.

6.  Cartesian-to-spherical conversion of the aligned vertices.

7.  Interpolation onto a regular `(grid_size, grid_size)` lat-lon grid.

8.  Spherical harmonic expansion via `pyshtools` (4pi normalization,
    zero-component normalized: all coefficients divided by `c(0, 0)`).

## Usage

``` r
spharm_from_meshes(
  stl_dir,
  lmax = 20,
  target_faces = 20000,
  grid_size = 256,
  smooth_iterations = 3,
  pre_decimate_threshold = 3e+06,
  pre_decimate_target = 5e+05,
  verbose = TRUE
)
```

## Arguments

- stl_dir:

  Character. Path to a directory containing `.stl` files. Each STL is
  treated as one specimen; the file basename (without extension) becomes
  its ID.

- lmax:

  Integer. Maximum spherical harmonic degree. Default 20.

- target_faces:

  Integer. Decimation target. Default 20000.

- grid_size:

  Integer. Latitude resolution of the interpolation grid; the longitude
  direction uses `2 * grid_size` points (Driscoll-Healy sampling 2).
  Default 256.

- smooth_iterations:

  Integer. Laplacian smoothing iterations after decimation. Default 3.
  Set to 0 to skip smoothing.

- pre_decimate_threshold:

  Integer. Face count above which streaming pre-decimation is triggered.
  Default 3000000.

- pre_decimate_target:

  Integer. Pre-decimation target face count. Default 500000.

- verbose:

  Logical. Print per-file progress. Default `TRUE`.

## Value

A list with one element per successfully-processed STL, named by file
basename. Each element is a list with `coefficients` and
`power_spectrum` (same structure as
[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)).
Failed specimens are reported via warnings and omitted from the result.

## Details

Requires the **mesh** Python extension. Run
`install_spharmlithic_python(mesh = TRUE)` once before first use.

**macOS:** not supported in a native installation. `open3d` ships its
own copy of the OpenMP runtime (`libomp`), which clashes with the copy
that comes with CRAN's R for macOS (used by R packages such as
`data.table`), and the R session aborts with `OMP: Error #15`. Use the
Docker image instead (see
[`vignette("Introduction", package = "spharmlithic")`](https://peiyuanxiao.github.io/spharmlithic/articles/Introduction.md)).
[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)
is not affected.

## References

Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working with
spherical harmonics. *Geochemistry, Geophysics, Geosystems*, **19**(8),
2574–2592.

## See also

[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md),
[`spharm_to_dataframe()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_to_dataframe.md),
[`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md),
[`install_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/install_spharmlithic_python.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires mesh extension
install_spharmlithic_python(mesh = TRUE)

result <- spharm_from_meshes(
  stl_dir      = "data/3D_models",
  lmax         = 20,
  target_faces = 20000
)

df <- spharm_to_dataframe(result)
write.csv(df, "spharm_meshes.csv", row.names = FALSE)
} # }
```
