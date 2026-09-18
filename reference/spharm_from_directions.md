# Spherical harmonic analysis from scar direction vectors

Pipeline (scar pattern analysis):

1.  Extract unit direction vectors `(d_x, d_y, d_z)` from each scar.

2.  Per-specimen von Mises-Fisher kernel density estimate on the sphere.

3.  Interpolate KDE values onto a Driscoll-Healy regular grid.

4.  Spherical harmonic expansion via `pyshtools` (4pi normalization).

## Usage

``` r
spharm_from_directions(
  data,
  lmax = 20,
  bandwidth = 0.35,
  n_bearing = 72,
  n_plunge = 36,
  dh_size = 64,
  id_col = "ID",
  dx_col = "d_x",
  dy_col = "d_y",
  dz_col = "d_z",
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame containing direction vectors. Must contain an ID column
  (default `"ID"`) and three direction-component columns (default
  `"d_x"`, `"d_y"`, `"d_z"`, matching the output of
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
  and
  [`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)).
  Rows where a direction vector is zero or NA are silently dropped.

- lmax:

  Integer. Maximum spherical harmonic degree. Default 20.

- bandwidth:

  Numeric. vMF KDE bandwidth (smaller -\> sharper peaks). Default 0.35.
  This is an empirical value; users are encouraged to conduct
  sensitivity analysis based on data characteristics (specimen scar
  count, directional concentration).

- n_bearing, n_plunge:

  Integers. Resolution of the intermediate evaluation grid for vMF KDE.
  Defaults 72 and 36 (5-degree spacing).

- dh_size:

  Integer. Driscoll-Healy grid latitude size; longitude is
  `2 * dh_size`. Default 64.

- id_col, dx_col, dy_col, dz_col:

  Column names. Defaults match
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
  output.

- verbose:

  Logical. Print per-specimen progress. Default `TRUE`.

## Value

A list with one element per specimen, named by ID. Each element is
itself a list with:

- coefficients:

  Numeric array of shape `(2, lmax+1, lmax+1)` - spherical harmonic
  coefficients in pyshtools 4pi normalization. Index `[1, l, m]` is the
  cosine coefficient (m \>= 0); index `[2, l, m]` is the sine
  coefficient (m \> 0).

- power_spectrum:

  Numeric vector of length `lmax+1` - raw power per degree.

Use
[`spharm_to_dataframe()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_to_dataframe.md)
to flatten the list into a wide-format data frame suitable for CSV
export.

## Details

This function calls the bundled Python backend via `reticulate`. Make
sure you have run
[`install_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/install_spharmlithic_python.md)
(or
[`use_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/use_spharmlithic_python.md))
at least once.

## References

Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working with
spherical harmonics. *Geochemistry, Geophysics, Geosystems*, **19**(8),
2574–2592.

## See also

[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md),
[`spharm_to_dataframe()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_to_dataframe.md),
[`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md),
[`install_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/install_spharmlithic_python.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Typical workflow
aligned <- align_scar_batch(my_scar_data)
result  <- spharm_from_directions(aligned, lmax = 20)

# Inspect one specimen
result[[1]]$power_spectrum

# Flatten for CSV / downstream analysis
df <- spharm_to_dataframe(result)
write.csv(df, "spharm_results.csv", row.names = FALSE)
} # }
```
