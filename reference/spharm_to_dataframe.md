# Flatten SPHARM result list into a wide-format data frame

Converts the per-specimen list returned by
[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)
or
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)
into a single wide-format data frame, with one row per specimen and one
column per power-spectrum value and per coefficient. Suitable for CSV
export and downstream multivariate analysis (PCA, UMAP, clustering) in
R.

## Usage

``` r
spharm_to_dataframe(x, include_coeffs = TRUE)
```

## Arguments

- x:

  A list returned by
  [`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)
  or
  [`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md).

- include_coeffs:

  Logical. If `TRUE` (default), include all flattened coefficients as
  columns named `coeff_0001`, `coeff_0002`, ... If `FALSE`, return only
  the power spectrum (much smaller).

## Value

A `tibble` with columns:

- ID:

  Specimen identifier.

- power_l0, power_l1, ..., power_lN:

  Raw power per spherical harmonic degree (N = lmax).

- coeff_0001, coeff_0002, ...:

  Flattened real-valued coefficients (only when
  `include_coeffs = TRUE`).

## See also

[`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md),
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md),
[`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md),
[`spharm_reconstruct()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_reconstruct.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- spharm_from_directions(my_aligned_data, lmax = 20)

# Full output (power + coefficients)
df_full <- spharm_to_dataframe(result)

# Power spectrum only (smaller, easier to inspect)
df_power <- spharm_to_dataframe(result, include_coeffs = FALSE)

write.csv(df_full, "spharm_results.csv", row.names = FALSE)
} # }
```
