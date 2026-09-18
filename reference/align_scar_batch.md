# Batch SVD alignment for all specimens

A convenience wrapper around `align_scar()` that accepts a complete
multi-specimen data frame and handles the grouping internally.
Equivalent to calling `group_by(ID) %>% group_modify(~ align_scar(.x))`.

## Usage

``` r
align_scar_batch(data, id_col = "ID")
```

## Arguments

- data:

  A data frame containing **all** specimens. Must contain columns `ID`,
  `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`.

- id_col:

  Character. Name of the specimen ID column (default `"ID"`).

## Value

The input data frame with nine additional columns:

- s_x, s_y, s_z:

  Aligned start-point coordinates.

- e_x, e_y, e_z:

  Aligned end-point coordinates.

- d_x, d_y, d_z:

  Aligned unit direction vectors.

## Details

Specimens are processed independently: the SVD plane and in-plane
rotation are estimated separately for each group defined by `id_col`.
Row order within each group is preserved; the output is ungrouped.

## See also

[`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
aligned <- align_scar_batch(raw_data)
head(aligned[, c("ID", "d_x", "d_y", "d_z")])
} # }
```
