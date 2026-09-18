# Batch Lin 2024 alignment for all specimens

A convenience wrapper around `align_morph()` that accepts a complete
multi-specimen data frame and handles the grouping internally.

## Usage

``` r
align_morph_batch(data, id_col = "ID")
```

## Arguments

- data:

  A data frame containing **all** specimens. Must contain columns `ID`,
  `Start_X`, `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`, `Norm_X`,
  `Norm_Y`, `Norm_Z`, and optionally `Length`.

- id_col:

  Character. Name of the specimen ID column (default `"ID"`).

## Value

The input data frame with nine additional columns (`s_x/y/z`, `e_x/y/z`,
`d_x/y/z`).

## Details

Specimens are processed independently: the morphological plane normal
and the longest-scar anchor are determined separately for each group
defined by `id_col`. Row order within each group is preserved; the
output is ungrouped.

## See also

[`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
aligned <- align_morph_batch(raw_data)
head(aligned[, c("ID", "d_x", "d_y", "d_z")])
} # }
```
