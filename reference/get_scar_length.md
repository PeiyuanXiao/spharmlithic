# Get scar lengths from a data frame

Returns the length of each scar. If a `Length` column already exists it
is used directly; otherwise lengths are computed from the start- and
end-point coordinates.

## Usage

``` r
get_scar_length(df)
```

## Arguments

- df:

  A data frame containing either a `Length` column or columns `Start_X`,
  `Start_Y`, `Start_Z`, `End_X`, `End_Y`, `End_Z`.

## Value

A numeric vector of scar lengths, one per row.

## Details

The `Length` column (if present) takes precedence over
coordinate-derived lengths. This allows pre-computed lengths from
external software to be used directly without recalculation.

## See also

[`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)

## Examples

``` r
if (FALSE) { # \dontrun{
lens <- get_scar_length(df_one)
longest <- df_one[which.max(lens), ]
} # }
```
