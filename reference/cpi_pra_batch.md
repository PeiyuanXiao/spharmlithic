# CPI and PRA for every peripherally exploited core

Computes the Centripetality Index (CPI) and the Plane-Relative Angle
(PRA) of every scar of each core that has an outline, and summarises
them per core by their mean and Gini's mean difference (GMD).

## Usage

``` r
cpi_pra_batch(
  data,
  outlines,
  id_col = "ID",
  start_cols = c("Start_X", "Start_Y", "Start_Z"),
  end_cols = c("End_X", "End_Y", "End_Z"),
  min_length = 0,
  keep = NULL,
  quiet = FALSE
)
```

## Arguments

- data:

  A data frame with one row per scar and the columns named by `id_col`,
  `start_cols` and `end_cols`. Scars of all cores are stacked.

- outlines:

  A named list with one element per core, named by core ID: the outline
  points (a three-column matrix or data frame ordered along the loop, as
  returned by
  [`read_outline_iges()`](https://peiyuanxiao.github.io/spharmlithic/reference/read_outline_iges.md))
  or a plane already fitted with
  [`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md).
  Cores without an outline are left out. The outlines must be in the
  same coordinate frame as the scar columns named by `start_cols` and
  `end_cols`.

- id_col:

  Name of the core ID column. Default `"ID"`.

- start_cols, end_cols:

  Names of the start- and end-point coordinate columns (X, Y, Z).
  Defaults `Start_X/Y/Z` and `End_X/Y/Z`, the columns the scars are
  recorded in. The aligned columns of
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
  (`s_x`... and `e_x`...) give the same results: these indices do not
  change when a core is rotated or moved.

- min_length:

  Scars shorter than this (in the units of the coordinates) are left
  out. Default 0 (keep all). Use it to repeat an analysis without small
  scars, for example `min_length = 10` for 10 mm.

- keep:

  Names of per-core columns (for example `"Typology"`) to carry into the
  output. Each must be constant within a core.

- quiet:

  Logical. Suppress the summary message? Default `FALSE`.

## Value

A tibble with one row per core that has an outline: the ID column, the
`keep` columns, `n_scar` (scars used), `n_cpi` (scars with a defined
CPI), `CPI_mean`, `CPI_GMD`, `PRA_mean`, `PRA_GMD` (degrees), and three
checks of the outline: `flatness` (see
[`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md)),
`start_dist` (the typical distance of the scar starts from the outline,
as a fraction of the outline's size) and `frame_ok` (`start_dist` at
most 0.03).

## Details

**The angles.** The plane of intersection is fitted to the outline and
its centre is the centroid of the enclosed area (see
[`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md)).
CPI is the angle, in the plane, between a scar's removal direction and
the direction from the centre to its start: 180 degrees for a scar
running towards the centre, 90 for one running along the edge, 0 for one
running away from the centre. PRA is the angle between the removal
direction and the plane: 0 for a scar running along the plane, 90 for
one driven into the core perpendicular to it. PRA does not distinguish
the two sides of the plane. Neither angle wraps around, so ordinary
means apply.

**The core summaries.** The mean summarises the overall centripetality
or inclination of a core; GMD, the mean absolute difference between the
angles of all pairs of its scars, summarises how much they vary. Both
are averages over single scars or pairs of scars, so their expectations
do not change with the scar count, and they can be compared directly
between cores and groups (for example with
[`compare_indices()`](https://peiyuanxiao.github.io/spharmlithic/reference/compare_indices.md));
the median, interquartile range and standard deviation do not have this
property (Xiao et al., in prep.).

**When to use.** CPI and PRA describe cores flaked from their periphery
whose plane of intersection can be identified, such as discoidal and
Levallois cores. Supply outlines only for such cores.

**The frame check.** Scar starts are recorded on the periphery, so they
lie on or near the outline when scars and outline share a coordinate
frame. `start_dist` is the median distance of the starts from the
outline (over the 85% of starts closest to it) divided by the outline's
bounding-box diagonal. Cores whose `start_dist` exceeds 0.03 are flagged
with a warning: their scars and outline are probably in different frames
or units, and their angles are not meaningful. Use
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md)
to look.

Rows with missing coordinates or zero length, and scars shorter than
`min_length`, are left out. Cores with fewer than three usable scars get
`NA`. A scar whose start lies on the centre, or which runs exactly along
the normal, has no CPI; `n_cpi` counts the scars that have one.

## References

Gini, C. (1912). *Variabilità e mutabilità*. Bologna: C. Cuppini.

Yitzhaki, S. (2003). Gini's mean difference: a superior measure of
variability for non-normal distributions. *Metron*, **61**(2), 285–316.

## See also

[`read_outline_iges()`](https://peiyuanxiao.github.io/spharmlithic/reference/read_outline_iges.md),
[`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md),
[`cpi_pra_per_scar()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_per_scar.md),
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md),
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md)

## Examples

``` r
# A disc-shaped core: twelve scars from the edge towards the centre
t <- seq(0, 2 * pi, length.out = 13)[-13]
start <- cbind(30 * cos(t), 30 * sin(t), 0)
end <- start * 0.4
end[, 3] <- -5
scars <- data.frame(ID = "disc", Start_X = start[, 1], Start_Y = start[, 2],
                    Start_Z = start[, 3], End_X = end[, 1], End_Y = end[, 2],
                    End_Z = end[, 3])
u <- seq(0, 2 * pi, length.out = 200)
outlines <- list(disc = cbind(30 * cos(u), 30 * sin(u), 0))
cpi_pra_batch(scars, outlines)
#> Warning: Outline of core disc: collapsing to unique 'x' values
#> Warning: Outline of core disc: collapsing to unique 'x' values
#> Warning: Outline of core disc: collapsing to unique 'x' values
#> CPI and PRA computed for 1 core(s) with an outline.
#> # A tibble: 1 × 10
#>   ID    n_scar n_cpi CPI_mean CPI_GMD PRA_mean  PRA_GMD flatness start_dist
#>   <chr>  <int> <int>    <dbl>   <dbl>    <dbl>    <dbl>    <dbl>      <dbl>
#> 1 disc      12    12      180       0     15.5 2.80e-15        0  0.0000288
#> # ℹ 1 more variable: frame_ok <lgl>
```
