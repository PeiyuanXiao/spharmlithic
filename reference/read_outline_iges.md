# Read core outlines from IGES files

Reads the peripheral outline of each core from an IGES file and returns
it as a closed loop of points, ready for
[`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)
and
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md).

## Usage

``` r
read_outline_iges(files, ids = NULL, n_per_seg = 300, gap_frac = 0.01)
```

## Arguments

- files:

  Paths to IGES files (`.igs` / `.iges`), one per core.

- ids:

  Core IDs, one per file, matching the ID column of the scar table.
  Default: the file names without extension.

- n_per_seg:

  Number of points sampled along each B-spline segment. Default 300.

- gap_frac:

  Largest gap between segment ends that is bridged with a straight
  chord, as a fraction of the outline's bounding-box diagonal. Default
  0.01.

## Value

A named list (names = `ids`) of numeric matrices with columns `X`, `Y`,
`Z`, one point per row, ordered along the loop with the first point
repeated as the last. Attributes record how the loop was assembled:
`n_segments`, `worst_gap` (largest bridged gap), `n_bridged`,
`n_discarded` (separate curves that were not part of the loop) and
`file`. Files that cannot be read are left out with a warning.

## Details

The outline must be exported as rational B-spline curves (IGES entity
type 126), the form Geomagic Wrap and most CAD programs use for curves.
Segments may come in any order and orientation; they are chained into
connected loops. If an export carries stray curves away from the
periphery, the longest closed loop is taken as the outline. Gaps between
segment ends up to `gap_frac` of the outline's size are closed with a
straight chord; larger gaps are an error, because a missing stretch of
outline must not be invented. Files with other curve entities (arcs,
polylines, offset curves) are rejected: export the outline as B-splines.

The declared IGES unit is ignored. The outline must be in the same
coordinate frame as the scar coordinates of its core;
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md)
checks this.

## See also

[`fit_intersection_plane()`](https://peiyuanxiao.github.io/spharmlithic/reference/fit_intersection_plane.md),
[`cpi_pra_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/cpi_pra_batch.md),
[`plot_core()`](https://peiyuanxiao.github.io/spharmlithic/reference/plot_core.md)

## Examples

``` r
if (FALSE) { # \dontrun{
files    <- list.files("outlines", "\\.(igs|iges)$", full.names = TRUE)
outlines <- read_outline_iges(files)
names(outlines)

# File names that differ from the core IDs
outlines <- read_outline_iges(files, ids = c("EXP01", "EXP02"))
} # }
```
