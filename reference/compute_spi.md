# Scar Pattern Index (Clarkson et al. 2006)

Computes the Scar Pattern Index of one core: the length of the resultant
of its scar directions divided by the number of scars, from 0
(directions cancel out) to 1 (all scars run the same way). By default
the function returns three readings:

## Usage

``` r
compute_spi(
  dx,
  dy,
  dz,
  lengths = NULL,
  raw = TRUE,
  rarefied = TRUE,
  random_interval = TRUE,
  m = 3,
  max_comb = 2000,
  seed = 2226
)
```

## Arguments

- dx, dy, dz:

  Numeric vectors of equal length: the X, Y and Z components of the scar
  direction vectors of **one** core, for example `End_X - Start_X` etc.,
  or the `d_x`, `d_y`, `d_z` columns returned by
  [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md).
  The vectors are normalised to unit length inside the function; rows of
  zero length or with missing values are left out.

- lengths:

  Optional numeric vector of scar lengths, one per vector. If given,
  each unit direction is weighted by its length, reproducing the
  length-weighted SPI of Clarkson et al. (2006). The rarefied value and
  the random interval are defined for unweighted SPI only and are `NA`
  when `lengths` is given.

- raw:

  Logical. Return the raw SPI (`SPI_raw`)? Default `TRUE`.

- rarefied:

  Logical. Return the rarefied SPI (`SPI_rarefied`)? Default `TRUE`.

- random_interval:

  Logical. Return the random interval at the core's scar count
  (`SPI_lower`, `SPI_upper`) and the position of the raw SPI relative to
  it (`SPI_vs_random`)? Default `TRUE`.

- m:

  Number of scars per combination for the rarefied value. Default 3.

- max_comb:

  Largest number of combinations averaged. Above it, `max_comb` distinct
  combinations are drawn at random (a Monte Carlo approximation of the
  mean over all combinations). Default 2000.

- seed:

  Seed for that random draw; the caller's random stream is restored
  afterwards. `NULL` draws from the current stream. Default 2226.

## Value

A one-row tibble with `n_scar` (the number of scars used) and the
columns selected by `raw`, `rarefied` and `random_interval`: `SPI_raw`,
`SPI_rarefied`, `SPI_lower`, `SPI_upper` and `SPI_vs_random` (`"within"`
or `"outside"`). All values are `NA` for cores with fewer than three
scars.

## Details

- **raw** SPI, computed from all scars together;

- **rarefied** SPI, the mean SPI over every combination of `m` scars
  (Xiao et al., in prep.). Its expectation does not change with the scar
  count, so it is the value to use when comparing cores or groups of
  cores;

- the **random interval** at the core's scar count: the central 95% of
  the SPI of cores whose scar directions are uniform on the sphere, and
  whether the raw SPI lies `"within"` or `"outside"` it.

**Which value for which question.** Raw SPI shifts with the scar count:
the fewer the scars, the less their directions cancel out, even when
they are random. A single core is therefore judged by comparing its raw
SPI with the random interval at its own scar count; cores and groups of
cores are compared with the rarefied SPI, whose expectation does not
depend on the scar count (Xiao et al., in prep.).

**Reading the random interval.** SPI above the upper bound indicates
scar directions more concentrated on one side than random directions
produce at this scar count. SPI below the lower bound is typical of
opposed (bidirectional) flaking, whose two sets of removals cancel out.
A value `"within"` the interval cannot be told apart from random
directions at this scar count. See
[`random_interval()`](https://peiyuanxiao.github.io/spharmlithic/reference/random_interval.md).

**Changes in 0.2.0.** The function returns a one-row tibble instead of a
single number, and normalises the input vectors itself: raw displacement
vectors (end minus start) now give the unweighted SPI, and the
length-weighted SPI is requested with `lengths`.

## References

Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
patterning on cores using 3D recording techniques. *Journal of
Archaeological Science*, **33**(1), 132–142.

Bretzke, K., & Conard, N. J. (2012). Evaluating morphological
variability in lithic assemblages using 3D models of stone artifacts.
*Journal of Archaeological Science*, **39**(12), 3741–3749.

Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique
and alternative parameters. *Ecology*, **52**(4), 577–586.

## See also

[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md),
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`random_interval()`](https://peiyuanxiao.github.io/spharmlithic/reference/random_interval.md),
[`compute_spi_angle()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi_angle.md)

## Examples

``` r
# Twelve scars running broadly the same way
set.seed(1)
d <- matrix(rnorm(36, sd = 0.4), 12, 3)
d[, 1] <- d[, 1] + 1
compute_spi(d[, 1], d[, 2], d[, 3])
#> # A tibble: 1 × 6
#>   n_scar SPI_raw SPI_rarefied SPI_lower SPI_upper SPI_vs_random
#>    <int>   <dbl>        <dbl>     <dbl>     <dbl> <chr>        
#> 1     12   0.902        0.929    0.0790     0.500 outside      

# Raw SPI only
compute_spi(d[, 1], d[, 2], d[, 3], rarefied = FALSE, random_interval = FALSE)
#> # A tibble: 1 × 2
#>   n_scar SPI_raw
#>    <int>   <dbl>
#> 1     12   0.902

# Length-weighted SPI (Clarkson et al. 2006)
compute_spi(d[, 1], d[, 2], d[, 3], lengths = runif(12, 5, 30),
            rarefied = FALSE, random_interval = FALSE)
#> # A tibble: 1 × 2
#>   n_scar SPI_raw
#>    <int>   <dbl>
#> 1     12   0.918
```
