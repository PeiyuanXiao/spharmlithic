# Elongation (ER) and Isotropy (IR) ratios from the orientation tensor

Constructs the 3 x 3 orientation (fabric) tensor of the scar directions
of one core, computes its eigenvalues \\\lambda_1 \ge \lambda_2 \ge
\lambda_3\\, and returns the elongation and isotropy ratios \$\$ER = 1 -
\lambda_2 / \lambda_1, \quad IR = \lambda_3 / \lambda_1.\$\$ As for
[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md),
each is returned raw, rarefied, and with its random interval at the
core's scar count.

## Usage

``` r
compute_ei(
  dx,
  dy,
  dz,
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

- raw:

  Logical. Return the raw ratios (`ER_raw`, `IR_raw`) and the
  eigenvalues? Default `TRUE`.

- rarefied:

  Logical. Return the rarefied ratios (`ER_rarefied`, `IR_rarefied`)?
  Default `TRUE`.

- random_interval:

  Logical. Return the random intervals (`ER_lower`, `ER_upper`,
  `IR_lower`, `IR_upper`) and the position of the raw ratios relative to
  them (`ER_vs_random`, `IR_vs_random`)? Default `TRUE`.

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

A one-row tibble with `n_scar`, then for ER and for IR the columns
selected by `raw`, `rarefied` and `random_interval` (`*_raw`,
`*_rarefied`, `*_lower`, `*_upper`, `*_vs_random`), and, with
`raw = TRUE`, the eigenvalues `lambda1`, `lambda2`, `lambda3` of the
orientation tensor of all scars. All values are `NA` for cores with
fewer than three scars.

## Details

ER measures how far the scars are aligned along one axis, whether they
run the same way or opposite ways; IR measures how far the scar
directions spread into three dimensions (low IR: directions concentrated
in a plane or along an axis). ER above the upper bound and IR below the
lower bound of the random interval indicate directions more concentrated
than random.

Raw ER and IR shift with the scar count; compare cores and groups with
the rarefied values (see
[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md)).
Negative eigenvalues caused by numerical error are clamped to zero, and
both ratios are `NA` when \\\lambda_1 \approx 0\\.

**Changes in 0.2.0.** The columns `E` and `I` are renamed `ER_raw` and
`IR_raw`; the input vectors are normalised inside the function.

## References

Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., Jatmiko, &
Sutikna, T. (2024). A new method for quantifying flake scar organisation
on cores using orientation statistics. *Journal of Archaeological
Science*, **167**, 105998.

McPherron, S. P. (2018). Additional statistical and graphical methods
for analyzing site formation processes using artifact orientations.
*PLOS ONE*, **13**(1), e0190195.

## See also

[`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md),
[`spi_ei_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/spi_ei_batch.md),
[`random_interval()`](https://peiyuanxiao.github.io/spharmlithic/reference/random_interval.md)

## Examples

``` r
set.seed(1)
d <- matrix(rnorm(45), 15, 3)
compute_ei(d[, 1], d[, 2], d[, 3])
#> # A tibble: 1 × 14
#>   n_scar ER_raw ER_rarefied ER_lower ER_upper ER_vs_random IR_raw IR_rarefied
#>    <int>  <dbl>       <dbl>    <dbl>    <dbl> <chr>         <dbl>       <dbl>
#> 1     15  0.507       0.622   0.0668    0.576 within        0.221      0.0582
#> # ℹ 6 more variables: IR_lower <dbl>, IR_upper <dbl>, IR_vs_random <chr>,
#> #   lambda1 <dbl>, lambda2 <dbl>, lambda3 <dbl>
```
