# Summarise power spectrum diagnostics by degree

Calculate per-degree coefficients of variation and cumulative mean power
from a power spectrum.

## Usage

``` r
degree_diagnostics(power_df, descriptor, max_degree = 20)
```

## Arguments

- power_df:

  A data frame containing columns named `"power_l1"`, `"power_l2"`, ...,
  `"power_lN"`.

- descriptor:

  Character string identifying the descriptor type (e.g., `"shape"`,
  `"amplitude"`, `"power"`).

- max_degree:

  Maximum harmonic degree to analyse.

## Value

A data frame with one row per degree containing:

- descriptor:

  Descriptor name supplied by the user.

- degree:

  Harmonic degree.

- mean_power:

  Mean power at that degree.

- cv_pct:

  Coefficient of variation (%).

- cumul_pct:

  Cumulative percentage of total mean power.

## Details

The coefficient of variation (CV) quantifies across-specimen variability
for each harmonic degree. The cumulative power indicates the proportion
of total mean power explained by successive degrees.

## Examples

``` r
x <- data.frame(
  power_l1 = c(0.4, 0.5, 0.6),
  power_l2 = c(0.3, 0.3, 0.2),
  power_l3 = c(0.3, 0.2, 0.2)
)

degree_diagnostics(x, "example", max_degree = 3)
#>   descriptor degree mean_power   cv_pct cumul_pct
#> 1    example      1  0.5000000 20.00000  50.00000
#> 2    example      2  0.2666667 21.65064  76.66667
#> 3    example      3  0.2333333 24.74358 100.00000
```
