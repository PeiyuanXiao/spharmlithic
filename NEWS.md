# spharmlithic 0.2.0

## New features

* Scar pattern measures that can be compared between cores with different scar counts (Xiao et al., submitted). They are written in R and do not need the Python back-end.
  * `spi_ei_batch()` computes raw and rarefied SPI, ER and IR for every core, with their random intervals; `random_interval()` and `plot_random_interval()` give and plot the random interval at a scar count.
  * `read_outline_iges()`, `fit_intersection_plane()`, `cpi_pra_batch()` and `cpi_pra_per_scar()` compute the Centripetality Index (CPI) and Plane-Relative Angle (PRA) of peripherally exploited cores.
  * `spharm_spatiodirection()` computes spatiodirectional SP-SPHARM features from pairs of different scars and from scar positions.
  * `compare_indices()` and `compare_groups()` compare groups of cores on single indices and on their overall scar patterns.
  * `check_scars()` checks a scar table before analysis; `plot_core()` and `export_outline_html()` show cores with their outlines.

# spharmlithic 0.1.0

* First release, accompanying Xiao et al. (2027), *Journal of Archaeological Method and Theory*.
