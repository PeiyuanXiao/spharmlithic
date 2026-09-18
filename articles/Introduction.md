# Introduction to spharmlithic

## Overview

**spharmlithic** is an R package aimed at supporting two complementary
lines of quantitative analysis of 3D shapes, with a focus on stone
artefacts:

- **Morphological analysis.** Spherical harmonic decomposition of 3D
  mesh surfaces (.stl), capturing overall artefact shape. This analysis
  can be adapted to any object that can be represented with a 3D mesh
  surface.
- **Scar pattern analysis.** Alignment, descriptive statistics, and
  spherical harmonic decomposition of flake scar orientation vectors,
  capturing how flaking removals are organised during core reduction.

This vignette is a getting-started guide: it walks through both
workflows using the small example data bundled with the package, taking
you from raw artefact data files to spherical harmonic coefficients and
an exported results table.

Once you have that table, our companion vignette
[`vignette("Statistics", package = "spharmlithic")`](https://peiyuanxiao.github.io/spharmlithic/articles/Statistics.md)
shows how to explore, analyse and visualize it with widely used
statistical methods. That companion vignette guides you through a
exploration of the data with PCA, group comparisons with PERMANOVA, and
tests of standardisation. These steps form the core of a typical
research project into 3D shapes.

------------------------------------------------------------------------

## Setup: Install the Package and the Python Environment

### Installing the R package

``` r

# install.packages("remotes")
remotes::install_github("PeiyuanXiao/spharmlithic")
```

### Installing the Python back-end

Spherical harmonic analysis requires a Python environment with
`pyshtools` and related packages. This step creates (or recreates) a
conda environment containing the Python packages required for spherical
harmonic analysis. This is a one-time setup step. Similar to
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
you only need to do this once per computer, not every time you work on
your project. If you only need alignment and descriptive statistics you
can skip this step entirely.

``` r

library(spharmlithic)

# Core install (scar pattern analysis only)
install_spharmlithic_python()

# Full install (also enables morphological analysis from meshes)
install_spharmlithic_python(mesh = TRUE)
```

### Optional: Docker

Another option for setting up the Python environment is to work within a
Docker container. This is convenient if you have an existing Python
install that you don’t want to modify, or are unable to.

Docker is also the supported way to run the mesh pipeline on **macOS**.
The direction-vector pipeline
([`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md))
works in a native macOS installation, but
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)
does not: `open3d` ships its own copy of the OpenMP runtime (`libomp`),
CRAN’s R for macOS ships another (used by R packages such as
`data.table`), and when both are loaded into the same R session the
runtime aborts it with `OMP: Error #15`. The pre-built Docker image runs
Linux, where this conflict does not arise:

``` bash
docker pull peiyuanxiao/spharmlithic
docker run -d -p 8787:8787 \
  -v /path/to/your/data:/home/rstudio/data \
  peiyuanxiao/spharmlithic
```

Open `http://localhost:8787` in your browser (user: `rstudio`, password:
`rstudio`). The container includes R, RStudio Server, and the full
Python environment pre-configured. Example data and a quick-start script
are in `~/examples/`. From this point on the instructions are the same
if you are working in the Docker container or directly on the desktop.

## Working with the Package: Activating the Python environment

With the Python environment installed, at the start of each R session,
all that’s required is to point `reticulate` to the conda environment:

``` r

library(spharmlithic)
library(readxl)

# point reticulate to the conda environment
use_spharmlithic_python("r-spharmlithic")
```

------------------------------------------------------------------------

## Inputting Data: Formats

### Scar orientation data (Excel / CSV)

The scar analysis pipeline expects a data frame with **one row per
scar**, where all scars from all specimens are stacked together.

**Required columns:**

| Column                          | Description                         |
|:--------------------------------|:------------------------------------|
| `ID`                            | Specimen identifier                 |
| `Start_X`, `Start_Y`, `Start_Z` | XYZ coordinates of scar start-point |
| `End_X`, `End_Y`, `End_Z`       | XYZ coordinates of scar end-point   |

**Additional columns required by
[`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)
(Lin 2024 pipeline):**

| Column | Description |
|:---|:---|
| `Norm_X`, `Norm_Y`, `Norm_Z` | Morphological best-fitted plane normal vector |
| `Pos_X`, `Pos_Y`, `Pos_Z` | Centroid of the morphological best-fitted plane |

**Optional columns:**

| Column | Description |
|:---|:---|
| `Scar_ID` | Individual scar identifier within a specimen |
| `Length` | Pre-computed scar length (if absent, computed from start/end coordinates) |
| `Typology` | Specimen type label (e.g. “Levallois”, “Discoid”); used for grouping in the interactive SPHARM viewer |

**Minimal example (3 scars from 2 specimens):**

| ID   | Start_X | Start_Y | Start_Z | End_X | End_Y | End_Z |  Typology |
|:-----|--------:|--------:|--------:|------:|------:|------:|----------:|
| S001 |    10.2 |     5.1 |     3.0 |  12.4 |   6.3 |   2.8 | Levallois |
| S001 |    11.0 |     4.8 |     3.1 |  13.5 |   5.9 |   2.5 | Levallois |
| S002 |    20.1 |     8.3 |     1.5 |  22.0 |   9.1 |   1.2 |   Discoid |

The bundled example data follows this format:

``` r

scar_path <- system.file("extdata", "example_scars.xlsx",
                          package = "spharmlithic")
example   <- readxl::read_excel(scar_path)
names(example)
#>  [1] "ID"       "Pos_X"    "Pos_Y"    "Pos_Z"    "Norm_X"   "Norm_Y"  
#>  [7] "Norm_Z"   "Scar_ID"  "Start_X"  "Start_Y"  "Start_Z"  "End_X"   
#> [13] "End_Y"    "End_Z"    "Typology"
```

### 3D mesh files

Place `.stl` mesh files in a single directory. Each file is treated as
one specimen; the filename (without extension) becomes the specimen ID.

**Mesh preparation recommendations:**

- Simplify meshes to 50,000–200,000 faces in advance using dedicated 3D
  software (e.g. MeshLab, CloudCompare). The package includes automatic
  pre-decimation for very high-resolution models (\>3M faces), but
  pre-simplification significantly reduces computation time.
- Meshes must be **watertight** (closed, manifold, genus-0). Repair any
  holes, non-manifold edges, or self-intersections before analysis.

``` r

stl_dir <- system.file("extdata", "meshes", package = "spharmlithic")
list.files(stl_dir, pattern = "\\.stl$", ignore.case = TRUE)
#> [1] "ClarksonEXP01_Levallois.stl" "ClarksonEXP02_Discoid.stl"
```

### Output columns after alignment

Both alignment pipelines
([`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
and
[`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md))
add nine columns to the input data frame:

| Column | Description |
|:---|:---|
| `s_x`, `s_y`, `s_z` | Aligned scar start-point coordinates |
| `e_x`, `e_y`, `e_z` | Aligned scar end-point coordinates |
| `d_x`, `d_y`, `d_z` | Aligned unit direction vectors — the primary input for [`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md), [`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md), and [`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md) |

------------------------------------------------------------------------

## Scar Pattern Analysis

This is the first major line of quantitative analysis that this package
supports. We draw on the pioneering work of Clarkson et al. (2006) and
Lin et al. (2024) to provide functions to streamline the analyses they
introduced.

### Load the example data

``` r

scar_path <- system.file("extdata", "example_scars.xlsx",
                          package = "spharmlithic")
raw_data  <- read_excel(scar_path)
head(raw_data)
#> # A tibble: 6 × 15
#>   ID      Pos_X Pos_Y Pos_Z Norm_X Norm_Y Norm_Z Scar_ID Start_X Start_Y Start_Z
#>   <chr>   <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <chr>     <dbl>   <dbl>   <dbl>
#> 1 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S1       -10.7    -21.3    20.7
#> 2 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S2         6.37   -34.2    26.7
#> 3 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S3        14.9    -38.7    21.4
#> 4 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S4        38.7    -33.3    23.1
#> 5 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S5        52.2    -22.7    21.8
#> 6 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S6        50.8     13.5    11.6
#> # ℹ 4 more variables: End_X <dbl>, End_Y <dbl>, End_Z <dbl>, Typology <chr>
```

### Alignment

Before comparing scar patterns across specimens, the orientation vectors
need to be rotated into a common coordinate frame. **spharmlithic**
provides two alignment pipelines.

#### SVD alignment (`align_scar_batch`)

The SVD method defines three principal axes from the full set of scar
direction vectors via singular value decomposition. It is
general-purpose and does not require a predefined morphological axis.

``` r

aligned_svd <- align_scar_batch(raw_data)
head(aligned_svd)
#> # A tibble: 6 × 27
#>   ID      Pos_X Pos_Y Pos_Z Norm_X Norm_Y Norm_Z Scar_ID Start_X Start_Y Start_Z
#>   <chr>   <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <chr>     <dbl>   <dbl>   <dbl>
#> 1 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S1       -10.7    -21.3    20.7
#> 2 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S2         6.37   -34.2    26.7
#> 3 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S3        14.9    -38.7    21.4
#> 4 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S4        38.7    -33.3    23.1
#> 5 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S5        52.2    -22.7    21.8
#> 6 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S6        50.8     13.5    11.6
#> # ℹ 16 more variables: End_X <dbl>, End_Y <dbl>, End_Z <dbl>, Typology <chr>,
#> #   Direct_X <dbl>, Direct_Y <dbl>, Direct_Z <dbl>, s_x <dbl>, s_y <dbl>,
#> #   s_z <dbl>, e_x <dbl>, e_y <dbl>, e_z <dbl>, d_x <dbl>, d_y <dbl>, d_z <dbl>
```

#### Lin et al. 2024 alignment (`align_morph_batch`)

The method of Lin et al. (2024) first orients the specimen along a
morphological long-axis and then rotates around that axis to align the
dominant flaking direction. It requires `Norm_X/Y/Z` and `Pos_X/Y/Z`
columns in the input data.

``` r

aligned_lin <- align_morph_batch(raw_data)
head(aligned_lin)
#> # A tibble: 6 × 27
#>   ID      Pos_X Pos_Y Pos_Z Norm_X Norm_Y Norm_Z Scar_ID Start_X Start_Y Start_Z
#>   <chr>   <dbl> <dbl> <dbl>  <dbl>  <dbl>  <dbl> <chr>     <dbl>   <dbl>   <dbl>
#> 1 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S1       -10.7    -21.3    20.7
#> 2 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S2         6.37   -34.2    26.7
#> 3 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S3        14.9    -38.7    21.4
#> 4 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S4        38.7    -33.3    23.1
#> 5 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S5        52.2    -22.7    21.8
#> 6 Clarks…  22.0 -8.18  18.0   0.08 -0.286 -0.955 S6        50.8     13.5    11.6
#> # ℹ 16 more variables: End_X <dbl>, End_Y <dbl>, End_Z <dbl>, Typology <chr>,
#> #   Direct_X <dbl>, Direct_Y <dbl>, Direct_Z <dbl>, s_x <dbl>, s_y <dbl>,
#> #   s_z <dbl>, e_x <dbl>, e_y <dbl>, e_z <dbl>, d_x <dbl>, d_y <dbl>, d_z <dbl>
```

Both pipelines produce output in the same format, so all downstream
functions work identically regardless of which alignment you choose. We
use the SVD-aligned data for the rest of this vignette.

#### Alignment visualisation

To export an interactive HTML report showing 3D visualisations of the
alignment steps:

``` r

export_alignment_html_svd(raw_data, 
                          out_path = "alignment_svd.html")

export_alignment_html_lin2024(raw_data, 
                              out_path = "alignment_lin2024.html")
```

### Descriptive statistics

#### Scar Pattern Index (SPI)

The SPI quantifies the overall directionality of flaking removals. The
default behaviour is unweighted, following Bretzke & Conard (2012),
where all scars contribute equally regardless of length. Values close to
1 indicate strong preferred orientation; values close to 0 indicate
isotropic / random patterning.

``` r

spi <- compute_spi(aligned_svd$d_x, 
                   aligned_svd$d_y, 
                   aligned_svd$d_z)
spi
#> [1] 0.2373782
```

For the original, length-weighted SPI (Clarkson et al. 2006), pass scar
lengths:

``` r

lens <- get_scar_length(aligned_svd)
spi_w <- compute_spi(aligned_svd$d_x, 
                     aligned_svd$d_y, 
                     aligned_svd$d_z,
                     lengths = lens)
spi_w
#> [1] 0.2665681
```

An angular variant is also available, this converts the Scar Pattern
Index into the expected pairwise angle between two randomly-selected
scars, following Clarkson et al.’s (2006) interpretation:

``` r

spi_angle <- compute_spi_angle(aligned_svd$d_x,
                               aligned_svd$d_y,
                               aligned_svd$d_z)
spi_angle
#> [1] 76.26815
```

#### Elongation and Isotropy

[`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md)
returns the elongation (E) and isotropy (I) ratios derived from the
orientation tensor eigenvalues (Lin et al. 2024), whose implementation
adapts the fabric analysis of McPherron (2018):

``` r

ei <- compute_ei(aligned_svd$d_x, 
                 aligned_svd$d_y, 
                 aligned_svd$d_z)
ei
#>           E         I   lambda1  lambda2   lambda3
#> 1 0.3699931 0.3725705 0.4993565 0.314598 0.1860455
```

#### Per-specimen statistics

Use
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
to compute the descriptive statistics for each specimen (row):

``` r

library(dplyr)
aligned_svd %>% 
  group_by(ID) %>%
  summarise(
    SPI       = compute_spi(d_x, d_y, d_z),
    SPI_angle = compute_spi_angle(d_x, d_y, d_z),
    compute_ei(d_x, d_y, d_z)
  )
```

### Spherical harmonic analysis

This step is a novel application of spherical harmonic analysis that
extends the work of Clarkson et al. (2006) and Lin et al. (2024). With
the Python back-end active, scar orientation vectors can be decomposed
into spherical harmonic coefficients. This function applies von
Mises–Fisher kernel density estimation on the sphere and then fits the
harmonic expansion:

``` r

sh_scar <- spharm_from_directions(aligned_svd, lmax = 20, bandwidth = 0.35)
sh_scar
```

The output here is a list, so to convert the results into a flat data
frame for downstream multivariate analysis (PCA, clustering, distance
matrices) we have a convenient fucntion:

``` r

df_scar <- spharm_to_dataframe(sh_scar)
head(df_scar)
```

------------------------------------------------------------------------

## Morphological Analysis

This is the second major line of quantitative analysis that this package
supports. Morphological analysis takes 3D mesh files (.stl) as input and
characterises overall artefact shape via spherical harmonics. This
requires the correctly configured Python environment and mesh extension:
`install_spharmlithic_python(mesh = TRUE)`. This
[`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)
function executes the morphological analysis pipline with the following
steps, per STL file:

- Read mesh via open3d. Files exceeding 3M faces are pre-decimated by
  streaming face sub-sampling to avoid memory pressure.  
- Quadric-error decimation to target_faces.  
- Laplacian smoothing (smooth_iterations passes).  
- Volume-centroid normalization to a unit sphere.  
- Area-weighted PCA alignment with energy-based sign convention.
- Cartesian-to-spherical conversion of the aligned vertices.
- Interpolation onto a regular grid.
- Spherical harmonic expansion via pyshtools (4pi normalization,
  zero-component normalized: all coefficients divided by c(0, 0)).

> **macOS users:**
> [`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)
> is not supported in a native macOS installation: the OpenMP runtime
> bundled with `open3d` clashes with the one that comes with R, and the
> R session aborts with `OMP: Error #15`. Please use the Docker image
> instead, see the Setup section above. The direction-vector pipeline is
> not affected.

``` r

stl_dir  <- system.file("extdata", "meshes", package = "spharmlithic")
sh_morph <- spharm_from_meshes(stl_dir, lmax = 20)
sh_morph
```

The output is a list with one element per successfully-processed STL, we
can convert to a flat data frame with a convenience function:

``` r

df_morph <- spharm_to_dataframe(sh_morph)
head(df_morph)
```

### Harmonic Degree (`lmax`)

The maximum spherical harmonic degree controls the level of detail in
the decomposition. The default is `lmax = 20`, which captures the
dominant features for most artefacts.

Morphological analysis can benefit from higher degrees if you need finer
geometric detail:

``` r

sh_fine <- spharm_from_meshes(stl_dir, lmax = 50)
```

Note that computation time and output size grow substantially with
higher degrees (the number of coefficients scales as `2 × (lmax + 1)²`).
When increasing `lmax` for morphology, you may also want to increase
`grid_size` (default 256) for adequate sampling.

Scar pattern analysis does not benefit from degrees above 20. The
underlying von Mises–Fisher KDE already smooths the directional density,
so increasing `lmax` beyond 20 adds noise fitting rather than meaningful
detail. If you do increase `lmax`, the DH grid size must satisfy
`dh_size >= 2 * (lmax + 1)`.

At this point you have all the outputs you need to use the statistical
methods we walk through in the next vignette. The remainder of this
vignette are optional visualisation steps.

------------------------------------------------------------------------

## Interactive SPHARM Viewer

It can be useful to visually inspect the 3D models and scar patterns. We
provide
[`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md)
to generate a self-contained HTML file with an interactive
Three.js-based 3D viewer for exploring spherical harmonic
reconstructions. The viewer synthesizes the harmonic expansion in real
time on the client side — only the coefficients are embedded, so file
sizes stay small.

You can pass results from either analysis, or both:

``` r

# Scar pattern only
export_spharm_html(
  scar     = sh_scar,
  out_path = "viewer_scar.html"
)
```

``` r

# Both morphology and scar pattern side-by-side
export_spharm_html(
  morph    = sh_morph,
  scar     = sh_scar,
  out_path = "spharm_viewer.html"
)
```

### Specimen metadata and grouping

An optional `meta` argument accepts a data frame with specimen metadata.
If a `Typology` column is present, the viewer groups the specimen
dropdown by type and enables a Type Mean toggle that averages
coefficients within each group:

``` r

meta <- data.frame(
  ID       = c("ClarksonEXP01_Levallois", "ClarksonEXP02_Discoid"),
  Typology = c("Levallois", "Discoid")
)
export_spharm_html(
  scar     = sh_scar,
  meta     = meta,
  out_path = "spharm_viewer.html",
  title    = "Example specimens"
)
```

### Viewer features

- **Dual viewports** — morphology and scar direction side-by-side, with
  synchronised rotation.
- **Degree slider and animation** — l = 1 → lmax, for visualising how
  successive harmonic degrees build up the reconstruction.
- **Six material presets** — Ceramic, Clay, Glass, Brushed Metal, X-Ray,
  Flat — plus a radial deviation colormap.
- **Type Mean** — average coefficients across specimens sharing the same
  Typology.
- **PNG screenshot** and **OBJ mesh export** at the current degree.

### Low-level Reconstruction

For users who need the reconstructed density as numerical data rather
than an interactive visualisation,
[`spharm_reconstruct()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_reconstruct.md)
performs the inverse spherical harmonic transform and returns:

- `density` — numeric matrix (grid_size × 2 × grid_size)
- `lon`, `lat` — longitude and latitude vectors (radians)
- `xyz` — unit-sphere Cartesian coordinates

``` r

recon <- spharm_reconstruct(
  sh_scar[["ClarksonEXP01_Levallois"]]$coefficients,
  grid_size = 64
)
```

This is useful for computing point-wise differences between specimens,
extracting density peaks, or building custom visualisations:

``` r

library(plotly)
nlat <- length(recon$lat)
nlon <- length(recon$lon)
x_mat <- matrix(recon$xyz[, "x"], nlat, nlon)
y_mat <- matrix(recon$xyz[, "y"], nlat, nlon)
z_mat <- matrix(recon$xyz[, "z"], nlat, nlon)

plot_ly(x = x_mat, y = y_mat, z = z_mat,
        surfacecolor = recon$density,
        type = "surface", colorscale = "Hot",
        showscale = TRUE)
```

------------------------------------------------------------------------

## References

Bretzke, K., & Conard, N. J. (2012). Evaluating morphological
variability in lithic assemblages using 3D models of stone artifacts.
*Journal of Archaeological Science*, 39(12), 3741–3749.

Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
patterning on cores using 3D recording techniques. *Journal of
Archaeological Science*, 33(1), 132–142.

Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., Jatmiko, &
Sutikna, T. (2024). A new method for quantifying flake scar organisation
on cores using orientation statistics. *Journal of Archaeological
Science*, 167, 105998.

McPherron, S. P. (2018). Additional statistical and graphical methods
for analyzing site formation processes using artifact orientations.
*PLoS ONE*, 13(1), e0190195.

Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working with
spherical harmonics. *Geochemistry, Geophysics, Geosystems*, 19(8),
2574–2592.

Ye, Z., Pei, S. W., Ma, D. D., Li, H., & Marwick, B. (2026). Spherical
harmonic analysis of faceted spheroids identifies shaping strategies and
standardisation at Qianshangying (North China). *Journal of
Archaeological Science*, 190, 106551.
