

<!-- README.md is generated from README.qmd. Please edit that file -->

<p align="center">

<br /> <samp>R PACKAGE</samp> <br />

<h1 align="center">

<b>spharmlithic</b>

</h1>

<p align="center">

<b>Spherical Harmonic Analysis of Lithic Morphology and Flaking Scar
Patterns</b>

</p>

<hr />

</p>

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R ≥
4.1](https://img.shields.io/badge/R-%E2%89%A5%204.1-276DC3?logo=r)](https://cran.r-project.org/)
[![R-CMD-check](https://github.com/PeiyuanXiao/spharmlithic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PeiyuanXiao/spharmlithic/actions/workflows/R-CMD-check.yaml)

**spharmlithic** is an R package for the quantitative analysis of 3D
stone artefacts. It supports two complementary lines of inquiry:
**artefact morphology** (overall 3D shape from mesh surfaces) and
**flaking scar patterns** (orientation and organisation of flaking
removals). The package brings together two coordinate-alignment
pipelines, a suite of descriptive statistics (e.g. SPI), interactive 3D
visualisation via Plotly, and an optional Python back-end for spherical
harmonic decomposition — all in a single, script-based workflow that
goes from raw data to publication-ready outputs.

The package is designed for **lithic analysts** who work with 3D-scanned
cores, shaped tools or other artifacts and want reproducible results
without requiring prior experience with package development or Python.

------------------------------------------------------------------------

### 📦 Installation

#### Step 1 — Install the R package

``` r
# install.packages("remotes")
remotes::install_github("PeiyuanXiao/spharmlithic")
```

#### Step 2 — Install the Python back-end

Spherical harmonic analysis is handled by a bundled Python module that
runs behind the scenes via **reticulate**. If you only need alignment
and descriptive statistics you can skip this step entirely.

``` r
library(spharmlithic)

# First-time setup — creates a dedicated conda environment "r-spharmlithic"
install_spharmlithic_python()

# If you also need mesh-based morphological analysis (Track A), install the mesh extras:
install_spharmlithic_python(mesh = TRUE)
```

On subsequent R sessions, activate the environment with:

``` r
use_spharmlithic_python("r-spharmlithic")
```

> **Note:** A working
> [conda](https://docs.conda.io/en/latest/miniconda.html) installation
> is required. If conda is not on your system PATH, point R to it before
> calling the install function — for example by adding
> `Sys.setenv(RETICULATE_CONDA = "path/to/conda")` to your
> `~/.Rprofile`.

------------------------------------------------------------------------

### 🚀 Quick Start

``` r
library(spharmlithic)
library(readxl)

# Activate the Python back-end (skip if you only need alignment / compute SPI & EI)
use_spharmlithic_python("r-spharmlithic")

# ── Track B: scar pattern analysis ──────────────────────────

# 1. Load the bundled example scar data
scar_path <- system.file("extdata", "example_scars.xlsx",
                          package = "spharmlithic")
raw_data  <- read_excel(scar_path)

# 2. Align all specimens (SVD pipeline)
aligned <- align_scar_batch(raw_data)

# 3. Descriptive statistics (per specimen)
library(dplyr)
aligned %>%
  group_by(ID) %>%
  summarise(
    SPI = compute_SPI(d_x, d_y, d_z),
    SPI_angle = compute_spi_angle(d_x, d_y, d_z),
    compute_EI(d_x, d_y, d_z)
  )

# 4. Spherical harmonic analysis on aligned directions
sh_scar <- spharm_from_directions(aligned, lmax = 20, bandwidth = 0.35)

# 5. Export an interactive HTML report for alignment
export_alignment_html_svd(aligned, out_path = "alignment_report.html")

# ── Track A: morphological analysis ─────────────────────────

# 6. Point to the bundled example meshes
stl_dir <- system.file("extdata", "meshes", package = "spharmlithic")

# 7. Spherical harmonic analysis on 3D shape
sh_morph <- spharm_from_meshes(stl_dir, lmax = 20)

# ── Interactive SPHARM viewer ───────────────────────────────

# 8. Export an interactive 3D viewer for spherical harmonic reconstructions
#    (supports morph only, scar only, or both)
export_spharm_html(
  morph    = sh_morph,
  scar     = sh_scar,
  out_path = "spharm_viewer.html"
)
```

------------------------------------------------------------------------

### 🔧 Main Functions

#### Alignment

| Function | Description |
|:---|:---|
| `align_scar()` / `align_scar_batch()` | SVD three-step alignment of scar orientation vectors |
| `align_morph()` / `align_morph_batch()` | Two-step alignment following Lin et al. (2024) |

#### Descriptive Statistics

| Function | Description |
|:---|:---|
| `compute_SPI()` | Scar Pattern Index (unweighted by default; pass `lengths` for Clarkson’s original version) |
| `compute_spi_angle()` | Scar Pattern Angle (angular conversion: `θ = arccos(SPI)`) |
| `compute_EI()` | Elongation and Isotropy ratio from the orientation tensor |
| `get_scar_length()` | Individual scar lengths from coordinate data |
| `get_rot_matrix()` | Rotation matrix between two 3D vectors (Rodrigues formula) |

#### Spherical Harmonic Analysis (Python back-end)

| Function | Description |
|:---|:---|
| `spharm_from_directions()` | Compute SH coefficients from scar orientation vectors (core scar patterns) |
| `spharm_from_meshes()` | Compute SH coefficients from 3D mesh surfaces (core morphology; requires mesh extras) |
| `spharm_to_dataframe()` | Convert results to a wide-format data frame for downstream analysis |
| `spharm_reconstruct()` | Reconstruct a density grid from SH coefficients (low-level utility) |
| `export_spharm_html()` | Export an interactive 3D viewer for SH reconstructions |

#### Visualisation & Export

| Function | Description |
|:---|:---|
| `add_scars_3d()` / `add_arrow_3d()` | Add scar vectors or arrows to a Plotly scene |
| `add_plane_3d()` / `add_tilted_plane_3d()` | Add reference planes to a 3D plot |
| `build_panel_scar()` / `build_panel_morph()` | Assemble multi-panel Plotly figures |
| `export_alignment_html_svd()` | Export SVD alignment as a self-contained HTML file |
| `export_alignment_html_lin2024()` | Export Lin 2024 alignment as a self-contained HTML |

------------------------------------------------------------------------

### ↔️ Two Alignment Pipelines

**spharmlithic** offers two ways to bring scar orientation data into a
common coordinate frame:

|  | SVD alignment | Lin 2024 alignment |
|:---|:---|:---|
| **Functions** | `align_scar()` / `align_scar_batch()` | `align_morph()` / `align_morph_batch()` |
| **Approach** | Singular value decomposition of the full direction-vector cloud | Two-step rotation along a morphological long-axis (Lin et al. 2024) |
| **Output** | Identical format — all downstream functions work with either | ← same |

------------------------------------------------------------------------

### 🌐 Spherical Harmonic Analysis

Spherical harmonic (SPHARM) decomposition expands a function defined on
the sphere into a series of orthogonal basis functions, analogous to a
Fourier transform on a circle. In lithic analysis this provides a
compact, quantitative representation that supports frequency-based
comparison across specimens or assemblages.

**spharmlithic** applies this framework to two complementary aspects of
stone artefact variability:

- **Track A — 3D morphology.** `spharm_from_meshes()` reads 3D mesh
  files (.stl/.ply), extracts surface normals, and maps them onto the
  unit sphere to characterise overall shape via spherical harmonic
  coefficients. This track requires the mesh extras
  (`install_spharmlithic_python(mesh = TRUE)`).

- **Track B — Scar patterns.** `spharm_from_directions()` takes
  pre-computed scar orientation vectors, applies von Mises–Fisher kernel
  density estimation on the sphere, and decomposes the resulting density
  into spherical harmonic coefficients. This captures how flaking
  removals are organised on core reduction.

Both functions return spherical harmonic coefficient arrays that can be
flattened into a wide-format data frame with `spharm_to_dataframe()` for
downstream multivariate analysis (PCA, clustering, distance matrices)
directly in R.

#### Interactive 3D viewer

`export_spharm_html()` generates a self-contained HTML file that
reconstructs the spherical harmonic surfaces in real time using
Three.js. The viewer supports:

- Dual viewports (morphology and scar direction side-by-side, with
  synchronised rotation)
- Degree-by-degree slider and animation (l = 1 → lmax)
- Six material presets and a radial deviation colormap
- Type-mean reconstruction (averaging coefficients within a Typology
  group)
- PNG screenshot and OBJ mesh export

#### Low-level reconstruction

For users who need the reconstructed density values as numerical data
(e.g. for computing point-wise differences between specimens, extracting
density peaks, or building custom visualisations),
`spharm_reconstruct()` performs the inverse spherical harmonic transform
and returns the density matrix, grid coordinates, and unit-sphere
Cartesian points.

#### Harmonic degree (`lmax`)

The default maximum degree is `lmax = 20`, which captures the dominant
features for most artefacts. For Track A (morphology), users who need
finer geometric detail can increase the degree — for example
`spharm_from_meshes(stl_dir, lmax = 50)` — but note that computation
time and output size grow substantially with higher degrees (the number
of coefficients scales as `2 × (lmax+1)²`).

Track B (scar patterns) does not benefit from higher degrees because the
underlying von Mises–Fisher KDE already smooths the directional density;
increasing `lmax` beyond 20 adds fitting noise rather than meaningful
detail, and is therefore not recommended.

------------------------------------------------------------------------

### 📊 Input Data Format

#### Track A: 3D mesh files

Place `.stl` or `.ply` mesh files in a single directory. Each file is
treated as one specimen; the filename (without extension) becomes the
specimen ID. No particular naming convention is required, but avoid
spaces and special characters in filenames.

**Mesh preparation:** Although the package includes automatic
pre-decimation for high-resolution models (\>3M faces), it is
recommended to simplify your meshes to a manageable face count
(e.g. 50,000–200,000 faces) in advance using dedicated 3D processing
software (e.g. MeshLab, CloudCompare). This significantly reduces
computation time without affecting the spherical harmonic results at
typical analysis degrees (lmax ≤ 20–50). In addition, input meshes must
be **watertight** (closed, manifold, genus-0), repair any holes,
non-manifold edges, or self-intersections before analysis.

``` r
stl_dir <- "path/to/my/meshes"
list.files(stl_dir, pattern = "\\.(stl|ply)$")
# [1] "Specimen_A.stl"  "Specimen_B.stl"
```

#### Track B: Scar orientation data (.xlsx / .csv)

The scar analysis pipeline expects a data frame with **one row per
scar**, where all scars from all specimens are stacked together. The
`ID` column identifies which specimen each scar belongs to.

**Required columns:**

| Column                          | Description                         |
|:--------------------------------|:------------------------------------|
| `ID`                            | Specimen identifier                 |
| `Start_X`, `Start_Y`, `Start_Z` | XYZ coordinates of scar start-point |
| `End_X`, `End_Y`, `End_Z`       | XYZ coordinates of scar end-point   |

**Additional columns required by `align_morph()` / `align_morph_batch()`
(Lin 2024 pipeline):**

| Column | Description |
|:---|:---|
| `Norm_X`, `Norm_Y`, `Norm_Z` | Morphological best-fitted plane normal vector |
| `Pos_X`, `Pos_Y`, `Pos_Z` | Centroid of the morphological best-fitted plane (used for visualisation in `build_panel_morph()`) |

**Optional columns:**

| Column | Description |
|:---|:---|
| `Scar_ID` | Individual scar identifier within a specimen |
| `Length` | Pre-computed scar length (if absent, computed from start/end coordinates) |
| `Typology` | Specimen type label (e.g. “Levallois”, “Discoid”, “Laminar”); used for grouping in the interactive SPHARM viewer |

**Minimal example (3 scars from 2 specimens):**

| ID   | Start_X | Start_Y | Start_Z | End_X | End_Y | End_Z |  Typology |
|:-----|--------:|--------:|--------:|------:|------:|------:|----------:|
| S001 |    10.2 |     5.1 |     3.0 |  12.4 |   6.3 |   2.8 | Levallois |
| S001 |    11.0 |     4.8 |     3.1 |  13.5 |   5.9 |   2.5 | Levallois |
| S002 |    20.1 |     8.3 |     1.5 |  22.0 |   9.1 |   1.2 |   Discoid |

The bundled example data (`inst/extdata/example_scars.xlsx`) follows
this format and can be used as a template:

``` r
scar_path <- system.file("extdata", "example_scars.xlsx", package = "spharmlithic")
example   <- readxl::read_excel(scar_path)
names(example)
```

#### Output columns after alignment

Both alignment pipelines (`align_scar_batch()` and
`align_morph_batch()`) add nine columns to the input data frame:

| Column | Description |
|:---|:---|
| `s_x`, `s_y`, `s_z` | Aligned scar start-point coordinates |
| `e_x`, `e_y`, `e_z` | Aligned scar end-point coordinates |
| `d_x`, `d_y`, `d_z` | Aligned unit direction vectors — the primary input for `compute_SPI()`, `compute_EI()`, and `spharm_from_directions()` |

------------------------------------------------------------------------

### 📖 References

Clarkson, C., Vinicius, L., & Lahr, M. M. (2006). Quantifying flake scar
patterning on cores using 3D recording techniques. *Journal of
Archaeological Science*, 33(1), 132–142.

Bretzke, K., & Conard, N. J. (2012). Evaluating morphological
variability in lithic assemblages using 3D models of stone artifacts.
*Journal of Archaeological Science*, 39(12), 3741–3749.

McPherron, S. P. (2018). Additional statistical and graphical methods
for analyzing site formation processes using artifact orientations.
*PLoS ONE*, 13(1), e0190195.

Wieczorek, M. A., & Meschede, M. (2018). SHTools: Tools for working with
spherical harmonics. *Geochemistry, Geophysics, Geosystems*, 19(8),
2574–2592.

Lin, S. C., Clarkson, C., Julianto, I. M. A., Ferdianto, A., & Sutikna,
T. (2024). A new method for quantifying flake scar organisation on cores
using orientation statistics. *Journal of Archaeological Science*, 167,
105998.

Ye, Z., Pei, S. W., Ma. D. D., Li, H., Marwick, B. (2026). Spherical
harmonic analysis of faceted spheroids identifies shaping strategies and
standardisation at Qianshangying (North China). *Journal of
Archaeological Science*, 190, 106551.

------------------------------------------------------------------------

### 📝 Citation

If you use **spharmlithic** in your research, please cite it as:

> Xiao, P. Y., Marwick, B. (2026). spharmlithic: Spherical Harmonic
> Analysis of Lithic Morphology and Flaking Scar Patterns. R package.
> https://github.com/PeiyuanXiao/spharmlithic

------------------------------------------------------------------------

### ⚖️ License

This project is licensed under the **MIT License** — see the
[LICENSE](LICENSE) file for details.
