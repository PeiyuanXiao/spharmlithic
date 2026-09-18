# spharmlithic: SPHARM framework for stone artifact shape and scar pattern analysis

`spharmlithic` provides a complete toolkit for the geometric alignment,
spherical harmonic decomposition for artifact shape and flaking scar
vector

### Main components

#### Alignment pipelines

- [`align_scar_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_scar_batch.md)
  — Three-step SVD pipeline (rotate -\> translate -\> in-plane rotate).

- [`align_morph_batch()`](https://peiyuanxiao.github.io/spharmlithic/reference/align_morph_batch.md)
  — Two-step Lin 2024 pipeline (morphological normal rotation -\>
  longest-scar translation).

#### Classical orientation statistics

- [`compute_spi()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi.md)
  — Scar Pattern Index (Clarkson et al. 2006); supports both unit-vector
  and length-weighted variants.

- [`compute_spi_angle()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_spi_angle.md)
  — SPI converted to expected pairwise angle (Clarkson et al. 2006).

- [`compute_ei()`](https://peiyuanxiao.github.io/spharmlithic/reference/compute_ei.md)
  — Elongation (E) and Isotropy (I) from the orientation tensor (Lin et
  al. 2024).

#### Spherical harmonic analysis (Python-backed)

- [`install_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/install_spharmlithic_python.md)
  — One-time Python backend setup.

- [`use_spharmlithic_python()`](https://peiyuanxiao.github.io/spharmlithic/reference/use_spharmlithic_python.md)
  — Point at an existing Python env.

- [`spharm_from_directions()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_directions.md)
  — Direction vectors -\> vMF KDE -\> SPHARM.

- [`spharm_from_meshes()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_from_meshes.md)
  — STL meshes -\> spherical interpolation -\> SPHARM.

- [`spharm_to_dataframe()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_to_dataframe.md)
  — Flatten SH results for CSV export and downstream multivariate
  analysis (PCA, clustering).

- [`spharm_reconstruct()`](https://peiyuanxiao.github.io/spharmlithic/reference/spharm_reconstruct.md)
  — Low-level inverse transform: coefficients -\> density grid (for
  custom analysis or plotting).

#### Interactive viewers and HTML export

- [`export_spharm_html()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_spharm_html.md)
  — Export an interactive Three.js-based 3D viewer for spherical
  harmonic reconstructions (morphology, scar direction, or both
  side-by-side).

- [`export_alignment_html_svd()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_svd.md)
  — Four-panel SVD alignment page.

- [`export_alignment_html_lin2024()`](https://peiyuanxiao.github.io/spharmlithic/reference/export_alignment_html_lin2024.md)
  — Three-panel Lin 2024 alignment page.

#### Geometry utilities

- [`get_rot_matrix()`](https://peiyuanxiao.github.io/spharmlithic/reference/get_rot_matrix.md)
  — Rodrigues rotation matrix between two unit vectors.

- [`get_scar_length()`](https://peiyuanxiao.github.io/spharmlithic/reference/get_scar_length.md)
  — Scar lengths from a data frame.

## See also

Useful links:

- <https://github.com/PeiyuanXiao/spharmlithic>

- <https://peiyuanxiao.github.io/spharmlithic/>

- Report bugs at <https://github.com/PeiyuanXiao/spharmlithic/issues>

## Author

**Maintainer**: Peiyuan Xiao <pyxiao@uw.edu>
([ORCID](https://orcid.org/0009-0000-9733-5875))

Authors:

- Peiyuan Xiao <pyxiao@uw.edu>
  ([ORCID](https://orcid.org/0009-0000-9733-5875))

- Ben Marwick <bmarwick@uw.edu>
  ([ORCID](https://orcid.org/0000-0001-7879-4531))
