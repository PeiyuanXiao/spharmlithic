# =============================================================================
# Dockerfile for spharmlithic
# Pre-built R + Python environment for spherical harmonic analysis of lithics.
# Solves the macOS open3d crash by running in a Linux container.
#
# Usage:
#   docker build -t peiyuanxiao/spharmlithic .
#   docker run -d -p 8787:8787 -v /path/to/your/data:/home/rstudio/data \
#     peiyuanxiao/spharmlithic
#   # Open http://localhost:8787  (user: rstudio, password: spharm)
#
# Push to Docker Hub:
#   docker push peiyuanxiao/spharmlithic
# =============================================================================

FROM rocker/rstudio:4.5.0

# --- 1. System dependencies --------------------------------------------------
#   libgl1 + libglu1-mesa: OpenGL for open3d (mesh pipeline)
#   xvfb: virtual framebuffer so open3d doesn't need a real display
#   wget, git: for Miniconda + GitHub installs
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    git \
    libgl1 \
    libglu1-mesa \
    xvfb \
    && rm -rf /var/lib/apt/lists/*

# --- 2. Install Miniconda ----------------------------------------------------
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
      -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# --- 3. Build Python environment (r-spharmlithic) ----------------------------
#   Accept conda TOS (required since conda 26.x), then create full env.
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda create -n r-spharmlithic -c conda-forge -y \
      python=3.10 pip \
      numpy scipy pandas trimesh && \
    /opt/conda/envs/r-spharmlithic/bin/pip install --no-cache-dir \
      pyshtools open3d && \
    conda clean -afy

# --- 4. Install spharmlithic R package from GitHub ----------------------------
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')" && \
    R -e "remotes::install_github('PeiyuanXiao/spharmlithic', dependencies = TRUE)"

# --- 5. Configure R to find conda + Python env -------------------------------
RUN echo 'Sys.setenv(RETICULATE_CONDA = "/opt/conda/bin/conda")' \
      >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'Sys.setenv(RETICULATE_PYTHON = "/opt/conda/envs/r-spharmlithic/bin/python")' \
      >> /usr/local/lib/R/etc/Rprofile.site

# --- 6. Copy example data into container -------------------------------------
RUN mkdir -p /home/rstudio/examples/meshes

# R is already installed with the package; copy extdata from the installed location
RUN cp /usr/local/lib/R/site-library/spharmlithic/extdata/example_scars.xlsx \
      /home/rstudio/examples/ && \
    cp /usr/local/lib/R/site-library/spharmlithic/extdata/meshes/*.stl \
      /home/rstudio/examples/meshes/

# --- 7. Create a startup script for quick-start ------------------------------
RUN cat > /home/rstudio/examples/quickstart.R << 'RSCRIPT'
# ============================================================
# spharmlithic Quick Start  (run line by line or source())
# ============================================================

library(spharmlithic)
use_spharmlithic_python("r-spharmlithic", check_mesh = TRUE)

# --- Scar Pattern Analysis ---
library(readxl)
raw_data <- read_excel("~/examples/example_scars.xlsx")
aligned  <- align_scar_batch(raw_data)

# Statistics
compute_SPI(aligned$d_x, aligned$d_y, aligned$d_z)
compute_EI(aligned$d_x, aligned$d_y, aligned$d_z)

# Spherical harmonic decomposition
scar_sh <- spharm_from_directions(aligned, lmax = 20)

# --- Morphological Analysis ---
morph_sh <- spharm_from_meshes("~/examples/meshes", lmax = 20)

# --- Interactive viewer ---
export_spharm_html(
  morph = morph_sh, scar = scar_sh,
  out_path = "~/spharm_viewer.html"
)
RSCRIPT

# --- 8. Fix ownership --------------------------------------------------------
RUN chown -R rstudio:rstudio /home/rstudio/examples

# --- 9. RStudio preferences --------------------------------------------------
RUN mkdir -p /home/rstudio/.config/rstudio && \
    echo '{"initial_working_directory": "/home/rstudio"}' \
      > /home/rstudio/.config/rstudio/rstudio-prefs.json && \
    chown -R rstudio:rstudio /home/rstudio/.config

# --- 10. Default environment variables ---------------------------------------
#   PASSWORD and DISABLE_AUTH are set at runtime via -e or docker-compose.
#   Defaults here for convenience; override with:
#     docker run -e PASSWORD=yourpass ...
ARG DEFAULT_PASSWORD=spharm
ENV PASSWORD=${DEFAULT_PASSWORD}
ENV DISABLE_AUTH=false

EXPOSE 8787