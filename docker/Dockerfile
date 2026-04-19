# STRING-GSEA Docker image
# Mixed Python 3.13 + R 4.5 + Quarto + Snakemake
# Multi-arch: linux/amd64, linux/arm64
#
# Build:  docker buildx build -f docker/Dockerfile -t string-gsea:local --load .
# Run:    docker run --rm -v $(pwd)/tests:/work \
#           -v ~/.config/string_gsea:/root/.config/string_gsea:ro \
#           string-gsea:local --cores 1

FROM rocker/r-ver:4.5.3

# ── System dependencies ─────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    cmake wget curl git && rm -rf /var/lib/apt/lists/*

# ── Quarto ───────────────────────────────────────────────────────────
ARG QUARTO_VERSION=1.9.36
ARG TARGETARCH
RUN wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb" \
    && dpkg -i quarto-*.deb && rm quarto-*.deb

# ── R / Bioconductor packages ───────────────────────────────────────
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Use Posit Package Manager for pre-built Linux binaries (CRAN + Bioc)
ENV PPM_CRAN="https://packagemanager.posit.co/cran/__linux__/noble/latest"
ENV PPM_BIOC="https://packagemanager.posit.co/bioconductor/3.22"
RUN Rscript -e "\
    repos <- c(CRAN='${PPM_CRAN}', BioCsoft='${PPM_BIOC}', BioCann='${PPM_BIOC}'); \
    options(repos=repos); \
    install.packages(c('BiocManager','remotes','Rcpp','yulab.utils','ggplot2','ggupset', \
      'patchwork','DT','jsonlite','knitr','rmarkdown','quarto','dplyr','tidyr', \
      'tibble','pillar','cli','rlang','vctrs','lifecycle','scales','reshape2', \
      'igraph','purrr','plyr','stringr','ggnewscale','ggrepel','ggfun','ggplotify', \
      'ggforce','gridGraphics','cowplot','aplot')); \
    install.packages(c('AnnotationDbi','GO.db','GOSemSim','BiocParallel','fgsea','DOSE','enrichplot')); \
    stopifnot(requireNamespace('DOSE'), requireNamespace('enrichplot'))"

# ── Python: uv + string-gsea from GitHub (includes snakemake) ────────
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/python
RUN uv python install 3.13 \
    && uv venv /opt/venv \
    && uv pip install --python /opt/venv "string-gsea @ git+https://github.com/wolski/stringdbpy.git" \
    && chmod -R a+rX /opt/python /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# ── stringGSEAplot R package — install from local source in build context.
#    Skip vignette compilation (fragile); copy QMD sources manually into doc/.  ──
COPY stringGSEAplot /tmp/stringGSEAplot
RUN R -e "remotes::install_local('/tmp/stringGSEAplot', build_vignettes=FALSE, dependencies=FALSE, force=TRUE)" \
    && inst=$(Rscript -e "cat(find.package('stringGSEAplot'))") \
    && mkdir -p "$inst/doc" \
    && cp /tmp/stringGSEAplot/vignettes/GSEA_report.qmd "$inst/doc/" \
    && rm -rf /tmp/stringGSEAplot

WORKDIR /work
ENTRYPOINT ["string_gsea_workflow"]
CMD ["--cores", "1"]
