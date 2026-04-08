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
# Override rocker's pinned CRAN snapshot with live CRAN to avoid
# stale dependency versions (e.g. yulab.utils 0.2.0 vs >= 0.2.1 needed by GOSemSim)
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /usr/local/lib/R/etc/Rprofile.site

RUN R -e "install.packages('BiocManager')" \
    && R -e "BiocManager::install(c('DOSE', 'enrichplot'), ask=FALSE, update=TRUE)" \
    && R -e "install.packages(c('ggplot2','ggupset','patchwork','DT','jsonlite','knitr','rmarkdown','quarto','remotes'))"

# ── stringGSEAplot R package (from GitHub, build_vignettes so doc/ has QMD sources) ──
RUN R -e "remotes::install_github('wolski/stringdbpy', subdir='stringGSEAplot', build_vignettes=TRUE, repos=BiocManager::repositories())"

# ── Python: uv + string-gsea from GitHub (includes snakemake) ────────
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/python
RUN uv python install 3.13 \
    && uv venv /opt/venv \
    && uv pip install --python /opt/venv "string-gsea @ git+https://github.com/wolski/stringdbpy.git" \
    && chmod -R a+rX /opt/python /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

WORKDIR /work
ENTRYPOINT ["string_gsea_workflow"]
CMD ["--cores", "1"]
