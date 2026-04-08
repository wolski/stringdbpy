# STRING-GSEA Docker image
# Mixed Python 3.13 + R 4.5 + Quarto + Snakemake
# Multi-arch: linux/amd64, linux/arm64
#
# Build:  docker buildx build -f docker/Dockerfile -t string-gsea:local --load .
# Run:    docker run --rm -v $(pwd)/tests:/work \
#           -v ~/.config/string_gsea:/root/.config/string_gsea:ro \
#           string-gsea:local --cores 1

FROM rocker/r-ver:4.5.0

# ── System dependencies ─────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    wget curl git && rm -rf /var/lib/apt/lists/*

# ── Quarto ───────────────────────────────────────────────────────────
ARG QUARTO_VERSION=1.9.36
ARG TARGETARCH
RUN wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb" \
    && dpkg -i quarto-*.deb && rm quarto-*.deb

# ── R / Bioconductor packages ───────────────────────────────────────
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('DOSE', 'enrichplot'), ask=FALSE, update=TRUE)" \
    && R -e "install.packages(c('ggplot2','ggupset','patchwork','DT','jsonlite','knitr','rmarkdown','quarto','remotes'), repos='https://cloud.r-project.org')"

# ── stringGSEAplot R package (from GitHub, build_vignettes so doc/ has QMD sources) ──
RUN R -e "options(repos = BiocManager::repositories()); remotes::install_github('wolski/stringdbpy', subdir='stringGSEAplot', build_vignettes=TRUE)"

# ── Python: uv + string-gsea from GitHub (includes snakemake) ────────
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
RUN uv python install 3.13 \
    && uv venv /opt/venv \
    && uv pip install --python /opt/venv "string-gsea @ git+https://github.com/wolski/stringdbpy.git"
ENV PATH="/opt/venv/bin:$PATH"

WORKDIR /work
ENTRYPOINT ["string_gsea_workflow"]
CMD ["--cores", "1"]
