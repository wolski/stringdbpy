# STRING-GSEA Docker image
# Mixed Python 3.13 + R 4.5 + Quarto + Snakemake
# Multi-arch: linux/amd64, linux/arm64
#
# Build:  docker buildx build -f docker/Dockerfile -t string-gsea:local --load .
# Run:    docker run --rm -v $(pwd)/tests:/work \
#           -v ~/.config/string_gsea:/root/.config/string_gsea:ro \
#           string-gsea:local --cores 1

FROM rocker/r-ver:4.5.0

# ── Layer 1: System dependencies ─────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# ── Layer 2: Quarto ──────────────────────────────────────────────────────────
ARG QUARTO_VERSION=1.9.36
ARG TARGETARCH
RUN wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb" \
    && dpkg -i quarto-*.deb \
    && rm quarto-*.deb

# ── Layer 3: R / Bioconductor packages (~20 min, cached aggressively) ────────
# Install everything to the system library
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('DOSE', 'enrichplot'), ask=FALSE, update=FALSE)" \
    && R -e "install.packages(c('ggplot2', 'ggupset', 'patchwork', 'DT', 'jsonlite', 'knitr', 'rmarkdown', 'quarto'), repos='https://cloud.r-project.org')"

# ── Layer 4: uv + Python venv with dependencies ─────────────────────────────
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
WORKDIR /build
RUN uv python install 3.13

# Install deps first (cached unless pyproject.toml/uv.lock change)
COPY pyproject.toml uv.lock README.md /build/
RUN uv sync --frozen --no-install-project

# ── Layer 5: Python source + install into venv ───────────────────────────────
COPY src/ /build/src/
RUN uv sync --frozen

# Install snakemake into the same venv
RUN uv pip install snakemake

# Put venv on PATH so all CLI entry points are available
ENV PATH="/build/.venv/bin:$PATH"

# ── Layer 6: R package (stringGSEAplot) ──────────────────────────────────────
COPY stringGSEAplot/ /tmp/stringGSEAplot/
# Copy vignette QMD into templates so system.file("templates", ...) finds it
RUN cp /tmp/stringGSEAplot/vignettes/GSEA_report.qmd /tmp/stringGSEAplot/inst/templates/ \
    && R CMD INSTALL /tmp/stringGSEAplot/ \
    && rm -rf /tmp/stringGSEAplot/

# ── Layer 7: Snakefile + entrypoint ──────────────────────────────────────────
COPY tests/Snakefile /app/Snakefile
COPY docker/entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

WORKDIR /work
ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["--cores", "1"]
