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
    libglpk-dev cmake wget curl git && rm -rf /var/lib/apt/lists/*

# ── Quarto ───────────────────────────────────────────────────────────
ARG QUARTO_VERSION=1.9.36
ARG TARGETARCH
RUN wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb" \
    && dpkg -i quarto-*.deb && rm quarto-*.deb

# ── R / Bioconductor packages ───────────────────────────────────────
ENV R_LIBS_USER=/usr/local/lib/R/site-library

# Use Posit Package Manager for pre-built Linux binaries (CRAN + Bioc)
COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R && rm /tmp/install_r_packages.R

# ── Python: uv + string-gsea from GitHub (includes snakemake) ────────
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/python
RUN uv python install 3.13 \
    && uv venv /opt/venv \
    && uv pip install --python /opt/venv "string-gsea @ git+https://github.com/wolski/stringdbpy.git" \
    && chmod -R a+rX /opt/python /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# ── FGCZ Quarto template — pinned, as prolfquapp pins it in its DESCRIPTION.
#    stringGSEAplot Imports it, and fgcz_render() stages the theme from here at
#    render time, so the styling of a report is fixed by this commit.  ──
ARG FGCZ_QUARTO_TEMPLATE_REF=dbeb852bc14f4374d59ff8b94350faf84f4513e1
RUN R -e "remotes::install_github('fgcz/fgczQuartoTemplate', ref='${FGCZ_QUARTO_TEMPLATE_REF}', dependencies=TRUE, upgrade='never')" \
    && Rscript -e "stopifnot(requireNamespace('fgczQuartoTemplate', quietly=TRUE))"

# ── stringGSEAplot R package — built WITH vignettes and installed from the
#    tarball, which is what activates `vignettes/.install_extras`: the qmd
#    sources and the FGCZ assets land in the installed package's doc/. Building
#    the vignette also renders the packaged example report, so a report that
#    cannot render fails the image build instead of a production run. Its
#    rendered HTML is then deleted: it embeds every figure and runs to ~100 MB,
#    and runtime renders its own report from the real results.  ──
COPY stringGSEAplot /tmp/stringGSEAplot
RUN cd /tmp \
    && R CMD build stringGSEAplot \
    && R CMD INSTALL stringGSEAplot_*.tar.gz \
    && inst=$(Rscript -e "cat(find.package('stringGSEAplot'))") \
    && for asset in GSEA_report.qmd _metadata.yml fgcz.scss fgcz_header_quarto.html fgcz-plot-finder.html; do \
         test -f "$inst/doc/$asset" || { echo "missing from installed doc/: $asset"; exit 1; }; \
       done \
    && rm -f "$inst"/doc/*.html "$inst"/doc/*.R \
    && rm -rf /tmp/stringGSEAplot /tmp/stringGSEAplot_*.tar.gz

WORKDIR /work
ENTRYPOINT ["string_gsea_workflow"]
CMD ["--cores", "1"]
