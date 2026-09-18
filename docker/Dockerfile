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

# ── Python: uv + string-gsea from THIS build context (includes snakemake).
#    Installed from the sources being built, not re-fetched from GitHub: a
#    `git+https://` install never changes this instruction's text, so Docker and
#    the GHA cache reused the layer and the image shipped whatever string-gsea
#    the cache happened to hold -- one release behind, for as long as the cache
#    lived. COPYing the sources ties the layer to their content instead.  ──
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/python
COPY pyproject.toml uv.lock README.md LICENSE.txt /src/
COPY src /src/src
RUN uv python install 3.13 \
    && uv venv /opt/venv \
    && uv pip install --python /opt/venv /src \
    && chmod -R a+rX /opt/python /opt/venv \
    && /opt/venv/bin/python -c "import importlib.metadata as meta, pathlib, tomllib; \
source = tomllib.loads(pathlib.Path('/src/pyproject.toml').read_text())['project']['version']; \
installed = meta.version('string-gsea'); \
assert installed == source, f'installed {installed}, sources say {source}'" \
    && rm -rf /src
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
#    rendered HTML is deleted BY NAME -- `doc/*.html` would also take
#    fgcz_header_quarto.html and fgcz-plot-finder.html, which are assets the
#    report includes, not output -- and the assets are checked afterwards, so
#    the check describes what the image actually ships.  ──
COPY stringGSEAplot /tmp/stringGSEAplot
RUN cd /tmp \
    && R CMD build stringGSEAplot \
    && R CMD INSTALL stringGSEAplot_*.tar.gz \
    && inst=$(Rscript -e "cat(find.package('stringGSEAplot'))") \
    && rm -f "$inst"/doc/GSEA_report.html "$inst"/doc/GSEA_report.R \
    && for asset in GSEA_report.qmd _metadata.yml fgcz.scss fgcz_header_quarto.html fgcz-plot-finder.html; do \
         test -f "$inst/doc/$asset" || { echo "missing from installed doc/: $asset"; exit 1; }; \
       done \
    && rm -rf /tmp/stringGSEAplot /tmp/stringGSEAplot_*.tar.gz

# ── Shipped-state check. The render step crosses from Python into R, so the
#    two halves have to agree in the image: a stale Python layer or a missing R
#    export is a broken image, not a broken run. The doc/ listing is printed so
#    a build log always records what the image ships.  ──
RUN python -c "import importlib.metadata as meta; print('string-gsea', meta.version('string-gsea'))" \
    && Rscript -e "stopifnot('render_gsea_reports' %in% getNamespaceExports('stringGSEAplot'))" \
    && ls -1 "$(Rscript -e "cat(find.package('stringGSEAplot'))")/doc"

WORKDIR /work
ENTRYPOINT ["string_gsea_workflow"]
CMD ["--cores", "1"]
