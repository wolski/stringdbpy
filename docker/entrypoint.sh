#!/bin/bash
set -euo pipefail

CONFIG_DIR="$HOME/.config/string_gsea"

# Dispatch on first argument
case "${1:-}" in
  config)
    # Generate config.toml — needs the config dir mounted as writable
    exec string_gsea_write_config
    ;;
  *)
    # Default: run snakemake
    if [ ! -f "$CONFIG_DIR/config.toml" ]; then
        echo "ERROR: config.toml not found."
        echo "Either mount it:  -v ~/.config/string_gsea:/root/.config/string_gsea:ro"
        echo "Or generate it:   docker run -it -v ~/.config/string_gsea:/root/.config/string_gsea IMAGE config"
        exit 1
    fi
    cd /work
    exec snakemake -s /app/Snakefile "$@"
    ;;
esac
