#!/bin/bash
# Run string_gsea_workflow in the STRING-GSEA Docker image.
#
# The current directory is mounted to /work and set as the working directory.
# All paths passed to the container must be relative to the current directory.
# The STRING-GSEA config (~/.config/string_gsea/config.toml) is mounted read-only.
# The container runs as the current user (UID/GID) so output files are owned by you.
#
# Prerequisites:
#   1. Docker (or Podman) installed
#   2. Config file at ~/.config/string_gsea/config.toml
#      (generate with: ./string_gsea_docker.sh config)
#
# Setup (first time):
#   ./string_gsea_docker.sh config
#
# Run pipeline on a single dataset:
#   ./docker/string_gsea_docker.sh data/input.zip WU123 results --cores 4
#
# RNK input (no XLSX):
#   ./docker/string_gsea_docker.sh data/input.zip WU123 results --which none --cores 2
#
# With custom FDR threshold:
#   ./docker/string_gsea_docker.sh data/input.zip WU123 results --fdr 0.1 --cores 4
#
# Dry-run (show what would execute without running):
#   ./docker/string_gsea_docker.sh data/input.zip WU123 results --dry-run
#
# Use a specific image version:
#   ./docker/string_gsea_docker.sh --image-version 0.1.0 -- data/input.zip WU123 results --cores 4
#
# Limitations:
#   - Paths must be relative to cwd and must not go outside of it.
#   - Quarto report rendering requires the stringGSEAplot R package (included in the image).

set -euo pipefail

# Default values
IMAGE_VERSION="latest"
IMAGE_REPO="docker.io/wolski/string-gsea"

usage() {
  echo "Usage: $0 [--image-version VERSION] [--image-repo REPO] [-- container_args]"
  echo ""
  echo "Options:"
  echo "  --image-version   Image tag to use (default: $IMAGE_VERSION)"
  echo "  --image-repo      Image repository (default: $IMAGE_REPO)"
  echo "  --help            Show this help message"
  echo ""
  echo "Everything after known options (or after --) is passed to string_gsea_workflow."
  echo ""
  echo "Examples:"
  echo "  $0 data/input.zip WU123 results --cores 4"
  echo "  $0 data/input.zip WU123 results --which none --cores 2"
  echo "  $0 data/input.zip WU123 results --dry-run"
  echo "  $0 config"
  exit 1
}

if [ "$#" -eq 0 ]; then
  usage
fi

# Parse our options, pass the rest through
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --image-version)
            IMAGE_VERSION="$2"
            shift 2
            ;;
        --image-repo)
            IMAGE_REPO="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        --)
            shift
            break
            ;;
        *)
            break
            ;;
    esac
done

CONTAINER_ARGS=("$@")

run() {
    local image_version="$1"
    local image_repo="$2"
    shift 2
    local container_args=("$@")

    # Prefer podman if available
    if command -v podman > /dev/null; then
        DOCKER="podman"
    else
        DOCKER="docker"
    fi

    local image="${image_repo}:${image_version}"
    echo "Using $DOCKER to run $image"

    # Pull if not available locally
    if ! $DOCKER image inspect "$image" > /dev/null 2>&1; then
        echo "Image $image not found locally. Pulling..."
        $DOCKER pull "$image"
    else
        echo "Image $image already exists locally."
        echo "To update: $DOCKER pull $image"
    fi

    # TTY detection
    if [ -t 0 ]; then
        docker_args="-it"
    else
        docker_args="-i"
    fi

    # Config mount: writable for "config" command, read-only otherwise
    local config_dir="${HOME}/.config/string_gsea"
    local config_mount_opts="type=bind,source=${config_dir},target=/root/.config/string_gsea"
    if [[ "${container_args[0]:-}" == "config" ]]; then
        mkdir -p "$config_dir"
    else
        config_mount_opts="${config_mount_opts},readonly"
    fi

    # Mount cwd as /work (read-write), run as current user
    $DOCKER run \
        --init \
        --user "$(id -u):$(id -g)" \
        --rm $docker_args \
        --mount "type=bind,source=$(pwd),target=/work" \
        --mount "$config_mount_opts" \
        -w /work \
        "$image" \
        "${container_args[@]}"
}

run "$IMAGE_VERSION" "$IMAGE_REPO" "${CONTAINER_ARGS[@]}"
