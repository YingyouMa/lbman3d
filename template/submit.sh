#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"
sbatch --job-name="$(basename "$PWD")" example.sh
