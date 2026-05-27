#!/bin/bash

#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH -t 5-00:00:00
#SBATCH --mail-user=yingyouma@hpcc.brandeis.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

set -euo pipefail

if command -v module >/dev/null 2>&1; then
    module purge || true
    # Edit these module loads if your HPCC toolchain differs.
    module load cmake || true
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OMP_PLACES=threads
export OMP_PROC_BIND=close
export OMP_NESTED=FALSE

cd "$SLURM_SUBMIT_DIR"
mkdir -p output
cmake -S . -B build-case
cmake --build build-case -j "${SLURM_CPUS_PER_TASK:-16}"
./build-case/main > output.resume 2> errors.resume
