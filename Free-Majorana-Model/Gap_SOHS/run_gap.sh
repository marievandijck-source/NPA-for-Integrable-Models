#!/bin/bash

#SBATCH --job-name=GapScaling
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/gap_%j.out
#SBATCH --error=logs/gap_%j.err

set -euo pipefail

module load env/software/doduo
module load Julia/1.11.6-linux-x86_64

cd "$SLURM_SUBMIT_DIR"
echo "Working directory: $(pwd)"

# Needed for Plots/GR on a compute node without a display
export GKSwstype=100

julia --project=. Dev.jl
julia --project=. Logs.jl
julia --project=. Power.jl

echo "Job finished successfully"