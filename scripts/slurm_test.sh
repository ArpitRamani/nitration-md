#!/usr/bin/env bash
#SBATCH --job-name=gpu_test
#SBATCH --partition=l4-4-gm96-c48-m192
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=gpu_test_%j.log
# Batch GPU smoke test — verifies env + CUDA without holding an interactive session.
#   sbatch slurm_test.sh   ; then:  cat gpu_test_<jobid>.log
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nitromd
echo "=== node: $(hostname) ==="
nvidia-smi
echo "=== openmm install check (expect CUDA platform, all pass) ==="
python -m openmm.testInstallation
