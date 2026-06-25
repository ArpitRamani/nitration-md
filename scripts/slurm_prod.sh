#!/usr/bin/env bash
#SBATCH --job-name=md_prod
#SBATCH --partition=l4-4-gm96-c48-m192
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-9%4
#SBATCH --output=/users/arama30/prod_%A_%a.log
# GPU production fan-out: each array task = one (system, replica), 1 GPU each.
# Override partition/array/concurrency on the CLI, e.g. for the 8-GPU RTX-Pro node:
#   sbatch --partition=rp6b-8-gm768-c192-m2048 --array=1-9%8 scripts/slurm_prod.sh
# jobs.txt line format:  <prmtop> <inpcrd> <out_prefix> <seed> <ns> <temp>  (paths relative to ROOT)
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nitromd

ROOT=/scratch/arama30/mol_dynamics       # data root: holds chi3l1/ lyz/ ttr/ and scripts/
cd "$ROOT"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ROOT/scripts/jobs.txt")
[ -z "$LINE" ] && { echo "no jobs.txt line ${SLURM_ARRAY_TASK_ID}"; exit 1; }
read -r PRM CRD OUT SEED NS TEMP <<< "$LINE"
echo "task ${SLURM_ARRAY_TASK_ID}: $OUT  (seed $SEED, $NS ns, $TEMP K)  GPU=$CUDA_VISIBLE_DEVICES  node=$(hostname)"

python "$ROOT/scripts/run_openmm.py" \
    --prmtop "$PRM" --inpcrd "$CRD" --out "$OUT" --seed "$SEED" --ns "$NS" --temp "$TEMP"
