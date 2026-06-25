#!/usr/bin/env bash
#SBATCH --job-name=md_analyze
#SBATCH --partition=c64-m512
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/users/arama30/analyze_%j.log
# Exp-1 theory test on every completed replica: PCET distance / thiol proximity /
# ring SASA, nitrated-vs-non. Run once all production runs finish (squeue empty).
#   sbatch scripts/slurm_analyze.sh
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nitromd

ROOT=/scratch/arama30/mol_dynamics
cd "$ROOT"
mkdir -p results
for sys in chi3l1 lyz ttr; do
  sites="$sys/sites_top.csv"
  top="$sys/sys.prmtop"
  [ -f "$sites" ] || { echo "skip $sys: no $sites"; continue; }
  [ -f "$top" ]   || { echo "skip $sys: no $top";   continue; }
  for dcd in "$sys"/*_r*.dcd; do
    [ -f "$dcd" ] || continue
    name=$(basename "$dcd" .dcd)
    echo "==================== $name ===================="
    python scripts/analyze_md_theory.py \
        --top "$top" --traj "$dcd" --sites "$sites" --stride 10 \
        | tee "results/${name}_theory.txt"
  done
done
echo "done -> results/*_theory.txt"
