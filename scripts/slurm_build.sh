#!/usr/bin/env bash
#SBATCH --job-name=md_build
#SBATCH --partition=c64-m512
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=build_%j.log
# CPU system build. Usage:
#   sbatch slurm_build.sh <raw.pdb> <tleap_input>    e.g. sbatch slurm_build.sh ttr_raw.pdb build_native.tleap
# Produces sys.prmtop / sys.inpcrd in the current directory.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate nitromd

RAW="$1"; TLEAP="$2"

# clean: add H (reduce), strip waters/heteroatoms (keep ligand manually if modeling holo),
# flag disulfides as CYX. Inspect clean.pdb before tleap (chain breaks, missing residues).
pdb4amber -i "$RAW" -o clean.pdb --reduce --dry

tleap -f "$TLEAP"
echo "build done: $(ls -1 sys.prmtop sys.inpcrd 2>/dev/null || echo 'MISSING — check tleap log')"
