#!/usr/bin/env bash
# One-time: create the MD conda env in YOUR HOME (the shared /opt/anaconda base is
# read-only; new envs go to ~/.conda/envs which you own). Run on a login/CPU node.
#   bash env_setup.sh        (or: nohup bash env_setup.sh > env_setup.log 2>&1 &)
set -euo pipefail

# force env + package cache into home (avoids the read-only system base entirely)
conda config --add pkgs_dirs "$HOME/.conda/pkgs" 2>/dev/null || true
conda config --add envs_dirs "$HOME/.conda/envs" 2>/dev/null || true

SOLVER=$(command -v mamba || command -v conda)   # mamba if present (faster), else conda

"$SOLVER" create -y -n nitromd -c conda-forge \
    python=3.11 openmm ambertools mdtraj numpy scipy parmed

echo
echo "created env in: $HOME/.conda/envs/nitromd"
echo "activate with:  conda activate nitromd"
echo "verify GPU (on a GPU node):  conda activate nitromd && python -m openmm.testInstallation"
