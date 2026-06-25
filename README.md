# nitration-md

Molecular dynamics pipeline for testing the **acidic + thiol microenvironment / PCET
hypothesis** of site-selective protein nitration (3-nitrotyrosine, nitrotryptophan) in
CSF proteins, using AMBER **ff19SB / OPC** built with AmberTools and run on **OpenMM**
(CUDA, 4 fs via hydrogen-mass repartitioning).

## The hypothesis

Peroxynitrite-driven tyrosyl/tryptophanyl-radical nitration proceeds by proton-coupled
electron transfer (PCET). If nitration is *functionally selective* rather than random
oxidative damage, the residues that get nitrated should sit in a distinctive local
microenvironment — a nearby carboxylate (Asp/Glu) to accept the phenolic/indole proton,
and a nearby thiol/thioether (Cys/Met) redox partner. This repo tests that prediction in
explicit-solvent dynamics.

## Two experiments

1. **Microenvironment test (native dynamics).** Simulate the *unmodified* protein and
   ask, within the same trajectory, whether the aromatics that are observed nitrated sit
   in PCET-competent geometry (short donor→carboxylate distance, high H-bond occupancy,
   close thiol) versus aromatics that are not nitrated. No PTM parameters required.
2. **Consequence test (nitrated dynamics).** Build the modified residue and compare
   native vs nitrated to probe the structural consequence (e.g. transthyretin tetramer-
   interface destabilization). Nitro-Tyr / nitro-Trp use **published** AMBER-compatible
   parameters (NIY / NIW, 6-nitrotryptophan as the dominant peroxynitrite isomer); no
   force field is defined here.

## Pipeline

| Stage | Script | What it does |
|---|---|---|
| Build | `scripts/build_system.py` | pdb4amber → strip to standard-AA heavy atoms → renumber → auto-detect disulfides (CYX) → tleap (ff19SB/OPC, neutralize + 20 mM salt) → `sys.prmtop`/`sys.inpcrd` |
| Assembly | `scripts/flatten_assembly.py` | Collapse a multi-MODEL biological assembly (`.pdb1`) to one model (e.g. TTR tetramer) |
| Map | `scripts/map_sites.py` | Align built structure to UniProt, emit topology-indexed sites file (`resindex,resname,nitrated,chain,unp_pos`) |
| Produce | `scripts/run_openmm.py` | min → NVT heat → NPT equil → NPT production; PME, HBonds + HMR (4 fs), Langevin + MonteCarlo barostat; DCD + checkpoints |
| Analyze | `scripts/analyze_md_theory.py` | Per-aromatic PCET donor→carboxylate distance, H-bond fraction, ring→thiol distance, ring SASA; nitrated-vs-non summary |

Slurm wrappers (`slurm_build.sh`, `slurm_prod.sh` GPU array, `slurm_analyze.sh`,
`slurm_test.sh`) and `env_setup.sh` (conda env) are provided for an HPC cluster.

## Requirements

`conda create -n nitromd -c conda-forge python=3.11 openmm ambertools mdtraj numpy scipy parmed`
(pin `cuda-version<=<your driver's max>` so OpenMM's CUDA build matches the GPU driver).

## Usage

```bash
# build one system (optionally keep a single chain)
cd chi3l1 && python ../scripts/build_system.py 1hjx.pdb A

# map MS-observed nitration sites onto the topology
python scripts/map_sites.py --built chi3l1/sys_built.pdb --uniprot P36222 \
    --sites md_sites/sites_CHI3L1.csv --out chi3l1/sites_top.csv

# production (1 GPU)
python scripts/run_openmm.py --prmtop chi3l1/sys.prmtop --inpcrd chi3l1/sys.inpcrd \
    --out chi3l1/chi3l1_r1 --seed 1 --ns 200 --temp 300

# theory test
python scripts/analyze_md_theory.py --top chi3l1/sys.prmtop \
    --traj chi3l1/chi3l1_r1.dcd --sites chi3l1/sites_top.csv --stride 10
```

## Data note

The per-residue nitration labels (`md_sites/*.csv`) are derived from an unpublished study
and are **not** distributed here. Provide your own site table (columns `resSeq,resname,
nitrated`) to `scripts/map_sites.py`. Starting structures are public PDB entries
(CHI3L1 `1HJX`, lysozyme `1LZ1`, transthyretin `1TTA`).

## License

MIT
