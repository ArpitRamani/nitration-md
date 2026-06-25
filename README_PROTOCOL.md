# Nitration MD — protocol

Tests the two-cohort microenvironment theory (acidic+thiol → nitration) on the OG-paper
proteins, and the TTR functional-consequence hypothesis. Toolchain: **AmberTools build
(tleap/antechamber) → OpenMM production**. No CHARMM-GUI. FF: ff19SB + OPC water, 0.15 M NaCl.

## Two experiments
- **Exp 1 — theory test (native MD, NO new params; start now).** Per aromatic, compare
  nitrated vs non-nitrated sites IN THE SAME protein on: (a) carboxylate→Tyr-OH / Trp-NE1
  proton-transfer distance (PCET geometry), (b) nearest Cys-SG/Met-SD approach, (c) ring
  SASA fluctuation (transient exposure). Targets (chosen for nitrated/non CONTRAST in our
  data): **CHI3L1 (11 nit / 13 non — best), human LYZ (3 / 8 — fast)**.
  NOT Hsp90 (Y33/Y56 are tissue sites, absent from our CSF data) and NOT TTR (7/7 nitrated,
  no contrast). Hsp90 N-domain optional as a literature-only positive control.
- **Exp 2 — TTR consequence (needs nitro params).** Tetramer ± nitration; W61 (PDB W41)
  is the validated AD-dominant site (ADNI rank 1/7, all channels) -> nitro-Trp (derive
  early, critical path). Y136 (PDB Y116) secondary -> nitro-Tyr via forcefield_PTM.
  Measure dimer–dimer interface separation, T4-pocket geometry, monomer-unfolding (does
  nitration destabilize the amyloid-protective tetramer?).

## Starting structures (fetch + clean; see scripts/01_prep_structure.md)
| protein | PDB | note |
|---|---|---|
| Hsp90α N-domain | 1YER or 3T0Z | contains Y33, Y56 (UniProt P07900 numbering) |
| Human lysozyme (LYZ) | 1LZ1 or 2NWD | small, fast; P61626 |
| CHI3L1 (YKL-40) | 1HJX / 1NWR | glyco — strip glycans for Exp 1 |
| **TTR tetramer** | **1F41 or 1TTA** | WT homotetramer; ± T4 in 1F41 |

UniProt/precursor numbering throughout (matches the MS data): TTR W61=mature W41, Y136=mature Y116.

## Pipeline (run order)
    # --- system build (HPC; module load ambertools) ---
    # 1. clean PDB: remove waters/hetero (keep T4 for TTR if modeling holo), add H, fix termini
    pdb4amber -i raw.pdb -o clean.pdb --reduce --dry
    # 2. (Exp 2 only) build nitro-Tyr / nitro-Trp residue libs once:  docs/nitro_params.md
    # 3. tleap build -> prmtop/inpcrd
    tleap -f scripts/build_native.tleap        # Exp 1
    tleap -f scripts/build_nitro_ttr.tleap     # Exp 2 (after params)

    # --- production (HPC GPU; module load openmm/cuda) ---
    python scripts/run_openmm.py --prmtop sys.prmtop --inpcrd sys.inpcrd \
        --out PROT_rep1 --ns 200 --temp 300
    # 3 replicas per system (vary --seed); TTR ± nitration each x3.

    # --- analysis ---
    python scripts/analyze_md_theory.py --top sys.prmtop --traj PROT_rep1.dcd \
        --sites sites_PROT.csv          # PCET / thiol / SASA per aromatic, nitrated vs non

## Simulation settings (in run_openmm.py)
PME, 1.0 nm cutoff, HBonds constraints, hydrogen-mass repartitioning → 4 fs dt; min →
NVT heat (10→300 K) → NPT equil (1 ns) → NPT production (default 200 ns). 3 replicas
(different seeds). On the cluster: 1 GPU per replica, embarrassingly parallel.

## Cluster (Emory AWS ParallelCluster, Slurm)
One-time env (conda-forge OpenMM+AmberTools+MDTraj, avoids module guessing):
    bash scripts/env_setup.sh        # creates env 'nitromd'

Partition mapping:
- build (tleap/antechamber/pdb4amber): `c64-m512` (CPU) -> `scripts/slurm_build.sh`
- production: `l4-4-gm96-c48-m192` (4x L4, 1 GPU/replica, %4 concurrent) -> `scripts/slurm_prod.sh`
- burst fan-out: `rp6b-8-gm768-c192-m2048` (8x L40S, always-on; set partition + `%8`)
- skip B200 (180 GB/GPU wasted on ~100k-atom systems; needs bleeding-edge CUDA/OpenMM)

Each system lives in its own subdir (hsp90/, lyz/, ttr/) with sys.prmtop/sys.inpcrd.
Workflow:
    # build each system (CPU)
    cd hsp90 && sbatch ../scripts/slurm_build.sh hsp90_raw.pdb ../scripts/build_native.tleap
    # fan out production (GPU array; edit jobs.txt + --array range)
    sbatch scripts/slurm_prod.sh
    # 3 replicas/system already in scripts/jobs.txt (9 lines = Exp-1 + TTR native)

Sizing: ~200 ns in <1 day (L40S) to ~1-2 days (L4) per replica; 7-day wall limit, with
checkpoints (.chk) written for restart. All systems fit easily in 24 GB (L4).

## sites_PROT.csv (you provide per protein; I generate from cohort-1 data)
    resid,resSeq,resname,nitrated     # resSeq = UniProt position; nitrated from MS data
The MD-target list (which aromatics are nitrated in each protein) comes from
`adni_replication` / cohort-1 — I'll emit these with the static features for target picking.
