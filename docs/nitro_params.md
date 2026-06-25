# Parameterizing 3-nitro-Tyr (NIY) and 6-nitro-Trp (NIW) for Amber/OpenMM

Needed only for Exp 2 (TTR ± nitration). This is the standard "non-standard amino acid"
workflow: derive sidechain charges by RESP with the BACKBONE atoms charge-constrained to
the standard ff19SB values, so the residue bonds correctly in-chain. Done once; reused.

## Option A (USE THIS for nitro-Tyr): forcefield_PTM
nitro-Tyr (2/3 of our sites, incl. TTR Y136 + Hsp90 Y33/Y56) does NOT need derivation.
- **forcefield_PTM** (Khoury et al., JCTC 2013) — AMBER-compatible RESP params for ~32
  PTMs; 3-nitrotyrosine is included. Ships as Amber lib/frcmod -> drops straight into
  tleap (loadoff + loadamberparams), no CHARMM-GUI. Verify current download (authors'
  site / paper SI / maintained mirror); confirm residue name + net charge = 0 before use.
- Alternatives if needed: SwissSidechain (NCAA params), or CHARMM36 + PyTM/CHARMM-GUI PTM
  Manager (only if going the GROMACS/CHARMM route).
Drop the nitro-Tyr lib/frcmod into `params/` and skip to tleap.

nitro-Trp (6-nitro-Trp, e.g. TTR W61) is rarely parameterized and likely NOT in these
sets -> derive it once via Option B. It is not on the critical path: start nitro-Tyr
(Y136) TTR runs with forcefield_PTM while deriving nitro-Trp in parallel.

## Option B: derive (AmberTools)
    # 1. Build a capped residue ACE-X-NME where X is the nitrated aromatic.
    #    Start from a Tyr/Trp dipeptide PDB; add the nitro group: an aromatic ring H is
    #    replaced by N(=O)O. Tyr -> NO2 on CE1 (ring position 3); Trp -> NO2 on CZ2 (pos 6).
    #    Build with a GUI (Avogadro) or PyMOL, save capped_niy.pdb / capped_niw.pdb.

    # 2. Geometry-optimize + ESP at HF/6-31G* (Gaussian or psi4), then RESP with backbone
    #    charge restraints (two-stage RESP; constrain ACE+backbone N,H,CA,HA,C,O to ff19SB):
    antechamber -i capped_niy.pdb -fi pdb -o niy.mol2 -fo mol2 -c resp -at amber -rn NIY
    parmchk2 -i niy.mol2 -f mol2 -o niy.frcmod      # fill missing nitro-group params

    # 3. Make an off library with correct head/tail for in-chain linkage:
    tleap -f - <<'EOF'
    source leaprc.protein.ff19SB
    NIY = loadmol2 niy.mol2
    set NIY head NIY.1.N
    set NIY tail NIY.1.C
    set NIY.1 connect0 NIY.1.N
    set NIY.1 connect1 NIY.1.C
    loadamberparams niy.frcmod
    saveoff NIY params/niy.lib
    quit
    EOF
    # repeat for NIW.

## Validation before using in TTR
- Net charge of NIY/NIW residue = 0 (zwitterion neutral in-chain); nitro group ~ -0.0.
- Short (5 ns) MD of a capped NIY in water: ring planar, nitro group stays conjugated
  (dihedral O-N-C-C oscillates near 0/180, not freely rotating off-plane).
- Compare backbone charges to standard TYR/TRP (should match within rounding).

Then `tleap -f scripts/build_nitro_ttr.tleap` substitutes NIY at 136 / NIW at 61 in the
chains being nitrated (one chain, two chains, or all four — run each stoichiometry).
