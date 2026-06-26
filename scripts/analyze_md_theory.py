#!/usr/bin/env python3
"""
analyze_md_theory.py - the Exp-1 test: does MD dynamically distinguish nitrated from
non-nitrated aromatics on the acidic+thiol microenvironment mechanism?

Per aromatic site, over the trajectory:
  pcet_mindist   - min distance donor -> nearest carboxylate O (Tyr OH / Trp NE1 to
                   Asp OD1/OD2, Glu OE1/OE2). Proton-acceptor geometry for the
                   tyrosyl/tryptophanyl radical (PCET).
  pcet_frac_HB   - fraction of frames with that distance < 0.35 nm (H-bond competent)
  thiol_mindist  - min distance ring centroid -> nearest Cys SG / Met SD
  sasa_ring      - mean solvent-accessible surface area of the ring atoms (nm^2)

Then nitrated-vs-non summary. Prediction (theory): nitrated sites have SMALLER
pcet_mindist / higher pcet_frac_HB, SMALLER thiol_mindist, and/or larger transient
sasa than non-nitrated aromatics in the same protein.

  python analyze_md_theory.py --top sys.prmtop --traj run.dcd --sites sites.csv [--stride 5]
sites.csv columns: resSeq,resname,nitrated   (resSeq = topology residue number)
"""
import argparse, csv
import numpy as np
import mdtraj as md

RING = {"TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]}
DONOR = {"TYR": "OH", "TRP": "NE1"}
HB = 0.35  # nm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--min_sasa", type=float, default=0.1,
                    help="ring SASA (nm^2) above which an aromatic counts as solvent-accessible")
    a = ap.parse_args()

    t = md.load(a.traj, top=a.top, stride=a.stride)
    # OPC water has a massless virtual site (element 'VS') with no SASA radius, and
    # explicit solvent would mask the protein surface — drop water before any geometry.
    t = t.atom_slice(t.topology.select("not water"))
    top = t.topology
    print(f"loaded {t.n_frames} frames, {t.n_atoms} atoms (water stripped)", flush=True)

    # acceptor/donor pools (protein-wide)
    carbox = top.select("(resname ASP and (name OD1 or name OD2)) or "
                        "(resname GLU and (name OE1 or name OE2))")
    thiol = top.select("(resname CYS and name SG) or (resname MET and name SD)")
    sasa = md.shrake_rupley(t, mode="atom")     # (n_frames, n_atoms) nm^2

    sites = list(csv.DictReader(open(a.sites)))
    use_index = "resindex" in (sites[0] if sites else {})   # prefer map_sites.py output
    rows = []
    for s in sites:
        rn = s["resname"].upper()
        nit = s["nitrated"].strip().lower() in ("1", "true", "yes")
        if use_index:
            res = top.residue(int(s["resindex"]))
            if res.name != rn:
                print(f"  WARN: resindex {s['resindex']} is {res.name}, expected {rn}"); continue
            tag = f"{rn}{s.get('unp_pos', s['resindex'])}"
        else:
            rs = int(s["resSeq"])
            res = next((r for r in top.residues if r.resSeq == rs and r.name == rn), None)
            tag = f"{rn}{rs}"
        if res is None:
            print(f"  WARN: {tag} not found in topology"); continue
        ring_idx = [a_.index for a_ in res.atoms if a_.name in RING[rn]]
        don = next((a_.index for a_ in res.atoms if a_.name == DONOR[rn]), None)

        # PCET: donor -> nearest carboxylate O per frame
        if don is not None and len(carbox):
            d = md.compute_distances(t, np.array([[don, c] for c in carbox]))  # (F, nC)
            pcet = d.min(axis=1)
            pcet_mn, pcet_hb = float(pcet.mean()), float((pcet < HB).mean())
        else:
            pcet_mn, pcet_hb = np.nan, np.nan

        # thiol: ring centroid -> nearest Cys SG / Met SD per frame
        if len(thiol):
            cen = t.xyz[:, ring_idx, :].mean(axis=1)                # (F,3)
            sx = t.xyz[:, thiol, :]                                  # (F,nT,3)
            thd = np.linalg.norm(sx - cen[:, None, :], axis=2).min(axis=1)
            thiol_mn = float(thd.mean())
        else:
            thiol_mn = np.nan

        sasa_ring = float(sasa[:, ring_idx].sum(axis=1).mean())
        rows.append(dict(site=tag, nitrated=nit, pcet_mindist=pcet_mn,
                         pcet_frac_HB=pcet_hb, thiol_mindist=thiol_mn, sasa_ring=sasa_ring))

    print(f"\n{'site':8} {'nit':>4} {'pcet_min(nm)':>13} {'pcet_HB':>8} {'thiol_min(nm)':>14} {'sasa(nm2)':>10}")
    for r in sorted(rows, key=lambda x: (not x["nitrated"], x["site"])):
        print(f"{r['site']:8} {str(r['nitrated']):>4} {r['pcet_mindist']:>13.3f} "
              f"{r['pcet_frac_HB']:>8.3f} {r['thiol_mindist']:>14.3f} {r['sasa_ring']:>10.3f}")

    def summarize(subset, label):
        nit = [r for r in subset if r["nitrated"]]
        non = [r for r in subset if not r["nitrated"]]
        if not (nit and non):
            print(f"\n{label}: need both groups (nitrated={len(nit)}, non={len(non)}) -- skipped")
            return
        print(f"\n=== {label}: nitrated (n={len(nit)}) vs non (n={len(non)}) means ===")
        for k in ["pcet_mindist", "pcet_frac_HB", "thiol_mindist", "sasa_ring"]:
            mn = np.nanmean([r[k] for r in nit]); mo = np.nanmean([r[k] for r in non])
            print(f"  {k:14} nitrated={mn:.3f}  non={mo:.3f}  diff={mn-mo:+.3f}")

    summarize(rows, "all aromatics")
    acc = [r for r in rows if r["sasa_ring"] >= a.min_sasa]
    summarize(acc, f"accessible only (sasa_ring >= {a.min_sasa} nm^2)")
    print("\nTheory predicts nitrated: smaller pcet_mindist, higher pcet_frac_HB,")
    print("smaller thiol_mindist, and/or larger sasa_ring than non-nitrated.")
    print("The accessible-only block is the honest test: among aromatics peroxynitrite can")
    print("actually reach, does the carboxylate/thiol geometry separate nitrated from non?")


if __name__ == "__main__":
    main()
