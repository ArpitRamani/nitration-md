#!/usr/bin/env python3
"""
make_md_sites.py - bridge from the MS data to MD: for each notable protein, emit the
aromatic inventory (nitrated vs non, from cohort-1) plus the static acidic/thiol/
accessibility features, so we (a) pick MD targets and (b) feed analyze_md_theory.py.

Reads cohort-1 results/site_features.csv. Writes one sites_<GENE>.csv per protein with
  resSeq,resname,nitrated   (the columns analyze_md_theory.py needs)
and prints a ranked preview (nitrated sites + their static features) for target picking.

  python make_md_sites.py            # writes to ./md_sites/
NOTE: resSeq = UniProt position. After building the system, confirm the topology keeps
UniProt numbering (pdb4amber --reduce + tleap preserve it); if it renumbers, remap.
"""
import csv, os

C1 = "/Users/ramani/Desktop/nitration-pipeline/results/site_features.csv"
NOTABLE = ["TTR", "C3", "CHI3L1", "LYZ", "HSP90AA1", "HSP90AB1"]
THREE = {"Y": "TYR", "W": "TRP"}
OUT = "md_sites"


def flank_chem(fl):
    if not fl:
        return (None, None)
    n = len(fl)
    acid = sum(c in "DE" for c in fl) / n
    cys = sum(c in "CM" for c in fl) / n
    return acid, cys


def main():
    os.makedirs(OUT, exist_ok=True)
    by_gene = {g: [] for g in NOTABLE}
    for r in csv.DictReader(open(C1)):
        g = r.get("gene", "")
        if g in by_gene and r["residue"] in THREE:
            acid, cys = flank_chem(r.get("flank", ""))
            by_gene[g].append(dict(
                resSeq=int(r["position"]), residue=r["residue"],
                nitrated=r["nitrated"] == "True",
                rsa=r.get("rsa", ""), flk_acidic=acid, flk_cysmet=cys))
    for g, sites in by_gene.items():
        if not sites:
            print(f"{g}: not in cohort-1 feature table"); continue
        sites.sort(key=lambda x: x["resSeq"])
        with open(f"{OUT}/sites_{g}.csv", "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["resSeq", "resname", "nitrated"])
            for s in sites:
                w.writerow([s["resSeq"], THREE[s["residue"]], s["nitrated"]])
        nit = [s for s in sites if s["nitrated"]]
        print(f"\n=== {g}: {len(sites)} aromatics, {len(nit)} nitrated -> {OUT}/sites_{g}.csv ===")
        print(f"  {'site':7} {'nit':>4} {'rsa':>6} {'flk_acidic':>11} {'flk_cysmet':>11}")
        for s in sites:
            ac = f"{s['flk_acidic']:.2f}" if s["flk_acidic"] is not None else "NA"
            cy = f"{s['flk_cysmet']:.2f}" if s["flk_cysmet"] is not None else "NA"
            rsa = s["rsa"] if s["rsa"] not in ("", "None") else "NA"
            print(f"  {s['residue']}{s['resSeq']:<5} {str(s['nitrated']):>4} {rsa:>6} {ac:>11} {cy:>11}")


if __name__ == "__main__":
    main()
