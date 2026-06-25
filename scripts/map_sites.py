#!/usr/bin/env python3
"""
map_sites.py - eliminate numbering errors. Aligns the BUILT system (sys_built.pdb,
which may be mature-numbered, multi-chain) to the UniProt sequence and emits a
topology-INDEXED sites file for analyze_md_theory.py. No manual signal-peptide offsets.

  python map_sites.py --built ttr/sys_built.pdb --uniprot P02766 \
      --sites ../adni... or md_sites/sites_TTR.csv --out ttr/sites_top.csv

Output columns: resindex,resname,nitrated,chain,unp_pos
  resindex = 0-based residue index in topology order (what mdtraj top.residue(i) uses).
sites_<GENE>.csv must have resSeq (UniProt position) + nitrated.
"""
import argparse, csv, urllib.request

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYX": "C",
             "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "HID": "H", "HIE": "H",
             "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
             "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def uniprot_seq(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url, timeout=60) as r:
        return "".join(l.strip() for l in r.read().decode().splitlines() if not l.startswith(">"))


def read_protein_residues(pdb):
    """ordered [(global_index, chain, resname3)] for protein residues, topology order."""
    out, seen, idx = [], None, -1
    for line in open(pdb):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname not in THREE2ONE:
            continue
        chain = line[21]
        key = (chain, line[22:27])               # resSeq + icode
        if key != seen:
            seen = key; idx += 1
            out.append([idx, chain, resname])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--built", required=True)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    seq = uniprot_seq(a.uniprot)
    nit = {}
    for r in csv.DictReader(open(a.sites)):
        nit[int(r["resSeq"])] = r["nitrated"].strip().lower() in ("1", "true", "yes")

    res = read_protein_residues(a.built)
    # per-chain offset = where the chain's sequence starts in UniProt
    chains = {}
    for gi, ch, rn in res:
        chains.setdefault(ch, []).append((gi, rn))
    rows = []
    for ch, lst in chains.items():
        cseq = "".join(THREE2ONE[rn] for _, rn in lst)
        off = seq.find(cseq[:40])
        if off < 0:
            off = seq.find(cseq[:20])
        if off < 0:
            print(f"  WARN chain {ch}: sequence not found in UniProt {a.uniprot} (mutations/gaps?)")
            continue
        for k, (gi, rn) in enumerate(lst):
            unp = off + k + 1
            if rn in ("TYR", "TRP"):
                rows.append((gi, rn, nit.get(unp, False), ch, unp))

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["resindex", "resname", "nitrated", "chain", "unp_pos"])
        for r in sorted(rows):
            w.writerow(r)
    n = sum(1 for r in rows if r[2])
    print(f"wrote {len(rows)} aromatics ({n} nitrated) -> {a.out}")


if __name__ == "__main__":
    main()
