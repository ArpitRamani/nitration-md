#!/usr/bin/env python3
"""
flatten_assembly.py - collapse a multi-MODEL biological-assembly PDB (e.g. an RCSB
.pdb1, like TTR 1tta.pdb1) into a SINGLE model with unique chain IDs, keeping only
altloc A/blank and protein ATOM records. Fixes the multi-model + alternate-conformation
parse error that breaks pdb4amber/parmed on biological assemblies.

  python flatten_assembly.py 1tta.pdb1 flat.pdb
  -> then: python build_system.py flat.pdb
"""
import sys

inp, out = sys.argv[1], sys.argv[2]
CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
cur_model, ci, seen, lines = 0, 0, {}, []
for line in open(inp):
    rec = line[:6]
    if rec == "MODEL ":
        cur_model += 1; continue
    if rec == "ENDMDL":
        continue
    if rec == "ATOM  ":                       # protein only; drop HETATM (waters/ligands)
        if line[16] not in (" ", "A"):        # keep one altloc
            continue
        key = (cur_model, line[21])
        if key not in seen:
            seen[key] = CHAINS[ci]; ci += 1
        line = line[:16] + " " + line[17:21] + seen[key] + line[22:]   # blank altloc, new chain
        lines.append(line)
with open(out, "w") as f:
    f.writelines(lines); f.write("END\n")
print(f"flattened {inp}: {len(seen)} chains -> {out}")
