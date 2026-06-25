#!/usr/bin/env python3
"""
build_system.py - robust automated build. Strips everything except standard-amino-acid
HEAVY atoms (removes glycans/ligands/ions/buffer + all H -> tleap rebuilds H), renumbers
residues sequentially (so disulfide bonds target the right atoms in multi-chain files),
auto-detects disulfides, neutralizes correctly (+ or - protein), runs tleap.

  cd lyz    && python ../scripts/build_system.py 1lz1.pdb          # all chains
  cd chi3l1 && python ../scripts/build_system.py 1hjx.pdb A        # keep only chain A
Produces sys.prmtop / sys.inpcrd / sys_built.pdb.
"""
import sys, subprocess, math

raw = sys.argv[1]
keep_chain = sys.argv[2] if len(sys.argv) > 2 else None
SALT = 20
AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HIS", "HID",
      "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
      "TYR", "VAL"}

subprocess.run(["pdb4amber", "-i", raw, "-o", "clean0.pdb", "--reduce", "--dry"], check=False)

# keep protein heavy atoms (+ optional chain), in order
records, cur = [], None
for line in open("clean0.pdb"):
    if line[:6] not in ("ATOM  ", "HETATM"):
        continue
    if line[17:20].strip() not in AA:
        continue
    chain = line[21]
    if keep_chain and chain != keep_chain:
        continue
    elem = line[76:78].strip()
    nm = line[12:16].strip().lstrip("0123456789")
    if elem == "H" or (not elem and nm[:1] == "H"):      # drop all hydrogens
        continue
    key = (chain, line[22:27])
    if cur is None or key != cur[0]:
        cur = [key, chain, []]
        records.append(cur)
    cur[2].append(line)

# renumber sequentially, TER between chains
out, n, prev = [], 0, None
for key, chain, atoms in records:
    n += 1
    if prev is not None and chain != prev:
        out.append("TER\n")
    prev = chain
    for a in atoms:
        out.append(a[:22] + f"{n:4d}" + " " + a[27:])
out.append("TER\nEND\n")

# detect disulfides on the renumbered residues
sg = []
for line in out:
    if line[:6] == "ATOM  " and line[12:16].strip() == "SG":
        sg.append((int(line[22:26]), float(line[30:38]), float(line[38:46]), float(line[46:54])))
bonds, used = [], set()
for i in range(len(sg)):
    for j in range(i + 1, len(sg)):
        if sg[i][0] in used or sg[j][0] in used:
            continue
        if math.dist(sg[i][1:], sg[j][1:]) < 2.5:
            bonds.append((sg[i][0], sg[j][0])); used.update((sg[i][0], sg[j][0]))

# rename bonded CYS -> CYX
final = []
for line in out:
    if line[:6] == "ATOM  " and line[22:26].strip():
        if int(line[22:26]) in used and line[17:20] == "CYS":
            line = line[:17] + "CYX" + line[20:]
    final.append(line)
open("clean.pdb", "w").writelines(final)
print(f"residues kept: {n} | chains: {len(set(c for _, c, _ in records))} | disulfides: {bonds}")

with open("build.in", "w") as f:
    f.write("source leaprc.protein.ff19SB\nsource leaprc.water.opc\n")
    f.write("mol = loadpdb clean.pdb\n")
    for a, b in bonds:
        f.write(f"bond mol.{a}.SG mol.{b}.SG\n")
    f.write("solvateBox mol OPCBOX 12.0\n")
    f.write("addIonsRand mol Na+ 0\naddIonsRand mol Cl- 0\n")
    f.write(f"addIonsRand mol Na+ {SALT} Cl- {SALT}\n")
    f.write("charge mol\nsaveamberparm mol sys.prmtop sys.inpcrd\nsavepdb mol sys_built.pdb\nquit\n")
subprocess.run(["tleap", "-f", "build.in"], check=False)
print("\n--> want 'Errors = 0' and net charge ~0 (Total perturbed charge: 0.000000).")
