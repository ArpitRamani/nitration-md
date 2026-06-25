#!/usr/bin/env python3
"""
run_openmm.py - reusable OpenMM production engine for the nitration MD project.
Reads Amber prmtop/inpcrd (from tleap), runs min -> NVT heat -> NPT equil -> NPT
production on GPU, writes DCD + periodic checkpoints (HPC-restartable). Same script
for native (Exp 1) and nitro (Exp 2) systems - the modified residue lives in the
prmtop, OpenMM does not need to know about it.

  module load cuda openmm
  python run_openmm.py --prmtop sys.prmtop --inpcrd sys.inpcrd --out PROT_rep1 \
      --ns 200 --temp 300 --seed 1

ff19SB + OPC, PME 1.0 nm, HBonds + HMR -> 4 fs. 3 replicas = run 3x with --seed 1/2/3.
"""
import argparse, sys
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmm import openmm as mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop", required=True)
    ap.add_argument("--inpcrd", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--ns", type=float, default=200.0, help="production length (ns)")
    ap.add_argument("--temp", type=float, default=300.0)
    ap.add_argument("--equil_ns", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--dt_fs", type=float, default=4.0, help="4 fs requires HMR (default on)")
    ap.add_argument("--report_ps", type=float, default=100.0, help="DCD/log interval (ps)")
    a = ap.parse_args()

    prm = app.AmberPrmtopFile(a.prmtop)
    crd = app.AmberInpcrdFile(a.inpcrd)
    system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True,
                              hydrogenMass=1.5 * unit.amu)          # HMR enables 4 fs
    system.addForce(MonteCarloBarostat(1.0 * unit.bar, a.temp * unit.kelvin, 25))

    integ = LangevinMiddleIntegrator(a.temp * unit.kelvin, 1.0 / unit.picosecond,
                                     a.dt_fs * unit.femtosecond)
    integ.setRandomNumberSeed(a.seed)
    try:
        plat = Platform.getPlatformByName("CUDA")
        props = {"Precision": "mixed", "DeterministicForces": "false"}
    except Exception:
        plat = Platform.getPlatformByName("CPU"); props = {}
        print("WARNING: CUDA not found, running on CPU (slow).", file=sys.stderr)

    sim = app.Simulation(prm.topology, system, integ, plat, props)
    sim.context.setPositions(crd.positions)
    if crd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*crd.boxVectors)

    print("minimizing...", flush=True)
    sim.minimizeEnergy(maxIterations=10000)

    # NVT heat 10 -> temp
    print("NVT heating...", flush=True)
    sim.context.setVelocitiesToTemperature(10 * unit.kelvin, a.seed)
    nsteps_heat = int(100 * unit.picosecond / (a.dt_fs * unit.femtosecond))
    for i in range(20):
        integ.setTemperature((10 + (a.temp - 10) * (i + 1) / 20) * unit.kelvin)
        sim.step(nsteps_heat // 20)
    integ.setTemperature(a.temp * unit.kelvin)

    # NPT equilibration
    print(f"NPT equilibration ({a.equil_ns} ns)...", flush=True)
    sim.step(int(a.equil_ns * unit.nanoseconds / (a.dt_fs * unit.femtosecond)))

    # production
    every = int(a.report_ps * unit.picosecond / (a.dt_fs * unit.femtosecond))
    nprod = int(a.ns * unit.nanoseconds / (a.dt_fs * unit.femtosecond))
    sim.reporters.append(app.DCDReporter(f"{a.out}.dcd", every))
    sim.reporters.append(app.StateDataReporter(f"{a.out}.log", every, step=True, time=True,
                         potentialEnergy=True, temperature=True, density=True,
                         progress=True, remainingTime=True, speed=True, totalSteps=nprod))
    sim.reporters.append(app.CheckpointReporter(f"{a.out}.chk", every * 10))
    print(f"production {a.ns} ns ({nprod} steps @ {a.dt_fs} fs)...", flush=True)
    sim.step(nprod)
    sim.saveState(f"{a.out}_final.xml")
    print("done.", flush=True)


if __name__ == "__main__":
    main()
