"""
Generate ORCA single-point input files for a constrained methyl rotation barrier scan.

All atoms are held fixed except the three methyl H atoms (atoms 67, 68, 69 in
1-based indexing), which are rotated around the C52-C66 bond axis.

Usage:
    python gen_methyl_scan.py /path/to/hfc.xyz [--nsteps N] [--ntasks N] [--mem STR]

    xyz        path to the geometry XYZ file (positional argument)
    --nsteps   number of steps over 360 degrees (default 36, i.e. 10 deg steps)
    --ntasks   number of tasks (CPUs) in the SLURM batch job (default 15)
               sets %pal nprocs directly — each ORCA job gets all ntasks cores
    --mem      total memory for the SLURM batch job, e.g. 300g, 128g, 64000m
               (default 300g); divided by ntasks to give %maxcore per job

Outputs:
    scan/step_000.inp, step_001.inp, ...
    scan/run_all.sh  (sequential runner; each job uses all allocated cores)
"""

import numpy as np
import os
import re
import argparse

# ── Settings ─────────────────────────────────────────────────
CHARGE     = 0
MULT       = 0
AXIS_ATOM1 = 0     # C bonded to the methyl C (fixed end of rotation axis)
AXIS_ATOM2 = 0     # methyl C
ROTOR_H    = []
ORCA_KEYWORDS = "! B3LYP def2-TZVP def2/J RIJCOSX TightSCF"
# ─────────────────────────────────────────────────────────────────────────────

def parse_mem_to_mb(mem_str):
    """Convert a SLURM-style memory string to total MB.
    Accepts: 300g, 300G, 300gb, 128000m, 128000M, 128000mb, 4t, 4T
    """
    m = re.fullmatch(r'([0-9]+(?:\.[0-9]+)?)\s*([gmtGMT][bB]?)', mem_str.strip())
    if not m:
        raise ValueError(f"Cannot parse memory string '{mem_str}'. "
                         "Use SLURM format: 300g, 128000m, 4t, etc.")
    value, unit = float(m.group(1)), m.group(2).lower().rstrip('b')
    return int(value * {'m': 1, 'g': 1024, 't': 1024**2}[unit])

def read_xyz(fname):
    with open(fname) as f:
        lines = f.readlines()
    natoms = int(lines[0])
    atoms = []
    for line in lines[2:2+natoms]:
        parts = line.split()
        atoms.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])
    return atoms

def rotation_matrix(axis, theta):
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(theta), np.sin(theta)
    t = 1 - c
    x, y, z = axis
    return np.array([
        [t*x*x + c,   t*x*y - s*z, t*x*z + s*y],
        [t*x*y + s*z, t*y*y + c,   t*y*z - s*x],
        [t*x*z - s*y, t*y*z + s*x, t*z*z + c  ]
    ])

def rotate_atoms(atoms, pivot, axis, theta, indices_1based):
    R = rotation_matrix(axis, theta)
    new_atoms = [list(a) for a in atoms]
    for idx in indices_1based:
        i = idx - 1
        pos = np.array(atoms[i][1:]) - pivot
        new_atoms[i][1:] = (R @ pos + pivot).tolist()
    return new_atoms

def write_orca_input(fname, atoms, charge, mult, nprocs, maxcore_mb, step, angle_deg):
    with open(fname, 'w') as f:
        f.write(f"{ORCA_KEYWORDS}\n")
        f.write(f"* xyz {charge} {mult}\n")
        for a in atoms:
            f.write(f"  {a[0]:<4s}  {a[1]:16.8f}  {a[2]:16.8f}  {a[3]:16.8f}\n")
        f.write("*\n")

def main():
    parser = argparse.ArgumentParser(description="Generate constrained methyl scan ORCA inputs")
    parser.add_argument("xyz",                 help="Path to the geometry XYZ file")
    parser.add_argument("--nsteps", type=int,  default=36,    help="Steps over 360 deg (default 36)")
    args = parser.parse_args()

    total_mem_mb = parse_mem_to_mb(args.mem)
    maxcore_mb   = max(256, total_mem_mb // args.ntasks)

    angles = np.linspace(0, 360, args.nsteps, endpoint=False)
    atoms  = read_xyz(args.xyz)

    p1    = np.array(atoms[AXIS_ATOM1 - 1][1:])
    p2    = np.array(atoms[AXIS_ATOM2 - 1][1:])
    axis  = p2 - p1
    pivot = p2

    os.makedirs(f"scan{AXIS_ATOM2}", exist_ok=True)

    for step, angle_deg in enumerate(angles):
        theta   = np.deg2rad(angle_deg)
        rotated = rotate_atoms(atoms, pivot, axis, theta, ROTOR_H)
        fname   = f"scan{AXIS_ATOM2}/step_{step:03d}.inp"
        write_orca_input(fname, rotated, CHARGE, MULT,
                         args.ntasks, maxcore_mb, step, angle_deg)

    print(f"Generated {args.nsteps} input files in ./scan/")
    print(f"  XYZ file      : {args.xyz}")
    print(f"  Angular step  : {360/args.nsteps:.1f} deg")
    print(f"\nTo run: cd scan && bash run_all.sh")
    print(f"To extract energies: grep 'FINAL SINGLE POINT ENERGY' scan/step_*.out")

if __name__ == "__main__":
    main()