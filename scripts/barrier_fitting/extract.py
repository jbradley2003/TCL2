"""
Extract ORCA single-point energies from a scan output directory into a CSV.

Usage:
    python extract_energies.py <scan_dir> [--out FILE] [--nsteps N]

    scan_dir   directory containing step_000.out, step_001.out, ...
    --out      output CSV filename (default: energies.csv in scan_dir)
    --nsteps   number of steps (default: auto-detect from *.out files)

Output CSV columns:
    step        step index (0-based)
    angle_deg   rotation angle in degrees
    energy_Eh   ORCA final single point energy in Hartree
    energy_rel_Eh   energy relative to the minimum (Hartree)
    energy_rel_kcal energy relative to the minimum (kcal/mol)
    energy_rel_kJ   energy relative to the minimum (kJ/mol)
    energy_rel_K    energy relative to the minimum (Kelvin, for barrier comparisons)
"""

import os
import re
import csv
import argparse
import glob

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_KJ   = 2625.49963
HARTREE_TO_K    = 315775.02  # 1 Eh = kB * T  ->  T in Kelvin

def find_energy(outfile):
    """Return the last 'FINAL SINGLE POINT ENERGY' value in an ORCA output file,
    or None if not found (e.g. job did not finish)."""
    energy = None
    pattern = re.compile(r"FINAL SINGLE POINT ENERGY\s+([-\d.]+)")
    with open(outfile) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                energy = float(m.group(1))
    return energy

def main():
    parser = argparse.ArgumentParser(description="Extract ORCA scan energies to CSV")
    parser.add_argument("scan_dir",         help="Directory containing step_NNN.out files")
    parser.add_argument("--out",  type=str, default=None,
                        help="Output CSV path (default: <scan_dir>/energies.csv)")
    parser.add_argument("--nsteps", type=int, default=None,
                        help="Number of steps (default: auto-detect)")
    args = parser.parse_args()

    scan_dir = args.scan_dir.rstrip("/")
    out_csv  = args.out or os.path.join(scan_dir, "energies.csv")

    # Find output files
    outfiles = sorted(glob.glob(os.path.join(scan_dir, "step_*.out")))
    if not outfiles:
        raise FileNotFoundError(f"No step_*.out files found in '{scan_dir}'")

    nsteps = args.nsteps or len(outfiles)
    angles = [360.0 * i / nsteps for i in range(nsteps)]

    # Parse energies
    rows = []
    missing = []
    for i, fpath in enumerate(outfiles):
        # Extract step index from filename
        m = re.search(r"step_(\d+)\.out$", fpath)
        step = int(m.group(1)) if m else i
        angle = angles[step] if step < len(angles) else float("nan")

        energy = find_energy(fpath)
        if energy is None:
            missing.append(os.path.basename(fpath))
        rows.append({"step": step, "angle_deg": angle, "energy_Eh": energy})

    if missing:
        print(f"Warning: no energy found in {len(missing)} file(s): {', '.join(missing)}")

    # Compute relative energies (skip None values)
    valid_energies = [r["energy_Eh"] for r in rows if r["energy_Eh"] is not None]
    if not valid_energies:
        raise ValueError("No energies could be parsed from any output file.")

    e_min = min(valid_energies)
    for r in rows:
        if r["energy_Eh"] is not None:
            delta = r["energy_Eh"] - e_min
            r["energy_rel_Eh"]    = delta
            r["energy_rel_kcal"]  = delta * HARTREE_TO_KCAL
            r["energy_rel_kJ"]    = delta * HARTREE_TO_KJ
            r["energy_rel_K"]     = delta * HARTREE_TO_K
        else:
            r["energy_rel_Eh"]    = None
            r["energy_rel_kcal"]  = None
            r["energy_rel_kJ"]    = None
            r["energy_rel_K"]     = None

    # Write CSV
    fieldnames = ["step", "angle_deg", "energy_Eh",
                  "energy_rel_Eh", "energy_rel_kcal", "energy_rel_kJ", "energy_rel_K"]
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_csv}")
    print(f"  Steps parsed  : {len(valid_energies)} / {len(rows)}")
    print(f"  Min energy    : {e_min:.10f} Eh")
    print(f"  Max barrier   : {(max(valid_energies) - e_min) * HARTREE_TO_K:.1f} K  "
          f"({(max(valid_energies) - e_min) * HARTREE_TO_KCAL:.3f} kcal/mol)")

if __name__ == "__main__":
    main()