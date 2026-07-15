#!/usr/bin/env python3
"""
gen_dummy_obs.py — 產生 PSI_D 所需的 DUMMY observation 檔案（paz_rad = 0.0）

bug 說明：
    paz_rad = π/4（0.7853981633974483）在 DUMMY 中是錯誤的。
    PSI_D 用 paz 計算 ζ = φ - BAZ - paz → SI = dt·sin(2ζ)。
    paz=π/4 → SI = -dt·cos(2(φ-BAZ))；paz=0 → SI = dt·sin(2(φ-BAZ)) ✓（Chevrot 2000）

用法（兩種模式）：

  模式 A — 修 PSI_DIR/psi_input/（uniform48/96/Kuo2018）：
      cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
      python3 scripts/gen_dummy_obs.py

  模式 B — 修某個 project 的 psi_input/（SinkingSlab 等）：
      cd /home/wl/work/ASPECT/0701_SinkingSlab_period
      python3 /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/scripts/gen_dummy_obs.py --project .

  模式 C — 直接把現有 DUMMY 檔的 paz 欄位從 π/4 改成 0（快速補丁，不重新生成）：
      python3 gen_dummy_obs.py --patch /path/to/psi_input/DUMMY_*.dat
"""

import sys
import argparse
from pathlib import Path
import math

PAZ_RAD   = 0.0                        # CORRECT for SKS (SV polarization along Q/radial)
PAZ_OLD   = "0.7853981633974483"       # buggy value that was in all DUMMY files
PAZ_NEW   = "0.0"

# ── Helpers ─────────────────────────────────────────────────────────────────

def read_ids(path):
    """Return list of ID strings from first column of a PSI CSV file."""
    ids = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        ids.append(line.split(",")[0].strip())
    return ids


def count_lines(path):
    n = 0
    for line in open(path):
        line = line.strip()
        if line and not line.startswith("#"):
            n += 1
    return n


def write_dummy_si(out_path, src_ids, rcv_ids, phase="SKS"):
    """SI, error, period_s, phase, src_id, rcv_id, ???, paz_rad"""
    n = 0
    with open(out_path, "w") as f:
        for sid in src_ids:
            for rid in rcv_ids:
                f.write(f"0.0, 0.1, 0.0, {phase}, {sid}, {rid}, ???, {PAZ_RAD}\n")
                n += 1
    print(f"  [{n:7d} lines] → {out_path}")


def write_dummy_sp(out_path, src_ids, rcv_ids, phase="SKS", period_s=50.0):
    """dt, phi, err_dt, err_phi, period_s, phase, src_id, rcv_id, quality, paz_rad
    NOTE: period_s must be > 0. PSI_D reads Phase.period from this column for SP
    observations (kernel_split_wavelet), and period<=0 throws 'Infinte frequency split!'.
    For SI observations, PSI_D uses dominant_period from the model toml instead.
    Set period_s to match dominant_period in the job toml (default 8.0 s).
    """
    n = 0
    with open(out_path, "w") as f:
        for sid in src_ids:
            for rid in rcv_ids:
                f.write(f"0.0, 0.0, 0.1, 0.1, {period_s}, {phase}, {sid}, {rid}, ???, {PAZ_RAD}\n")
                n += 1
    print(f"  [{n:7d} lines] → {out_path}  (period_s={period_s})")


def patch_paz_in_file(path):
    """Replace paz_rad = π/4 → 0.0 in-place (last CSV column)."""
    with open(path) as f:
        content = f.read()
    n = content.count(PAZ_OLD)
    if n == 0:
        print(f"  [SKIP] no π/4 found in {path}")
        return
    content2 = content.replace(PAZ_OLD, PAZ_NEW)
    with open(path, "w") as f:
        f.write(content2)
    print(f"  [patched {n:6d} lines] {path}")


# ── Mode A: PSI_DIR/psi_input (uniform48/96 + Kuo2018 real stations) ────────

def gen_psidir_dummies(psi_dir, sp_period=8.0):
    indir = psi_dir / "psi_input"
    print(f"\n=== Mode A: PSI_DIR/psi_input = {indir} ===")

    rcv_file = indir / "Receivers.dat"
    if not rcv_file.exists():
        print(f"  ERROR: {rcv_file} not found"); return
    rcv_ids = read_ids(rcv_file)
    print(f"  Receivers: {len(rcv_ids)}")

    # uniform48
    src48 = indir / "Sources_uniform48.dat"
    if src48.exists():
        sid48 = list(range(1, count_lines(src48) + 1))
        print(f"  uniform48: {len(sid48)} sources")
        write_dummy_si(indir / "DUMMY_SI_uniform48.dat", sid48, rcv_ids, phase="SKS")
        write_dummy_sp(indir / "DUMMY_SP_uniform48.dat", sid48, rcv_ids, phase="SKS", period_s=sp_period)
    else:
        print(f"  [SKIP] {src48} not found")

    # uniform96
    src96 = indir / "Sources_uniform96.dat"
    if src96.exists():
        sid96 = list(range(1, count_lines(src96) + 1))
        print(f"  uniform96: {len(sid96)} sources")
        write_dummy_si(indir / "DUMMY_SI_uniform96.dat", sid96, rcv_ids, phase="SKS")
        write_dummy_sp(indir / "DUMMY_SP_uniform96.dat", sid96, rcv_ids, phase="SKS", period_s=sp_period)
    else:
        print(f"  [SKIP] {src96} not found")

    # Kuo2018 / main real-station DUMMY
    src_kuo = indir / "Sources.dat"
    if src_kuo.exists():
        sids = list(range(1, count_lines(src_kuo) + 1))
        print(f"  Kuo2018 real: {len(sids)} sources")
        write_dummy_si(indir / "DUMMY_SplittingIntensity_ShearWave.dat",
                       sids, rcv_ids, phase="SKS")
    else:
        print(f"  [SKIP] {src_kuo} not found")


# ── Mode B: project/psi_input (SinkingSlab, any synthetic project) ───────────

def gen_project_dummies(project_dir, sp_period=8.0):
    indir = project_dir / "psi_input"
    print(f"\n=== Mode B: project/psi_input = {indir} ===")

    src_file = indir / "Sources.dat"
    rcv_file = indir / "Receivers.dat"
    if not src_file.exists() or not rcv_file.exists():
        print(f"  ERROR: Sources.dat or Receivers.dat not found in {indir}"); return

    src_ids = list(range(1, count_lines(src_file) + 1))
    rcv_ids = read_ids(rcv_file)
    print(f"  Sources: {len(src_ids)}, Receivers: {len(rcv_ids)}")

    # Determine phase from existing DUMMY if present, default S
    phase = "S"
    existing = indir / "DUMMY_SplittingIntensity_ShearWave.dat"
    if existing.exists():
        first = open(existing).readline().strip()
        if first:
            cols = [c.strip() for c in first.split(",")]
            if len(cols) >= 4:
                phase = cols[3]
    print(f"  Phase: {phase}")

    write_dummy_si(indir / "DUMMY_SplittingIntensity_ShearWave.dat",
                   src_ids, rcv_ids, phase=phase)
    write_dummy_sp(indir / "DUMMY_SplittingParameters_ShearWave.dat",
                   src_ids, rcv_ids, phase=phase, period_s=sp_period)


# ── Mode C: patch existing files in-place ────────────────────────────────────

def patch_files(paths):
    print(f"\n=== Mode C: patching {len(paths)} file(s) ===")
    for p in paths:
        patch_paz_in_file(Path(p))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Generate/fix PSI_D DUMMY obs files with paz=0")
    parser.add_argument("--project", metavar="DIR",
                        help="Mode B: project directory (has psi_input/Sources.dat)")
    parser.add_argument("--patch",   metavar="FILE", nargs="+",
                        help="Mode C: patch paz in existing DUMMY files in-place")
    parser.add_argument("--sp-period", metavar="PERIOD", type=float, default=50.0,
                        help="Period (s) for DUMMY_SP files (default: 50.0). "
                             "SP observations read Phase.period from this column for wavelet "
                             "construction; period=0 crashes PSI_D. Use large T>>δt (~50s) "
                             "to approximate infinite-frequency (ray theory) SP.")
    args = parser.parse_args()

    print(f"paz_rad will be set to: {PAZ_RAD}  (0.0 = correct for SKS/SV)")

    if args.patch:
        patch_files(args.patch)
    elif args.project:
        gen_project_dummies(Path(args.project).resolve(), sp_period=args.sp_period)
    else:
        # Mode A: PSI_DIR = parent of this script's scripts/ dir
        psi_dir = Path(__file__).parent.parent
        gen_psidir_dummies(psi_dir, sp_period=args.sp_period)

    print("\n✓ Done. Rerun all PSI_D jobs to get correct SI outputs.")
    print("  jobA (Ray SP) + jobC (HFFK SI uniform48/96) + jobD (HFFK SI Kuo2018)")


if __name__ == "__main__":
    main()
