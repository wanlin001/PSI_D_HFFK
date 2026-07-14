#!/usr/bin/env python3
"""
gen_dummy_obs.py — 產生 PSI_D 所需的 DUMMY observation 檔案

paz_rad 必須為 0.0（SKS 偏振沿徑向 Q 方向）。
過去使用 π/4 是 bug，導致 HFFK SI = −dt·cos(2(φ−BAZ))，
與 Ray 理論 SI = +dt·sin(2(φ−BAZ)) 完全相反。

用法（在 PSI_DIR 或 project 目錄執行）：
    python3 gen_dummy_obs.py

會在 PSI_DIR/psi_input/ 下產生：
    DUMMY_SI_uniform48.dat    (48 src × N_rcv)
    DUMMY_SI_uniform96.dat    (96 src × N_rcv)
    DUMMY_SP_uniform48.dat    (48 src × N_rcv)
    DUMMY_SplittingIntensity_ShearWave.dat  (Kuo2018 real stations, 24 src)

Usage on wl2:
    cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
    python3 scripts/gen_dummy_obs.py
"""

import os
from pathlib import Path

PSI_DIR = Path(__file__).parent.parent
INDIR   = PSI_DIR / "psi_input"

PAZ_RAD = 0.0   # CORRECT: SKS polarization along Q (radial) direction

# ── Reader helpers ──────────────────────────────────────────────────────────

def read_ids(path):
    """Return list of (id_str,) from first column of a CSV file."""
    ids = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        ids.append(line.split(",")[0].strip())
    return ids


def read_source_count(path):
    """Count non-comment, non-empty lines in Sources.dat."""
    n = 0
    for line in open(path):
        line = line.strip()
        if line and not line.startswith("#"):
            n += 1
    return n


# ── Writer helpers ──────────────────────────────────────────────────────────

def write_dummy_si(out_path, src_ids, rcv_ids, phase="SKS"):
    """
    SplittingIntensity format:
      SI, error, period_s, phase, src_id, rcv_id, ???, paz_rad
    """
    n = 0
    with open(out_path, "w") as f:
        for sid in src_ids:
            for rid in rcv_ids:
                f.write(f"0.0, 0.1, 0.0, {phase}, {sid}, {rid}, ???, {PAZ_RAD}\n")
                n += 1
    print(f"  Wrote {n} lines → {out_path}")


def write_dummy_sp(out_path, src_ids, rcv_ids, phase="SKS"):
    """
    SplittingParameters format:
      dt, phi, err_dt, err_phi, period_s, phase, src_id, rcv_id, quality, paz_rad
    """
    n = 0
    with open(out_path, "w") as f:
        for sid in src_ids:
            for rid in rcv_ids:
                f.write(f"0.0, 0.0, 0.1, 0.1, 0.0, {phase}, {sid}, {rid}, ???, {PAZ_RAD}\n")
                n += 1
    print(f"  Wrote {n} lines → {out_path}")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print(f"PSI_DIR  : {PSI_DIR}")
    print(f"psi_input: {INDIR}")
    print(f"paz_rad  : {PAZ_RAD}  (0.0 = correct for SKS)\n")

    # ── Receivers: same for all uniform* runs ────────────────────────────────
    rcv_file = INDIR / "Receivers.dat"
    if not rcv_file.exists():
        print(f"ERROR: {rcv_file} not found — copy Receivers.dat here first")
        return
    rcv_ids = read_ids(rcv_file)
    print(f"Receivers: {len(rcv_ids)} stations from {rcv_file}")

    # ── uniform48 (16 dir × 3 depth = 48 sources) ───────────────────────────
    src48 = INDIR / "Sources_uniform48.dat"
    if src48.exists():
        sid48 = list(range(1, read_source_count(src48) + 1))
        print(f"\nuniform48: {len(sid48)} sources")
        write_dummy_si(INDIR / "DUMMY_SI_uniform48.dat",  sid48, rcv_ids, phase="SKS")
        write_dummy_sp(INDIR / "DUMMY_SP_uniform48.dat",  sid48, rcv_ids, phase="SKS")
    else:
        print(f"\n[SKIP] {src48} not found")

    # ── uniform96 (32 dir × 3 depth = 96 sources) ───────────────────────────
    src96 = INDIR / "Sources_uniform96.dat"
    if src96.exists():
        sid96 = list(range(1, read_source_count(src96) + 1))
        print(f"\nuniform96: {len(sid96)} sources")
        write_dummy_si(INDIR / "DUMMY_SI_uniform96.dat",  sid96, rcv_ids, phase="SKS")
    else:
        print(f"\n[SKIP] {src96} not found")

    # ── Kuo2018 real stations (DUMMY_SplittingIntensity_ShearWave.dat) ───────
    kuo_src = INDIR / "Sources.dat"
    if kuo_src.exists():
        kuo_sids = list(range(1, read_source_count(kuo_src) + 1))
        print(f"\nKuo2018 real: {len(kuo_sids)} sources")
        write_dummy_si(INDIR / "DUMMY_SplittingIntensity_ShearWave.dat",
                       kuo_sids, rcv_ids, phase="SKS")
    else:
        print(f"\n[SKIP] {kuo_src} not found")

    print("\n✓ All DUMMY files generated with paz_rad = 0.0")
    print("  → Rerun all HFFK and Ray SP jobs on wl2 to get correct SI.")


if __name__ == "__main__":
    main()
