#!/usr/bin/env python3
"""
plot_benchmark_tier1.py — Tier 1 單元正確性圖（獨立輸出）

Col A: 格網解析度（1L_A vs 1L_B，均勻場）
Col B: 單層橫向邊界（bench_1L_lat_B）— Ray SI + HFFK T=4–50s
Col C: 雙層下層橫向邊界（bench_2L_lat_B）— Ray SI + HFFK T=4–50s
Col D: 全深度橫向邊界（bench_lateral_B）— Ray SI + HFFK T=4–50s

用法：
    python3 validation/plot_benchmark_tier1.py \
        --bench-out /lfs/wl/bench_psi/bench_output
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

PERIODS = [4., 8., 16., 20., 25., 33., 50.]
T_COLORS = cm.plasma(np.linspace(0.1, 0.9, len(PERIODS)))
RAY_SI_LABEL = "Ray SI (inf. freq.)"
PSI_DIR_DEFAULT = "/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK"
RCV_CENTER_LON, RCV_CENTER_LAT = 123.0, 24.0


def _baz(src_lon, src_lat,
         rcv_lon=RCV_CENTER_LON, rcv_lat=RCV_CENTER_LAT) -> float:
    phi_r = np.radians(rcv_lat)
    phi_s = np.radians(src_lat)
    dlam = np.radians(src_lon - rcv_lon)
    y = np.sin(dlam) * np.cos(phi_s)
    x = (np.cos(phi_r) * np.sin(phi_s)
         - np.sin(phi_r) * np.cos(phi_s) * np.cos(dlam))
    return float(np.degrees(np.arctan2(y, x)) % 360)


def read_sources_baz(sources_dat):
    baz = {}
    with open(sources_dat) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            if len(p) >= 3:
                baz[p[0]] = _baz(float(p[1]), float(p[2]))
    return baz


def read_dummy_si_srcids(dummy_si):
    ids = []
    with open(dummy_si) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            ids.append(p[4])
    return ids


def read_si_output(path):
    vals = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            vals.append(float(ln.split(',')[0].strip()))
    return np.array(vals)


def si_per_baz(si_arr, src_ids, src_baz_map):
    from collections import defaultdict
    groups = defaultdict(list)
    for si, sid in zip(si_arr, src_ids):
        if sid in src_baz_map:
            groups[src_baz_map[sid]].append(si)
    if not groups:
        return np.array([]), np.array([])
    baz_vals = np.array(sorted(groups.keys()))
    si_mean = np.array([np.mean(groups[b]) for b in baz_vals])
    return baz_vals, si_mean


def load_ray_si(bench_out, model, si_src_ids, src_baz):
    dat = bench_out / f"{model}_ray_si" / "SYN_SplittingIntensity_ShearWave.dat"
    if not dat.exists():
        return None, None
    si = read_si_output(dat)
    return si_per_baz(si, si_src_ids[:len(si)], src_baz)


def load_hffk(bench_out, model, period_s, si_src_ids, src_baz):
    TT = f"{int(period_s)}s"
    dat = bench_out / f"{model}_hffk_T{TT}" / "SYN_SplittingIntensity_ShearWave.dat"
    if not dat.exists():
        return None, None
    si = read_si_output(dat)
    return si_per_baz(si, si_src_ids[:len(si)], src_baz)


def rms_diff(baz_a, si_a, baz_b, si_b):
    if baz_a is None or baz_b is None or len(baz_a) == 0:
        return np.nan
    si_a_at_b = np.interp(baz_b, baz_a, si_a)
    d = si_b - si_a_at_b
    return float(np.sqrt(np.mean(d ** 2)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi-dir", default=PSI_DIR_DEFAULT)
    parser.add_argument("--bench-out", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    psi = Path(args.psi_dir)
    bench_out = Path(args.bench_out) if args.bench_out else psi / "validation" / "bench_output"
    bench_inp = psi / "validation" / "bench_psi_input"
    out_p = Path(args.out) if args.out else bench_out / "benchmark_tier1"
    out_p.parent.mkdir(parents=True, exist_ok=True)

    src_baz = read_sources_baz(bench_inp / "Sources.dat")
    si_src_ids = read_dummy_si_srcids(bench_inp / "DUMMY_SI.dat")

    if not bench_out.exists():
        sys.exit(f"ERROR: {bench_out} not found")

    models = ["bench_1L_A", "bench_1L_B", "bench_1L_lat_B",
              "bench_2L_lat_B", "bench_lateral_B"]
    global_absmax = 0.05
    for mod in models:
        bR, sR = load_ray_si(bench_out, mod, si_src_ids, src_baz)
        if sR is not None:
            global_absmax = max(global_absmax, float(np.abs(sR).max()))
        for T in PERIODS:
            bH, sH = load_hffk(bench_out, mod, T, si_src_ids, src_baz)
            if sH is not None:
                global_absmax = max(global_absmax, float(np.abs(sH).max()))
    y_si = global_absmax * 1.15

    fig, axes = plt.subplots(2, 4, figsize=(22, 10), constrained_layout=True)
    fig.suptitle(
        "Tier 1 — Unit correctness (synthetic benchmark)\n"
        "A: Grid resolution  |  B: 1L lateral φ  |  C: 2L lower lateral φ  |  D: Full lateral φ",
        fontsize=12, fontweight="bold",
    )

    def plot_si(ax, baz, si, label, color, lw=1.5, ls="-", marker=None, zorder=2):
        if baz is None or len(baz) == 0:
            return
        ax.plot(baz, si, color=color, lw=lw, ls=ls,
                marker=marker, ms=4, label=label, zorder=zorder)

    def plot_residual(ax, baz_ref, si_ref, baz_h, si_h, label, color):
        if baz_ref is None or baz_h is None or len(baz_ref) == 0:
            return
        si_ref_at_h = np.interp(baz_h, baz_ref, si_ref)
        ax.plot(baz_h, si_h - si_ref_at_h, color=color, lw=1.2, label=label)

    # Col A — resolution only
    axA0, axA1 = axes[0, 0], axes[1, 0]
    axA0.set_title("A — Grid resolution\n1L_A (5×5) vs 1L_B (9×9), uniform φ=45°", fontsize=10)
    axA1.set_title("A — 1L_B − 1L_A  HFFK T=8s\n(uniform: → 0)", fontsize=9)

    for mod, color, ls in [("bench_1L_A", "tab:blue", "--"),
                           ("bench_1L_B", "tab:orange", "-")]:
        label = mod.replace("bench_", "")
        bR, sR = load_ray_si(bench_out, mod, si_src_ids, src_baz)
        b8, s8 = load_hffk(bench_out, mod, 8., si_src_ids, src_baz)
        plot_si(axA0, bR, sR, f"{label} {RAY_SI_LABEL}", color, ls="--")
        plot_si(axA0, b8, s8, f"{label} HFFK T=8s", color, ls=ls)

    b8A, s8A = load_hffk(bench_out, "bench_1L_A", 8., si_src_ids, src_baz)
    b8B, s8B = load_hffk(bench_out, "bench_1L_B", 8., si_src_ids, src_baz)
    if b8A is not None and b8B is not None:
        plot_residual(axA1, b8A, s8A, b8B, s8B, "1L_B − 1L_A", "tab:green")

    # Col B/C/D — period-sensitive lateral models
    panels = [
        (1, "bench_1L_lat_B", "B — Single-layer lateral\nφ_L=0°/φ_R=90° @ 123°E"),
        (2, "bench_2L_lat_B", "C — Two-layer: upper uniform φ=0°\nlower layer lateral boundary"),
        (3, "bench_lateral_B", "D — Full-depth lateral\nφ_L=0°/φ_R=90° @ 123°E"),
    ]

    for col, model, title in panels:
        ax0, ax1 = axes[0, col], axes[1, col]
        ax0.set_title(title + "\nRay SI + HFFK T=4–50 s", fontsize=10)
        ax1.set_title(f"{title.split()[0]} — HFFK(T) − Ray SI", fontsize=9)

        bRay, sRay = load_ray_si(bench_out, model, si_src_ids, src_baz)
        plot_si(ax0, bRay, sRay, RAY_SI_LABEL, "red", lw=2, marker="o", zorder=4)

        for T, col_c in zip(PERIODS, T_COLORS):
            bH, sH = load_hffk(bench_out, model, T, si_src_ids, src_baz)
            lw = 2.0 if T in (4., 50.) else 1.0
            plot_si(ax0, bH, sH, f"HFFK T={int(T)}s", col_c, lw=lw)
            if bRay is not None:
                plot_residual(ax1, bRay, sRay, bH, sH, f"T={int(T)}s", col_c)

    for ax0, ax1 in zip(axes[0], axes[1]):
        ax0.set_xlim(0, 360)
        ax0.set_ylim(-y_si, y_si)
        ax0.axhline(0, color="gray", lw=0.5)
        ax0.set_ylabel("SI (s)", fontsize=9)
        ax0.legend(fontsize=6.5, ncol=1, loc="upper right")
        ax0.set_xlabel("BAZ (°)", fontsize=9)

        ax1.set_xlim(0, 360)
        ax1.set_ylim(-y_si, y_si)
        ax1.axhline(0, color="gray", lw=0.5, ls="--")
        ax1.set_ylabel("Δ SI (s)", fontsize=9)
        ax1.legend(fontsize=6.5, ncol=1, loc="upper right")
        ax1.set_xlabel("BAZ (°)", fontsize=9)

    fig_path = str(out_p) + ".png"
    plt.savefig(fig_path, dpi=150)
    print(f"Saved: {fig_path}")
    plt.close()


if __name__ == "__main__":
    main()
