#!/usr/bin/env python3
"""
plot_benchmark_tier2.py — Tier 2 方法學收斂圖（2×3）

用法：
    python3 validation/plot_benchmark_tier2.py \
        --bench-out /lfs/wl/bench_psi/bench_output \
        --nrings-dir /lfs/wl/bench_psi/bench_output/nrings_conv \
        --nazimuth-dir /lfs/wl/bench_psi/bench_output/nazimuth_conv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_convergence_common import convergence_stats, rms_for_log_plot, mark_reference

PERIODS = [4., 8., 16., 20., 25., 33., 50.]
NRINGS_LIST = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24]
NAZ_LIST = [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64]
REF_RINGS = 16
REF_AZ = 48
PSI_DIR_DEFAULT = "/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK"
RCV_CENTER_LON, RCV_CENTER_LAT = 123.0, 24.0


def _baz(src_lon, src_lat, rcv_lon=RCV_CENTER_LON, rcv_lat=RCV_CENTER_LAT):
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


def corr_and_rms(baz_ref, si_ref, baz_h, si_h):
    if baz_ref is None or baz_h is None or len(baz_ref) == 0:
        return np.nan, np.nan
    si_ref_at_h = np.interp(baz_h, baz_ref, si_ref)
    valid = ~(np.isnan(si_h) | np.isnan(si_ref_at_h))
    if valid.sum() < 3:
        return np.nan, np.nan
    d = si_h[valid] - si_ref_at_h[valid]
    r = np.corrcoef(si_h[valid], si_ref_at_h[valid])[0, 1]
    return r, float(np.sqrt(np.mean(d ** 2)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi-dir", default=PSI_DIR_DEFAULT)
    parser.add_argument("--bench-out", default=None)
    parser.add_argument("--nrings-dir", default=None)
    parser.add_argument("--nazimuth-dir", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    psi = Path(args.psi_dir)
    bench_out = Path(args.bench_out) if args.bench_out else psi / "validation" / "bench_output"
    nrings_dir = Path(args.nrings_dir) if args.nrings_dir else bench_out / "nrings_conv"
    naz_dir = Path(args.nazimuth_dir) if args.nazimuth_dir else bench_out / "nazimuth_conv"
    out_p = Path(args.out) if args.out else bench_out / "benchmark_tier2"
    out_p.parent.mkdir(parents=True, exist_ok=True)

    bench_inp = psi / "validation" / "bench_psi_input"
    src_baz = read_sources_baz(bench_inp / "Sources.dat")
    si_src_ids = read_dummy_si_srcids(bench_inp / "DUMMY_SI.dat")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    fig.suptitle("Tier 2 — Methodological convergence", fontsize=13, fontweight="bold")

    models_r = {
        "1L_lat_B":  ("bench_1L_lat_B",  "tab:blue",   "o-"),
        "2L_lat_B":  ("bench_2L_lat_B",  "tab:orange", "s-"),
        "lateral_B": ("bench_lateral_B", "tab:red",    "^-"),
        "lateral_C": ("bench_lateral_C", "tab:purple", "v--"),
        "lateral_D": ("bench_lateral_D", "tab:green",  "d-"),
    }

    ax_r = axes[0, 0]
    ax_r.set_title("T2-1 — corr(HFFK, Ray SI) vs Period", fontsize=10)
    for label, (mod, color, style) in models_r.items():
        bR, sR = load_ray_si(bench_out, mod, si_src_ids, src_baz)
        if bR is None:
            continue
        r_vals = [corr_and_rms(bR, sR, *load_hffk(bench_out, mod, T, si_src_ids, src_baz))[0]
                  for T in PERIODS]
        ax_r.plot(PERIODS, r_vals, style, color=color, lw=2, ms=6, label=label)
    ax_r.set_xscale("log"); ax_r.set_xlim(3, 60); ax_r.set_ylim(0.98, 1.005)
    ax_r.set_xticks(PERIODS); ax_r.set_xticklabels([f"{int(T)}s" for T in PERIODS])
    ax_r.set_xlabel("Period T (s)"); ax_r.set_ylabel("corr")
    ax_r.legend(fontsize=7); ax_r.grid(True, alpha=0.3)

    ax_rms = axes[0, 1]
    ax_rms.set_title("T2-1 — RMS(HFFK − Ray SI) vs Period", fontsize=10)
    for label, (mod, color, style) in models_r.items():
        bR, sR = load_ray_si(bench_out, mod, si_src_ids, src_baz)
        if bR is None:
            continue
        rms_vals = [corr_and_rms(bR, sR, *load_hffk(bench_out, mod, T, si_src_ids, src_baz))[1]
                    for T in PERIODS]
        ax_rms.plot(PERIODS, rms_vals, style, color=color, lw=2, ms=6, label=label)
    ax_rms.set_xscale("log"); ax_rms.set_xlim(3, 60)
    ax_rms.set_xticks(PERIODS); ax_rms.set_xticklabels([f"{int(T)}s" for T in PERIODS])
    ax_rms.set_xlabel("Period T (s)"); ax_rms.set_ylabel("RMS (s)")
    ax_rms.legend(fontsize=7); ax_rms.grid(True, alpha=0.3)

    ax_grid = axes[0, 2]
    ax_grid.set_title("T2-4 — Grid resolution (HFFK T=8s)", fontsize=10)
    pairs = [
        ("1L_A vs 1L_B\nuniform", "bench_1L_A", "bench_1L_B", "tab:blue"),
        ("B vs C\n9×9 vs 17×17", "bench_lateral_B", "bench_lateral_C", "tab:orange"),
        ("C vs D\n17×17 vs 25×25", "bench_lateral_C", "bench_lateral_D", "tab:red"),
        ("B vs D\n9×9 vs 25×25", "bench_lateral_B", "bench_lateral_D", "tab:green"),
    ]
    labels, rms_g, colors = [], [], []
    for lbl, mA, mB, c in pairs:
        bA, sA = load_hffk(bench_out, mA, 8., si_src_ids, src_baz)
        bB, sB = load_hffk(bench_out, mB, 8., si_src_ids, src_baz)
        _, rms = corr_and_rms(bA, sA, bB, sB)
        if np.isnan(rms):
            continue
        labels.append(lbl); rms_g.append(rms); colors.append(c)
    bars = ax_grid.bar(range(len(labels)), rms_g, color=colors, alpha=0.8)
    ax_grid.set_xticks(range(len(labels)))
    ax_grid.set_xticklabels(labels, fontsize=7)
    ax_grid.set_ylabel("RMS (s)")
    for bar, val in zip(bars, rms_g):
        ax_grid.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                     f"{val:.2e}" if val < 0.01 else f"{val:.3f}",
                     ha="center", va="bottom", fontsize=7)
    ax_grid.set_yscale("symlog", linthresh=1e-10)

    ax_nr = axes[1, 0]
    ax_nr.set_title(f"T2-2 — n_rings (ref={REF_RINGS} ★)", fontsize=10)
    nr_result = convergence_stats(nrings_dir, "lateral_B_hffk_T25s_rings{}", NRINGS_LIST, REF_RINGS)
    if nr_result:
        _, rms_nr, _, _ = nr_result
        rms_log = rms_for_log_plot(NRINGS_LIST, rms_nr, REF_RINGS)
        ax_nr.semilogy(NRINGS_LIST, rms_log, "ko-", lw=2, ms=6)
        mark_reference(ax_nr, REF_RINGS, REF_RINGS, for_log=True)
        ax_nr.axvline(3, color="tab:green", ls="--", lw=1, label="preset=3")
        ax_nr.axhline(0.01, color="green", ls="--", lw=1, label="0.01 s")
        ax_nr.set_xlabel("n_rings"); ax_nr.set_ylabel("RMS vs ref (s)")
        ax_nr.set_xticks(NRINGS_LIST); ax_nr.legend(fontsize=7); ax_nr.grid(True, alpha=0.3)
    else:
        ax_nr.text(0.5, 0.5, "n_rings data pending", ha="center", va="center",
                   transform=ax_nr.transAxes)

    ax_naz = axes[1, 1]
    ax_naz.set_title(f"T2-3 — n_azimuth (ref={REF_AZ} ★)", fontsize=10)
    az_result = convergence_stats(naz_dir, "lateral_B_hffk_T25s_az{}", NAZ_LIST, REF_AZ)
    if az_result:
        _, rms_naz, _, _ = az_result
        rms_log = rms_for_log_plot(NAZ_LIST, rms_naz, REF_AZ)
        ax_naz.semilogy(NAZ_LIST, rms_log, "bo-", lw=2, ms=6)
        mark_reference(ax_naz, REF_AZ, REF_AZ, for_log=True)
        ax_naz.axvline(8, color="tab:green", ls="--", lw=1, label="preset=8")
        ax_naz.axhline(0.01, color="green", ls="--", lw=1, label="0.01 s")
        ax_naz.set_xlabel("n_azimuth"); ax_naz.set_ylabel("RMS vs ref (s)")
        ax_naz.set_xticks(NAZ_LIST); ax_naz.legend(fontsize=7); ax_naz.grid(True, alpha=0.3)
    else:
        ax_naz.text(0.5, 0.5, "n_azimuth data pending", ha="center", va="center",
                    transform=ax_naz.transAxes)

    ax_ref = axes[1, 2]
    ax_ref.set_title("T2-4 — lateral grid chain (vs lateral_D)", fontsize=10)
    bD, sD = load_hffk(bench_out, "bench_lateral_D", 8., si_src_ids, src_baz)
    chain = [
        ("B (9×9)", "bench_lateral_B"),
        ("C (17×17)", "bench_lateral_C"),
    ]
    if bD is not None:
        chain_rms = []
        for lbl, mod in chain:
            bH, sH = load_hffk(bench_out, mod, 8., si_src_ids, src_baz)
            _, rms = corr_and_rms(bD, sD, bH, sH)
            chain_rms.append((lbl, rms))
        lbls = [x[0] for x in chain_rms]
        vals = [x[1] for x in chain_rms]
        bars = ax_ref.bar(range(len(lbls)), vals, color=["tab:orange", "tab:red"], alpha=0.8)
        ax_ref.set_xticks(range(len(lbls))); ax_ref.set_xticklabels(lbls, fontsize=9)
        ax_ref.set_ylabel("RMS vs lateral_D (s)")
        for bar, val in zip(bars, vals):
            ax_ref.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                        f"{val:.3f}", ha="center", va="bottom", fontsize=9)
    else:
        ax_ref.text(0.5, 0.5, "lateral_D pending", ha="center", va="center",
                    transform=ax_ref.transAxes)

    fig_path = str(out_p) + ".png"
    plt.savefig(fig_path, dpi=150)
    print(f"Saved: {fig_path}")
    plt.close()


if __name__ == "__main__":
    main()
