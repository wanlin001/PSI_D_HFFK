#!/usr/bin/env python3
"""
plot_benchmark_groups.py — PSI_D Benchmark 分組比較圖

4 組比較（2 rows × 4 cols）：
  Row 0  SI vs BAZ:
    Col A — Resolution:   1L_A vs 1L_B (Ray SI + HFFK T=8s)
    Col B — Single layer: 1L_B  analytical + Ray SI + HFFK T=4–50s
    Col C — Two layers:   2L_B  Silver-Savage(ref) + Ray SI + HFFK T=4–50s
    Col D — Period effect:lateral_B Ray SI + HFFK T=4–50s

  Row 1  HFFK(T) − Ray SI residuals vs BAZ:
    Col A — 1L_A vs 1L_B HFFK T=8s 殘差（應趨近 0）
    Col B — 1L_B  各 T 殘差（均勻場：應趨近 0）
    Col C — 2L_B  各 T 殘差（垂直分層：應趨近 0）
    Col D — lateral_B 各 T 殘差（T↑ → Fresnel 平均效應）

用法（在 PSI_DIR 執行）：
    cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
    python3 validation/plot_benchmark_groups.py --bench-out /lfs/wl/bench_psi/bench_output

需要：
    bench_output/  (job_bench.sh 跑完後，含 *_ray_si/)
    validation/bench_models/analytical_si.csv
    validation/bench_psi_input/Sources.dat
    validation/bench_psi_input/DUMMY_SI.dat
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ── 常數 ─────────────────────────────────────────────────────────
PERIODS_ALL = [4., 8., 16., 20., 25., 33., 50.]
PERIODS     = PERIODS_ALL
T_COLORS    = cm.plasma(np.linspace(0.1, 0.9, len(PERIODS)))
RAY_SI_LABEL = "Ray SI (inf. freq.)"

PSI_DIR_DEFAULT = "/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK"

RCV_CENTER_LON = 123.0
RCV_CENTER_LAT = 24.0


# ══════════════════════════════════════════════════════════════════
# 資料讀取
# ══════════════════════════════════════════════════════════════════

def _baz(src_lon, src_lat,
         rcv_lon=RCV_CENTER_LON, rcv_lat=RCV_CENTER_LAT) -> float:
    """計算 BAZ（接收站 → 震源方位角，度，0–360°）。"""
    phi_r = np.radians(rcv_lat)
    phi_s = np.radians(src_lat)
    dlam  = np.radians(src_lon - rcv_lon)
    y = np.sin(dlam) * np.cos(phi_s)
    x = (np.cos(phi_r)*np.sin(phi_s)
         - np.sin(phi_r)*np.cos(phi_s)*np.cos(dlam))
    return float(np.degrees(np.arctan2(y, x)) % 360)


def read_sources_baz(sources_dat) -> dict:
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


def read_dummy_si_srcids(dummy_si) -> list:
    ids = []
    with open(dummy_si) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            ids.append(p[4])
    return ids


def read_si_output(path) -> np.ndarray:
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
    si_mean  = np.array([np.mean(groups[b]) for b in baz_vals])
    return baz_vals, si_mean


def load_analytical(csv_path, col_name):
    data = np.genfromtxt(csv_path, delimiter=',', names=True)
    names_clean = {n.strip(): n for n in data.dtype.names}
    col = names_clean.get(col_name.strip(), col_name)
    if col not in data.dtype.names:
        return None, None
    return data['BAZ_deg'], data[col]


# ══════════════════════════════════════════════════════════════════
# 診斷
# ══════════════════════════════════════════════════════════════════

def _stats(arr):
    if arr is None or len(arr) == 0:
        return (0, 0., 0., 0., 0., 0)
    a = np.asarray(arr, float)
    nz = int(np.sum(np.abs(a) > 1e-9))
    return (len(a), float(a.min()), float(a.max()),
            float(a.mean()), float(np.abs(a).max()), nz)


def diagnose_outputs(bench_out, all_models, periods):
    print("\n" + "=" * 78)
    print("BENCHMARK 輸出診斷（畫圖前先檢查實際數值）")
    print("=" * 78)
    print(f"{'model':<18}{'kind':<14}{'n':>5}{'min':>10}{'max':>10}"
          f"{'|max|':>10}{'nonzero':>9}")
    print("-" * 78)

    global_absmax = 0.0
    for model in all_models:
        ray_si_dat = bench_out / f"{model}_ray_si" / "SYN_SplittingIntensity_ShearWave.dat"
        if ray_si_dat.exists():
            si = read_si_output(ray_si_dat)
            n, mn, mx, me, am, nz = _stats(si)
            print(f"{model:<18}{'ray SI':<14}{n:>5}{mn:>10.4f}{mx:>10.4f}"
                  f"{am:>10.4f}{nz:>9}")
            global_absmax = max(global_absmax, am)
        else:
            print(f"{model:<18}{'ray SI':<14}{'--- MISSING ---'}")

        for T in periods:
            TT = f"{int(T)}s"
            si_dat = (bench_out / f"{model}_hffk_T{TT}"
                      / "SYN_SplittingIntensity_ShearWave.dat")
            if not si_dat.exists():
                print(f"{'':<18}{'hffk '+TT:<14}{'--- MISSING ---'}")
                continue
            si = read_si_output(si_dat)
            n, mn, mx, me, am, nz = _stats(si)
            print(f"{'':<18}{'hffk '+TT:<14}{n:>5}{mn:>10.4f}{mx:>10.4f}"
                  f"{am:>10.4f}{nz:>9}")
            global_absmax = max(global_absmax, am)
        print("-" * 78)

    print(f"全體 |max| = {global_absmax:.4f} s")
    print("=" * 78 + "\n")
    return global_absmax


# ══════════════════════════════════════════════════════════════════
# 單模型 SI 讀取入口
# ══════════════════════════════════════════════════════════════════

def load_ray_si(bench_out, model, si_src_ids, src_baz):
    """讀 Ray SI（TauP + SI kernel，無限高頻極限）。"""
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


# ══════════════════════════════════════════════════════════════════
# 主程式
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi-dir", default=PSI_DIR_DEFAULT)
    parser.add_argument("--bench-out", default=None,
                        help="job_bench.sh 輸出目錄")
    parser.add_argument("--out", default=None,
                        help="輸出圖片前綴")
    args = parser.parse_args()

    psi       = Path(args.psi_dir)
    bench_out = (Path(args.bench_out) if args.bench_out
                 else psi / "validation" / "bench_output")
    bench_inp = psi / "validation" / "bench_psi_input"
    csv_path  = psi / "validation" / "bench_models" / "analytical_si.csv"
    src_file  = bench_inp / "Sources.dat"
    dummy_si  = bench_inp / "DUMMY_SI.dat"
    out_p     = Path(args.out) if args.out else bench_out / "benchmark_groups"
    out_p.parent.mkdir(parents=True, exist_ok=True)

    for f in [src_file, dummy_si]:
        if not f.exists():
            sys.exit(f"ERROR: {f} not found")

    src_baz    = read_sources_baz(src_file)
    si_src_ids = read_dummy_si_srcids(dummy_si)

    if not bench_out.exists():
        sys.exit(f"ERROR: {bench_out} not found — run job_bench.sh first")

    all_models = ["bench_1L_A", "bench_1L_B", "bench_2L_A", "bench_2L_B",
                  "bench_lateral_B", "bench_lateral_C"]
    global_absmax = diagnose_outputs(bench_out, all_models, PERIODS)
    y_si  = max(0.05, global_absmax * 1.15)

    fig, axes = plt.subplots(2, 4, figsize=(22, 10), constrained_layout=True)
    fig.suptitle(
        "PSI_D Benchmark — grouped comparison\n"
        "Col A: Resolution  |  Col B: Single layer (1L_B)  |"
        "  Col C: Two layers (2L_B)  |  Col D: Period effect (lateral_B)",
        fontsize=12, fontweight='bold'
    )

    def plot_si(ax, baz, si, label, color, lw=1.5, ls='-', marker=None, zorder=2):
        if baz is None or len(baz) == 0:
            return
        ax.plot(baz, si, color=color, lw=lw, ls=ls,
                marker=marker, ms=4, label=label, zorder=zorder)

    def plot_residual(ax, baz_ref, si_ref, baz_h, si_h, label, color):
        if baz_ref is None or baz_h is None or len(baz_ref) == 0:
            return
        si_ref_at_h = np.interp(baz_h, baz_ref, si_ref)
        ax.plot(baz_h, si_h - si_ref_at_h, color=color, lw=1.2, label=label)

    # Col A — Resolution
    axA0, axA1 = axes[0, 0], axes[1, 0]
    axA0.set_title("A — Resolution\n1L_A (coarse) vs 1L_B (fine)", fontsize=10)
    axA1.set_title("A — Residual 1L_B - 1L_A\n(HFFK T=8s, should -> 0)", fontsize=9)

    for mod, color, ls in [("bench_1L_A", "tab:blue", "--"),
                            ("bench_1L_B", "tab:orange", "-")]:
        label = mod.replace("bench_", "")
        bR, sR = load_ray_si(bench_out, mod, si_src_ids, src_baz)
        b8, s8 = load_hffk(bench_out, mod, 8., si_src_ids, src_baz)
        plot_si(axA0, bR, sR,  f"{label} {RAY_SI_LABEL}", color, ls="--")
        plot_si(axA0, b8, s8,  f"{label} HFFK T=8s", color, ls=ls)

    b8A, s8A = load_hffk(bench_out, "bench_1L_A", 8., si_src_ids, src_baz)
    b8B, s8B = load_hffk(bench_out, "bench_1L_B", 8., si_src_ids, src_baz)
    if b8A is not None and b8B is not None:
        plot_residual(axA1, b8A, s8A, b8B, s8B, "1L_B − 1L_A  HFFK T=8s", "tab:green")

    bLB, sLB = load_hffk(bench_out, "bench_lateral_B", 8., si_src_ids, src_baz)
    bLC, sLC = load_hffk(bench_out, "bench_lateral_C", 8., si_src_ids, src_baz)
    plot_si(axA0, bLB, sLB, "lateral_B HFFK T=8s", "tab:red",   ls="-.", lw=1)
    plot_si(axA0, bLC, sLC, "lateral_C HFFK T=8s", "tab:purple", ls=":", lw=1)
    if bLB is not None and bLC is not None:
        plot_residual(axA1, bLB, sLB, bLC, sLC, "lateral_C − lateral_B  HFFK T=8s", "tab:red")

    # Col B — Single layer
    axB0, axB1 = axes[0, 1], axes[1, 1]
    axB0.set_title("B — Single layer (1L_B)\nanalytical + Ray SI + HFFK T=4–50 s", fontsize=10)
    axB1.set_title("B — HFFK(T) − Ray SI  (uniform: ~0)", fontsize=9)

    baz_an, si_an = load_analytical(csv_path, "SI_1L_B")
    if baz_an is not None:
        plot_si(axB0, baz_an, si_an, "Analytical", "black", lw=2, ls="--", zorder=5)

    bRayB, sRayB = load_ray_si(bench_out, "bench_1L_B", si_src_ids, src_baz)
    plot_si(axB0, bRayB, sRayB, RAY_SI_LABEL, "red", lw=2, marker='o', zorder=4)

    for T, col in zip(PERIODS, T_COLORS):
        bH, sH = load_hffk(bench_out, "bench_1L_B", T, si_src_ids, src_baz)
        lw = 2.0 if T in [4., 50.] else 1.0
        plot_si(axB0, bH, sH, f"HFFK T={int(T)}s", col, lw=lw)
        if bRayB is not None:
            plot_residual(axB1, bRayB, sRayB, bH, sH, f"T={int(T)}s", col)

    # Col C — Two layers
    axC0, axC1 = axes[0, 2], axes[1, 2]
    axC0.set_title("C — Two layers (2L_B)\nSilver-Savage(ref) + Ray SI + HFFK T=4–50 s", fontsize=10)
    axC1.set_title("C — HFFK(T) − Ray SI  (2L vertical: ~0)", fontsize=9)

    baz_an2, si_an2 = load_analytical(csv_path, "SI_2L_B_T8s")
    if baz_an2 is not None:
        plot_si(axC0, baz_an2, si_an2, "Silver-Savage T=8s (ref)",
                "black", lw=2, ls="--", zorder=5)

    bRay2B, sRay2B = load_ray_si(bench_out, "bench_2L_B", si_src_ids, src_baz)
    plot_si(axC0, bRay2B, sRay2B, RAY_SI_LABEL, "red", lw=2, marker='o', zorder=4)

    for T, col in zip(PERIODS, T_COLORS):
        bH, sH = load_hffk(bench_out, "bench_2L_B", T, si_src_ids, src_baz)
        lw = 2.0 if T in [4., 50.] else 1.0
        plot_si(axC0, bH, sH, f"HFFK T={int(T)}s", col, lw=lw)
        if bRay2B is not None:
            plot_residual(axC1, bRay2B, sRay2B, bH, sH, f"T={int(T)}s", col)

    # Col D — Lateral boundary (MAIN period plot)
    axD0, axD1 = axes[0, 3], axes[1, 3]
    axD0.set_title(
        "D — Lateral boundary (lateral_B)\nφ_L=0°/φ_R=90°  Ray SI + HFFK T=4–50 s",
        fontsize=10
    )
    axD1.set_title(
        "D — HFFK(T) − Ray SI\n(T↑ → Fresnel averages both sides)", fontsize=9
    )

    bRayLB, sRayLB = load_ray_si(bench_out, "bench_lateral_B", si_src_ids, src_baz)
    plot_si(axD0, bRayLB, sRayLB, RAY_SI_LABEL, "red", lw=2, marker='o', zorder=4)

    for T, col in zip(PERIODS, T_COLORS):
        bH, sH = load_hffk(bench_out, "bench_lateral_B", T, si_src_ids, src_baz)
        lw = 2.0 if T in [4., 50.] else 1.0
        plot_si(axD0, bH, sH, f"HFFK T={int(T)}s", col, lw=lw)
        if bRayLB is not None:
            plot_residual(axD1, bRayLB, sRayLB, bH, sH, f"T={int(T)}s", col)

    for ax0, ax1 in zip(axes[0], axes[1]):
        ax0.set_xlim(0, 360)
        ax0.set_ylim(-y_si, y_si)
        ax0.axhline(0, color='gray', lw=0.5)
        ax0.set_ylabel("SI (s)", fontsize=9)
        ax0.legend(fontsize=7, ncol=1, loc="upper right")
        ax0.set_xlabel("BAZ (°)", fontsize=9)

        ax1.set_xlim(0, 360)
        ax1.set_ylim(-y_si, y_si)
        ax1.axhline(0, color='gray', lw=0.5, ls='--')
        ax1.set_ylabel("Δ SI (s)  (see panel title)", fontsize=9)
        ax1.legend(fontsize=7, ncol=1, loc="upper right")
        ax1.set_xlabel("BAZ (°)", fontsize=9)

    fig_path = str(out_p) + ".png"
    plt.savefig(fig_path, dpi=150)
    print(f"Saved: {fig_path}")
    plt.close()

    # 第二張圖：HFFK(T) vs Ray SI 相關係數
    fig2, ax = plt.subplots(figsize=(9, 5))
    ax.set_title("r vs Period — corr( HFFK(T), Ray SI )\n"
                 "uniform/2L: ~1 all T;  lateral: drops as T grows", fontsize=11)

    models_r = {
        "1L_B":      ("bench_1L_B",      "tab:blue",   "o-"),
        "2L_B":      ("bench_2L_B",      "tab:orange",  "s-"),
        "lateral_B": ("bench_lateral_B", "tab:red",    "^-"),
        "lateral_C": ("bench_lateral_C", "tab:purple", "v--"),
    }

    for label, (mod, color, style) in models_r.items():
        r_vals = []
        bR, sR = load_ray_si(bench_out, mod, si_src_ids, src_baz)
        if bR is None or len(bR) == 0:
            continue
        for T in PERIODS:
            bH, sH = load_hffk(bench_out, mod, T, si_src_ids, src_baz)
            if bH is None or len(bH) == 0:
                r_vals.append(np.nan)
                continue
            sR_at_H = np.interp(bH, bR, sR)
            valid = ~(np.isnan(sH) | np.isnan(sR_at_H))
            if valid.sum() > 2:
                r = np.corrcoef(sH[valid], sR_at_H[valid])[0, 1]
            else:
                r = np.nan
            r_vals.append(r)
        ax.plot(PERIODS, r_vals, style, color=color, lw=2, ms=8, label=label)

    ax.set_xlabel("Period T (s)", fontsize=11)
    ax.set_ylabel("r  ( HFFK(T) vs Ray SI )", fontsize=11)
    ax.set_xscale("log")
    ax.set_xlim(3, 60)
    ax.set_ylim(0, 1.05)
    ax.axhline(1.0, color='gray', lw=0.5, ls='--')
    ax.set_xticks(PERIODS)
    ax.set_xticklabels([f"{int(T)}s" for T in PERIODS])
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig2.tight_layout()

    fig2_path = str(out_p) + "_r_vs_period.png"
    plt.savefig(fig2_path, dpi=150)
    print(f"Saved: {fig2_path}")
    plt.close()


if __name__ == "__main__":
    main()
