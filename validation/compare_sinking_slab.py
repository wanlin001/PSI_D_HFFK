#!/usr/bin/env python3
"""
compare_sinking_slab.py — T3-3 SinkingSlab Ray vs HFFK 比較

比較：
  1. HFFK Ray SI vs 官方 PSI_D SinkingSlab ray 輸出（應機器精度一致）
  2. HFFK(T) vs Ray SI（period 效應）

用法：
    module load anaconda/2022.05
    python3 validation/compare_sinking_slab.py \
        --sink-out /lfs/wl/bench_psi/bench_output/sinking_slab
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OFFICIAL_RAY = Path(
    "/home/wl/software/ECOMAN2.0-seismology.PSI_D/examples/SinkingSlab"
    "/psi_output/SYN_SinkingBlock/SYN_SplittingIntensity_ShearWave.dat"
)
PERIODS = [3., 5., 8., 15., 25.]


def read_si(path):
    vals = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                vals.append(float(ln.split(",")[0]))
    return np.array(vals)


def stats(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    d = b - a
    r = np.corrcoef(a, b)[0, 1] if n > 2 else np.nan
    return {
        "n": n,
        "rms": float(np.sqrt(np.mean(d ** 2))),
        "max": float(np.max(np.abs(d))),
        "corr": float(r),
        "mean_a": float(a.mean()),
        "mean_b": float(b.mean()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sink-out", required=True)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    sink = Path(args.sink_out)
    out_p = Path(args.out) if args.out else sink / "sinking_slab_comparison"
    out_p.parent.mkdir(parents=True, exist_ok=True)

    ray_si_path = sink / "SinkingSlab_ray_si" / "SYN_SplittingIntensity_ShearWave.dat"
    if not ray_si_path.exists():
        sys.exit(f"ERROR: {ray_si_path} not found — run job_sinking_slab.sh first")

    ray_si = read_si(ray_si_path)

    print("\n" + "=" * 70)
    print("T3-3 — SinkingSlab validation")
    print("=" * 70)

    if OFFICIAL_RAY.exists():
        official = read_si(OFFICIAL_RAY)
        s = stats(official, ray_si)
        print("\n[1a] Ray SI (HFFK repo) vs Official PSI_D SinkingSlab")
        print(f"     n={s['n']}  RMS={s['rms']:.3e}  max|diff|={s['max']:.3e}  corr={s['corr']:.6f}")
        if s["rms"] < 1e-10:
            print("     PASS — 與官方 Ray 結果一致（機器精度）")
        elif s["corr"] > 0.9999:
            print("     PASS — 高度一致（corr>0.9999）")
        else:
            print("     WARN — HFFK fork 的 TauP Ray SI 與原版 PSI_D 有系統差異")
            print("            → T3-3 應以官方 Ray 為基準比較 HFFK（見 [1b]）")
    else:
        print(f"\n[WARN] Official reference not found: {OFFICIAL_RAY}")
        official = None

    if official is not None:
        print("\n[1b] HFFK(T) vs Official PSI_D Ray（T3-3 主比較）")
        print(f"{'Period':>8} {'corr':>10} {'RMS(s)':>12} {'max|diff|':>12}")
        print("-" * 46)
        for T in PERIODS:
            p = sink / f"SinkingSlab_hffk_T{int(T)}s" / "SYN_SplittingIntensity_ShearWave.dat"
            if not p.exists():
                print(f"{int(T):>6}s  --- MISSING ---")
                continue
            hffk = read_si(p)
            s = stats(official, hffk)
            print(f"{int(T):>6}s  {s['corr']:>10.5f} {s['rms']:>12.5f} {s['max']:>12.5f}")

    print("\n[2] HFFK(T) vs Ray SI（HFFK repo 內部一致性）")
    print(f"{'Period':>8} {'corr':>10} {'RMS(s)':>12} {'max|diff|':>12}")
    print("-" * 46)

    corrs, rms_vals = [], []
    for T in PERIODS:
        p = sink / f"SinkingSlab_hffk_T{int(T)}s" / "SYN_SplittingIntensity_ShearWave.dat"
        if not p.exists():
            print(f"{int(T):>6}s  --- MISSING ---")
            corrs.append(np.nan)
            rms_vals.append(np.nan)
            continue
        hffk = read_si(p)
        s = stats(ray_si, hffk)
        print(f"{int(T):>6}s  {s['corr']:>10.5f} {s['rms']:>12.5f} {s['max']:>12.5f}")
        corrs.append(s["corr"])
        rms_vals.append(s["rms"])

    # 圖
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    fig.suptitle("T3-3 — SinkingSlab: Ray SI vs HFFK", fontsize=12, fontweight="bold")

    ax = axes[0]
    valid = [not np.isnan(c) for c in corrs]
    T_valid = [T for T, v in zip(PERIODS, valid) if v]
    c_valid = [c for c, v in zip(corrs, valid) if v]
    ax.plot(T_valid, c_valid, "o-", lw=2, ms=8, color="tab:blue")
    ax.set_xlabel("Period T (s)")
    ax.set_ylabel("corr(HFFK, Ray SI)")
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    r_valid = [r for r, v in zip(rms_vals, valid) if v]
    ax.plot(T_valid, r_valid, "s-", lw=2, ms=8, color="tab:red")
    ax.set_xlabel("Period T (s)")
    ax.set_ylabel("RMS(HFFK − Ray SI) (s)")
    ax.grid(True, alpha=0.3)

    fig.savefig(str(out_p) + ".png", dpi=150)
    print(f"\nSaved: {out_p}.png")
    plt.close()

    print("\n預期（VanderBeek 2021）：")
    print("  - Ray SI = 官方 PSI_D（應完全一致）")
    print("  - HFFK 趨勢與 Ray 相同，slab 邊界附近更平滑")
    print("  - corr 不必接近 1.0（真實 3D 結構的 Fresnel 平均效應）")


if __name__ == "__main__":
    main()
