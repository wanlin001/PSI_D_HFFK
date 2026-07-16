#!/usr/bin/env python3
"""
plot_nazimuth_convergence.py — HFFK n_azimuth 收斂圖 (T2-3)

重點：RMS 主圖 + monotonic envelope，避免 6→8 非單調造成「沒收斂」錯覺。

用法：
    python3 validation/plot_nazimuth_convergence.py \
        --conv-dir /lfs/wl/bench_psi/bench_output/nazimuth_conv
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_convergence_common import convergence_stats

NAZ_LIST = [4, 6, 8, 12, 16, 24, 32, 48]
PATTERN = "lateral_B_hffk_T25s_az{}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--conv-dir", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--ref-az", type=int, default=48)
    args = parser.parse_args()

    conv = Path(args.conv_dir)
    result = convergence_stats(conv, PATTERN, NAZ_LIST, args.ref_az)
    if result is None:
        raise SystemExit(f"Reference not found in {conv}")

    naz_list, rms_vals, max_vals, corrs = result
    mono = np.minimum.accumulate(
        [v if not np.isnan(v) else np.inf for v in rms_vals]
    )
    out = Path(args.out) if args.out else conv / "nazimuth_convergence"
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    fig.suptitle(
        "T2-3 — HFFK n_azimuth convergence (raw SI, 144 obs)\n"
        f"lateral_B, T=25s, n_rings=3, reference: n_azimuth={args.ref_az}",
        fontsize=11, fontweight="bold",
    )

    # 左：全範圍 semilogy + monotonic envelope
    ax = axes[0]
    ax.semilogy(naz_list, rms_vals, "o-", color="tab:blue", lw=2, ms=9,
                label="RMS per n_azimuth", zorder=3)
    ax.semilogy(naz_list, mono, "s--", color="tab:green", lw=1.5, ms=6,
                label="monotonic best-so-far", zorder=2)
    ax.axhline(0.01, color="red", ls=":", lw=1.2, label="pass: 0.01 s")
    for x, y in zip(naz_list, rms_vals):
        if not np.isnan(y):
            ax.annotate(f"{y:.4f}", (x, y), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=7)
    ax.set_xlabel("n_azimuth")
    ax.set_ylabel("RMS diff vs ref (s)")
    ax.set_xticks(naz_list)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3, which="both")

    # 右：az>=6 線性尺度放大（看 6→48 收斂斜率）
    ax2 = axes[1]
    mask = [n >= 6 for n in naz_list]
    x6 = [n for n, m in zip(naz_list, mask) if m]
    r6 = [r for r, m in zip(rms_vals, mask) if m and not np.isnan(r)]
    m6 = [m for m, mm in zip(mono, mask) if mm and not np.isnan(m)]
    ax2.plot(x6, r6, "o-", color="tab:blue", lw=2, ms=9, label="RMS")
    ax2.plot(x6, m6, "s--", color="tab:green", lw=1.5, ms=6, label="best-so-far")
    ax2.axhline(0.01, color="red", ls=":", lw=1.2, label="0.01 s")
    ax2.axhline(0.001, color="orange", ls=":", lw=1.0, label="0.001 s")
    ax2.set_xlabel("n_azimuth (≥6)")
    ax2.set_ylabel("RMS diff vs ref (s)")
    ax2.set_xticks(x6)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_title("Zoom: az ≥ 6 (linear y)", fontsize=10)

    fig.savefig(str(out) + ".png", dpi=150)
    print(f"Saved: {out}.png")

    print("\n=== n_azimuth convergence (raw SI) ===")
    print(f"{'n_az':>8} {'RMS(s)':>10} {'mono':>10} {'corr':>8}")
    for naz, r, m, c in zip(naz_list, rms_vals, mono, corrs):
        print(f"{naz:>8} {r:>10.6f} {m:>10.6f} {c:>8.5f}")
    print(f"\nPreset n_azimuth=8: RMS={rms_vals[2]:.6f} s  (pass < 0.01 s)")


if __name__ == "__main__":
    main()
