#!/usr/bin/env python3
"""
plot_nrings_convergence.py — HFFK n_rings 收斂圖 (T2-2)

左：semilogy（跨 3+ 數量級時看清 1→16 收斂）
右：linear（0.001 s 附近細節 + ref ★ 在 y=0）

用法：
    python3 validation/plot_nrings_convergence.py \
        --conv-dir /lfs/wl/bench_psi/bench_output/nrings_conv
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_convergence_common import (
    convergence_stats,
    plot_rms_convergence_axes,
)

RINGS_LIST = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24]
PATTERN = "lateral_B_hffk_T25s_rings{}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--conv-dir", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--ref-rings", type=int, default=16)
    args = parser.parse_args()

    conv = Path(args.conv_dir)
    result = convergence_stats(conv, PATTERN, RINGS_LIST, args.ref_rings)
    if result is None:
        raise SystemExit(f"Reference not found in {conv}")

    rings_list, rms_vals, max_vals, corrs = result
    out = Path(args.out) if args.out else conv / "nrings_convergence"
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    fig.suptitle(
        "T2-2 — HFFK n_rings convergence (raw SI, 144 obs)\n"
        f"lateral_B, T=25s, n_azimuth=8, reference: n_rings={args.ref_rings}",
        fontsize=11, fontweight="bold",
    )

    plot_rms_convergence_axes(
        axes[0], axes[1], rings_list, rms_vals, args.ref_rings,
        preset_val=3, preset_label="preset",
        xlabel="n_rings",
        title_log="Log y — full dynamic range",
        title_lin="Linear y — fine structure + ref ★",
    )

    fig.savefig(str(out) + ".png", dpi=150)
    print(f"Saved: {out}.png")

    print("\n=== n_rings convergence (raw SI) ===")
    print(f"{'n_rings':>8} {'RMS(s)':>10} {'max(s)':>10} {'corr':>8}")
    for nr, r, m, c in zip(rings_list, rms_vals, max_vals, corrs):
        tag = "  ← ref" if nr == args.ref_rings else ""
        print(f"{nr:>8} {r:>10.6f} {m:>10.6f} {c:>8.5f}{tag}")


if __name__ == "__main__":
    main()
