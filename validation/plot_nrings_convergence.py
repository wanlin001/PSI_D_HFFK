#!/usr/bin/env python3
"""
plot_nrings_convergence.py — HFFK n_rings 收斂圖 (T2-2)

用法：
    module load anaconda/2022.05
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
from plot_convergence_common import convergence_stats

RINGS_LIST = [1, 2, 3, 4, 5, 6, 8, 12, 16]
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

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    fig.suptitle(
        "T2-2 — HFFK n_rings convergence (raw SI, 144 obs)\n"
        f"lateral_B, T=25s, n_azimuth=8, reference: n_rings={args.ref_rings}",
        fontsize=11, fontweight="bold",
    )

    axes[0].semilogy(rings_list, rms_vals, "o-", color="tab:red", lw=2, ms=8)
    axes[0].axhline(0.01, color="green", ls="--", lw=1)
    axes[0].set_xlabel("n_rings")
    axes[0].set_ylabel("RMS diff vs ref (s)")
    axes[0].set_xticks(rings_list)
    axes[0].grid(True, alpha=0.3)

    axes[1].semilogy(rings_list, max_vals, "s-", color="tab:orange", lw=2, ms=8)
    axes[1].set_xlabel("n_rings")
    axes[1].set_ylabel("max |diff| vs ref (s)")
    axes[1].set_xticks(rings_list)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(rings_list, corrs, "^-", color="tab:blue", lw=2, ms=8)
    axes[2].set_xlabel("n_rings")
    axes[2].set_ylabel(f"corr( SI(nr), SI({args.ref_rings}) )")
    axes[2].set_ylim(0.9, 1.001)
    axes[2].set_xticks(rings_list)
    axes[2].axhline(1.0, color="gray", ls="--", lw=0.5)
    axes[2].grid(True, alpha=0.3)

    fig.savefig(str(out) + ".png", dpi=150)
    print(f"Saved: {out}.png")

    print("\n=== n_rings convergence (raw SI) ===")
    print(f"{'n_rings':>8} {'RMS(s)':>10} {'max(s)':>10} {'corr':>8}")
    for nr, r, m, c in zip(rings_list, rms_vals, max_vals, corrs):
        print(f"{nr:>8} {r:>10.6f} {m:>10.6f} {c:>8.5f}")


if __name__ == "__main__":
    main()
