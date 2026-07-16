#!/usr/bin/env python3
"""
plot_nazimuth_convergence.py — HFFK n_azimuth 收斂圖 (T2-3)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_convergence_common import (
    convergence_stats,
    plot_rms_convergence_axes,
)

NAZ_LIST = [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64]
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

    # 左：log + monotonic envelope
    from plot_convergence_common import rms_for_log_plot, mark_reference, RMS_FLOOR
    rms_log = rms_for_log_plot(naz_list, rms_vals, args.ref_az)
    mono_log = [m if not np.isinf(m) and m > 0 else RMS_FLOOR for m in mono]
    ax = axes[0]
    ax.semilogy(naz_list, rms_log, "o-", color="tab:blue", lw=2, ms=8, label="RMS", zorder=3)
    ax.semilogy(naz_list, mono_log, "s--", color="tab:green", lw=1.5, ms=6,
                label="best-so-far", zorder=2)
    mark_reference(ax, args.ref_az, args.ref_az, for_log=True)
    ax.axhline(0.01, color="red", ls=":", lw=1.2, label="pass 0.01 s")
    ax.axvline(8, color="tab:green", ls="--", lw=1.2, label="preset=8")
    ax.set_xlabel("n_azimuth")
    ax.set_ylabel("RMS vs ref (s)")
    ax.set_xticks(naz_list)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_title("Log y + monotonic envelope", fontsize=10)

    # 右：linear
    ax2 = axes[1]
    valid = [(n, r) for n, r in zip(naz_list, rms_vals) if not np.isnan(r)]
    xs = [n for n, _ in valid]
    ys = [r for _, r in valid]
    ax2.plot(xs, ys, "o-", color="tab:blue", lw=2, ms=8)
    mark_reference(ax2, args.ref_az, args.ref_az, for_log=False)
    ax2.axhline(0.01, color="red", ls=":", lw=1.2)
    ax2.axvline(8, color="tab:green", ls="--", lw=1.2, label="preset=8")
    for x, y in zip(xs, ys):
        ax2.annotate(f"{y:.4f}", (x, y), textcoords="offset points",
                     xytext=(0, 6), ha="center", fontsize=6)
    ax2.set_xlabel("n_azimuth")
    ax2.set_ylabel("RMS vs ref (s)")
    ax2.set_xticks(xs)
    ax2.legend(fontsize=7)
    ax2.grid(True, alpha=0.3)
    ax2.set_title("Linear y — ref ★ at y=0", fontsize=10)

    fig.savefig(str(out) + ".png", dpi=150)
    print(f"Saved: {out}.png")

    print("\n=== n_azimuth convergence (raw SI) ===")
    print(f"{'n_az':>8} {'RMS(s)':>10} {'mono':>10} {'corr':>8}")
    for naz, r, m, c in zip(naz_list, rms_vals, mono, corrs):
        tag = "  ← ref" if naz == args.ref_az else ""
        print(f"{naz:>8} {r:>10.6f} {m:>10.6f} {c:>8.5f}{tag}")


if __name__ == "__main__":
    main()
