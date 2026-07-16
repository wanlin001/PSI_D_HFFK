#!/usr/bin/env python3
"""
plot_convergence_common.py — 收斂測試共用工具

供 plot_nrings_convergence.py / plot_nazimuth_convergence.py /
plot_benchmark_tier2.py 共用，確保 RMS 計算方式一致。

規則：一律用 **raw SI**（全觀測點，不經 BAZ 平均）。
"""

from pathlib import Path

import numpy as np

# ref 點 RMS=0 無法上 semilogy；用 floor 標出參考位置
RMS_FLOOR = 1e-16


def read_si(path):
    vals = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                vals.append(float(ln.split(",")[0]))
    return np.array(vals)


def convergence_stats(conv_dir, pattern_fmt, param_list, ref_val):
    """
    回傳 (param_list, rms_vals, max_vals, corrs) 相對 ref_val 的收斂統計。

    pattern_fmt: 例如 "lateral_B_hffk_T25s_rings{}" 或 "..._az{}"
    """
    ref_path = Path(conv_dir) / pattern_fmt.format(ref_val) / "SYN_SplittingIntensity_ShearWave.dat"
    if not ref_path.exists():
        return None

    ref = read_si(ref_path)
    rms_vals, max_vals, corrs = [], [], []

    for p in param_list:
        path = Path(conv_dir) / pattern_fmt.format(p) / "SYN_SplittingIntensity_ShearWave.dat"
        if not path.exists():
            rms_vals.append(np.nan)
            max_vals.append(np.nan)
            corrs.append(np.nan)
            continue
        si = read_si(path)
        d = si - ref
        rms_vals.append(0.0 if p == ref_val else float(np.sqrt(np.mean(d ** 2))))
        max_vals.append(0.0 if p == ref_val else float(np.max(np.abs(d))))
        corrs.append(1.0 if p == ref_val else float(np.corrcoef(si, ref)[0, 1]))

    return param_list, rms_vals, max_vals, corrs


def rms_for_log_plot(param_list, rms_vals, ref_val):
    """ref 點 RMS=0 → floor，供 semilogy 顯示 ★ ref。"""
    out = []
    for p, r in zip(param_list, rms_vals):
        if np.isnan(r):
            out.append(np.nan)
        elif p == ref_val or r <= 0.0:
            out.append(RMS_FLOOR)
        else:
            out.append(r)
    return out


def mark_reference(ax, ref_x, ref_label, *, for_log=False, color="gold"):
    """畫 ref 垂直線 + ★ 標記。"""
    ax.axvline(ref_x, color=color, ls=":", lw=1.5, alpha=0.85, zorder=1)
    y_ref = RMS_FLOOR if for_log else 0.0
    ax.scatter(
        [ref_x], [y_ref], marker="*", s=220, color=color,
        edgecolors="black", linewidths=0.6, zorder=6,
        label=f"★ ref={ref_label}",
    )


def plot_rms_convergence_axes(
    ax_log, ax_lin, param_list, rms_vals, ref_val, *,
    preset_val=None, preset_label="preset",
    xlabel="parameter", title_log="log y", title_lin="linear y",
):
    """左 log、右 linear 的 RMS 收斂雙面板。"""
    rms_log = rms_for_log_plot(param_list, rms_vals, ref_val)
    valid = [(p, r) for p, r in zip(param_list, rms_vals) if not np.isnan(r)]

    ax_log.semilogy(param_list, rms_log, "o-", color="tab:blue", lw=2, ms=7, zorder=3)
    mark_reference(ax_log, ref_x=ref_val, ref_label=ref_val, for_log=True)
    ax_log.axhline(0.01, color="red", ls=":", lw=1.2, label="pass 0.01 s")
    ax_log.axhline(0.001, color="orange", ls=":", lw=1.0, label="0.001 s")
    if preset_val is not None:
        ax_log.axvline(preset_val, color="tab:green", ls="--", lw=1.2,
                       label=f"{preset_label}={preset_val}")
    ax_log.set_xlabel(xlabel)
    ax_log.set_ylabel("RMS vs ref (s)")
    ax_log.set_title(title_log, fontsize=10)
    ax_log.set_xticks(param_list)
    ax_log.legend(fontsize=7, loc="upper right")
    ax_log.grid(True, alpha=0.3, which="both")

    xs = [p for p, r in valid]
    ys = [r for p, r in valid]
    ax_lin.plot(xs, ys, "o-", color="tab:blue", lw=2, ms=7, zorder=3)
    mark_reference(ax_lin, ref_x=ref_val, ref_label=ref_val, for_log=False)
    ax_lin.axhline(0.01, color="red", ls=":", lw=1.2)
    ax_lin.axhline(0.001, color="orange", ls=":", lw=1.0)
    if preset_val is not None:
        ax_lin.axvline(preset_val, color="tab:green", ls="--", lw=1.2)
    for x, y in zip(xs, ys):
        ax_lin.annotate(f"{y:.4f}", (x, y), textcoords="offset points",
                        xytext=(0, 6), ha="center", fontsize=6)
    ax_lin.set_xlabel(xlabel)
    ax_lin.set_ylabel("RMS vs ref (s)")
    ax_lin.set_title(title_lin, fontsize=10)
    ax_lin.set_xticks(xs)
    ax_lin.grid(True, alpha=0.3)
