#!/usr/bin/env python3
"""
plot_convergence_common.py — 收斂測試共用工具

供 plot_nrings_convergence.py / plot_nazimuth_convergence.py /
plot_benchmark_tier2.py 共用，確保 RMS 計算方式一致。

規則：一律用 **raw SI**（全觀測點，不經 BAZ 平均）。
"""

from pathlib import Path

import numpy as np


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
