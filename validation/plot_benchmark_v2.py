#!/usr/bin/env python3
"""
plot_benchmark_v2.py — PSI_D Benchmark 結果可視化

用法（在 PSI_DIR 執行）：
    cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
    python3 validation/plot_benchmark_v2.py
    # 輸出在 validation/bench_output/benchmark_v2.png

輸出圖（4欄×3列）：
  Row 0: SI vs BAZ — 解析解 vs Ray SP→SI vs HFFK T=4/8/16/20/25/33/50s
  Row 1: (HFFK_T − Ray)/Ray 相對差，顯示頻率依賴
  Row 2: 各方法 r（與解析解 Pearson 相關係數）vs period

4 欄對應 bench_1L_A / bench_1L_B / bench_2L_A / bench_2L_B
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm

PERIODS   = [4., 8., 16., 20., 25., 33., 50.]
MODELS    = ["bench_1L_A", "bench_1L_B", "bench_2L_A", "bench_2L_B"]
MODEL_LAB = ["1L-Coarse\n5×5×11, 60–360 km",
             "1L-Fine\n9×9×21, 30–630 km",
             "2L-Coarse\n5×5×11, 60–360 km",
             "2L-Fine\n9×9×21, 30–630 km"]

COLORS_HFFK = cm.plasma(np.linspace(0.1, 0.9, len(PERIODS)))


# ══════════════════════════════════════════════════════════════════
# I/O helpers
# ══════════════════════════════════════════════════════════════════

def read_sp(path) -> tuple:
    """讀 SplittingParameters: (dt_arr, phi_arr)"""
    dt_v, phi_v = [], []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            dt_v.append(float(p[0]))
            phi_v.append(float(p[1]))
    return np.array(dt_v), np.array(phi_v)


def read_si(path) -> np.ndarray:
    """讀 SplittingIntensity: SI 值陣列"""
    si = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            si.append(float(ln.split(',')[0].strip()))
    return np.array(si)


def read_obs(dummy_sp_path) -> tuple:
    """從 DUMMY_SP 讀 (baz_arr, src_id_arr, rcv_id_arr)"""
    baz_list, src_list, rcv_list = [], [], []
    with open(dummy_sp_path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            baz_list.append(float(p[-1]))   # paz_rad ≈ 0 → BAZ ≈ p[-1]?
            src_list.append(p[6])
            rcv_list.append(p[7])
    return np.array(baz_list), src_list, rcv_list


def read_sources_baz(sources_dat) -> dict:
    """讀 Sources.dat，返回 {src_id: baz_deg}"""
    baz = {}
    with open(sources_dat) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            sid = p[0].strip()
            baz[sid] = float(p[3])   # col 4 = BAZ (deg), 1-indexed columns
    return baz


def sp_to_si(dt_arr, phi_arr, baz_arr) -> np.ndarray:
    """SP (δt, φ) → SI = δt·sin(2(φ−BAZ))"""
    return dt_arr * np.sin(2 * (phi_arr - baz_arr))


def read_analytical_csv(csv_path, model_key, period_s=None):
    """
    從 analytical_si.csv 讀解析解。
    model_key: '1L_A', '1L_B', '2L_A', '2L_B'
    period_s: 若為 None 且是 1L → 與 period 無關
    """
    data = np.genfromtxt(csv_path, delimiter=',', names=True)
    baz_deg = data['BAZ_deg']
    if '1L' in model_key:
        col = f"SI_1L_{model_key[-1]}"
        return baz_deg, data[col]
    else:
        T = int(period_s) if period_s else 8
        col = f"SI_2L_{model_key[-1]}_T{T}s"
        # handle 33s special case
        if col not in data.dtype.names:
            col = f"SI_2L_{model_key[-1]}_T33s"
        return baz_deg, data[col]


# ══════════════════════════════════════════════════════════════════
# 主程式
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    PSI_DIR_DEFAULT = "/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK"
    parser.add_argument("--psi-dir", default=PSI_DIR_DEFAULT)
    parser.add_argument("--src", type=int, default=None, help="Only plot this source ID")
    parser.add_argument("--out", default=None,
                        help="Output figure prefix (default: validation/bench_output/benchmark_v2)")
    args = parser.parse_args()

    psi   = Path(args.psi_dir)
    out_p = Path(args.out) if args.out else psi / "validation" / "bench_output" / "benchmark_v2"
    out_p.parent.mkdir(parents=True, exist_ok=True)

    csv_path    = psi / "validation" / "bench_models" / "analytical_si.csv"
    dummy_sp    = psi / "psi_input" / "DUMMY_SP_uniform48.dat"
    dummy_si    = psi / "psi_input" / "DUMMY_SI_uniform48.dat"
    sources_dat = psi / "psi_input" / "Sources_uniform48.dat"

    if not csv_path.exists():
        sys.exit(f"ERROR: {csv_path} not found — run gen_benchmark_psitomo.py first")

    # 讀 BAZ 對應表（source ID → BAZ deg）
    src_baz = read_sources_baz(sources_dat)

    # 建 (src_id, rcv_id) → baz 查找（從 DUMMY_SI 讀 src_id）
    si_src_ids = []
    with open(dummy_si) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            si_src_ids.append(str(p[3]))  # SI format: SI, err, 0, phase, src, rcv, ???, paz

    sp_src_ids = []
    sp_baz_arr = []
    with open(dummy_sp) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            sp_src_ids.append(str(p[6]))
            sp_baz_arr.append(float(p[-1]))

    si_baz_arr = np.array([float(src_baz.get(sid, 0)) * np.pi / 180
                           for sid in si_src_ids])
    sp_baz_arr = np.array(sp_baz_arr)

    # ── 繪圖 ──────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 4, figsize=(20, 14), sharex='row')
    fig.suptitle("PSI_D Benchmark — Resolution × Layer × Period",
                 fontsize=14, fontweight='bold')

    r_matrix = np.full((len(MODELS), len(PERIODS)+1), np.nan)  # [model, period_idx]

    for mi, (mod, lab) in enumerate(zip(MODELS, MODEL_LAB)):
        ax0 = axes[0, mi]  # SI vs BAZ
        ax1 = axes[1, mi]  # relative diff
        ax2 = axes[2, mi]  # r vs period

        mod_key = mod.replace("bench_", "")   # e.g. '1L_A'
        is_2L   = "2L" in mod

        # 解析解 (用 T=8s for 2L default)
        baz_an, si_an = read_analytical_csv(csv_path, mod_key,
                                             period_s=8. if is_2L else None)
        ax0.plot(baz_an, si_an, 'k--', lw=1.5, label="Analytical", zorder=5)

        # Ray SP → SI
        ray_sp_dir = psi / "validation" / "bench_output" / f"bench_{mod}_ray"
        ray_sp_dat = ray_sp_dir / "SYN_SplittingParameters_ShearWave.dat"
        if ray_sp_dat.exists():
            dt_r, phi_r = read_sp(ray_sp_dat)
            si_ray = sp_to_si(dt_r, phi_r, sp_baz_arr[:len(dt_r)])
            ax0.plot(np.degrees(sp_baz_arr[:len(si_ray)]), si_ray,
                     'r-o', ms=3, lw=1.5, label="Ray SP→SI", zorder=4)
            r_ray = np.corrcoef(si_ray, np.interp(
                np.degrees(sp_baz_arr[:len(si_ray)]), baz_an, si_an))[0, 1]
            r_matrix[mi, 0] = r_ray

        # HFFK SI 各頻率
        for pi, T in enumerate(PERIODS):
            TT = f"{int(T)}s"
            hffk_dir = psi / "validation" / "bench_output" / f"bench_{mod}_hffk_T{TT}"
            hffk_dat = hffk_dir / "SYN_SplittingIntensity_ShearWave.dat"
            if not hffk_dat.exists():
                continue
            si_h = read_si(hffk_dat)
            baz_deg_h = np.degrees(si_baz_arr[:len(si_h)])
            an_at_h   = np.interp(baz_deg_h, baz_an,
                                   read_analytical_csv(csv_path, mod_key, T)[1]
                                   if is_2L else si_an)
            color = COLORS_HFFK[pi]
            ax0.plot(baz_deg_h, si_h, color=color, lw=1,
                     alpha=0.7, label=f"HFFK T={TT}")
            diff = (si_h - np.interp(baz_deg_h, baz_an, si_an))
            ax1.plot(baz_deg_h, diff, color=color, lw=1, alpha=0.7,
                     label=f"T={TT}")
            r_h = np.corrcoef(si_h, an_at_h)[0, 1]
            r_matrix[mi, pi+1] = r_h

        ax0.set_title(f"{mod}\n{lab}", fontsize=9)
        ax0.set_ylabel("SI (s)" if mi == 0 else "")
        ax0.legend(fontsize=6, ncol=2)
        ax0.set_ylim(-2.5, 2.5)
        ax0.axhline(0, color='gray', lw=0.5)

        ax1.set_ylabel("HFFK−Analytical (s)" if mi == 0 else "")
        ax1.axhline(0, color='gray', lw=0.5)
        ax1.legend(fontsize=6)

        # r vs period
        ax2.axhline(r_matrix[mi, 0], color='r', ls='--', label=f"Ray r={r_matrix[mi,0]:.3f}")
        ax2.plot(PERIODS, r_matrix[mi, 1:], 'o-', color='navy', label="HFFK r")
        ax2.set_xlabel("Period (s)")
        ax2.set_ylabel("r (vs Analytical)" if mi == 0 else "")
        ax2.set_ylim(0, 1.05)
        ax2.set_xscale('log')
        ax2.legend(fontsize=7)

    axes[0, 0].set_xlabel("BAZ (°)")
    axes[1, 0].set_xlabel("BAZ (°)")

    plt.tight_layout()
    fig_path = str(out_p) + ".png"
    plt.savefig(fig_path, dpi=150)
    print(f"Saved: {fig_path}")
    plt.close()


if __name__ == "__main__":
    main()
