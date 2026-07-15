#!/usr/bin/env python3
"""
plot_real_period_map.py — 真實 Kuo2018 模型的 period 差異地圖

在台灣網格上畫 HFFK SplittingIntensity 隨 period 的變化：
  每個接收點 (21×25 網格) 的 ΔSI = SI(T_hi) − SI(T_lo)，
  對所有事件取 RMS，顯示「哪裡的 period 效應最強」。
  （period 效應 = 有限頻率 Fresnel zone 掃過橫向異向性梯度的結果；
    均勻區看不到差異，橫向 φ 變化大的地方差異最大。）

輸入（jobD 的輸出，psi_output/hffk_Kuo2018_<PHASE>_T{N}s/）：
  SYN_SplittingIntensity_ShearWave.dat  每個 period 一個資料夾
資料列順序 = DUMMY 順序：事件(outer) → 接收點(inner)。

用法（在 project dir，即 jobD 的 SLURM_SUBMIT_DIR 執行）：
  python3 /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/scripts/plot_real_period_map.py \
      --project . --t-lo 4 --t-hi 50
  # 輸出 psi_output/period_map_T4_vs_T50.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PSI_DIR_DEFAULT = "/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK"


def read_receiver_grid(rcv_file):
    """讀 Receivers.dat → (rcv_ids, lon_arr, lat_arr, nlon, nlat)。
    假設規則網格，lon 變化最快（每列 nlon 個），再換 lat。"""
    ids, lons, lats = [], [], []
    with open(rcv_file) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            p = [x.strip() for x in ln.split(',')]
            ids.append(p[0]); lons.append(float(p[1])); lats.append(float(p[2]))
    lons = np.array(lons); lats = np.array(lats)
    ulon = np.unique(np.round(lons, 4))
    ulat = np.unique(np.round(lats, 4))
    return ids, lons, lats, ulon, ulat


def read_si(path):
    """讀 SYN SplittingIntensity 第 1 欄。"""
    vals = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            vals.append(float(ln.split(',')[0].strip()))
    return np.array(vals)


def si_dir(project, phase, period):
    return (project / "psi_output"
            / f"hffk_Kuo2018_{phase}_T{int(period)}s"
            / "SYN_SplittingIntensity_ShearWave.dat")


def per_receiver_rms(si_lo, si_hi, n_rcv):
    """ΔSI = si_hi − si_lo，reshape (n_event, n_rcv)，對事件取 RMS → (n_rcv,)。"""
    n = min(len(si_lo), len(si_hi))
    n_ev = n // n_rcv
    if n_ev == 0:
        return None, 0
    d = (si_hi[:n_ev*n_rcv] - si_lo[:n_ev*n_rcv]).reshape(n_ev, n_rcv)
    rms = np.sqrt(np.nanmean(d**2, axis=0))
    return rms, n_ev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default=".", help="jobD 的 project dir（含 psi_output/）")
    ap.add_argument("--psi-dir", default=PSI_DIR_DEFAULT)
    ap.add_argument("--phases", nargs="+", default=["SKS", "SKKS"])
    ap.add_argument("--t-lo", type=float, default=4.0, help="短週期（ray 極限）")
    ap.add_argument("--t-hi", type=float, default=50.0, help="長週期")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    project = Path(args.project).resolve()
    rcv_file = Path(args.psi_dir) / "psi_input" / "Receivers.dat"
    if not rcv_file.exists():
        sys.exit(f"ERROR: {rcv_file} not found")

    ids, lons, lats, ulon, ulat = read_receiver_grid(rcv_file)
    n_rcv = len(ids)
    nlon, nlat = len(ulon), len(ulat)
    print(f"Receivers: {n_rcv} ({nlon} lon × {nlat} lat), "
          f"lon {ulon[0]:.2f}-{ulon[-1]:.2f}, lat {ulat[0]:.2f}-{ulat[-1]:.2f}")

    phases = [ph for ph in args.phases
              if si_dir(project, ph, args.t_lo).exists()
              and si_dir(project, ph, args.t_hi).exists()]
    if not phases:
        sys.exit(f"ERROR: no SI output found for T={args.t_lo}/{args.t_hi}s under "
                 f"{project}/psi_output/hffk_Kuo2018_*")

    fig, axes = plt.subplots(1, len(phases), figsize=(7*len(phases), 6),
                             squeeze=False)
    fig.suptitle(
        f"Kuo2018 real model — period sensitivity of SI\n"
        f"per-receiver RMS over events of  SI(T={int(args.t_hi)}s) - SI(T={int(args.t_lo)}s)",
        fontsize=12, fontweight='bold')

    # pcolormesh 網格邊界
    def edges(u):
        d = np.diff(u).mean()
        return np.concatenate([u - d/2, [u[-1] + d/2]])
    LON, LAT = np.meshgrid(edges(ulon), edges(ulat))

    vmax = 0.0
    grids = {}
    for ph in phases:
        si_lo = read_si(si_dir(project, ph, args.t_lo))
        si_hi = read_si(si_dir(project, ph, args.t_hi))
        rms, n_ev = per_receiver_rms(si_lo, si_hi, n_rcv)
        if rms is None:
            print(f"  [WARN] {ph}: obs 數不足 ({len(si_lo)} < {n_rcv})"); continue
        # reshape 到 (nlat, nlon)；接收點順序 = lon 最快
        grid = rms[:nlat*nlon].reshape(nlat, nlon)
        grids[ph] = (grid, n_ev)
        vmax = max(vmax, np.nanpercentile(grid, 98))

    for ax, ph in zip(axes[0], phases):
        if ph not in grids:
            continue
        grid, n_ev = grids[ph]
        pcm = ax.pcolormesh(LON, LAT, grid, cmap="magma",
                            vmin=0, vmax=vmax, shading="auto")
        cb = fig.colorbar(pcm, ax=ax, shrink=0.8)
        cb.set_label("RMS ΔSI over events (s)", fontsize=9)
        ax.set_title(f"{ph}  ({n_ev} events)", fontsize=11)
        ax.set_xlabel("Longitude (°E)"); ax.set_ylabel("Latitude (°N)")
        ax.set_aspect("equal")
        # 台灣本島大致範圍參考框
        ax.plot([120.0, 122.0, 122.0, 120.0, 120.0],
                [21.9, 21.9, 25.3, 25.3, 21.9],
                color="cyan", lw=0.8, ls="--", alpha=0.7)

    out = (Path(args.out) if args.out
           else project / "psi_output"
           / f"period_map_T{int(args.t_lo)}_vs_T{int(args.t_hi)}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    plt.savefig(out, dpi=150)
    print(f"Saved: {out}")

    # 也印每個 phase 的整體統計
    for ph, (grid, n_ev) in grids.items():
        print(f"  {ph}: RMS ΔSI  median={np.nanmedian(grid):.4f}s  "
              f"max={np.nanmax(grid):.4f}s  "
              f"@ (lon={ulon[np.nanargmax(grid)%nlon]:.2f}, "
              f"lat={ulat[np.nanargmax(grid)//nlon]:.2f})")


if __name__ == "__main__":
    main()
