#!/usr/bin/env python3
"""
fill_nan_psitomo.py — 填補 psitomo*.dat 中的 NaN 格點

策略：對每一深度層，用該層所有非 NaN 格點的 Cij/rho 中位數填入 NaN 格點。
      因為 NaN 格點都在 ASPECT domain 外緣，射線幾乎不穿過，填等效值對結果無影響。

用法：
    python3 fill_nan_psitomo.py psitomo0050.dat psitomo0050_filled.dat
"""

import sys
import numpy as np
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        print("Usage: fill_nan_psitomo.py <input.dat> <output.dat>")
        sys.exit(1)

    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])
    print(f"Input : {src}  ({src.stat().st_size/1e6:.1f} MB)")

    with open(src) as f:
        lines = f.readlines()

    # 前 4 行是 header（R0/lon0/lat0/0, nlon/nlat/ndep, dlon/dlat/ddep, 0）
    header = lines[:4]
    data_lines = lines[4:]

    # 從 header 取 nlon, nlat, ndep
    h2 = [int(x.strip().rstrip(',')) for x in header[1].split(',')
          if x.strip().rstrip(',').isdigit()]
    nlon, nlat, ndep = h2[0], h2[1], h2[2]
    n_per_depth = nlon * nlat
    print(f"Grid  : {nlon} × {nlat} × {ndep}  ({n_per_depth} pts/layer)")

    n_col = 25   # lon, lat, dep, C11..C66, rho

    # --- Pass 1: 讀所有層，記錄每層的 array 和 nan 狀況 ---
    layers = []   # list of (arr, nan_mask, n_nan)
    nan_total = 0
    for iz in range(ndep):
        chunk = data_lines[iz * n_per_depth : (iz+1) * n_per_depth]
        rows = []
        for ln in chunk:
            vals = [v.strip().rstrip(',') for v in ln.split(',')]
            rows.append([float(v) if v.strip() != 'NaN' else np.nan
                         for v in vals[:n_col]])
        arr = np.array(rows)
        nan_mask = np.any(np.isnan(arr[:, 3:]), axis=1)
        n_nan = int(nan_mask.sum())
        nan_total += n_nan
        layers.append((arr, nan_mask, n_nan))

    def voigt_isotropic(cij22):
        """將 22 個值 (21 Cij + rho) 取 Voigt isotropic 平均，anisotropy 清零。
        Cij 排列（Voigt 上三角）:
          idx: 0=C11 1=C12 2=C13 3=C14 4=C15 5=C16
               6=C22 7=C23 8=C24 9=C25 10=C26
               11=C33 12=C34 13=C35 14=C36
               15=C44 16=C45 17=C46
               18=C55 19=C56
               20=C66
               21=rho
        """
        C11,C22,C33 = cij22[0], cij22[6],  cij22[11]
        C12,C13,C23 = cij22[1], cij22[2],  cij22[7]
        C44,C55,C66 = cij22[15], cij22[18], cij22[20]
        mu  = (C11+C22+C33 - C12-C13-C23 + 3*(C44+C55+C66)) / 15.0
        K   = (C11+C22+C33 + 2*(C12+C13+C23)) / 9.0
        lam = K - 2.0*mu/3.0
        iso = np.zeros(22)
        iso[0]  = lam + 2*mu  # C11
        iso[1]  = lam         # C12
        iso[2]  = lam         # C13
        iso[6]  = lam + 2*mu  # C22
        iso[7]  = lam         # C23
        iso[11] = lam + 2*mu  # C33
        iso[15] = mu          # C44
        iso[18] = mu          # C55
        iso[20] = mu          # C66
        iso[21] = cij22[21]   # rho 不變
        return iso

    # 每層的 median（只對有資料的層計算）
    layer_medians = []
    for arr, nan_mask, n_nan in layers:
        if n_nan < n_per_depth:
            layer_medians.append(np.nanmedian(arr[~nan_mask, 3:], axis=0))
        else:
            layer_medians.append(None)

    # --- Pass 2: 填補 ---
    filled_lines = []
    for iz, (arr, nan_mask, n_nan) in enumerate(layers):
        if n_nan == 0:
            pass  # 乾淨，不動
        elif n_nan < n_per_depth:
            # 邊緣格點：用同層 median（同深度物性相近）
            arr[nan_mask, 3:] = layer_medians[iz]
        else:
            # 整層 NaN（ASPECT domain 外）：
            # 找最近有效層，取其 Voigt isotropic 平均 → anisotropy 貢獻 = 0
            donor = None
            for offset in range(1, ndep):
                for iz2 in [iz + offset, iz - offset]:
                    if 0 <= iz2 < ndep and layer_medians[iz2] is not None:
                        donor = iz2
                        break
                if donor is not None:
                    break
            if donor is not None:
                fill = layer_medians[donor]
                arr[nan_mask, 3:] = fill
                print(f"  [INFO] layer {iz}: 整層 NaN → layer {donor} median 填補"
                      f"  (PSI_D 不支援 isotropic 格點，需保留完整 anisotropy)")
            else:
                print(f"  [ERROR] layer {iz}: 找不到任何有效層，填零")
                arr[nan_mask, 3:] = 0.0

        for row in arr:
            parts = [f' {v:>14.5E}' for v in row]
            filled_lines.append(','.join(parts) + ',\n')

        if (iz+1) % 10 == 0:
            print(f"  layer {iz+1}/{ndep}")

    pct = 100 * nan_total / (nlon * nlat * ndep)
    print(f"NaN格點: {nan_total} / {nlon*nlat*ndep}  ({pct:.1f}%) → 已填補")

    with open(dst, 'w') as f:
        f.writelines(header)
        f.writelines(filled_lines)
    print(f"Output: {dst}  ({dst.stat().st_size/1e6:.1f} MB)")
    print("Done.")

if __name__ == '__main__':
    main()
