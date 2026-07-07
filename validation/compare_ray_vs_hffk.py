#!/usr/bin/env python3
"""
compare_ray_vs_hffk.py
比較 Ray Theory vs HFFK 的 splitting intensity 預測
用法: python compare_ray_vs_hffk.py
"""
import numpy as np
import sys
import os

def read_si_dat(path):
    """讀取 PSI_D 輸出的 SplittingIntensity .dat 檔"""
    si_vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [p.strip() for p in line.split(',')]
            si_vals.append(float(parts[0]))
    return np.array(si_vals)

# ── 路徑設定 ──
base = os.path.dirname(os.path.abspath(__file__))
ray_dir  = os.path.join(base, "../psi_output/SYN_SinkingBlock")
hffk_dir = os.path.join(base, "psi_output/HFFK_SinkingSlab")

si_file = "SYN_SplittingIntensity_ShearWave.dat"

ray_path  = os.path.join(ray_dir, si_file)
hffk_path = os.path.join(hffk_dir, si_file)

if not os.path.exists(ray_path):
    print(f"ERROR: Ray theory output not found: {ray_path}")
    print("請先執行 run_ray_forward.jl")
    sys.exit(1)
if not os.path.exists(hffk_path):
    print(f"ERROR: HFFK output not found: {hffk_path}")
    print("請先執行 run_hffk_forward.jl")
    sys.exit(1)

ray  = read_si_dat(ray_path)
hffk = read_si_dat(hffk_path)

if len(ray) != len(hffk):
    print(f"WARNING: 筆數不同 ray={len(ray)}, hffk={len(hffk)}")
    n = min(len(ray), len(hffk))
    ray, hffk = ray[:n], hffk[:n]

diff = hffk - ray

print("\n" + "="*55)
print("  Ray Theory vs HFFK  —  Splitting Intensity 比較")
print("="*55)
print(f"  樣本數           : {len(ray)}")
print(f"  Ray theory SI    : mean={ray.mean():.4f}, std={ray.std():.4f}")
print(f"  HFFK SI          : mean={hffk.mean():.4f}, std={hffk.std():.4f}")
print(f"  差值 (HFFK-Ray)  : mean={diff.mean():.4f}, std={diff.std():.4f}")
print(f"  RMS 差異         : {np.sqrt(np.mean(diff**2)):.4f}")
print(f"  最大正差         : {diff.max():.4f}")
print(f"  最大負差         : {diff.min():.4f}")
print("="*55)
print("\n預期結果：HFFK 與 Ray theory 整體趨勢相同，")
print("但在橫向速度梯度大的區域（sinking slab 邊界）HFFK 會更平滑。")
print("\n引用: VanderBeek & Faccenda (2021, GJI) — HFFK RMSE ≈ Born kernel")
