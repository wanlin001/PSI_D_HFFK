#!/usr/bin/env python3
"""
gen_benchmark_psitomo.py — PSI_D benchmark 合成各向異性模型生成器

輸出 4 個 psitomo 檔 + 1 個解析解 CSV：
  bench_models/bench_1L_A.dat  — 單層，粗網格 5×5×11, dep 60–360 km
  bench_models/bench_1L_B.dat  — 單層，細網格 9×9×21, dep 30–630 km
  bench_models/bench_2L_A.dat  — 雙層，粗網格 5×5×11, dep 60–360 km
  bench_models/bench_2L_B.dat  — 雙層，細網格 9×9×21, dep 30–630 km
  bench_models/analytical_si.csv

異向性設計（均一 δVs = 0.101 km/s = 2.25% of Vs0）：
  單層: φ = 45°（NE 方向），整個深度範圍
  雙層: 上半層 φ₁ = 0°（N），下半層 φ₂ = 45°（NE）— 差 45°，有 BAZ 依賴性

用法：
  cd /Users/wanlin/Documents/ASPECT/PSI_D_HFFK
  python3 validation/gen_benchmark_psitomo.py

wl 上：
  cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
  git pull && python3 validation/gen_benchmark_psitomo.py
  cp -r validation/bench_models /home/wl/work/ASPECT/<PROJECT>/
"""

import numpy as np
from pathlib import Path

# ── 背景物性（地幔典型值）──────────────────────────────────────────
VP0 = 8.0    # km/s
VS0 = 4.5    # km/s
RHO = 3.3    # g/cm³   [Cij unit = ρV² → GPa when ρ in g/cm³, V in km/s]
DVS = 0.101  # km/s — 各向異性大小 (δVs/Vs0 ≈ 2.25%)

# 快軸方位角定義（從 North 順時針，°）
PHI_1L       = 45.0   # 單層
PHI_2L_UPPER = 0.0    # 雙層上半層 (North)
PHI_2L_LOWER = 45.0   # 雙層下半層 (NE)   — 差 45°，不互相抵消


# ══════════════════════════════════════════════════════════════════
# Cij 構造
# ══════════════════════════════════════════════════════════════════

def cij_hti(phi_deg: float, dVs: float = DVS) -> np.ndarray:
    """
    HTI 各向異性彈性張量（Voigt 上三角 21 值 + rho = 22 值，GPa / g/cm³）

    快軸在 (East, North) 水平面，方位角 phi_deg 從 North 順時針。
    快軸方向向量 ê_fast = (sin φ, cos φ) in (E, N)

    對垂直傳播（n = (0,0,1)），Christoffel 矩陣 Γ_ij = C_i3j3/ρ：
      Γ₁₁ = C₅₅/ρ  (East 偏振)
      Γ₂₂ = C₄₄/ρ  (North 偏振)
      Γ₁₂ = C₄₅/ρ  (非零 → 快軸非正好 E 或 N 時)

    本徵值 = Vs_fast², Vs_slow²，本徵向量 = ê_fast, ê_slow：
      C₅₅ = ρ(Vs_fast²·sin²φ + Vs_slow²·cos²φ)
      C₄₄ = ρ(Vs_fast²·cos²φ + Vs_slow²·sin²φ)
      C₄₅ = ρ(Vs_fast² − Vs_slow²)·sinφ·cosφ

    其餘 Cij 維持等向性值。
    """
    lam  = RHO * (VP0**2 - 2*VS0**2)   # Lamé λ ≈ 77.6 GPa
    mu   = RHO * VS0**2                  # Lamé μ ≈ 66.8 GPa
    C11  = lam + 2*mu                    # ≈ 211.2 GPa

    Vf = VS0 + dVs/2
    Vs = VS0 - dVs/2

    ph      = np.radians(phi_deg)
    sf, cf  = np.sin(ph), np.cos(ph)

    C44 = RHO * (Vf**2 * cf**2 + Vs**2 * sf**2)
    C55 = RHO * (Vf**2 * sf**2 + Vs**2 * cf**2)
    C45 = RHO * (Vf**2 - Vs**2) * sf * cf

    # Voigt 上三角（PSI_D psitomo 順序）：
    #   idx 0–5:   C11 C12 C13 C14 C15 C16
    #   idx 6–10:  C22 C23 C24 C25 C26
    #   idx 11–14: C33 C34 C35 C36
    #   idx 15–17: C44 C45 C46
    #   idx 18–19: C55 C56
    #   idx 20:    C66
    #   idx 21:    rho
    cij = [
        C11, lam, lam,  0.,  0.,  0.,   # C11 C12 C13 C14 C15 C16
             C11, lam,  0.,  0.,  0.,   # C22 C23 C24 C25 C26
                  C11,  0.,  0.,  0.,   # C33 C34 C35 C36
                        C44, C45, 0.,   # C44 C45 C46
                             C55,  0.,  # C55 C56
                                   mu,  # C66
        RHO                             # rho
    ]
    return np.array(cij, dtype=float)


def delta_t_km(thick_km: float, dVs: float = DVS) -> float:
    """垂直射線穿 thick_km 公里的分裂時間 δt (s)。
    一階近似：δt ≈ thick_km · δVs / Vs0²"""
    return thick_km * dVs / VS0**2


# ══════════════════════════════════════════════════════════════════
# psitomo 格式輸出
# ══════════════════════════════════════════════════════════════════

def write_psitomo(out_path, lon_arr, lat_arr, dep_arr, cij_func):
    """
    cij_func(dep_km) → ndarray(22) [C11..C66, rho]
    迴圈順序：depth-outer → lat-mid → lon-inner（與 VIZTOMO 一致）
    """
    nlon, nlat, ndep = len(lon_arr), len(lat_arr), len(dep_arr)
    dlon = float(lon_arr[1] - lon_arr[0]) if nlon > 1 else 1.0
    dlat = float(lat_arr[1] - lat_arr[0]) if nlat > 1 else 1.0
    ddep = float(dep_arr[1] - dep_arr[0]) if ndep > 1 else 1.0
    clon = float(lon_arr[nlon // 2])
    clat = float(lat_arr[nlat // 2])

    with open(out_path, 'w') as f:
        f.write(f"  6370.0000, {clon:9.4f}, {clat:9.4f},   0.0000,\n")
        f.write(f"    {nlon},    {nlat},    {ndep},\n")
        f.write(f"  {dlon:.4f},   {dlat:.4f},  {ddep:.4f},\n")
        f.write("   0,\n")
        for dep in dep_arr:
            cij_rho = cij_func(float(dep))
            for lat in lat_arr:
                for lon in lon_arr:
                    row = np.concatenate([[lon, lat, dep], cij_rho])
                    f.write(','.join(f' {v:>14.5E}' for v in row) + ',\n')

    sz = Path(out_path).stat().st_size / 1e6
    print(f"  → {out_path}  ({nlon}×{nlat}×{ndep}={nlon*nlat*ndep:,} pts, {sz:.1f} MB)")


# ══════════════════════════════════════════════════════════════════
# 解析解
# ══════════════════════════════════════════════════════════════════

def si_1layer(phi_deg, dt_s, baz_arr):
    """單層解析 SI。SI = δt·sin(2(φ−BAZ))"""
    return dt_s * np.sin(2 * (np.radians(phi_deg) - np.radians(baz_arr)))


def si_2layer(phi1_deg, dt1_s, phi2_deg, dt2_s, baz_arr, period_s=8.0):
    """
    雙層解析 SI — Silver & Savage (1994) 精確公式。
    Layer 1 = 淺層（波最後穿過），Layer 2 = 深層（波先穿過）。

    複數分裂：
      H = e^{2iφ₁}·sin(πδt₁/T)·cos(πδt₂/T)
        + e^{2iφ₂}·sin(πδt₂/T)·cos(πδt₁/T)

    視分裂參數（頻率依賴）：
      δt_app(T) = (T/π)·arcsin(|H|)
      φ_app(T)  = arg(H)/2

    SI(BAZ, T) = δt_app·sin(2(φ_app − BAZ))
    """
    omega = 2 * np.pi / period_s
    ph1   = np.radians(phi1_deg)
    ph2   = np.radians(phi2_deg)
    baz   = np.radians(baz_arr)

    H = (np.exp(2j * ph1) * np.sin(omega * dt1_s / 2) * np.cos(omega * dt2_s / 2) +
         np.exp(2j * ph2) * np.sin(omega * dt2_s / 2) * np.cos(omega * dt1_s / 2))

    absH   = np.clip(np.abs(H), 0., 1.)
    dt_app = (2. / omega) * np.arcsin(absH)
    ph_app = 0.5 * np.angle(H)

    return dt_app * np.sin(2 * (ph_app - baz))


# ══════════════════════════════════════════════════════════════════
# 主程式
# ══════════════════════════════════════════════════════════════════

def main():
    out_dir = Path(__file__).parent / "bench_models"
    out_dir.mkdir(exist_ok=True)

    # ── 網格定義（台灣為中心）────────────────────────────────────
    # Grid A: 粗, 5×5×11, dep 60-360 km, ddep=30 km
    lon_A = np.linspace(121., 125., 5)
    lat_A = np.linspace(22.,  26.,  5)
    dep_A = np.linspace(60.,  360., 11)
    thk_A = float(dep_A[-1] - dep_A[0])   # 300 km

    # Grid B: 細, 9×9×21, dep 30-630 km, ddep=30 km
    lon_B = np.linspace(119., 127., 9)
    lat_B = np.linspace(20.,  28.,  9)
    dep_B = np.linspace(30.,  630., 21)
    thk_B = float(dep_B[-1] - dep_B[0])   # 600 km

    # 各種組合的 δt（均一 DVs = 0.101 km/s）
    dt_1A  = delta_t_km(thk_A)          # 單層 Grid A: 1.50 s
    dt_1B  = delta_t_km(thk_B)          # 單層 Grid B: 3.00 s
    dt_2A  = delta_t_km(thk_A / 2)      # 雙層各半 Grid A: 0.75 s each
    dt_2B  = delta_t_km(thk_B / 2)      # 雙層各半 Grid B: 1.50 s each

    print("=" * 60)
    print("PSI_D Benchmark — 合成各向異性模型生成")
    print("=" * 60)
    print(f"δVs = {DVS} km/s  ({100*DVS/VS0:.2f}% of Vs0={VS0} km/s)")
    print()
    print(f"  bench_1L_A: 單層 φ={PHI_1L}°, δt={dt_1A:.3f}s  (300 km / {VS0}²)")
    print(f"  bench_1L_B: 單層 φ={PHI_1L}°, δt={dt_1B:.3f}s  (600 km / {VS0}²)")
    print(f"  bench_2L_A: 雙層 φ₁={PHI_2L_UPPER}°+φ₂={PHI_2L_LOWER}°, δt₁=δt₂={dt_2A:.3f}s each")
    print(f"  bench_2L_B: 雙層 φ₁={PHI_2L_UPPER}°+φ₂={PHI_2L_LOWER}°, δt₁=δt₂={dt_2B:.3f}s each")
    print()

    # ── 單層模型 ──────────────────────────────────────────────────
    print("=== 單層模型 ===")
    cij_1L = cij_hti(PHI_1L)
    write_psitomo(out_dir / "bench_1L_A.dat", lon_A, lat_A, dep_A,
                  lambda dep: cij_1L)
    write_psitomo(out_dir / "bench_1L_B.dat", lon_B, lat_B, dep_B,
                  lambda dep: cij_1L)

    # ── 雙層模型 ──────────────────────────────────────────────────
    print("\n=== 雙層模型 ===")
    cij_upper = cij_hti(PHI_2L_UPPER)
    cij_lower = cij_hti(PHI_2L_LOWER)

    def make_2L(dep_arr):
        mid = (dep_arr[0] + dep_arr[-1]) / 2.0
        def f(dep):
            return cij_upper if dep <= mid else cij_lower
        return f

    write_psitomo(out_dir / "bench_2L_A.dat", lon_A, lat_A, dep_A, make_2L(dep_A))
    write_psitomo(out_dir / "bench_2L_B.dat", lon_B, lat_B, dep_B, make_2L(dep_B))

    # ── 解析解 CSV ────────────────────────────────────────────────
    print("\n=== 解析解 CSV ===")
    baz     = np.arange(0., 361., 5.)
    periods = [4., 8., 16., 20., 25., 33., 50.]

    rows = []
    for b in baz:
        row = [b,
               si_1layer(PHI_1L, dt_1A, b),
               si_1layer(PHI_1L, dt_1B, b)]
        for T in periods:
            row.append(si_2layer(PHI_2L_UPPER, dt_2A,
                                 PHI_2L_LOWER,  dt_2A, b, T))
        for T in periods:
            row.append(si_2layer(PHI_2L_UPPER, dt_2B,
                                 PHI_2L_LOWER,  dt_2B, b, T))
        rows.append(row)

    hdr = (["BAZ_deg", "SI_1L_A", "SI_1L_B"] +
           [f"SI_2L_A_T{int(T)}s" for T in periods] +
           [f"SI_2L_B_T{int(T)}s" for T in periods])
    np.savetxt(out_dir / "analytical_si.csv", rows,
               delimiter=',', header=','.join(hdr), comments='', fmt='%.6f')
    print(f"  → {out_dir}/analytical_si.csv  ({len(rows)} BAZ steps × {len(hdr)} columns)")

    print("\n=== Done ===")
    print(f"  模型在: {out_dir}/")
    print(f"  wl: cp -r validation/bench_models /home/wl/work/ASPECT/<PROJECT>/")
    print(f"       sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench.sh")


if __name__ == "__main__":
    main()
