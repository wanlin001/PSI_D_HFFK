#!/usr/bin/env python3
"""
gen_benchmark_psitomo.py — PSI_D benchmark 合成各向異性模型生成器

輸出 6 個 psitomo 檔 + 1 個解析解 CSV：
  bench_1L_A.dat    — 單層均勻，粗 5×5×11, 119-127°E 20-28°N dep 60-360 km (2°spacing)
  bench_1L_B.dat    — 單層均勻，細 9×9×21, 119-127°E 20-28°N dep 30-630 km (1°spacing)
  bench_2L_A.dat    — 雙層均勻，粗 5×5×11 (同 A)
  bench_2L_B.dat    — 雙層均勻，細 9×9×21 (同 B)           (two-layer interference)
  bench_lateral_B.dat — 橫向邊界，細 9×9×21                (HFFK period 差異最明顯)
  bench_lateral_B_fine.dat — 橫向邊界，超細 17×17×21       (確認 resolution 無影響)
  analytical_si.csv

橫向模型 (bench_lateral_B)：
  左半邊 (lon < 123°E): φ_L=0°  (North)
  右半邊 (lon ≥ 123°E): φ_R=90° (East)
  接收站在中心 (123°E, 24°N) — 正好在邊界上

  HFFK period 效應機制：
    Fresnel zone 半徑 R_F ≈ √(Vs·T·(L-z)·z/L)
    T=4s,  dep=330km: R_F ≈ √(4.5×4×330×300/630) ≈  90 km  → 接近 ray theory
    T=50s, dep=330km: R_F ≈ √(4.5×50×330×300/630) ≈ 333 km → 左右兩側都在 zone 內
  → T 越大，φ_eff 越往 45° 靠（兩側平均）

用法：
  python3 validation/gen_benchmark_psitomo.py

wl 上：
  cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
  git pull && python3 validation/gen_benchmark_psitomo.py
  mkdir -p /lfs/wl/bench_psi && cd /lfs/wl/bench_psi
  sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench.sh
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
PHI_LAT_L    = 0.0    # 橫向左半 (North)
PHI_LAT_R    = 90.0   # 橫向右半 (East)   — 互相垂直，Fresnel zone 效應最明顯
LON_BOUNDARY = 123.0  # 橫向邊界經度


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

def write_psitomo_3d(out_path, lon_arr, lat_arr, dep_arr, cij_func3d):
    """
    cij_func3d(lon, lat, dep) → ndarray(22) [C11..C66, rho]
    用於橫向非均勻模型。
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
            for lat in lat_arr:
                for lon in lon_arr:
                    cij_rho = cij_func3d(float(lon), float(lat), float(dep))
                    row = np.concatenate([[lon, lat, dep], cij_rho])
                    f.write(','.join(f' {v:>14.5E}' for v in row) + ',\n')

    sz = Path(out_path).stat().st_size / 1e6
    print(f"  → {out_path}  ({nlon}×{nlat}×{ndep}={nlon*nlat*ndep:,} pts, {sz:.1f} MB)")


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
# 獨立 psi_input 生成（benchmark 專用，不依賴真實台灣網絡）
# ══════════════════════════════════════════════════════════════════

CENTER_LON = 123.0
CENTER_LAT = 24.0
N_BAZ      = 16
DELTA_DEG  = 100.0   # 震央距（度）


def _src_position(baz_deg, delta_deg, rcv_lon, rcv_lat):
    """大圓弧：由接收站位置 + BAZ + 震央距 → 震源 (lon, lat)。"""
    d   = np.radians(delta_deg)
    phi = np.radians(rcv_lat)
    lam = np.radians(rcv_lon)
    az  = np.radians(baz_deg)
    phi_s = np.arcsin(np.sin(phi)*np.cos(d) + np.cos(phi)*np.sin(d)*np.cos(az))
    lam_s = lam + np.arctan2(np.sin(az)*np.sin(d)*np.cos(phi),
                              np.cos(d) - np.sin(phi)*np.sin(phi_s))
    return float(np.degrees(lam_s)), float(np.degrees(phi_s))


def generate_bench_psi_input(out_dir, sp_period=50.0):
    """
    產生 benchmark 專用 psi_input（與真實台灣測站無關）：
      Receivers.dat  — 9 虛擬接收站 3×3 grid ±0.1° 圍繞 (123°E,24°N), dep=0
      Sources.dat    — 16 均勻 BAZ 震源，Δ=100°
      DUMMY_SI.dat   — SI 觀測（paz=0）
      DUMMY_SP.dat   — SP 觀測（paz=0, period={sp_period}s）
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 接收站：3×3 in ±0.1°（都在模型範圍內，dep=0）
    rcv = []
    for i, dy in enumerate([-0.1, 0.0, 0.1]):
        for j, dx in enumerate([-0.1, 0.0, 0.1]):
            rcv.append((f"B{i*3+j+1:02d}",
                        CENTER_LON + dx, CENTER_LAT + dy, 0.0))
    with open(out_dir / "Receivers.dat", "w") as f:
        for code, lo, la, dep in rcv:
            f.write(f"{code}, {lo:.6f}, {la:.6f}, {dep:.1f}\n")
    print(f"  [{len(rcv):3d} receivers] → {out_dir}/Receivers.dat")

    # 震源：16 均勻 BAZ
    baz_arr = np.arange(0., 360., 360. / N_BAZ)
    src = []
    for i, baz in enumerate(baz_arr):
        lo, la = _src_position(baz, DELTA_DEG, CENTER_LON, CENTER_LAT)
        src.append((i + 1, lo, la, -100.0))
    with open(out_dir / "Sources.dat", "w") as f:
        for sid, lo, la, dep in src:
            f.write(f"{sid}, {lo:.6f}, {la:.6f}, {dep:.1f}\n")
    print(f"  [{len(src):3d} sources   ] → {out_dir}/Sources.dat")

    # DUMMY 觀測檔
    rcv_ids = [r[0] for r in rcv]
    src_ids = [s[0] for s in src]
    n_obs   = len(src_ids) * len(rcv_ids)
    with open(out_dir / "DUMMY_SI.dat", "w") as f:
        for sid in src_ids:
            for rid in rcv_ids:
                f.write(f"0.0, 0.1, 0.0, SKS, {sid}, {rid}, ???, 0.0\n")
    print(f"  [{n_obs:5d} obs     ] → {out_dir}/DUMMY_SI.dat")
    with open(out_dir / "DUMMY_SP.dat", "w") as f:
        for sid in src_ids:
            for rid in rcv_ids:
                f.write(f"0.0, 0.0, 0.1, 0.1, {sp_period}, SKS, {sid}, {rid}, ???, 0.0\n")
    print(f"  [{n_obs:5d} obs     ] → {out_dir}/DUMMY_SP.dat  (period_s={sp_period})")


# ══════════════════════════════════════════════════════════════════
# 主程式
# ══════════════════════════════════════════════════════════════════

def main():
    out_dir = Path(__file__).parent / "bench_models"
    out_dir.mkdir(exist_ok=True)

    # ── 網格定義（台灣為中心 123°E, 24°N）──────────────────────
    # 全部模型從 dep=0 開始，確保地表接收站在模型範圍內
    # Grid A/B 覆蓋「相同」地理+深度範圍，只差網格疏密 → 真正的解析度測試
    #   → 1L_A 與 1L_B 是同一物理模型，SI 必須重疊（δt 相同）
    #
    # Grid A: 粗, 5×5×13,  dep 0-600 km (step 50km), lon/lat 2° spacing
    # Grid B: 細, 9×9×21,  dep 0-600 km (step 30km), lon/lat 1° spacing
    lon_A = np.linspace(119., 127., 5)
    lat_A = np.linspace(20.,  28.,  5)
    dep_A = np.linspace(0., 600., 13)           # 13 pts: 0-600km, step 50km
    thk_A = float(dep_A[-1] - dep_A[0])         # 600 km

    lon_B = np.linspace(119., 127., 9)
    lat_B = np.linspace(20.,  28.,  9)
    dep_B = np.arange(0., 601., 30.)            # 21 pts: 0-600km, step 30km
    thk_B = float(dep_B[-1] - dep_B[0])         # 600 km

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

    # ── 橫向非均勻模型（period 效應診斷）────────────────────────
    # 接收站在邊界 LON_BOUNDARY=123°E
    # 左半 lon < 123°: φ=0°(N)；右半 lon ≥ 123°: φ=90°(E)
    # 只用 Grid B（細網格）+ 超細 17×17 做解析度驗證
    print("\n=== 橫向非均勻模型（HFFK period 效應診斷）===")
    cij_lat_L = cij_hti(PHI_LAT_L)
    cij_lat_R = cij_hti(PHI_LAT_R)

    def make_lateral(lon_arr):
        def f_lat(lon, lat, dep):
            return cij_lat_L if lon < LON_BOUNDARY else cij_lat_R
        return f_lat

    # Grid B: 9×9×21 細網格
    write_psitomo_3d(out_dir / "bench_lateral_B.dat",
                     lon_B, lat_B, dep_B, make_lateral(lon_B))

    # Grid C: 17×17×21 超細格，間距 0.5°（解析度驗證）
    lon_C = np.linspace(119., 127., 17)
    lat_C = np.linspace(20.,  28.,  17)
    dep_C = dep_B.copy()
    write_psitomo_3d(out_dir / "bench_lateral_C.dat",
                     lon_C, lat_C, dep_C, make_lateral(lon_C))

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

    # ── 獨立 psi_input ────────────────────────────────────────────
    print("\n=== Benchmark psi_input（獨立，不依賴台灣真實網絡）===")
    inp_dir = Path(__file__).parent / "bench_psi_input"
    generate_bench_psi_input(inp_dir, sp_period=50.0)

    print("\n=== Done ===")
    print(f"  模型在: {out_dir}/")
    print(f"  psi_input: {inp_dir}/")
    print(f"  wl: cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK")
    print(f"       git pull && python3 validation/gen_benchmark_psitomo.py")
    print(f"       mkdir -p /lfs/wl/bench_psi && cd /lfs/wl/bench_psi")
    print(f"       sbatch /home/wl/.../validation/job_bench.sh")


if __name__ == "__main__":
    main()
