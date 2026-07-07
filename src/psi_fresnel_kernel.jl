
# =============================================================================
# HFFK — Heuristic Finite-Frequency slowness Kernel
#
# 方法來源:
#   Schmandt & Humphreys (2010b, Geochem. Geophys. Geosyst., 11, Q05004,
#     doi:10.1029/2010GC003042)  ← eq.(2) 提出 sin 型 kernel 公式
#   VanderBeek & Faccenda (2021, Geophys. J. Int., 225, 2097,
#     doi:10.1093/gji/ggab081)  ← eq.(4)–(5) 推廣至各向異性介質 P 波
#   Faccenda & VanderBeek (2023, J. Geodynamics, 158, 102003,
#     doi:10.1016/j.jog.2023.102003)  ← eq.(2) 同公式擴展至 S 波
#
# 本檔案被 PSI_D.jl 以 include("psi_fresnel_kernel.jl") 載入，
# 函數由 psi_forward.jl 的 psi_forward_hffk() 呼叫。
#
# ── 公式 (VanderBeek & Faccenda 2021, eq.4 / Faccenda & VanderBeek 2023, eq.2) ──
#
#   K(x, r) = Q / (π R_f²(x)) * sin(π r²/R_f²(x))
#
#   其中 Q 為正規化常數，使 ∫ K dV = total ray length L。
#   推導：令 ∫₀^{R_f} K·2πr dr = 1  →  Q = π/2。
#   本實作 Q=1（省略 π/2 因子），caller 應用加權平均（÷ sum of weights），
#   Q 自動消去，結果不受影響。
#
# ── Fresnel zone 半徑 (eq.5) ──
#
#   R_f(x) = sqrt(T * x*(L-x) / (L * u(x)))
#
#   在射線兩端 x=0, x=L 時 R_f=0（donut hole 最大、Fresnel zone 消失）。
#
# ── Donut hole（環形孔） ──
#
#   K(x, r=0) = sin(0)/(π R_f²) = 0
#   射線本身（r=0）的敏感度嚴格為零，與 Born kernel 一致。
#   psi_forward.jl 的 fresnel_sample_offsets 中心點 weight=0，
#   psi_forward_hffk() 以 `continue` 跳過中心點，donut hole 已正確實作。
#
# ── 本檔案匯出的函數 ──
#   fresnel_radius(x, L, period, u_ref)          → R_f (km)
#   hffk_weight(r, R_f)                          → K (1/km², Q=1)
#   fresnel_sample_offsets(R_f, n_rings, n_az)   → [(Δe, Δn, w), ...]
#   offset_to_geographic(lon, lat, az, el, Δe, Δn) → (lon, lat)
#
# 變數：
#   x   : 沿射線距離 (km)
#   r   : 垂直射線距離 (km)
#   L   : 射線總長度 (km)
#   T   : dominant period (s)
#   u   : 1D 參考 slowness (s/km)  [= 1/v_ref]
# =============================================================================

# -----------------------------------------------------------------------------
# Fresnel zone 半徑
# -----------------------------------------------------------------------------
"""
    fresnel_radius(x, L, T, u_ref) -> R_f (km)

計算射線上距源點 x km 處的第一 Fresnel zone 半徑。
  x      : 沿射線距離 (km), 0 <= x <= L
  L      : 射線總長度 (km)
  T      : dominant period (s)
  u_ref  : 該點的 1D 參考慢度 (s/km), = 1/v_ref
"""
function fresnel_radius(x::T, L::T, period::T, u_ref::T) where {T <: AbstractFloat}
    (x <= 0.0 || x >= L) && return zero(T)   # 括號必要：&&優先序高於||
    return sqrt(period * x * (L - x) / (L * u_ref))
end

# -----------------------------------------------------------------------------
# HFFK 權重
# -----------------------------------------------------------------------------
"""
    hffk_weight(r, R_f) -> weight (1/km²)

在某深度 x 的橫截面上，距射線 r km 處的 HFFK 敏感度密度。
歸零條件: r > R_f（Fresnel zone 外敏感度為零）
Donut hole: K(r=0) = sin(0)/(πR_f²) = 0，射線本身無一階敏感度（Fermat）
正規化: Q=1，∫₀^{R_f} K·2πr dr = 2/π ≈ 0.637（嚴格歸一需 Q=π/2）
        caller 用加權平均（÷ sum of weights），Q 自動消去
"""
function hffk_weight(r::T, R_f::T) where {T <: AbstractFloat}
    R_f <= 0.0 && return zero(T)
    r >= R_f   && return zero(T)
    return sin(π * r^2 / R_f^2) / (π * R_f^2)
end

# -----------------------------------------------------------------------------
# 給定射線點，生成 Fresnel zone 內的取樣點
# -----------------------------------------------------------------------------
"""
    fresnel_sample_offsets(R_f, n_rings, n_azimuth) -> offsets::Vector{Tuple{Float64,Float64,Float64}}

在 Fresnel zone 截面（垂直射線的平面）產生取樣點的二維偏移量。
傳回: [(Δeast_km, Δnorth_km, weight), ...]
  中心點 (r=0) 有 weight = 0（doughnut hole 近似中心亦包含，但 HFFK 在 r→0 時 sin→0 自然趨零）
  n_rings   : 徑向環數（預設 3）
  n_azimuth : 每環的方位取樣數（預設 8）
"""
function fresnel_sample_offsets(R_f::T, n_rings::Int=3, n_azimuth::Int=8) where {T <: AbstractFloat}
    offsets = Tuple{T, T, T}[]

    # 中心點（幾何射線本身，HFFK sin(0)=0，weight→0 但仍保留以保持向後相容）
    push!(offsets, (zero(T), zero(T), zero(T)))

    # 同心環
    for ir in 1:n_rings
        r = R_f * ir / (n_rings + 0.5)   # 不取到 R_f 邊界（sin→0）
        w = hffk_weight(r, R_f)
        for ia in 1:n_azimuth
            az = 2π * (ia - 1) / n_azimuth
            Δe = r * cos(az)
            Δn = r * sin(az)
            push!(offsets, (Δe, Δn, w))
        end
    end

    return offsets
end

# -----------------------------------------------------------------------------
# 輔助：把垂直射線的偏移量轉換到地理座標
# -----------------------------------------------------------------------------
"""
    offset_to_geographic(lon_ray, lat_ray, azimuth_deg, elevation_deg, Δe_km, Δn_km)
    -> (lon, lat)

把射線點 (lon_ray, lat_ray) 附近、垂直射線方向的偏移量 (Δe, Δn) km
轉換為地理座標。azimuth_deg = 射線方位角（從北順時針），elevation_deg = 射線仰角。
"""
function offset_to_geographic(lon_ray::T, lat_ray::T,
                               azimuth_deg::T, elevation_deg::T,
                               Δe_km::T, Δn_km::T) where {T <: AbstractFloat}
    # 地球半徑 (km)
    R_earth = 6371.0
    km_per_deg = π * R_earth / 180.0

    # 射線水平方位角決定「垂直射線平面」的方向
    # 垂直射線的兩個正交方向：
    #   ê₁ = 水平面內垂直射線方位 (perpendicular-horizontal)
    #   ê₂ = 垂直方向（但此處我們只做水平方向偏移近似，對幾乎垂直入射的 SKS 合理）
    az_rad = deg2rad(azimuth_deg)

    # 水平平面內垂直射線的方向（方位 + 90°）
    perp_az = az_rad + π/2.0
    Δlon_km = Δe_km * sin(perp_az) + Δn_km * cos(perp_az)
    Δlat_km = Δe_km * cos(perp_az) - Δn_km * sin(perp_az)

    lon_new = lon_ray + Δlon_km / (km_per_deg * cos(deg2rad(lat_ray)))
    lat_new = lat_ray + Δlat_km / km_per_deg

    return (lon_new, lat_new)
end
