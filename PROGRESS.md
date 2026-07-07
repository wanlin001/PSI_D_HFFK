# PSI_D_HFFK 實作進度記錄

建立日期：2026-06-24  
作者：Wanlin  
目標：在 PSI_D 中實作 HFFK（Heuristic Finite-Frequency slowness Kernel），  
用於 Taiwan double subduction SKS splitting forward modeling。

---

## 專案概覽

| 項目 | 路徑 |
|------|------|
| 原版 PSI_D（本地） | `/Users/wanlin/Documents/ASPECT/PSI_D/` |
| HFFK 修改版（本地） | `/Users/wanlin/Documents/ASPECT/PSI_D_HFFK/` |
| 原版 PSI_D（HPC） | `/home/wl/software/ECOMAN2.0-seismology.PSI_D/` |
| HFFK 修改版（HPC） | `/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/` |
| 驗證工作目錄（HPC） | `/home/wl/work/ASPECT/0624_HFFK_validation/` |
| Taiwan SKS 工作目錄（HPC） | `/home/wl/work/ASPECT/0622_PSI_D_TaiwanRyukyuManila/` |

**原則：原版 PSI_D 任何檔案均未修改。**

---

## 方法理論：HFFK

來源：VanderBeek & Faccenda (2021, GJI) eq.4  
延伸自：Schmandt & Humphreys (2010, Nature)

**Kernel 公式：**

```
K(x, r) = Q/(π R_f²(x)) · sin(π r²/R_f²(x))
```

**Fresnel zone 半徑：**

```
R_f(x) = sqrt(T · x(L-x) / (L · u(x)))
```

變數說明：
- `x`：沿射線距離 (km)
- `r`：垂直射線距離 (km)  
- `L`：射線總長度 (km)
- `T`：dominant period (s)
- `u(x)`：1D 參考 slowness (s/km) = 1/v_ref

射線路徑仍用 TauP（1D reference model），HFFK 只替換 kernel 的積分方式。

---

## 修改的檔案（僅 PSI_D_HFFK copy）

### 新增

#### `src/psi_fresnel_kernel.jl`
HFFK 數學核心：
- `fresnel_radius(x, L, period, u_ref)` → R_f (km)
- `hffk_weight(r, R_f)` → 歸一化 HFFK 敏感度 (1/km²)
- `fresnel_sample_offsets(R_f, n_rings, n_azimuth)` → [(Δe, Δn, weight), ...]
- `offset_to_geographic(lon, lat, azimuth, elevation, Δe, Δn)` → (lon_new, lat_new)

#### `validation/psi_parameters_hffk_forward.toml`
HFFK 的 TOML 設定格式（SinkingSlab 驗證用）：
```toml
[Model.Methods.TauP]
reference_model = "psi_input/ak135_no_crust.tvel"
DL = 10.0

[Model.Methods.FiniteFrequency]
dominant_period = 15.0
n_rings = 3
n_azimuth = 8
```

#### `validation/job_ray_validation.sh` / `job_hffk_validation.sh`
HPC sbatch 腳本，partition `8358`。

#### `validation/compare_ray_vs_hffk.py`
比較 ray theory 與 HFFK 輸出的 splitting intensity。

### 修改（相對於原版）

#### `src/PSI_D.jl`
加入：
```julia
include("psi_fresnel_kernel.jl")   # HFFK: fresnel_radius, hffk_weight, fresnel_sample_offsets
```

#### `src/psi_forward.jl`

1. **`ForwardFiniteFrequency` struct 加入參數：**
```julia
struct ForwardFiniteFrequency <: ForwardMethod
    dominant_period::Float64
    n_rings::Int
    n_azimuth::Int
end
```

2. **dispatch block 加入 HFFK 分支：**（約 line 416–430）

3. **新增 `psi_forward_hffk(Observation, Model)`：**
   - 呼叫 `return_kernel` 取得 TauP 射線路徑
   - 沿射線每個取樣點計算 R_f
   - 在 Fresnel zone 截面生成 n_rings × n_azimuth 偏移取樣點
   - 每個偏移點用 trilinear 插值模型，計算 splitting intensity
   - HFFK weight 加權累積，最後積分

4. **新增 `interpolate_si_at_point`：**
   - 任意地理座標的 splitting intensity 計算
   - 使用正確 Thomsen 參數（ratio_ϵ, ratio_η, ratio_γ）
   - 與 `kernel_splitting_intensity` 相同公式：`SI = 0.5*(1/vqs2 - 1/vqs1)*sin(2ζ)`

#### `src/psi_buildinputs.jl`

1. **TOML 解析加入 HFFK 參數：**
```julia
elseif haskey(D["Model"]["Methods"], "FiniteFrequency")
    ff = D["Model"]["Methods"]["FiniteFrequency"]
    FwdInstance = ForwardFiniteFrequency(
        get(ff, "dominant_period", 12.0),
        get(ff, "n_rings", 3),
        get(ff, "n_azimuth", 8)
    )
    FwdType = ForwardFiniteFrequency
```

2. **`build_observations` 加入 `FwdInstance` 參數：**
   避免對 `ForwardFiniteFrequency`（有必填欄位）呼叫 `FwdType()` 失敗。

3. **HFFK 允許 TauP + FiniteFrequency 同時存在：**
   解除 `length(Methods) > 1` 的限制（HFFK 需要 TauP 提供射線路徑）。

---

## 第三個 session 修復記錄（2026-06-26）

### Bug 1：`UndefVarError: SeismicDijkstra` (Job 33582)
- **原因：** `wrapper_seismic_dijkstra.jl` 被 include，`try-catch using SeismicDijkstra` 雖然不報錯，但 wrapper 裡直接呼叫 `SeismicDijkstra.xxx` 仍失敗
- **修法：** `src/PSI_D.jl` 中 comment 掉該 include

### Bug 2：Taiwan 近對蹠距離問題 (Job 33583)
- **原因：** 141 對觀測中，南美 → 台灣距離 ~160-178°，標為 "SKS" 但 SKS 視窗只到 ~130°；TauP 找不到 SKS → 空陣列 → `maximum()` crash
- **修法：** 用 Lin 2014a 條件（88-115°、depth>50km）過濾，產生 `SplittingParameters_filtered.dat`（14 對）

### Bug 3：`cannot parse "" as Float64` (Job 33586)
- **原因：** psitomo0020.dat 第 3 行有 trailing comma，`split(",")` 後最後一元素是 `""`；HFFK 版 `psi_buildinputs.jl` 的 `dz_ext = length(line)>3 ? parse(T, line[4]) : 0.0` 試圖 parse 空字串
- **修法：** 改為 `(length(line) > 3 && !isempty(strip(line[4]))) ? parse(T, strip(line[4])) : 0.0`

### Bug 4：`matrix contains Infs or NaNs` (Job 33587)
- **原因：** `psitomo0050.dat` 有 ~610 萬個 NaN（~15% 格點在 ASPECT domain 外 VIZTOMO 無輸出）
- **修法：** 純 Python 腳本（`/tmp/fill_nan_pure.py`）對每個 Voigt 分量用欄位中位數填充，產生 `psitomo0050_filled.dat`
- **狀態：** 背景執行中（2026-06-26）

---

## 驗證計畫

### 用 SinkingSlab 範例驗證

SinkingSlab 是 PSI_D 官方範例，有：
- 解析式的 HexagonalVectoralVelocity 模型（含 sinking anisotropic slab）
- 預先計算的 DUMMY 觀測資料
- 已知 ray theory 輸出可對比

**驗證步驟：**

| 步驟 | Job ID | 狀態 | 說明 |
|------|--------|------|------|
| A. VIZTOMO bug fix | 33565 | ✅ 完成 | psitomo0050.dat 標頭正確（lon=122.7, lat=24.76） |
| B. Ray theory forward (SinkingSlab) | 33566 | ✅ 完成 | 原版 PSI_D，15000 筆 SI |
| C. HFFK forward (SinkingSlab) | 33591 | 🔄 在跑 | PSI_D_HFFK，psi_buildinputs dz_ext 修復後提交 |
| D. 比較 Ray vs HFFK | — | 待 C 完成 | `python validation/compare_ray_vs_hffk.py` |
| E. Taiwan ray theory | 33587 | ❌ 失敗→重跑 | NaN 問題，等 psitomo0050_filled.dat 完成 |

**預期結果（根據 VanderBeek 2021 推斷）：**
- HFFK 與 ray theory 趨勢相同（相關係數 > 0.95）
- slab 邊界附近 HFFK 更平滑（Fresnel zone 平均效應）
- 不應有數量級差異；若出現 NaN 或極大值 → 有 bug

---

## VIZTOMO / psitomo 進度

### 問題
VIZTOMO Fortran 中 `rad2deg` 未初始化（`comvar` module 宣告但未賦值）→ Fortran 零初始化 → 所有輸出 lon=0, lat=90。

### 修復
在 `/home/wl/software/ECOMAN2.0-geodynamics_org/VIZTOMO/VIZTOMO.f90` 第 162 行加入：
```fortran
rad2deg = 180.0d0/pi
```

### 狀態
Job 33565（VIZTOMO）仍在佇列 PD。  
執行完成後，需確認 `psitomo0050.dat` 標頭第一行：
```
 0.63710E+04,  1.22700E+02,  2.47600E+01,  0.00000E+00,   ← 正確（lon≈122.7, lat≈24.76）
```
（現有輸出：`0.63710E+04, 0.00000E+00, 0.90000E+02, 0.00000E+00` — 仍是舊的壞資料）

---

## Taiwan SKS Pipeline 整體狀態

```
ASPECT → D-Rex_M → VIZTOMO → psitomo0050.dat → PSI_D → SKS SI predictions
                      ↑                            ↑
               job 33565 (PD)              job 33566 (PD) [ray theory]
                  FIXING                           比較基準
```

待辦事項（依序）：
1. 等 job 33565 完成 → 確認 psitomo0050.dat 標頭 → 若正確繼續
2. 更新 Taiwan TOML 的 `theModel` 路徑（psitomo0050.dat）
3. 提交 PSI_D ray theory forward job（Taiwan real data）
4. 等 job 33566 (SinkingSlab) 完成 → 確認輸出 → 提交 HFFK job (33567)
5. 比較 Ray vs HFFK（validation/compare_ray_vs_hffk.py）
6. 若 HFFK 通過 → 對 Taiwan 資料跑 HFFK forward
7. 套用 event 過濾（88°–115°，depth > 50 km）：`python validation/filter_events.py`

---

## 引用

- VanderBeek, B. P., & Faccenda, M. (2021). Imaging upper mantle anisotropy with teleseismic P-wave delays. *GJI*, 225(3), 2097–2119. https://doi.org/10.1093/gji/ggab081
- Schmandt, B., & Humphreys, E. (2010). Complex subduction and small-scale convection revealed by body-wave tomography. *Nature Geoscience*, 3, 55–59.
- Dahlen, F. A., Hung, S.-H., & Nolet, G. (2000). Fréchet kernels for finite-frequency traveltimes. *GJI*, 141, 157–174.
