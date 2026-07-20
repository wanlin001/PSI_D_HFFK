# PSI_D_HFFK 驗證文件 — Tier 2 方法學收斂測試

**專案**：ECOMAN2.0-seismology.PSI_D_HFFK  
**作者**：Wanlin Hu  
**最後更新**：2026-07-16  
**平台**：HPC `wl`（SLURM partition `8358`）

> 本文件供 GitHub 公開檢視，說明 HFFK（Heuristic Finite-Frequency Kernel）實作的方法學收斂測試設計、執行方式與目前結果。  
> 完整驗證框架分三層：Tier 1（單元正確性）→ **Tier 2（方法學收斂）** → Tier 3（應用有效性）。

---

## 1. 背景：Ray vs HFFK

| 方法 | 物理意義 | PSI_D 實作 |
|------|----------|------------|
| **Ray SI** | 無限高頻極限；敏感度集中在幾何射線（δ-function） | `SplittingIntensity` + `ForwardTauP`：TauP 射線 + 沿射線 SI kernel 積分 |
| **HFFK SI** | 有限頻率 banana-donut 敏感核 | `SplittingIntensity` + `ForwardFiniteFrequency`：同射線 + Fresnel zone 加權平均 |

**HFFK kernel**（VanderBeek & Faccenda 2021, GJI, eq.4）：

```
K(x, r) = Q / (π · R_f²) · sin(π · r² / R_f²)

R_f(x) = √( T · x · (L − x) / (L · u(x)) )
```

- `T`：dominant period (s)  
- `r`：垂直射線距離 (km)  
- `L`：射線總長度 (km)  
- `u(x)`：參考 S-wave slowness (s/km)

**比較原則**：Ray 與 HFFK 必須在**同一 observable**（`SplittingIntensity`）上比較。  
不可使用 `SplittingParameters`（δt, φ）事後轉 SI——PSI_D 內部 φ 與 SI 符號慣例不一致，會造成假反相。

---

## 2. 驗證框架總覽

```
Tier 1 — 單元正確性     「算得對不對？」     synthetic benchmark（均勻 / 雙層 / 橫向）
Tier 2 — 方法學收斂   「離散近似夠不夠準？」  period / n_rings / n_azimuth / 格網 / 極限
Tier 3 — 應用有效性   「科學上有沒有說服力？」  SinkingSlab / Kuo2018 / Taiwan SKS
```

**Tier 2 通過標準**（建議）：

- 收斂測試：參數加倍後 RMS 殘差下降 ≥ 10×，或已 < 0.01 s  
- 極限測試：HFFK T=4s 與 Ray SI 的 RMS < 0.05 s（橫向邊界模型）  
- 格網測試：均勻場粗/細格網 RMS < 1e-10 s；橫向邊界模型需記錄格網依賴性（見 T2-4）

---

## 3. Tier 2 測試清單

### T2-1. 頻率極限：HFFK(T) → Ray SI

| 項目 | 內容 |
|------|------|
| **目的** | 驗證 T→小 時 HFFK 趨近無限高頻 Ray 極限 |
| **模型** | `bench_lateral_B`（橫向 φ 邊界，period 效應最明顯） |
| **比較** | `HFFK(T)` vs `Ray SI (inf. freq.)`，T = 4, 8, 16, 20, 25, 33, 50 s |
| **通過標準** | T=4s：corr > 0.99，RMS < 0.05 s；T↑ 時 RMS 單調增大（Fresnel 平均效應） |
| **腳本** | `validation/job_bench.sh` + `validation/plot_benchmark_groups.py` |
| **狀態** | ✅ 已完成（Job 33970, 2026-07-15） |

**結果（`bench_lateral_B`，144 obs，按 BAZ 分組平均）**：

| Period T (s) | corr(HFFK, Ray SI) | RMS (s) |
|:------------:|:------------------:|:-------:|
| 4 | 0.99980 | 0.023 |
| 8 | 0.99924 | 0.048 |
| 16 | 0.99763 | 0.092 |
| 20 | 0.99678 | 0.111 |
| 25 | 0.99579 | 0.134 |
| 33 | 0.99441 | 0.165 |
| 50 | 0.99216 | 0.222 |

**解讀**：T=4s 已接近 Ray SI；T=50s 偏差 ~0.22 s 為 Fresnel zone 跨邊界平均的**預期物理效應**，非 bug。

---

### T2-2. 徑向離散收斂：n_rings

| 項目 | 內容 |
|------|------|
| **目的** | 驗證 Fresnel zone 徑向環數離散取樣是否收斂 |
| **模型** | `bench_lateral_B` |
| **設定** | T = 25 s，n_azimuth = 8（固定），n_rings = 1, 3, 5, 8 |
| **參考** | n_rings = 8 |
| **通過標準** | n_rings = 3 時 RMS < 0.01 s（生產預設值） |
| **腳本** | `validation/job_bench_nrings.sh` + `validation/plot_nrings_convergence.py` |
| **狀態** | ✅ 已完成（Job 33971, 2026-07-15） |

**結果（vs n_rings = 8）**：

| n_rings | RMS (s) | max \|diff\| (s) | corr |
|:-------:|:-------:|:----------------:|:----:|
| 1 | 0.00719 | 0.0188 | 0.99994 |
| **3**（**預設**） | **0.00097** | **0.0024** | **1.00000** |
| 5 | 0.00031 | 0.0008 | 1.00000 |
| 8 | 0（參考） | 0 | 1.00000 |

**結論**：**n_rings = 3 已足夠**（RMS < 0.001 s）。生產環境維持預設即可；高精度需求可改 5。

**輸出圖**：`/lfs/wl/bench_psi/bench_output/nrings_conv/nrings_convergence.png`

---

### T2-3. 方位離散收斂：n_azimuth

| 項目 | 內容 |
|------|------|
| **目的** | 驗證 Fresnel zone 方位角取樣是否收斂 |
| **模型** | `bench_lateral_B` |
| **設定** | T = 25 s，n_rings = 3（固定），n_azimuth = 4, 8, 16, 32 |
| **參考** | n_azimuth = 32 |
| **通過標準** | n_azimuth = 8 時 RMS < 0.01 s |
| **腳本** | 待建：`validation/job_bench_nazimuth.sh`（仿 `job_bench_nrings.sh`） |
| **狀態** | ⏳ 待執行 |

**執行方式（草案）**：

```bash
cd /lfs/wl/bench_psi
# 建立 job_bench_nazimuth.sh 後：
sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench_nazimuth.sh
```

---

### T2-4. 格網解析度收斂（模型空間）

| 項目 | 內容 |
|------|------|
| **目的** | 區分「均勻場格網無影響」vs「尖銳橫向邊界對 HFFK 的格網依賴」 |
| **測試 A** | 均勻單層：`bench_1L_A`（5×5×13, 2°）vs `bench_1L_B`（9×9×21, 1°） |
| **測試 B** | 橫向邊界：`bench_lateral_B`（9×9×21, 1°）vs `bench_lateral_C`（17×17×21, 0.5°） |
| **通過標準** | 測試 A：RMS < 1e-10 s；測試 B：記錄 RMS 並在論文中說明格網需求 |
| **腳本** | `validation/job_bench.sh`（已含全部模型） |
| **狀態** | ✅ 已完成 |

**結果（HFFK T = 8 s）**：

| 比較 | 格網 | RMS (s) | 解讀 |
|------|------|:-------:|------|
| 1L_A vs 1L_B | 5×5 vs 9×9 | **1.8×10⁻¹⁶** | ✅ 均勻場完全收斂 |
| lateral_B vs lateral_C | 9×9 vs 17×17 | **0.258** | ⚠️ 尖銳邊界 + HFFK Fresnel 取樣對格網敏感 |

**重要說明**：Col A benchmark 圖下排紅色大振幅殘差是 **`lateral_C − lateral_B`**，不是 `1L_B − 1L_A`。  
均勻場 resolution 測試已通過；橫向邊界模型建議使用 ≥ 1° 格網，或接受 ~0.26 s 的格網誤差。

---

### T2-5. 均勻 / 垂直分層場的 period 獨立性

| 項目 | 內容 |
|------|------|
| **目的** | 確認均勻或純垂直分層模型上，HFFK 不應有 period 效應 |
| **模型** | `bench_1L_B`（單層）、`bench_2L_B`（雙層垂直分層） |
| **比較** | HFFK T = 4–50 s 之間互相比較，及 vs Ray SI |
| **通過標準** | 所有 T 的 SI 差異 < 1e-10 s |
| **狀態** | ✅ 已完成 |

**結果**：

| 模型 | HFFK T=4s vs T=50s max diff | Ray SI vs HFFK（所有 T） |
|------|:---------------------------:|:------------------------:|
| bench_1L_B | ~1.3×10⁻¹⁵ | corr = 1.000, RMS = 0 |
| bench_2L_B | ~2.9×10⁻¹⁶ | corr = 1.000, RMS = 0 |

---

### T2-6. SinkingSlab 官方範例（跨模型驗證）

| 項目 | 內容 |
|------|------|
| **目的** | 在 PSI_D 官方各向異性幾何（HexagonalVectoralVelocity）上驗證 HFFK |
| **模型** | SinkingSlab（`psi_output/SinkingSlab_*`） |
| **比較** | Ray vs HFFK T = 3, 5, 8, 15, 25 s |
| **預期** | 趨勢一致但 HFFK 振幅壓縮（Fresnel 平均板內/板外）；corr ~0.3–0.5 為物理效應 |
| **腳本** | `validation/job_ray_validation.sh`, `validation/job_hffk_validation.sh`, `validation/plot_si_3d.py` |
| **狀態** | ✅ 歷史結果已有（wl2 Jobs 44189, 44255–44260）；待整理進本文件附錄 |

---

### T2-7. HFFK SP 路徑（SplittingParameters + tf_hffk_sp）

| 項目 | 內容 |
|------|------|
| **目的** | 驗證 HFFK 對 `SplittingParameters` 的 Fresnel Cij 平滑路徑 |
| **機制** | `tf_hffk_sp = true`：Fresnel zone 加權平均 Cij → waveform propagation |
| **模型** | Taiwan SKS 或 SinkingSlab |
| **狀態** | ⏳ 待系統測試（Taiwan 14 對 SKS） |

---

## 4. 執行指令速查

### 4.1 環境

```bash
ssh wl
cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
module load anaconda/2022.05   # Python 繪圖需要
```

### 4.2 Tier 1 + Tier 2 主 benchmark（含 Ray SI）

```bash
# 生成合成模型（首次或模型更新時）
python3 validation/gen_benchmark_psitomo.py

# 提交（必須從 /lfs 目錄）
mkdir -p /lfs/wl/bench_psi && cd /lfs/wl/bench_psi
sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench.sh
```

**輸出**：`/lfs/wl/bench_psi/bench_output/`（54 runs：6 模型 × 9 cases）

### 4.3 n_rings 收斂（T2-2）

```bash
cd /lfs/wl/bench_psi
sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench_nrings.sh
```

### 4.4 繪圖

```bash
python3 validation/plot_benchmark_groups.py \
  --bench-out /lfs/wl/bench_psi/bench_output \
  --out /lfs/wl/bench_psi/bench_output/benchmark_groups

python3 validation/plot_nrings_convergence.py \
  --conv-dir /lfs/wl/bench_psi/bench_output/nrings_conv
```

---

## 5. 輸出檔案對照

| 路徑 | 內容 |
|------|------|
| `bench_output/bench_{model}_ray_si/` | Ray SI（無限高頻） |
| `bench_output/bench_{model}_hffk_T{N}s/` | HFFK SI，period = N s |
| `bench_output/benchmark_groups.png` | 四欄比較圖 |
| `bench_output/benchmark_groups_r_vs_period.png` | corr(HFFK, Ray SI) vs T |
| `bench_output/nrings_conv/` | n_rings 收斂輸出 |
| `bench_output/nrings_conv/nrings_convergence.png` | n_rings 收斂圖 |

**模型代號**：

| 代號 | 說明 | 格網 |
|------|------|------|
| `bench_1L_A` | 單層均勻 φ=45° | 5×5×13 |
| `bench_1L_B` | 單層均勻 φ=45° | 9×9×21 |
| `bench_2L_A/B` | 雙層垂直 φ=0°/45° | 粗/細 |
| `bench_lateral_B` | 橫向邊界 φ_L=0°/φ_R=90° | 9×9×21 |
| `bench_lateral_C` | 同上（格網收斂用） | 17×17×21 |

---

## 6. Tier 2 完成度摘要

| 編號 | 測試 | 狀態 | 結論 |
|:----:|------|:----:|------|
| T2-1 | HFFK(T) → Ray SI 極限 | ✅ | T=4s RMS=0.023s；T=50s RMS=0.222s（物理正確） |
| T2-2 | n_rings 收斂 | ✅ | n_rings=3 足夠（RMS=0.001s） |
| T2-3 | n_azimuth 收斂 | ⏳ | 待執行 |
| T2-4 | 格網解析度 | ✅ | 均勻場通過；橫向邊界有格網依賴（需 ≥1° 或接受 ~0.26s） |
| T2-5 | 均勻/2L period 獨立 | ✅ | 所有 T 完全相同（機器精度） |
| T2-6 | SinkingSlab | ✅ | 歷史結果可用，待整理 |
| T2-7 | HFFK SP 路徑 | ⏳ | 待 Taiwan SKS 測試 |

**Tier 2 整體判定**：核心收斂測試（T2-1, T2-2, T2-4, T2-5）已通過。  
**待補**：T2-3（n_azimuth）、T2-7（HFFK SP），方可宣稱「數值方法完全收斂」。

---

## 7. 相關腳本與設定檔

```
validation/
├── gen_benchmark_psitomo.py          # 合成模型生成
├── job_bench.sh                      # 主 benchmark（Ray SI + HFFK T=4–50s）
├── job_bench_nrings.sh               # n_rings 收斂
├── plot_benchmark_groups.py          # 四欄比較圖
├── plot_nrings_convergence.py        # n_rings 收斂圖
├── psi_config_template_bench.toml      # HFFK / Ray SP 設定
├── psi_config_template_bench_ray_si.toml  # Ray SI 設定
├── bench_models/                     # 合成 psitomo .dat
└── bench_psi_input/                  # Sources, Receivers, DUMMY_SI/SP

src/
├── psi_fresnel_kernel.jl             # HFFK 數學核心
└── psi_forward.jl                    # psi_forward_hffk()
```

---

## 8. 引用

- VanderBeek, B. P., & Faccenda, M. (2021). Imaging upper mantle anisotropy with teleseismic P-wave delays. *GJI*, 225(3), 2097–2119. https://doi.org/10.1093/gji/ggab081
- Schmandt, B., & Humphreys, E. (2010). Complex subduction and small-scale convection revealed by body-wave tomography. *Nature Geoscience*, 3, 55–59.
- Faccenda, M., & VanderBeek, B. P. (2023). S-wave splitting intensity imaging. *J. Geodynamics*, 158, 102003.
- Chevrot, S. (2000). Multiplet waveform cross correlation on the western United States. *GJI*, 140, 480–496.

---

## 9. 附錄：Tier 1 與 Tier 3 簡表

### Tier 1 — 單元正確性（已完成）

- 均勻 1L：Ray SI = HFFK（所有 T），corr = 1.0  
- 橫向 lateral_B：HFFK period 效應符合 Fresnel 平均預期  
- 解析解 analytical_si.csv 與 HFFK 一致（SI 符號慣例已修正，2026-07-15）

### Tier 3 — 應用有效性（進行中）

- Kuo2018 真實 3D 模型 period sensitivity map（`scripts/plot_real_period_map.py`）  
- Taiwan SKS：14 對有效觀測（Lin 2014a 條件，88°–115°，depth > 50 km）  
- ASPECT → VIZTOMO → psitomo → PSI_D 完整 pipeline

---

*本文件將隨測試進度更新。GitHub 路徑建議：`docs/validation/Tier2_Methodological_Convergence.md`*
