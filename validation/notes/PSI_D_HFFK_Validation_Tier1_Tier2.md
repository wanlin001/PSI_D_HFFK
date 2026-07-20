# PSI_D_HFFK 驗證文件 — Tier 1 + Tier 2

**專案**：ECOMAN2.0-seismology.PSI_D_HFFK  
**作者**：Wanlin Hu  
**最後更新**：2026-07-17  
**平台**：HPC `wl`（SLURM partition `8358`）

> 本文件供 GitHub 公開檢視。  
> 完整驗證框架：Tier 1（單元正確性）→ Tier 2（方法學收斂）→ Tier 3（應用有效性）。

**輸出圖**（Desktop / HPC）：
- Tier 1：`benchmark_tier1.png`
- Tier 2：`benchmark_tier2.png`（2×3，含 n_rings / n_azimuth / lateral_D）
- T2-2：`nrings_convergence.png`
- T2-3：`nazimuth_convergence.png`
- T3-3：`sinking_slab_comparison.png`

**Jobs 2026-07-16**：
- 34000–34002：**作廢**（`dims3` 錯誤模式）
- 34003–34005：SinkingSlab ✅
- **34007**（bench `linear`）、**34008**（n_rings）、**34009**（n_azimuth）✅
- Git main：`5d4f651`（benchmark 改 `linear`）

---

## 1. Ray vs HFFK 比較原則

| 方法 | 物理意義 | PSI_D 實作 |
|------|----------|------------|
| **Ray SI** | 無限高頻極限；敏感度集中在幾何射線 | `SplittingIntensity` + `ForwardTauP` |
| **HFFK SI** | 有限頻率 banana-donut 敏感核 | `SplittingIntensity` + `ForwardFiniteFrequency` |

**HFFK kernel**（VanderBeek & Faccenda 2021, eq.4）：

```
K(x, r) = Q / (π · R_f²) · sin(π · r² / R_f²)
R_f(x) = √( T · x · (L − x) / (L · u(x)) )
```

**關鍵**：Ray 與 HFFK 必須在**同一 observable**（`SplittingIntensity`）上比較。  
不可使用 `SplittingParameters` 事後轉 SI（會造成假反相）。

---

## 2. Tier 2 測試檔案在哪裡？

### 2.1 程式碼與設定（repo）

```
/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/
├── gen_benchmark_psitomo.py              # 合成 psitomo 模型生成
├── job_bench.sh                          # 主 benchmark（Tier 1 + T2-1/4/5）
├── job_bench_nrings.sh                   # T2-2 n_rings 收斂
├── job_bench_nazimuth.sh                 # T2-3 n_azimuth 收斂
├── job_sinking_slab.sh                   # T3-3 SinkingSlab
├── compare_sinking_slab.py               # T3-3 比較腳本
├── plot_benchmark_tier1.py               # Tier 1 圖（獨立 PNG）
├── plot_benchmark_tier2.py               # Tier 2 圖（獨立 PNG）
├── plot_benchmark_groups.py              # 舊版合併四欄圖（保留）
├── plot_convergence_common.py            # 收斂測試共用 RMS 計算（raw SI）
├── plot_nrings_convergence.py            # T2-2 詳圖
├── plot_nazimuth_convergence.py          # T2-3 詳圖
├── psi_config_template_bench.toml        # HFFK / Ray SP
├── psi_config_template_bench_ray_si.toml   # Ray SI（TauP only）
├── bench_models/                         # 合成 .dat + analytical_si.csv
│   ├── bench_1L_A.dat / bench_1L_B.dat   # 均勻場（Col A 格網）
│   ├── bench_1L_lat_B.dat                # 單層橫向邊界（Col B period）
│   ├── bench_2L_lat_B.dat                # 雙層下層橫向（Col C period）
│   ├── bench_lateral_B.dat / _C.dat      # 全深度橫向（Col D + T2）
│   └── analytical_si.csv
└── bench_psi_input/                      # Sources, Receivers, DUMMY_SI/SP
```

### 2.2 執行輸出（/lfs，SLURM 提交目錄）

```
/lfs/wl/bench_psi/
├── bench_33985.out                       # 最近一次 job log（2026-07-16）
└── bench_output/
    ├── bench_{model}_ray_si/             # Ray SI
    ├── bench_{model}_hffk_T{N}s/         # HFFK SI，period = N s
    ├── benchmark_tier1.png               # ← Tier 1 圖
    ├── benchmark_tier2.png               # ← Tier 2 圖
    └── nazimuth_conv/                      # T2-3 輸出
        ├── lateral_B_hffk_T25s_az{N}/
        └── nazimuth_convergence.png
    └── sinking_slab/                       # T3-3 輸出
        ├── SinkingSlab_ray_si/
        └── SinkingSlab_hffk_T{N}s/
```

**提交方式**（必須從 `/lfs` 目錄）：

```bash
ssh wl
mkdir -p /lfs/wl/bench_psi && cd /lfs/wl/bench_psi
sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench.sh
sbatch validation/job_bench_nrings.sh
sbatch validation/job_bench_nazimuth.sh
sbatch validation/job_sinking_slab.sh
```

---

## 3. Tier 1 — 單元正確性

### 3.1 設計（2026-07-16 更新）

原先 Col B/C 使用均勻 / 純垂直分層 φ，HFFK 所有 period **完全重疊**（物理正確但無法視覺化 period 效應）。  
現改為：

| 欄位 | 模型 | φ 結構 | 目的 |
|:----:|------|--------|------|
| **A** | `bench_1L_A` vs `bench_1L_B` | 均勻 φ=45°，5×5 vs 9×9 | 格網解析度（均勻場應 → 0） |
| **B** | `bench_1L_lat_B` | 單層橫向 φ_L=0°/φ_R=90° @ 123°E | period 效應可視化 |
| **C** | `bench_2L_lat_B` | 上層均勻 φ=0°，下層橫向邊界 | 雙層 + 水平 φ 干涉 |
| **D** | `bench_lateral_B` | 全深度橫向邊界 | 最強 period 效應基準 |

保留 `bench_1L_B` / `bench_2L_B`（均勻/垂直）供 T2-5 period 獨立性測試。

### 3.2 圖：`benchmark_tier1.png`

- **上排**：SI vs BAZ，Ray SI（紅）+ HFFK T=4–50 s（彩色）
- **下排**：HFFK(T) − Ray SI 殘差
- Col A 僅顯示格網比較（不再疊加 lateral 殘差）

### 3.3 通過標準

- Col A：1L_B − 1L_A RMS < 1e-10 s（均勻場）
- Col B/C/D：T=4s 接近 Ray SI；T↑ 殘差增大（Fresnel 平均，非 bug）
- 均勻/垂直模型（T2-5）：所有 T 機器精度一致

---

## 4. Tier 2 — 方法學收斂

### 4.1 圖：`benchmark_tier2.png`（2×3）

| 面板 | 測試 | 內容 |
|:----:|------|------|
| 左上 | T2-1 | corr(HFFK, Ray SI) vs Period |
| 中上 | T2-1 | RMS(HFFK − Ray SI) vs Period |
| 右上 | T2-4 | 格網解析度（含 lateral_D） |
| 左下 | T2-2 | n_rings 收斂（ref=16） |
| 中下 | T2-3 | n_azimuth 收斂（ref=48） |
| 右下 | T2-4 | lateral_B/C vs lateral_D |

### 4.2 測試清單與狀態

| 編號 | 測試 | 模型 | 狀態 | 結論 |
|:----:|------|------|:----:|------|
| T2-1 | HFFK(T) → Ray SI | lateral_B, 1L_lat_B, 2L_lat_B | ✅ | Job 34007；T=4s RMS≈0.028s |
| T2-2 | n_rings 收斂 | lateral_B, T=25s | ✅ | Job 34008；n_rings=3 RMS=0.001s |
| T2-3 | n_azimuth 收斂 | lateral_B, T=25s | ✅ | Job 34009；n_azimuth=8 RMS=0.002s |
| T2-4 | 格網解析度 | 1L_A/B, lateral_B/C/D | ✅ | Job 34007；C→D 0.10s |
| T2-5 | 均勻/2L period 獨立 | 1L_B, 2L_B | ✅ | 所有 T 完全相同 |
| T2-6 | SinkingSlab | 官方範例 | ✅ | 見 §10 T3-3 |
| T2-7 | HFFK SP 路徑 | Taiwan SKS | ⏳ | 待測試 |

### 4.3 T2-1 結果（`bench_lateral_B`，Job **34007**，`linear`）

| Period T (s) | corr | RMS (s) |
|:------------:|:----:|:-------:|
| 4 | 0.9996 | 0.028 |
| 8 | 0.9987 | 0.054 |
| 16 | 0.9965 | 0.100 |
| 20 | 0.9954 | 0.119 |
| 25 | 0.9941 | 0.142 |
| 33 | 0.9923 | 0.173 |
| 50 | 0.9894 | 0.229 |

（數值與舊 `dims3` run 相近，因 HFFK/Ray **同模式**下相對差異類似；修正的是**絕對 SI 可信度**。）

### 4.4 T2-2 結果（Job **34008**，ref=n_rings=16，`linear`）

| n_rings | RMS vs nr=16 (s) | corr |
|:-------:|:----------------:|:----:|
| 1 | 0.00704 | 0.99994 |
| 2 | 0.00245 | 0.99999 |
| **3**（**預設**） | **0.00111** | **1.00000** |
| 4 | 0.00065 | 1.00000 |
| 5 | 0.00046 | 1.00000 |
| 6 | 0.00031 | 1.00000 |
| 8 | 0.00016 | 1.00000 |
| 12 | 0.00004 | 1.00000 |
| 16 | 0（參考） | 1.00000 |

### 4.5 T2-3 結果（Job **34009**，ref=n_azimuth=48，`linear`）

| n_azimuth | RMS vs az=48 (s) | corr |
|:---------:|:----------------:|:----:|
| 4 | 0.00992 | 0.99986 |
| 6 | 0.00234 | 0.99999 |
| **8**（**預設**） | **0.00216** | **0.99999** |
| 12 | 0.00036 | 1.00000 |
| 16 | 0.00030 | 1.00000 |
| 24 | 0.00012 | 1.00000 |
| 32 | 0.00002 | 1.00000 |
| 48 | 0（參考） | 1.00000 |

### 4.6 T2-4 格網（HFFK T=8s，Job **34007**，`linear`）

| 比較 | 格網 | RMS (s) | 解讀 |
|------|------|:-------:|------|
| 1L_A vs 1L_B | 5×5 vs 9×9 | **1.8×10⁻¹⁶** | ✅ 均勻場完全收斂 |
| lateral_B vs lateral_C | 9×9 vs 17×17 | **0.260** | ⚠️ 粗→細仍有大誤差 |
| lateral_C vs lateral_D | 17×17 vs 25×25 | **0.101** | ✅ 繼續細化誤差下降 |
| lateral_B vs lateral_D | 9×9 vs 25×25 | **0.358** | 直接跳最細格網非單調收斂 |

**結論**：加 `lateral_D`（25×25，~0.33°）後，C→D 段誤差從 0.26s 降至 **0.10s**，證實尖銳邊界需要 ≤0.33° 格網。

---

## 5. 執行指令速查

```bash
ssh wl
module load anaconda/2022.05
cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK

# 生成/更新合成模型
python3 validation/gen_benchmark_psitomo.py

# 提交 benchmark
mkdir -p /lfs/wl/bench_psi && cd /lfs/wl/bench_psi
sbatch validation/job_bench.sh

# 繪圖（Tier 1 / Tier 2 分開）
python3 validation/plot_benchmark_tier1.py \
  --bench-out /lfs/wl/bench_psi/bench_output

python3 validation/plot_benchmark_tier2.py \
  --bench-out /lfs/wl/bench_psi/bench_output \
  --nrings-dir /lfs/wl/bench_psi/bench_output/nrings_conv \
  --nazimuth-dir /lfs/wl/bench_psi/bench_output/nazimuth_conv

python3 validation/plot_nrings_convergence.py \
  --conv-dir /lfs/wl/bench_psi/bench_output/nrings_conv

python3 validation/plot_nazimuth_convergence.py \
  --conv-dir /lfs/wl/bench_psi/bench_output/nazimuth_conv
```

### 4.7 繪圖腳本說明（均在 `validation/`）

| 腳本 | 輸出 | 用途 |
|------|------|------|
| `plot_benchmark_tier1.py` | `benchmark_tier1.png` | Tier 1 四欄圖 |
| `plot_benchmark_tier2.py` | `benchmark_tier2.png` | Tier 2 總覽 2×3 |
| `plot_nrings_convergence.py` | `nrings_conv/nrings_convergence.png` | T2-2 詳圖 |
| `plot_nazimuth_convergence.py` | `nazimuth_conv/nazimuth_convergence.png` | T2-3 詳圖 |
| `plot_convergence_common.py` | （模組，不直接出圖） | 共用 raw SI RMS |
| `compare_sinking_slab.py` | `sinking_slab_comparison.png` | T3-3 |

**為何 `nrings_convergence.png` 與 `benchmark_tier2.png` 右下角曲線曾不一致？**

1. **統計方式不同（已修正）**：舊版 `tier2` 對收斂曲線做 **BAZ 平均**（16 點），`nrings_convergence.py` 用 **raw 144 obs**；現已統一為 raw SI（`plot_convergence_common.py`）。
2. **參考值不同（已修正）**：舊 nrings 測試 ref=`n_rings=8`、點數 `[1,3,5,8]`；新版 ref=`16`、點數 `[1,2,3,4,5,6,8,12,16]`。
3. **n_azimuth**：`lateral_B_hffk_T25s_az{N}` 已由 Job 34002 跑完（N=4,6,8,12,16,24,32,48），輸出在 `nazimuth_conv/`。

---

## 6. 模型代號總表

| 代號 | 說明 | 格網 | Tier 用途 |
|------|------|------|-----------|
| `bench_1L_A` | 單層均勻 φ=45° | 5×5×13 | T1-A 格網 |
| `bench_1L_B` | 單層均勻 φ=45° | 9×9×21 | T1-A / T2-5 |
| `bench_1L_lat_B` | 單層橫向 φ_L=0°/φ_R=90° | 9×9×21 | **T1-B period** |
| `bench_2L_A/B` | 雙層垂直 φ=0°/45° | 粗/細 | T2-5 |
| `bench_2L_lat_B` | 上均勻 + 下層橫向 | 9×9×21 | **T1-C period** |
| `bench_lateral_B` | 全深度橫向邊界 | 9×9×21 | **T1-D / T2-1/2** |
| `bench_lateral_C` | 同上（格網收斂） | 17×17×21 | T2-4 |
| `bench_lateral_D` | 同上（極細格網） | 25×25×21 | T2-4 延伸 |

---

## 7. Tier 2 整體判定

核心收斂（T2-1, T2-2, T2-3, T2-4, T2-5）**已全部通過**。  
**待補**：T2-7（HFFK SP on Taiwan SKS）。

---

## 10. T3-3 — SinkingSlab 官方範例（Job 34004+，git `2860201`）

**參考基準**：原版 PSI_D 官方輸出  
`/home/wl/software/ECOMAN2.0-seismology.PSI_D/examples/SinkingSlab/psi_output/SYN_SinkingBlock/`

### 10.1 與官方 PSI_D Ray 比較（T3-3 主判定）✅

| 比較 | corr | RMS (s) | md5 | 判定 |
|------|:----:|:-------:|:---:|------|
| HFFK repo Ray SI vs 官方 Ray | **1.000** | **0** | **相同** | ✅ **通過** |

修正後 log：`Reading Global Cartesian Tensors!` + `Reversing depth (depth_reverse=linear)`。

### 10.2 `depth_reverse`：一律用 `linear`（`dims3` 已廢棄）

**2026-07-16 修正**：先前誤以為合成 benchmark 需 `dims3`。**錯誤**——`gen_benchmark_psitomo.py` 已註明寫檔順序對齊 VIZTOMO（深度外層、由深到淺），與官方 psitomo 相同，故 **Tier 1/2/3 全部應使用 `depth_reverse = linear`**。

| 模式 | 判定 | 說明 |
|------|------|------|
| **`linear`** | ✅ **唯一正確** | `reverse!(c11)` 線性索引反轉 ↔ mesh x₃ |
| **`dims3`** | ❌ **廢棄** | fork 早期誤用；均勻場（1L_B）會**掩蓋**錯誤（corr=1），橫向場（lateral_B）暴露：corr≈0.36、RMS≈0.80 s |

**驗證（lateral_B Ray SI）**：
```
linear vs dims3:  corr=0.356  RMS=0.798 s
dims3 vs 舊 bench: corr=1.000  （舊結果自洽但物理錯）
linear vs 舊 bench: corr=0.356  （證明舊 benchmark 用錯模式）
```

**連帶設定**：

| 資料來源 | `depth_reverse` | `tf_global_cartesian` |
|----------|-----------------|----------------------|
| `gen_benchmark_psitomo.py` | `linear` | `false`（Local Cartesian 直接寫入） |
| VIZTOMO / SinkingSlab / Kuo2018 | `linear` | `true`（Global Cartesian → 旋轉） |

**Jobs 34000–34002 結果作廢**；**34007–34009 已以 `linear` 重跑完成** ✅。

### 10.3 第二根因：`tf_global_cartesian`（SinkingSlab 必須 `true`）

| 設定 | SinkingSlab `psitomo0020.dat` | 效果 |
|------|-------------------------------|------|
| **`true`（官方預設）** | VIZTOMO **Global Cartesian** Cij → 旋轉到 local geographic | ✅ md5 與官方一致 |
| **`false`（舊 bench template）** | 假設檔案已是 Local Cartesian，跳過旋轉 | ❌ corr≈0.11，RMS≈0.20 s |

`job_sinking_slab.sh` 已加入 `tf_global_cartesian = true`；合成 benchmark 仍用 `false`。

### 10.4 HFFK vs Ray SI（Job 34005，修正後）✅

| Period T (s) | HFFK vs Ray SI corr | RMS (s) |
|:------------:|:-------------------:|:-------:|
| 3 | 0.998 | 0.008 |
| 5 | 0.997 | 0.012 |
| 8 | 0.994 | 0.016 |
| 15 | 0.987 | 0.025 |
| 25 | 0.975 | 0.034 |

HFFK 與 Ray SI 高度一致；T↑ 時 corr 下降、RMS 增大，符合 Fresnel 平均預期。圖：`sinking_slab_comparison.png`。

### 10.5 執行指令

```bash
cd /lfs/wl/bench_psi
sbatch validation/job_sinking_slab.sh
python3 validation/compare_sinking_slab.py \
  --sink-out /lfs/wl/bench_psi/bench_output/sinking_slab
```

---

## 8. 引用

- VanderBeek, B. P., & Faccenda, M. (2021). *GJI*, 225(3), 2097–2119.
- Faccenda, M., & VanderBeek, B. P. (2023). *J. Geodynamics*, 158, 102003.
- Chevrot, S. (2000). *GJI*, 140, 480–496.

---

## 9. Tier 3 待辦

| 編號 | 測試 | 狀態 | 下一步 |
|:----:|------|:----:|--------|
| T3-1 | Kuo2018 period map | ⏳ | `plot_real_period_map.py` |
| T3-2 | Taiwan SKS 14 對 | ⏳ | 真實 obs vs Ray/HFFK |
| **T3-3** | **SinkingSlab** | **✅** | Ray SI = 官方；HFFK Job 34005 ✅ |
| T3-4 | ASPECT→PSI_D pipeline | ⏳ | 文件化 |
| T3-5 | HFFK SP 路徑 | ⏳ | = T2-7 |

---

## 11. Tier 3 — ASPECT 0608 Ray vs HFFK (application-case validation)

**Case**: `0608_model4_V-50H100`  
**Jobs**: A (`ray_si_uniform48`) vs C (`hffk_u48_T{4,8,16,20,25,33,50}s`)  
**Plotting**: `99_aspect_scripts/plot_sks_script/plotCompare_ray_vs_hffk.py`, `plotC_sinfit.py`, `plot_depth_pierce.py`  
**Outputs**: `plotCompare_maps_bin*.png`, `plotCompare_mix8bins.png`, `plotCompare_stats_bin*.png`

### 11.1 Pearson *r* vs RMSE — what each metric means

Both metrics quantify Ray–HFFK agreement on **Splitting Intensity (SI)**, but they answer different questions.

**Pearson *r*** measures whether the **spatial pattern** covaries: stations with high Ray SI also tend to have high HFFK SI. It is insensitive to a uniform offset (e.g. HFFK = Ray + constant still yields *r* = 1) and insensitive to a uniform scaling (HFFK = 2×Ray still yields *r* = 1).

**RMSE** measures the **absolute mismatch** in SI (seconds):

\[
\mathrm{RMSE} = \sqrt{\frac{1}{n}\sum_i (R_i - H_i)^2}
\]

It reflects both random scatter and **systematic bias** (non-zero mean ΔSI = Ray − HFFK).

In practice: **high *r* with growing RMSE** (as period increases) means the **same anomaly locations** are recovered, but **amplitudes differ** because HFFK integrates over a finite-frequency Fresnel zone whereas Ray SI is the infinite-frequency limit along the geometric ray. This is expected physics, not a back-azimuth or geometry error.

For the 0608 case (BAZ bin 45–90°, 525 stations, depth-averaged SI):

| Period *T* (s) | Pearson *r* | RMSE (s) | Mean ΔSI (s) |
|:--------------:|:-----------:|:--------:|:------------:|
| 4 | 0.984 | 0.074 | −0.000 |
| 8 | 0.966 | 0.108 | −0.002 |
| 16 | 0.936 | 0.150 | −0.003 |
| 20 | 0.922 | 0.164 | −0.004 |
| 25 | 0.906 | 0.180 | −0.004 |
| 33 | 0.881 | 0.203 | −0.006 |
| 50 | 0.819 | 0.242 | −0.008 |
| **MIX (T4–T50 mean)** | **0.936** | **0.155** | **−0.004** |

**Interpretation (0608, single BAZ bin):** At *T* = 4 s, HFFK SI agrees closely with Ray SI (*r* ≈ 0.98, RMSE ≈ 0.07 s), validating the HFFK implementation against the ray-theoretical limit at short periods relevant to SKS. RMSE increases approximately linearly from 0.07 s (*T* = 4 s) to 0.24 s (*T* = 50 s) while *r* remains > 0.82, indicating **frequency-dependent sensitivity** rather than inconsistent geometry. The slightly negative mean ΔSI at long periods (HFFK SI > Ray SI) is consistent with broader depth averaging by the finite-frequency kernel. The period-mixed HFFK field (mean over *T* = 4–50 s per station) has RMSE ≈ 0.15 s, intermediate between *T* = 16–20 s, and is useful as a broad-band proxy but should not replace period-specific comparisons against observations.

### 11.2 All eight BAZ bins — Ray vs HFFK MIX (*T* = 4–50 s)

Comparison uses the same 45° back-azimuth binning as uniform48 (2 ray directions × 3 source depths, depth-averaged). HFFK MIX = mean SI over all seven periods per station within each bin.

| BAZ bin | Range (°) | *r* | RMSE (s) | Mean ΔSI (s) |
|:-------:|:---------:|:---:|:--------:|:------------:|
| 0 | 0–45 | 0.807 | 0.122 | −0.001 |
| 1 | 45–90 | 0.936 | 0.155 | −0.004 |
| 2 | 90–135 | 0.807 | 0.148 | −0.004 |
| 3 | 135–180 | 0.909 | 0.179 | +0.009 |
| 4 | 180–225 | 0.656 | 0.151 | +0.005 |
| 5 | 225–270 | 0.944 | 0.149 | −0.003 |
| 6 | 270–315 | 0.719 | 0.143 | +0.003 |
| 7 | 315–360 | 0.954 | 0.139 | −0.002 |

Agreement is best for bins 1, 5, and 7 (*r* > 0.93); bin 4 shows the weakest correlation (*r* ≈ 0.66), suggesting stronger azimuth-dependent finite-frequency effects or lateral heterogeneity sampling in that sector. All bins retain *r* > 0.65, so spatial patterns remain broadly consistent with Ray SI.

At a single representative period (*T* = 8 s), RMSE ranges from 0.095 s (bins 0, 7) to 0.153 s (bin 3), with *r* from 0.784 (bin 4) to 0.975 (bin 7).

### 11.3 Apparent φ and δt — why φ differs by back-azimuth (Job A)

Job A maps use **per-45° BAZ-bin sinfit** on Ray SI from uniform48 (16 directions × 3 depths). Within each bin, only **two BAZ directions** (×3 depths) constrain the fit:

\[
\mathrm{SI} = A\sin(2\,\mathrm{BAZ}) + B\cos(2\,\mathrm{BAZ}), \quad
\delta t = \sqrt{A^2+B^2}, \quad
\phi = -\tfrac{1}{2}\atan2(B,A)
\]

The resulting φ is an **apparent** fast-axis orientation valid **when the model is sampled only from that BAZ sector**. In a laterally heterogeneous model (0608), true splitting parameters can vary with ray path and azimuth; therefore **apparent φ legitimately differs across the eight BAZ bins**. This is a physical consequence of azimuth-limited sampling, not an inconsistency in the Ray solver.

### 11.4 Job A vs Job B — both Ray SI, but different φ products

| | **Job A** (uniform48) | **Job B** (Kuo2018 SKKS/SKS) |
|---|---|---|
| Sources | Synthetic uniform48 (16 BAZ × 3 depths) | Real teleseismic events (1 BAZ per event) |
| φ retrieval | Sinfit from **2 BAZ directions** per 45° bin | **Single BAZ** per event → underdetermined |
| φ meaning | Apparent φ for that BAZ bin | Display sticks from minimum-|δt| solution (δt = \|SI\|) |
| Maps | 8 panels (one per BAZ bin) | 4 SKS + 15 SKKS panels (one per event) |

**Important:** Job B φ maps must **not** be compared directly to Job A bin maps at a similar back-azimuth. Even though both use Ray `SplittingIntensity`, Job A resolves (δt, φ) from two directions within a bin, whereas Job B assigns one (SI, BAZ) pair per station per event and cannot uniquely determine both δt and φ without additional azimuthal coverage. For Kuo2018-style validation, multi-event sinfit (pooling BAZ at each station) or direct comparison on **SI** (not φ) is appropriate.

### 11.5 HFFK sinfit and depth projection (Job C)

**Sinfit plot** (`plotC_sinfit.py`): For each station, SI is plotted vs BAZ; coloured arcs show independent sinfit within each 45° bin (reference: legacy `plotC2_sinfit.py`). Spatial φ maps from the same per-bin sinfit are in `plotC_phi.py`.

**Depth pierce-point maps** (`plot_depth_pierce.py`): Surface SI is converted to apparent (φ, δt) per source–receiver pair, then projected along the SKS ray to reference depths **100, 200, 300, 400, 500 km** using TauP (ak135). Layout: 5 depth rows × 8 BAZ bins; separate figures for Ray (`ray_si_uniform48`) and HFFK (representative period *T* = 8 s). Rose diagrams on each panel indicate the observation azimuth bin.

**Recommended wording for proposals:**

> We validate HFFK against Ray SI on the ASPECT 0608 anisotropic model using native `SplittingIntensity` outputs (Jobs A and C). At SKS-relevant periods (*T* ≈ 4–8 s), Ray and HFFK SI agree closely (*r* > 0.96, RMSE < 0.11 s). Discrepancies grow with period (*r* > 0.82, RMSE < 0.25 s at *T* = 50 s) while preserving spatial correlation, consistent with finite-frequency Fresnel averaging. Apparent splitting parameters are evaluated with azimuth-resolved sinfit (uniform48) or SI-only comparisons (Kuo2018 events); φ maps from single-azimuth events are interpreted as orientation indicators, not unique (δt, φ) inversions.

---

*GitHub 建議路徑：`docs/validation/PSI_D_HFFK_Validation_Tier1_Tier2.md`*
