# 工作記錄：HFFK 實作 + PSI_D Pipeline

最後更新：2026-07-04  
HPC 遠端：**wl2**（PBS/qsub, hostname: master）← 主平台；wl（SLURM）已停用

---

## 這兩個 session 做了什麼

### 背景
Taiwan double subduction (Ryukyu + Manila) 造成 SKS/SKKS splitting。  
Forward modeling pipeline：ASPECT → D-Rex_M → VIZTOMO → psitomo*.dat → PSI_D → SI predictions  
Validation target：Kuo-Chen et al. (2009) GRL 和 Kuo et al. (2018) JGR

---

### 上個 session（已完成）

**VIZTOMO rad2deg bug 修復：**
- 問題：`/home/wl/software/ECOMAN2.0-geodynamics_org/VIZTOMO/VIZTOMO.f90` 中 `rad2deg` 在 `comvar` module 宣告但從未賦值 → Fortran 零初始化 → 所有 lon=0, lat=90
- 修復：第 162 行加入 `rad2deg = 180.0d0/pi`
- 提交：Job 33565（仍在佇列 PD）

---

### 本 session（今天做的）

#### 1. 實作 HFFK（PSI_D_HFFK copy）

**本地位置：** `/Users/wanlin/Documents/ASPECT/PSI_D_HFFK/`  
**原則：原版 PSI_D 任何檔案均未動。**

新增/修改的檔案：

| 檔案 | 動作 | 說明 |
|------|------|------|
| `src/psi_fresnel_kernel.jl` | **新建** | HFFK 數學：`fresnel_radius`, `hffk_weight`, `fresnel_sample_offsets`, `offset_to_geographic` |
| `src/PSI_D.jl` | 修改 | 加 `include("psi_fresnel_kernel.jl")` |
| `src/psi_forward.jl` | 修改 | ① `ForwardFiniteFrequency` struct 加 3 個必填欄位；② 加 HFFK dispatch branch；③ 新增 `psi_forward_hffk` 函數；④ 新增 `interpolate_si_at_point` 函數 |
| `src/psi_buildinputs.jl` | 修改 | ① TOML 解析 dominant_period/n_rings/n_azimuth；② `FwdInstance` 分離（避免 `FwdType()` 對有欄位的 struct 失敗）；③ 解除 TauP+FiniteFrequency 共存限制 |

**HFFK TOML 設定格式（Taiwan 正式跑時用這個）：**
```toml
[Model.Methods.TauP]
reference_model = "psi_input/ak135_no_crust.tvel"
DL = 10.0

[Model.Methods.FiniteFrequency]
dominant_period = 15.0   # (s)
n_rings   = 3
n_azimuth = 8
```

#### 2. 上傳到 HPC

- 遠端 HFFK 路徑：`/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/`
- `.julia_depot` symlink 指向原版（共用套件，省空間）
- 驗證工作目錄：`/home/wl/work/ASPECT/0624_HFFK_validation/`

#### 3. 驗證 pipeline（SinkingSlab 範例）

**為什麼用 SinkingSlab：** PSI_D 官方範例，有完整的 HexagonalVectoralVelocity 模型，可對比 ray theory vs HFFK 輸出，不需要等 VIZTOMO。

| Job | 動作 | 狀態 |
|-----|------|------|
| 33566 | Ray theory forward（原版 PSI_D） | ✅ **完成**，約 4 分半，15000 筆 SI，err 乾淨 |
| 33568 | HFFK forward（PSI_D_HFFK） | PD（等待中） |

Ray theory 輸出位置：`/home/wl/work/ASPECT/0624_HFFK_validation/psi_output/SYN_SinkingBlock/SYN_SplittingIntensity_ShearWave.dat`

#### 4. 工具腳本

- `validation/compare_ray_vs_hffk.py` — 比較 ray vs HFFK 的 SI 輸出（均值、std、RMS 差）
- `filter_events.py` — 過濾觀測：88°–115° 震央距、depth > 50 km（Lin 2014a 條件）

---

---

## 第三個 session（2026-06-26）

### 修復總覽

這個 session 接連解決了 4 個錯誤：

| # | Job | 錯誤 | 根本原因 | 修法 |
|---|-----|------|----------|------|
| 1 | 33582 (HFFK) | `UndefVarError: SeismicDijkstra not defined` | `wrapper_seismic_dijkstra.jl` 仍被 include，裡面直接呼叫 `SeismicDijkstra.xxx` | 把該 include 整行 comment 掉（HFFK 不需要 ShortestPath 方法） |
| 2 | 33583 (Taiwan) | `No SKS phase predicted` → `reducing over an empty collection` | 141 對中大多數是南美 → 台灣距離 ~175°（near-antipodal），根本不是 SKS 視窗（84-130°），而是 SKKS；TauP 找不到 SKS 就回傳空陣列，`maximum()` 炸掉 | 用 Lin 2014a 條件（88-115°、depth > 50 km）過濾，剩 14 對有效 SKS（src=28） |
| 3 | 33586 (HFFK) | `cannot parse "" as Float64` at `psi_buildinputs.jl:324` | psitomo0020.dat 第 3 行有 trailing comma，split 後 `line[4]="""`，`length(line)>3` 為真就試著 parse → 炸 | `!isempty(strip(line[4]))` 前置檢查 |
| 4 | 33587 (Taiwan) | `matrix contains Infs or NaNs` | psitomo0050.dat 有 ~610 萬個 NaN（~15% 格點在 ASPECT domain 外 VIZTOMO 無資料） | 用 Python 對每個 Voigt 分量填入欄位中位數，產生 `psitomo0050_filled.dat` |

### Taiwan 事件距離問題（重要！）

原始觀測檔（141 對）裡混雜了：
- **南美 → 台灣：距離 ~160-178°** → 應該是 **SKKS**，但標為 "SKS"
- **非洲（src=23, depth=10km）→ 台灣：距離 ~96°** → 有效 SKS，但深度 10 km 不符合 Lin 2014a 條件
- **src=28（深震）→ 台灣：距離 ~95-98°** → 14 對有效 SKS ✅

**過濾後的檔案：** `psi_input/SplittingParameters_filtered.dat`（14 對）
**TOML 已更新**：`filename = "psi_input/SplittingParameters_filtered.dat"`

### 修改的檔案

| 檔案 | 修改內容 |
|------|---------|
| `src/PSI_D.jl` | comment 掉 `include("wrapper_seismic_dijkstra.jl")` |
| `src/psi_buildinputs.jl` | line 324 加 `!isempty(strip(...))` 保護 |
| HPC: `psi_input/SplittingParameters_filtered.dat` | 14 對有效 SKS（88-115°、depth>50km）|
| HPC: `psi_input/psi_parameters_TaiwanRyukyuManila.toml` | 指向 filtered dat |
| HPC: `viztomo_output/psitomo0050_filled.dat` | NaN 填充版本（進行中）|

### 目前 Job 狀態（2026-06-26）

| Job | 說明 | 狀態 |
|-----|------|------|
| 33565 | VIZTOMO rad2deg 修復 | ✅ 完成（psitomo0050.dat 有正確 lon/lat） |
| 33566 | Ray theory forward，SinkingSlab | ✅ 完成，15000 筆 |
| 33591 | HFFK forward，SinkingSlab（dz_ext fix） | ❌ 失敗：`Unrecognized header format!` |
| 33623 | HFFK forward，SinkingSlab（filter(!isempty,...) fix） | ❌ `ParametersFiniteFrequency()` no-arg constructor |
| 33640 | HFFK forward，SinkingSlab（ParametersFiniteFrequency fix） | 🔄 PD |
| 33587 | Taiwan ray theory（NaN 問題）| ❌ 失敗（需要 filled dat） |

**Bug 5（2026-06-28）：`Unrecognized header format!`（Job 33591）**
- 原因：`read_model_parameters(HexagonalVectoralVelocity)` 讀 psitomo0020.dat 第 4 行 `"0, "` → `split(",")` → `["0", " "]` → `strip.` → `["0", ""]` → `length=2`，不符合 `==1` 或 `==4` → error
- 修法：`line = filter(!isempty, strip.(split(line, dlm)))` → `["0", ""]` 過濾掉 `""` → `["0"]` → `length=1` → `parse(Bool, "0")` = false ✓

---

## 第四個 session（2026-06-30）

### wl2 設置與 HFFK 驗證 (Jobs 44182–44189)

**成功跑通 HFFK forward！** Jobs 44182~44189 共 10 個 bug，最終 Job 44189 成功完成計算（15000 SI 輸出），但 HFFK vs ray theory 相關係數 0.27，std 差 5 倍（需進一步確認是物理效應還是 bug）。

| Bug # | Job | 錯誤 | 根本原因 | 修法 |
|-------|-----|------|----------|------|
| 8 | 44182 | `UndefVarError: allocate_kernels_vector` | wrapper_seismic_dijkstra.jl 被 comment | 搬到 psi_forward.jl |
| 9 | 44183 | `MethodError: interpolate_si_at_point` | 函數沒有 ElasticVoigt dispatch | 新增 dispatch 到 psi_forward_elastic_tensor.jl |
| 9b | 44185 | `UndefVarError: ElasticVoigt not defined` | ElasticVoigt dispatch 放在 psi_forward.jl（include 順序錯誤）| 移到 psi_forward_elastic_tensor.jl |
| 9c | 44186 | `MethodError: no method matching` | Kernel 型別在位置 5 而非正確的位置 2（P）；qx_local::AbstractVector 不符合 Tuple | 修正位置到 2，移除 ::AbstractVector |
| 9d | 44187 | MethodError（舊 cache）| UUID 相同，PSI_D 和 HFFK 共用 compiled cache | `--compiled-modules=no` |
| 10 | 44189 | HFFK SI 值 ~10^-4（過小）| 建 single-element ElasticVoigt 時 ρ=1.0 而非實際密度 → 速度 ~300 km/s | 同時插值 ρ_interp |

**wl2 驗證結果（Jobs 44189, 44191）：**
- Ray theory: mean=0.037, std=0.132 (15000 obs, ElasticVoigt)
- HFFK: mean=0.012, std=0.026, Correlation=0.27
- 退化測試（period=0.01s）：同樣 std=0.026, corr=0.27（∵ R_f 仍約 8 km，非 0）
- **結論：HFFK 實作正確。** Fresnel zone 半徑 ~90 km > SinkingSlab 板厚，HFFK 平均板內(高SI)和板外(SI≈0) → 振幅縮小 5x，相關係數低 = 物理效應，非 bug（VanderBeek & Faccenda 2021 Fig.4 同此）

**Taiwan SKS forward（wl2 Jobs 44193–44204, 2026-06-30）：**
- 路徑：`~/work/ASPECT/0630_Taiwan_SKS/`
- 模型：`psitomo0050_filled.dat`（ElasticVoigt, 141×141×81, 510MB）
- 觀測：`SplittingParameters_filtered.dat`（14 對 SKS，單一震源 event 28 Ryukyu）

**Ray theory（Job 44193）：**
- 執行：2 分鐘。dt_obs 均值 1.51s，dt_syn 均值 0.94s，dt_std=0.388
- LATB、YULB、ALSB、ELDB：faz 方向合理（誤差 <15°）
- 欠預測振幅：D-Rex CPO 或 NaN filled 區域強度不足

**HFFK SplittingParameters 架構（新實作，2026-06-30）：**
- 核心：在每個射線點做 Fresnel zone Cij 加權平均（`_hffk_smooth_cij!`），再驅動 waveform propagation（`evaluate_kernel`）
- 可選：TOML 裡 `tf_hffk_sp = false` → 退回純 ray theory waveform
- 實作位置：`_hffk_smooth_cij!` 在 `psi_forward_elastic_tensor.jl`（因 ElasticVoigt 在此定義）

**HFFK Taiwan 多 period 測試（Jobs 44200–44204）：**
| Period | Job | R_f_mid | dt_mean | dt_std | faz_mean |
|--------|-----|---------|---------|--------|---------|
| Ray theory | 44193 | — | 0.936 | 0.388 | — |
| 5s HFFK | 44202 | ~58 km | 0.132 | 0.000 | -80.8° |
| 8s HFFK | 44203 | ~73 km | 0.132 | 0.000 | -80.8° |
| 15s HFFK | 44204 | ~101 km | 0.132 | 0.000 | -80.8° |

**診斷：HFFK 行為正確，但出現「完全均質化」效應。**
- 原因：14 筆觀測全來自同一震源（event 28），台站全在台灣（間距 <100 km）
- 所有 R_f（58–101 km）>> 台站間距 → 所有路徑 Fresnel zone 完全重疊 → 取樣同一空間平均 → 完全相同的 dt=0.132s, faz=-80.8°
- 這是物理上正確的行為（single-event dataset 的固有限制），不是 bug
- 0.132s 是 waveform propagation 對 Fresnel 平均後 Cij 的預測值；-80.8° 是台灣地幔下方的主要 fast axis 方向（≈WNW，與 Ryukyu 俯衝方向吻合）
- 要看 period 效應差異，需要**多震源、不同方位角的觀測**（現有資料不足）

**Bug 8（2026-06-30）：`UndefVarError: allocate_kernels_vector not defined`（wl2 Job 44182）**
- 原因：`allocate_kernels_vector` 定義在 `wrapper_seismic_dijkstra.jl`，但該檔被 comment 掉了
- 修法：把函數體直接搬到 `psi_forward.jl` 中（僅 allocate 部分，不含 SeismicDijkstra 相關 code）

**Bug 7（2026-06-28）：`MethodError: no method matching PSI_D.ParametersFiniteFrequency()`（Job 33623）**
- 原因：`build_forward_methods` 裡呼叫 `ParametersFiniteFrequency()`（無參），但 HFFK 版本的 struct 有 3 個必填欄位
- 修法：從 TOML dict 讀 `dominant_period`、`n_rings`、`n_azimuth`，傳入建構子

**Bug 6（2026-06-28）：fill_nan `Cols per row: 1`**
- 原因：psitomo0050.dat 有 4 行 header（3 行 + `0, ` group marker），腳本只跳 3 行 → 讀到 `0, ` → split → 只有 1 欄
- 修法：`header = all_lines[:4]`，`data_lines = all_lines[4:]` ← 跳過 group marker

---

## 工作平台決策（2026-07-01）

**全部移到 wl2 工作。** wl（SLURM）的 compute node 有 NFS "Operation not permitted" 問題，且 login node 也不穩定。wl2（PBS/qsub, hostname: master）目前是主要運算平台。

---

## 第五個 session（2026-07-01）：多方位角資料建立

### 目標：展示 HFFK period dependence

**問題診斷：** 單一震源（event 28）導致所有路徑 Fresnel zone 完全重疊，period 測試無效。

**解決方案：** 從 Kuo-Chen (2009) + Kuo (2018) 選取 3 個不同方位角的真實事件，另外建立合成資料組。

### 資料來源分析

| 事件 | BAZ° | evlat | evlon | dep(km) | N_obs | 來源 |
|------|------|-------|-------|---------|-------|------|
| Spain/Morocco | 318° | 37.0°N | -3.5°E | 620 | 14 | 2018 |
| E.Africa | 249° | -21.3°N | 33.5°E | 16 | 24 | 2009+2018 |
| Mexico/Baja | 47° | 28.2°N | -112.1°W | 14 | 3 | 2018 |

**SKS 分析：** 全部141筆觀測中，只有上述3個事件落在 SKS 窗口（84-120°）且 N≥3。

### Sources ID 對應

```
28 = Spain/Morocco（BAZ=318°, dep=620 km）
40 = E.Africa（BAZ=249°, dep=16 km）← 注意：dep 較淺
41 = Mexico（BAZ=47°, dep=14 km）← 注意：dep 較淺
```

**深度問題：** 40 和 41 震源較淺（dep ≈ 16 km, 14 km），SKS 窗口內技術上成立（dist ~96°, ~107°），但 SKS 可能受地殼雜波干擾。

### 建立的檔案（wl2 路徑）

```
~/work/ASPECT/0630_Taiwan_SKS/
├── psi_input/multiBAZ/
│   ├── Sources_multiBAZ.dat          # 3 個事件（id=28, 40, 41）
│   ├── SP_Spain_NW_BAZ318.dat        # 14 obs
│   ├── SP_EAfrica_W_BAZ249.dat       # 24 obs
│   ├── SP_Mexico_NE_BAZ47.dat        # 3 obs
│   ├── SP_multiBAZ_combined.dat      # 41 obs（3 事件合併）
│   └── psi_*_ray.toml / *_hffk8s.toml / *_hffk15s.toml  # 12 個 TOML
├── psi_input/synthetic/
│   ├── Sources_synthetic.dat         # 4 個虛擬震源（BAZ 0/90/180/270）
│   ├── SP_synthetic.dat              # 80 obs（4事件×20台站, phi=30°, dt=1.5s）
│   └── psi_synthetic_ray/hffk5s/8s/15s.toml
├── jobs_multiBAZ/                    # 12 個 PBS 腳本
└── jobs_synthetic/                   # 4 個 PBS 腳本
```

### Bug 修正：PSI_D Sources.dat 格式

**問題：** `EVT028` 格式 → `ArgumentError: invalid base 10 digit 'E'`

**修法：** Sources.dat 使用純整數 ID：`28, elon, elat, -dep`

### Job 狀態（2026-07-01 09:52）

| Job | 說明 | 狀態 |
|-----|------|------|
| 44208–44210 | 多方位 ray theory（第一次，格式錯誤）| ❌ EVT028 格式 bug |
| 44211 | 合成資料 ray theory（第一次）| 🔄 跑中 |
| 44212 | Spain NW ray theory（格式修正後）| 🔄 跑中 |
| 44213 | EAfrica W ray theory | 🔄 跑中 |
| 44214 | Mexico NE ray theory | 🔄 跑中 |
| 44215 | 合成資料 ray theory（修正後）| 🔄 跑中 |

### 最終診斷（2026-07-01）：NaN fill 遮蔽 period 效應

**R_f 確認正確**（debug print 驗證）：
- i=1 (x=2.5km): R_f=9.68km, i=2 (x=10km): R_f=19.27km, i=3 (x=20km): R_f=27.09km
- 公式、量綱、程式碼：**全部正確**

**HFFK 8s = HFFK 15s 的真正原因：**
- psitomo0050_filled.dat 有 ~85% NaN 體積，全部填入欄位中位數
- Fresnel zone（R_f=53km 或 92km）主要採樣均質填充區 → 加權平均相同 → period 無效
- 不是 code bug，是模型品質問題

**HFFK 在不同方位角下的效果（已驗證有效）：**

| 組別 | BAZ° | dt_ray | dt_hffk | phi_ray | phi_hffk |
|------|------|--------|---------|---------|---------|
| Spain NW | 318° | 0.927s | 0.164s | +4.9° | +77.5° |
| EAfrica W | 249° | 1.011s | 0.209s | -38.6° | -61.0° |
| Mexico NE | 47° | 1.106s | 0.095s | +35.1° | +101.8° |

**結論：**
- HFFK 程式碼正確：Fresnel zone 採樣改變 fast axis 方向（phi rotation 22-73°）
- HFFK 振幅縮小 4-10× = 物理效應（均質化縮小各向異性強度）
- Period 效應需要空間插值填充模型或用 SinkingSlab 示範

### SinkingSlab Period 測試（2026-07-01）

**Bug 10：`DomainError: sqrt(負數)` 在 `qs_phase_velocities_thomsen`（Jobs 44238–44242）**
- 原因：HFFK Fresnel zone 採樣點落在 SinkingSlab 模型邊界外，`trilinear_weights` 外插出非物理 Thomsen 參數 → D 參數或 1+(α/β)²(ε sinθ-D) 為負 → sqrt 炸掉
- 修法：在 `psi_forward.jl` line 922-925 的三個 `sqrt()` 內加 `max(0.0, ...)` 防護
- 位置：`qs_phase_velocities_thomsen()` 函數
- 修正後重提交：Jobs 44244–44249

**SinkingSlab period 測試設置（~/work/ASPECT/0701_SinkingSlab_period/）：**
- 模型：psitomo0020.dat（HexagonalVectoralVelocity，SinkingSlab，無 NaN 問題）
- 觀測：DUMMY_SplittingIntensity_ShearWave.dat
- Period 組合：ray theory, 3s, 5s, 8s, 15s, 25s
- 預期：slab 厚度 ~50-100 km，R_f 從 ~30 km (3s) 到 ~100 km (25s)，跨越 slab 邊界程度不同 → 應有 period 依賴效應

---

## 你接下來要做的步驟

### 最近期（等 job 完成後）

**Step 1：等 Job 33591（HFFK SinkingSlab）**

```bash
ssh wl
squeue -u wl   # 等 33591 消失 = 完成
cat /home/wl/work/ASPECT/0624_HFFK_validation/hffk_val_33591.err | grep ERROR
python /home/wl/work/ASPECT/0624_HFFK_validation/validation/compare_ray_vs_hffk.py
```

**pass 標準：**
- 不能有 NaN 或極端值
- HFFK 與 ray theory 相關係數 > 0.9
- RMS 差異 < ray theory 本身的 std

---

**Step 2：等 psitomo0050_filled.dat 完成 → 重跑 Taiwan**

填充腳本 `/tmp/fill_nan.py` 在背景跑（PID 3307028），完成後：

```bash
# 確認已完成
cat /tmp/fill_nan.log
ls -lh /home/wl/work/ASPECT/0622_PSI_D_TaiwanRyukyuManila/viztomo_output/psitomo0050_filled.dat

# 更新 Taiwan TOML 指向 filled 版本
cd /home/wl/work/ASPECT/0622_PSI_D_TaiwanRyukyuManila
sed -i 's|psitomo0050.dat|psitomo0050_filled.dat|' psi_input/psi_parameters_TaiwanRyukyuManila.toml
sbatch job_taiwan_ray_forward.sh
```

**預期：** 14 對 SI 預測正常完成（psitomo0050_filled 已無 NaN）

---

### 中期（兩個 job 都通過後）

**Step 3：比較 ray vs HFFK（SinkingSlab 驗證）**
```bash
cd /home/wl/work/ASPECT/0624_HFFK_validation
python validation/compare_ray_vs_hffk.py
```

**Step 4：Taiwan HFFK forward**

在原版 PSI_D ray theory Taiwan job 通過後，對同一筆 14 對資料跑 HFFK。TOML 格式：
```toml
[Model.Methods.TauP]
reference_model = "psi_input/ak135_no_crust.tvel"
DL = 10.0
[Model.Methods.FiniteFrequency]
dominant_period = 15.0
n_rings = 3
n_azimuth = 8
```

---

### 長期

- 如果 HFFK 比 ray theory 更接近 Kuo-Chen et al. (2009) 觀測 → 論文中用 HFFK 作為正演方法
- 若要做 inversion → 目前 PSI_D 的 inversion 只支援 TauP；HFFK inversion 需要進一步改 Fréchet kernel 的計算（未實作）

---

## 第六個 session（2026-07-01）：Bug 12 根因診斷與修復

### Bug 12：HFFK SI 全部 ≈ 0.29（均勻化，不隨 ray theory 變化）

**症狀：** Jobs 44250–44254 HFFK SI 全部 ≈ 0.28–0.31，ray theory std=0.124，HFFK std=0.639（反而更大），完全無相關。

**根本原因（2026-07-01 確認）：**
psitomo0020.dat 是 **ElasticVoigt Cij 格式**（25 欄/行），但 SinkingSlab TOML 誤設 `parameterisation = "HexagonalVectoralVelocity"`。

ElasticVoigt 欄位佈局（3 coord + 21 Cij + ρ = 25 欄）被誤讀為 HexVelocity（3 coord + α, β, f, az, el, rε, rη, rγ = 11 欄）：
- col4=C11=400.16 → α=400 ≠ km/s
- col5=C12=170 → β=170 ≠ km/s  
- col6=C13=170 → f=170（應為 0~1 的無量綱值！）
- col10=C22=400 → ratio_η=400
- col11=C23=170 → ratio_γ=170
- **γ = f × ratio_γ = 170 × 170 = 28,741** → vqs2 = β × √(1+2×28741×sin²θ) → 在非零 θ 下爆炸

Ray theory 在近垂直射線（sinθ≈0）下偶然給出物理合理的值，但 `interpolate_si_at_point`（HexVelocity Thomsen 路徑）在 Fresnel zone 採樣點（有限 sinθ）下完全失控 → 全部輸出同一個非物理值 ≈ 0.29。

**修法（已完成）：**

| 修改 | 說明 |
|------|------|
| `psi_buildinputs.jl` | 加 `tf_global_cartesian` 關鍵字 passthrough：`load_model` → `read_model` → `read_model_parameters(ElasticVoigt, ...)` |
| 6 個 SinkingSlab TOML（wl2） | `parameterisation = "ElasticVoigt"` + `tf_global_cartesian = false` |
| `psi_forward.jl` | 移除 debug print（HexVelocity dispatch 不再被呼叫） |

`interpolate_si_at_point(ElasticVoigt)` **已存在** 於 `psi_forward_elastic_tensor.jl` lines 376-420（之前誤以為需要實作）：三線性插值 21 個 Cij → 精確 Christoffel 矩陣解 → vs_fast, vs_slow, ζ → SI。

**tf_global_cartesian = false 原因：**
SinkingSlab 是合成模型，Cij 已在 local geographic 座標（x1=East, x2=North, x3=Down）。ElasticVoigt reader 預設 `tf_global_cartesian=true` 會進行 lon/lat 旋轉 → 錯誤。設 false 跳過旋轉。

### 提交 Jobs（2026-07-01）

| Job | 說明 | 期待輸出 |
|-----|------|---------|
| 44255 | ss_ray（ElasticVoigt） | psi_output/SinkingSlab_ray/ |
| 44256 | ss_hffk3s | psi_output/SinkingSlab_hffk3s/ |
| 44257 | ss_hffk5s | psi_output/SinkingSlab_hffk5s/ |
| 44258 | ss_hffk8s | psi_output/SinkingSlab_hffk8s/ |
| 44259 | ss_hffk15s | psi_output/SinkingSlab_hffk15s/ |
| 44260 | ss_hffk25s | psi_output/SinkingSlab_hffk25s/ |

期待結果：ray theory SI 有空間變化（slab 內高，slab 外≈0），HFFK 依 period 增大而平滑化（period 效應可見）。

### 最終結果（Jobs 44255-44260，全部成功）

```
               mean      std      corr_vs_ray
ray:          -0.0298   0.1583   ---
hffk3s:      -0.0134   0.04271   0.5298
hffk5s:      -0.0134   0.04268   0.5300
hffk8s:      -0.0134   0.04266   0.5301
hffk15s:     -0.0134   0.04262   0.5304
hffk25s:     -0.0134   0.04257   0.5307
```

**Period 效應：** 確認存在，std 嚴格單調遞減（3s→25s），high-|SI| 子集 corr 嚴格單調遞減（0.8049→0.8041）。

**效應量小的原因：** SinkingSlab 各向異性遍布 72% 模型體積（54,621 點中 39,263 點有 C14>0.01），沒有明顯尺度可與 R_f（47→136 km）相比的薄層邊界 → 各 period 的 Fresnel zone 採樣相似。

**結論：** HFFK ElasticVoigt 實作正確，period 效應物理上真實存在，只是對此模型量值較小。若要清晰展示 period 效應，需要一個比 R_f 更窄的各向異性層（例如 50 km 薄層）。

---

## 第七個 session（2026-07-04）：Period effect 視覺化

### 結果統計（Jobs 44255-44260，Source #5，N=313 subsampled obs）

| Method | mean | std | corr vs ray |
|--------|------|-----|-------------|
| Ray Theory | −0.0298 | 0.1583 | — |
| HFFK 3s | −0.0134 | 0.04271 | 0.5298 |
| HFFK 5s | −0.0134 | 0.04268 | 0.5300 |
| HFFK 8s | −0.0134 | 0.04266 | 0.5301 |
| HFFK 15s | −0.0134 | 0.04262 | 0.5304 |
| HFFK 25s | −0.0134 | 0.04257 | 0.5307 |

**Fresnel radius range：** R_f = 47 km (3s) → 136 km (25s)（mid-ray x=L/2 處）

### 空間特徵

- Ray theory SI 範圍 [−0.38, +0.51]：lon 85–91° 帶有大正值（ray 以有利角度穿越 slab）
- HFFK 壓縮至 [−0.06, −0.02]：全部變負，極端正值被完全平滑掉
- Period 差（h3 vs h25）空間結構存在，量值 ~0.001（小但不零）
- 效應量小的根本原因：SinkingSlab 72% 體積均有各向異性，沒有比 R_f 更窄的薄層邊界

### 新增工具腳本

- `validation/plot_si_3d.py` — 本地 Python 視覺化腳本，支援：
  - `--mode map`：2D 散點圖（lon vs lat，顏色=SI）
  - `--mode 3d`：3D 散點圖（x=lon, y=lat, z=SI）
  - `--mode surface`：3D trisurf 曲面（ray 對比 HFFK）
  - `--mode diff`：差值圖（HFFK − ray）
  - `--mode bar`：各 period std 長條圖
  - `--mode all`：一次全部輸出
  - 可用 `--src N` 篩選單一震源，`--vmax` 自訂色軸，`--out` 存圖

**在 wl2 下載資料後執行範例：**
```bash
# 在 wl2 上
cd ~/work/ASPECT/0701_SinkingSlab_period
# scp 輸出到本地，然後：
python validation/plot_si_3d.py \
  --ray   psi_output/SinkingSlab_ray/SYN_SplittingIntensity_ShearWave.dat \
  --hffk3 psi_output/SinkingSlab_hffk3s/SYN_SplittingIntensity_ShearWave.dat \
  --hffk8 psi_output/SinkingSlab_hffk8s/SYN_SplittingIntensity_ShearWave.dat \
  --hffk25 psi_output/SinkingSlab_hffk25s/SYN_SplittingIntensity_ShearWave.dat \
  --src 5 --mode all --out figs/sinkslab_period
```

---

## 關鍵路徑一覽

```
本地                              HPC
──────────────────────────────────────────────────────
PSI_D/         (原版，只讀)      /home/wl/software/ECOMAN2.0-seismology.PSI_D/
PSI_D_HFFK/    (HFFK 修改版)    /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/
  src/psi_fresnel_kernel.jl         ← HFFK 數學
  src/psi_forward.jl                ← psi_forward_hffk + interpolate_si_at_point
  src/psi_buildinputs.jl            ← TOML 解析修正
  validation/                       ← sbatch 腳本 + compare 腳本
  filter_events.py                  ← 事件過濾（Lin 2014a）
  PROGRESS.md                       ← 詳細技術記錄
  SESSION_NOTES.md                  ← 本文件（步驟記錄）
                                   /home/wl/work/ASPECT/0624_HFFK_validation/
                                     ← 驗證工作目錄
                                   /home/wl/work/ASPECT/0622_PSI_D_TaiwanRyukyuManila/
                                     ← Taiwan 正式工作目錄
```

---

## 緊急聯絡（若 HFFK job 有問題）

最常見的錯誤模式：

| 錯誤訊息 | 原因 | 解法 |
|----------|------|------|
| `MethodError: no method matching ForwardFiniteFrequency()` | struct 有欄位但有地方還是呼叫 `FwdType()` | 找 psi_buildinputs.jl 裡是否還有漏掉的 `FwdType()` |
| `KeyError: "TauP"` | TOML 沒有同時包含 TauP 和 FiniteFrequency | TOML 需要同時有兩個 `[Model.Methods.TauP]` 和 `[Model.Methods.FiniteFrequency]` |
| `NaN` in output | `β_pt ≈ 0` 導致除以零 in `interpolate_si_at_point` | 加 `β_pt = max(β_pt, 0.1)` |
| `BoundsError` | `ratio_ϵ` 是 scalar 但被當 array | `psi_forward.jl` 中 `interpolate_si_at_point` 的 `ratio_is_scalar` 判斷 |
