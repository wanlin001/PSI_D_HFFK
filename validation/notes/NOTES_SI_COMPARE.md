# SKS/SI 比較筆記（plot_sks_script）

持續記錄：定義差異、已修 bug、待重跑事項。日期為本地工作日。

---

## 2026-07-20

### SI 定義：為何 PSI 有 ½、觀測圖沒有

| 來源 | 公式 | 備註 |
|------|------|------|
| **觀測圖** (`plotOBS_Kuo2018_SI`, `psi_common.obs_si_from_dt_phi`) | `SI_obs = δt · sin(2·(BAZ−φ))` | 用測得的 δt、φ 轉成「參數型」SI，**無 ½** |
| **PSI_D kernel** (`psi_forward.jl` ~812) | `si = ½·(1/v_slow − 1/v_fast)·sin(2ζ)` | 沿射線積分；單層時 ≈ `½·δt·sin(2β)` |
| **Benchmark 註解** (`validation/gen_benchmark_psitomo.py`) | `SI_HALF = 0.5` | 明寫：要與 PSI_D 對齊必須乘 0.5；Chevrot 2000 慣例 |

**為什麼要 ½？**  
δt = ∫(1/v_slow − 1/v_fast) dr。Chevrot / PSI 的 splitting intensity 是從橫向能量投影來的振幅，小分裂近似下振幅是 **½·δt·sin(2ζ)**，不是 δt·sin(2ζ)。這是 kernel 設計，不是漏乘。

**`MODEL_SI_TO_OBS_SCALE = 2.0` 的意義**  
只用於 **obs vs model 比較圖**：把 PSI 的 SYN SI ×2，換成與觀測圖同一套「δt·sin2」刻度。  
**不改** PSI 輸出檔、也不改 Job A/C 內部 Ray vs HFFK 比較（兩邊都是 PSI 慣例，不必 ×2）。

### paz = π/4 bug（相位）

- **問題**：DUMMY SI 若 `paz=π/4`，forward 用 `ζ = φ − BAZ − paz` → 預測近似 −cos 而非 sin，**圖形旋轉**，振幅包絡可仍在。
- **已修（DUMMY 輸入）**：
  - 較早：`DUMMY_SI_uniform48` / `uniform96` 已是 `paz=0`（Job A/C 正確）。
  - 2026-07-20：patch 所有 `DUMMY_SI_Kuo2018*.dat`（0608 project + PSI_D_HFFK `psi_input/`）→ `paz=0`。
  - 腳本：`PSI_D_HFFK/scripts/gen_dummy_obs.py`（說明 + `--patch`）。
- **尚未生效（SYN 輸出）**：現有  
  `ray_si_Kuo2018_SKS/SKKS`、`hffk_Kuo2018_*`  
  檔內仍寫 `paz=0.785…`，因 **Job B/D 尚未用新 DUMMY 重跑**。  
  → 修 DUMMY ≠ 修舊 SYN；要重跑 Job B/D 才改相位。

核對指令：
```bash
# DUMMY 應為 0
awk -F',' 'NR==1{print $NF}' psi_input/DUMMY_SI_Kuo2018_SKS.dat
# SYN 重跑前仍為 π/4；重跑後應為 0
awk -F',' 'NR==2{print $NF}' psi_output/ray_si_Kuo2018_SKS/SYN_SplittingIntensity_ShearWave.dat
```

### 為何 `plotCompare_ray_vs_hffk` 看起來 SI「比較強」？

那支腳本比的是 **Job A vs Job C（uniform48）**，不是 vs Kuo 觀測：

1. **paz=0**（已正確）→ 相位對、場看起來乾淨。
2. **16 方位 × 3 深度**，總有接近最大分裂角的方向 → TW 區 |SI| p90 ~0.7 s、max ~1.6 s（raw PSI）。
3. 色階 `SI_VLIM=1.0`，點用 |SI| 放大，視覺上很「滿」。
4. **沒有觀測 δt·sin2（無 ½）當對照**，所以不會覺得模型「矮一截」。

Obs vs model（Job B/D + Kuo）則是：觀測 |SI|~1.3 s（無 ½）對 raw 模型 ~0.2 s，再 ×2 後仍可能因 **舊 SYN 的 paz=π/4** 與幾何而對不齊。

### 輸出檔名（obs vs model）

- `obs_vs_{model}_SKS_*.png` / `*_SKKS_*.png`（field / misfit / fit / scatter）
- SKS：震央距離配對；SKKS：proxy → **BAZ 配對**（`max_baz_deg=15`）

### 待辦

1. wl 上用已 patch DUMMY **重跑 Job B + D**（SKS/SKKS SI），確認 SYN `paz=0`。
2. 重跑後再出 obs vs model（可先 `MODEL_SI_TO_OBS_SCALE=2`）。
3. 可選：觀測圖也改畫 `0.5·δt·sin2`（與 PSI raw 直接比，則 scale=1）。

---

## wl 上要做什麼（2026-07-20）

Repo: `/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK`

```bash
ssh wl
cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
git pull origin main

# 確認 DUMMY paz=0
awk -F',' 'NR==1{print $NF}' psi_input/DUMMY_SI_Kuo2018_SKS.dat

# 若 0608 專案用自己的 psi_input，也要同步 DUMMY（或從 repo 複製）
# cp psi_input/DUMMY_SI_Kuo2018_SKS.dat  /lfs/wl/.../0608.../psi_input/
# cp psi_input/DUMMY_SI_Kuo2018_SKKS.dat /lfs/wl/.../0608.../psi_input/

# 然後重跑 Job B (Ray SI) + Job D (HFFK) — SKS 與 SKKS
# 重跑後檢查 SYN:
# awk -F',' 'NR==2{print $NF}' psi_output/ray_si_Kuo2018_SKS/SYN_...dat   # 應為 0.0
```

筆記目錄：`validation/notes/`

### 2026-07-20 — Job A/B/C/D 一律用 software `psi_input`

- **規則**：Sources / DUMMY / templates / Receivers / TauP → `${PSI_DIR}/psi_input/`
- **專案只保留**：`viztomo_output/`（MODEL_DAT）與 `psi_output/`
- Templates 改為含 `__TF_GLOBAL_CARTESIAN__` / `__DEPTH_REVERSE__`（VIZTOMO 預設 true / linear）
- 共用：`scripts/psi_project_common.sh`、`scripts/slurm/jobA.sh`–`jobD.sh`
- 專案內 `jobA.sh`–`jobD.sh` 與 long-name jobs 已去掉 `PROJECT_DIR/psi_input` 優先邏輯
