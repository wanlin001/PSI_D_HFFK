# PSI_D_HFFK — Operations Guide

Last checked: 2026-07-22  
Canonical software tree on HPC (`wl`):

```text
/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
```

GitHub remote: `git@github.com:wanlin001/PSI_D_HFFK.git`  
Branch of record: `main`

---

## 1. Layout (what lives where)

| Location | Contents |
|----------|----------|
| `PSI_DIR` = software clone above | Julia package, **shared** `psi_input/`, SLURM/PBS jobs, validation |
| `PSI_DIR/psi_input/` | Sources, DUMMY obs, Receivers, TauP, TOML templates (**do not copy into each ASPECT project**) |
| `PSI_DIR/scripts/slurm/` | Canonical SI jobs only (`jobA.sh` … `jobD_perturbed_*.sh`) |
| `PSI_DIR/scripts/slurm/bk/` | Archived duplicate / SP-only scripts (do not submit) |
| `PSI_DIR/validation/` | Benchmark SLURM jobs + plotting scripts |
| ASPECT **project** folder | `viztomo_output/psitomo0050_filled.dat` (model) + `psi_output/` (results only) |

Shared inputs sit next to the code (`psi_input/` beside `src` / `scripts`), not inside every model directory.

Default model file resolved by `scripts/psi_project_common.sh`:

```text
${PROJECT_DIR}/viztomo_output/psitomo${MODEL_STEP}_filled.dat
# MODEL_STEP defaults to 0050
```

---

## 2. Keep Mac and `wl` git identical

On Mac (or any clone that should lead):

```bash
cd /path/to/PSI_D_HFFK
git checkout main
git pull origin main
git status   # should be clean on main, matching origin/main
```

On `wl`:

```bash
ssh wl
cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
git fetch origin
git branch backup/wl-main-$(date +%Y%m%d)   # safety backup of current tip
git reset --hard origin/main
git status
git rev-parse HEAD   # must equal Mac / origin/main
```

Do **not** keep divergent uncommitted job copies inside ASPECT projects as the source of truth. Always refresh from:

```bash
cp /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/scripts/slurm/jobA.sh .
# … etc
```

---

## 3. Shared `psi_input` catalogue

Path: `${PSI_DIR}/psi_input/`

### Synthetic event catalogues

| File | Events | Notes |
|------|--------|-------|
| `Sources_uniform48.dat` | 48 | 16 BAZ × 3 depths |
| `Sources_uniform96.dat` | 96 | 32 BAZ × 3 depths |
| `DUMMY_SI_uniform48.dat` | — | SI dummy / paz for u48 |
| `DUMMY_SI_uniform96.dat` | — | SI dummy / paz for u96 |
| `Receivers.dat` | TW grid | Shared receivers |
| `ak135_no_crust.tvel` | — | TauP reference |
| `psi_config_template_ray_si.toml` | — | Ray SI template |
| `psi_config_template_v2.toml` | — | HFFK (and SP) template |

### Kuo2018 (real-event geometry)

| File | Role |
|------|------|
| `Sources_Kuo2018_SKS.dat` / `_SKKS.dat` | Event sources |
| `DUMMY_SI_Kuo2018_SKS.dat` / `_SKKS.dat` | SI dummy (**paz must be 0**) |

Check paz:

```bash
awk -F',' 'NR==1{print $NF}' psi_input/DUMMY_SI_Kuo2018_SKS.dat
# expect: 0.0
```

---

## 4. Periods (HFFK)

All HFFK SI jobs use:

```text
T = 4, 8, 16, 20, 25, 33, 50 s
```

Output tags look like `hffk_u48_T8s`, `hffk_u96_T50s`, `hffk_Kuo2018_SKS_T16s`.

---

## 5. Production jobs on an ASPECT model (`wl`)

### Prerequisites

1. ASPECT → D-Rex → viztomo finished  
2. Filled model exists:

```bash
ls -lh $PROJ/viztomo_output/psitomo0050_filled.dat
```

### Copy canonical SLURM scripts into the project

```bash
ssh wl
SOFT=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJ=/lfs/wl/ASPECT/YOUR_MODEL_DIR   # <-- edit

cd "$PROJ"
cp "$SOFT/scripts/slurm"/job{A,A96,B,B_perturbed_ray_si,C,C96,D,D_perturbed_hffk}.sh .
```

### Submit (must `cd` into `$PROJ` first)

```bash
cd "$PROJ"
sbatch jobA.sh                    # ray_si_uniform48
sbatch jobA96.sh                  # ray_si_uniform96
sbatch jobB.sh                    # ray_si_Kuo2018_SKS / SKKS
sbatch jobB_perturbed_ray_si.sh   # ray_si_Kuo2018_SKS_perturbed
sbatch jobC.sh                    # hffk_u48_T{4…50}s
sbatch jobC96.sh                  # hffk_u96_T{4…50}s
sbatch jobD.sh                    # hffk_Kuo2018_*_T{4…50}s
sbatch jobD_perturbed_hffk.sh     # hffk_Kuo2018_SKS_perturbed_T*
```

Jobs are independent. Partition: `8358`.

### Expected outputs under `$PROJ/psi_output/`

| Job | Tag(s) |
|-----|--------|
| A | `ray_si_uniform48/` |
| A96 | `ray_si_uniform96/` |
| B | `ray_si_Kuo2018_SKS/`, `ray_si_Kuo2018_SKKS/` |
| B-pert | `ray_si_Kuo2018_SKS_perturbed/` |
| C | `hffk_u48_T4s` … `hffk_u48_T50s` |
| C96 | `hffk_u96_T4s` … `hffk_u96_T50s` |
| D | `hffk_Kuo2018_SKS_T*s/`, `hffk_Kuo2018_SKKS_T*s/` |
| D-pert | `hffk_Kuo2018_SKS_perturbed_T{4…50}s/` |

Each complete SI file: `SYN_SplittingIntensity_ShearWave.dat`.

### Non-default timestep

```bash
export MODEL_STEP=0040   # then sbatch …
# or: export MODEL_DAT=/absolute/path/to/psitomoXXXX_filled.dat
```

### Force re-run

Delete the tag directory (jobs SKIP if line counts look complete):

```bash
rm -rf psi_output/ray_si_uniform48
rm -rf psi_output/hffk_u48_T8s
```

### Receivers for local plotting

```bash
mkdir -p "$PROJ/psi_output"
cp "$SOFT/psi_input/Receivers.dat" "$PROJ/psi_output/"
```

---

## 6. Job cheat sheet (canonical SI only)

Active scripts in `scripts/slurm/` (everything else is under `scripts/slurm/bk/`):

| Script | Catalogue | Method | Periods | Purpose |
|--------|-----------|--------|---------|---------|
| `jobA.sh` | u48 | Ray SI | ∞-freq | Synthetic Ray baseline |
| `jobA96.sh` | u96 | Ray SI | ∞-freq | Denser BAZ synthetic Ray |
| `jobB.sh` | Kuo SKS+SKKS | Ray SI | ∞-freq | Real-event Ray (obs compare) |
| `jobB_perturbed_ray_si.sh` | Kuo SKS perturbed | Ray SI | ∞-freq | BAZ/depth sensitivity (Ray) |
| `jobC.sh` | u48 | HFFK SI | 4–50 s | Synthetic FF vs Job A |
| `jobC96.sh` | u96 | HFFK SI | 4–50 s | Synthetic FF vs Job A96 |
| `jobD.sh` | Kuo SKS+SKKS | HFFK SI | 4–50 s | Real-event FF vs Job B |
| `jobD_perturbed_hffk.sh` | Kuo SKS perturbed | HFFK SI | 4–50 s | FF vs Job B-pert |

Archived duplicates / SP jobs: see `scripts/slurm/bk/README.md`.

---

## 7. Validation benchmarks + plotting

Software path: `validation/`

### Run benchmarks on `wl`

```bash
cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
# or from a bench work dir that points at this tree

sbatch validation/job_bench.sh            # main Ray SI vs HFFK suite
sbatch validation/job_bench_nrings.sh     # n_rings convergence
sbatch validation/job_bench_nazimuth.sh   # n_azimuth convergence
sbatch validation/job_sinking_slab.sh     # optional real-model check
```

Typical bench output (example): `/lfs/wl/bench_psi/bench_output/`

### Plot on `wl` (or Mac after download)

```bash
module load anaconda/2022.05   # on wl, if needed

python3 validation/plot_benchmark_groups.py \
  --bench-out /lfs/wl/bench_psi/bench_output

python3 validation/plot_nrings_convergence.py \
  --conv-dir /lfs/wl/bench_psi/bench_output/nrings_conv

python3 validation/plot_nazimuth_convergence.py \
  --conv-dir /lfs/wl/bench_psi/bench_output/nazimuth_conv

python3 validation/plot_benchmark_tier1.py   # tier-1 panels
python3 validation/plot_benchmark_tier2.py   # tier-2 panels
```

Other helpers:

| Script | Purpose |
|--------|---------|
| `validation/compare_ray_vs_hffk.py` | Pairwise Ray vs HFFK |
| `validation/compare_sinking_slab.py` | Sinking-slab case |
| `validation/plot_si_3d.py` | 3-D SI view |
| `validation/gen_benchmark_psitomo.py` | Build analytic bench models |
| `scripts/plot_real_period_map.py` | Period map on a real model |
| `scripts/gen_dummy_obs.py` | Build / patch DUMMY obs (`--patch` for paz=0) |

Notes / problem log: `validation/notes/` (`NOTES_SI_COMPARE.md`, `PROBLEMS.md`, Tier1/2 docs).

---

## 8. Obs vs model (Kuo2018) on the Mac

Plotting lives under ASPECT scripts (not required inside this repo):

```text
/Users/wanlin/Documents/ASPECT/99_aspect_scripts/plot_sks_script/plotOBS_vs_model_Kuo.py
```

Rules used there:

1. Same event (`src_id`) only  
2. Apply `drift_xy` (km) to model receivers before misfit  
3. Compare **HFFK MIX** = mean over T = 4…50 s  
4. `MODEL_SI_TO_OBS_SCALE = 2` (PSI ½-kernel vs obs δt·sin2)

Edit `MODELS` → set `project` + `drift_xy`, then run in Spyder.

---

## 9. Quick health checks

```bash
# Software HEAD matches GitHub
git rev-parse HEAD
git rev-parse origin/main

# Shared inputs present
ls psi_input/Sources_uniform{48,96}.dat \
   psi_input/DUMMY_SI_uniform{48,96}.dat \
   psi_input/Receivers.dat

# Project model present
ls $PROJ/viztomo_output/psitomo0050_filled.dat

# Job finished (example)
wc -l $PROJ/psi_output/ray_si_uniform48/SYN_SplittingIntensity_ShearWave.dat
# expect ≥ 25201 lines for u48
```

---

## 10. Design rules (do not break)

1. **Always** read Sources / DUMMY / Receivers / templates from `${PSI_DIR}/psi_input/`.  
2. **Never** prefer a stale `PROJECT_DIR/psi_input` for Job A–D.  
3. Project folder owns only `viztomo_output/` (input model) and `psi_output/` (SYN results).  
4. Submit with `cd $PROJ && sbatch jobX.sh` so `SLURM_SUBMIT_DIR` is the project.  
5. For SI comparisons use Ray SI ↔ HFFK SI (same observable), not SP→SI conversions.
