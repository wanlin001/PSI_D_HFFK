#!/bin/bash
# job_bench.sh — PSI_D 全套 benchmark
#   4 均勻模型 (1L/2L × coarse/fine) × ray SP + HFFK SI T=4-50s
#   2 橫向模型 (lateral_B/C)         × ray SP + HFFK SI T=4-50s
#
# 提交（必須從 /lfs 下的目錄提交）：
#   mkdir -p /lfs/wl/bench_psi
#   cd /lfs/wl/bench_psi
#   sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench.sh
#
# 輸出在：/lfs/wl/bench_psi/bench_output/
#
#SBATCH -J psi_bench
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=72:00:00
#SBATCH -o bench_%j.out
#SBATCH -e bench_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench.toml"  # 自足，含 __RECEIVERS_DAT__/__REFMODEL__ placeholder
REFMODEL="${PSI_DIR}/psi_input/ak135_no_crust.tvel"
BENCH_INP="${PSI_DIR}/validation/bench_psi_input"   # 獨立 benchmark 測站/震源
SRC="${BENCH_INP}/Sources.dat"
OBS_SP="${BENCH_INP}/DUMMY_SP.dat"
OBS_SI="${BENCH_INP}/DUMMY_SI.dat"
BENCH_DIR="${PSI_DIR}/validation/bench_models"
OUT_BASE="${SLURM_SUBMIT_DIR}/bench_output"   # /lfs/wl/bench_psi/bench_output

PERIODS=(4.0 8.0 16.0 20.0 25.0 33.0 50.0)

# 均勻模型（1L/2L × A/B）
UNIFORM_MODELS=(bench_1L_A bench_1L_B bench_2L_A bench_2L_B)

# 橫向非均勻模型（HFFK period 差異診斷）
LATERAL_MODELS=(bench_lateral_B bench_lateral_C)

MIN_LINES=145   # 16 sources × 9 receivers + 1
export JULIA_COPY_STACKS=1

# ── 函式 ─────────────────────────────────────────────────────────
run_psi() {
    local TAG="$1" MODEL_FILE="$2" OBS="$3" OBS_TYPE="$4" PERIOD="$5" HFFK_SP="$6"
    local OUTDIR="${OUT_BASE}/${TAG}"
    local TMPDIR="${OUT_BASE}/.tmp_${TAG}_$$"

    if [[ "${OBS_TYPE}" == "SplittingParameters" ]]; then
        OUTDAT="${OUTDIR}/SYN_SplittingParameters_ShearWave.dat"
    else
        OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
    fi

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG}"; return 0
    fi
    rm -rf "${OUTDIR}"
    echo ""; echo "===== ${TAG}  T=${PERIOD}s  $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_FILE}|g" \
        -e "s|__OBS_TYPE__|${OBS_TYPE}|g" \
        -e "s|__OBS_FILE__|${OBS}|g" \
        -e "s|__HFFK_SP__|${HFFK_SP}|g" \
        -e "s|__PERIOD__|${PERIOD}|g" \
        -e "s|__N_RINGS__|3|g" \
        -e "s|__N_AZIMUTH__|8|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__RECEIVERS_DAT__|${BENCH_INP}/Receivers.dat|g" \
        -e "s|__REFMODEL__|${REFMODEL}|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    # 保險：若曾用寫死 receiver_data 的舊 template，也強制指向 benchmark 測站
    # （benchmark 模型格網 119-127°E；真實 TW0001 在 116.5°E 會出界）。
    sed -i "s|^receiver_data.*|receiver_data = \"${BENCH_INP}/Receivers.dat\"|g" \
        "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

run_all_periods() {
    local MOD="$1" MODEL_FILE="$2"
    # Ray SP (high-freq approx)
    run_psi "${MOD}_ray" "${MODEL_FILE}" \
            "${OBS_SP}" "SplittingParameters" "8.0" "tf_hffk_sp = false"
    # HFFK SI 各頻率
    for T in "${PERIODS[@]}"; do
        run_psi "${MOD}_hffk_T${T%.*}s" "${MODEL_FILE}" \
                "${OBS_SI}" "SplittingIntensity" "${T}" ""
    done
}

# ── 主程式 ──────────────────────────────────────────────────────
mkdir -p "${OUT_BASE}"
echo "===== Benchmark start $(date) ====="
echo "PSI_DIR : ${PSI_DIR}"
echo "Models  : ${BENCH_DIR}"
echo "Output  : ${OUT_BASE}"

echo ""
echo "--- 均勻模型（resolution + 雙層干涉）---"
for MOD in "${UNIFORM_MODELS[@]}"; do
    MODEL_FILE="${BENCH_DIR}/${MOD}.dat"
    [ -f "${MODEL_FILE}" ] || { echo "[WARN] ${MODEL_FILE} not found"; continue; }
    run_all_periods "${MOD}" "${MODEL_FILE}"
done

echo ""
echo "--- 橫向模型（HFFK period 差異診斷）---"
for MOD in "${LATERAL_MODELS[@]}"; do
    MODEL_FILE="${BENCH_DIR}/${MOD}.dat"
    [ -f "${MODEL_FILE}" ] || { echo "[WARN] ${MODEL_FILE} not found"; continue; }
    run_all_periods "${MOD}" "${MODEL_FILE}"
done

echo ""; echo "===== All benchmark done $(date) ====="
