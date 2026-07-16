#!/bin/bash
# job_bench.sh — PSI_D 全套 benchmark
#   4 均勻模型 (1L/2L × coarse/fine) × ray SP + ray SI + HFFK SI T=4-50s
#   2 橫向 period 模型 (1L_lat_B, 2L_lat_B) × ray SP + ray SI + HFFK SI T=4-50s
#   2 橫向模型 (lateral_B/C)         × ray SP + ray SI + HFFK SI T=4-50s
#   ray SI = TauP 射線 + SI kernel 積分（無限高頻極限，與 HFFK 同 observable）
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
TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench.toml"
RAY_SI_TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench_ray_si.toml"
REFMODEL="${PSI_DIR}/psi_input/ak135_no_crust.tvel"
BENCH_INP="${PSI_DIR}/validation/bench_psi_input"
SRC="${BENCH_INP}/Sources.dat"
OBS_SP="${BENCH_INP}/DUMMY_SP.dat"
OBS_SI="${BENCH_INP}/DUMMY_SI.dat"
BENCH_DIR="${PSI_DIR}/validation/bench_models"
OUT_BASE="${SLURM_SUBMIT_DIR}/bench_output"

PERIODS=(4.0 8.0 16.0 20.0 25.0 33.0 50.0)

UNIFORM_MODELS=(bench_1L_A bench_1L_B bench_2L_A bench_2L_B)
PERIOD_TEST_MODELS=(bench_1L_lat_B bench_2L_lat_B)
LATERAL_MODELS=(bench_lateral_B bench_lateral_C bench_lateral_D)

MIN_LINES=145
export JULIA_COPY_STACKS=1

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
        -e "s|__DEPTH_REVERSE__|dims3|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    sed -i "s|^receiver_data.*|receiver_data = \"${BENCH_INP}/Receivers.dat\"|g" \
        "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

run_ray_si() {
    local TAG="$1" MODEL_FILE="$2"
    local OUTDIR="${OUT_BASE}/${TAG}"
    local TMPDIR="${OUT_BASE}/.tmp_${TAG}_$$"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG}"; return 0
    fi
    rm -rf "${OUTDIR}"
    echo ""; echo "===== ${TAG}  Ray SI (inf. freq.)  $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_FILE}|g" \
        -e "s|__OBS_FILE__|${OBS_SI}|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__RECEIVERS_DAT__|${BENCH_INP}/Receivers.dat|g" \
        -e "s|__REFMODEL__|${REFMODEL}|g" \
        -e "s|__DEPTH_REVERSE__|dims3|g" \
        "${RAY_SI_TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    sed -i "s|^receiver_data.*|receiver_data = \"${BENCH_INP}/Receivers.dat\"|g" \
        "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

run_all_periods() {
    local MOD="$1" MODEL_FILE="$2"
    run_psi "${MOD}_ray" "${MODEL_FILE}" \
            "${OBS_SP}" "SplittingParameters" "8.0" "tf_hffk_sp = false"
    run_ray_si "${MOD}_ray_si" "${MODEL_FILE}"
    for T in "${PERIODS[@]}"; do
        run_psi "${MOD}_hffk_T${T%.*}s" "${MODEL_FILE}" \
                "${OBS_SI}" "SplittingIntensity" "${T}" ""
    done
}

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
echo "--- Col B/C period 診斷（含水平 φ 變化）---"
for MOD in "${PERIOD_TEST_MODELS[@]}"; do
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
