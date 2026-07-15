#!/bin/bash
# job_bench.sh — PSI_D 全套 benchmark（4 模型 × ray SP + HFFK SI T=4-50s）
#
# 前置：
#   git pull && python3 validation/gen_benchmark_psitomo.py
#
# 提交（直接在 software 目錄提交，不需要 project 目錄）：
#   cd /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
#   sbatch validation/job_bench.sh
#
# 輸出全部在：
#   /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/bench_output/
#
#SBATCH -J psi_bench
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -o validation/bench_%j.out
#SBATCH -e validation/bench_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_v2.toml"
SRC="${PSI_DIR}/psi_input/Sources_uniform48.dat"
OBS_SP="${PSI_DIR}/psi_input/DUMMY_SP_uniform48.dat"   # SplittingParameters, period=8s
OBS_SI="${PSI_DIR}/psi_input/DUMMY_SI_uniform48.dat"   # SplittingIntensity,  period=0 (from toml)
BENCH_DIR="${PSI_DIR}/validation/bench_models"
OUT_BASE="${PSI_DIR}/validation/bench_output"
PERIODS=(4.0 8.0 16.0 20.0 25.0 33.0 50.0)
MODELS=(bench_1L_A bench_1L_B bench_2L_A bench_2L_B)
export JULIA_COPY_STACKS=1

# 48 sources × 525 receivers + 1 header
MIN_LINES=25201

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
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

mkdir -p "${OUT_BASE}"

echo "===== Benchmark start $(date) ====="
echo "PSI_DIR: ${PSI_DIR}"
echo "Models:  ${BENCH_DIR}"
echo "Output:  ${OUT_BASE}"

for MOD in "${MODELS[@]}"; do
    MODEL_FILE="${BENCH_DIR}/${MOD}.dat"
    if [ ! -f "${MODEL_FILE}" ]; then
        echo "[WARN] ${MODEL_FILE} not found — skipping ${MOD}"; continue
    fi

    # Ray SP (high-freq approx, tf_hffk_sp=false, period=8s in DUMMY_SP)
    run_psi "bench_${MOD}_ray" "${MODEL_FILE}" \
            "${OBS_SP}" "SplittingParameters" "8.0" "tf_hffk_sp = false"

    # HFFK SI 各頻率
    for T in "${PERIODS[@]}"; do
        TT="${T%.*}s"
        run_psi "bench_${MOD}_hffk_T${TT}" "${MODEL_FILE}" \
                "${OBS_SI}" "SplittingIntensity" "${T}" ""
    done
done

echo ""; echo "===== All benchmark done $(date) ====="
