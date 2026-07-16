#!/bin/bash
# job_bench_nazimuth.sh — HFFK n_azimuth 收斂測試 (T2-3)
#   模型：bench_lateral_B
#   Period：T=25s
#   n_rings：3（固定）
#   n_azimuth：4, 6, 8, 12, 16, 24, 32, 48
#   參考：n_azimuth=48
#
# 提交：
#   cd /lfs/wl/bench_psi
#   sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench_nazimuth.sh
#
# 輸出：/lfs/wl/bench_psi/bench_output/nazimuth_conv/
#
#SBATCH -J psi_naz
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o nazimuth_%j.out
#SBATCH -e nazimuth_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench.toml"
REFMODEL="${PSI_DIR}/psi_input/ak135_no_crust.tvel"
BENCH_INP="${PSI_DIR}/validation/bench_psi_input"
SRC="${BENCH_INP}/Sources.dat"
OBS_SI="${BENCH_INP}/DUMMY_SI.dat"
MODEL_FILE="${PSI_DIR}/validation/bench_models/bench_lateral_B.dat"
OUT_BASE="${SLURM_SUBMIT_DIR}/bench_output/nazimuth_conv"

PERIOD=25.0
N_RINGS=3
NAZ_LIST=(4 6 8 10 12 16 20 24 32 40 48 64)
MIN_LINES=145
export JULIA_COPY_STACKS=1

run_hffk_nazimuth() {
    local NAZ="$1"
    local TAG="lateral_B_hffk_T25s_az${NAZ}"
    local OUTDIR="${OUT_BASE}/${TAG}"
    local TMPDIR="${OUT_BASE}/.tmp_${TAG}_$$"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG}"; return 0
    fi
    rm -rf "${OUTDIR}"
    echo ""; echo "===== ${TAG}  T=${PERIOD}s  n_azimuth=${NAZ}  $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_FILE}|g" \
        -e "s|__OBS_TYPE__|SplittingIntensity|g" \
        -e "s|__OBS_FILE__|${OBS_SI}|g" \
        -e "s|__HFFK_SP__||g" \
        -e "s|__PERIOD__|${PERIOD}|g" \
        -e "s|__N_RINGS__|${N_RINGS}|g" \
        -e "s|__N_AZIMUTH__|${NAZ}|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__RECEIVERS_DAT__|${BENCH_INP}/Receivers.dat|g" \
        -e "s|__REFMODEL__|${REFMODEL}|g" \
        -e "s|__DEPTH_REVERSE__|linear|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    sed -i "s|^receiver_data.*|receiver_data = \"${BENCH_INP}/Receivers.dat\"|g" \
        "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

mkdir -p "${OUT_BASE}"
echo "===== n_azimuth convergence start $(date) ====="
echo "Model    : bench_lateral_B"
echo "Period   : T=${PERIOD}s"
echo "n_rings  : ${N_RINGS}"
echo "Output   : ${OUT_BASE}"

for NAZ in "${NAZ_LIST[@]}"; do
    run_hffk_nazimuth "${NAZ}"
done

echo ""; echo "===== n_azimuth convergence done $(date) ====="

REF="${OUT_BASE}/lateral_B_hffk_T25s_az48/SYN_SplittingIntensity_ShearWave.dat"
if [ -f "${REF}" ]; then
    echo ""
    echo "=== Convergence vs n_azimuth=48 (RMS diff, max diff) ==="
    for NAZ in "${NAZ_LIST[@]}"; do
        F="${OUT_BASE}/lateral_B_hffk_T25s_az${NAZ}/SYN_SplittingIntensity_ShearWave.dat"
        if [ -f "${F}" ]; then
            paste <(awk -F',' '{print $1}' "${REF}") <(awk -F',' '{print $1}' "${F}") | \
                awk '{d=$1-$2; if(d<0)d=-d; s+=d*d; if(d>m)m=d; c++} END{printf "n_azimuth=%s  n=%d  RMS=%.6f  max=%.6f\n", naz, c, sqrt(s/c), m}' naz="${NAZ}"
        fi
    done
fi
