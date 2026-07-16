#!/bin/bash
# job_bench_nrings.sh — HFFK n_rings 收斂測試
#   模型：bench_lateral_B（橫向邊界，period 效應最明顯）
#   Period：T=25s
#   n_rings：1, 2, 3, 4, 5, 6, 8, 12, 16（n_azimuth 固定 8）
#   參考：n_rings=16
#
# 提交：
#   cd /lfs/wl/bench_psi
#   sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_bench_nrings.sh
#
# 輸出：/lfs/wl/bench_psi/bench_output/nrings_conv/
#
#SBATCH -J psi_nrings
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH -o nrings_%j.out
#SBATCH -e nrings_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench.toml"
REFMODEL="${PSI_DIR}/psi_input/ak135_no_crust.tvel"
BENCH_INP="${PSI_DIR}/validation/bench_psi_input"
SRC="${BENCH_INP}/Sources.dat"
OBS_SI="${BENCH_INP}/DUMMY_SI.dat"
MODEL_FILE="${PSI_DIR}/validation/bench_models/bench_lateral_B.dat"
OUT_BASE="${SLURM_SUBMIT_DIR}/bench_output/nrings_conv"

PERIOD=25.0
N_AZIMUTH=8
NRINGS_LIST=(1 2 3 4 5 6 8 12 16)
MIN_LINES=145
export JULIA_COPY_STACKS=1

run_hffk_nrings() {
    local NR="$1"
    local TAG="lateral_B_hffk_T25s_rings${NR}"
    local OUTDIR="${OUT_BASE}/${TAG}"
    local TMPDIR="${OUT_BASE}/.tmp_${TAG}_$$"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG}"; return 0
    fi
    rm -rf "${OUTDIR}"
    echo ""; echo "===== ${TAG}  T=${PERIOD}s  n_rings=${NR}  $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_FILE}|g" \
        -e "s|__OBS_TYPE__|SplittingIntensity|g" \
        -e "s|__OBS_FILE__|${OBS_SI}|g" \
        -e "s|__HFFK_SP__||g" \
        -e "s|__PERIOD__|${PERIOD}|g" \
        -e "s|__N_RINGS__|${NR}|g" \
        -e "s|__N_AZIMUTH__|${N_AZIMUTH}|g" \
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

mkdir -p "${OUT_BASE}"
echo "===== n_rings convergence start $(date) ====="
echo "Model   : bench_lateral_B"
echo "Period  : T=${PERIOD}s"
echo "n_azimuth: ${N_AZIMUTH}"
echo "Output  : ${OUT_BASE}"

for NR in "${NRINGS_LIST[@]}"; do
    run_hffk_nrings "${NR}"
done

echo ""; echo "===== n_rings convergence done $(date) ====="

# 內建快速統計（不需 numpy）
REF="${OUT_BASE}/lateral_B_hffk_T25s_rings16/SYN_SplittingIntensity_ShearWave.dat"
if [ -f "${REF}" ]; then
    echo ""
    echo "=== Convergence vs n_rings=16 (RMS diff, max diff) ==="
    for NR in "${NRINGS_LIST[@]}"; do
        F="${OUT_BASE}/lateral_B_hffk_T25s_rings${NR}/SYN_SplittingIntensity_ShearWave.dat"
        if [ -f "${F}" ]; then
            paste <(awk -F',' '{print $1}' "${REF}") <(awk -F',' '{print $1}' "${F}") | \
                awk '{d=$1-$2; if(d<0)d=-d; s+=d*d; if(d>m)m=d; c++} END{printf "n_rings=%s  n=%d  RMS=%.6f  max=%.6f\n", nr, c, sqrt(s/c), m}' nr="${NR}"
        fi
    done
fi
