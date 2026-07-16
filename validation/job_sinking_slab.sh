#!/bin/bash
# job_sinking_slab.sh — T3-3 SinkingSlab 官方範例驗證
#   Ray SI（TauP）+ HFFK T=3,5,8,15,25 s
#   參考：ECOMAN2.0-seismology.PSI_D/examples/SinkingSlab/psi_output/SYN_SinkingBlock/
#
# 提交：
#   cd /lfs/wl/bench_psi
#   sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/validation/job_sinking_slab.sh
#
# 輸出：/lfs/wl/bench_psi/bench_output/sinking_slab/
#
#SBATCH -J psi_sink
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH -o sinking_%j.out
#SBATCH -e sinking_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
JULIA=/home/wl/software/julia-1.10.0/bin/julia
SINK_INP="/home/wl/software/ECOMAN2.0-seismology.PSI_D/examples/SinkingSlab/psi_input"
RAY_SI_TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench_ray_si.toml"
HFFK_TEMPLATE="${PSI_DIR}/validation/psi_config_template_bench.toml"
MODEL_FILE="${SINK_INP}/psitomo0020.dat"
OBS_SI="${SINK_INP}/DUMMY_SplittingIntensity_ShearWave.dat"
SRC="${SINK_INP}/Sources.dat"
RCV="${SINK_INP}/Receivers.dat"
REFMODEL="${SINK_INP}/ak135_no_crust.tvel"
OUT_BASE="${SLURM_SUBMIT_DIR}/bench_output/sinking_slab"

PERIODS=(3.0 5.0 8.0 15.0 25.0)
MIN_LINES=14000
export JULIA_COPY_STACKS=1

run_ray_si() {
    local TAG="SinkingSlab_ray_si"
    local OUTDIR="${OUT_BASE}/${TAG}"
    local TMPDIR="${OUT_BASE}/.tmp_${TAG}_$$"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG}"; return 0
    fi
    rm -rf "${OUTDIR}"
    echo ""; echo "===== ${TAG}  Ray SI  $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_FILE}|g" \
        -e "s|__OBS_FILE__|${OBS_SI}|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__RECEIVERS_DAT__|${RCV}|g" \
        -e "s|__REFMODEL__|${REFMODEL}|g" \
        -e "s|__DEPTH_REVERSE__|linear|g" \
        -e "s|tf_global_cartesian = false|tf_global_cartesian = true|g" \
        "${RAY_SI_TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

run_hffk() {
    local T="$1"
    local TAG="SinkingSlab_hffk_T${T%.*}s"
    local OUTDIR="${OUT_BASE}/${TAG}"
    local TMPDIR="${OUT_BASE}/.tmp_${TAG}_$$"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG}"; return 0
    fi
    rm -rf "${OUTDIR}"
    echo ""; echo "===== ${TAG}  T=${T}s  $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_FILE}|g" \
        -e "s|__OBS_TYPE__|SplittingIntensity|g" \
        -e "s|__OBS_FILE__|${OBS_SI}|g" \
        -e "s|__HFFK_SP__||g" \
        -e "s|__PERIOD__|${T}|g" \
        -e "s|__N_RINGS__|3|g" \
        -e "s|__N_AZIMUTH__|8|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__RECEIVERS_DAT__|${RCV}|g" \
        -e "s|__REFMODEL__|${REFMODEL}|g" \
        -e "s|__DEPTH_REVERSE__|linear|g" \
        -e "s|tf_global_cartesian = false|tf_global_cartesian = true|g" \
        "${HFFK_TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    echo "===== ${TAG} done ====="
}

mkdir -p "${OUT_BASE}"
echo "===== SinkingSlab T3-3 start $(date) ====="
echo "Input   : ${SINK_INP}"
echo "Output  : ${OUT_BASE}"

run_ray_si
for T in "${PERIODS[@]}"; do
    run_hffk "${T}"
done

echo ""; echo "===== SinkingSlab T3-3 done $(date) ====="

# 快速比對官方 Ray 輸出
OFFICIAL="/home/wl/software/ECOMAN2.0-seismology.PSI_D/examples/SinkingSlab/psi_output/SYN_SinkingBlock/SYN_SplittingIntensity_ShearWave.dat"
OURS="${OUT_BASE}/SinkingSlab_ray_si/SYN_SplittingIntensity_ShearWave.dat"
if [ -f "${OFFICIAL}" ] && [ -f "${OURS}" ]; then
    echo ""
    echo "=== Ray SI vs Official PSI_D SinkingSlab ==="
    paste <(awk -F',' '{print $1}' "${OFFICIAL}") <(awk -F',' '{print $1}' "${OURS}") | \
        awk '{d=$1-$2; if(d<0)d=-d; s+=d*d; if(d>m)m=d; c++} END{printf "n=%d  RMS=%.6e  max=%.6e\n", c, sqrt(s/c), m}'
fi
