#!/bin/bash
# Job B — Ray SplittingIntensity, Kuo2018 SKS + SKKS
# psi_input: always software PSI_DIR. Output: project psi_output/
# Submit : cd <PROJECT> && sbatch jobB.sh
#SBATCH -J psi_B_raysi_kuo
#SBATCH -p 8358
#SBATCH --time=04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o B_ray_si_Kuo2018_%j.out
#SBATCH -e B_ray_si_Kuo2018_%j.err

set -e
PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
source "${PSI_DIR}/scripts/psi_project_common.sh"
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_ray_si.toml"

cd "${PROJECT_DIR}"
run_ray_si() {
    local TAG="$1" SRC="$2" OBS="$3" MIN_LINES="$4"
    local OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
    local TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"
    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG} already complete"; return 0
    fi
    echo ""; echo "===== ${TAG} $(date) ====="
    echo "  SRC=${SRC}"; echo "  OBS=${OBS}"
    mkdir -p "${OUTDIR}" "${TMPDIR}"
    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_DAT}|g" \
        -e "s|__OBS_FILE__|${OBS}|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__TF_GLOBAL_CARTESIAN__|${TF_GLOBAL_CARTESIAN}|g" \
        -e "s|__DEPTH_REVERSE__|${DEPTH_REVERSE}|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"
    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"
    { printf "# SI, error, period_s, phase, src_id, rcv_id, unknown, paz_rad
"; cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
    wc -l "${OUTDAT}"
}
run_ray_si "ray_si_Kuo2018_SKS"  "${PSI_DIR}/psi_input/Sources_Kuo2018_SKS.dat"  "${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKS.dat"  2101
run_ray_si "ray_si_Kuo2018_SKKS" "${PSI_DIR}/psi_input/Sources_Kuo2018_SKKS.dat" "${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKKS.dat" 7876
echo "===== Job B done $(date) ====="
