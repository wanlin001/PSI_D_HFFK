#!/bin/bash
# Job D — HFFK SplittingIntensity, Kuo2018 SKS + SKKS, T=4/8/16/20/25/33/50 s
# psi_input: always software PSI_DIR. Output: project psi_output/
# Submit : cd <PROJECT> && sbatch jobD.sh
#SBATCH -J psi_D_hffk_kuo
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH -o D_hffk_Kuo2018_%j.out
#SBATCH -e D_hffk_Kuo2018_%j.err

set -e
PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
source "${PSI_DIR}/scripts/psi_project_common.sh"
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_v2.toml"
SRC_SKS="${PSI_DIR}/psi_input/Sources_Kuo2018_SKS.dat"
OBS_SKS="${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKS.dat"
SRC_SKKS="${PSI_DIR}/psi_input/Sources_Kuo2018_SKKS.dat"
OBS_SKKS="${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKKS.dat"
export JULIA_COPY_STACKS=1
cd "${PROJECT_DIR}"

run_hffk() {
    local TAG="$1" PERIOD="$2" SRC="$3" OBS="$4" MIN_LINES="$5"
    local OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
    local TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"
    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG} already complete"; return 0
    fi
    echo ""; echo "===== ${TAG} T=${PERIOD}s n_rings=${N_RINGS} n_azimuth=${N_AZIMUTH} $(date) ====="
    echo "  SRC=${SRC}"; echo "  OBS=${OBS}"
    mkdir -p "${OUTDIR}" "${TMPDIR}"
    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_DAT}|g" \
        -e "s|__OBS_TYPE__|SplittingIntensity|g" \
        -e "s|__OBS_FILE__|${OBS}|g" \
        -e "s|__HFFK_SP__||g" \
        -e "s|__PERIOD__|${PERIOD}|g" \
        -e "s|__N_RINGS__|${N_RINGS}|g" \
        -e "s|__N_AZIMUTH__|${N_AZIMUTH}|g" \
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
for T in 4.0 8.0 16.0 20.0 25.0 33.0 50.0; do
    run_hffk "hffk_Kuo2018_SKS_T${T%.*}s"  "$T" "${SRC_SKS}"  "${OBS_SKS}"  2101
    run_hffk "hffk_Kuo2018_SKKS_T${T%.*}s" "$T" "${SRC_SKKS}" "${OBS_SKKS}" 7876
done
echo "===== Job D done $(date) ====="
