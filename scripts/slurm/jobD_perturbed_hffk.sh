#!/bin/bash
# Job D-pert — HFFK SI, Kuo2018 SKS perturbed (all DEFAULT periods)
# Prerequisite: gen_Sources_Kuo2018_perturbed.py
# Output : psi_output/hffk_Kuo2018_SKS_perturbed_T{4,8,16,20,25,33,50}s/
# Submit : cd <PROJECT> && sbatch jobD_perturbed_hffk.sh
#SBATCH -J psi_Dpert_kuo
#SBATCH -p 8358
#SBATCH --time=48:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o Dpert_hffk_%j.out
#SBATCH -e Dpert_hffk_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
source "${PSI_DIR}/scripts/psi_project_common.sh"
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_v2.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_filled.dat"

SRC="${PSI_DIR}/psi_input/Sources_Kuo2018_SKS_perturbed.dat"
OBS="${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKS_perturbed.dat"
MIN_LINES=31501

export JULIA_COPY_STACKS=1
cd "${PROJECT_DIR}"

run_hffk() {
    local TAG="$1" PERIOD="$2"
    local OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
    local TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"
    echo ""; echo "===== ${TAG} T=${PERIOD}s $(date) ====="
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
"; \
      cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
    wc -l "${OUTDAT}"
}

for T in 4 8 16 20 25 33 50; do
  run_hffk "hffk_Kuo2018_SKS_perturbed_T${T}s" "${T}.0"
done
echo ""; echo "===== Job D-pert done $(date) ====="
