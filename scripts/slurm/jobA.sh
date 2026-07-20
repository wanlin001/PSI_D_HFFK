#!/bin/bash
# Job A — Ray SplittingIntensity, synthetic uniform48 events
# psi_input: always software PSI_DIR (shared). Output: project psi_output/
# Submit : cd <PROJECT> && sbatch jobA.sh
#SBATCH -J psi_A_raysi_u48
#SBATCH -p 8358
#SBATCH --time=06:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o A_ray_si_uniform48_%j.out
#SBATCH -e A_ray_si_uniform48_%j.err

set -e
PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
source "${PSI_DIR}/scripts/psi_project_common.sh"
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_ray_si.toml"
SRC="${PSI_DIR}/psi_input/Sources_uniform48.dat"
OBS="${PSI_DIR}/psi_input/DUMMY_SI_uniform48.dat"
TAG="ray_si_uniform48"
OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"
MIN_LINES=25201

cd "${PROJECT_DIR}"
if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
    echo "[SKIP] ${TAG} already complete"; exit 0
fi
echo "===== ${TAG} n_rings=${N_RINGS} n_azimuth=${N_AZIMUTH} $(date) ====="
echo "  TEMPLATE=${TEMPLATE}"
echo "  SRC=${SRC}  OBS=${OBS}"
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
echo "===== ${TAG} done $(date) ====="
