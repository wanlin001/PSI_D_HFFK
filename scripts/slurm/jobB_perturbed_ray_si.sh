#!/bin/bash
# Job B-pert — Ray SI, Kuo2018 SKS perturbed catalogue (±BAZ × ±depth)
# Prerequisite: python gen_Sources_Kuo2018_perturbed.py  (in plot_sks_script)
# Output : psi_output/ray_si_Kuo2018_SKS_perturbed/
# Submit : cd <PROJECT> && sbatch jobB_perturbed_ray_si.sh
#SBATCH -J psi_Bpert_kuo
#SBATCH -p 8358
#SBATCH --time=08:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o Bpert_ray_si_%j.out
#SBATCH -e Bpert_ray_si_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
source "${PSI_DIR}/scripts/psi_project_common.sh"
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_ray_si.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_filled.dat"

SRC="${PSI_DIR}/psi_input/Sources_Kuo2018_SKS_perturbed.dat"
OBS="${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKS_perturbed.dat"
TAG="ray_si_Kuo2018_SKS_perturbed"
# 60 src × 525 rcv + header ≈ 31501 (4×5×3=60 sources)
MIN_LINES=31501

export JULIA_COPY_STACKS=1
cd "${PROJECT_DIR}"

OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"

if [ ! -f "${SRC}" ]; then
  echo "MISSING ${SRC} — run gen_Sources_Kuo2018_perturbed.py first"
  exit 1
fi

echo "===== ${TAG} $(date) ====="
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
"; \
  cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
wc -l "${OUTDAT}"
echo "===== ${TAG} done ====="
