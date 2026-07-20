#!/bin/bash
# Job A2 — Ray theory SplittingIntensity, uniform48 (no HFFK)
# Output : psi_output/ray_si_uniform48/
# Paired with jobA (SP) and jobC (HFFK SI) for direct comparison.
# Submit : cd /home/wl/work/ASPECT/<PROJECT> && qsub jobA2_ray_si_uniform48.sh
#PBS -N psi_A2_raysi_u48
#PBS -q q25g
#PBS -l walltime=06:00:00
#PBS -l nodes=1:ncpus=4
#PBS -l mem=16gb
#PBS -o A2_ray_si_uniform48.out
#PBS -e A2_ray_si_uniform48.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${PBS_O_WORKDIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_ray_si.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_fix2.dat"
SRC="${PSI_DIR}/psi_input/Sources_uniform48.dat"
OBS="${PSI_DIR}/psi_input/DUMMY_SI_uniform48.dat"

TAG="ray_si_uniform48"
OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"
MIN_LINES=25201   # 48 × 525 + 1

cd "${PROJECT_DIR}"

if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
    echo "[SKIP] ${TAG} already complete ($(wc -l < "${OUTDAT}") lines)"; exit 0
fi

echo "===== ${TAG} $(date) ====="
mkdir -p "${OUTDIR}" "${TMPDIR}"

sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
    -e "s|__MODEL_DAT__|${MODEL_DAT}|g" \
    -e "s|__OBS_FILE__|${OBS}|g" \
    -e "s|__SOURCES_DAT__|${SRC}|g" \
    -e "s|__TF_GLOBAL_CARTESIAN__|${TF_GLOBAL_CARTESIAN:-true}|g" \
    -e "s|__DEPTH_REVERSE__|${DEPTH_REVERSE:-linear}|g" \
    "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
rm -rf "${TMPDIR}"

{ printf "# SI, error, period_s, phase, src_id, rcv_id, unknown, paz_rad
"; \
  cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
wc -l "${OUTDAT}"
echo "===== ${TAG} done $(date) ====="
