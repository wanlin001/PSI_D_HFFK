#!/bin/bash
# Job B2 — Ray theory SplittingIntensity, Kuo2018 events SKS+SKKS (no HFFK)
# Output : psi_output/ray_si_Kuo2018_SKS/
#          psi_output/ray_si_Kuo2018_SKKS/
# Paired with jobB (SP) and jobD (HFFK SI) for direct comparison.
# Submit : cd /home/wl/work/ASPECT/<PROJECT> && qsub jobB2_ray_si_Kuo2018.sh
#PBS -N psi_B2_raysi_kuo
#PBS -q q25g
#PBS -l walltime=04:00:00
#PBS -l nodes=1:ncpus=4
#PBS -l mem=16gb
#PBS -o B2_ray_si_Kuo2018.out
#PBS -e B2_ray_si_Kuo2018.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${PBS_O_WORKDIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_ray_si.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_fix2.dat"

cd "${PROJECT_DIR}"

run_ray_si() {
    local TAG="$1" SRC="$2" OBS="$3" MIN_LINES="$4"
    local OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
    local OUTDAT="${OUTDIR}/SYN_SplittingIntensity_ShearWave.dat"
    local TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG} already complete ($(wc -l < "${OUTDAT}") lines)"; return 0
    fi

    echo ""; echo "===== ${TAG} $(date) ====="
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
    wc -l "${OUTDAT}"; echo "===== ${TAG} done ====="
}

# SKS: 4 events → 4 × 525 + 1 = 2101 lines
run_ray_si "ray_si_Kuo2018_SKS" \
    "${PSI_DIR}/psi_input/Sources_Kuo2018_SKS.dat" \
    "${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKS.dat" \
    2101

# SKKS: 15 proxy events → 15 × 525 + 1 = 7876 lines
run_ray_si "ray_si_Kuo2018_SKKS" \
    "${PSI_DIR}/psi_input/Sources_Kuo2018_SKKS.dat" \
    "${PSI_DIR}/psi_input/DUMMY_SI_Kuo2018_SKKS.dat" \
    7876

echo ""; echo "===== Job B2 done $(date) ====="
