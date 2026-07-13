#!/bin/bash
# Job D — HFFK Splitting Intensity, Kuo2018 events (SKS + SKKS), T=4/8/16/20/25/33/50 s
# Output : psi_output/hffk_Kuo2018_SKS_T{4,8,16,20,25,33,50}s/
#          psi_output/hffk_Kuo2018_SKKS_T{4,8,16,20,25,33,50}s/
# Submit : cd /home/wl/work/ASPECT/<PROJECT> && sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/scripts/slurm/jobD_hffk_Kuo2018.sh
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
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_v2.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_fix2.dat"

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
        echo "[SKIP] ${TAG} already complete ($(wc -l < "${OUTDAT}") lines)"; return 0
    fi

    echo ""; echo "===== ${TAG} T=${PERIOD}s $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"
    rm -f "${OUTDIR}/psi_config.toml"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_DAT}|g" \
        -e "s|__OBS_TYPE__|SplittingIntensity|g" \
        -e "s|__OBS_FILE__|${OBS}|g" \
        -e "s|__HFFK_SP__||g" \
        -e "s|__PERIOD__|${PERIOD}|g" \
        -e "s|__N_RINGS__|3|g" \
        -e "s|__N_AZIMUTH__|8|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"

    { printf "# SI, error, period_s, phase, src_id, rcv_id, unknown, paz_rad\n"; \
      cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
    wc -l "${OUTDAT}"; echo "===== ${TAG} done ====="
}

# SKS: 4 real events → 4 × 525 + 1 = 2101 lines
run_hffk "hffk_Kuo2018_SKS_T4s"  4.0  "${SRC_SKS}"  "${OBS_SKS}"  2101
run_hffk "hffk_Kuo2018_SKS_T8s"  8.0  "${SRC_SKS}"  "${OBS_SKS}"  2101
run_hffk "hffk_Kuo2018_SKS_T16s" 16.0 "${SRC_SKS}"  "${OBS_SKS}"  2101
run_hffk "hffk_Kuo2018_SKS_T20s" 20.0 "${SRC_SKS}"  "${OBS_SKS}"  2101
run_hffk "hffk_Kuo2018_SKS_T25s" 25.0 "${SRC_SKS}"  "${OBS_SKS}"  2101
run_hffk "hffk_Kuo2018_SKS_T33s" 33.0 "${SRC_SKS}"  "${OBS_SKS}"  2101
run_hffk "hffk_Kuo2018_SKS_T50s" 50.0 "${SRC_SKS}"  "${OBS_SKS}"  2101

# SKKS: 15 proxy events → 15 × 525 + 1 = 7876 lines
run_hffk "hffk_Kuo2018_SKKS_T4s"  4.0  "${SRC_SKKS}" "${OBS_SKKS}" 7876
run_hffk "hffk_Kuo2018_SKKS_T8s"  8.0  "${SRC_SKKS}" "${OBS_SKKS}" 7876
run_hffk "hffk_Kuo2018_SKKS_T16s" 16.0 "${SRC_SKKS}" "${OBS_SKKS}" 7876
run_hffk "hffk_Kuo2018_SKKS_T20s" 20.0 "${SRC_SKKS}" "${OBS_SKKS}" 7876
run_hffk "hffk_Kuo2018_SKKS_T25s" 25.0 "${SRC_SKKS}" "${OBS_SKKS}" 7876
run_hffk "hffk_Kuo2018_SKKS_T33s" 33.0 "${SRC_SKKS}" "${OBS_SKKS}" 7876
run_hffk "hffk_Kuo2018_SKKS_T50s" 50.0 "${SRC_SKKS}" "${OBS_SKKS}" 7876

echo ""; echo "===== Job D done $(date) ====="
