#!/bin/bash
# Job B — Ray, Kuo2018 real events (SKS 4 evts + SKKS 15 evts proxy)
# Output : psi_output/ray_Kuo2018_SKS/
#          psi_output/ray_Kuo2018_SKKS/
# Submit : cd /home/wl/work/ASPECT/<PROJECT> && sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/scripts/slurm/jobB_ray_Kuo2018.sh
#SBATCH -J psi_B_ray_kuo
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH -o B_ray_Kuo2018_%j.out
#SBATCH -e B_ray_Kuo2018_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_v2.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_filled.dat"

export JULIA_COPY_STACKS=1

cd "${PROJECT_DIR}"

run_ray() {
    local TAG="$1" SRC="$2" OBS="$3" MIN_LINES="$4"
    local OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
    local OUTDAT="${OUTDIR}/SYN_SplittingParameters_ShearWave.dat"
    local TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG} already complete ($(wc -l < "${OUTDAT}") lines)"; return 0
    fi

    echo ""; echo "===== ${TAG} $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"
    rm -f "${OUTDIR}/psi_config.toml"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_DAT}|g" \
        -e "s|__OBS_TYPE__|SplittingParameters|g" \
        -e "s|__OBS_FILE__|${OBS}|g" \
        -e "s|__HFFK_SP__|tf_hffk_sp = false|g" \
        -e "s|__PERIOD__|8.0|g" \
        -e "s|__N_RINGS__|3|g" \
        -e "s|__N_AZIMUTH__|8|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"

    { printf "# dt_s, phi_rad, err_dt, err_phi, period_s, phase, src_id, rcv_id, quality, paz_rad\n"; \
      cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
    wc -l "${OUTDAT}"; echo "===== ${TAG} done ====="
}

# SKS: 4 real events, Δ≈91–108°
run_ray "ray_Kuo2018_SKS" \
    "${PSI_DIR}/psi_input/Sources_Kuo2018_SKS.dat" \
    "${PSI_DIR}/psi_input/DUMMY_SP_Kuo2018_SKS.dat" \
    2101

# SKKS: 15 proxy events, Δ≈150° same BAZ
run_ray "ray_Kuo2018_SKKS" \
    "${PSI_DIR}/psi_input/Sources_Kuo2018_SKKS.dat" \
    "${PSI_DIR}/psi_input/DUMMY_SP_Kuo2018_SKKS.dat" \
    7876

echo ""; echo "===== Job B done $(date) ====="
