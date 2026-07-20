#!/bin/bash
# Job C-SP — HFFK SplittingParameters, uniform48 SYN events, T=4/8/16/20/25/33/50 s
# HFFK finite-frequency kernels applied to SP (tf_hffk_sp = true)
# Output : psi_output/hffksp_u48_T4s/
#          psi_output/hffksp_u48_T8s/
#          psi_output/hffksp_u48_T16s/
#          psi_output/hffksp_u48_T20s/
#          psi_output/hffksp_u48_T25s/
#          psi_output/hffksp_u48_T33s/
#          psi_output/hffksp_u48_T50s/
# Submit : cd /home/wl/work/ASPECT/<PROJECT> && sbatch /home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/scripts/slurm/jobC_hffksp_uniform48.sh
#SBATCH -J psi_Csp_hffk_u48
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH -o Csp_hffk_uniform48_%j.out
#SBATCH -e Csp_hffk_uniform48_%j.err

set -e

PSI_DIR=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
JULIA=/home/wl/software/julia-1.10.0/bin/julia
TEMPLATE="${PSI_DIR}/psi_input/psi_config_template_v2.toml"
MODEL_DAT="${PROJECT_DIR}/viztomo_output/psitomo0050_filled.dat"
SRC="${PSI_DIR}/psi_input/Sources_uniform48.dat"
OBS="${PSI_DIR}/psi_input/DUMMY_SP_uniform48.dat"
MIN_LINES=25201   # 48 × 525 + 1

export JULIA_COPY_STACKS=1

cd "${PROJECT_DIR}"

run_hffk() {
    local TAG="$1" PERIOD="$2"
    local OUTDIR="${PROJECT_DIR}/psi_output/${TAG}"
    local OUTDAT="${OUTDIR}/SYN_SplittingParameters_ShearWave.dat"
    local TMPDIR="${PROJECT_DIR}/.psi_tmp_${TAG}_$$"

    if [ -f "${OUTDAT}" ] && [ "$(wc -l < "${OUTDAT}")" -ge "${MIN_LINES}" ]; then
        echo "[SKIP] ${TAG} already complete ($(wc -l < "${OUTDAT}") lines)"; return 0
    fi
    rm -rf "${OUTDIR}"

    echo ""; echo "===== ${TAG} T=${PERIOD}s $(date) ====="
    mkdir -p "${OUTDIR}" "${TMPDIR}"
    rm -f "${OUTDIR}/psi_config.toml"

    sed -e "s|__OUTPUT_DIR__|${OUTDIR}|g" \
        -e "s|__MODEL_DAT__|${MODEL_DAT}|g" \
        -e "s|__OBS_TYPE__|SplittingParameters|g" \
        -e "s|__OBS_FILE__|${OBS}|g" \
        -e "s|__HFFK_SP__|tf_hffk_sp = true|g" \
        -e "s|__PERIOD__|${PERIOD}|g" \
        -e "s|__N_RINGS__|3|g" \
        -e "s|__N_AZIMUTH__|8|g" \
        -e "s|__SOURCES_DAT__|${SRC}|g" \
        -e "s|__TF_GLOBAL_CARTESIAN__|${TF_GLOBAL_CARTESIAN:-true}|g" \
        -e "s|__DEPTH_REVERSE__|${DEPTH_REVERSE:-linear}|g" \
        "${TEMPLATE}" > "${TMPDIR}/psi_config.toml"

    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
    rm -rf "${TMPDIR}"

    { printf "# dt_s, phi_rad, err_dt, err_phi, period_s, phase, src_id, rcv_id, quality, paz_rad
"; \
      cat "${OUTDAT}"; } > "${OUTDAT}.tmp" && mv "${OUTDAT}.tmp" "${OUTDAT}"
    wc -l "${OUTDAT}"; echo "===== ${TAG} done ====="
}

run_hffk "hffksp_u48_T4s"  4.0
run_hffk "hffksp_u48_T8s"  8.0
run_hffk "hffksp_u48_T16s" 16.0
run_hffk "hffksp_u48_T20s" 20.0
run_hffk "hffksp_u48_T25s" 25.0
run_hffk "hffksp_u48_T33s" 33.0
run_hffk "hffksp_u48_T50s" 50.0

echo ""; echo "===== Job C-SP (HFFK SplittingParameters) done $(date) ====="
