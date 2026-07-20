# Shared defaults for project job scripts (source from software, not project folder).
#   source "${PSI_DIR}/scripts/psi_project_common.sh"
#
# psi_input (Sources / DUMMY / templates / Receivers / TauP) → always PSI_DIR
# MODEL_DAT / psi_output → project folder only

: "${PSI_DIR:=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK}"
: "${JULIA:=/home/wl/software/julia-1.10.0/bin/julia}"
: "${PROJECT_DIR:=${SLURM_SUBMIT_DIR:-$(pwd)}}"
: "${MODEL_STEP:=0050}"
: "${MODEL_DAT:=${PROJECT_DIR}/viztomo_output/psitomo${MODEL_STEP}_filled.dat}"
: "${TF_GLOBAL_CARTESIAN:=true}"
: "${DEPTH_REVERSE:=linear}"
: "${N_RINGS:=3}"
: "${N_AZIMUTH:=8}"
: "${PSI_INPUT:=${PSI_DIR}/psi_input}"
: "${RECEIVERS_DAT:=${PSI_INPUT}/Receivers.dat}"
: "${REFMODEL:=${PSI_INPUT}/ak135_no_crust.tvel}"

psi_sed_model_flags() {
    # Append to sed:  sed ... $(psi_sed_model_flags)
    echo "-e s|__TF_GLOBAL_CARTESIAN__|${TF_GLOBAL_CARTESIAN}|g -e s|__DEPTH_REVERSE__|${DEPTH_REVERSE}|g"
}

psi_run() {
    local TMPDIR="$1"
    ${JULIA} --project=${PSI_DIR} ${PSI_DIR}/scripts/run_psi.jl "${TMPDIR}"
}
