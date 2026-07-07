#!/bin/bash
#PBS -N psi_hffk_val
#PBS -o hffk_val.out
#PBS -e hffk_val.err
#PBS -l nodes=1:ppn=4
#PBS -l walltime=02:00:00
#PBS -q q25g2

# ============================================================
# STEP B: HFFK Forward — SinkingSlab validation (PBS)
# wl2 (master cluster), qsub system
# ============================================================

set -e
echo "=== HFFK Forward: $(date) ==="

export PSI_D_ROOT="$HOME/software/ECOMAN2.0-seismology.PSI_D_HFFK"
export JULIA_BIN="$HOME/software/julia-1.10.0/bin/julia"
export JULIA_DEPOT_PATH="$PSI_D_ROOT/.julia_depot:$HOME/.julia"
export JULIA_PROJECT="$PSI_D_ROOT"
export JULIA_COPY_STACKS="yes"

WORK_DIR="$HOME/work/ASPECT/0624_HFFK_validation"
cd "${WORK_DIR}"

echo "Working dir: $(pwd)"
echo "PSI_D_HFFK: ${PSI_D_ROOT}"
echo "Julia: $($JULIA_BIN --version)"

$JULIA_BIN --compiled-modules=no -e "
    ENV[\"JULIA_COPY_STACKS\"] = \"yes\"
    delete!(ENV, \"TAUP_JAR\")
    using PSI_D
    psi_forward(\"validation/psi_parameters_hffk_forward.toml\")
"

echo "=== HFFK Done: $(date) ==="
ls -la "${WORK_DIR}/psi_output/HFFK_SinkingSlab/" 2>/dev/null
