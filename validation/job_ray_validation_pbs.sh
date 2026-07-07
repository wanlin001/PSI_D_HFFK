#!/bin/bash
#PBS -N psi_ray_val
#PBS -o ray_val.out
#PBS -e ray_val.err
#PBS -l nodes=1:ppn=4
#PBS -l walltime=02:00:00
#PBS -q q25g2

# ============================================================
# STEP A: Ray Theory Forward — SinkingSlab validation (PBS)
# wl2 (master cluster), qsub system
# ============================================================

set -e
echo "=== Ray Theory Forward: $(date) ==="

export PSI_D_ROOT="$HOME/software/ECOMAN2.0-seismology.PSI_D"
export JULIA_BIN="$HOME/software/julia-1.10.0/bin/julia"
export JULIA_DEPOT_PATH="$PSI_D_ROOT/.julia_depot"
export JULIA_PROJECT="$PSI_D_ROOT"
export JULIA_COPY_STACKS="yes"

WORK_DIR="$HOME/work/ASPECT/0624_HFFK_validation"
cd "${WORK_DIR}"

echo "Working dir: $(pwd)"
echo "PSI_D: ${PSI_D_ROOT}"
echo "Julia: $($JULIA_BIN --version)"

$JULIA_BIN -e "
    ENV[\"JULIA_COPY_STACKS\"] = \"yes\"
    delete!(ENV, \"TAUP_JAR\")
    using PSI_D
    psi_forward(\"validation/psi_parameters_ray_forward.toml\")
"

echo "=== Ray Done: $(date) ==="
ls -la "${WORK_DIR}/psi_output/SYN_SinkingBlock/" 2>/dev/null
