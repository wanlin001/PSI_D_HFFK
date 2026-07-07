#!/bin/bash
#SBATCH -J psi_ray_val
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -o ray_val_%j.out
#SBATCH -e ray_val_%j.err
#SBATCH --time=02:00:00

# ============================================================
# STEP A: Ray Theory Forward — SinkingSlab validation
# 目的: 用原版 PSI_D 跑 SinkingSlab 範例的 SI forward prediction
#       作為 HFFK 的比較基準
# ============================================================

set -e
echo "=== Ray Theory Forward: $(date) ==="

# --- 環境 ---
export PSI_D_ROOT="/home/wl/software/ECOMAN2.0-seismology.PSI_D"
export PATH="/home/wl/software/julia-1.10.0/bin:$PATH"
export JULIA_DEPOT_PATH="$PSI_D_ROOT/.julia_depot"
export JULIA_PROJECT="$PSI_D_ROOT"
export JULIA_COPY_STACKS="yes"

WORK_DIR="/home/wl/work/ASPECT/0624_HFFK_validation"
cd "${WORK_DIR}"

echo "Working dir: $(pwd)"
echo "Julia: $(julia --version)"

# --- Forward run ---
julia -e "
    ENV[\"JULIA_COPY_STACKS\"] = \"yes\"
    delete!(ENV, \"TAUP_JAR\")
    using PSI_D
    psi_forward(\"psi_input/psi_parameters_synthetic.toml\")
"

echo "=== Ray Theory Done: $(date) ==="
echo "Output: ${WORK_DIR}/psi_output/SYN_SinkingBlock/"
ls -la "${WORK_DIR}/psi_output/SYN_SinkingBlock/" 2>/dev/null || echo "Output dir not yet created"
