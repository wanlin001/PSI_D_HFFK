#!/bin/bash
#SBATCH -J psi_hffk_val
#SBATCH -p 8358
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -o hffk_val_%j.out
#SBATCH -e hffk_val_%j.err
#SBATCH --time=04:00:00

# ============================================================
# STEP B: HFFK Forward — SinkingSlab validation
# 目的: 用 PSI_D_HFFK 跑相同的 SinkingSlab 範例
#       與 STEP A 的 ray theory 輸出比較，驗證 HFFK kernel 正確性
# 預期: HFFK 與 ray theory 趨勢相同，但結構邊界附近更平滑
# ============================================================

set -e
echo "=== HFFK Forward: $(date) ==="

# --- 環境：指向 HFFK copy ---
export PSI_D_ROOT="/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK"
export PATH="/home/wl/software/julia-1.10.0/bin:$PATH"
export JULIA_DEPOT_PATH="$PSI_D_ROOT/.julia_depot:$HOME/.julia"   # 需要 ~/.julia 提供 Plots, SeismicDijkstra
export JULIA_PROJECT="$PSI_D_ROOT"
export JULIA_COPY_STACKS="yes"

WORK_DIR="/home/wl/work/ASPECT/0624_HFFK_validation"
cd "${WORK_DIR}"

echo "Working dir: $(pwd)"
echo "PSI_D version: ${PSI_D_ROOT}"
echo "Julia: $(julia --version)"

# --- Forward run with HFFK ---
julia -e "
    ENV[\"JULIA_COPY_STACKS\"] = \"yes\"
    delete!(ENV, \"TAUP_JAR\")
    using PSI_D
    psi_forward(\"validation/psi_parameters_hffk_forward.toml\")
"

echo "=== HFFK Done: $(date) ==="
echo "Output: ${WORK_DIR}/psi_output/HFFK_SinkingSlab/"
ls -la "${WORK_DIR}/psi_output/HFFK_SinkingSlab/" 2>/dev/null || echo "Output dir not yet created"
