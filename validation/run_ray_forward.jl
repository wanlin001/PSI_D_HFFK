# VALIDATION STEP 1: Ray theory forward (original PSI_D, SinkingSlab example)
# 目的：產生 ray theory SI 預測，作為 HFFK 比較的基準
# 對應: VanderBeek (2021) Table 2 ray theory RMSE = 254–325 ms (for P-wave)

ENV["JULIA_COPY_STACKS"] = "yes"
delete!(ENV, "TAUP_JAR")

# 用 HFFK 版本的 PSI_D（包含 wrapper_seismic_dijkstra.jl）
push!(LOAD_PATH, "/home/wl/software/ECOMAN2.0-seismology.PSI_D/src")
using PSI_D

toml = joinpath(@__DIR__, "../psi_input/psi_parameters_synthetic.toml")
println("=== Ray Theory Forward (SinkingSlab) ===")
println("TOML: ", toml)
psi_forward(toml)
println("=== Done ===")
