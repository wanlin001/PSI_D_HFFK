# VALIDATION STEP 2: HFFK forward (PSI_D_HFFK, SinkingSlab example)
# 目的：用 HFFK kernel 計算同一組觀測的 SI 預測，與 ray theory 輸出比較
# VanderBeek (2021): HFFK RMSE 148–165 ms vs ray 254–325 ms (P-wave)
# SKS 的比較結果預期方向相同但數值不同（SKS 更深，Fresnel zone 更窄）

ENV["JULIA_COPY_STACKS"] = "yes"
delete!(ENV, "TAUP_JAR")

push!(LOAD_PATH, "/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK/src")
using PSI_D

toml = joinpath(@__DIR__, "psi_parameters_hffk_forward.toml")
println("=== HFFK Forward (SinkingSlab) ===")
println("TOML: ", toml)
psi_forward(toml)
println("=== Done ===")
