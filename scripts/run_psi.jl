#!/usr/bin/env julia
# Usage: julia --project=<psi_dir> scripts/run_psi.jl <project_dir>
# Example: julia --project=/home/wl/software/ECOMAN2.0-seismology.PSI_D_HFFK \
#            scripts/run_psi.jl /home/wl/work/ASPECT/0608_model4_V-50H100

ENV["JULIA_COPY_STACKS"] = "yes"
delete!(ENV, "TAUP_JAR")

using PSI_D

project_dir = length(ARGS) >= 1 ? ARGS[1] : error("Usage: run_psi.jl <project_dir>")
config = joinpath(project_dir, "psi_config.toml")
isfile(config) || error("Config not found: $config")

println("Running PSI_D: $config")
psi_forward(config)
println("Done: $project_dir")
