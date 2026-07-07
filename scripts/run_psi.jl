#!/usr/bin/env julia
# Usage: julia --project=<repo_root> scripts/run_psi.jl <project_dir>
# Example: julia --project=. scripts/run_psi.jl projects/0608_model4_V-50H100

using PSI_D

project_dir = length(ARGS) >= 1 ? ARGS[1] : error("Usage: run_psi.jl <project_dir>")
config = joinpath(project_dir, "psi_config.toml")
isfile(config) || error("Config not found: $config")

println("Running PSI_D: $config")
PSI_D.run(config)
println("Done: $project_dir")
