module PSI_D
println("PSI_D: Plateform for Seismic Imaging - Deterministic")

using TauP
using StaticArrays
using SparseArrays
using IterativeSolvers
using Distributions
using FFTW
using LinearAlgebra
using TOML
using WriteVTK
using Dates
# Plots / SeismicDijkstra: optional, not needed for HFFK forward modeling
# Use try-catch to avoid hard dependency; these may not be in HFFK depot
try; @eval using Plots; catch; end
try; @eval using SeismicDijkstra; catch; end

include("psi_coordinate_systems.jl")
include("psi_fresnel_kernel.jl")   # HFFK: fresnel_radius, hffk_weight, fresnel_sample_offsets
include("psi_forward.jl")
include("psi_forward_elastic_tensor.jl")
include("psi_forward_splitting_parameters.jl")
include("psi_inverse.jl")
include("psi_output.jl")
include("psi_buildinputs.jl")
include("utilities.jl")

# wrapper_seismic_dijkstra.jl requires SeismicDijkstra (ShortestPath method)
# Not needed for HFFK forward modeling — skip to avoid missing dependency
# include("wrapper_seismic_dijkstra.jl")

export build_inputs, psi_forward, psi_inverse, psi_inverse!
export CompressionalWave, ShearWave, TravelTime, SplittingIntensity, SplittingParameters

end
