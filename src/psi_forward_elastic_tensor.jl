
struct ElasticVoigt{T} <: SeismicVelocity
    ρ::Union{T, Nothing}
    c11::T
    c12::T
    c13::T
    c14::T
    c15::T
    c16::T
    c22::T
    c23::T
    c24::T
    c25::T
    c26::T
    c33::T
    c34::T
    c35::T
    c36::T
    c44::T
    c45::T
    c46::T
    c55::T
    c56::T
    c66::T
end
function ElasticVoigt(T::DataType, ndims; tf_density_normalized = false)
    ρ = tf_density_normalized ? nothing : zeros(T, ndims)

    return ElasticVoigt(ρ, zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims),
    zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims),
    zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims), zeros(T, ndims))
end
# Build ElasticVoigt (density-normalized) from HexagonalVectoralVelocity
function ElasticVoigt(Parameters::HexagonalVectoralVelocity)
    # Allocate parameters
    T = eltype(Parameters.f)
    ndims = size(Parameters.f)
    V = ElasticVoigt(T, ndims; tf_density_normalized = true)
    # Derive elastic tensor components
    for i in eachindex(Parameters.f)
        # Build tensor in principal coordinates
        α, β, ϵ, δ, γ = return_thomsen_parameters(Parameters, i)
        C = elastic_from_thomsen(1000.0*α, 1000.0*β, ϵ, δ, γ, Parameters.tf_exact)
        # Rotate the tensor
        R = rotation_matrix((0.5*π - Parameters.elevation[i], Parameters.azimuth[i]), (2,3))
        C = rotate_voigt(C, R)
        # Store result
        V.c11[i], V.c12[i], V.c13[i], V.c14[i], V.c15[i], V.c16[i],
        V.c22[i], V.c23[i], V.c24[i], V.c25[i], V.c26[i],
        V.c33[i], V.c34[i], V.c35[i], V.c36[i],
        V.c44[i], V.c45[i], V.c46[i],
        V.c55[i], V.c56[i], V.c66[i] = 
        (C[1,1], C[1,2], C[1,3], C[1,4], C[1,5], C[1,6],
        C[2,2], C[2,3], C[2,4], C[2,5], C[2,6],
        C[3,3], C[3,4], C[3,5], C[3,6],
        C[4,4], C[4,5], C[4,6],
        C[5,5], C[5,6], C[6,6])
    end

    return V
end

function read_model_parameters(io, parameterisation::Type{ElasticVoigt}, Mesh; dlm = ",", T = Float64,
    tf_global_cartesian = true, depth_reverse::String = "linear")
    # depth_reverse:
    #   "linear" — 官方 PSI_D；VIZTOMO 與 gen_benchmark_psitomo 相同寫檔順序（深度外層遞減）
    #   "dims3"  — 已廢棄（fork 早期誤用；非均勻場會錯置 Cij，勿用）
    # Note! There is currently a mix of psitomo models in local vs global coordinates
    tf_global_cartesian ? println("Reading Global Cartesian Tensors!") : println("Reading Local Cartesian Tensors!")
    # Read header information 
    tf_density_normalized = split(readline(io), dlm)
    tf_density_normalized = parse(Bool, strip(tf_density_normalized[1]))
    # Scalar to convert elastic coefficients from GPa to Pa
    # Density-normalized elastic coefficients are assumed to be in m²/s²
    a = tf_density_normalized ? 1.0 : 1.0e9
    # Allocate parameter structure
    Parameters = ElasticVoigt(T, size(Mesh); tf_density_normalized = tf_density_normalized)
    # Read and populate parameters
    k = 0
    x3_last = zeros(T,1)
    for line in readlines(io)
        k += 1
        line = split(line, dlm)
        x3_last[1] = parse(T, line[3])
        # Store elastic parameters (Pa or m²/s²)
        Parameters.c11[k] = a*parse(T, line[4])
        Parameters.c12[k] = a*parse(T, line[5])
        Parameters.c13[k] = a*parse(T, line[6])
        Parameters.c14[k] = a*parse(T, line[7])
        Parameters.c15[k] = a*parse(T, line[8])
        Parameters.c16[k] = a*parse(T, line[9])
        Parameters.c22[k] = a*parse(T, line[10])
        Parameters.c23[k] = a*parse(T, line[11])
        Parameters.c24[k] = a*parse(T, line[12])
        Parameters.c25[k] = a*parse(T, line[13])
        Parameters.c26[k] = a*parse(T, line[14])
        Parameters.c33[k] = a*parse(T, line[15])
        Parameters.c34[k] = a*parse(T, line[16])
        Parameters.c35[k] = a*parse(T, line[17])
        Parameters.c36[k] = a*parse(T, line[18])
        Parameters.c44[k] = a*parse(T, line[19])
        Parameters.c45[k] = a*parse(T, line[20])
        Parameters.c46[k] = a*parse(T, line[21])
        Parameters.c55[k] = a*parse(T, line[22])
        Parameters.c56[k] = a*parse(T, line[23])
        Parameters.c66[k] = a*parse(T, line[24])
        # Store density (kg/m³)
        tf_density_normalized ? nothing : Parameters.ρ[k] = parse(T, line[25])
        # Rotate to local geographic coordinates
        if tf_global_cartesian
            # Rotate the Voigt tensor from global cartesian to the local geographic coordinate system
            lon, lat = ( (π/180.0)*parse(T, line[1]),  (π/180.0)*parse(T, line[2]) )
            R = rotation_matrix((-lon, lat), (3, 2))
            C = return_voigt_matrix(Parameters, k)
            C = rotate_voigt(C, R)
            # Re-arrange the tensor such that c₁₁ = East = y_global, c₂₂ = North = z_global, and c₃₃ = Radial = x_global
            # after applying the above rotation
            Parameters.c11[k], Parameters.c12[k], Parameters.c13[k], Parameters.c14[k], Parameters.c15[k],
            Parameters.c16[k], Parameters.c22[k], Parameters.c23[k], Parameters.c24[k], Parameters.c25[k],
            Parameters.c26[k], Parameters.c33[k], Parameters.c34[k], Parameters.c35[k], Parameters.c36[k],
            Parameters.c44[k], Parameters.c45[k], Parameters.c46[k], Parameters.c55[k], Parameters.c56[k],
            Parameters.c66[k] = (C[2,2], C[2,3], C[1,2], C[2,5], C[2,6], C[2,4], C[3,3],
            C[1,3], C[3,5], C[3,6], C[3,4], C[1,1], C[1,5], C[1,6],
            C[1,4], C[5,5], C[5,6], C[4,5], C[6,6], C[4,6], C[4,4])
        end
    end
    close(io)

    # VIZTOMO elastic models are ordered by decreasing depth
    # PSI_D models are ordered by increasing depth
    should_reverse = if depth_reverse == "linear"
        true
    elseif depth_reverse == "dims3"
        x3_last[1] == Mesh.x[3][1]
    else
        error("Unknown depth_reverse = \"$depth_reverse\" (use \"linear\" or \"dims3\")")
    end
    if should_reverse
        println("Reversing depth (depth_reverse=$depth_reverse)")
        _reverse_elasticvoigt_depth!(Parameters, tf_density_normalized; mode = depth_reverse)
    end

    return Parameters
end


"""Reverse ElasticVoigt depth ordering after reading psitomo."""
function _reverse_elasticvoigt_depth!(Parameters::ElasticVoigt, tf_density_normalized::Bool; mode::String)
    if mode == "linear"
        reverse!(Parameters.c11)
        reverse!(Parameters.c12)
        reverse!(Parameters.c13)
        reverse!(Parameters.c14)
        reverse!(Parameters.c15)
        reverse!(Parameters.c16)
        reverse!(Parameters.c22)
        reverse!(Parameters.c23)
        reverse!(Parameters.c24)
        reverse!(Parameters.c25)
        reverse!(Parameters.c26)
        reverse!(Parameters.c33)
        reverse!(Parameters.c34)
        reverse!(Parameters.c35)
        reverse!(Parameters.c36)
        reverse!(Parameters.c44)
        reverse!(Parameters.c45)
        reverse!(Parameters.c46)
        reverse!(Parameters.c55)
        reverse!(Parameters.c56)
        reverse!(Parameters.c66)
        tf_density_normalized ? nothing : reverse!(Parameters.ρ)
    elseif mode == "dims3"
        reverse!(Parameters.c11, dims = 3)
        reverse!(Parameters.c12, dims = 3)
        reverse!(Parameters.c13, dims = 3)
        reverse!(Parameters.c14, dims = 3)
        reverse!(Parameters.c15, dims = 3)
        reverse!(Parameters.c16, dims = 3)
        reverse!(Parameters.c22, dims = 3)
        reverse!(Parameters.c23, dims = 3)
        reverse!(Parameters.c24, dims = 3)
        reverse!(Parameters.c25, dims = 3)
        reverse!(Parameters.c26, dims = 3)
        reverse!(Parameters.c33, dims = 3)
        reverse!(Parameters.c34, dims = 3)
        reverse!(Parameters.c35, dims = 3)
        reverse!(Parameters.c36, dims = 3)
        reverse!(Parameters.c44, dims = 3)
        reverse!(Parameters.c45, dims = 3)
        reverse!(Parameters.c46, dims = 3)
        reverse!(Parameters.c55, dims = 3)
        reverse!(Parameters.c56, dims = 3)
        reverse!(Parameters.c66, dims = 3)
        tf_density_normalized ? nothing : reverse!(Parameters.ρ, dims = 3)
    else
        error("Unknown depth reverse mode: $mode")
    end
    return Parameters
end


function return_kernel_parameter_type(::Type{ElasticVoigt{T}}, V) where {T <: Array}
    return ElasticVoigt{Vector{V}}
end

# RETURN KERNEL PARAMETERS: Elastic Voigt
function return_kernel_parameters(::SeismicObservable, Model::PsiModel{<:ElasticVoigt}, kernel_coordinates)
    # Define length and element type for kernel parameters
    n = length(kernel_coordinates)
    T = eltype(Model.Parameters.c11)
    # Initialise new parameter structure
    KernelParameters = ElasticVoigt(T, n; tf_density_normalized = isnothing(Model.Parameters.ρ))
    # Linearly interpolate model to the kernel
    interpolate_kernel_parameters!(KernelParameters, kernel_coordinates, Model)
    
    return KernelParameters
end
# INTERPOLATE KERNEL PARAMETERS: Elastic Voigt
function interpolate_kernel_parameters!(KernelParameters::ElasticVoigt, kernel_coordinates::Vector{NTuple{3, T}}, Model::PsiModel{<:ElasticVoigt}) where {T}
    # Loop to interpolate parameters to Kernel
    tf_density_normalized = isnothing(KernelParameters.ρ)
    for (i, qx_global) in enumerate(kernel_coordinates)
        # Convert the global kernel coordinates into the local coordinate system
        qx = global_to_local(qx_global[1], qx_global[2], qx_global[3], Model.Mesh.Geometry)
        # Get trilinear interpolation weights
        wind, wval = trilinear_weights(Model.Mesh.x, qx; tf_extrapolate = false, scale = 1.0)
        # Interpolate fields
        for (j, wj) in enumerate(wval)
            k = wind[j]
            tf_density_normalized ? nothing : KernelParameters.ρ[i] += wj*Model.Parameters.ρ[k]
            KernelParameters.c11[i] += wj*Model.Parameters.c11[k]
            KernelParameters.c12[i] += wj*Model.Parameters.c12[k]
            KernelParameters.c13[i] += wj*Model.Parameters.c13[k]
            KernelParameters.c14[i] += wj*Model.Parameters.c14[k]
            KernelParameters.c15[i] += wj*Model.Parameters.c15[k]
            KernelParameters.c16[i] += wj*Model.Parameters.c16[k]
            KernelParameters.c22[i] += wj*Model.Parameters.c22[k]
            KernelParameters.c23[i] += wj*Model.Parameters.c23[k]
            KernelParameters.c24[i] += wj*Model.Parameters.c24[k]
            KernelParameters.c25[i] += wj*Model.Parameters.c25[k]
            KernelParameters.c26[i] += wj*Model.Parameters.c26[k]
            KernelParameters.c33[i] += wj*Model.Parameters.c33[k]
            KernelParameters.c34[i] += wj*Model.Parameters.c34[k]
            KernelParameters.c35[i] += wj*Model.Parameters.c35[k]
            KernelParameters.c36[i] += wj*Model.Parameters.c36[k]
            KernelParameters.c44[i] += wj*Model.Parameters.c44[k]
            KernelParameters.c45[i] += wj*Model.Parameters.c45[k]
            KernelParameters.c46[i] += wj*Model.Parameters.c46[k]
            KernelParameters.c55[i] += wj*Model.Parameters.c55[k]
            KernelParameters.c56[i] += wj*Model.Parameters.c56[k]
            KernelParameters.c66[i] += wj*Model.Parameters.c66[k]

        end
    end

    return nothing
end

function kernel_phase_velocity(::CompressionalWave, Kernel::ObservableKernel{T1, T2}, index) where {T1, T2 <: ElasticVoigt}
    # Propagation vector
    k = spherical_to_cartesian(Kernel.weights[index].azimuth, Kernel.weights[index].elevation, 1.0)
    # Solve Christoffel equations. Because we do not need any polarization information, we can more efficiently
    # compute just the eigenvalues instead of using qp_phase_velocity.
    G = build_christoffel_matrix(k, Kernel.Parameters, index)
    # Compute phase velocities
    vq = eigvals(G)

    return 0.001*sqrt(maximum(vq)) # km/s
end
function kernel_phase_velocity(Phase::ShearWave, Kernel::ObservableKernel{T1, T2}, index) where {T1, T2 <: ElasticVoigt}
    # Compute phase velocities
    vs_slow, vs_fast, _, _, ps_fast, _ = phase_velocities(Kernel.weights[index].azimuth, Kernel.weights[index].elevation, Kernel.Parameters, index)
    # Angles
    fast_azimuth, fast_elevation, _ = cartesian_to_spherical(ps_fast[1], ps_fast[2], ps_fast[3])
    _, ζ = symmetry_axis_cosine(fast_azimuth, fast_elevation, Kernel.weights[index].azimuth, Kernel.weights[index].elevation, Phase.paz)
    # Effective anisotropic shear slowness
    uq = (1.0/vs_slow) - ((1.0/vs_slow) - (1.0/vs_fast))*(cos(ζ)^2)

    return 1.0/uq
end
function kernel_splitting_intensity(Phase::ShearWave, Kernel::ObservableKernel{T1, T2}, index) where {T1, T2 <: ElasticVoigt}
    # Compute phase velocities
    vs_slow, vs_fast, _, _, ps_fast, _ = phase_velocities(Kernel.weights[index].azimuth, Kernel.weights[index].elevation, Kernel.Parameters, index)
    # Angles
    fast_azimuth, fast_elevation, _ = cartesian_to_spherical(ps_fast[1], ps_fast[2], ps_fast[3])
    _, ζ = symmetry_axis_cosine(fast_azimuth, fast_elevation, Kernel.weights[index].azimuth, Kernel.weights[index].elevation, Phase.paz)
    # Splitting Intensity
    si = 0.5*((1.0/vs_slow) - (1.0/vs_fast))*sin(2.0*ζ)

    return si
end


function qp_phase_velocity(propagation_azimuth, propagation_elevation, Parameters::ElasticVoigt, index)
    # Return qP phase velocity and polarization computed via the Christoffel equations
    _, _, vqp, _, _, pqp = phase_velocities(propagation_azimuth, propagation_elevation, Parameters, index)
    qp_azimuth, qp_elevation, _ = cartesian_to_spherical(pqp[1], pqp[2], pqp[3])
    cosθ = symmetry_axis_cosine(qp_azimuth, qp_elevation, propagation_azimuth, propagation_elevation)
    # Note! The cosθ value returned is with respect to the qP-wave polarization and not the true symmetry axis.

    return vqp, cosθ
end
function qs_phase_velocities(propagation_azimuth, propagation_elevation, qt_polarization, Parameters::ElasticVoigt, index)
    # Compute qS-velocities
    vs_slow, vs_fast, _, _, ps_fast, _ = phase_velocities(propagation_azimuth, propagation_elevation, Parameters, index)
    # Polarization angles. For general anisotropic medium, there may not be a single symmetry axis and so the fast-polarisation
    # direction is taken as the reference.
    fast_azimuth, fast_elevation, _ = cartesian_to_spherical(ps_fast[1], ps_fast[2], ps_fast[3])
    cosθ, ζ = symmetry_axis_cosine(fast_azimuth, fast_elevation, propagation_azimuth, propagation_elevation, qt_polarization)
    # Note! The cosθ value returned is with respect to the fast-S-wave polarization and not the true symmetry axis.
    return vs_fast, vs_slow, cosθ, ζ
end
function phase_velocities(propagation_azimuth, propagation_elevation, V::ElasticVoigt, index)
    # Propagation vector
    k = spherical_to_cartesian(propagation_azimuth, propagation_elevation, 1.0)
    # Solve the Christoffel equations
    G = build_christoffel_matrix(k, V, index)
    vq, pq = eigen(G) # Returns static arrays because G::StaticArray
    vq = 0.001*sqrt.(vq) # km/s
    # Sort ascending (vs_slow, vs_fast, vp)
    n = sortperm(vq)
    vs_slow, vs_fast, vp = (vq[n[1]], vq[n[2]], vq[n[3]])
    ps_slow = (pq[1,n[1]], pq[2,n[1]], pq[3,n[1]])
    ps_fast = (pq[1,n[2]], pq[2,n[2]], pq[3,n[2]])
    pp = (pq[1,n[3]], pq[2,n[3]], pq[3,n[3]])

    return vs_slow, vs_fast, vp, ps_slow, ps_fast, pp
end



function build_christoffel_matrix(k::NTuple{3,T}, V::ElasticVoigt, index) where {T}
    # Unpack variables
    kx, ky, kz = k
    c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66 =
    (V.c11[index], V.c12[index], V.c13[index], V.c14[index], V.c15[index], V.c16[index],
    V.c22[index], V.c23[index], V.c24[index], V.c25[index], V.c26[index],
    V.c33[index], V.c34[index], V.c35[index], V.c36[index],
    V.c44[index], V.c45[index], V.c46[index],
    V.c55[index], V.c56[index], V.c66[index])
    # Elements of the 3x3 symmetric Christoffel matrix
    g11 = kx*(c11*kx + c16*ky + c15*kz) + ky*(c16*kx + c66*ky + c56*kz) + kz*(c15*kx + c56*ky + c55*kz)
    g12 = kx*(c16*kx + c66*ky + c56*kz) + ky*(c12*kx + c26*ky + c25*kz) + kz*(c14*kx + c46*ky + c45*kz)
    g13 = kx*(c15*kx + c56*ky + c55*kz) + ky*(c14*kx + c46*ky + c45*kz) + kz*(c13*kx + c36*ky + c35*kz)
    g22 = kx*(c66*kx + c26*ky + c46*kz) + ky*(c26*kx + c22*ky + c24*kz) + kz*(c46*kx + c24*ky + c44*kz)
    g23 = kx*(c56*kx + c25*ky + c45*kz) + ky*(c46*kx + c24*ky + c44*kz) + kz*(c36*kx + c23*ky + c34*kz)
    g33 = kx*(c55*kx + c45*ky + c35*kz) + ky*(c45*kx + c44*ky + c34*kz) + kz*(c35*kx + c34*ky + c33*kz)
    # Static arrays to avoid memory allocations
    G = @SMatrix [g11 g12 g13; g12 g22 g23; g13 g23 g33]
    # Scale inversely by density
    G = isnothing(V.ρ) ? G : (1.0/V.ρ[index])*G

    return G
end

function elastic_from_thomsen(α, β, ϵ, δ, γ, tf_exact; ρ = 1.0)
    # Diagonal elements
    c33 = ρ*α^2
    c44 = ρ*β^2
    c11 = (1.0 + 2.0*ϵ)*c33
    c66 = (1.0 + 2.0*γ)*c44
    # Off-diagonal elements
    c12 = c11 - 2.0*c66
    if tf_exact
        c13 = sqrt(δ*(c33^2) + 0.5*(c33 - c44)*(c11 + c33 - 2.0*c44)) - c44
    else
        c13 = sqrt(2.0*δ*c33*(c33 - c44) + ((c33 - c44)^2)) - c44
    end
    # Construct Voigt tensor
    C = @SMatrix [
        c11 c12 c13 0.0 0.0 0.0;
        c12 c11 c13 0.0 0.0 0.0;
        c13 c13 c33 0.0 0.0 0.0;
        0.0 0.0 0.0 c44 0.0 0.0;
        0.0 0.0 0.0 0.0 c44 0.0;
        0.0 0.0 0.0 0.0 0.0 c66
    ]

    return C
end
function return_voigt_matrix(V::ElasticVoigt, index)
    c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66 =
    (V.c11[index], V.c12[index], V.c13[index], V.c14[index], V.c15[index], V.c16[index],
    V.c22[index], V.c23[index], V.c24[index], V.c25[index], V.c26[index],
    V.c33[index], V.c34[index], V.c35[index], V.c36[index],
    V.c44[index], V.c45[index], V.c46[index],
    V.c55[index], V.c56[index], V.c66[index])
    C = @SMatrix [
        c11 c12 c13 c14 c15 c16;
        c12 c22 c23 c24 c25 c26;
        c13 c23 c33 c34 c35 c36;
        c14 c24 c34 c44 c45 c46;
        c15 c25 c35 c45 c55 c56;
        c16 c26 c36 c46 c56 c66
    ]

    return C
end

"""
    rotate_voigt(C, R)

    Rotates 6x6 Voigt tensor `C` given 3x3 rotation matrix `R` using
    the method of Bower (2009) 'Applied Mechanics of Solids'; Chapter 3
    that avoids trransformation to 4ᵗʰ-order tensor.
"""
function rotate_voigt(C, R)
    # Build 6x6 rotation matrix
    K = @SMatrix [
          R[1,1]^2      R[1,2]^2      R[1,3]^2       2.0*R[1,2]*R[1,3]               2.0*R[1,1]*R[1,3]               2.0*R[1,1]*R[1,2];
          R[2,1]^2      R[2,2]^2      R[2,3]^2       2.0*R[2,2]*R[2,3]               2.0*R[2,1]*R[2,3]               2.0*R[2,1]*R[2,2];
          R[3,1]^2      R[3,2]^2      R[3,3]^2       2.0*R[3,2]*R[3,3]               2.0*R[3,1]*R[3,3]               2.0*R[3,1]*R[3,2];
          R[2,1]*R[3,1] R[2,2]*R[3,2] R[2,3]*R[3,3] (R[2,2]*R[3,3] + R[2,3]*R[3,2]) (R[2,1]*R[3,3] + R[2,3]*R[3,1]) (R[2,1]*R[3,2] + R[2,2]*R[3,1]);
          R[1,1]*R[3,1] R[1,2]*R[3,2] R[1,3]*R[3,3] (R[1,2]*R[3,3] + R[1,3]*R[3,2]) (R[1,1]*R[3,3] + R[1,3]*R[3,1]) (R[1,1]*R[3,2] + R[1,2]*R[3,1]);
          R[1,1]*R[2,1] R[1,2]*R[2,2] R[1,3]*R[2,3] (R[1,2]*R[2,3] + R[1,3]*R[2,2]) (R[1,1]*R[2,3] + R[1,3]*R[2,1]) (R[1,1]*R[2,2] + R[1,2]*R[2,1])
          ]
    return K*C*transpose(K)
end

function interpolate_si_at_point(Phase::ShearWave,
                                  Kernel::ObservableKernel{<:SplittingIntensity, <:ElasticVoigt, <:Any, <:Any, <:Any},
                                  ray_idx::Int,
                                  qx_local,
                                  Model::PsiModel{<:ElasticVoigt},
                                  w)
    wind, wval = trilinear_weights(Model.Mesh.x, qx_local; tf_extrapolate = true, scale = 1.0)
    P = Model.Parameters

    c11=0.0; c12=0.0; c13=0.0; c14=0.0; c15=0.0; c16=0.0
    c22=0.0; c23=0.0; c24=0.0; c25=0.0; c26=0.0
    c33=0.0; c34=0.0; c35=0.0; c36=0.0
    c44=0.0; c45=0.0; c46=0.0
    c55=0.0; c56=0.0; c66=0.0
    ρ_interp = 0.0
    for (j, wj) in enumerate(wval)
        k = wind[j]
        c11 += wj*P.c11[k]; c12 += wj*P.c12[k]; c13 += wj*P.c13[k]
        c14 += wj*P.c14[k]; c15 += wj*P.c15[k]; c16 += wj*P.c16[k]
        c22 += wj*P.c22[k]; c23 += wj*P.c23[k]; c24 += wj*P.c24[k]
        c25 += wj*P.c25[k]; c26 += wj*P.c26[k]
        c33 += wj*P.c33[k]; c34 += wj*P.c34[k]; c35 += wj*P.c35[k]; c36 += wj*P.c36[k]
        c44 += wj*P.c44[k]; c45 += wj*P.c45[k]; c46 += wj*P.c46[k]
        c55 += wj*P.c55[k]; c56 += wj*P.c56[k]; c66 += wj*P.c66[k]
        if !isnothing(P.ρ)
            ρ_interp += wj * P.ρ[k]
        end
    end

    ρ_field = isnothing(P.ρ) ? nothing : [max(ρ_interp, 1.0)]
    ev = ElasticVoigt(ρ_field,
        [c11],[c12],[c13],[c14],[c15],[c16],
        [c22],[c23],[c24],[c25],[c26],
        [c33],[c34],[c35],[c36],
        [c44],[c45],[c46],
        [c55],[c56],[c66])

    vs_slow, vs_fast, _, _, ps_fast, _ = phase_velocities(w.azimuth, w.elevation, ev, 1)
    if vs_slow <= 0.0 || vs_fast <= 0.0
        return 0.0
    end
    fast_azimuth, fast_elevation, _ = cartesian_to_spherical(ps_fast[1], ps_fast[2], ps_fast[3])
    _, ζ = symmetry_axis_cosine(fast_azimuth, fast_elevation, w.azimuth, w.elevation, Phase.paz)
    return 0.5 * ((1.0/vs_slow) - (1.0/vs_fast)) * sin(2.0 * ζ)
end

# _hffk_smooth_cij! — defined here because ElasticVoigt is defined in this file
# Replaces Kernel.Parameters Cij at each ray point with the Fresnel-zone-weighted average,
# producing a banana-donut-smoothed elastic tensor for use by evaluate_kernel (waveform propagation).
function _hffk_smooth_cij!(Kernel::ObservableKernel{<:SplittingParameters, <:ElasticVoigt},
                            Observation::SplittingParameters,
                            Model::PsiModel{<:ElasticVoigt})
    period    = Observation.Forward.dominant_period
    n_rings   = Observation.Forward.n_rings
    n_azimuth = Observation.Forward.n_azimuth
    L         = sum(w.dr for w in Kernel.weights)
    x_cumul   = 0.0
    P         = Kernel.Parameters

    for (i, w) in enumerate(Kernel.weights)
        dr    = w.dr
        x_mid = x_cumul + dr / 2.0

        vs_ref = kernel_phase_velocity(Observation.Phase, Kernel, i)
        u_ref  = 1.0 / max(vs_ref, 0.1)
        R_f    = fresnel_radius(x_mid, L, period, u_ref)

        if R_f < 0.5
            x_cumul += dr
            continue
        end

        qx_global     = Kernel.coordinates[i]
        qx_local      = global_to_local(qx_global[1], qx_global[2], qx_global[3], Model.Mesh.Geometry)
        lon_ray       = qx_local[1]
        lat_ray       = qx_local[2]
        elv_ray       = qx_local[3]
        azimuth_ray   = w.azimuth
        elevation_ray = w.elevation

        offsets = fresnel_sample_offsets(R_f, n_rings, n_azimuth)

        c11=0.0;c12=0.0;c13=0.0;c14=0.0;c15=0.0;c16=0.0
        c22=0.0;c23=0.0;c24=0.0;c25=0.0;c26=0.0
        c33=0.0;c34=0.0;c35=0.0;c36=0.0
        c44=0.0;c45=0.0;c46=0.0
        c55=0.0;c56=0.0;c66=0.0
        ρv=0.0; total_w=0.0

        for (Δe, Δn, hffk_w) in offsets
            if Δe == 0.0 && Δn == 0.0
                continue
            end
            lon_off, lat_off = offset_to_geographic(lon_ray, lat_ray,
                                                     rad2deg(azimuth_ray),
                                                     rad2deg(elevation_ray),
                                                     Δe, Δn)
            # lon_off/lat_off 已在 local frame (degrees)，elv_ray = elevation (km)
            qx_off_local  = (lon_off, lat_off, elv_ray)
            wind, wval = trilinear_weights(Model.Mesh.x, qx_off_local;
                                           tf_extrapolate = true, scale = 1.0)
            MP = Model.Parameters
            for (j, wj) in enumerate(wval)
                k = wind[j]
                c11+=hffk_w*wj*MP.c11[k]; c12+=hffk_w*wj*MP.c12[k]; c13+=hffk_w*wj*MP.c13[k]
                c14+=hffk_w*wj*MP.c14[k]; c15+=hffk_w*wj*MP.c15[k]; c16+=hffk_w*wj*MP.c16[k]
                c22+=hffk_w*wj*MP.c22[k]; c23+=hffk_w*wj*MP.c23[k]; c24+=hffk_w*wj*MP.c24[k]
                c25+=hffk_w*wj*MP.c25[k]; c26+=hffk_w*wj*MP.c26[k]
                c33+=hffk_w*wj*MP.c33[k]; c34+=hffk_w*wj*MP.c34[k]; c35+=hffk_w*wj*MP.c35[k]
                c36+=hffk_w*wj*MP.c36[k]
                c44+=hffk_w*wj*MP.c44[k]; c45+=hffk_w*wj*MP.c45[k]; c46+=hffk_w*wj*MP.c46[k]
                c55+=hffk_w*wj*MP.c55[k]; c56+=hffk_w*wj*MP.c56[k]; c66+=hffk_w*wj*MP.c66[k]
                isnothing(MP.ρ) || (ρv += hffk_w*wj*MP.ρ[k])
            end
            total_w += hffk_w
        end

        if total_w > 0.0
            inv_w = 1.0 / total_w
            P.c11[i]=c11*inv_w; P.c12[i]=c12*inv_w; P.c13[i]=c13*inv_w
            P.c14[i]=c14*inv_w; P.c15[i]=c15*inv_w; P.c16[i]=c16*inv_w
            P.c22[i]=c22*inv_w; P.c23[i]=c23*inv_w; P.c24[i]=c24*inv_w
            P.c25[i]=c25*inv_w; P.c26[i]=c26*inv_w
            P.c33[i]=c33*inv_w; P.c34[i]=c34*inv_w; P.c35[i]=c35*inv_w; P.c36[i]=c36*inv_w
            P.c44[i]=c44*inv_w; P.c45[i]=c45*inv_w; P.c46[i]=c46*inv_w
            P.c55[i]=c55*inv_w; P.c56[i]=c56*inv_w; P.c66[i]=c66*inv_w
            isnothing(P.ρ) || (P.ρ[i] = ρv*inv_w)
        end

        x_cumul += dr
    end
    return nothing
end