using LinearAlgebra

"""
    ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns r, z, and ψ interpolant
"""
function ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    z = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    PSI_interpolant = Interpolations.cubic_spline_interpolation((r, z), eqt2d.psi; extrapolation_bc=Interpolations.Line())
    return r, z, PSI_interpolant
end

"""
    ψ_interpolant(eqt2dv::IDSvector{IMAS.equilibrium__time_slice___profiles_2d})

Returns r, z, and ψ interpolant automatically choosing from available equilibrium profiles_2d grids
"""
function ψ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return ψ_interpolant(eqt2d)
end

"""
    ψ_interpolant(dd::IMAS.dd)

Returns r, z, and ψ interpolant automatically choosing from dd
"""
function ψ_interpolant(dd::IMAS.dd)
    return ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
end

"""
    Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::Array{T}, z::Array{T}) where {T<:Real}

Returns Br and Bz tuple evaluated at r and z starting from ψ interpolant
"""
function Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::Array{T}, z::Array{T}) where {T<:Real}
    # Check that r and z are the same size
    @assert size(r) == size(z)
    grad = [Interpolations.gradient(PSI_interpolant, r[idx], z[idx]) for idx in CartesianIndices(r)]
    Br = [grad[idx][2] / r[idx] / (2π) for idx in CartesianIndices(r)]
    Bz = [-grad[idx][1] / r[idx] / (2π) for idx in CartesianIndices(r)]
    return Br, Bz
end

function Br_Bz(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    Z, R = meshgrid(z, r)
    return Br_Bz(PSI_interpolant, R, Z)
end

function Br_Bz(eqt2dv::IDSvector{IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return Br_Bz(eqt2d)
end

"""
    Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{T}, z::Vector{T}) where {T<:Real}

Returns Bp evaluated at r and z starting from ψ interpolant
"""
function Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{T}, z::Vector{T}) where {T<:Real}
    Br, Bz = Br_Bz(PSI_interpolant, r, z)
    return sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)
end

function Bp(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    Z, R = meshgrid(z, r)
    return Bp(PSI_interpolant, R, Z)
end

function Bp(eqt2dv::IDSvector{IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return Bp(eqt2d)
end

"""
    find_psi_boundary(eqt::IMAS.equilibrium__time_slice; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)

Find psi value of the last closed flux surface
"""
function find_psi_boundary(eqt::IMAS.equilibrium__time_slice; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    if eqt2d !== nothing
        dimR = eqt2d.grid.dim1
        dimZ = eqt2d.grid.dim2
        PSI = eqt2d.psi
        psi = eqt.profiles_1d.psi
        R0 = eqt.global_quantities.magnetic_axis.r
        Z0 = eqt.global_quantities.magnetic_axis.z
        return find_psi_boundary(dimR, dimZ, PSI, psi, R0, Z0; precision, raise_error_on_not_open, raise_error_on_not_closed)
    else
        return eqt.profiles_1d.psi[end]
    end
end

"""
    find_psi_boundary(
        dimR::Union{AbstractVector{T},AbstractRange{T}},
        dimZ::Union{AbstractVector{T},AbstractRange{T}},
        PSI::Matrix{T},
        psi::Union{AbstractVector{T},AbstractRange{T}},
        R0::T,
        Z0::T;
        precision::Float64=1e-6,
        raise_error_on_not_open::Bool,
        raise_error_on_not_closed::Bool) where {T<:Real}

Find psi value of the last closed flux surface based on the cartesian grid flux surface tracing
"""
function find_psi_boundary(
    dimR::Union{AbstractVector{T},AbstractRange{T}},
    dimZ::Union{AbstractVector{T},AbstractRange{T}},
    PSI::Matrix{T},
    psi::Union{AbstractVector{T},AbstractRange{T}},
    R0::T,
    Z0::T;
    precision::Float64=1e-6,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool) where {T<:Real}

    psirange_init = [psi[1] * 0.9 + psi[end] * 0.1, psi[end] + 0.5 * (psi[end] - psi[1])]

    # innermost tentative flux surface (which should be closed!)
    surface, _ = flux_surface(dimR, dimZ, PSI, psi, R0, Z0, psirange_init[1], :closed)
    if isempty(surface)
        if raise_error_on_not_closed
            error("Flux surface at ψ=$(psirange_init[1]) is not closed; ψ=[$(psi[1])...$(psi[end])]")
        else
            return nothing
        end
    end

    # outermost tentative flux surface (which should be open!)
    surface, _ = flux_surface(dimR, dimZ, PSI, psi, R0, Z0, psirange_init[end], :closed)
    if !isempty(surface)
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open; ψ=[$(psi[1])...$(psi[end])]")
        else
            return nothing
        end
    end

    δd = sqrt((dimR[2] - dimR[1])^2 + (dimZ[2] - dimZ[1])^2)
    psirange = deepcopy(psirange_init)
    for k in 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        surface, _ = flux_surface(dimR, dimZ, PSI, psi, R0, Z0, psimid, :closed)
        # closed flux surface
        if !isempty(surface)
            ((pr, pz),) = surface
            psirange[1] = psimid
            if (abs(psirange[end] - psirange[1]) / abs(psirange[end] + psirange[1]) / 2.0) < precision
                if any(abs.([(minimum(pr) - minimum(dimR)), (maximum(pr) - maximum(dimR)), (minimum(pz) - minimum(dimZ)), (maximum(pz) - maximum(dimZ))]) .< δd)
                    return psi[end], psi[end]
                else
                    return psimid, psirange[end]
                end
            end
            # open flux surface
        else
            psirange[end] = psimid
        end
    end

    return error("Could not find closed boundary between ψ=$(psirange_init[1]) and ψ=$(psirange_init[end])")
end

"""
    find_psi_boundary(dd::IMAS.dd; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)

Find psi value of the last closed flux surface
"""
function find_psi_boundary(dd::IMAS.dd; precision::Float64=1e-6, raise_error_on_not_open::Bool=true, raise_error_on_not_closed::Bool=true)
    return find_psi_boundary(dd.equilibrium.time_slice[]; precision, raise_error_on_not_open, raise_error_on_not_closed)
end


"""
    find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation) 

Returns psi of the second magentic separatrix. This relies on the fact that find_x_points! saves the x points in such a way
that the last one is the null with Z opposite to the first x point which is the closest in psi to the lcfs.
"""
function find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation)
    psi2nd = PSI_interpolant.(eqt.boundary.x_point[end].r, eqt.boundary.x_point[end].z)
    return psi2nd * 0.999 .+ eqt.profiles_1d.psi[end] * 0.001
end

"""
    find_psi_2nd_separatrix(dd::IMAS.dd) 

Returns psi of the second magentic separatrix. This relies on the fact that find_x_points! saves the x points in such a way 
that the last one is the null with Z opposite to the first x point which is the closest in psi to the lcfs
"""

function find_psi_2nd_separatrix(dd::IMAS.dd)
    r, z, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
    return find_psi_2nd_separatrix(dd.equilibrium.time_slice[], PSI_interpolant)
end

"""
    find_psi_last_diverted(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::Vector{<:Real},
        wall_z::Vector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation;
        precision::Float64=1e-7)

Returns `psi_last_lfs, `psi_first_lfs_far`, and `null_within_wall`

psi_first_lfs_far will be the first surface inside OFL[:lfs_far]; psi_last_lfs will be the last surface inside OFL[:lfs]

Precision between the two is defined on the poloidal crossection area at the OMP (Psol*precision = power flowing between psi_first_lfs_far and psi_last_lfs ~ 0)
"""
function find_psi_last_diverted(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::Vector{<:Real},
    wall_z::Vector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation;
    precision::Float64=1e-7)

    # if no wall in dd, psi_last diverted not defined
    if isempty(wall_r) || isempty(wall_z)
        return (psi_last_lfs=NaN, psi_first_lfs_far=NaN, null_within_wall=true)
    end

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    Xpoint2 = [eqt.boundary.x_point[end].r, eqt.boundary.x_point[end].z]
    sign_z = sign(Xpoint2[2]) # sign of Z coordinate of 2nd null

    _, psi_separatrix = find_psi_boundary(eqt; raise_error_on_not_open=true) # psi LCFS
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt, PSI_interpolant) # psi second magnetic separatrix
    crossings = intersection([RA, maximum(wall_r)], [ZA, ZA], wall_r, wall_z)[2] # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
    psi_wall_midplane = PSI_interpolant.(r_wall_midplane, ZA)[1] # psi at the intersection between wall and omp

    # intersect 2nd separatrix with wall, and look 
    surface, _ = flux_surface(eqt, psi_2ndseparatrix, :open)
    r_intersect = Float64[]
    z_intersect = Float64[]
    r_max = 0.0
    for (r, z) in surface
        indexes, crossings = intersection(r, z, wall_r, wall_z) # find where flux surface crosses wall ("strike points" of surface)
        if isempty(crossings)
            continue
        end
        rr = (cr[1] for cr in crossings) # R coordiante of intersections btw 2nd separatrix and wall (could be more than 2)
        zz = (cr[2] for cr in crossings) # Z coordiante of intersections btw 2nd separatrix and wall (could be more than 2)

        # save all intersections with wall 
        append!(r_intersect, rr)
        append!(z_intersect, zz)
        r_max = max(r_max, maximum(r))
    end

    # r_mid(ψ) interpolator for region of interest
    r_mid_of_interest = 10.0 .^ range(log10(maximum(eqt.boundary.outline.r) * 0.99), log10(r_max), 1000)
    r_mid_itp = interp_rmid_at_psi(PSI_interpolant, r_mid_of_interest, ZA)

    if length(r_intersect) > 2
        # pick only the two intersections closest to 2nd Xpoint
        dist = (r_intersect .- Xpoint2[1]) .^ 2 + (z_intersect .- Xpoint2[2]) .^ 2
        r_intersect = r_intersect[partialsortperm(dist, 1:2)]
        z_intersect = z_intersect[partialsortperm(dist, 1:2)]
    end

    # order intersection clockwise (needed for null_within_wall)
    angle = 2 * π .- mod.(atan.(z_intersect .- ZA, r_intersect .- RA), 2 * π) # clockwise angle from midplane
    r_intersect = r_intersect[sortperm(angle)]
    z_intersect = z_intersect[sortperm(angle)]

    @assert length(r_intersect) == 2 # for safety, and to simplify eventual debugging

    # check if upper null is inside the wall, by checking if upper null is left/right of the vector between the 2 (ordered) intersections
    # This is an approximation (should work except for exotic walls)
    vec = [diff(r_intersect)[1], diff(z_intersect)[1]]
    vec2 = Xpoint2 - [r_intersect[1], z_intersect[1]]
    if vec[1] * vec2[2] - vec[2] * vec2[1] > 0
        # upper null on the left = is outside wall
        null_within_wall = false
    else
        # upper null on the right = is inside wall
        null_within_wall = true
    end

    # find the two surfaces `psi_first_lfs_far` and `psi_last_lfs` around the last diverted flux surface
    counter_max = 50
    counter = 0
    psi_axis_level = eqt.profiles_1d.psi[1] # psi value on axis 
    psi_sign = sign(psi_separatrix - psi_axis_level) # sign of the poloidal flux taking psi_axis = 0
    psi_first_lfs_far = psi_wall_midplane
    psi_last_lfs = psi_separatrix
    psi = (psi_first_lfs_far + psi_last_lfs) / 2
    err = Inf
    while abs(err) > precision && counter < counter_max
        surface, _ = flux_surface(eqt, psi, :open)
        for (r, z) in surface
            rr, zz, strike_angles, wall_index = line_wall_2_wall(r, z, wall_r, wall_z, RA, ZA)

            if isempty(rr) || all(zz .> ZA) || all(zz .< ZA)
                continue
            end

            if sign_z * zz[1] > 0 || sign_z * zz[end] > 0
                # psi intersects FW -> update upper bound
                psi_first_lfs_far = psi
            else
                # psi intersects divertor -> update lower bound
                psi_last_lfs = psi
            end
        end

        # better to compute error on poloidal area between [psi_last_lfs, psi_first_lfs_far] (needed for accurate power balance)
        r_up = r_mid_itp(psi_first_lfs_far)
        r_low = r_mid_itp(psi_last_lfs)

        A = π * (r_up^2 - r_low^2) # annular area between r_up and r_low [m^2] 
        err = abs(A)
        psi = (psi_first_lfs_far + psi_last_lfs) / 2

        counter = counter + 1
    end

    return (psi_last_lfs=psi_last_lfs, psi_first_lfs_far=psi_first_lfs_far, null_within_wall=null_within_wall)
end

"""
    find_psi_last_diverted(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, PSI_interpolant::Interpolations.AbstractInterpolation; precision::Float64=1e-7)

Returns `psi_last_lfs, `psi_first_lfs_far`, and `null_within_wall`

psi_first_lfs_far will be the first surface inside OFL[:lfs_far]; psi_last_lfs will be the last surface inside OFL[:lfs]

Precision between the two is defined on the poloidal crossection area at the OMP (Psol*precision = power flowing between psi_first_lfs_far and psi_last_lfs ~ 0)
"""

function find_psi_last_diverted(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, PSI_interpolant::Interpolations.AbstractInterpolation; precision::Float64=1e-7)
    fw = first_wall(wall)
    return find_psi_last_diverted(eqt, fw.r, fw.z, PSI_interpolant; precision)
end

"""
    find_psi_last_diverted(dd.IMAS.dd; precision::Float64=1e-7)

Returns `psi_last_lfs, `psi_first_lfs_far`, and `null_within_wall`

psi_first_lfs_far will be the first surface inside OFL[:lfs_far]; psi_last_lfs will be the last surface inside OFL[:lfs]

Precision between the two is defined on the poloidal crossection area at the OMP (Psol*precision = power flowing between psi_first_lfs_far and psi_last_lfs ~ 0)
"""

function find_psi_last_diverted(dd::IMAS.dd; precision::Float64=1e-7)
    rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
    return find_psi_last_diverted(dd.equilibrium.time_slice[], dd.wall, PSI_interpolant; precision)
end


"""
    interp_rmid_at_psi(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation, R::AbstractVector{<:Real})

Returns the interpolant r_mid(ψ) to compute the r at the midplane of the flux surface identified by ψ

The vector `R` defines the sampling of interest for thie interpolation
"""
function interp_rmid_at_psi(PSI_interpolant::Interpolations.AbstractInterpolation, R::AbstractVector{T}, ZA::T) where {T<:Real}
    return interp1d(PSI_interpolant.(R, R .* 0.0 .+ ZA), R, :cubic)
end

"""
    flux_surfaces(eq::equilibrium; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities in the equilibrium IDS
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eq::equilibrium; upsample_factor::Int=1)
    for time_index in eachindex(eq.time_slice)
        flux_surfaces(eq.time_slice[time_index]; upsample_factor)
    end
    return eq
end

"""
    flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice; upsample_factor::Int=1)
    R0, B0 = vacuum_r0_b0(eqt)
    return flux_surfaces(eqt, B0, R0; upsample_factor)
end

function find_magnetic_axis!(r::AbstractVector{<:Real}, z::AbstractVector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation, psi_sign::Real)
    res = Optim.optimize(
        x -> begin
            try
                PSI_interpolant(x[1], x[2]) * psi_sign
            catch e
                if typeof(e) <: BoundsError
                    return Inf
                else
                    rethrow(e)
                end
            end
        end,
        [r[Int(round(length(r) / 2))], z[Int(round(length(z) / 2))]],
        Optim.Newton(),
        Optim.Options(; g_tol=1E-8);
        autodiff=:forward
    )
    return res.minimizer[1], res.minimizer[2]
end

"""
    flux_surfaces(eqt::equilibrium__time_slice{T}, B0::T, R0::T; upsample_factor::Int=1) where {T<:Real}

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice, B0 and R0
The original psi grid can be upsampled by a `upsample_factor` to get higher resolution flux surfaces
"""
function flux_surfaces(eqt::equilibrium__time_slice{T}, B0::T, R0::T; upsample_factor::Int=1) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    PSI = eqt2d.psi

    # upsampling for high-resolution r,z flux surface coordinates
    if upsample_factor > 1
        r = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1) * upsample_factor)
        z = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2) * upsample_factor)
        PSI = PSI_interpolant(r, z)
    end

    psi_sign = sign(eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1])

    # find magnetic axis
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = find_magnetic_axis!(r, z, PSI_interpolant, psi_sign)
    psi_axis = PSI_interpolant(eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z)
    eqt.profiles_1d.psi =
        (eqt.profiles_1d.psi .- eqt.profiles_1d.psi[1]) ./ (eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1]) .* (eqt.profiles_1d.psi[end] - psi_axis) .+ psi_axis

    for item in (
        :b_field_average,
        :b_field_max,
        :b_field_min,
        :elongation,
        :triangularity_lower,
        :triangularity_upper,
        :squareness_lower_inner,
        :squareness_lower_outer,
        :squareness_upper_inner,
        :squareness_upper_outer,
        :r_inboard,
        :r_outboard,
        :q,
        :surface,
        :dvolume_dpsi,
        :j_tor,
        :gm1,
        :gm2,
        :gm4,
        :gm5,
        :gm8,
        :gm9,
        :fsa_bp,
        :trapped_fraction
    )
        setproperty!(eqt.profiles_1d, item, zeros(eltype(eqt.profiles_1d.psi), size(eqt.profiles_1d.psi)))
    end

    PR = Vector{T}[]
    PZ = Vector{T}[]
    LL = Vector{T}[]
    FLUXEXPANSION = Vector{T}[]
    INT_FLUXEXPANSION_DL = zeros(T, length(eqt.profiles_1d.psi))
    BPL = zeros(T, length(eqt.profiles_1d.psi))
    for (k, psi_level0) in reverse!(collect(enumerate(eqt.profiles_1d.psi)))

        if k == 1 # on axis flux surface is a synthetic one
            eqt.profiles_1d.elongation[1] = eqt.profiles_1d.elongation[2] - (eqt.profiles_1d.elongation[3] - eqt.profiles_1d.elongation[2])
            eqt.profiles_1d.triangularity_upper[1] = 0.0
            eqt.profiles_1d.triangularity_lower[1] = 0.0

            a = (eqt.profiles_1d.r_outboard[2] - eqt.profiles_1d.r_inboard[2]) / 100.0
            b = eqt.profiles_1d.elongation[1] * a

            t = range(0, 2π, 17)
            pr = cos.(t) .* a .+ eqt.global_quantities.magnetic_axis.r
            pz = sin.(t) .* b .+ eqt.global_quantities.magnetic_axis.z

            # Extrema on array indices
            (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)

        else  # other flux surfaces
            # trace flux surface
            ((pr, pz),), psi_level =
                flux_surface(r, z, PSI, eqt.profiles_1d.psi, eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, psi_level0, :closed)
            if isempty(pr)
                # p = heatmap(r, z, PSI'; colorbar=true, aspect_ratio=:equal)
                # contour!(r, z, PSI'; color=:white, levels=100)
                # contour!(r, z, PSI'; levels=[eqt.profiles_1d.psi[end]], color=:white, lw=2)
                # display(p)
                error("IMAS: Could not trace closed flux surface $k out of $(length(eqt.profiles_1d.psi)) at ψ = $(psi_level)")
            end

            # Extrema on array indices
            (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)

            # accurate geometric quantities by finding geometric extrema as optimization problem
            w = 1E-4 # push away from magnetic axis
            function fx(x::AbstractVector{<:Real}, psi_level::Float64, eqt::IMAS.equilibrium__time_slice, w::Float64)
                try
                    (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[1] - eqt.global_quantities.magnetic_axis.r)^2 * w
                catch
                    return 100
                end
            end
            function fz(x::AbstractVector{<:Real}, psi_level::Float64, eqt::IMAS.equilibrium__time_slice, w::Float64)
                try
                    (PSI_interpolant(x[1], x[2]) - psi_level)^2 - (x[2] - eqt.global_quantities.magnetic_axis.z)^2 * w
                catch
                    return 100
                end
            end
            res = Optim.optimize(x -> fx(x, psi_level, eqt, w), [max_r, z_at_max_r], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
            (max_r, z_at_max_r) = (res.minimizer[1], res.minimizer[2])
            res = Optim.optimize(x -> fx(x, psi_level, eqt, w), [min_r, z_at_min_r], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
            (min_r, z_at_min_r) = (res.minimizer[1], res.minimizer[2])
            if psi_level0 != eqt.profiles_1d.psi[end]
                res = Optim.optimize(x -> fz(x, psi_level, eqt, w), [r_at_max_z, max_z], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
                (r_at_max_z, max_z) = (res.minimizer[1], res.minimizer[2])
                res = Optim.optimize(x -> fz(x, psi_level, eqt, w), [r_at_min_z, min_z], Optim.Newton(), Optim.Options(; g_tol=1E-8); autodiff=:forward)
                (r_at_min_z, min_z) = (res.minimizer[1], res.minimizer[2])
            end
            # p = plot(pr, pz, label = "")
            # plot!([max_r], [z_at_max_r], marker = :cicle)
            # plot!([min_r], [z_at_min_r], marker = :cicle)
            # plot!([r_at_max_z], [max_z], marker = :cicle)
            # plot!([r_at_min_z], [min_z], marker = :cicle)
            # display(p)

            # plasma boundary information
            if k == length(eqt.profiles_1d.psi)
                eqt.boundary.outline.r = pr
                eqt.boundary.outline.z = pz
            end
        end

        eqt.profiles_1d.r_outboard[k] = max_r
        eqt.profiles_1d.r_inboard[k] = min_r

        # miller geometric coefficients
        _, _, κ, δu, δl, ζou, ζol, ζil, ζiu = miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
        eqt.profiles_1d.elongation[k] = κ
        eqt.profiles_1d.triangularity_upper[k] = δu
        eqt.profiles_1d.triangularity_lower[k] = δl
        eqt.profiles_1d.squareness_lower_outer[k] = ζol
        eqt.profiles_1d.squareness_upper_outer[k] = ζou
        eqt.profiles_1d.squareness_lower_inner[k] = ζil
        eqt.profiles_1d.squareness_upper_inner[k] = ζiu

        # poloidal magnetic field (with sign)
        Br, Bz = Br_Bz(PSI_interpolant, pr, pz)
        Bp2 = Br .^ 2.0 .+ Bz .^ 2.0
        Bp_abs = sqrt.(Bp2)
        Bp = (
            Bp_abs .*
            sign.((pz .- eqt.global_quantities.magnetic_axis.z) .* Br .- (pr .- eqt.global_quantities.magnetic_axis.r) .* Bz)
        )

        # flux expansion
        dl = vcat(0.0, sqrt.(diff(pr) .^ 2 + diff(pz) .^ 2))
        ll = cumsum(dl)
        fluxexpansion = 1.0 ./ Bp_abs
        int_fluxexpansion_dl = integrate(ll, fluxexpansion)
        Bpl = integrate(ll, Bp)

        # save flux surface coordinates for later use
        pushfirst!(PR, pr)
        pushfirst!(PZ, pz)
        pushfirst!(LL, ll)
        pushfirst!(FLUXEXPANSION, fluxexpansion)
        INT_FLUXEXPANSION_DL[k] = int_fluxexpansion_dl
        BPL[k] = Bpl

        # trapped fraction
        Bt = eqt.profiles_1d.f[k] ./ pr
        Btot = sqrt.(Bp2 .+ Bt .^ 2)
        Bmin = minimum(Btot)
        Bmax = maximum(Btot)
        Bratio = Btot ./ Bmax
        avg_Btot = flxAvg(Btot, ll, fluxexpansion, int_fluxexpansion_dl)
        avg_Btot2 = flxAvg(Btot .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)
        hf = flxAvg((1.0 .- sqrt.(1.0 .- Bratio) .* (1.0 .+ Bratio ./ 2.0)) ./ Bratio .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)
        h = avg_Btot / Bmax
        h2 = avg_Btot2 / Bmax^2
        ftu = 1.0 - h2 / (h^2) * (1.0 - sqrt(1.0 - h) * (1.0 + 0.5 * h))
        ftl = 1.0 - h2 * hf
        eqt.profiles_1d.trapped_fraction[k] = 0.75 * ftu + 0.25 * ftl

        # Bavg
        eqt.profiles_1d.b_field_average[k] = avg_Btot

        # Bmax
        eqt.profiles_1d.b_field_max[k] = Bmax

        # Bmin
        eqt.profiles_1d.b_field_min[k] = Bmin

        # gm1 = <1/R^2>
        eqt.profiles_1d.gm1[k] = flxAvg(1.0 ./ pr .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)

        # gm4 = <1/B^2>
        eqt.profiles_1d.gm4[k] = flxAvg(1.0 ./ Btot .^ 2, ll, fluxexpansion, int_fluxexpansion_dl)

        # gm5 = <B^2>
        eqt.profiles_1d.gm5[k] = avg_Btot2

        # gm8 = <R>
        eqt.profiles_1d.gm8[k] = flxAvg(pr, ll, fluxexpansion, int_fluxexpansion_dl)

        # gm9 = <1/R>
        eqt.profiles_1d.gm9[k] = flxAvg(1.0 ./ pr, ll, fluxexpansion, int_fluxexpansion_dl)

        # fsa_bp = <Bp>
        eqt.profiles_1d.fsa_bp[k] = flxAvg(Bp, ll, fluxexpansion, int_fluxexpansion_dl)

        # j_tor = <j_tor/R> / <1/R> [A/m²]
        eqt.profiles_1d.j_tor[k] =
            (
                -(eqt.profiles_1d.dpressure_dpsi[k] + eqt.profiles_1d.f_df_dpsi[k] * eqt.profiles_1d.gm1[k] / constants.μ_0) *
                (2π)
            ) / eqt.profiles_1d.gm9[k]

        # dvolume_dpsi
        eqt.profiles_1d.dvolume_dpsi[k] = sign(eqt.profiles_1d.fsa_bp[k]) * int_fluxexpansion_dl

        # surface area
        eqt.profiles_1d.surface[k] = 2π * sum(pr .* dl)

        # q
        eqt.profiles_1d.q[k] = eqt.profiles_1d.dvolume_dpsi[k] .* eqt.profiles_1d.f[k] .* eqt.profiles_1d.gm1[k] ./ (2π)

        # quantities calculated on the last closed flux surface
        if k == length(eqt.profiles_1d.psi)
            # perimeter
            eqt.global_quantities.length_pol = ll[end]
        end
    end

    # area
    eqt.profiles_1d.area = cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* eqt.profiles_1d.gm9) ./ 2π

    # volume
    eqt.profiles_1d.volume = cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi)

    # phi
    eqt.profiles_1d.phi = cumul_integrate(eqt.profiles_1d.volume, eqt.profiles_1d.f .* eqt.profiles_1d.gm1) / (2π)

    # rho_tor_norm
    rho = sqrt.(abs.(eqt.profiles_1d.phi ./ (π * B0)))
    rho_meters = rho[end]
    eqt.profiles_1d.rho_tor = rho
    eqt.profiles_1d.rho_tor_norm = rho ./ rho_meters

    # phi 2D
    eqt2d.phi =
        Interpolations.cubic_spline_interpolation(
            to_range(eqt.profiles_1d.psi) * psi_sign,
            eqt.profiles_1d.phi;
            extrapolation_bc=Interpolations.Line()
        ).(eqt2d.psi * psi_sign)

    # rho 2D in meters
    RHO = sqrt.(abs.(eqt2d.phi ./ (π * B0)))

    # gm2: <∇ρ²/R²>
    if false
        RHO_interpolant = Interpolations.cubic_spline_interpolation((r, z), RHO)
        for k in 1:length(eqt.profiles_1d.psi)
            tmp = [Interpolations.gradient(RHO_interpolant, PR[k][j], PZ[k][j]) for j in 1:length(PR[k])]
            dPHI2 = [j[1] .^ 2.0 .+ j[2] .^ 2.0 for j in tmp]
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k] .^ 2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    else
        dRHOdR, dRHOdZ = gradient(collect(r), collect(z), RHO)
        dPHI2_interpolant = Interpolations.cubic_spline_interpolation((r, z), dRHOdR .^ 2.0 .+ dRHOdZ .^ 2.0)
        for k in 1:length(eqt.profiles_1d.psi)
            dPHI2 = dPHI2_interpolant.(PR[k], PZ[k])
            eqt.profiles_1d.gm2[k] = flxAvg(dPHI2 ./ PR[k] .^ 2.0, LL[k], FLUXEXPANSION[k], INT_FLUXEXPANSION_DL[k])
        end
    end
    eqt.profiles_1d.gm2[1] =
        Interpolations.cubic_spline_interpolation(
            to_range(eqt.profiles_1d.psi[2:end]) * psi_sign,
            eqt.profiles_1d.gm2[2:end];
            extrapolation_bc=Interpolations.Line()
        ).(eqt.profiles_1d.psi[1] * psi_sign)

    # ip
    eqt.global_quantities.ip = IMAS.integrate(eqt.profiles_1d.area, eqt.profiles_1d.j_tor)
    # eqt.global_quantities.ip = IMAS.integrate(eqt.profiles_1d.volume, eqt.profiles_1d.j_tor.*eqt.profiles_1d.gm9) / (2π) # equivalent

    # Geometric major and minor radii
    Rgeo = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
    a = (eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]) / 2.0

    # vacuum magnetic field at the geometric center
    Btvac = B0 * R0 / Rgeo

    # average poloidal magnetic field
    Bpave = eqt.global_quantities.ip * constants.μ_0 / eqt.global_quantities.length_pol

    # li
    Bp2v = integrate(eqt.profiles_1d.psi, BPL)
    eqt.global_quantities.li_3 = 2.0 * Bp2v / Rgeo / (eqt.global_quantities.ip * constants.μ_0)^2

    # beta_tor
    avg_press = volume_integrate(eqt, eqt.profiles_1d.pressure) / eqt.profiles_1d.volume[end]
    eqt.global_quantities.beta_tor = abs(avg_press / (Btvac^2 / 2.0 / constants.μ_0))

    # beta_pol
    eqt.global_quantities.beta_pol = abs(avg_press / (Bpave^2 / 2.0 / constants.μ_0))

    # beta_normal
    ip = eqt.global_quantities.ip / 1e6
    eqt.global_quantities.beta_normal = eqt.global_quantities.beta_tor / abs(ip / a / Btvac) * 100

    # find quantities on separatrix
    find_x_point!(eqt)
    find_strike_points!(eqt)

    # secondary separatrix
    if length(eqt.boundary.x_point) > 1
        psi2nd = find_psi_2nd_separatrix(eqt, PSI_interpolant)
        tmp, _ = flux_surface(r, z, PSI, eqt.profiles_1d.psi, eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, psi2nd, :encircling)
        if !isempty(tmp)
            (pr2nd, pz2nd) = tmp[1]
            eqt.boundary_secondary_separatrix.outline.r = pr2nd
            eqt.boundary_secondary_separatrix.outline.z = pz2nd
        end
    end

    return eqt
end

"""
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real)

Returns r,z coordiates of closed flux surface at given psi_level
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real)
    return flux_surface(eqt, psi_level, :closed)
end

"""
    flux_surface(eqt::equilibrium__time_slice, psi_level::Real, type::Symbol)

Returns tuple of (r,z) coordiates of flux surface at given psi_level, and `psi_level` actually used

  - psi[1] returns psi[2]
  - psi[end] triggers accurate finding of lcfs value

The `type` parameter:

  - :any, return all contours
  - :closed, all closed flux-surface that encircle the magnetic axis
  - :open, all open flux-surfaces
  - :encircling, open flux-surfaces encircling the magnetic axis
"""
function flux_surface(eqt::equilibrium__time_slice, psi_level::Real, type::Symbol)
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    dim1 = eqt2d.grid.dim1
    dim2 = eqt2d.grid.dim2
    PSI = eqt2d.psi
    psi = eqt.profiles_1d.psi
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z
    return flux_surface(dim1, dim2, PSI, psi, R0, Z0, psi_level, type)
end

function flux_surface(
    dim1::Union{AbstractVector{T},AbstractRange{T}},
    dim2::Union{AbstractVector{T},AbstractRange{T}},
    PSI::AbstractArray{T},
    psi::Union{AbstractVector{T},AbstractRange{T}},
    R0::T,
    Z0::T,
    psi_level::T,
    type::Symbol) where {T<:Real}

    if psi_level == psi[1]
        # handle on axis value as the first flux surface
        psi_level = psi[2]

    elseif psi_level == psi[end]
        # handle boundary by finding accurate lcfs psi
        psi__boundary_level, _ = find_psi_boundary(dim1, dim2, PSI, psi, R0, Z0; raise_error_on_not_open=false, raise_error_on_not_closed=false)
        if psi__boundary_level !== nothing
            if abs(psi__boundary_level - psi_level) < abs(psi[end] - psi[end-1])
                psi_level = psi__boundary_level
            end
        end
    end

    # contouring routine
    cl = Contour.contour(dim1, dim2, PSI, psi_level)

    prpz = Tuple{Vector{T},Vector{T}}[]
    if type == :any
        # if no open/closed check, then return all contours
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
            push!(prpz, (pr, pz))
        end

    elseif type == :closed
        # look for closed flux-surface
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surface that closes and contains magnetic axis
            if (pr[1] == pr[end]) && (pz[1] == pz[end]) && (PolygonOps.inpolygon((R0, Z0), collect(zip(pr, pz))) == 1)
                reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
                push!(prpz, (pr, pz))
                break
            end
        end

    elseif type == :open
        # look for open flux-surfaces
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surfaces that do not close
            if (pr[1] != pr[end]) || (pz[1] != pz[end])
                reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
                push!(prpz, (pr, pz))
            end
        end

    elseif type == :encircling
        # look for open flux-surfaces that encircle the magnetic axis
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surfaces that do not close
            if (pr[1] != pr[end]) || (pz[1] != pz[end])
                tmp = collect(zip(pr, pz))
                push!(tmp, tmp[1]) # close it
                if PolygonOps.inpolygon((R0, Z0), tmp) == 1
                    reorder_flux_surface!(pr, pz, R0, Z0; force_close=false)
                    push!(prpz, (pr, pz))
                end
            end
        end
    else
        error("flux_surface type `$type` is not recognized: can be [:any, :closed, :open, :encircling]")
    end

    return Tuple(prpz), psi_level
end

"""
    tweak_psi_to_match_psilcfs!(eqt::IMAS.equilibrium.time_slice{D}; ψbound::Union{Nothing,D}=nothing) where {D<:Real}

Tweak `eqt.profiles_2d[:].psi` and `eqt.profiles_1d.psi` so that contouring finds LCFS at `eqt.profiles_1d.psi[end]`.

If `ψbound !== nothing` this routine also shifts `eqt.profiles_2d[:].psi` and `eqt.profiles_1d.psi` so that `eqt.profiles_1d.psi[end] == ψbound`.
"""
function tweak_psi_to_match_psilcfs!(eqt::IMAS.equilibrium__time_slice{D}; ψbound::Union{Nothing,D}=nothing) where {D<:Real}
    eqt1d = eqt.profiles_1d
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)

    psia = eqt1d.psi[1]
    psib = eqt1d.psi[end]
    if ψbound === nothing
        delta_psib = 0.0
    else
        delta_psib = ψbound - psib
    end

    # retrace the last closed flux surface
    true_psib, _ = IMAS.find_psi_boundary(eqt)
    if true_psib !== nothing
        # scale psirz so to match original psi bounds (also add delta_psib to get desired ψbound)
        @. eqt2d.psi = (eqt2d.psi - psia) * (psib - psia) / (true_psib - psia) + psia + delta_psib
        @. eqt1d.psi = eqt1d.psi + delta_psib
    end

    return nothing
end

function flxAvg(input::AbstractVector{T}, ll::AbstractVector{T}, fluxexpansion::AbstractVector{T}, int_fluxexpansion_dl::T)::T where {T<:Real}
    return integrate(ll, input .* fluxexpansion) / int_fluxexpansion_dl
end

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over volume
"""
function volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    return integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what)
end

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over volume
"""
function cumlul_volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    return cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what)
end

"""
    surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over surface
"""
function surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    return integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what .* eqt.profiles_1d.gm9) ./ 2π
end

"""
    surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over surface
"""
function cumlul_surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    return cumul_integrate(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi .* what .* eqt.profiles_1d.gm9) ./ 2π
end

"""
    find_x_point!(eqt::IMAS.equilibrium__time_slice)::IDSvector{<:IMAS.equilibrium__time_slice___boundary__x_point}

Find the `n` X-points that are closest to the separatrix
"""
function find_x_point!(eqt::IMAS.equilibrium__time_slice)::IDSvector{<:IMAS.equilibrium__time_slice___boundary__x_point}
    ((rlcfs, zlcfs),), _ = flux_surface(eqt, eqt.profiles_1d.psi[end], :closed)
    private, _ = flux_surface(eqt, eqt.profiles_1d.psi[end], :open)
    Z0 = sum(zlcfs) / length(zlcfs)
    empty!(eqt.boundary.x_point)

    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        elseif (sum(pz) < Z0)
            index = argmax(pz)
        elseif (sum(pz) > Z0)
            # upper private region
            index = argmin(pz)
        else
            continue
        end

        # add x-point info to the data structure
        resize!(eqt.boundary.x_point, length(eqt.boundary.x_point) + 1)
        eqt.boundary.x_point[end].r = pr[index]
        eqt.boundary.x_point[end].z = pz[index]

    end

    # if only one x-point is found, then look on the other side of the magnetic axis to find the other one
    if length(eqt.boundary.x_point) == 1
        push!(eqt.boundary.x_point, deepcopy(eqt.boundary.x_point[1]))
        eqt.boundary.x_point[2].z = -(eqt.boundary.x_point[2].z - eqt.global_quantities.magnetic_axis.z) + eqt.global_quantities.magnetic_axis.z
    end

    psi_separatrix, _ = find_psi_boundary(eqt; raise_error_on_not_open=true) # psi at LCFS
    psi_axis_level = eqt.profiles_1d.psi[1] # psi value on axis 
    psi_sign = sign(psi_separatrix - psi_axis_level) # +1 if psi increases / -1 if psi decreases

    if !isempty(eqt.boundary.x_point)

        # refine x-points location and re-sort
        psidist_lcfs_xpoints = Float64[]
        r, z, PSI_interpolant = ψ_interpolant(eqt.profiles_2d)
        for (k, x_point) in enumerate(eqt.boundary.x_point)
            res = Optim.optimize(
                x -> Bp(PSI_interpolant, [x_point.r + x[1]], [x_point.z + x[2]])[1],
                [0.0, 0.0],
                Optim.NelderMead(),
                Optim.Options(; g_tol=1E-8)
            )
            x_point.r += res.minimizer[1]
            x_point.z += res.minimizer[2]

            # record the distance from this x-point to the separatrix
            push!(psidist_lcfs_xpoints, PSI_interpolant(x_point.r, x_point.z)[1] - psi_separatrix)

        end

        # find distances among pairs of x-points (d_x) and record which one is closest to each (i_x)
        r_x = [x_point.r for x_point in eqt.boundary.x_point]
        z_x = [x_point.z for x_point in eqt.boundary.x_point]
        d_x = r_x .* Inf
        i_x = zeros(Int, length(eqt.boundary.x_point))
        for (k1, (r1, z1)) in enumerate(zip(r_x, z_x))
            for (k2, (r2, z2)) in enumerate(zip(r_x, z_x))
                if k2 == k1
                    continue
                end
                d = sqrt((r1 - r2)^2 + (z1 - z2)^2)
                if d < d_x[k1]
                    d_x[k1] = d
                    i_x[k1] = k2
                end
            end
        end

        # NOTE: this handles cases with two x-points in the same spot
        i_x = i_x[d_x.<1E-3]
        d_x = d_x[d_x.<1E-3]
        index = sortperm(d_x)
        d_x = d_x[index[1:2:end]]
        i_x = i_x[index[1:2:end]]
        for k in reverse!(sort(i_x))
            deleteat!(eqt.boundary.x_point, k)
            deleteat!(psidist_lcfs_xpoints, k)
            deleteat!(z_x, k)
        end


        # remove x-points that have fallen on the magnetic axis 
        sign_closest = sign(psidist_lcfs_xpoints[argmin(abs.(psidist_lcfs_xpoints))])# sign of psi of closest X-point in psi to LCFS
        index = psidist_lcfs_xpoints .* psi_sign .>= (psi_sign - sign_closest * 1E-5) * psidist_lcfs_xpoints[argmin(abs.(psidist_lcfs_xpoints))]
        psidist_lcfs_xpoints = psidist_lcfs_xpoints[index]
        eqt.boundary.x_point = eqt.boundary.x_point[index]
        z_x = z_x[index]

        # sort a second time now by distance in psi
        index = sortperm(abs.(psidist_lcfs_xpoints))
        eqt.boundary.x_point = eqt.boundary.x_point[index]
        z_x = z_x[index]
        # save up to the x_point with Z coordinate opposite to first x point
        # save up to first index where z_x.*z_x[1].<0 is 1 
        eqt.boundary.x_point = eqt.boundary.x_point[1:argmax(z_x .* z_x[1] .< 0)]
    end

    return eqt.boundary.x_point
end

"""
    miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)

Returns R0, a, κ, δu, δl, ζou, ζol, ζil, ζiu of a contour
"""
function miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}, r_at_max_z::T, max_z::T, r_at_min_z::T, min_z::T, z_at_max_r::T, max_r::T, z_at_min_r::T, min_r::T) where {T<:Real}
    R0 = 0.5 * (max_r + min_r)
    a = 0.5 * (max_r - min_r)
    b = 0.5 * (max_z - min_z)
    κ = b / a
    δu = (R0 - r_at_max_z) / a
    δl = (R0 - r_at_min_z) / a
    ζou, ζol, ζil, ζiu = luce_squareness(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
    return R0, a, κ, δu, δl, ζou, ζol, ζil, ζiu
end

function miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}) where {T<:Real}
    (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)
    return miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
end

"""
    fluxsurface_extrema(pr::Vector{T}, pz::Vector{T}) where {T<:Real}

Returns extrema indexes and values of R,Z flux surfaces vectors:

    imaxr, iminr,
    imaxz, iminz,
    r_at_max_z, max_z,
    r_at_min_z, min_z,
    z_at_max_r, max_r,
    z_at_min_r, min_r
"""
function fluxsurface_extrema(pr::Vector{T}, pz::Vector{T}) where {T<:Real}
    _, imaxr = findmax(pr)
    _, iminr = findmin(pr)
    _, imaxz = findmax(pz)
    _, iminz = findmin(pz)
    r_at_max_z, max_z = pr[imaxz], pz[imaxz]
    r_at_min_z, min_z = pr[iminz], pz[iminz]
    z_at_max_r, max_r = pz[imaxr], pr[imaxr]
    z_at_min_r, min_r = pz[iminr], pr[iminr]
    return (imaxr, iminr, imaxz, iminz,
        r_at_max_z, max_z,
        r_at_min_z, min_z,
        z_at_max_r, max_r,
        z_at_min_r, min_r)
end

"""
    luce_squareness(pr::AbstractVector{T}, pz::AbstractVector{T}, r_at_max_z::T, max_z::T, r_at_min_z::T, min_z::T, z_at_max_r::T, max_r::T, z_at_min_r::T, min_r::T) where {T<:Real}

Squareness from: "An analytic functional form for characterization and generation of axisymmetric plasma boundaries"
T.C. Luce, Plasma Phys. Control. Fusion 55 (2013) http://dx.doi.org/10.1088/0741-3335/55/9/095009

Returns: zetaou, zetaol, zetail, zetaiu
"""
function luce_squareness(
    pr::AbstractVector{T}, pz::AbstractVector{T},
    r_at_max_z::T, max_z::T,
    r_at_min_z::T, min_z::T,
    z_at_max_r::T, max_r::T,
    z_at_min_r::T, min_r::T) where {T<:Real}

    # zetaou
    PO = (r_at_max_z, z_at_max_r)
    PE = (max_r, max_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaou = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetaol
    PO = (r_at_min_z, z_at_max_r)
    PE = (max_r, min_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaol = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetaiu
    PO = (r_at_max_z, z_at_min_r)
    PE = (min_r, max_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaiu = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetail
    PO = (r_at_min_z, z_at_min_r)
    PE = (min_r, min_z)
    PD = IMAS.intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz)[2][1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetail = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    return zetaou, zetaol, zetail, zetaiu
end

