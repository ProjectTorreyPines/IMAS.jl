using LinearAlgebra
using Plots

"""
    ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns r, z, and ψ interpolant named tuple
"""
function ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    z = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    return ψ_interpolant(r, z, eqt2d.psi)
end

function ψ_interpolant(r::AbstractRange{T}, z::AbstractRange{T}, psi::Matrix{T}) where {T<:Real}
    PSI_interpolant = Interpolations.cubic_spline_interpolation((r, z), psi; extrapolation_bc=Interpolations.Line())
    return (r=r, z=z, PSI_interpolant=PSI_interpolant)
end

function ψ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return ψ_interpolant(eqt2d)
end

function ψ_interpolant(eqt::IMAS.equilibrium__time_slice)
    return ψ_interpolant(eqt.profiles_2d)
end

function ψ_interpolant(dd::IMAS.dd)
    return ψ_interpolant(dd.equilibrium.time_slice[])
end

"""
    ρ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d{T}, phi_norm::T) where {T<:Real}

Returns r, z, and ρ interpolant named tuple
"""
function ρ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d{T}, phi_norm::T) where {T<:Real}
    r = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    z = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    return ρ_interpolant(r, z, sqrt.(abs.(eqt2d.phi ./ phi_norm)))
end

function ρ_interpolant(r::AbstractRange{T}, z::AbstractRange{T}, rho::Matrix{T}) where {T<:Real}
    RHO_interpolant = Interpolations.cubic_spline_interpolation((r, z), rho; extrapolation_bc=Interpolations.Line())
    return (r=r, z=z, RHO_interpolant=RHO_interpolant)
end

function ρ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d{T}}, phi_norm::T) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt2dv)
    return ρ_interpolant(eqt2d, phi_norm)
end

function ρ_interpolant(eqt::IMAS.equilibrium__time_slice)
    return ρ_interpolant(eqt.profiles_2d, eqt.profiles_1d.phi[end])
end

function ρ_interpolant(dd::IMAS.dd)
    return ρ_interpolant(dd.equilibrium.time_slice[])
end

"""
    Br_Bz(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns Br and Bz named tuple evaluated at r and z starting from ψ interpolant
"""
function Br_Bz(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    Z, R = meshgrid(z, r)
    return Br_Bz(PSI_interpolant, R, Z)
end

function Br_Bz(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return Br_Bz(eqt2d)
end

function Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::T, z::T) where {T<:Real}
    grad = Interpolations.gradient(PSI_interpolant, r, z)
    inv_twopi_r = 1.0 / (2π * r)
    Br = grad[2] * inv_twopi_r
    Bz = -grad[1] * inv_twopi_r
    return (Br=Br, Bz=Bz)
end

function Br_Bz(PSI_interpolant::Interpolations.AbstractInterpolation, r::Array{T}, z::Array{T}) where {T<:Real}
    # Check that r and z are the same size
    @assert size(r) == size(z)
    Br, Bz = similar(r), similar(r)
    for k in eachindex(r)
        Br[k], Bz[k] = Br_Bz(PSI_interpolant, r[k], z[k])
    end
    return (Br=Br, Bz=Bz)
end

"""
    Bp(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns Bp evaluated at r and z starting from ψ interpolant
"""
function Bp(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    Z, R = meshgrid(z, r)
    return Bp.(Ref(PSI_interpolant), R, Z)
end

function Bp(PSI_interpolant::Interpolations.AbstractInterpolation, r::T, z::T) where {T<:Real}
    Br, Bz = Br_Bz(PSI_interpolant, r, z)
    return sqrt(Br^2.0 + Bz^2.0)
end

function Bp(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return Bp(eqt2d)
end

"""
    find_psi_boundary(
        eqt::IMAS.equilibrium__time_slice{T},
        wall_r::AbstractVector{T},
        wall_z::AbstractVector{T};
        precision::Float64=1e-6,
        raise_error_on_not_open::Bool=true,
        raise_error_on_not_closed::Bool=true
    ) where {T<:Real}

Find psi value of the last-closed and first-open flux surface

Results are returned as a named tuple `(last_closed=..., first_open=...)`
"""
function find_psi_boundary(
    eqt::IMAS.equilibrium__time_slice{T},
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T};
    precision::Float64=1e-6,
    raise_error_on_not_open::Bool=true,
    raise_error_on_not_closed::Bool=true
) where {T<:Real}

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    original_psi_boundary = eqt.profiles_1d.psi[end]
    if eqt2d !== nothing
        dimR = IMAS.to_range(eqt2d.grid.dim1)
        dimZ = IMAS.to_range(eqt2d.grid.dim2)
        psi_axis = eqt.profiles_1d.psi[1]
        RA = eqt.global_quantities.magnetic_axis.r
        ZA = eqt.global_quantities.magnetic_axis.z
        return find_psi_boundary(dimR, dimZ, eqt2d.psi, psi_axis, original_psi_boundary, RA, ZA, wall_r, wall_z; precision, raise_error_on_not_open, raise_error_on_not_closed)
    else
        # closed boundary equilibrium should end up here
        return (last_closed=original_psi_boundary, first_open=nothing)
    end
end

function find_psi_boundary(
    dimR::Union{AbstractVector{T},AbstractRange{T}},
    dimZ::Union{AbstractVector{T},AbstractRange{T}},
    PSI::Matrix{T},
    psi_axis::T,
    original_psi_boundary::T,
    RA::T,
    ZA::T,
    fw_r::AbstractVector{T},
    fw_z::AbstractVector{T};
    PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
    precision::Float64=1e-6,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool) where {T<:Real}

    verbose = false

    # here we figure out the range of psi to use to find the psi boundary
    if !isempty(fw_r)
        psi_edge = PSI_interpolant.(fw_r, fw_z)
    else
        psi_edge = [PSI[1, :]; PSI[end, :]; PSI[:, 1]; PSI[:, end]]
    end
    if psi_axis < original_psi_boundary
        psi_edge0 = maximum(psi_edge)
    else
        psi_edge0 = minimum(psi_edge)
    end
    psirange_init = [psi_axis + (psi_edge0 - psi_axis) / 100.0, psi_edge0]

    if verbose
        @show psirange_init
        plot(; aspect_ratio=:equal)
        contour!(dimR, dimZ, transpose(PSI); color=:gray, clim=(min(psirange_init...), max(psirange_init...)))
        if !isempty(fw_r)
            plot!(fw_r, fw_z; color=:black, lw=2, label="")
        end
        display(plot!())
    end

    # innermost tentative flux surface (which should be closed!)
    surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psirange_init[1], :closed)
    if isempty(surface)
        if raise_error_on_not_closed
            error("Flux surface at ψ=$(psirange_init[1]) is not closed; ψ=[$(psirange_init[1])...$(psirange_init[end])]")
        else
            return (last_closed=nothing, first_open=nothing)
        end
    end
    if verbose
        for surf in surface
            plot!(surf.r, surf.z; color=:blue, label="")
        end
        display(plot!())
    end

    # outermost tentative flux surface (which should be open!)
    surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psirange_init[end], :closed)
    if !isempty(surface)
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open; ψ=[$(psirange_init[1])...$(psirange_init[end])]")
        else
            return (last_closed=nothing, first_open=nothing)
        end
    end
    if verbose
        surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psirange_init[end], :open)
        for surf in surface
            plot!(surf.r, surf.z; color=:red, label="")
        end
        display(plot!())
    end

    psirange = deepcopy(psirange_init)
    for k in 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psimid, :closed)
        # closed flux surface
        if !isempty(surface)
            ((pr, pz),) = surface
            if verbose
                display(plot!(pr, pz; label="", color=:green))
            end
            psirange[1] = psimid
            if (abs(psirange[end] - psirange[1]) / abs(psirange[end] + psirange[1]) / 2.0) < precision
                return (last_closed=psimid, first_open=psirange[end])
            end
            # open flux surface
        else
            psirange[end] = psimid
        end
    end

    return error("Could not find closed boundary between ψ=$(psirange_init[1]) and ψ=$(psirange_init[end])")
end

function find_psi_boundary(
    dimR::Union{AbstractVector{T},AbstractRange{T}},
    dimZ::Union{AbstractVector{T},AbstractRange{T}},
    PSI::Matrix{T},
    psi_axis::T,
    axis2bnd::Symbol,
    RA::T,
    ZA::T,
    fw_r::AbstractVector{T}=T[],
    fw_z::AbstractVector{T}=T[];
    PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
    precision::Float64=1e-6,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool) where {T<:Real}

    @assert axis2bnd in (:increasing, :decreasing)
    verbose = false

    # here we figure out the range of psi to use to find the psi boundary
    f_ext = (axis2bnd === :increasing) ? maximum : minimum

    if !isempty(fw_r)
        psi_edge0 = f_ext(PSI_interpolant.(fw_r[k], fw_z[k]) for k in eachindex(fw_r))
    else
        @views psi_edge0 = f_ext((f_ext(PSI[1, :]), f_ext(PSI[end, :]), f_ext(PSI[:, 1]), f_ext(PSI[:, end])))
    end
    psirange_init = StaticArrays.@MVector[psi_axis + (psi_edge0 - psi_axis) / 100.0, psi_edge0]

    if verbose
        @show psirange_init
        plot(; aspect_ratio=:equal)
        contour!(dimR, dimZ, transpose(PSI); color=:gray, clim=(min(psirange_init...), max(psirange_init...)))
        if !isempty(fw_r)
            plot!(fw_r, fw_z; color=:black, lw=2, label="")
        end
        display(plot!())
    end

    # innermost tentative flux surface (which should be closed!)
    surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psirange_init[1], :closed)
    if isempty(surface)
        if raise_error_on_not_closed
            error("Flux surface at ψ=$(psirange_init[1]) is not closed; ψ=[$(psirange_init[1])...$(psirange_init[end])]")
        else
            return (last_closed=nothing, first_open=nothing)
        end
    end
    if verbose
        for surf in surface
            plot!(surf.r, surf.z; color=:blue, label="")
        end
        display(plot!())
    end

    # outermost tentative flux surface (which should be open!)
    surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psirange_init[end], :closed)
    if !isempty(surface)
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open; ψ=[$(psirange_init[1])...$(psirange_init[end])]")
        else
            return (last_closed=nothing, first_open=nothing)
        end
    end
    if verbose
        surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psirange_init[end], :open)
        for surf in surface
            plot!(surf.r, surf.z; color=:red, label="")
        end
        display(plot!())
    end

    psirange = deepcopy(psirange_init)
    for k in 1:100
        psimid = (psirange[1] + psirange[end]) / 2.0
        surface = flux_surface(dimR, dimZ, PSI, RA, ZA, fw_r, fw_z, psimid, :closed)
        # closed flux surface
        if !isempty(surface)
            ((pr, pz),) = surface
            if verbose
                display(plot!(pr, pz; label="", color=:green))
            end
            psirange[1] = psimid
            if (abs(psirange[end] - psirange[1]) / abs(psirange[end] + psirange[1]) / 2.0) < precision
                return (last_closed=psimid, first_open=psirange[end])
            end
            # open flux surface
        else
            psirange[end] = psimid
        end
    end

    return error("Could not find closed boundary between ψ=$(psirange_init[1]) and ψ=$(psirange_init[end])")
end

"""
    find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice; type::Symbol=:not_diverted, precision::Float64=1E-7)

Returns psi of the second magentic separatrix. This relies only on eqt and finds the 2nd sep geometrically.
"""
function find_psi_2nd_separatrix(
    eqt::IMAS.equilibrium__time_slice{T},
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T};
    type::Symbol=:not_diverted,
    precision::Float64=1E-7
) where {T<:Real}
    psi_separatrix = eqt.profiles_1d.psi[end]
    surface = flux_surface(eqt, psi_separatrix, :open, wall_r, wall_z)

    # First check if we are in a double null configuration
    ZA = eqt.global_quantities.magnetic_axis.z
    for (r, z) in surface
        if isempty(r) || all(z .> ZA) || all(z .< ZA)
            continue
        end
        if (z[end] - ZA) * (z[1] - ZA) < 0
            #if double null, all open surfaces in the SOL start and finish in opposite sides of the midplane
            return psi_separatrix
        end
    end

    # Single null case
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    psi_axis = eqt.profiles_1d.psi[1] # psi value on axis
    psi_sign = sign(psi_separatrix - psi_axis) # +1 for increasing psi / -1 for decreasing psi
    psi_low = psi_separatrix
    if psi_sign > 0
        # increasing psi
        psi_up = maximum(eqt2d.psi)
    else
        # decresing psi
        psi_up = minimum(eqt2d.psi)
    end
    psi = (psi_up + psi_low) / 2.0
    err = Inf
    counter = 0
    counter_max = 50
    while abs(err) > precision && counter < counter_max
        surface = flux_surface(eqt, psi, :open, wall_r, wall_z)
        for (r, z) in surface
            if isempty(r) || all(z .> ZA) || all(z .< ZA)
                continue
            end

            if (z[end] - ZA) * (z[1] - ZA) > 0
                # the surface starts and finishes on the same side of the midplane; aka is diverted
                psi_low = psi
            else
                # the surface starts and finishes on opposite sides of the midplane
                psi_up = psi
            end
        end

        if err == abs(psi_up - psi_low) / abs(psi_low)
            # neither psi_up nor psi_low were updated; all surfaces in surface were discarded
            # update psi_up, because psi_low is intrinsically safe for the algorithm
            psi_up = psi
        end

        # update new error
        err = abs(psi_up - psi_low) / abs(psi_low)
        psi = (psi_up + psi_low) / 2.0

        counter = counter + 1
    end

    if type == :not_diverted
        return psi_up
    elseif type == :diverted
        return psi_low
    end
end

"""
    find_psi_last_diverted(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation;
        precision::Float64=1e-7)

Returns `psi_last_lfs, `psi_first_lfs_far`, and `null_within_wall`

psi_first_lfs_far will be the first surface inside OFL[:lfs_far]; psi_last_lfs will be the last surface inside OFL[:lfs]

Precision between the two is defined on the poloidal crossection area at the OMP (Psol*precision = power flowing between psi_first_lfs_far and psi_last_lfs ~ 0)
"""
function find_psi_last_diverted(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation;
    precision::Float64=1e-7)

    # if no wall in dd, psi_last diverted not defined
    if isempty(wall_r) || isempty(wall_z) || isempty(eqt.boundary.x_point)
        return (psi_last_lfs=NaN, psi_first_lfs_far=NaN, null_within_wall=true)
    end

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    Xpoint2 = [eqt.boundary.x_point[end].r, eqt.boundary.x_point[end].z]

    psi_axis = eqt.profiles_1d.psi[1] # psi value on axis
    psi_separatrix = eqt.profiles_1d.psi[end] # psi value last closed
    psi_first_open = find_psi_boundary(eqt, wall_r, wall_z; raise_error_on_not_open=true).first_open
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt, wall_r, wall_z) # psi second magnetic separatrix
    psi_sign = sign(psi_separatrix - psi_axis) # +1 incresing psi / -1 decreasing psi

    # intersect 2nd separatrix with wall, and look
    surface = flux_surface(eqt, psi_2ndseparatrix, :open, wall_r, wall_z)
    rz_intersects = StaticArrays.SVector{2,Float64}[]
    r_max = 0.0
    for (r, z) in surface
        crossings = intersection(r, z, wall_r, wall_z, 1E-6).crossings # find where flux surface crosses wall ("strike points" of surface)
        if isempty(crossings)
            continue
        end

        # save all intersections with wall
        append!(rz_intersects, crossings)

        r_max = max(r_max, maximum(r))
    end

    r_intersect = Float64[rr for (rr, zz) in rz_intersects]
    z_intersect = Float64[zz for (rr, zz) in rz_intersects]

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

    # for safety, and to simplify eventual debugging
    if length(r_intersect) != 2
        plot(wall_r, wall_z)
        for (r, z) in surface
            @show point_to_path_distance(r[end], z[end], wall_r, wall_z)
            plot!(r, z)
        end
        scatter!(r_intersect, z_intersect)
        display(plot!())
        @assert length(r_intersect) == 2
    end

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

    if isempty(eqt.boundary.strike_point)
        find_strike_points!(eqt, wall_r, wall_z, psi_first_open)
    end

    # First we treat the double null case:
    # if we have 4 strike points, it is a double null
    if length(eqt.boundary.strike_point) == 4
        # LDFS is the separatrix
        # psi_first_lfs_far is determined using the precision condition
        psi_first_lfs_far = psi_separatrix + psi_sign * precision * abs(psi_separatrix)
        return (psi_last_lfs=psi_separatrix, psi_first_open=psi_first_open, psi_first_lfs_far=psi_first_lfs_far, null_within_wall=null_within_wall)
    end

    # Single_null case
    # find the two surfaces `psi_first_lfs_far` and `psi_last_lfs` around the last diverted flux surface
    psi_2ndseparatrix_notdiverted = psi_2ndseparatrix
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt, wall_r, wall_z; type=:diverted)
    psi_first_lfs_far = psi_2ndseparatrix
    psi_last_lfs = psi_separatrix
    psi = (psi_first_lfs_far + psi_last_lfs) / 2

    counter_max = 50
    err = Inf
    for counter in 1:counter_max
        surface = flux_surface(eqt, psi, :open, wall_r, wall_z)

        for (r, z) in surface
            if isempty(r) || all(z .> ZA) || all(z .< ZA)
                continue
            end

            crossings = intersection(r, z, wall_r, wall_z).crossings # find where flux surface crosses wall ("strike points" of surface)

            if isempty(crossings)
                continue
            end

            # r and z coordiante of intersections with wall
            r_intersect = (cr[1] for cr in crossings)
            z_intersect = (cr[2] for cr in crossings)

            closest_to_strike_points = Int64[]
            for point in eqt.boundary.strike_point
                dist = (r_intersect .- point.r) .^ 2 .+ (z_intersect .- point.z) .^ 2
                push!(closest_to_strike_points, argmin(dist))
            end
            sort!(closest_to_strike_points)

            if !isempty(closest_to_strike_points)
                # discard crossings occuring after second strike point
                if closest_to_strike_points[end] < length(crossings)
                    for k in reverse(closest_to_strike_points[end]+1:length(crossings))
                        deleteat!(crossings, k)
                    end
                end
                # discard crossings occuring before first strike point
                if closest_to_strike_points[1] > 1
                    for k in reverse(1:closest_to_strike_points[1]-1)
                        deleteat!(crossings, k)
                    end
                end
            end

            if length(crossings) == 2
                # psi corresonds to a diverted surface
                psi_last_lfs = psi
            else
                psi_first_lfs_far = psi
            end
        end

        # better to compute error on poloidal area between [psi_last_lfs, psi_first_lfs_far] (needed for accurate power balance)
        r_up = r_mid_itp(psi_first_lfs_far)
        r_low = r_mid_itp(psi_last_lfs)

        A = π * (r_up^2 - r_low^2) # annular area between r_up and r_low [m²]

        if err == abs(A)
            # no psi was updated; update psi_first_lfs_far because psi_last_lfs is safe
            psi_first_lfs_far = psi
        end

        err = abs(A)
        psi = (psi_first_lfs_far + psi_last_lfs) / 2

        if abs(err) < precision
            break
        end
    end

    # if LDFS is the 2nd separatrix, be consistent with the find_psi_2nd_separatrix function
    if psi_first_lfs_far == psi_2ndseparatrix
        psi_last_lfs = psi_first_lfs_far
        psi_first_lfs_far = psi_2ndseparatrix_notdiverted
    end

    return (psi_last_lfs=psi_last_lfs, psi_first_open=psi_first_open, psi_first_lfs_far=psi_first_lfs_far, null_within_wall=null_within_wall)
end

"""
    find_psi_tangent_omp(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation;
        precision::Float64=1e-7)

Returns the psi of the magnetic surface in the SOL which is tangent to the wall near the outer midplane
"""
function find_psi_tangent_omp(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation;
    precision::Float64=1e-7)

    # if no wall in dd, psi_last diverted not defined
    if isempty(wall_r) || isempty(wall_z)
        return (psi_tangent_in=NaN, psi_tangent_out=NaN)
    end

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    # find intersection of the wall with the outer midlpane
    crossings2 = intersection([RA, maximum(wall_r) * 2], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_omp = [cr[1] for cr in crossings2] # R coordinate of the wall at OMP
    r_wall_omp = r_wall_omp[1] # make it float

    psi_axis = eqt.profiles_1d.psi[1] # psi value on axis
    psi_separatrix = eqt.profiles_1d.psi[end] # psi LCFS
    psi_sign = sign(psi_separatrix - psi_axis) # +1 incresing psi / -1 decreasing psi

    ((_, zsep),) = flux_surface(eqt, psi_separatrix, :closed, wall_r, wall_z)

    b = maximum(abs.([minimum(zsep), maximum(zsep)]))

    wallr = wall_r[wall_r.>RA.&&abs.(wall_z).<b/2]
    wallz = wall_z[wall_r.>RA.&&abs.(wall_z).<b/2]

    # add intersection at +/- b/2 meter
    crossings = intersection([RA, 3 * RA], [1.0, 1.0] * b / 2, wall_r, wall_z).crossings # find where flux surface crosses wall around MP
    wallr = vcat(crossings[1][1], wallr)
    wallz = vcat(crossings[1][2], wallz)
    crossings = intersection([RA, 3 * RA], [-1.0, -1.0] * b / 2, wall_r, wall_z).crossings # find where flux surface crosses wall around MP
    push!(wallr, crossings[1][1])
    push!(wallz, crossings[1][2])

    # find the two surfaces `psi_wall_up` and `psi_wall_low` around the tangent surface to the wall at OMP
    if psi_sign == 1
        psi_tangent_out = maximum(PSI_interpolant.(wallr, wallz)) # psi increases
    else
        psi_tangent_out = minimum(PSI_interpolant.(wallr, wallz)) # psi decreases
    end
    psi_tangent_in = psi_separatrix
    psi = (psi_tangent_out + psi_tangent_in) / 2

    counter_max = 50
    err = Inf
    for counter in 1:counter_max
        surface = flux_surface(eqt, psi, :open, wall_r, wall_z)

        for (r, z) in surface
            # exclude empty vectors, surfaces that do not cross the midplane ans surfaces that cross the midplane at the HFS
            if isempty(r) || all(z .> ZA) || all(z .< ZA) || sum(r .> RA) == 0
                continue
            end

            if maximum(r) > r_wall_omp
                psi_tangent_out = psi
                psi = (psi_tangent_out + psi_tangent_in) / 2
            else
                psi_tangent_in = psi
                psi = (psi_tangent_out + psi_tangent_in) / 2
            end
        end

        err = abs(psi_tangent_out - psi_tangent_in) / abs(psi_tangent_in)
        if abs(err) < precision
            break
        end
    end

    return (psi_tangent_in=psi_tangent_in, psi_tangent_out=psi_tangent_out)
end

"""
    find_psi_max(
        eqt::IMAS.equilibrium__time_slice{T},
        wall_r::AbstractVector{T},
        wall_z::AbstractVector{T};
        precision::Float64=1e-2) where {T<:Real}

Returns the max psi useful for an ofl in the SOL with no wall.
"""
function find_psi_max(
    eqt::IMAS.equilibrium__time_slice{T},
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T};
    precision::Float64=1e-2) where {T<:Real}

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt, wall_r, wall_z)
    psi_axis = eqt.profiles_1d.psi[1] # psi value on axis
    psi_sign = sign(psi_2ndseparatrix - psi_axis) # +1 for increasing psi / -1 for decreasing psi

    # find the two surfaces `psi_first_lfs_far` and `psi_last_lfs` around the last diverted flux surface
    if psi_sign > 0
        # increasing psi
        psi_up = maximum(eqt2d.psi)
    else
        # decresing psi
        psi_up = minimum(eqt2d.psi)
    end

    # check first if psi_up is already the psi we want
    surface = flux_surface(eqt, psi_up, :open, wall_r, wall_z)
    for (r, z) in surface

        if isempty(r) || all(z .> ZA) || all(z .< ZA)
            continue
        end

        crossings = intersection(r, z, [1, 10] * RA, [1, 1] * ZA).crossings # find intersection with midplane

        if isempty(crossings)
            continue
        end

        rr = crossings[1][1] # R coordiante of intersection
        zz = crossings[1][2] # Z coordiante of intersection

        if rr < RA
            # surface is in :hfs
            continue
        end

        if z[1] * z[end] < 0
            #surface starts and finishes in opposite sides of the midplane is the surface we want
            return psi_up
        end

    end

    # psi_up corresponds to a surface that does not start and finish on opposite sides of the midplane
    psi_low = psi_2ndseparatrix # psi_low starts and finishes on opposite sides of the midlpane
    psi = (psi_up + psi_low) / 2.0
    zmax = maximum(eqt2d.grid.dim2)
    zmin = minimum(eqt2d.grid.dim2)

    counter_max = 50
    err = Inf
    for counter in 1:counter_max
        surface = flux_surface(eqt, psi, :open, wall_r, wall_z)
        for (r, z) in surface

            if isempty(r) || all(z .> ZA) || all(z .< ZA)
                continue
            end

            crossings = intersection(r, z, [1, 10] * RA, [1, 1] * ZA).crossings # find intersection of surface with midplane

            if isempty(crossings)
                continue
            end

            rr = crossings[1][1] # R coordiante of intersection
            zz = crossings[1][2] # Z coordiante of intersection

            if rr < RA
                # surface is in :hfs
                continue
            end

            if z[end] * z[1] < 0
                # psi starts and finishes on opposite sides of the midplane
                if abs(z[end] * z[1] - zmax * zmin) / (zmax^2) < 0.1
                    # make sure that indeed the surface goes from zmin to zmax, or close to that
                    psi_low = psi
                else
                    psi_up = psi
                end
            else
                # psi does not start and finish on opposite sides of the midplane
                psi_up = psi
            end
        end

        if err == abs(psi_up - psi_low) / abs(psi_low)
            # neither psi_up nor psi_low were updated; all surfaces in surface were discarded
            # update psi_up, because psi_low is intrinsically safe for the algorithm
            psi_up = psi
        end

        err = abs(psi_up - psi_low) / abs(psi_low)
        psi = (psi_up + psi_low) / 2.0

        if abs(err) < precision
            break
        end
    end


    return psi_low
end

"""
    find_psi_wall_omp(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real})

Returns the psi of the magnetic surface in the SOL which intersects the wall at the outer midplane
"""
function find_psi_wall_omp(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real})

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)

    return find_psi_wall_omp(PSI_interpolant, RA, ZA, wall_r, wall_z)
end

function find_psi_wall_omp(
    PSI_interpolant::Interpolations.AbstractInterpolation,
    RA::T1,
    ZA::T1,
    wall_r::AbstractVector{T2},
    wall_z::AbstractVector{T2}) where {T1<:Real,T2<:Real}

    crossings = intersection([RA, maximum(wall_r) * 1.1], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
    return PSI_interpolant.(r_wall_midplane, ZA)[1] # psi at the intersection between wall and omp
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
    flux_surfaces(eq::equilibrium{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Update flux surface averaged and geometric quantities in the equilibrium IDS
"""
function flux_surfaces(eq::equilibrium{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    for time_index in eachindex(eq.time_slice)
        flux_surfaces(eq.time_slice[time_index], wall_r, wall_z)
    end
    return eq
end

function find_magnetic_axis(r::AbstractVector{<:Real}, z::AbstractVector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation, psi_sign::Real;
    rguess::Real=r[Int(round(length(r) / 2))], zguess::Real=z[Int(round(length(z) / 2))])
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
        [rguess, zguess],
        Optim.Newton();
        autodiff=:forward
    )
    return res.minimizer[1], res.minimizer[2]
end

mutable struct FluxSurface{T}
    psi::T
    r::Vector{T}
    z::Vector{T}
    r_at_max_z::T
    max_z::T
    r_at_min_z::T
    min_z::T
    z_at_max_r::T
    max_r::T
    z_at_min_r::T
    min_r::T
    Br::Vector{T}
    Bz::Vector{T}
    Bp::Vector{T}
    Btot::Vector{T}
    ll::Vector{T}
    fluxexpansion::Vector{T}
    int_fluxexpansion_dl::T
end

@recipe function plot_FluxSurface(surface::FluxSurface)
    @series begin
        aspect_ratio := :equal
        surface.r, surface.z
    end
    @series begin
        primary := false
        aspect_ratio --> :equal
        seriestype := :scatter
        markerstrokewidth := 0
        markercolor --> :black
        markersize := 1
        [surface.r_at_max_z, surface.r_at_min_z, surface.max_r, surface.min_r], [surface.max_z, surface.min_z, surface.z_at_max_r, surface.z_at_min_r]
    end
end

@recipe function plot_FluxSurfaces(surfaces::Vector{FluxSurface})
    for k in eachindex(surfaces)
        @series begin
            label --> ""
            primary := k == 1
            surfaces[k]
        end
    end
end

function trace_surfaces(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z
    return trace_surfaces(eqt.profiles_1d.psi, eqt.profiles_1d.f, r, z, eqt2d.psi, eqt2d.b_field_r, eqt2d.b_field_z, PSI_interpolant, RA, ZA, wall_r, wall_z)
end

function trace_surfaces(
    psi::AbstractVector{T},
    f::AbstractVector{T},
    r::AbstractVector{T},
    z::AbstractVector{T},
    PSI::Matrix{T},
    BR::Matrix{T},
    BZ::Matrix{T},
    PSI_interpolant::Interpolations.AbstractInterpolation,
    RA::T,
    ZA::T,
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T}) where {T<:Real}

    surfaces = Dict{Int,FluxSurface}()
    N = length(psi)
    for k in N:-1:1
        psi_level = psi[k]

        # trace surfaces
        if k == 1 # on axis flux surface is a artificial one, generated from the second surface
            pr = (surfaces[2].r .- RA) ./ 100.0 .+ RA
            pz = (surfaces[2].z .- ZA) ./ 100.0 .+ ZA

        else  # other flux surfaces
            # trace flux surface
            tmp = flux_surface(r, z, PSI, RA, ZA, wall_r, wall_z, psi_level, :closed)
            if isempty(tmp)
                # p = heatmap(r, z, PSI'; colorbar=true, aspect_ratio=:equal)
                # contour!(r, z, PSI'; color=:white, levels=100)
                # contour!(r, z, PSI'; levels=[psi[end]], color=:white, lw=2)
                # display(p)
                error("IMAS: Could not trace closed flux surface $k out of $(N) at ψ = $(psi_level)")
            end
            (pr, pz) = tmp[1]
        end

        # surface length
        dl = vcat(0.0, sqrt.(diff(pr) .^ 2 + diff(pz) .^ 2))
        ll = cumsum(dl)

        # poloidal magnetic field (with sign)
        Br, Bz = Br_Bz(PSI_interpolant, pr, pz)
        tmp = Br .^ 2.0 .+ Bz .^ 2.0 #Bp2
        Bp_abs = sqrt.(tmp)
        Bp = Bp_abs .* sign.((pz .- ZA) .* Br .- (pr .- RA) .* Bz)
        Btot = sqrt.(tmp .+ (f[k] ./ pr) .^ 2)

        # flux expansion
        fluxexpansion = 1.0 ./ Bp_abs
        int_fluxexpansion_dl = trapz(ll, fluxexpansion)

        # Extrema on array indices
        (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)

        # create 
        surfaces[k] = FluxSurface(
            psi_level,
            pr,
            pz,
            r_at_max_z,
            max_z,
            r_at_min_z,
            min_z,
            z_at_max_r,
            max_r,
            z_at_min_r,
            min_r,
            Br,
            Bz,
            Bp,
            Btot,
            ll,
            fluxexpansion,
            int_fluxexpansion_dl)
    end

    N2 = Int(ceil(N / 2))
    algorithm = Optim.NelderMead()

    # extrema in R
    lines = Contour.lines(Contour.contour(r, z, BR, 0.0))
    k = 0
    d = Inf
    for (kk, line) in enumerate(lines)
        pr, pz = Contour.coordinates(line)
        dd = minimum(sqrt.((pr .- surfaces[N2].max_r) .^ 2 .+ (pz .- surfaces[N2].z_at_max_r) .^ 2))
        if dd < d
            d = dd
            k = kk
        end
    end
    leftright_r, leftright_z = Contour.coordinates(lines[k])
    for k in [N2:-1:2; N2+1:N]
        cost = x -> extrema_cost(x, psi[k], PSI_interpolant, leftright_r, leftright_z, RA, ZA, :right)
        surfaces[k].max_r, surfaces[k].z_at_max_r = Optim.optimize(cost, [surfaces[k].max_r, surfaces[k].z_at_max_r], algorithm; autodiff=:forward).minimizer
        cost = x -> extrema_cost(x, psi[k], PSI_interpolant, leftright_r, leftright_z, RA, ZA, :left)
        surfaces[k].min_r, surfaces[k].z_at_min_r = Optim.optimize(cost, [surfaces[k].min_r, surfaces[k].z_at_min_r], algorithm; autodiff=:forward).minimizer
        if 2 < k <= N2
            surfaces[k-1].z_at_max_r = surfaces[k].z_at_max_r
            surfaces[k-1].z_at_min_r = surfaces[k].z_at_min_r
        elseif k < N
            surfaces[k+1].z_at_max_r = surfaces[k].z_at_max_r
            surfaces[k+1].z_at_min_r = surfaces[k].z_at_min_r
        end
    end

    # extrema in Z
    lines = Contour.lines(Contour.contour(r, z, BZ, 0.0))
    k = 0
    d = Inf
    for (kk, line) in enumerate(lines)
        pr, pz = Contour.coordinates(line)
        dd = minimum(sqrt.((pr .- surfaces[N2].r_at_max_z) .^ 2 .+ (pz .- surfaces[N2].max_z) .^ 2))
        if dd < d
            d = dd
            k = kk
        end
    end
    updown_r, updown_z = Contour.coordinates(lines[k])
    for k in [N2:-1:2; N2+1:N]
        cost = x -> extrema_cost(x, psi[k], PSI_interpolant, updown_r, updown_z, RA, ZA, :up)
        surfaces[k].r_at_max_z, surfaces[k].max_z = Optim.optimize(cost, [surfaces[k].r_at_max_z, surfaces[k].max_z], algorithm; autodiff=:forward).minimizer
        cost = x -> extrema_cost(x, psi[k], PSI_interpolant, updown_r, updown_z, RA, ZA, :down)
        surfaces[k].r_at_min_z, surfaces[k].min_z = Optim.optimize(cost, [surfaces[k].r_at_min_z, surfaces[k].min_z], algorithm; autodiff=:forward).minimizer
        if 2 < k <= N2
            surfaces[k-1].r_at_max_z = surfaces[k].r_at_max_z
            surfaces[k-1].r_at_min_z = surfaces[k].r_at_min_z
        elseif k < N
            surfaces[k+1].r_at_max_z = surfaces[k].r_at_max_z
            surfaces[k+1].r_at_min_z = surfaces[k].r_at_min_z
        end
    end

    # first flux surface just a scaled down version of the second one
    surfaces[1].r_at_max_z = (surfaces[2].r_at_max_z .- RA) ./ 100.0 .+ RA
    surfaces[1].max_z = (surfaces[2].max_z .- ZA) ./ 100.0 .+ ZA
    surfaces[1].r_at_min_z = (surfaces[2].r_at_min_z .- RA) ./ 100.0 .+ RA
    surfaces[1].min_z = (surfaces[2].min_z .- ZA) ./ 100.0 .+ ZA
    surfaces[1].z_at_max_r = (surfaces[2].z_at_max_r .- ZA) ./ 100.0 .+ ZA
    surfaces[1].max_r = (surfaces[2].max_r - RA) / 100.0 + RA
    surfaces[1].z_at_min_r = (surfaces[2].z_at_min_r .- ZA) ./ 100.0 .+ ZA
    surfaces[1].min_r = (surfaces[2].min_r - RA) / 100.0 + RA

    return FluxSurface[surfaces[k] for k in 1:N]
end

# accurate geometric quantities by finding geometric extrema as optimization problem
function extrema_cost(
    x::AbstractVector{<:Real},
    psi_level::T,
    PSI_interpolant,
    r_rail::AbstractVector{T},
    z_rail::AbstractVector{T},
    RA::T,
    ZA::T,
    direction::Symbol
) where {T<:Real}
    d = point_to_path_distance(x[1], x[2], r_rail, z_rail)
    if direction == :right
        cost = (x[1] - RA)^2 * (x[1] < RA)
    elseif direction == :left
        cost = (x[1] - RA)^2 * (x[1] > RA)
    elseif direction == :up
        cost = (x[2] - ZA)^2 * (x[2] < ZA)
    elseif direction == :down
        cost = (x[2] - ZA)^2 * (x[2] > ZA)
    end
    cost += (PSI_interpolant(x[1], x[2]) - psi_level)^2 + 0.1 * d^2
    return cost
end

"""
    flux_surfaces(eqt::equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice.
"""
function flux_surfaces(eqt::equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    PSI = eqt2d.psi

    psi_sign = sign(eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1])
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    # ensure certain global quantities are consistent with 1d profiles by making them expressions
    for field in (:psi_boundary, :psi_axis, :q_95, :q_axis, :q_min)
        empty!(eqt.global_quantities, field)
    end

    # accurately find magnetic axis and lcfs and scale psi accordingly
    RA, ZA = find_magnetic_axis(r, z, PSI_interpolant, psi_sign)
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = RA, ZA
    psi_axis = PSI_interpolant(RA, ZA)
    original_psi_boundary = eqt.profiles_1d.psi[end]
    psi_boundaries =
        find_psi_boundary(r, z, eqt2d.psi, psi_axis, original_psi_boundary, RA, ZA, wall_r, wall_z; PSI_interpolant, raise_error_on_not_open=false, raise_error_on_not_closed=false)

    find_strike_points!(eqt, wall_r, wall_z, psi_boundaries.first_open)
    eqt.profiles_1d.psi =
        (eqt.profiles_1d.psi .- eqt.profiles_1d.psi[1]) ./ (eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1]) .* (psi_boundaries.last_closed - psi_axis) .+ psi_axis

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
        :gm10,
        :fsa_bp,
        :trapped_fraction
    )
        setproperty!(eqt.profiles_1d, item, zeros(eltype(eqt.profiles_1d.psi), size(eqt.profiles_1d.psi)))
    end

    # trace flux surfaces
    surfaces = trace_surfaces(eqt.profiles_1d.psi, eqt.profiles_1d.f, r, z, eqt2d.psi, eqt2d.b_field_r, eqt2d.b_field_z, PSI_interpolant, RA, ZA, wall_r, wall_z)

    # calculate flux surface averaged and geometric quantities
    N = length(eqt.profiles_1d.psi)

    Np = maximum(length(surface.r) for surface in surfaces)
    shared_tmp = Vector{T}(undef, Np)
    for k in N:-1:1
        surface = surfaces[k]
        pr = surface.r
        pz = surface.z
        Btot = surface.Btot
        Bp = surface.Bp

        tmp = @views shared_tmp[1:length(pr)]

        eqt.profiles_1d.r_outboard[k] = surface.max_r
        eqt.profiles_1d.r_inboard[k] = surface.min_r

        # miller geometric coefficients
        _, _, κ, δu, δl, ζou, ζol, ζil, ζiu =
            miller_R_a_κ_δ_ζ(pr, pz, surface.r_at_max_z, surface.max_z, surface.r_at_min_z, surface.min_z, surface.z_at_max_r, surface.max_r, surface.z_at_min_r, surface.min_r)
        eqt.profiles_1d.elongation[k] = κ
        eqt.profiles_1d.triangularity_upper[k] = δu
        eqt.profiles_1d.triangularity_lower[k] = δl
        eqt.profiles_1d.squareness_lower_outer[k] = ζol
        eqt.profiles_1d.squareness_upper_outer[k] = ζou
        eqt.profiles_1d.squareness_lower_inner[k] = ζil
        eqt.profiles_1d.squareness_upper_inner[k] = ζiu

        # trapped fraction
        Bmin = minimum(Btot)
        Bmax = maximum(Btot)
        avg_Btot = flux_surface_avg(Btot, surface)
        tmp .= Btot .^ 2
        avg_Btot2 = flux_surface_avg(tmp, surface)
        tmp .= Btot ./ Bmax
        tmp .= (1.0 .- sqrt.(1.0 .- tmp) .* (1.0 .+ tmp ./ 2.0)) ./ tmp .^ 2
        hf = flux_surface_avg(tmp, surface)
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
        tmp .= 1.0 ./ pr .^ 2
        eqt.profiles_1d.gm1[k] = flux_surface_avg(tmp, surface)

        # gm4 = <1/B^2>
        tmp .= 1.0 ./ Btot .^ 2
        eqt.profiles_1d.gm4[k] = flux_surface_avg(tmp, surface)

        # gm5 = <B^2>
        eqt.profiles_1d.gm5[k] = avg_Btot2

        # gm8 = <R>
        eqt.profiles_1d.gm8[k] = flux_surface_avg(pr, surface)

        # gm9 = <1/R>
        tmp .= 1.0 ./ pr
        eqt.profiles_1d.gm9[k] = flux_surface_avg(tmp, surface)

        # gm10 = <R^2>
        tmp .= pr .^ 2
        eqt.profiles_1d.gm10[k] = flux_surface_avg(tmp, surface)

        # fsa_bp = <Bp>
        eqt.profiles_1d.fsa_bp[k] = flux_surface_avg(Bp, surface)

        # j_tor = <j_tor/R> / <1/R> [A/m²]
        eqt.profiles_1d.j_tor[k] =
            (
                -(eqt.profiles_1d.dpressure_dpsi[k] + eqt.profiles_1d.f_df_dpsi[k] * eqt.profiles_1d.gm1[k] / constants.μ_0) *
                (2π)
            ) / eqt.profiles_1d.gm9[k]

        # dvolume_dpsi
        eqt.profiles_1d.dvolume_dpsi[k] = sign(eqt.profiles_1d.fsa_bp[k]) * surface.int_fluxexpansion_dl

        # surface area
        eqt.profiles_1d.surface[k] = 2π * trapz(surface.ll, pr)

        # q
        eqt.profiles_1d.q[k] = eqt.profiles_1d.dvolume_dpsi[k] * eqt.profiles_1d.f[k] * eqt.profiles_1d.gm1[k] / (2π)

        # quantities calculated on the last closed flux surface
        if k == N
            # perimeter
            eqt.global_quantities.length_pol = surface.ll[end]

            # boundary
            eqt.boundary.outline.r = pr
            eqt.boundary.outline.z = pz
        end
    end

    tmp = similar(eqt.profiles_1d.psi)

    # area
    tmp .= eqt.profiles_1d.dvolume_dpsi .* eqt.profiles_1d.gm9 ./ 2π
    eqt.profiles_1d.area = cumtrapz(eqt.profiles_1d.psi, tmp)
    eqt.global_quantities.area = eqt.profiles_1d.area[end]

    # volume
    eqt.profiles_1d.volume = cumtrapz(eqt.profiles_1d.psi, eqt.profiles_1d.dvolume_dpsi)
    eqt.global_quantities.volume = eqt.profiles_1d.volume[end]

    # phi
    tmp .= eqt.profiles_1d.f .* eqt.profiles_1d.gm1 / (2π)
    eqt.profiles_1d.phi = cumtrapz(eqt.profiles_1d.volume, tmp)

    # rho_tor_norm
    rho = sqrt.(abs.(eqt.profiles_1d.phi ./ (π * B0)))
    rho_meters = rho[end]
    eqt.profiles_1d.rho_tor = rho
    eqt.profiles_1d.rho_tor_norm = rho ./ rho_meters

    # phi 2D
    tmp .= eqt.profiles_1d.psi .* psi_sign
    phi_itp = interp1d(tmp, eqt.profiles_1d.phi, :cubic)
    eqt2d.phi = phi_itp.(psi_sign .* eqt2d.psi)

    # rho 2D in meters
    RHO = sqrt.(abs.(eqt2d.phi ./ (π * B0)))

    # gm2: <∇ρ²/R²>
    dRHOdR, dRHOdZ = gradient(collect(r), collect(z), RHO)
    dPHI2_interpolant = Interpolations.cubic_spline_interpolation((r, z), dRHOdR .^ 2.0 .+ dRHOdZ .^ 2.0)
    dPHI2_R2 = Vector{T}(undef, Np)
    for k in eachindex(surfaces)
        surface = surfaces[k]
        n = length(surface.r)
        dPHI2_R2[1:n] .= dPHI2_interpolant.(surface.r, surface.z) ./ surface.r .^ 2.0
        @views eqt.profiles_1d.gm2[k] = flux_surface_avg(dPHI2_R2[1:n], surface)
    end
    @views gm2_itp = interp1d(tmp[2:end], eqt.profiles_1d.gm2[2:end], :cubic)
    eqt.profiles_1d.gm2[1] = gm2_itp(tmp[1])

    # ip
    eqt.global_quantities.ip = trapz(eqt.profiles_1d.area, eqt.profiles_1d.j_tor)
    # eqt.global_quantities.ip = trapz(eqt.profiles_1d.volume, eqt.profiles_1d.j_tor.*eqt.profiles_1d.gm9) / (2π) # equivalent

    # Geometric major and minor radii
    Rgeo = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
    a = (eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]) / 2.0

    # vacuum magnetic field at the geometric center
    Btvac = B0 * R0 / Rgeo

    # average poloidal magnetic field
    Bpave = eqt.global_quantities.ip * constants.μ_0 / eqt.global_quantities.length_pol

    # li
    Bp2v = trapz(eqt.profiles_1d.psi, T[trapz(surface.ll, surface.Bp) for surface in surfaces])
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
    find_x_point!(eqt, wall_r, wall_z)

    # secondary separatrix
    if length(eqt.boundary.x_point) > 1
        psi2nd = find_psi_2nd_separatrix(eqt, wall_r, wall_z)
        pts = flux_surface(r, z, eqt2d.psi, RA, ZA, wall_r, wall_z, psi2nd, :encircling)
        if !isempty(pts)
            (pr2nd, pz2nd) = pts[1]
            eqt.boundary_secondary_separatrix.outline.r = pr2nd
            eqt.boundary_secondary_separatrix.outline.z = pz2nd
        end
    end

    return eqt
end

function flux_surfaces(eqt::equilibrium__time_slice{T}, wall::IMAS.wall) where {T<:Real}
    fw = IMAS.first_wall(wall)
    return flux_surfaces(eqt, fw.r, fw.z)
end

"""
    flux_surface(eqt::equilibrium__time_slice{T}, psi_level::Real, type::Symbol, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Returns a vector with the (r,z) coordiates of flux surface at given psi_level

The `type` parameter:

  - :any, return all contours
  - :closed, all closed flux-surface that encircle the magnetic axis and do not cross the wall
  - :open, all open flux-surfaces (considerning open even closed flux surfaces that hit the first wall)
  - :open_no_wall, all open flux-surfaces independently of wall
  - :encircling, open flux-surfaces encircling the magnetic axis
"""
function flux_surface(eqt::equilibrium__time_slice{T}, psi_level::Real, type::Symbol, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    dim1 = IMAS.to_range(eqt2d.grid.dim1)
    dim2 = IMAS.to_range(eqt2d.grid.dim2)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z
    return flux_surface(dim1, dim2, eqt2d.psi, RA, ZA, wall_r, wall_z, psi_level, type)
end

function flux_surface(
    dim1::AbstractVector{T},
    dim2::AbstractVector{T},
    PSI::AbstractArray{T},
    RA::T,
    ZA::T,
    fw_r::AbstractVector{T},
    fw_z::AbstractVector{T},
    psi_level::T,
    type::Symbol) where {T<:Real}

    # contouring routine
    cl = Contour.contour(dim1, dim2, PSI, psi_level; VT=NTuple{2,T})

    return flux_surface(dim1, dim2, cl, RA, ZA, fw_r, fw_z, type)
end

function flux_surface(
    dim1::Union{AbstractVector{T},AbstractRange{T}},
    dim2::Union{AbstractVector{T},AbstractRange{T}},
    cl::Contour.ContourLevel,
    RA::T,
    ZA::T,
    fw_r::AbstractVector{T},
    fw_z::AbstractVector{T},
    type::Symbol) where {T<:Real}

    prpz = NamedTuple{(:r, :z),Tuple{Vector{T},Vector{T}}}[]
    if type == :any
        # if no open/closed check, then return all contours
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            reorder_flux_surface!(pr, pz, RA, ZA; force_close=false)
            push!(prpz, (r=pr, z=pz))
        end

    elseif type == :closed
        # look for closed flux-surface
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surface that closes, contains magnetic axis, and does not intersect any wall element
            if (is_closed_polygon(pr, pz) && (PolygonOps.inpolygon((RA, ZA), collect(zip(pr, pz))) == 1) && !IMAS.intersects(pr, pz, fw_r, fw_z))
                push!(prpz, (r=pr, z=pz))
                break
            end
        end

    elseif type == :open_no_wall || (type == :open && isempty(fw_r))
        # look for open flux-surfaces
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # pick flux surfaces that do not close
            if is_open_polygon(pr, pz)
                reorder_flux_surface!(pr, pz, RA, ZA; force_close=false)
                push!(prpz, (r=pr, z=pz))
            end
        end

    elseif type == :open && !isempty(fw_r)
        # look for open flux-surfaces with wall
        fw = collect(zip(fw_r, fw_z))
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            # only lines that intersect with the wall are open
            if IMAS.intersects(pr, pz, fw_r, fw_z)
                reorder_flux_surface!(pr, pz, RA, ZA; force_close=false)
                segments = intersection_split(pr, pz, fw_r, fw_z)
                for segment in segments
                    # we retain only segments that are within the wall (disregard the extrema)
                    if length(segment.r) > 2 && PolygonOps.inpolygon((segment.r[2], segment.z[2]), fw) == 1
                        push!(prpz, segment)
                    end
                end
            end
        end

    elseif type == :encircling
        # look for open flux-surfaces that encircle the magnetic axis
        for (pr, pz) in flux_surface(dim1, dim2, cl, RA, ZA, fw_r, fw_z, :open)
            # close open flux surface
            tmp = collect(zip(pr, pz))
            push!(tmp, tmp[1])
            # and see if it contains the magnetic axis
            if PolygonOps.inpolygon((RA, ZA), tmp) == 1
                reorder_flux_surface!(pr, pz, RA, ZA; force_close=false)
                push!(prpz, (r=pr, z=pz))
            end
        end

    else
        error("flux_surface type `$type` is not recognized. It can be one of [:any, :closed, :open, :open_no_wall, :encircling]")
    end

    return prpz
end

"""
    flux_surface_avg(quantity::AbstractVector{T}, surface::FluxSurface{T}) where {T<:Real}

Flux surface averaging of a quantity
"""
function flux_surface_avg(quantity::AbstractVector{T}, surface::FluxSurface{T}) where {T<:Real}
    f = (k, xx) -> quantity[k] * surface.fluxexpansion[k]
    return trapz(surface.ll, f) / surface.int_fluxexpansion_dl
end

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over volume
"""
function volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    dv = eqt.profiles_1d.dvolume_dpsi
    f = (k, xx) -> dv[k] * what[k]
    return trapz(eqt.profiles_1d.psi, f)
end

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over volume
"""
function cumlul_volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    dv = eqt.profiles_1d.dvolume_dpsi
    f = (k, xx) -> dv[k] * what[k]
    return cumtrapz(eqt.profiles_1d.psi, f)
end

"""
    surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over surface
"""
function surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    dv = eqt.profiles_1d.dvolume_dpsi
    gm9 = eqt.profiles_1d.gm9
    f = (k, xx) -> dv[k] * what[k] * gm9[k]
    return trapz(eqt.profiles_1d.psi, f) / 2π
end

"""
    surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over surface
"""
function cumlul_surface_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    dv = eqt.profiles_1d.dvolume_dpsi
    gm9 = eqt.profiles_1d.gm9
    f = (k, xx) -> dv[k] * what[k] * gm9[k] / (2π)
    return cumtrapz(eqt.profiles_1d.psi, f)
end

"""
    find_x_point!(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Find the `n` X-points that are closest to the separatrix
"""
function find_x_point!(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    ((rlcfs, zlcfs),) = flux_surface(eqt, eqt.profiles_1d.psi[end], :closed, wall_r, wall_z)
    private = flux_surface(eqt, eqt.profiles_1d.psi[end], :open, wall_r, wall_z)
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
        resize!(eqt.boundary.x_point, 2)
        eqt.boundary.x_point[2].r = eqt.boundary.x_point[1].r
        eqt.boundary.x_point[2].z = -(eqt.boundary.x_point[1].z - eqt.global_quantities.magnetic_axis.z) + eqt.global_quantities.magnetic_axis.z
    end

    psi_separatrix = eqt.profiles_1d.psi[end] # psi value at LCFS
    psi_axis_level = eqt.profiles_1d.psi[1] # psi value on axis
    psi_sign = sign(psi_separatrix - psi_axis_level) # +1 if psi increases / -1 if psi decreases

    if !isempty(eqt.boundary.x_point)
        # refine x-points location and re-sort
        psidist_lcfs_xpoints = Float64[]
        r, z, PSI_interpolant = ψ_interpolant(eqt.profiles_2d)
        for (k, x_point) in enumerate(eqt.boundary.x_point)
            res = Optim.optimize(
                x -> Bp.(Ref(PSI_interpolant), [x_point.r + x[1]], [x_point.z + x[2]])[1],
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
        min_psidist = psidist_lcfs_xpoints[argmin(abs(x) for x in psidist_lcfs_xpoints)]

        sign_closest = sign(min_psidist)# sign of psi of closest X-point in psi to LCFS
        index = psidist_lcfs_xpoints .* psi_sign .>= (psi_sign - sign_closest * 1E-5) * min_psidist
        psidist_lcfs_xpoints = psidist_lcfs_xpoints[index]
        eqt.boundary.x_point = eqt.boundary.x_point[index]
        z_x = z_x[index]

        # sort a second time now by distance in psi
        index = sortperm(psidist_lcfs_xpoints; by=abs)
        eqt.boundary.x_point = eqt.boundary.x_point[index]
        z_x = z_x[index]
        # save up to the x_point with Z coordinate opposite to first x point
        # save up to first index where z_x.*z_x[1].<0 is 1
        eqt.boundary.x_point = eqt.boundary.x_point[1:argmax(z_x .* z_x[1] .< 0)]
    end

    return eqt.boundary.x_point
end

"""
    x_points_inside_wall(x_points::IDSvector{<:IMAS.equilibrium__time_slice___boundary__x_point}, wall::IMAS.wall)

Returns vector of x_points that are inside of the first wall
"""
function x_points_inside_wall(x_points::IDSvector{T}, wall::IMAS.wall) where {T<:IMAS.equilibrium__time_slice___boundary__x_point}
    outline = first_wall(wall)
    if isempty(outline.r)
        return T[x_point for x_point in x_points]
    else
        outline = collect(zip(outline.r, outline.z))
        indexes = findall(x_point -> PolygonOps.inpolygon((x_point.r, x_point.z), outline) == 1, x_points)
        return T[x_points[index] for index in indexes]
    end
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
    PD = intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz).crossings[1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaou = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetaol
    PO = (r_at_min_z, z_at_max_r)
    PE = (max_r, min_z)
    PD = intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz).crossings[1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaol = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetaiu
    PO = (r_at_max_z, z_at_min_r)
    PE = (min_r, max_z)
    PD = intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz).crossings[1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetaiu = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    # zetail
    PO = (r_at_min_z, z_at_min_r)
    PE = (min_r, min_z)
    PD = intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz).crossings[1]
    PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
    zetail = (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC)

    return zetaou, zetaol, zetail, zetaiu
end
