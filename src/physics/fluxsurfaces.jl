document[Symbol("Physics flux-surfaces")] = Symbol[]

using LinearAlgebra
using Plots

# setting global variable for precision for computing relevant surfaces (i.e. lcfs, first_open, ldfs, 1st sep, 2nd sep,..)
const flux_surfaces_precision::Float64 = 1E-6

"""
    find_psi_boundary(
        dimR::Union{AbstractVector{T1},AbstractRange{T1}},
        dimZ::Union{AbstractVector{T1},AbstractRange{T1}},
        PSI::Matrix{T2},
        psi_axis::Real,
        original_psi_boundary::Real,
        RA::T3,
        ZA::T3,
        fw_r::AbstractVector{T4},
        fw_z::AbstractVector{T5};
        PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
        precision::Float64=flux_surfaces_precision,
        raise_error_on_not_open::Bool,
        raise_error_on_not_closed::Bool,
        verbose::Bool=false
    ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}
"""
function find_psi_boundary(
    dimR::Union{AbstractVector{T1},AbstractRange{T1}},
    dimZ::Union{AbstractVector{T1},AbstractRange{T1}},
    PSI::Matrix{T2},
    psi_axis::Real,
    original_psi_boundary::Real,
    RA::T3,
    ZA::T3,
    fw_r::AbstractVector{T4},
    fw_z::AbstractVector{T5},
    r_cache::AbstractVector{T1}=T1[],
    z_cache::AbstractVector{T1}=T1[];
    PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
    precision::Float64=flux_surfaces_precision,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool,
    verbose::Bool=false
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}

    # here we figure out the range of psi to use to find the psi boundary
    if isempty(fw_r)
        @views psi_edge = [PSI[1, :]; PSI[end, :]; PSI[:, 1]; PSI[:, end]]
    else
        psi_edge = PSI_interpolant.(fw_r, fw_z)
    end
    psi_edge0 = original_psi_boundary + (original_psi_boundary - psi_axis)
    # BCL 11/22/24: Need to handle local extrema (coils) inside domain
    if psi_axis < original_psi_boundary
        psi_edge1 = maximum(psi_edge)
        psi_edge_guess = min(psi_edge1, psi_edge0)
    else
        psi_edge1 = minimum(psi_edge)
        psi_edge_guess = max(psi_edge1, psi_edge0)
    end
    psirange_init = [psi_axis + (original_psi_boundary - psi_axis) / 100.0, psi_edge_guess]

    if verbose
        @show psi_axis
        @show original_psi_boundary
    end

    psi_boundaries = find_psi_boundary(
        dimR,
        dimZ,
        PSI,
        psirange_init,
        RA,
        ZA,
        fw_r,
        fw_z,
        r_cache,
        z_cache;
        PSI_interpolant,
        precision,
        raise_error_on_not_open,
        raise_error_on_not_closed,
        verbose
    )

    return psi_boundaries
end

"""
    find_psi_boundary(
        dimR::Union{AbstractVector{T1},AbstractRange{T1}},
        dimZ::Union{AbstractVector{T1},AbstractRange{T1}},
        PSI::Matrix{T2},
        psi_axis::Real,
        axis2bnd::Symbol,
        RA::T3,
        ZA::T3,
        fw_r::AbstractVector{T4}=T1[],
        fw_z::AbstractVector{T5}=T1[],
        r_cache::AbstractVector{T1}=T1[],
        z_cache::AbstractVector{T1}=T1[];
        PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
        precision::Float64=flux_surfaces_precision,
        raise_error_on_not_open::Bool,
        raise_error_on_not_closed::Bool,
        verbose::Bool=false
    ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}
"""
function find_psi_boundary(
    dimR::Union{AbstractVector{T1},AbstractRange{T1}},
    dimZ::Union{AbstractVector{T1},AbstractRange{T1}},
    PSI::Matrix{T2},
    psi_axis::Real,
    axis2bnd::Symbol,
    RA::T3,
    ZA::T3,
    fw_r::AbstractVector{T4}=T1[],
    fw_z::AbstractVector{T5}=T1[],
    r_cache::AbstractVector{<:Real}=T1[],
    z_cache::AbstractVector{<:Real}=T1[];
    PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
    precision::Float64=flux_surfaces_precision,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool,
    verbose::Bool=false
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}

    @assert axis2bnd in (:increasing, :decreasing)

    # here we figure out the range of psi to use to find the psi boundary
    f_ext = (axis2bnd === :increasing) ? maximum : minimum

    if !isempty(fw_r)
        psi_edge0 = f_ext(PSI_interpolant(fw_r[k], fw_z[k]) for k in eachindex(fw_r))
    else
        @views psi_edge0 = f_ext((f_ext(PSI[1, :]), f_ext(PSI[end, :]), f_ext(PSI[:, 1]), f_ext(PSI[:, end])))
    end
    psirange_init = [psi_axis + (psi_edge0 - psi_axis) / 100.0, psi_edge0]

    psi_boundaries = find_psi_boundary(
        dimR,
        dimZ,
        PSI,
        psirange_init,
        RA,
        ZA,
        fw_r,
        fw_z,
        r_cache,
        z_cache;
        PSI_interpolant,
        precision,
        raise_error_on_not_open,
        raise_error_on_not_closed,
        verbose
    )

    return psi_boundaries
end

"""
    find_psi_boundary(
        dimR::Union{AbstractVector{T1},AbstractRange{T1}},
        dimZ::Union{AbstractVector{T1},AbstractRange{T1}},
        PSI::Matrix{T2},
        psirange_init::AbstractVector,
        RA::T3,
        ZA::T3,
        fw_r::AbstractVector{T4},
        fw_z::AbstractVector{T5},
        r_cache::AbstractVector{T1}=T1[],
        z_cache::AbstractVector{T1}=T1[];
        PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
        precision::Float64=flux_surfaces_precision,
        raise_error_on_not_open::Bool,
        raise_error_on_not_closed::Bool,
        verbose::Bool=false
    )
"""
function find_psi_boundary(
    dimR::Union{AbstractVector{T1},AbstractRange{T1}},
    dimZ::Union{AbstractVector{T1},AbstractRange{T1}},
    PSI::Matrix{T2},
    psirange_init::AbstractVector,
    RA::T3,
    ZA::T3,
    fw_r::AbstractVector{T4},
    fw_z::AbstractVector{T5},
    r_cache::AbstractVector{<:Real}=T1[],
    z_cache::AbstractVector{<:Real}=T1[];
    PSI_interpolant=IMAS.ψ_interpolant(dimR, dimZ, PSI).PSI_interpolant,
    precision::Float64=flux_surfaces_precision,
    raise_error_on_not_open::Bool,
    raise_error_on_not_closed::Bool,
    verbose::Bool=false
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}

    if verbose
        @show psirange_init
        plot(; aspect_ratio=:equal)
        contour!(dimR, dimZ, transpose(PSI); color=:gray, clim=(min(psirange_init...), max(psirange_init...)))
        if !isempty(fw_r)
            plot!(fw_r, fw_z; color=:black, lw=2, label="")
        end
        display(plot!())
    end

    if isempty(r_cache)
        r_cache, z_cache = IMASutils.contour_cache(dimR, dimZ)
    end

    # innermost tentative flux surface (which should be closed!)
    psi_axis = PSI_interpolant(RA, ZA)
    pr, pz = IMASutils.contour_from_midplane!(r_cache, z_cache, PSI, dimR, dimZ, psirange_init[1], RA, ZA, psi_axis)
    if isempty(pr) || is_open_polygon(pr, pz)
        if raise_error_on_not_closed
            error("Flux surface at ψ=$(psirange_init[1]) is not closed; ψ=[$(psirange_init[1])...$(psirange_init[end])]")
        else
            return (last_closed=nothing, first_open=nothing)
        end
    end
    if verbose
        plot!(pr, pz; color=:blue, label="")
        display(plot!())
    end

    # outermost tentative flux surface (which should be open!)
    pr, pz = IMASutils.contour_from_midplane!(r_cache, z_cache, PSI, dimR, dimZ, psirange_init[end], RA, ZA, psi_axis)
    if is_closed_surface(pr, pz, fw_r, fw_z)
        if raise_error_on_not_open
            error("Flux surface at ψ=$(psirange_init[end]) is not open; ψ=[$(psirange_init[1])...$(psirange_init[end])]")
        else
            return (last_closed=nothing, first_open=nothing)
        end
    end
    if verbose
        plot!(pr, pz; color=:red, label="")
        display(plot!())
    end

    psibeg = psirange_init[1]
    psiend = psirange_init[end]
    for k in 1:100
        psimid = (psibeg + psiend) / 2.0
        pr, pz = IMASutils.contour_from_midplane!(r_cache, z_cache, PSI, dimR, dimZ, psimid, RA, ZA, psi_axis)
        # closed flux surface
        if is_closed_surface(pr, pz, fw_r, fw_z)
            if verbose
                display(plot!(pr, pz; label="", color=:green))
            end
            psibeg = psimid
            if (abs(psiend - psibeg) / (abs(psiend + psibeg) / 2.0)) < precision
                return (last_closed=psimid, first_open=psiend)
            end
            # open flux surface
        else
            psiend = psimid
        end
    end

    return error("Could not find closed boundary between ψ=$(psirange_init[1]) and ψ=$(psirange_init[end])")
end

@compat public find_psi_boundary
push!(document[Symbol("Physics flux-surfaces")], :find_psi_boundary)

function is_closed_surface(pr::AbstractVector{T1}, pz::AbstractVector{T1}, fw_r::AbstractVector{T2}=T1[], fw_z::AbstractVector{T2}=T1[]) where {T1<:Real, T2<:Real}
    @assert length(pr) == length(pz)
    closed = !isempty(pr) && is_closed_polygon(pr, pz)
    if closed
        @assert length(fw_r) == length(fw_z)
        if !isempty(fw_r)
            closed = !intersects(pr, pz, fw_r, fw_z)
        end
    end
    return closed
end

"""
    find_psi_separatrix(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=flux_surfaces_precision) where {T<:Real}

Returns psi of the first magentic separatrix

Note: The first separatrix is the LCFS only in diverted plasmas
"""
function find_psi_separatrix(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=flux_surfaces_precision) where {T<:Real}
    psi_up = find_psi_2nd_separatrix(eqt).diverted
    psi_low = eqt.profiles_1d.psi[1]

    psi = (psi_up + psi_low) / 2.0
    err = Inf
    counter_max = 50

    for k in 1:counter_max
        surface = flux_surface(eqt, psi, :encircling, Float64[], Float64[])
        if length(surface) > 1
            error("find_psi_separatrix: more than one :encircling flux surfaces")
        end

        if abs(surface[1].r[1] - surface[1].r[end]) == 0 && abs(surface[1].z[1] - surface[1].z[end]) == 0
            # the surface is inside the separatrix (= "closed" with no wall)
            psi_low = psi
        else
            # the surface is outside the first separatrix
            psi_up = psi
        end

        err = abs(psi_up - psi_low) / abs(psi_low)
        psi = (psi_up + psi_low) / 2.0
        if abs(err) < precision
            break
        end
    end

    return (psi_sep_closed=psi_up, psi_sep_open=psi_low)
end

@compat public find_psi_separatrix
push!(document[Symbol("Physics flux-surfaces")], :find_psi_separatrix)

"""
    find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=flux_surfaces_precision) where {T<:Real}

Returns psi of the second magentic separatrix

returns `diverted` and `not_diverted` surface in a NamedTuple

  - `diverted` is when the surface starts and finishes on the same side of the midplane

  - `not_diverted` is when the surface starts and finishes on opposite sides of the midplane
"""
function find_psi_2nd_separatrix(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=flux_surfaces_precision) where {T<:Real}

    psi_separatrix = eqt.profiles_1d.psi[end]
    surface = flux_surface(eqt, psi_separatrix, :open, Float64[], Float64[]) # last closed flux surface

    # First check if we are in a double null configuration
    ZA = eqt.global_quantities.magnetic_axis.z
    # retrieve b = elongation * minor radius
    b = eqt.boundary.elongation * eqt.boundary.minor_radius

    for (r, z) in surface
        # exclude empty vectors, private region and surfaces starting and finishing in the same z range as the plasma
        if isempty(r) || all(z .> ZA) || all(z .< ZA) || all(z .> (ZA - b) .&& z .< (ZA + b))
            continue
        end
        if (z[end] - ZA) * (z[1] - ZA) < 0
            # if perfect double null, all open surfaces in the SOL start and finish in opposite sides of the midplane
            psi_axis = eqt.profiles_1d.psi[1] # psi value on axis
            psi_sign = sign(psi_separatrix - psi_axis) # +1 for increasing psi / -1 for decreasing psi
            if psi_sign > 0
                #increasing psi
                return (diverted=psi_separatrix, not_diverted=psi_separatrix * (1 + flux_surfaces_precision))
            else
                #decreasing psi
                return (diverted=psi_separatrix, not_diverted=psi_separatrix * (1 - flux_surfaces_precision))
            end
        end
    end

    # Single null case
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    psi_axis = eqt.profiles_1d.psi[1] # psi value on axis
    psi_sign = sign(psi_separatrix - psi_axis) # +1 for increasing psi / -1 for decreasing psi
    psi_low = psi_separatrix # this is closed
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
        surface = flux_surface(eqt, psi, :open, Float64[], Float64[])
        for (r, z) in surface
            if isempty(r) || all(z .> ZA) || all(z .< ZA) || all(z .> (ZA - b) .&& z .< (ZA + b))
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

    return (diverted=psi_low, not_diverted=psi_up)
end

@compat public find_psi_2nd_separatrix
push!(document[Symbol("Physics flux-surfaces")], :find_psi_2nd_separatrix)

"""
    find_psi_last_diverted(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation;
        precision::Float64=flux_surfaces_precision)

Returns `psi_last_lfs, `psi_first_lfs_far`, and `null_within_wall`

  - `psi_first_lfs_far` will be the first surface inside `OFL[:lfs_far]`

  - `psi_last_lfs` will be the last surface inside `OFL[:lfs]`

Precision between the two is defined on the poloidal crossection area at the OMP (`Psol*precision` = power flowing between `psi_first_lfs_far` and `psi_last_lfs` ~ 0)
"""
function find_psi_last_diverted(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation;
    precision::Float64=flux_surfaces_precision)

    # if no wall in dd, psi_last diverted not defined
    if isempty(wall_r) || isempty(wall_z) || isempty(eqt.boundary.x_point)
        return (psi_last_lfs=NaN, psi_first_open=NaN, psi_first_lfs_far=NaN, null_within_wall=true)
    end

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    psi_boundaries = (last_closed=eqt.boundary.psi, first_open=eqt.boundary_separatrix.psi)
    if psi_boundaries.last_closed == psi_boundaries.first_open
        # the plasma is limited and not diverted
        limited = true
        Xpoint2 = [eqt.boundary.x_point[1].r, eqt.boundary.x_point[1].z]
        psi_2ndseparatrix = psi_sep_open # psi first magnetic separatrix
    else
        limited = false
        Xpoint2 = [eqt.boundary.x_point[end].r, eqt.boundary.x_point[end].z]
        psi_2ndseparatrix = find_psi_2nd_separatrix(eqt).not_diverted # psi second magnetic separatrix
    end

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
    @assert r_max > 0.0 "Could not trace an open flux surface for 2nd separatrix"

    r_intersect = Float64[rr for (rr, zz) in rz_intersects]
    z_intersect = Float64[zz for (rr, zz) in rz_intersects]

    # r_mid(ψ) interpolator for region of interest
    r_mid_of_interest = 10.0 .^ range(log10(maximum(eqt.boundary.outline.r) * 0.99), log10(r_max * 1.1), 1000)
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
        find_strike_points!(eqt, wall_r, wall_z, psi_boundaries.last_closed, psi_boundaries.first_open)
    end

    # find the two surfaces `psi_first_lfs_far` and `psi_last_lfs` around the last diverted flux surface
    psi_2ndseparatrix_notdiverted = psi_2ndseparatrix
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt).diverted
    if limited
        psi_first_lfs_far = psi_boundaries.last_closed
    else
        psi_first_lfs_far = psi_2ndseparatrix
    end

    psi_last_lfs = psi_boundaries.first_open
    psi = (psi_first_lfs_far + psi_last_lfs) / 2

    counter_max = 100
    err = Inf
    for counter in 1:counter_max
        surface = flux_surface(eqt, psi, :any, Float64[], Float64[])

        for (r, z) in surface
            if isempty(r) || all(z .> ZA) || all(z .< ZA)
                continue
            end

            crossings = intersection(r, z, wall_r, wall_z).crossings # find where flux surface crosses wall ("strike points" of surface)

            if isempty(crossings)
                continue
            end
            l = length(crossings) # save number of intersections before filtering

            # r and z coordiante of intersections with wall
            r_intersect = (cr[1] for cr in crossings)
            z_intersect = (cr[2] for cr in crossings)

            # find which intersection is the closest to each strike-point, to remove them
            closest_to_strike_points = Int64[]
            for point in eqt.boundary.strike_point
                dist = (r_intersect .- point.r) .^ 2 .+ (z_intersect .- point.z) .^ 2
                push!(closest_to_strike_points, argmin(dist))
            end
            # order closest intersections with the same oder as crossings ([1] is closest to [1] of surface; [end] is closest to [end] of surface)
            sort!(closest_to_strike_points)

            # discard all points in crossings before/after the closest intersections,
            # only if they are not crossings[1] and crossings[end] which by construction are the closest to the ends of the magnetic surface
            if !isempty(closest_to_strike_points)
                # discard crossings occuring after second strike point (clockwise)
                if closest_to_strike_points[end] < length(crossings)
                    for k in reverse(closest_to_strike_points[end]+1:length(crossings))
                        deleteat!(crossings, k)
                    end
                end
                # discard crossings occuring before first strike point (clockwise)
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
                if l == 2
                    # if orginally psi had 2 intersections, move up psi_low instead of lowering psi_up
                    # needed for limited case, probably cannot be done before filtering intersections (check)
                    psi_last_lfs = psi
                else
                    psi_first_lfs_far = psi
                end
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

    # if LDFS is the 1st separatrix, be consistent with the find_psi_separatrix function
    if psi_first_lfs_far == psi_boundaries.last_closed
        psi_last_lfs = psi_boundaries.last_closed
        psi_first_lfs_far = psi_boundaries.first_open
    end

    return (psi_last_lfs=psi_last_lfs, psi_first_open=psi_boundaries.first_open, psi_first_lfs_far=psi_first_lfs_far, null_within_wall=null_within_wall)
end

@compat public find_psi_last_diverted
push!(document[Symbol("Physics flux-surfaces")], :find_psi_last_diverted)

"""
    find_psi_tangent_omp(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation;
        precision::Float64=flux_surfaces_precision)

Returns the psi of the magnetic surface in the SOL which is tangent to the wall near the outer midplane
"""
function find_psi_tangent_omp(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation;
    precision::Float64=flux_surfaces_precision)

    psi_max = find_psi_max(eqt)

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
        surface = flux_surface(eqt, psi, :open, Float64[], Float64[])

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
    # if tangent flux surface is outside the maximum useful value (psi_max), use psi_max and throw a warning
    if psi_sign > 0
        # psi is increasing
        psi_tangent_in = minimum([psi_max, psi_tangent_in])
        psi_tangent_out = minimum([psi_max, psi_tangent_out])
    else
        # psi is decreasing
        psi_tangent_in = maximum([psi_max, psi_tangent_in])
        psi_tangent_out = maximum([psi_max, psi_tangent_out])
    end
    if psi_tangent_in == psi_max || psi_tangent_out == psi_max
        @warn "The tangent flux surface to the omp is outside the last useful surface (psi_max) [point of contact: G.Dose]"
    end
    return (psi_tangent_in=psi_tangent_in, psi_tangent_out=psi_tangent_out)
end

@compat public find_psi_tangent_omp
push!(document[Symbol("Physics flux-surfaces")], :find_psi_tangent_omp)

"""
    find_psi_max(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=1e-2) where {T<:Real}

Returns the max psi useful for an ofl in the SOL with no wall.
"""
function find_psi_max(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=1e-2) where {T<:Real}

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt).not_diverted
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
    surface = flux_surface(eqt, psi_up, :open, Float64[], Float64[])
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
        surface = flux_surface(eqt, psi, :open, Float64[], Float64[])
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

@compat public find_psi_max
push!(document[Symbol("Physics flux-surfaces")], :find_psi_max)

"""
    find_psi_wall_omp(eqt::IMAS.equilibrium__time_slice, wall_r::AbstractVector{<:Real}, wall_z::AbstractVector{<:Real})

Returns the psi of the magnetic surface in the SOL which intersects the wall at the outer midplane
"""
function find_psi_wall_omp(eqt::IMAS.equilibrium__time_slice, wall_r::AbstractVector{<:Real}, wall_z::AbstractVector{<:Real})
    @assert length(wall_r) == length(wall_z)

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    psi_max = find_psi_max(eqt)
    psi_sign = sign(psi_max - eqt.profiles_1d.psi[1])

    return find_psi_wall_omp(PSI_interpolant, RA, ZA, wall_r, wall_z, psi_max, psi_sign)
end

"""
    find_psi_wall_omp(
        PSI_interpolant::Interpolations.AbstractInterpolation,
        RA::T1,
        ZA::T1,
        wall_r::AbstractVector{T2},
        wall_z::AbstractVector{T2},
        psi_max::T1,
        psi_sign::T1
    ) where {T1<:Real,T2<:Real}
"""
function find_psi_wall_omp(
    PSI_interpolant::Interpolations.AbstractInterpolation,
    RA::T1,
    ZA::T1,
    wall_r::AbstractVector{T2},
    wall_z::AbstractVector{T2},
    psi_max::T1,
    psi_sign::T1
) where {T1<:Real,T2<:Real}
    crossings = intersection([RA, maximum(wall_r) * 1.1], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
    # if flux surface passing through the wall at the omp is outside the maximum useful value (psi_max), use psi_max and throw a warning
    if psi_sign > 0
        # psi is increasing
        psi_wall_omp = minimum([psi_max, PSI_interpolant.(r_wall_midplane, ZA)[1]])
    else
        # psi is decreasing
        psi_wall_omp = maximum([psi_max, PSI_interpolant.(r_wall_midplane, ZA)[1]])
    end
    if psi_wall_omp == psi_max
        @warn "The flux surface passing through the wall at the omp is outside the last useful surface (psi_max) [point of contact: G.Dose]"
    end
    return psi_wall_omp
end

@compat public find_psi_wall_omp
push!(document[Symbol("Physics flux-surfaces")], :find_psi_wall_omp)

"""
    interp_rmid_at_psi(eqt::IMAS.equilibrium__time_slice, PSI_interpolant::Interpolations.AbstractInterpolation, R::AbstractVector{<:Real})

Returns the interpolant r_mid(ψ) to compute the r at the midplane of the flux surface identified by ψ

The vector `R` defines the sampling of interest for thie interpolation
"""
function interp_rmid_at_psi(PSI_interpolant::Interpolations.AbstractInterpolation, R::AbstractVector{T}, ZA::T) where {T<:Real}
    return cubic_interp1d(PSI_interpolant.(R, R .* 0.0 .+ ZA), R)
end

@compat public interp_rmid_at_psi
push!(document[Symbol("Physics flux-surfaces")], :interp_rmid_at_psi)

function find_magnetic_axis(
    r::AbstractVector{<:T1},
    z::AbstractVector{<:T1},
    PSI_interpolant::Interpolations.AbstractInterpolation,
    psi_sign::T2;
    rguess::T1=r[round(Int, length(r) / 2)],
    zguess::T1=z[round(Int, length(z) / 2)]) where {T1<:Real,T2<:Real}

    res = Optim.optimize(
        x -> begin
            try
                PSI_interpolant(x[1], x[2]) * psi_sign
            catch e
                if typeof(e) <: BoundsError
                    return T2(Inf)
                else
                    rethrow(e)
                end
            end
        end,
        [rguess, zguess],
        Optim.Newton()#;
        #autodiff=:forward
    )

    return (RA=res.minimizer[1], ZA=res.minimizer[2])
end


abstract type AbstractFluxSurface{T<:Real} end

"""
    struct SimpleSurface{T<:Real} <: AbstractFluxSurface{T}
        psi::T
        r::Vector{T}
        z::Vector{T}
        ll::Vector{T}
        fluxexpansion::Vector{T}
        int_fluxexpansion_dl::T
    end

A simplified version of FluxSurface that only has the contour points and
what is needed to compute flux surface averages
"""
mutable struct SimpleSurface{T<:Real} <: AbstractFluxSurface{T}
    psi::T
    r::Vector{T}
    z::Vector{T}
    ll::Vector{T}
    fluxexpansion::Vector{T}
    int_fluxexpansion_dl::T
end

@compat public SimpleSurface
push!(document[Symbol("Physics flux-surfaces")], :SimpleSurface)

@recipe function plot_SimpleSurface(surface::SimpleSurface)
    @series begin
        aspect_ratio := :equal
        surface.r, surface.z
    end
end

"""
    struct FluxSurface{T<:Real} <: AbstractFluxSurface{T}
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
"""
mutable struct FluxSurface{T<:Real} <: AbstractFluxSurface{T}
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

@compat public FluxSurface
push!(document[Symbol("Physics flux-surfaces")], :FluxSurface)

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

@recipe function plot_FluxSurfaces(surfaces::Vector{<:AbstractFluxSurface})
    for k in eachindex(surfaces)
        @series begin
            label --> ""
            primary := k == 1
            surfaces[k]
        end
    end
end

"""
    trace_simple_surfaces(
        psi::AbstractVector{T},
        r::AbstractVector{T},
        z::AbstractVector{T},
        PSI::Matrix{T},
        PSI_interpolant::Interpolations.AbstractInterpolation,
        RA::T,
        ZA::T,
        wall_r::AbstractVector{T},
        wall_z::AbstractVector{T}
    ) where {T<:Real}

Trace flux surfaces and returns vector of SimpleSurface structures. The result
contains only the contours and what is needed to perform flux-surface averaging.
"""
function trace_simple_surfaces(
    psi::AbstractVector{T},
    r::AbstractVector{T},
    z::AbstractVector{T},
    PSI::Matrix{T},
    PSI_interpolant::Interpolations.AbstractInterpolation,
    RA::T,
    ZA::T,
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T}) where {T<:Real}

    surfaces = Vector{SimpleSurface{T}}(undef, length(psi))
    r_cache, z_cache = IMASutils.contour_cache(r, z)
    return trace_simple_surfaces!(surfaces, psi, r, z, PSI, PSI_interpolant, RA, ZA, wall_r, wall_z, r_cache, z_cache)
end

@compat public trace_simple_surfaces
push!(document[Symbol("Physics flux-surfaces")], :trace_simple_surfaces)

"""
    trace_simple_surfaces!(
        surfaces::Vector{SimpleSurface{T}},
        psi::AbstractVector{T},
        r::AbstractVector{T},
        z::AbstractVector{T},
        PSI::Matrix{T},
        PSI_interpolant::Interpolations.AbstractInterpolation,
        RA::T,
        ZA::T,
        wall_r::AbstractVector{T},
        wall_z::AbstractVector{T},
        r_cache::AbstractVector{T}=T[],
        z_cache::AbstractVector{T}=T[])
    ) where {T<:Real}

Trace flux surfaces and store in `surfaces` vector of SimpleSurface structures.
The result contains only the contours and what is needed to perform flux-surface averaging.
"""
function trace_simple_surfaces!(
    surfaces::Vector{SimpleSurface{T}},
    psi::AbstractVector{<:Real},
    r::AbstractVector{<:Real},
    z::AbstractVector{<:Real},
    PSI::Matrix{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation,
    RA::Real,
    ZA::Real,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    r_cache::AbstractVector{T}=T[],
    z_cache::AbstractVector{T}=T[]) where {T<:Real}

    N = length(psi)

    if isempty(r_cache) || isempty(z_cache)
        r_cache, z_cache = IMASutils.contour_cache(r, z)
    end
    PSIA = PSI_interpolant(RA, ZA)

    for k in N:-1:1
        psi_level = psi[k]

        # trace surfaces
        if k == 1 # on axis flux surface is a artificial one, generated from the second surface
            pr = (surfaces[2].r .- RA) ./ 100.0 .+ RA
            pz = (surfaces[2].z .- ZA) ./ 100.0 .+ ZA

        else  # other flux surfaces
            # trace flux surface
            pr, pz = IMASutils.contour_from_midplane!(r_cache, z_cache, PSI, r, z, psi_level, RA, ZA, PSIA)

            if !is_closed_surface(pr, pz, wall_r, wall_z)
                error("IMAS: Could not trace closed flux surface $k out of $(N) at ψ = $(psi_level)")
            end
        end

        ll = zero(pr)
        fluxexpansion = similar(pr)
        for i in eachindex(pr)
            if i > 1
                # surface length
                dl = sqrt((pr[i] - pr[i-1])^2 + (pz[i] - pz[i-1])^2)
                ll[i] = ll[i-1] + dl
            end

            # flux expansion = 1 / abs(Bp)
            Br, Bz = Br_Bz(PSI_interpolant, pr[i], pz[i])
            fluxexpansion[i] = 1.0 / sqrt(Br^2.0 + Bz^2.0)
        end
        int_fluxexpansion_dl = trapz(ll, fluxexpansion)

        # create
        surfaces[k] = SimpleSurface(psi_level, collect(pr), collect(pz), ll, fluxexpansion, int_fluxexpansion_dl)
    end

    return surfaces
end

@compat public trace_simple_surfaces!
push!(document[Symbol("Physics flux-surfaces")], :trace_simple_surfaces!)

"""
    trace_surfaces(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}; refine_extrema::Bool=true) where {T<:Real}

Trace flux surfaces and returns vector of FluxSurface structures
"""
function trace_surfaces(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}; refine_extrema::Bool=true) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z
    Br, Bz = Br_Bz(eqt2d)
    return trace_surfaces(eqt.profiles_1d.psi, eqt.profiles_1d.f, r, z, eqt2d.psi, Br, Bz, PSI_interpolant, RA, ZA, wall_r, wall_z; refine_extrema)
end

"""
    trace_surfaces(
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
        wall_z::AbstractVector{T};
        refine_extrema::Bool=true
    ) where {T<:Real}
"""
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
    wall_z::AbstractVector{T};
    refine_extrema::Bool=true
) where {T<:Real}

    N = length(psi)
    surfaces = Vector{FluxSurface{T}}(undef, N)
    r_cache, z_cache = IMASutils.contour_cache(r, z)
    PSIA = PSI_interpolant(RA, ZA)
    for k in N:-1:1
        psi_level = psi[k]

        # trace surfaces
        if k == 1 # on axis flux surface is a artificial one, generated from the second surface
            pr = (surfaces[2].r .- RA) ./ 100.0 .+ RA
            pz = (surfaces[2].z .- ZA) ./ 100.0 .+ ZA

        else  # other flux surfaces
            tmp = IMASutils.contour_from_midplane!(r_cache, z_cache, PSI, r, z, psi_level, RA, ZA, PSIA)
            pr, pz = collect(tmp[1]), collect(tmp[2])
            if isempty(pr) && k == 2
                # If there are too many flux surfaces for a given grid resolution, the inner-most flux surface may fail to trace
                pr = (surfaces[k+1].r .- RA) ./ 2 .+ RA
                pz = (surfaces[k+1].z .- ZA) ./ 2 .+ ZA
            end
            if k == N && !is_closed_surface(pr, pz, wall_r, wall_z) || isempty(pr)
                # contour(r, z, PSI'; levels=psi)
                # display(plot!(wall_r, wall_z; aspect_ratio=:equal, color=:black))
                error("IMAS: Could not trace closed flux surface $k out of $(N) at ψ = $(psi_level)")
            end
        end

        # surface length
        ll = arc_length(pr, pz)

        # poloidal magnetic field (with sign)
        Br, Bz = Br_Bz(PSI_interpolant, pr, pz)
        Bp2 = Br .^ 2.0 .+ Bz .^ 2.0
        Bp_abs = sqrt.(Bp2)
        Bp = Bp_abs .* sign.((pz .- ZA) .* Br .- (pr .- RA) .* Bz)
        Btot = sqrt.(Bp2 .+ (f[k] ./ pr) .^ 2)

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

    if refine_extrema
        N2 = round(Int, N / 2, RoundUp)
        psi_norm = abs(psi[end] - psi[1]) / N
        space_norm = (surfaces[N].max_r - surfaces[N].min_r) / 2 / N

        # Find where Br changes sign
        lines = Contour.lines(Contour.contour(r, z, BR, 0.0))
        k = 0
        d = Inf
        for (kk, line) in enumerate(lines)
            pr, pz = Contour.coordinates(line)
            # plot!(pr, pz)
            dd = minimum(filter(!isnan, sqrt.((pr .- surfaces[N2].max_r) .^ 2 .+ (pz .- surfaces[N2].z_at_max_r) .^ 2)))
            if dd < d
                d = dd
                k = kk
            end
        end
        leftright_r, leftright_z = Contour.coordinates(lines[k])
        index = .!(isnan.(leftright_r) .|| isnan.(leftright_z))
        leftright_r = @view leftright_r[index]
        leftright_z = @view leftright_z[index]

        # extrema in R
        interp_r = interp1d(1:length(leftright_r), leftright_r)
        interp_z = interp1d(1:length(leftright_z), leftright_z)
        for k in 1:N
            #@show "R", k
            index = _extrema_index(leftright_r, leftright_z, surfaces[k].max_r, surfaces[k].z_at_max_r, :right)
            cost =
                x -> _extrema_cost(
                    interp_r(x),
                    interp_z(x),
                    psi[k],
                    PSI_interpolant,
                    surfaces[k].max_r,
                    surfaces[k].z_at_max_r,
                    RA,
                    ZA,
                    psi_norm,
                    space_norm,
                    :right
                )
            x = Optim.optimize(cost, index[1], index[end], Optim.Brent()).minimizer
            surfaces[k].max_r, surfaces[k].z_at_max_r = interp_r(x), interp_z(x)
            index = _extrema_index(leftright_r, leftright_z, surfaces[k].min_r, surfaces[k].z_at_min_r, :left)
            cost =
                x -> _extrema_cost(
                    interp_r(x),
                    interp_z(x),
                    psi[k],
                    PSI_interpolant,
                    surfaces[k].min_r,
                    surfaces[k].z_at_min_r,
                    RA,
                    ZA,
                    psi_norm,
                    space_norm,
                    :left
                )
            #plot!(leftright_r[index], leftright_z[index]; label="")
            x = Optim.optimize(cost, index[1], index[end], Optim.Brent()).minimizer
            min_r, z_at_min_r = interp_r(x), interp_z(x)
            surfaces[k].min_r, surfaces[k].z_at_min_r = min_r, z_at_min_r
            @assert surfaces[k].min_r < surfaces[k].max_r
            if k < 3
                surfaces[k+1].z_at_max_r = ZA
                surfaces[k+1].z_at_min_r = ZA
            elseif k < N
                surfaces[k+1].z_at_max_r = (surfaces[k+1].z_at_max_r + surfaces[k].z_at_max_r) / 2.0
                surfaces[k+1].z_at_min_r = (surfaces[k].z_at_min_r + surfaces[k+1].z_at_min_r) / 2.0
            end
        end

        # Find where Bz changes sign
        lines = Contour.lines(Contour.contour(r, z, BZ, 0.0))
        k = 0
        d = Inf
        for (kk, line) in enumerate(lines)
            pr, pz = Contour.coordinates(line)
            # plot!(pr, pz)
            dd = minimum(filter(!isnan, sqrt.((pr .- surfaces[N2].r_at_max_z) .^ 2 .+ (pz .- surfaces[N2].max_z) .^ 2)))
            if dd < d
                d = dd
                k = kk
            end
        end
        updown_r, updown_z = Contour.coordinates(lines[k])
        index = .!(isnan.(updown_r) .|| isnan.(updown_z))
        updown_r = @view updown_r[index]
        updown_z = @view updown_z[index]

        # extrema in Z
        interp_r = interp1d(1:length(updown_r), updown_r)
        interp_z = interp1d(1:length(updown_z), updown_z)
        for k in 1:N
            #@show "Z", k
            index = _extrema_index(updown_r, updown_z, surfaces[k].r_at_max_z, surfaces[k].max_z, :up)
            cost =
                x -> _extrema_cost(
                    interp_r(x),
                    interp_z(x),
                    psi[k],
                    PSI_interpolant,
                    surfaces[k].r_at_max_z,
                    surfaces[k].max_z,
                    RA,
                    ZA,
                    psi_norm,
                    space_norm,
                    :up
                )
            x = Optim.optimize(cost, index[1], index[end], Optim.Brent()).minimizer
            surfaces[k].r_at_max_z, surfaces[k].max_z = interp_r(x), interp_z(x)
            index = _extrema_index(updown_r, updown_z, surfaces[k].r_at_min_z, surfaces[k].min_z, :down)
            cost =
                x -> _extrema_cost(
                    interp_r(x),
                    interp_z(x),
                    psi[k],
                    PSI_interpolant,
                    surfaces[k].r_at_min_z,
                    surfaces[k].min_z,
                    RA,
                    ZA,
                    psi_norm,
                    space_norm,
                    :down
                )
            # plot!(updown_r[index], updown_z[index]; label="")
            x = Optim.optimize(cost, index[1], index[end], Optim.Brent()).minimizer
            surfaces[k].r_at_min_z, surfaces[k].min_z = interp_r(x), interp_z(x)
            @assert surfaces[k].min_z < surfaces[k].max_z
            if k < 3
                surfaces[k+1].r_at_max_z = RA
                surfaces[k+1].r_at_min_z = RA
            elseif k < N
                surfaces[k+1].r_at_max_z = (surfaces[k+1].r_at_max_z + surfaces[k].r_at_max_z) / 2.0
                surfaces[k+1].r_at_min_z = (surfaces[k+1].r_at_min_z + surfaces[k].r_at_min_z) / 2.0
            end
        end

        # first flux surface just a scaled down version of the second one
        frac = 0.01
        surfaces[1].r_at_max_z = (surfaces[2].r_at_max_z .- RA) .* frac .+ RA
        surfaces[1].max_z = (surfaces[2].max_z .- ZA) .* frac .+ ZA
        surfaces[1].r_at_min_z = (surfaces[2].r_at_min_z .- RA) .* frac .+ RA
        surfaces[1].min_z = (surfaces[2].min_z .- ZA) .* frac .+ ZA
        surfaces[1].z_at_max_r = (surfaces[2].z_at_max_r .- ZA) .* frac .+ ZA
        surfaces[1].max_r = (surfaces[2].max_r - RA) * frac + RA
        surfaces[1].z_at_min_r = (surfaces[2].z_at_min_r .- ZA) .* frac .+ ZA
        surfaces[1].min_r = (surfaces[2].min_r - RA) * frac + RA
    end

    return FluxSurface[surfaces[k] for k in 1:N]
end

@compat public trace_surfaces
push!(document[Symbol("Physics flux-surfaces")], :trace_surfaces)

function _extrema_index(r::AbstractVector{T}, z::AbstractVector{T}, r0::T, Z0::T, direction::Symbol) where {T<:Real}
    i = argmin((r .- r0) .^ 2 .+ (z .- Z0) .^ 2)
    N = length(z)
    j = round(Int, N / 2, RoundUp)
    n = 3
    if direction == :right
        if (r[j] - r[j-1]) > 0 # oriented right
            return max(1, i - 1):min(N, i + n)
        else # opposite orientation
            return max(1, i - n):min(N, i + 1)
        end
    elseif direction == :left
        if (r[j] - r[j-1]) < 0 # oriented left
            return max(1, i - 1):min(N, i + n)
        else
            return max(1, i - n):min(N, i + 1)
        end
    elseif direction == :up
        if (z[j] - z[j-1]) > 0 # oriented up
            return max(1, i - 1):min(N, i + n)
        else
            return max(1, i - n):min(N, i + 1)
        end
    elseif direction == :down
        if (z[j] - z[j-1]) < 0 # oriented down
            return max(1, i - 1):min(N, i + n)
        else
            return max(1, i - n):min(N, i + 1)
        end
    else
        error("_extrema_index(..., direction::Symbol) can only be (:left, :right, :up, :down)")
    end
end

# accurate geometric quantities by finding geometric extrema as optimization problem
function _extrema_cost(
    r::T,
    z::T,
    psi_level::T,
    PSI_interpolant,
    r_orig::T,
    z_orig::T,
    RA::T,
    ZA::T,
    psi_norm::T,
    space_norm::T,
    direction::Symbol
) where {T<:Real}
    cost_psi = (PSI_interpolant(r, z) - psi_level) / psi_norm * space_norm # convert psi cost into spatial units
    if direction == :right
        cost_dir = abs(r - r_orig) * (r > r_orig)
        cost_dir += abs(r - RA) * (r < RA) + (r < RA) # extra penalty needed for flux surfaces near axis
    elseif direction == :left
        cost_dir = abs(r - r_orig) * (r < r_orig)
        cost_dir += abs(r - RA) * (r > RA) + (r > RA)
    elseif direction == :up
        cost_dir = abs(z - z_orig) * (z > z_orig)
        cost_dir += abs(z - ZA) * (z < ZA) + (z < ZA)
    elseif direction == :down
        cost_dir = abs(z - z_orig) * (z < z_orig)
        cost_dir += abs(z - ZA) * (z > ZA) + (z > ZA)
    else
        error("_extrema_cost(..., direction::Symbol) can only be (:left, :right, :up, :down)")
    end
    cost = norm((cost_psi, cost_dir^2)) # linear in psi since it already varies quadratically in space
    return cost
end

"""
    flux_surfaces(eq::equilibrium{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Update flux surface averaged and geometric quantities for all time_slices in the equilibrium IDS
"""
function flux_surfaces(eq::equilibrium{T1}, wall_r::AbstractVector{T2}, wall_z::AbstractVector{T2}) where {T1<:Real,T2<:Real}
    Threads.@threads for time_index in eachindex(eq.time_slice)
        eqt = eq.time_slice[time_index]
        eqt2d = findfirst(:rectangular, eqt.profiles_2d)
        if eqt2d !== nothing
            flux_surfaces(eqt, wall_r, wall_z)
        end
    end
    return eq
end

"""
    flux_surfaces(eqt::equilibrium__time_slice{T1}, wall_r::AbstractVector{T2}, wall_z::AbstractVector{T2}) where {T1<:Real,T2<:Real}

Update flux surface averaged and geometric quantities for a given equilibrum IDS time slice.
"""
function flux_surfaces(eqt::equilibrium__time_slice{T1}, wall_r::AbstractVector{T2}, wall_z::AbstractVector{T2}) where {T1<:Real,T2<:Real}
    eqt1d = eqt.profiles_1d
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)
    PSI = eqt2d.psi

    psi_sign = sign(eqt1d.psi[end] - eqt1d.psi[1])
    R0 = eqt.global_quantities.vacuum_toroidal_field.r0
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0

    # ensure certain global quantities are consistent with 1d profiles by making them expressions
    for field in (:psi_boundary, :psi_axis, :q_95, :q_axis, :q_min)
        empty!(eqt.global_quantities, field)
    end

    # accurately find magnetic axis
    # If there's a wall, try the following:
    # 1. guess from the previous time slice
    # 2. guess by  finding all interior extrema that could be near axes then exclude points outside the first wall (i.e., likely coils)
    rguess = r[round(Int, length(r) / 2)]
    zguess = z[round(Int, length(z) / 2)]
    eqts = parent(eqt)
    if eqts !== nothing && index(eqts, eqt) > 1 && hasdata(eqts[index(eqts, eqt)-1].global_quantities.magnetic_axis, :r)
        rguess = eqts[index(eqts, eqt)-1].global_quantities.magnetic_axis.r
        zguess = eqts[index(eqts, eqt)-1].global_quantities.magnetic_axis.z
    elseif !isempty(wall_r)
        pts = (psi_sign > 0.0) ? IMASutils.findall_interior_argmin(PSI) : IMASutils.findall_interior_argmax(PSI)
        min_wall_r = minimum(wall_r)
        max_wall_r = maximum(wall_r)
        min_wall_z = minimum(wall_z)
        max_wall_z = maximum(wall_z)
        pts = [pt for pt in pts if (r[pt[1]] > min_wall_r && r[pt[1]] < max_wall_r && z[pt[2]] > min_wall_z && z[pt[2]] < max_wall_z)]
        if length(pts) == 1
            rguess = r[pts[1][1]]
            zguess = z[pts[1][2]]
        else
            wall = collect(zip(wall_r, wall_z))
            pts_in_wall = [PolygonOps.inpolygon((r[i], z[j]), wall) == 1 for (i, j) in pts] # 1 if in wall
            if length(pts_in_wall) == 1
                ig, jg = pts[argmax(pts_in_wall)] # this are the indices for the only guess inside the wall
                rguess = r[ig]
                zguess = z[jg]
            end
        end
    end
    RA, ZA = find_magnetic_axis(r, z, PSI_interpolant, psi_sign; rguess, zguess)
    eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z = RA, ZA
    psi_axis = PSI_interpolant(RA, ZA)

    # accurately find the lcfs and scale psi accordingly
    if !isempty(wall_r) || eqt.global_quantities.free_boundary == 1
        psi_boundaries = find_psi_boundary(r, z, eqt2d.psi, psi_axis, eqt1d.psi[end], RA, ZA, wall_r, wall_z;
            PSI_interpolant, raise_error_on_not_open=false, raise_error_on_not_closed=false)
    else
        psi_boundaries = (first_open=eqt1d.psi[end], last_closed=eqt1d.psi[end])
    end
    if psi_boundaries.first_open !== psi_boundaries.last_closed
        eqt1d.psi =
            (eqt1d.psi .- eqt1d.psi[1]) ./ (eqt1d.psi[end] - eqt1d.psi[1]) .* (psi_boundaries.last_closed - psi_axis) .+ psi_axis
    end
    eqt.boundary.psi = psi_boundaries.last_closed
    eqt.boundary_separatrix.psi = psi_boundaries.first_open

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
        setproperty!(eqt1d, item, zeros(eltype(eqt1d.psi), size(eqt1d.psi)))
    end

    # trace flux surfaces
    Br, Bz = Br_Bz(eqt2d)
    surfaces = trace_surfaces(eqt1d.psi, eqt1d.f, r, z, eqt2d.psi, Br, Bz, PSI_interpolant, RA, ZA, wall_r, wall_z)

    # calculate flux surface averaged and geometric quantities
    N = length(eqt1d.psi)

    Np = maximum(length(surface.r) for surface in surfaces)
    shared_tmp = Vector{T1}(undef, Np)
    for k in N:-1:1
        surface = surfaces[k]
        pr = surface.r
        pz = surface.z
        Btot = surface.Btot
        Bp = surface.Bp

        tmp = @views shared_tmp[1:length(pr)]

        eqt1d.r_outboard[k] = surface.max_r
        eqt1d.r_inboard[k] = surface.min_r

        # miller geometric coefficients
        _, _, κ, δu, δl, ζou, ζol, ζil, ζiu =
            miller_R_a_κ_δ_ζ(pr, pz, surface.r_at_max_z, surface.max_z, surface.r_at_min_z, surface.min_z, surface.z_at_max_r, surface.max_r, surface.z_at_min_r, surface.min_r)
        eqt1d.elongation[k] = κ
        eqt1d.triangularity_upper[k] = δu
        eqt1d.triangularity_lower[k] = δl
        eqt1d.squareness_lower_outer[k] = ζol
        eqt1d.squareness_upper_outer[k] = ζou
        eqt1d.squareness_lower_inner[k] = ζil
        eqt1d.squareness_upper_inner[k] = ζiu

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
        eqt1d.trapped_fraction[k] = 0.75 * ftu + 0.25 * ftl

        # Bavg
        eqt1d.b_field_average[k] = avg_Btot

        # Bmax
        eqt1d.b_field_max[k] = Bmax

        # Bmin
        eqt1d.b_field_min[k] = Bmin

        # gm1 = <1/R^2>
        tmp .= 1.0 ./ pr .^ 2
        eqt1d.gm1[k] = flux_surface_avg(tmp, surface)

        # gm4 = <1/B^2>
        tmp .= 1.0 ./ Btot .^ 2
        eqt1d.gm4[k] = flux_surface_avg(tmp, surface)

        # gm5 = <B^2>
        eqt1d.gm5[k] = avg_Btot2

        # gm8 = <R>
        eqt1d.gm8[k] = flux_surface_avg(pr, surface)

        # gm9 = <1/R>
        tmp .= 1.0 ./ pr
        eqt1d.gm9[k] = flux_surface_avg(tmp, surface)

        # gm10 = <R^2>
        tmp .= pr .^ 2
        eqt1d.gm10[k] = flux_surface_avg(tmp, surface)

        # fsa_bp = <Bp>
        eqt1d.fsa_bp[k] = flux_surface_avg(Bp, surface)

        # j_tor = <j_tor/R> / <1/R> [A/m²]
        eqt1d.j_tor[k] =
            (
                -(eqt1d.dpressure_dpsi[k] + eqt1d.f_df_dpsi[k] * eqt1d.gm1[k] / mks.μ_0) *
                (2π)
            ) / eqt1d.gm9[k]

        # dvolume_dpsi
        eqt1d.dvolume_dpsi[k] = sign(eqt1d.fsa_bp[k]) * surface.int_fluxexpansion_dl

        # surface area
        eqt1d.surface[k] = 2π * trapz(surface.ll, pr)

        # q
        eqt1d.q[k] = eqt1d.dvolume_dpsi[k] * eqt1d.f[k] * eqt1d.gm1[k] / (2π)

        # quantities calculated on the last closed flux surface
        if k == N
            # perimeter
            eqt.global_quantities.length_pol = surface.ll[end]

            # boundary
            eqt.boundary.outline.r = pr
            eqt.boundary.outline.z = pz
        end
    end

    tmp = similar(eqt1d.psi)

    # area
    tmp .= eqt1d.dvolume_dpsi .* eqt1d.gm9 ./ 2π
    eqt1d.area = cumtrapz(eqt1d.psi, tmp)
    eqt.global_quantities.area = eqt1d.area[end]

    # volume
    eqt1d.volume = cumtrapz(eqt1d.psi, eqt1d.dvolume_dpsi)
    eqt.global_quantities.volume = eqt1d.volume[end]

    # phi
    tmp .= eqt1d.f .* eqt1d.gm1 / (2π)
    eqt1d.phi = cumtrapz(eqt1d.volume, tmp)

    # rho_tor_norm
    rho = sqrt.(abs.(eqt1d.phi ./ (π * B0)))
    rho_meters = rho[end]
    eqt1d.rho_tor = rho
    eqt1d.rho_tor_norm = rho ./ rho_meters

    # phi 2D
    tmp .= eqt1d.psi .* psi_sign
    phi_itp = interp1d(tmp, eqt1d.phi, :linear) # must be linear
    eqt2d.phi = phi_itp.(psi_sign .* eqt2d.psi)

    # ip
    eqt.global_quantities.ip = Ip(eqt)

    # gm2: <∇ρ²/R²>
    cumtrapz!(tmp, eqt1d.area, eqt1d.j_tor) # It(psi)
    eqt1d.gm2 = (mks.μ_0 * (2π)^2) .* tmp ./ (eqt1d.dvolume_dpsi .* (eqt1d.dpsi_drho_tor .^ 2))
    @views gm2_itp = cubic_interp1d(eqt1d.rho_tor_norm[2:5], eqt1d.gm2[2:5])
    eqt.profiles_1d.gm2[1] = gm2_itp(0.0) # extrapolate to axis due to zero / zero division

    # Geometric major and minor radii
    Rgeo = (eqt1d.r_outboard[end] + eqt1d.r_inboard[end]) / 2.0
    a = (eqt1d.r_outboard[end] - eqt1d.r_inboard[end]) / 2.0

    # vacuum magnetic field at the geometric center
    Btvac = B0 * R0 / Rgeo

    # average poloidal magnetic field
    Bpave = eqt.global_quantities.ip * mks.μ_0 / eqt.global_quantities.length_pol

    # li
    Bp2v = trapz(eqt1d.psi, T1[trapz(surface.ll, surface.Bp) for surface in surfaces])
    eqt.global_quantities.li_3 = 2.0 * Bp2v / Rgeo / (eqt.global_quantities.ip * mks.μ_0)^2

    # beta_tor
    avg_press = volume_integrate(eqt, eqt1d.pressure) / eqt1d.volume[end]
    eqt.global_quantities.beta_tor = abs(avg_press / (Btvac^2 / 2.0 / mks.μ_0))

    # beta_pol
    eqt.global_quantities.beta_pol = abs(avg_press / (Bpave^2 / 2.0 / mks.μ_0))

    # beta_normal
    ip = eqt.global_quantities.ip / 1e6
    eqt.global_quantities.beta_normal = eqt.global_quantities.beta_tor / abs(ip / a / Btvac) * 100

    # find quantities on separatrix
    find_x_point!(eqt, wall_r, wall_z)

    # secondary separatrix
    if length(eqt.boundary.x_point) > 1
        psi2nd = find_psi_2nd_separatrix(eqt).diverted
        pts = flux_surface(r, z, eqt2d.psi, RA, ZA, wall_r, wall_z, psi2nd, :encircling)
        if !isempty(pts)
            eqt.boundary_secondary_separatrix.outline.r = pts[1][1]
            eqt.boundary_secondary_separatrix.outline.z = pts[1][2]
        end
    end

    # find strike points
    find_strike_points!(eqt, wall_r, wall_z, psi_boundaries.last_closed, psi_boundaries.first_open)

    return eqt
end

"""
    flux_surfaces(eqt::equilibrium__time_slice, wall::IMAS.wall)
"""
function flux_surfaces(eqt::equilibrium__time_slice, wall::IMAS.wall)
    fw = first_wall(wall)
    return flux_surfaces(eqt, fw.r, fw.z)
end

@compat public flux_surfaces
push!(document[Symbol("Physics flux-surfaces")], :flux_surfaces)

"""
    flux_surface(eqt::equilibrium__time_slice{T}, psi_level::Real, type::Symbol, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Returns a vector with the (r,z) coordiates of flux surface at given psi_level

The `type` parameter:

  - `:any`: return all contours
  - `:closed`: all closed flux-surface that encircle the magnetic axis and do not cross the wall
  - `:open`: all open flux-surfaces (considerning open even closed flux surfaces that hit the first wall)
  - `:open_no_wall`: all open flux-surfaces independently of wall
  - `:encircling`: open flux-surfaces encircling the magnetic axis
"""
function flux_surface(eqt::equilibrium__time_slice{T}, psi_level::Real, type::Symbol, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    if eqt2d === nothing
        error("Equilibrium at $(eqt.time) [s] does not have a rectangular grid for tracing flux surfaces")
    end
    dim1 = to_range(eqt2d.grid.dim1)
    dim2 = to_range(eqt2d.grid.dim2)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z
    return flux_surface(dim1, dim2, eqt2d.psi, RA, ZA, wall_r, wall_z, psi_level, type)
end

"""
    flux_surface(
        dim1::AbstractVector{T1},
        dim2::AbstractVector{T1},
        PSI::AbstractArray{T2},
        RA::T3,
        ZA::T3,
        fw_r::AbstractVector{T4},
        fw_z::AbstractVector{T4},
        psi_level::T5,
        type::Symbol
    ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}
"""
function flux_surface(
    dim1::AbstractVector{T1},
    dim2::AbstractVector{T1},
    PSI::AbstractArray{T2},
    RA::T3,
    ZA::T3,
    fw_r::AbstractVector{T4},
    fw_z::AbstractVector{T4},
    psi_level::T5,
    type::Symbol
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}

    # contouring routine
    cl = Contour.contour(dim1, dim2, PSI, psi_level; VT=NTuple{2,T2})

    return flux_surface(dim1, dim2, cl, RA, ZA, fw_r, fw_z, type)
end

"""
    flux_surface(
        dim1::Union{AbstractVector{T},AbstractRange{T}},
        dim2::Union{AbstractVector{T},AbstractRange{T}},
        cl::Contour.ContourLevel,
        RA::T,
        ZA::T,
        fw_r::AbstractVector{T},
        fw_z::AbstractVector{T},
        type::Symbol) where {T<:Real}
"""
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
            if (is_closed_polygon(pr, pz) && (PolygonOps.inpolygon((RA, ZA), collect(zip(pr, pz))) == 1) && !intersects(pr, pz, fw_r, fw_z))
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
            if intersects(pr, pz, fw_r, fw_z)
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
        for (pr, pz) in flux_surface(dim1, dim2, cl, RA, ZA, fw_r, fw_z, :any)
            if all(pz .> ZA) || all(pz .< ZA)
                continue
            end
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

@compat public flux_surface
push!(document[Symbol("Physics flux-surfaces")], :flux_surface)

"""
    flux_surface_avg(quantity::AbstractVector{T}, surface::FluxSurface{T}) where {T<:Real}

Flux surface averaging of a quantity
"""
function flux_surface_avg(quantity::AbstractVector{T}, surface::AbstractFluxSurface{T}) where {T<:Real}
    f = (k, xx) -> quantity[k] * surface.fluxexpansion[k]
    return flux_surface_avg(f, surface)
end

"""
    flux_surface_avg(f::F1, surface::FluxSurface{T}) where {F1<:Function, T<:Real}

Flux surface averaging of a function
"""
@inline function flux_surface_avg(f::F1, surface::AbstractFluxSurface{T}) where {F1<:Function,T<:Real}
    return trapz(surface.ll, f) / surface.int_fluxexpansion_dl
end

@compat public flux_surface_avg
push!(document[Symbol("Physics flux-surfaces")], :flux_surface_avg)

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Integrate quantity over volume
"""
function volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::T where {T<:Real}
    dv = eqt.profiles_1d.dvolume_dpsi
    f = (k, xx) -> dv[k] * what[k]
    return trapz(eqt.profiles_1d.psi, f)
end

@compat public volume_integrate
push!(document[Symbol("Physics flux-surfaces")], :volume_integrate)

"""
    volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}

Cumulative integrate quantity over volume
"""
function cumlul_volume_integrate(eqt::IMAS.equilibrium__time_slice, what::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    dv = eqt.profiles_1d.dvolume_dpsi
    f = (k, xx) -> dv[k] * what[k]
    return cumtrapz(eqt.profiles_1d.psi, f)
end

@compat public cumlul_volume_integrate
push!(document[Symbol("Physics flux-surfaces")], :cumlul_volume_integrate)

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

@compat public surface_integrate
push!(document[Symbol("Physics flux-surfaces")], :surface_integrate)

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

@compat public cumlul_surface_integrate
push!(document[Symbol("Physics flux-surfaces")], :cumlul_surface_integrate)

"""
    find_x_point!(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}

Find the `n` X-points that are closest to the separatrix
"""
function find_x_point!(eqt::IMAS.equilibrium__time_slice{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}) where {T<:Real}
    ((rlcfs, zlcfs),) = flux_surface(eqt, eqt.profiles_1d.psi[end], :closed, wall_r, wall_z)
    private = flux_surface(eqt, eqt.profiles_1d.psi[end], :open, wall_r, wall_z)
    Z0 = sum(extrema(zlcfs)) / 2.0
    empty!(eqt.boundary.x_point)

    for (pr, pz) in private
        if (sum(pz) < Z0)
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
        for k in sort(collect(Set(i_x)); rev=true)
            deleteat!(eqt.boundary.x_point, k)
            deleteat!(psidist_lcfs_xpoints, k)
            deleteat!(z_x, k)
        end

        # sort a second time now by distance in psi
        index = sortperm(psidist_lcfs_xpoints; by=abs)
        eqt.boundary.x_point = eqt.boundary.x_point[index]
        z_x = z_x[index]
        # save up to the x_point with Z coordinate opposite to first x point
        # save up to first index where z_x.*z_x[1].<0 is 1
        eqt.boundary.x_point = eqt.boundary.x_point[1:argmax(z_x .* z_x[1] .< 0)]

        #check if primary x-point consistent with strike points
        if isempty(eqt.boundary.strike_point)
            #case with case not already saved  - look at first open surface
            psi_first_open = eqt.boundary_separatrix.psi
            (_, z_first_open) = flux_surface(eqt, psi_first_open, :encircling, Float64[], Float64[])[1]
            if sign(-z_first_open[1]) !== sign(eqt.boundary.x_point[end].z)
                # x_point are not oredered correctly
                if length(eqt.boundary.x_point) <= 2
                    # if just 2, it is enough to flip them
                    reverse!(eqt.boundary.x_point)
                else
                    # warning for exotic divertor config (SuperX, SF,..)
                    @warn("X-points are more than 2 and not ordered correctly following our convention - Point of contact: G. Dose")
                end
            end
        else
            # case with strike_points already saved
            if sign(-eqt.boundary.strike_point[1].z) !== sign(eqt.boundary.x_point[end].z)
                # x_point are not oredered correctly
                if length(eqt.boundary.x_point) <= 2
                    reverse!(eqt.boundary.x_point)
                else
                    # warning for exotic divertor config (SuperX, SF,..)
                    @warn("X-points are more than 2 and not ordered correctly following our convention - Point of contact: G. Dose")
                end
            end
        end

    end

    return eqt.boundary.x_point
end

@compat public find_x_point!
push!(document[Symbol("Physics flux-surfaces")], :find_x_point!)

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

@compat public x_points_inside_wall
push!(document[Symbol("Physics flux-surfaces")], :x_points_inside_wall)

"""
    miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)

Returns named tuple with `R0`, `a`, `κ`, `δu`, `δl`, `ζou`, `ζol`, `ζil`, `ζiu` of a contour
"""
function miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}, r_at_max_z::T, max_z::T, r_at_min_z::T, min_z::T, z_at_max_r::T, max_r::T, z_at_min_r::T, min_r::T) where {T<:Real}
    R0 = 0.5 * (max_r + min_r)
    a = 0.5 * (max_r - min_r)
    b = 0.5 * (max_z - min_z)
    κ = b / a
    δu = (R0 - r_at_max_z) / a
    δl = (R0 - r_at_min_z) / a
    ζou, ζol, ζil, ζiu = luce_squareness(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
    return (R0=R0, a=a, κ=κ, δu=δu, δl=δl, ζou=ζou, ζol=ζol, ζil=ζil, ζiu=ζiu)
end

"""
    miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}) where {T<:Real}
"""
function miller_R_a_κ_δ_ζ(pr::Vector{T}, pz::Vector{T}) where {T<:Real}
    (imaxr, iminr, imaxz, iminz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r) = fluxsurface_extrema(pr, pz)
    return miller_R_a_κ_δ_ζ(pr, pz, r_at_max_z, max_z, r_at_min_z, min_z, z_at_max_r, max_r, z_at_min_r, min_r)
end

@compat public miller_R_a_κ_δ_ζ
push!(document[Symbol("Physics flux-surfaces")], :miller_R_a_κ_δ_ζ)

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

@compat public fluxsurface_extrema
push!(document[Symbol("Physics flux-surfaces")], :fluxsurface_extrema)

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

    names = (:zetaou, :zetaol, :zetaiu, :zetail)
    POs = ((r_at_max_z, z_at_max_r), (r_at_min_z, z_at_max_r), (r_at_max_z, z_at_min_r), (r_at_min_z, z_at_min_r))
    PEs = ((max_r, max_z), (max_r, min_z), (min_r, max_z), (min_r, min_z))

    z = T[]
    for (name, PO, PE) in zip(names, POs, PEs)
        try
            PD = intersection([PO[1], PE[1]], [PO[2], PE[2]], pr, pz).crossings[1]
            PC = (cos(π / 4.0) * (PE[1] - PO[1]) + PO[1], sin(π / 4.0) * (PE[2] - PO[2]) + PO[2])
            push!(z, (norm(PD .- PO) - norm(PC .- PO)) / norm(PE .- PC))
        catch e
            push!(z, T(0.0))
            # plot(pr,pz;aspect_ratio=:equal,title=name)
            # display(plot!([PO[1], PE[1]], [PO[2], PE[2]]))
            # rethrow(e)
        end
    end
    zetaou, zetaol, zetail, zetaiu = z
    return zetaou, zetaol, zetail, zetaiu
end

@compat public luce_squareness
push!(document[Symbol("Physics flux-surfaces")], :luce_squareness)

"""
    areal_elongation(eqt::IMAS.equilibrium__time_slice)

A measure of the plasma elongation based on the averaged cross-sectional area of the plasma, most notably used in the H98y2 scaling

See: https://iopscience.iop.org/article/10.1088/0029-5515/48/9/099801/pdf
"""
function areal_elongation(eqt::IMAS.equilibrium__time_slice)
    eqt1d = eqt.profiles_1d
    a = 0.5 * (eqt1d.r_outboard[end] - eqt1d.r_inboard[end])
    R = 0.5 * (eqt1d.r_outboard[end] + eqt1d.r_inboard[end])
    return eqt1d.volume[end] / (2π * R) / (π * a^2)
end

@compat public areal_elongation
push!(document[Symbol("Physics flux-surfaces")], :areal_elongation)
