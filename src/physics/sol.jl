document[Symbol("Physics sol")] = Symbol[]

"""
    r::Vector{Float64}
    z::Vector{Float64}
    Br::Vector{Float64}
    Bz::Vector{Float64}
    Bp::Vector{Float64}
    Bt::Vector{Float64}
    pitch::Vector{Float64}
    s::Vector{Float64}
    midplane_index::Int
    strike_angles::Vector{Float64}           # Angle in radiants between flux surface and the wall; poloidal angle
    pitch_angles::Vector{Float64}            # Angle in radiants between B and Btoroidal; atan(Bp/Bt)
    grazing_angles::Vector{Float64}          # Angle in radiants between B and the wall; grazing angle
    total_flux_expansion::Vector{Float64}    # Total flux expansion
    poloidal_flux_expansion::Vector{Float64} # Poloidal flux expansion
    wall_index::Vector{Int}                  # index in dd.wall where strike points intersect
"""
struct OpenFieldLine
    r::Vector{Float64}
    z::Vector{Float64}
    Br::Vector{Float64}
    Bz::Vector{Float64}
    Bp::Vector{Float64}
    Bt::Vector{Float64}
    pitch::Vector{Float64}
    s::Vector{Float64}
    midplane_index::Int
    strike_angles::Vector{Float64}           # Angle in radiants between flux surface and the wall; poloidal angle
    pitch_angles::Vector{Float64}            # Angle in radiants between B and Btoroidal; atan(Bp/Bt)
    grazing_angles::Vector{Float64}          # Angle in radiants between B and the wall; grazing angle
    total_flux_expansion::Vector{Float64}    # Total flux expansion
    poloidal_flux_expansion::Vector{Float64} # Poloidal flux expansion
    wall_index::Vector{Int}                  # index in dd.wall where strike points intersect
    psi::Float64                             # Psi level at that open flux surface
end

"""
    OpenFieldLine(
        PSI_interpolant,
        r::AbstractVector{T},
        z::AbstractVector{T},
        wall_r::AbstractVector{T},
        wall_z::AbstractVector{T},
        B0::T,
        R0::T,
        RA::T,
        ZA::T,
        level::T)

OpenFieldLine constructor
"""
function OpenFieldLine(
    PSI_interpolant,
    r::AbstractVector{T},
    z::AbstractVector{T},
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T},
    B0::T,
    R0::T,
    RA::T,
    ZA::T,
    level::T
) where {T<:Real}
    @assert length(wall_r) == length(wall_z)
    if isempty(wall_r)
        # SOL without wall
        rr = r
        zz = z
        strike_angles = [NaN, NaN]
        wall_index = Int64[]
    else
        # SOL with wall
        # returns poloidal angles of each surface
        # rr and zz are clockwise
        # crossing points with wall = (rr[1], zz[1]) (rr[end], zz[end])
        # this is the order at which angles are computed (strike, pitch and grazing)
        # Example - OFL[2].[1<n<length(levels)].strike_angle[1] is computed at (rr[1], zz[1])
        rr, zz, strike_angles, wall_index = line_wall_2_wall(r, z, wall_r, wall_z, RA, ZA)
    end

    if isempty(rr) || all(zz .> ZA) || all(zz .< ZA)
        return nothing
    end

    # add a point exactly at the (preferably outer) midplane
    crossing_index, crossings = intersection([0.0, 1000.0], [ZA, ZA], rr, zz)
    r_midplane = [cr[1] for cr in crossings] # R coordinate of points in SOL surface at MP (inner and outer)
    outer_index = argmax(r_midplane)  #index of point @ MP: this is OMP (for OFL_lfs); IMP for OFL_hfs
    crossing_index = crossing_index[outer_index] #indexes of point at MP in SOL surface
    r_midplane = r_midplane[outer_index] # R coordinate of point at MP in SOL surface
    rr = [rr[1:crossing_index[2]]; r_midplane; rr[crossing_index[2]+1:end]] #Insert in r of SOL surface a point @ MP
    zz = [zz[1:crossing_index[2]]; ZA; zz[crossing_index[2]+1:end]] #Insert in z of SOL surface a point @ MP
    midplane_index = crossing_index[2] + 1 #index at which point @ MP

    # calculate quantities along field line
    Br, Bz = Br_Bz(PSI_interpolant, rr, zz) # r and z component of B for each point in (r,z)
    Bp = sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)     # poloidal component of B for each point in (r,z)
    Bt = abs.(B0 .* R0 ./ rr)              # toroidal component of B for each point in (r,z)
    B = sqrt.(Bp .^ 2 + Bt .^ 2)           # total magnetic field B for each point in (r,z)
    dp = sqrt.(gradient(rr) .^ 2.0 .+ gradient(zz) .^ 2.0) # curvilinear abscissa increments of poloidal projection of SOL surface
    pitch = sqrt.(1.0 .+ (Bt ./ Bp) .^ 2) # ds = dp*sqrt(1 + (Bt/Bp)^2) (pythagora)
    s = cumsum(pitch .* dp) # s = integral(ds)
    s = abs.(s .- s[midplane_index]) # fix 0 at outer midplane

    # Parameters to map heat flux from OMP to wall
    pitch_angles = atan.(Bp, Bt)
    grazing_angles = asin.(sin.([pitch_angles[1], pitch_angles[end]]) .* sin.(strike_angles))
    total_flux_expansion = B[midplane_index] ./ B # total flux expansion(r,z) =  Bomp / B(r,z) [magentic flux conservation]
    poloidal_flux_expansion = total_flux_expansion .* rr[midplane_index] ./ rr .* sin(pitch_angles[midplane_index]) ./ sin.(pitch_angles) # poloidal flux expansion

    return OpenFieldLine(
        rr,
        zz,
        Br,
        Bz,
        Bp,
        Bt,
        pitch,
        s,
        midplane_index,
        strike_angles,
        pitch_angles,
        grazing_angles,
        total_flux_expansion,
        poloidal_flux_expansion,
        wall_index,
        level
    )
end

@compat public OpenFieldLine
push!(document[Symbol("Physics sol")], :OpenFieldLine)

@recipe function plot_ofl(ofl::OpenFieldLine)
    @series begin
        aspect_ratio --> :equal
        label --> ""
        colorbar_title := "log₁₀(Connection length [m] + 1.0)"
        line_z --> log10.(ofl.s .+ 1)
        ofl.r, ofl.z
    end
end

@recipe function plot_ofl(OFL_hfs_lfs_lfsfar::OrderedCollections.OrderedDict{Symbol,Vector{IMAS.OpenFieldLine}})
    @series begin
        OFL_hfs_lfs_lfsfar[:hfs]
    end
    @series begin
        OFL_hfs_lfs_lfsfar[:lfs]
    end
    @series begin
        OFL_hfs_lfs_lfsfar[:lfs_far]
    end
end

@recipe function plot_ofl(OFL::Vector{OpenFieldLine})
    for ofl in OFL
        @series begin
            ofl
        end
    end
end

"""
    sol(eqt::IMAS.equilibrium__time_slice, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true) where {T<:Real}

Returns vectors of hfs and lfs OpenFieldLine

If levels is a vector, it has the values of psi from 0 to max psi_wall_midplane. The function will modify levels of psi to introduce relevant sol surfaces
"""
function sol(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{T},
    wall_z::AbstractVector{T};
    levels::Union{Int,AbstractVector}=20,
    use_wall::Bool=!isempty(wall_r)
) where {T<:Real}
    @assert length(wall_r) == length(wall_z)
    OFL = OrderedCollections.OrderedDict(:hfs => OpenFieldLine[], :lfs => OpenFieldLine[], :lfs_far => OpenFieldLine[])

    # empty wall if use_wall = false
    if use_wall == false
        wall_r = T[]
        wall_z = T[]
    end

    ############
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    ############
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)
    psi__axis_level = eqt.profiles_1d.psi[1] # psi value on axis

    psi_boundaries = (last_closed=eqt.boundary.psi, first_open=eqt.boundary_separatrix.psi)
    if psi_boundaries.last_closed == psi_boundaries.first_open
        return OFL
    end

    # find psi at second magnetic separatrix
    psi__2nd_separatix = find_psi_2nd_separatrix(eqt).not_diverted # find psi at 2nd magnetic separatrix
    psi_sign = sign(psi_boundaries.first_open - psi__axis_level) # sign of the poloidal flux taking psi_axis = 0
    psi_max = find_psi_max(eqt)

    if !isempty(wall_r)
        psi_wall_midplane = find_psi_wall_omp(PSI_interpolant, RA, ZA, wall_r, wall_z, psi_max, psi_sign) # psi at the intersection between wall and omp
        psi_last_lfs, _, psi_first_lfs_far, _ = find_psi_last_diverted(eqt, wall_r, wall_z, PSI_interpolant) # find psi at LDFS, NaN if not a diverted plasma
        threshold = (psi_last_lfs + psi_first_lfs_far) / 2.0
        # limited plasma
        if isnan(psi_last_lfs) && isnan(psi_first_lfs_far)
            error("IMAS.sol cannot yet handle limited plasmas")
        end
    else
        # SOL without wall
        psi_wall_midplane = psi_max
        psi_last_lfs = psi_boundaries.first_open
        psi_first_lfs_far = psi_boundaries.first_open .+ 1E-5
        threshold = psi__2nd_separatix
    end
    ############

    # smart picking of psi levels
    if levels == 1
        levels = [psi_boundaries.first_open]

    elseif levels == 2
        levels = [psi_boundaries.first_open, psi__2nd_separatix]

    elseif typeof(levels) <: Int
        levels = psi_boundaries.first_open .+ psi_sign .* 10.0 .^ LinRange(-9, log10(abs(psi_wall_midplane - psi_boundaries.first_open)), levels)

        indexx = argmin_abs(levels, psi_last_lfs)
        levels = vcat(levels[1:indexx-1], psi_last_lfs, psi_first_lfs_far, levels[indexx+1:end]) # remove closest point + add last_lfs and first_lfs_far
        # add 2nd sep, sort in increasing order and remove doubles (it could happen that psi_boundaries.first_open = psi_last_lfs = psi_2ndseparatrix in DN)
        levels = unique!(sort(vcat(levels, psi__2nd_separatix)))

        # Case of limited plasma - Also add first magnetic separatrix
        psi_sep_closed, psi_sep_open = find_psi_separatrix(eqt)

        if psi_sign * psi_sep_open > psi_sign * psi_boundaries.first_open
            # if psi_boundary is inside separatrix -> limited case
            levels = unique!(sort(vcat(levels, psi_sep_closed, psi_sep_open)))
        end

        if psi_sign == -1
            # if psi is decreasing we must sort in decreasing order
            levels = reverse!(levels)
        end

    else
        #levels is a vector of psi_levels for the discretization of the SOL
        @assert psi_sign * levels[1] >= psi_sign * psi_boundaries.first_open "psi_boundaries.first_open = $(psi_boundaries.first_open) , psi_levels[1] = $(levels[1]) "
        @assert psi_sign * levels[end] <= psi_sign * psi_wall_midplane "psi_wall_midplane = $psi_wall_midplane , psi_levels[end] = $(levels[end]) "
        levels_is_not_monotonic_in_Ip_direction = all(psi_sign * diff(levels) .>= 0)
        @assert levels_is_not_monotonic_in_Ip_direction # levels must be monotonic according to plasma current direction
        # make sure levels includes separatrix and wall
        levels[1] = psi_boundaries.first_open
        # push!(levels,psi_wall_midplane - psi_sign * 1E-3 * abs(psi_wall_midplane))
        levels = unique!(sort!(levels))

        if psi_sign == -1
            # if psi is decreasing we must sort in decreasing order
            levels = reverse!(levels)
        end
    end

    # TO DO for the future: insert private flux regions (upper and lower)
    for level in levels
        lines = flux_surface(eqt, level, :open, wall_r, wall_z)

        for (r, z) in lines
            ofl = OpenFieldLine(PSI_interpolant, r, z, wall_r, wall_z, B0, R0, RA, ZA, level)
            if ofl === nothing
                continue
            end

            # identify field line as `hfs`, `lfs`, or `lfs_far`
            if ofl.r[ofl.midplane_index] < RA
                # Add SOL surface in OFL_hfs
                ofl_type = :hfs
            else
                if psi_sign * level < psi_sign * threshold
                    # if z[1] and z[end] have same sign, check psi
                    # Add SOL surface in OFL_lfs_far
                    ofl_type = :lfs
                else
                    ofl_type = :lfs_far
                end

            end

            push!(OFL[ofl_type], ofl) # add result
        end
    end

    # In the case of a perfectly balanced double null OFL_lfs will be empty and all of the flux surfaces will appear in OFL_lfs_far
    # Here we switch the two
    if isempty(OFL[:lfs])
        tmp = OFL[:lfs_far]
        OFL[:lfs_far] = OFL[:lfs]
        OFL[:lfs] = tmp
    end
    @assert !isempty(OFL[:lfs]) "No magnetic surfaces are found in the SOL. Try checking IMAS.line_wall_2_wall is working properly."
    return OFL
end

"""
    sol(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true)
"""
function sol(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall; levels::Union{Int,AbstractVector}=20, use_wall::Bool=!isempty(first_wall(wall).r))
    fw = first_wall(wall)
    return sol(eqt, fw.r, fw.z; levels, use_wall)
end

"""
    sol(dd::IMAS.dd; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true)
"""
function sol(dd::IMAS.dd; levels::Union{Int,AbstractVector}=20, use_wall::Bool=!isempty(first_wall(dd.wall).r))
    return sol(dd.equilibrium.time_slice[], dd.wall; levels, use_wall)
end

@compat public sol
push!(document[Symbol("Physics sol")], :sol)

"""
    find_levels_from_P(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation,
        r::Vector{<:Real},
        q::Vector{<:Real},
        levels::Int
    )

Function for the discretization of the poloidal flux ψ on the SOL, based on an hypotesis of OMP radial transport through arbitrary q(r)
returns vector with level of ψ, vector with matching r_midplane and q.
Discretization with even steps of P = integral_sep^wal 2πrq(r)dr (same power in each flux tube)
"""
function find_levels_from_P(
    eqt::IMAS.equilibrium__time_slice,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation,
    r::Vector{<:Real},
    q::Vector{<:Real},
    levels::Int
)
    ################### Housekeeping on function q(r) ###################
    @assert length(wall_r) == length(wall_z)
    @assert length(r) == length(q)
    @assert all(q .>= 0) # q is all positive
    @assert all(r .>= 0) # r is all positive
    @assert r[1] == 0

    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis
    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    # r_mid(ψ) interpolator for region of interest
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    rmax = eqt2d.grid.dim1[end]
    r_mid_of_interest = 10.0 .^ range(log10(maximum(eqt.boundary.outline.r) * 0.99), log10(rmax), 1000)
    r_mid = interp_rmid_at_psi(PSI_interpolant, r_mid_of_interest, ZA)
    psi_mid = PSI_interpolant.(r_mid_of_interest, r_mid_of_interest .* 0.0 .+ ZA)
    psi_sign = sign(psi_mid[end] - psi_mid[1])
    if psi_sign < 0
        r_mid = DataInterpolations.CubicSpline(r_mid_of_interest, -psi_mid; extrapolation=ExtrapolationType.Extension)
    end

    psi_boundaries = (last_closed=eqt.boundary.psi, first_open=eqt.boundary_separatrix.psi)
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt).not_diverted # psi of the second magnetic separatrix
    if psi_sign > 0
        r_separatrix_midplane = r_mid(psi_boundaries.first_open) # R OMP at separatrix
        r_2ndseparatrix_midplane = r_mid(psi_2ndseparatrix) # R coordinate at OMP of 2nd magnetic separatrix
    else
        r_separatrix_midplane = r_mid(-psi_boundaries.first_open) # R OMP at separatrix
        r_2ndseparatrix_midplane = r_mid(-psi_2ndseparatrix) # R coordinate at OMP of 2nd magnetic separatrix
    end

    if isempty(wall_r) .|| isempty(wall_z)
        # no wall
        psi_wall_midplane = psi_sign * maximum(psi_sign .* eqt2d.psi) - psi_sign # if no wall, upper bound of psi is maximum value in eqt -1 (safe)
        r_wall_midplane = rmax
        null_within_wall = true
        r_last_diverted = [1, 1] * r_2ndseparatrix_midplane
    else
        # there is a wall
        crossings = intersection([RA, maximum(wall_r)], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
        r_wall_midplane = crossings[1][1] # R coordinate of the wall at OMP
        psi_wall_midplane = find_psi_wall_omp(eqt, wall_r, wall_z)
        psi_last_lfs, _, psi_first_lfs_far, null_within_wall = find_psi_last_diverted(eqt, wall_r, wall_z, PSI_interpolant) # psi of grazing surface

        if psi_sign > 0
            r_last_diverted = r_mid.([psi_last_lfs, psi_first_lfs_far]) # R coordinate at OMP of grazing surface
        else

            r_last_diverted = r_mid.([-psi_last_lfs, -psi_first_lfs_far]) # R coordinate at OMP of grazing surface
        end
    end
    r = r .+ r_separatrix_midplane
    order = sortperm(r)
    q = q[order]
    r = r[order] # r is now increasing monotonically

    if r[1] >= r_wall_midplane
        # all the vector r is inside the wall
        r = [r_separatrix_midplane, r_wall_midplane] # constant value of q as closest point to SOL
        q = [1, 0.999] * q[1]  # keep preference for monotonic decrease

    end
    if r[end] <= r_separatrix_midplane
        # all of r is inside the separatrix
        r = [r_separatrix_midplane, r_wall_midplane] # constant value of q as closest point to SOL
        q = [1.001, 1] * q[end] # keep preference for monotonic decrease
    end

    # r[1] can be either >, = or < than r_separatrix_midplane; add r_separatrix_midplane
    if r[1] > r_separatrix_midplane
        # r starts from inside the sol
        r = append!([r_separatrix_midplane], r) # add a point at R_OMP
        q = append!([q[1] * 1.001], q)            # Repeat first value of q with slight increment, favoring monotonic decrease
    end
    # if r[1] == r_separatrix_midplane do nothing
    if r[1] < r_separatrix_midplane
        # r starts from inside the separatrix, q(r) must be cut
        index = argmin_abs(r, r_separatrix_midplane) # closest point
        # index2 is the position in r, such that r_separatrix_midplane is between r[index2] and r[index]
        if r[index] > r_separatrix_midplane
            index2 = index - 1
            #interp linearly value at r_separatrix_midplane between r[index2] and r[index]
            qq = q[index] + (q[index2] - q[index]) / (r[index2] - r[index]) * (r_separatrix_midplane - r[index])
            r = append!([r_separatrix_midplane], r[index:end]) # cut r and q
            q = append!([qq], q[index:end])
        else
            index2 = index + 1
            #interp linearly value at r_separatrix_midplane between r[index2] and r[index]
            qq = q[index] + (q[index2] - q[index]) / (r[index2] - r[index]) * (r_separatrix_midplane - r[index])
            r = append!([r_separatrix_midplane], r[index2:end]) # cut r and q
            q = append!([qq], q[index2:end])
        end
    end

    # r[end] can be either >, = < than r_wall_midplane; add r_wall_midplane
    if r[end] < r_wall_midplane
        # r ends inside sol
        r = push!(r, r_wall_midplane) # add a point at r_wall_midplane
        q = push!(q, q[end] * 0.999)    # Repeat last value of q with slight reduction, favoring monotonic decrease
    end
    # if r[end]==r_wall_midplane do nothing
    if r[end] > r_wall_midplane
        # r ends inside the wall, q(r) must be cut
        index = argmin_abs(r, r_wall_midplane) # closest point
        if r[index] > r_wall_midplane
            index2 = index - 1
            #interp linearly value at r_wall_midplane between r[index2] and r[index]
            qq = q[index] + (q[index2] - q[index]) / (r[index2] - r[index]) * (r_wall_midplane - r[index])
            r = push!(r[1:index-1], r_wall_midplane) # cut + add point between index and index2
            q = push!(q[1:index-1], qq)
        else
            index2 = index + 1
            #interp linearly value at r_wall_midplane between r[index2] and r[index]
            qq = q[index] + (q[index2] - q[index]) / (r[index2] - r[index]) * (r_wall_midplane - r[index])
            r = push!(r[1:index], r_wall_midplane) # cut + add a point between index and index
            q = push!(q[1:index], qq)
        end
    end

    ############### finished houskeeping of q(r) ####################
    #################################################################

    # build P(r) = integral_sep^r q(ρ)2πρdρ
    P = q * 0
    for index in 2:length(q)
        P[index] = P[index-1] + trapz(r[index-1:index], 2 * π .* r[index-1:index] .* q[index-1:index])
    end
    # being 2πr q(r) positive-definite, P(r) is strictly monotonic, therefore also injective (one-to-one)
    # P(r) is always invertible for every q(r)>0
    r = Interpolations.deduplicate_knots!(r)
    interp_P = cubic_interp1d(r, P) # interpolant of P(r)

    p_levels = collect(LinRange(0, maximum(P), levels))   # levels to interpolate P
    # add flux surfaces of interest: last diverted surface and 2nd magnetic separatrix
    # NOTE: number of level increases
    P_low = interp_P(r_last_diverted[1]) # last diverted surface (up to precision); inside OFL[:lfs]
    P_up = interp_P(r_last_diverted[2]) # first surface crossing the top first wall; inside OFL[:lfs_far]
    if null_within_wall
        # last diverted surface is the 2nd separatrix; add only 2nd separatrix_up and 2nd separatix_low
        if !(P_low in p_levels)
            p_levels = push!(p_levels, P_low)
        end
        if !(P_up in p_levels)
            p_levels = push!(p_levels, P_up)
        end
        # TO BE CHECKED: in case of double null, power balance must work
    else
        # last diverted surface is different from the 2nd separatix: add 2nd magnetic separatix and last diverted surface
        P_2ndseparatrix = interp_P(r_2ndseparatrix_midplane) # interp value at 2nd separatrix
        # add all surfaces, but avoid repetition
        if !(P_2ndseparatrix in p_levels)
            p_levels = push!(p_levels, P_2ndseparatrix)
        end
        if !(P_low in p_levels)
            p_levels = push!(p_levels, P_low)
        end
        if !(P_up in p_levels)
            p_levels = push!(p_levels, P_up)
        end
    end

    p_levels = sort!(p_levels)
    P = Interpolations.deduplicate_knots!(P)
    interp_inverseP = cubic_interp1d(P, r) # interpolant of inverse function of r(P)
    R = interp_inverseP.(p_levels)

    # using ψ(R), go from discretization in R to discretization in ψ
    psi_levels = PSI_interpolant.(R, R .* 0.0 .+ ZA) # ψ(R)

    #force psi_sep and psi_wall_midplane
    psi_levels[1] = psi_boundaries.first_open
    psi_levels[end] = psi_wall_midplane
    # filter possible errors
    if psi_sign > 0
        psi_levels = psi_levels[psi_levels.>=psi_boundaries.first_open]
    else
        psi_levels = psi_levels[psi_levels.<=psi_boundaries.first_open]
    end
    return psi_levels, R, p_levels
end

"""
    find_levels_from_P(
        eqt::IMAS.equilibrium__time_slice,
        wall::IMAS.wall,
        PSI_interpolant::Interpolations.AbstractInterpolation,
        r::Vector{<:Real},
        q::Vector{<:Real},
        levels::Int
    )
"""
function find_levels_from_P(
    eqt::IMAS.equilibrium__time_slice,
    wall::IMAS.wall,
    PSI_interpolant::Interpolations.AbstractInterpolation,
    r::Vector{<:Real},
    q::Vector{<:Real},
    levels::Int
)
    return find_levels_from_P(eqt, first_wall(wall).r, first_wall(wall).z, PSI_interpolant, q, r, levels)
end

"""
    find_levels_from_P(dd::IMAS.dd, r::Vector{<:Real}, q::Vector{<:Real}, levels::Int)
"""
function find_levels_from_P(dd::IMAS.dd, r::Vector{<:Real}, q::Vector{<:Real}, levels::Int)
    _, _, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
    return find_levels_from_P(dd.equilibrium.time_slice[], dd.wall, PSI_interpolant, q, r, levels)
end

@compat public find_levels_from_P
push!(document[Symbol("Physics sol")], :find_levels_from_P)

"""
    find_levels_from_wall(
        eqt::IMAS.equilibrium__time_slice,
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        PSI_interpolant::Interpolations.AbstractInterpolation
    )

Function for that computes the value of psi at the points of the wall mesh in dd
"""
function find_levels_from_wall(
    eqt::IMAS.equilibrium__time_slice{T},
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    PSI_interpolant::Interpolations.AbstractInterpolation
) where {T<:Real}

    @assert length(wall_r) == length(wall_z)
    if isempty(wall_r)
        # no wall
        return T[]
    end

    psi_boundaries = (last_closed=eqt.boundary.psi, first_open=eqt.boundary_separatrix.psi)
    psi_wall_midplane = find_psi_wall_omp(eqt, wall_r, wall_z)

    levels = PSI_interpolant.(wall_r, wall_z)
    psi_tangent_in = find_psi_tangent_omp(eqt, wall_r, wall_z, PSI_interpolant).psi_tangent_in
    push!(levels, psi_tangent_in)

    index = levels .>= psi_boundaries.first_open .&& levels .<= psi_wall_midplane
    levels = levels[index]
    sort!(levels)

    return levels
end

"""
    find_levels_from_wall(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, PSI_interpolant::Interpolations.AbstractInterpolation)
"""
function find_levels_from_wall(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, PSI_interpolant::Interpolations.AbstractInterpolation)
    return find_levels_from_wall(eqt, first_wall(wall).r, first_wall(wall).z, PSI_interpolant)
end

"""
    find_levels_from_wall(dd::IMAS.dd)
"""
function find_levels_from_wall(dd::IMAS.dd)
    _, _, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
    return find_levels_from_wall(dd.equilibrium.time_slice[], dd.wall, PSI_interpolant)
end

@compat public find_levels_from_wall
push!(document[Symbol("Physics sol")], :find_levels_from_wall)

"""
    line_wall_2_wall(r::AbstractVector{T}, z::AbstractVector{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}, RA::Real, ZA::Real) where {T<:Real}

Returns r, z coordinates of open field line contained within wall, as well as angles of incidence at the strike locations

RA and ZA are the coordinate of the magnetic axis
"""
function line_wall_2_wall(r::AbstractVector{T}, z::AbstractVector{T}, wall_r::AbstractVector{T}, wall_z::AbstractVector{T}, RA::Real, ZA::Real) where {T<:Real}
    @assert length(wall_r) == length(wall_z)

    indexes, crossings = intersection(r, z, wall_r, wall_z, 1E-6) # find where flux surface crosses wall ("strike points" of surface)
    # crossings -  Vector{Tuple{Float64, Float64}} - crossings[1] contains (r,z) of first "strike point"
    # indexes   -  Vector{Tuple{Float64, Float64}} - indexes[1] contains indexes of (r,z) and (wall_r, wall_z) of first "strike point"
    r_z_index = [k[1] for k in indexes] #index of vectors (r,z) of all crossing point
    wall_index = [k[2] for k in indexes] #index of vectors (wall_r, wall_z) of all crossing point

    crossings2 = intersection([0, RA], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw inner midplane (IMP) with wall
    r_wall_imp = crossings2[1][1] # R coordinate of the wall at IMP

    crossings2 = intersection([RA, 2 * maximum(wall_r)], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_omp = crossings2[1][1] # R coordinate of the wall at OMP

    if isempty(r_z_index) # if the flux surface does not cross the wall return empty vector (it is not a surf in SOL)
        return Float64[], Float64[], Float64[], Int64[]

    elseif length(r_z_index) == 1
        error("""line_wall_2_wall: open field line should intersect wall at least twice.
                 If it does not it's likely because the extent of the equilibrium grid is too small.
                 Suggestion: plot dd.wall + eqt.profiles_2d to debug.""")
    end

    # angle of incidence
    strike_angles = intersection_angles(r, z, wall_r, wall_z, indexes) # find poloidal angle btw SOL magnetic surface and wall

    if length(r_z_index) == 2 #it is a magnetic surface in the SOL
    # pass

    else
        # more than 2 intersections with wall

        #  index in (r,z) of closest midplane point (favoring low field side)
        j0 = argmin(abs.(z .- ZA) .+ (r .< RA)) # min of abs(z-ZA) + 1 meter only on hfs

        # the closest intersection point (in steps) to z=ZA
        i1 = sortperm(abs.(r_z_index .- j0))[1] # identifies which crossing point is closest to OMP (outer "strike point")

        # the intersection on the other side of the midplane
        j1 = r_z_index[i1] #index of outer "strike point"; point closer to the OMP

        if j0 < j1
            # (r,z) is ordered such that the outer "strike point" comes after OMP
            i2 = i1 - 1 # inner "strike point" is the point before in r_z_index
            if i2 == 0
                # if closest intersection of line with wall is below the midplane (j0<j1), i1 = 1 and i2 = 0
                i2 = 2 # fix that such that there is no index = 0
            end
        else
            # (r,z) is ordered such that the OMP comes after the outer "strike point"
            i2 = i1 + 1 #  inner "strike point" is the second point in r_z_index
        end
        if abs(j0 - j1) == 1
            # intersection with wall has same index than closest point to the OMP
            # the intersection occurs at the wall OMP, within 1 index
            if r[j1+1] > r_wall_omp
                # the next point is outside the wall at omp
                # take PREVIOUS intersection
                i2 = i1 - 1
            end

            if r[j1-1] > r_wall_omp
                #the previous point is outside the wall at omp
                # take NEXT intersection
                i2 = i1 + 1
            end

            # NOTE: if r[j1 + 1] <= r_wall_omp ||  r[j1 - 1] <= r_wall_omp, do nothing!

        end
        if i1 == length(r_z_index)
            i2 = i1 - 1
        end

        i = sort([i1, i2])
        r_z_index = r_z_index[i]
        wall_index = wall_index[i]
        crossings = crossings[i]
        strike_angles = strike_angles[i]
    end

    rr = vcat(crossings[1][1], r[r_z_index[1]+1:r_z_index[2]], crossings[2][1]) # r coordinate of magnetic surface between one "strike point" and the other
    zz = vcat(crossings[1][2], z[r_z_index[1]+1:r_z_index[2]], crossings[2][2]) # z coordinate of magnetic surface between one "strike point" and the other

    # remove surfaces that cross midplane outiside the wall
    crossings2 = intersection([0, RA], [ZA, ZA], rr, zz).crossings # (r,z) point of intersection btw inner midplane (IMP) with magnetic surface in the SOL
    r_imp = [cr[1] for cr in crossings2] # R coordinate of the wall at IMP
    if !isempty(r_imp)
        r_imp = r_imp[1] # make it float
    else
        r_imp = RA # the surface crosses only at the OMP
    end

    crossings2 = intersection([RA, 2 * maximum(wall_r)], [ZA, ZA], rr, zz).crossings # (r,z) point of intersection btw outer midplane (OMP) with magnetic surface in the SOL
    r_omp = [cr[1] for cr in crossings2] # R coordinate of the wall at OMP
    if !isempty(r_omp)
        r_omp = r_omp[1] # make it float
    else
        r_omp = RA # the surface crosses only at the IMP
    end

    if r_imp .< r_wall_imp || r_omp .> r_wall_omp
        return Float64[], Float64[], Float64[], Int64[]
    end
    # sort clockwise (COCOS 11)
    angle = mod.(atan.(zz .- ZA, rr .- RA), 2 * π) # counterclockwise angle from midplane
    angle_is_monotonic = all(abs.(diff(angle)) .< π) # this finds if the field line crosses the OMP
    if angle_is_monotonic
        if angle[1] < angle[end]
            rr = reverse(rr)
            zz = reverse(zz)
            strike_angles = reverse(strike_angles)
            wall_index = reverse(wall_index)
        end
    else
        if angle[1] > angle[end]
            rr = reverse(rr)
            zz = reverse(zz)
            strike_angles = reverse(strike_angles)
            wall_index = reverse(wall_index)
        end
    end

    return rr, zz, strike_angles, wall_index
end

@compat public line_wall_2_wall
push!(document[Symbol("Physics sol")], :line_wall_2_wall)

"""
    identify_strike_surface(ofl::OpenFieldLine, divertors::IMAS.divertors)

Returns vector of two tuples with three integers each, identifying the indexes of the divertor/target/tile that the field line intersections

When a field line does not intersect a divertor target, then the tuple returned is (0, 0, 0)
"""
function identify_strike_surface(ofl::OpenFieldLine, divertors::IMAS.divertors)
    identifiers = Tuple{Int,Int,Int}[]
    for strike_index in (1, length(ofl.r))
        distances = OrderedCollections.OrderedDict()
        for (k_divertor, divertor) in enumerate(divertors.divertor)
            for (k_target, target) in enumerate(divertor.target)
                for (k_tile, tile) in enumerate(target.tile)
                    id = (k_divertor, k_target, k_tile)
                    d = point_to_path_distance(ofl.r[strike_index], ofl.z[strike_index], tile.surface_outline.r, tile.surface_outline.z)
                    distances[id] = d
                end
            end
        end
        d = minimum(values(distances))
        if d < 1E-2
            k = argmin(collect(values(distances)))
            push!(identifiers, collect(keys(distances))[k])
        else
            push!(identifiers, (0, 0, 0))
        end
    end
    return identifiers
end

@compat public identify_strike_surface
push!(document[Symbol("Physics sol")], :identify_strike_surface)

"""
    divertor_totals_from_targets!(divertor::IMAS.divertors__divertor)

Sets total values of properties by summing over the corresponding value on the targets

  - power_black_body
  - power_conducted
  - power_convected
  - power_currents
  - power_incident
  - power_neutrals
  - power_radiated
  - power_recombination_neutrals
  - power_recombination_plasma
"""
function divertor_totals_from_targets!(divertor::IMAS.divertors__divertor{T}) where {T<:Real}
    time0 = global_time(divertor)
    for field in (
        :power_black_body,
        :power_conducted,
        :power_convected,
        :power_currents,
        :power_incident,
        :power_neutrals,
        :power_radiated,
        :power_recombination_neutrals,
        :power_recombination_plasma
    )
        total = zero(T)
        for target in divertor.target
            property = getproperty(target, field)
            if !ismissing(property, :data)
                value = get_time_array(property, :data, time0)
                total += value
            end
        end
        total_property = getproperty(divertor, field)
        @ddtime(total_property.data = total)
    end
end

@compat public divertor_totals_from_targets!
push!(document[Symbol("Physics sol")], :divertor_totals_from_targets!)

"""
    Bpol(a::T, κ::T, Ip::T) where {T<:Real}

Average poloidal magnetic field magnitude
"""
function Bpol(a::T, κ::T, Ip::T) where {T<:Real}
    return (mks.μ_0 * Ip) / (2π * a * sqrt((1.0 + κ^2) / 2.0))
end

@compat public Bpol
push!(document[Symbol("Physics sol")], :Bpol)

"""
    Bpol_omp(eqt::IMAS.equilibrium__time_slice)

Poloidal magnetic field magnitude evaluated at the outer midplane
"""
function Bpol_omp(eqt::IMAS.equilibrium__time_slice)
    _, _, PSI_interpolant = ψ_interpolant(eqt.profiles_2d)
    eqt1d = eqt.profiles_1d
    R_omp = eqt1d.r_outboard[end]
    Z_omp = eqt.global_quantities.magnetic_axis.z
    return Bp(PSI_interpolant, R_omp, Z_omp)
end

@compat public Bpol_omp
push!(document[Symbol("Physics sol")], :Bpol_omp)

"""
    power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; time0::Float64=global_time(cp1d))

Total power coming out of the SOL [W]

NOTE: This function returns 1.0 [W] if power is less than that so that SOL quantities remain finite
"""
function power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; time0::Float64=global_time(cp1d))
    tot_pow_in = total_power_inside(core_sources, cp1d; time0, include_radiation=true, include_time_derivative=true)
    tot_pow_in = max(1.0, tot_pow_in)
    return tot_pow_in
end

"""
    power_sol(dd::IMAS.dd; time0::Float64=dd.global_time)
"""
function power_sol(dd::IMAS.dd; time0::Float64=dd.global_time)
    return power_sol(dd.core_sources, dd.core_profiles.profiles_1d[time0]; time0)
end

@compat public power_sol
push!(document[Symbol("Physics sol")], :power_sol)

# ====== #
# Loarte #
# ====== #
"""
    widthSOL_loarte(B0::T, q95::T, Psol::T) where {T<:Real}

Returns midplane power decay length λ_q in meters
"""
function widthSOL_loarte(B0::T, q95::T, Psol::T) where {T<:Real}
    return 0.00265 * Psol^0.38 * B0^(-0.71) * q95^0.3
end

"""
    widthSOL_loarte(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
"""
function widthSOL_loarte(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    B0 = B0_geo(eqt)
    q95 = eqt.global_quantities.q_95
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_loarte(B0, q95, Psol)
end

"""
    widthSOL_loarte(dd::IMAS.dd)
"""
function widthSOL_loarte(dd::IMAS.dd)
    return widthSOL_loarte(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

@compat public widthSOL_loarte
push!(document[Symbol("Physics sol")], :widthSOL_loarte)

# ======= #
# Sieglin #
# ======= #
"""
    widthSOL_sieglin(R0::T, a::T, Bpol_omp::T, Psol::T, ne_ped::T) where {T<:Real}

Returns integral power decay length λ_int in meters
Eich scaling(NF 53 093031) & B. Sieglin PPCF 55 (2013) 124039
"""
function widthSOL_sieglin(R0::T, a::T, Bpol_omp::T, Psol::T, ne_ped::T) where {T<:Real}
    λ_q = widthSOL_eich(R0, a, Bpol_omp, Psol)

    # From B. Sieglin PPCF 55 (2013) 124039
    # S includes the geometrical effects of the divertor assembly itself
    ne_ped /= 1.e19
    S = 0.09 * 1E-3 * ne_ped^1.02 * Bpol_omp^-1.01

    # extrapolate S from ASDEX results
    S *= (R0 / 1.65)

    # This is a valid approximation when S/λ_q < 10
    if S / λ_q > 10
        @warn "S/λ_q = $(S/λ_q) > 10 integral power decay approximation is inaccurate"
    end
    λ_int = λ_q + 1.64 * S

    return λ_int
end

"""
    widthSOL_sieglin(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
"""
function widthSOL_sieglin(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eqt1d = eqt.profiles_1d
    R0 = (eqt1d.r_outboard[end] .+ eqt1d.r_inboard[end]) / 2.0
    a = (eqt1d.r_outboard[end] .- eqt1d.r_inboard[end])
    Psol = power_sol(core_sources, cp1d)
    ne_ped = interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(0.95)
    return widthSOL_sieglin(R0, a, Bpol_omp(eqt), Psol, ne_ped)
end

"""
    widthSOL_sieglin(dd::IMAS.dd)
"""
function widthSOL_sieglin(dd::IMAS.dd)
    return widthSOL_sieglin(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

@compat public widthSOL_sieglin
push!(document[Symbol("Physics sol")], :widthSOL_sieglin)

# ==== #
# Eich #
# ==== #
"""
    widthSOL_eich(R0::T, a::T, Bpol_omp::T, Psol::T) where {T<:Real}

Returns midplane power decay length λ_q in meters

Eich scaling (NF 53 093031)
"""
function widthSOL_eich(R0::T, a::T, Bpol_omp::T, Psol::T) where {T<:Real}
    return 1.35 * 1E-3 * (Psol / 1E6)^-0.02 * R0^0.04 * Bpol_omp^-0.92 * (a / R0)^0.42
end

"""
    widthSOL_eich(eqt::IMAS.equilibrium__time_slice, Psol::Real)
"""
function widthSOL_eich(eqt::IMAS.equilibrium__time_slice, Psol::Real)
    eqt1d = eqt.profiles_1d
    R0 = (eqt1d.r_outboard[end] .+ eqt1d.r_inboard[end]) / 2.0
    a = (eqt1d.r_outboard[end] .- eqt1d.r_inboard[end])
    return widthSOL_eich(R0, a, Bpol_omp(eqt), Psol)
end

"""
    widthSOL_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
"""
function widthSOL_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_eich(eqt, Psol)
end

"""
    widthSOL_eich(dd::IMAS.dd)
"""
function widthSOL_eich(dd::IMAS.dd)
    return widthSOL_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

@compat public widthSOL_eich
push!(document[Symbol("Physics sol")], :widthSOL_eich)

"""
    q_pol_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)

Poloidal heat flux [W/m²] at the outer midplane based on Eich λ_q
"""
function q_pol_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eqt1d = eqt.profiles_1d
    R_omp = eqt1d.r_outboard[end]
    Psol = power_sol(core_sources, cp1d)
    channel_area = 2π * R_omp * widthSOL_eich(eqt, cp1d, core_sources)
    return Psol / channel_area
end

"""
    q_pol_omp_eich(dd::IMAS.dd)
"""
function q_pol_omp_eich(dd::IMAS.dd)
    return q_pol_omp_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

@compat public q_pol_omp_eich
push!(document[Symbol("Physics sol")], :q_pol_omp_eich)

"""
    q_par_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)

Parallel heat flux [W/m²] at the outer midplane based on Eich λ_q
"""
function q_par_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eqt1d = eqt.profiles_1d
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    R_omp = eqt1d.r_outboard[end]
    Bt_omp = B0 * R0 / R_omp
    return q_pol_omp_eich(eqt, cp1d, core_sources) / sin(atan(Bpol_omp(eqt) / Bt_omp))
end

"""
    q_par_omp_eich(dd::IMAS.dd)
"""
function q_par_omp_eich(dd::IMAS.dd)
    return q_par_omp_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

@compat public q_par_omp_eich
push!(document[Symbol("Physics sol")], :q_par_omp_eich)

# ==== #
# Zohm #
# ==== #
"""
    zohm_divertor_figure_of_merit(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Computes a figure of merit for the divertor (Zohm) PB/R/q/A [W T/m]
"""
function zohm_divertor_figure_of_merit(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    R0 = eqt.boundary.geometric_axis.r
    a = eqt.boundary.minor_radius
    A = R0 / a
    q95 = eqt.global_quantities.q_95
    B0 = eqt.global_quantities.vacuum_toroidal_field.b0
    Psol = power_sol(core_sources, cp1d)

    zohm = Psol * B0 / R0 / A / q95 # W T/m
    return zohm
end

"""
    zohm_divertor_figure_of_merit(dd::IMAS.dd)
"""
function zohm_divertor_figure_of_merit(dd::IMAS.dd)
    return zohm_divertor_figure_of_merit(dd.core_sources, dd.core_profiles.profiles_1d[], dd.equilibrium.time_slice[])
end

@compat public q_par_omp_eich
push!(document[Symbol("Physics sol")], :q_par_omp_eich)

# ============= #
# Strike points #
# ============= #
"""
    find_strike_points(pr::AbstractVector{T1}, pz::AbstractVector{T1}, wall_r::AbstractVector{T2}, wall_z::AbstractVector{T2}) where {T1<:Real, T2<:Real}Real}

Finds strike points and angles of incidence between two paths
"""
function find_strike_points(pr::AbstractVector{T1}, pz::AbstractVector{T1}, wall_r::AbstractVector{T2}, wall_z::AbstractVector{T2}) where {T1<:Real,T2<:Real}
    indexes, crossings = intersection(wall_r, wall_z, pr, pz)
    Rxx = [cr[1] for cr in crossings]
    Zxx = [cr[2] for cr in crossings]
    return (Rxx=Rxx, Zxx=Zxx)
end

"""
    find_strike_points(
        eqt::IMAS.equilibrium__time_slice{T1},
        wall_r::AbstractVector{T2},
        wall_z::AbstractVector{T2},
        psi_last_closed::Real,
        psi_first_open::Real;
        strike_surfaces_r::AbstractVector{T3}=wall_r,
        strike_surfaces_z::AbstractVector{T3}=wall_z,
        private_flux_regions::Bool=true
    ) where {T1<:Real,T2<:Real,T3<:Real}

Finds equilibrium strike points, angle of incidence between wall and strike leg, and
the minimum distance between the surface where a strike point is located (both private and encircling) and the last closed flux surface (encircling)
"""
function find_strike_points(
    eqt::IMAS.equilibrium__time_slice{T1},
    wall_r::AbstractVector{T2},
    wall_z::AbstractVector{T2},
    psi_last_closed::Real,
    psi_first_open::Real;
    strike_surfaces_r::AbstractVector{T3}=wall_r,
    strike_surfaces_z::AbstractVector{T3}=wall_z,
    private_flux_regions::Bool=true
) where {T1<:Real,T2<:Real,T3<:Real}

    Rxx = Float64[]
    Zxx = Float64[]
    dxx = Float64[] # minimum distance between the surface where a strike point is located (both private and encircling) and the last closed surface (encircling)

    if !isempty(wall_r)
        # find separatrix as first surface in SOL, not in private region
        if psi_first_open !== nothing
            bnd = flux_surface(eqt, psi_last_closed, :closed, wall_r, wall_z)[1]
            sep = flux_surface(eqt, psi_first_open, :any, Float64[], Float64[])
            raxis = eqt.boundary.geometric_axis.r
            zaxis = eqt.boundary.geometric_axis.z
            for (pr, pz) in sep
                if isempty(pr)
                    continue
                end
                if private_flux_regions && (all(z < zaxis for z in pz) || all(z > zaxis for z in pz))
                    #pass, private flux region
                elseif any(z > zaxis for z in pz) && any(z < zaxis for z in pz) &&
                       sign(pz[1] - zaxis) == sign(pz[end] - zaxis)
                    #pass, going around the confined plasma
                else
                    continue
                end
                Rx_, Zx_ = find_strike_points(pr, pz, strike_surfaces_r, strike_surfaces_z)
                # compute dxx
                dx_, k1, _ = minimum_distance_polygons_vertices(pr, pz, bnd.r, bnd.z)

                # We save 2 strike points per surface
                # if there are more than 2 strike points per surface, filter shadowed strike-points
                if length(Rx_) > 2
                    # the surface identified by (pr,pz) itersects the wall more than twice,
                    # Retrieve only the segments of (pr,pz) inside the wall
                    ps = Tuple{Float64,Float64}[]
                    fw = collect(zip(wall_r, wall_z))
                    segments = intersection_split(pr, pz, wall_r, wall_z)
                    for segment in segments
                        # we retain only segments that are within the wall (disregard the extrema)
                        if length(segment.r) > 2 && PolygonOps.inpolygon((segment.r[2], segment.z[2]), fw) == 1
                            # save intersections of each segment
                            push!(ps, (segment.r[1], segment.z[1]))
                            push!(ps, (segment.r[end], segment.z[end]))
                        end
                    end
                    # how many segments inside the wall are found?
                    L = length(ps) / 2

                    if L == 1
                        # only one segment inside wall: strike points found.
                        Rx_ = [p[1] for p in ps]
                        Zx_ = [p[2] for p in ps]
                    else
                        # pick point on (pr,pz) closest to boundary (lcfs)
                        X = [pr[k1], pz[k1]]
                        # pick segment with intersections closest to X
                        dist = Vector{Float64}(undef, Int(2 * L))
                        for (k, point) in enumerate(ps)
                            dist[k] = sqrt((point[1] - X[1])^2 + (point[2] - X[2])^2)
                        end
                        # odd positions in ps are starting points of the segment, even position are the end
                        indx = argmin(dist) # index of closest point in ps to X
                        if iseven(indx)
                            Rx_ = [ps[indx-1][1], ps[indx][1]]
                            Zx_ = [ps[indx-1][2], ps[indx][2]]
                        else
                            Rx_ = [ps[indx][1], ps[indx+1][1]]
                            Zx_ = [ps[indx][2], ps[indx+1][2]]
                        end
                    end
                end

                # save strike-points in counter-clockwise order. Note: from here on, length(Rx_) = 2
                angle = mod.(atan.(Zx_ .- zaxis, Rx_ .- raxis), 2 * π) # counter-clockwise angle form geom axis

                append!(Rxx, Rx_[sortperm(angle)])
                append!(Zxx, Zx_[sortperm(angle)])
                append!(dxx, dx_ .* ones(length(Rx_)))
            end
        end
    end
    indexx = sortperm(dxx) # save in order by dxx
    return (Rxx=Rxx[indexx], Zxx=Zxx[indexx], dxx=dxx[indexx])
end

@compat public find_strike_points
push!(document[Symbol("Physics sol")], :find_strike_points)

"""
    find_strike_points!(
        eqt::IMAS.equilibrium__time_slice{T1},
        wall_r::AbstractVector{T2},
        wall_z::AbstractVector{T2},
        psi_last_closed::Real,
        psi_first_open::Real,
        dv::IMAS.divertors{T1};
        private_flux_regions::Bool=true,
        in_place::Bool=true
    ) where {T1<:Real,T2<:Real}

Adds strike points location to equilibrium IDS
"""
function find_strike_points!(
    eqt::IMAS.equilibrium__time_slice{T1},
    wall_r::AbstractVector{T2},
    wall_z::AbstractVector{T2},
    psi_last_closed::Real,
    psi_first_open::Real,
    dv::IMAS.divertors{T1};
    private_flux_regions::Bool=true,
    in_place::Bool=true
) where {T1<:Real,T2<:Real}

    Rxx = Float64[]
    Zxx = Float64[]
    dxx = Float64[]

    for divertor in dv.divertor
        for target in divertor.target
            Rxx0, Zxx0, dxx0 = find_strike_points(
                eqt,
                wall_r,
                wall_z,
                psi_last_closed,
                psi_first_open;
                private_flux_regions,
                strike_surfaces_r=target.tile[1].surface_outline.r,
                strike_surfaces_z=target.tile[1].surface_outline.z
            )
            # allow for strike points to miss the divertors
            if isempty(Rxx0)
                continue
            end
            push!(Rxx, Rxx0[1])
            push!(Zxx, Zxx0[1])
            push!(dxx, dxx0[1])
        end
    end

    if in_place
        resize!(eqt.boundary.strike_point, length(Rxx))
        for (k, strike_point) in enumerate(eqt.boundary.strike_point)
            strike_point.r = Rxx[k]
            strike_point.z = Zxx[k]
            strike_point.last_closed_flux_surface_gap = dxx[k]
        end
    end

    return (Rxx=Rxx, Zxx=Zxx, dxx=dxx)
end

"""
    find_strike_points(
        eqt::IMAS.equilibrium__time_slice{T1},
        wall_r::AbstractVector{T2},
        wall_z::AbstractVector{T2},
        psi_last_closed::Real,
        psi_first_open::Real,
        dv::IMAS.divertors{T1};
        private_flux_regions::Bool=true
    ) where {T1<:Real,T2<:Real}

Return strike points location in the divertors
"""
function find_strike_points(
    eqt::IMAS.equilibrium__time_slice{T1},
    wall_r::AbstractVector{T2},
    wall_z::AbstractVector{T2},
    psi_last_closed::Real,
    psi_first_open::Real,
    dv::IMAS.divertors{T1};
    private_flux_regions::Bool=true
) where {T1<:Real,T2<:Real}

    return find_strike_points!(eqt, wall_r, wall_z, psi_last_closed, psi_first_open, dv; private_flux_regions, in_place=false)
end

"""
    find_strike_points!(
        eqt::IMAS.equilibrium__time_slice{T1},
        wall_r::AbstractVector{T2},
        wall_z::AbstractVector{T2},
        psi_last_closed::Real,
        psi_first_open::Real
    ) where {T1<:Real,T2<:Real}
"""
function find_strike_points!(
    eqt::IMAS.equilibrium__time_slice{T1},
    wall_r::AbstractVector{T2},
    wall_z::AbstractVector{T2},
    psi_last_closed::Real,
    psi_first_open::Real
) where {T1<:Real,T2<:Real}

    if psi_last_closed == psi_first_open
        Rxx=Float64[]
        Zxx=Float64[]
        dxx=Float64[]

    else
        Rxx, Zxx, dxx = find_strike_points(eqt, wall_r, wall_z, psi_last_closed, psi_first_open)
        resize!(eqt.boundary.strike_point, length(Rxx))
        for (k, strike_point) in enumerate(eqt.boundary.strike_point)
            strike_point.r = Rxx[k]
            strike_point.z = Zxx[k]
            strike_point.last_closed_flux_surface_gap = dxx[k]
        end
    end

    return (Rxx=Rxx, Zxx=Zxx, dxx=dxx)
end

@compat public find_strike_points!
push!(document[Symbol("Physics sol")], :find_strike_points!)