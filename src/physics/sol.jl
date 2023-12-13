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
end

@recipe function plot_ofl(ofl::OpenFieldLine)
    @series begin
        aspect_ratio --> :equal
        label --> ""
        colorbar_title := "log₁₀(Connection length [m] + 1.0)"
        line_z := log10.(ofl.s .+ 1)
        ofl.r, ofl.z
    end
end

@recipe function plot_OFL(OFL_hfs_lfs_lfsfar::OrderedCollections.OrderedDict{Symbol,Vector{IMAS.OpenFieldLine}})
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

@recipe function plot_OFL(OFL::Vector{OpenFieldLine})
    for ofl in OFL
        @series begin
            ofl
        end
    end
end

"""
    sol(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{T}, wall_z::Vector{T}; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true) where {T<:Real}

Returns vectors of hfs and lfs OpenFieldLine

If levels is a vector, it has the values of psi from 0 to max psi_wall_midplane. The function will modify levels of psi to introduce relevant sol surfaces
"""
function sol(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{T}, wall_z::Vector{T}; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true) where {T<:Real}
    ############ 
    R0, B0 = vacuum_r0_b0(eqt)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    ############
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)
    psi__axis_level = eqt.profiles_1d.psi[1] # psi value on axis 
    psi__boundary_level = find_psi_boundary(eqt; raise_error_on_not_open=true) # find psi at LCFS
    # find psi at second magnetic separatrix 
    psi__2nd_separatix = find_psi_2nd_separatrix(eqt, PSI_interpolant) # find psi at 2nd magnetic separatrix
    psi_sign = sign(psi__boundary_level - psi__axis_level) # sign of the poloidal flux taking psi_axis = 0
    if use_wall
        # SOL with wall
        crossings = intersection([RA, maximum(wall_r)], [ZA, ZA], wall_r, wall_z)[2] # (r,z) point of intersection btw outer midplane (OMP) with wall
        r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
        psi_wall_midplane = PSI_interpolant.(r_wall_midplane, ZA)[1] # psi at the intersection between wall and omp
        psi_last_diverted, null_is_inside = find_psi_last_diverted(eqt, wall_r, wall_z, PSI_interpolant) # find psi at LDFS
    else
        # SOL without wall
        psi_wall_midplane = maximum(psi_sign .* eqt2d.psi) - psi_sign # if no wall, upper bound of psi is maximum value in eqt -1 (safe)
        r_wall_midplane = eqt2d.grid.dim1[end] # if no wall, take max R in psi grid
        psi_last_diverted = [0,1].*0.00001.* abs(psi__boundary_level)
        null_is_inside = true
    end
    ############

    # smart picking of psi levels
    if levels == 1
        levels = [psi__boundary_level + psi_sign .* 1E-3]

    elseif levels == 2
        levels = [psi__boundary_level + psi_sign .* 1E-3, psi__2nd_separatix]

    elseif typeof(levels) <: Int
        levels = psi__boundary_level .+ psi_sign .* 10.0 .^ LinRange(-3, log10(abs(psi_wall_midplane - psi_sign * 0.001 * abs(psi_wall_midplane) - psi__boundary_level)), levels)
        if null_is_inside
            levels[argmin(abs.(levels .- psi__2nd_separatix))] = psi__2nd_separatix # make sure 2nd separatrix is in levels
        else
            indexx = argmin(abs.(levels .- psi_last_diverted[1]))
            levels = vcat(levels[1:indexx-1], psi_last_diverted, levels[indexx+1:end]) # remove closest point + add LDFS (it is a vector)
            levels = sort(vcat(levels, psi__2nd_separatix))
        end

    else
        #levels is a vector of psi_levels for the discretization of the SOL
        @assert levels[1] >= psi__boundary_level
        @assert levels[end] <= psi_wall_midplane
        levels_is_not_monotonic_in_Ip_direction = all(psi_sign * diff(levels) .>= 0)
        @assert levels_is_not_monotonic_in_Ip_direction # levels must be monotonic according to plasma current direction
        # make sure levels includes separatrix and wall
        levels[1] = psi__boundary_level + psi_sign * abs.(diff(psi_last_diverted))[1] # if psi = psi__boundary_level, flux_surface does not work
        levels[end] = psi_wall_midplane - psi_sign * 0.001 * abs(psi_wall_midplane)
    end

    OFL_hfs = OpenFieldLine[] # field lines magnetically isolated from OMP
    OFL_lfs = OpenFieldLine[] # field lines magnetically connected to OMP inside  last diverted flux surface
    OFL_lfs_far = OpenFieldLine[] # field lines magnetically connected to OMP outside last diverted flux surface
    # TO DO for the future: insert private flux regions (upper and lower)

    # r_mid(ψ) interpolator for region of interest
    r_mid_of_interest = 10.0 .^ range(log10(maximum(eqt.boundary.outline.r) * 0.99), log10(maximum(r_wall_midplane)), 1000)
    r_mid_itp = interp_rmid_at_psi(PSI_interpolant, r_mid_of_interest, ZA)

    for level in levels
        lines, _ = flux_surface(eqt, level, :open) #returns (r,z) of surfaces with psi = level
        for (r, z) in lines
            if use_wall
                # SOL with wall
                # returns poloidal angles of each surface
                # rr and zz are clockwise
                # crossing points with wall = (rr[1], zz[1]) (rr[end], zz[end])
                # this is the order at which angles are computed (strike, pitch and grazing)
                # Example - OFL[2].[1<n<length(levels)].strike_angle[1] is computed at (rr[1], zz[1])
                rr, zz, strike_angles, wall_index = line_wall_2_wall(r, z, wall_r, wall_z, RA, ZA)
            else
                # SOL without wall
                rr = r
                zz = z
                strike_angles = [NaN, NaN]
                wall_index    = [0, 0]
            end
            if isempty(rr) || all(zz .> ZA) || all(zz .< ZA)
                continue
            end

            # add a point exactly at the (preferably outer) midplane
            crossing_index, crossings = intersection([0, maximum(wall_r) * 1.5], [ZA, ZA], rr, zz)
            r_midplane = [cr[1] for cr in crossings] # R coordinate of points in SOL surface at MP (inner and outer)
            z_midplane = [cr[2] for cr in crossings] # Z coordinate of points in SOL surface at MP (inner and outer)
            outer_index = argmax(r_midplane)  #index of point @ MP: this is OMP (for OFL_lfs); IMP for OFL_hfs
            crossing_index = crossing_index[outer_index] #indexes of point at MP in SOL surface
            r_midplane = r_midplane[outer_index] # R coordinate of point at MP in SOL surface
            z_midplane = z_midplane[outer_index] # Z coordinate of point at MP in SOL surface
            rr = [rr[1:crossing_index[2]]; r_midplane; rr[crossing_index[2]+1:end]] #Insert in r of SOL surface a point @ MP
            zz = [zz[1:crossing_index[2]]; z_midplane; zz[crossing_index[2]+1:end]] #Insert in z of SOL surface a point @ MP
            midplane_index = crossing_index[2] + 1 #index at which point @ MP

            # calculate quantities along field line
            Br, Bz = Br_Bz(PSI_interpolant, rr, zz) #r and z component of B for each point in (r,z)
            Bp = sqrt.(Br .^ 2.0 .+ Bz .^ 2.0)     #poloidal component of B for each point in (r,z)
            Bt = abs.(B0 .* R0 ./ rr)              #toroidal component of B for each point in (r,z)
            B = sqrt.(Bp .^ 2 + Bt .^ 2)                 # total magnetic field B for each point in (r,z)
            dp = sqrt.(gradient(rr) .^ 2.0 .+ gradient(zz) .^ 2.0) # curvilinear abscissa increments of poloidal projection of SOL surface
            pitch = sqrt.(1.0 .+ (Bt ./ Bp) .^ 2) # ds = dp*sqrt(1 + (Bt/Bp)^2) (pythagora)
            s = cumsum(pitch .* dp) # s = integral(ds)
            s = abs.(s .- s[midplane_index]) # fix 0 at outer midplane

            # Parameters to map heat flux from OMP to wall
            pitch_angles = atan.(Bp, Bt)
            grazing_angles = asin.(sin.([pitch_angles[1], pitch_angles[end]]) .* sin.(strike_angles))
            total_flux_expansion = B[midplane_index] ./ B # total flux expansion(r,z) =  Bomp / B(r,z) [magentic flux conservation]
            poloidal_flux_expansion = total_flux_expansion .* rr[midplane_index] ./ rr .* sin(pitch_angles[midplane_index]) ./ sin.(pitch_angles) # poloidal flux expansion

            # select HFS or LFS and add line to the list 
            if rr[midplane_index] < RA
                # Add SOL surface in OFL_hfs
                OFL = OFL_hfs # surfaces magnetically isolated from OMP
            else
                # update R coordinate of point at OMP in SOL surface, such that PSI_interpolant(rr[midplane_index],ZA) == level
                rr[midplane_index] = r_mid_itp(level)
                if level <= sum(psi_last_diverted)/length(psi_last_diverted)
                    # Add SOL surface in OFL_lfs
                    OFL = OFL_lfs
                else
                    # Add SOL surface in OFL_lfs_far
                    OFL = OFL_lfs_far
                end
            end
            push!(OFL, OpenFieldLine(rr, zz, Br, Bz, Bp, Bt, pitch, s, midplane_index, strike_angles, pitch_angles, grazing_angles, total_flux_expansion, poloidal_flux_expansion,wall_index)) # add result
        end
    end

    # In the case of a perfectly balanced double null OFL_lfs will be empty and all of the flux surfaces will appear in OFL_lfs_far
    # Here we switch the two
    if isempty(OFL_lfs)
        tmp = OFL_lfs_far
        OFL_lfs_far = OFL_lfs
        OFL_lfs = tmp
    end

    return OrderedCollections.OrderedDict(:hfs => OFL_hfs, :lfs => OFL_lfs, :lfs_far => OFL_lfs_far)
end

function sol(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true)
    return sol(eqt, first_wall(wall).r, first_wall(wall).z; levels, use_wall)
end

function sol(dd::IMAS.dd; levels::Union{Int,AbstractVector}=20, use_wall::Bool=true)
    return sol(dd.equilibrium.time_slice[], dd.wall; levels, use_wall)
end


"""

function find_levels_from_P(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{<:Real}, wall_z::Vector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{<:Real}, q::Vector{<:Real}, levels::Int) 

    Function for the discretization of the poloidal flux ψ on the SOL, based on an hypotesis of OMP radial transport through arbitrary q(r)
    returns vector with level of ψ, vector with matching r_midplane and q.
    Discretization with even steps of P = integral_sep^wal 2πrq(r)dr (same power in each flux tube)
"""
function find_levels_from_P(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{<:Real}, wall_z::Vector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{<:Real}, q::Vector{<:Real}, levels::Int) 
    ################### Housekeeping on function q(r) ###################
    @assert length(r) == length(q)
    @assert all(q .>=0) # q is all positive
    @assert all(r .>=0) # r is all positive  
    @assert r[1] == 0
    
    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis
    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    # r_mid(ψ) interpolator for region of interest
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    rmax = eqt2d.grid.dim1[end]
    r_mid_of_interest = 10.0 .^ range(log10(maximum(eqt.boundary.outline.r) * 0.99), log10(rmax), 1000)
    r_mid = interp_rmid_at_psi(PSI_interpolant, r_mid_of_interest, ZA)

    psi__boundary_level   = find_psi_boundary(eqt; raise_error_on_not_open=true) # psi at LCFS
    r_separatrix_midplane = r_mid(psi__boundary_level)      # R OMP at separatrix 
    crossings = intersection([RA, maximum(wall_r)], [ZA, ZA], wall_r, wall_z)[2] # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
    r_wall_midplane = r_wall_midplane[1] # make it float
    psi_wall_midplane = PSI_interpolant(r_wall_midplane,ZA)[1]
    psi_2ndseparatrix = find_psi_2nd_separatrix(eqt,PSI_interpolant) # psi of the second magnetic separatrix
    r_2ndseparatrix_midplane = r_mid(psi_2ndseparatrix) # R coordinate at OMP of 2nd magnetic separatrix
    psi_last_diverted, null_is_inside = find_psi_last_diverted(eqt,wall_r,wall_z,PSI_interpolant) # psi of grazing surface
    r_last_diverted = r_mid.(psi_last_diverted) # R coordinate at OMP of grazing surface

    r = r .+ r_separatrix_midplane
    order = sortperm(r)
    q = q[order] 
    r = r[order] # r is now increasing monotonically

    if r[1] >= r_wall_midplane
        # all the vector r is inside the wall
        r = [r_separatrix_midplane, r_wall_midplane] # constant value of q as closest point to SOL
        q = [1,0.999]*q[1]  # keep preference for monotonic decrease

    end
    if r[end] <= r_separatrix_midplane
        # all of r is inside the separatrix
        r = [r_separatrix_midplane, r_wall_midplane] # constant value of q as closest point to SOL
        q = [1.001,1]*q[end] # keep preference for monotonic decrease
    end

    # r[1] can be either >, = or < than r_separatrix_midplane; add r_separatrix_midplane
    if r[1] > r_separatrix_midplane 
        # r starts from inside the sol
        r = append!([r_separatrix_midplane], r) # add a point at R_OMP
        q = append!([q[1]*1.001], q)            # Repeat first value of q with slight increment, favoring monotonic decrease
    end
    # if r[1] == r_separatrix_midplane do nothing
    if r[1] < r_separatrix_midplane
        # r starts from inside the separatrix, q(r) must be cut
        index = argmin(abs.(r .- r_separatrix_midplane)) # closest point
        # index2 is the position in r, such that r_separatrix_midplane is between r[index2] and r[index]
        if r[index] > r_separatrix_midplane
            index2 = index - 1 
            #interp linearly value at r_separatrix_midplane between r[index2] and r[index]
            qq = q[index]+ (q[index2]-q[index])/(r[index2]-r[index])*(r_separatrix_midplane - r[index]) 
            r = append!([r_separatrix_midplane], r[index:end]) # cut r and q
            q = append!([qq]                   , q[index:end])
        else
            index2 = index + 1
            #interp linearly value at r_separatrix_midplane between r[index2] and r[index]
            qq = q[index]+ (q[index2]-q[index])/(r[index2]-r[index])*(r_separatrix_midplane - r[index]) 
            r = append!([r_separatrix_midplane], r[index2:end]) # cut r and q
            q = append!([qq]                   , q[index2:end])
        end
    end

    # r[end] can be either >, = < than r_wall_midplane; add r_wall_midplane
    if r[end] < r_wall_midplane
        # r ends inside sol
        r = push!(r, r_wall_midplane) # add a point at r_wall_midplane
        q = push!(q, q[end]*0.999)    # Repeat last value of q with slight reduction, favoring monotonic decrease
    end
    # if r[end]==r_wall_midplane do nothing
    if r[end] > r_wall_midplane
        # r ends inside the wall, q(r) must be cut
        index = argmin(abs.(r .- r_wall_midplane)) # closest point
        if r[index] > r_wall_midplane 
            index2 = index-1
            #interp linearly value at r_wall_midplane between r[index2] and r[index]
            qq = q[index]+ (q[index2]-q[index])/(r[index2]-r[index])*(r_wall_midplane - r[index]) 
            r = push!(r[1:index-1],r_wall_midplane) # cut + add point between index and index2
            q = push!(q[1:index-1],qq)
        else
            index2 = index+1
            #interp linearly value at r_wall_midplane between r[index2] and r[index]
            qq = q[index]+ (q[index2]-q[index])/(r[index2]-r[index])*(r_wall_midplane - r[index]) 
            r = push!(r[1:index],r_wall_midplane) # cut + add a point between index and index 
            q = push!(q[1:index],qq)
        end
    end

    ############### finished houskeeping of q(r) ####################
    #################################################################

    # build P(r) = integral_sep^r q(ρ)2πρdρ
    P = q*0;
    for index in 2:length(q)
        P[index] = P[index-1] + integrate(r[index-1:index], 2*π.*r[index-1:index].*q[index-1:index])
    end
    # being 2πr q(r) positive-definite, P(r) is strictly monotonic, therefore also injective (one-to-one)
    # P(r) is always invertible for every q(r)>0
    r = Interpolations.deduplicate_knots!(r) 
    interp_P = interp1d(r, P, :cubic) # interpolant of P(r)

    p_levels = collect(LinRange(0,maximum(P),levels))   # levels to interpolate P
    # add flux surfaces of interest: last diverted surface and 2nd magnetic separatrix
    # NOTE: number of level increases
    P_low = interp_P(r_last_diverted[1]) # last diverted surface (up to precision); inside OFL[:lfs]
    P_up  = interp_P(r_last_diverted[2]) # first surface crossing the top first wall; inside OFL[:lfs_far]
    if null_is_inside
        # last diverted surface is the 2nd separatrix; add only 2nd separatrix_up and 2nd separatix_low
        if !(P_low in p_levels)
            p_levels = push!(p_levels,P_low)
        end
        if !(P_up in p_levels)
            p_levels = push!(p_levels,P_up)
        end
        # TO BE CHECKED: in case of double null, power balance must work
    else
        # last diverted surface is different from the 2nd separatix: add 2nd magnetic separatix and last diverted surface
        P_2ndseparatrix = interp_P(r_2ndseparatrix_midplane) # interp value at 2nd separatrix
        # add all surfaces, but avoid repetition
        if !(P_2ndseparatrix in p_levels)
            p_levels = push!(p_levels,P_2ndseparatrix)
        end
        if !(P_low in p_levels)
            p_levels = push!(p_levels,P_low)
        end
        if !(P_up in p_levels)
            p_levels = push!(p_levels,P_up)
        end
    end 

    p_levels = sort!(p_levels)
    P = Interpolations.deduplicate_knots!(P) 
    interp_inverseP = interp1d(P, r, :cubic) # interpolant of inverse function of r(P)
    R = interp_inverseP.(p_levels)

    # using ψ(R), go from discretization in R to discretization in ψ
    psi_levels = PSI_interpolant.(R, R.*0.0 .+ZA) # ψ(R)

    #force psi_sep and psi_wall_midplane
    psi_levels[1] = psi__boundary_level
    psi_levels[end] = psi_wall_midplane
    return psi_levels, R, p_levels
    
end

function find_levels_from_P(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall, PSI_interpolant::Interpolations.AbstractInterpolation, r::Vector{<:Real}, q::Vector{<:Real}, levels::Int) 
    return find_levels_from_P(eqt, first_wall(wall).r,first_wall(wall).z, PSI_interpolant, q, r, levels)
end

function find_levels_from_P(dd::IMAS.dd, r::Vector{<:Real}, q::Vector{<:Real}, levels::Int) 
    rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
    return find_levels_from_P(dd.equilibrium.time_slice[], dd.wall, PSI_interpolant, q, r, levels) 
end


"""

function find_levels_from_wall(wall_r::Vector{<:Real}, wall_z::Vector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation) 

    Function for that computes the value of psi at the points of the wall mesh in dd
"""
function find_levels_from_wall(eqt::IMAS.equilibrium__time_slice, wall_r::Vector{<:Real}, wall_z::Vector{<:Real}, PSI_interpolant::Interpolations.AbstractInterpolation) 
    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    psi_separatrix = find_psi_boundary(eqt; raise_error_on_not_open=true) #psi on separatrix
    crossings = intersection([(minimum(wall_r)+maximum(wall_r))/2, maximum(wall_r)*1.05], [ZA, ZA], wall_r, wall_z)[2] # (r,z) point of intersection btw outer midplane (OMP) with wall
    r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
    r_wall_midplane = r_wall_midplane[1];
    psi_wall_midplane = PSI_interpolant(r_wall_midplane,ZA); 

    levels =  PSI_interpolant.(wall_r,wall_z)
    levels = levels[levels .>= psi_separatrix .&& levels.<= psi_wall_midplane]
    return sort!(levels)
end

function find_levels_from_wall(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall,PSI_interpolant::Interpolations.AbstractInterpolation) 
    return find_levels_from_wall(eqt,first_wall(wall).r,first_wall(wall).z, PSI_interpolant)
end

function find_levels_from_wall(dd::IMAS.dd) 
    rr, zz, PSI_interpolant = ψ_interpolant(dd.equilibrium.time_slice[].profiles_2d)
return find_levels_from_wall(dd.equilibrium.time_slice[], dd.wall, PSI_interpolant)
end

"""
    line_wall_2_wall(r::T, z::T, wall_r::T, wall_z::T, RA::Real, ZA::Real) where {T<:AbstractVector{<:Real}}

Returns r, z coordinates of open field line contained within wall, as well as angles of incidence at the strike locations

RA and ZA are the coordinate of the magnetic axis
"""
function line_wall_2_wall(r::T, z::T, wall_r::T, wall_z::T, RA::Real, ZA::Real) where {T<:AbstractVector{<:Real}}
    indexes, crossings = intersection(r, z, wall_r, wall_z) # find where flux surface crosses wall ("strike points" of surface)
    # crossings -  Vector{Tuple{Float64, Float64}} - crossings[1] contains (r,z) of first "strike point"
    # indexes   -  Vector{Tuple{Float64, Float64}} - indexes[1] contains indexes of (r,z) and (wall_r, wall_z) of first "strike point"
    r_z_index = [k[1] for k in indexes] #index of vectors (r,z) of all crossing point
    wall_index = [k[2] for k in indexes] #index of vectors (wall_r, wall_z) of all crossing point

    if length(r_z_index) == 0  # if the flux surface does not cross the wall return empty vector (it is not a surf in SOL)
        return Float64[], Float64[], Float64[], Int64[]

    elseif length(r_z_index) == 1
        error("line_wall_2_wall: open field line should intersect wall at least twice.
            If it does not it's likely because the equilibrium grid was too small.
            Suggestion: plot dd.wall + eqt.profiles_2d to debug.")
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
        i = sort([i1, i2])
        r_z_index = r_z_index[i]
        wall_index = wall_index[i]
        crossings = crossings[i]
        strike_angles = strike_angles[i]
    end

    rr = vcat(crossings[1][1], r[r_z_index[1]+1:r_z_index[2]], crossings[2][1]) # r coordinate of magnetic surface between one "strike point" and the other
    zz = vcat(crossings[1][2], z[r_z_index[1]+1:r_z_index[2]], crossings[2][2]) # z coordinate of magnetic surface between one "strike point" and the other
    if sum(rr .< minimum(wall_r))>0 
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

"""
    divertor_totals_from_targets(divertor::IMAS.divertors__divertor, field::Symbol)

Returns time dependent vectors of :field summed over all divertor targets
"""
function divertor_totals_from_targets(divertor::IMAS.divertors__divertor, field::Symbol)
    total = []
    time = []
    for target in divertor.target
        value = getproperty(target, field)
        if !ismissing(value, :data)
            push!(total, value.data)
            push!(time, value.time)
        end
    end
    return time[1], reduce(+, total)
end

"""
    Bpol(a::T, κ::T, Ip::T) where {T<:Real}

Average poloidal magnetic field magnitude
"""
function Bpol(a::T, κ::T, Ip::T) where {T<:Real}
    return (constants.μ_0 * Ip) / (2π * a * sqrt((1.0 + κ^2) / 2.0))
end

"""
    Bpol_omp(eqt::IMAS.equilibrium__time_slice)

Poloidal magnetic field magnitude evaluated at the outer midplane
"""
function Bpol_omp(eqt::IMAS.equilibrium__time_slice)
    r, z, PSI_interpolant = ψ_interpolant(eqt.profiles_2d)
    eq1d = eqt.profiles_1d
    R_omp = eq1d.r_outboard[end]
    Z_omp = eqt.global_quantities.magnetic_axis.z
    return Bp(PSI_interpolant, [R_omp], [Z_omp])[1]
end

"""
    power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d)

Total power coming out of the SOL [W]

NOTE: This function returns 1.0 [W] if power is less than that so that SOL quantities remain finite
"""
function power_sol(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d)
    Psol = total_power_source(total_sources(core_sources, cp1d; fields=[:power_inside, :total_ion_power_inside]))
    if Psol < 1.0
        return one(Psol)
    else
        return Psol
    end
end

function power_sol(dd::IMAS.dd)
    return power_sol(dd.core_sources, dd.core_profiles.profiles_1d[])
end

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

function widthSOL_loarte(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    B0 = B0_geo(eqt)
    q95 = eqt.global_quantities.q_95
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_loarte(B0, q95, Psol)
end

function widthSOL_loarte(dd::IMAS.dd)
    return widthSOL_loarte(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

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

function widthSOL_sieglin(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eq1d = eqt.profiles_1d
    R0 = (eq1d.r_outboard[end] .+ eq1d.r_inboard[end]) / 2.0
    a = (eq1d.r_outboard[end] .- eq1d.r_inboard[end])
    Psol = power_sol(core_sources, cp1d)
    ne_ped = interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(0.95)
    return widthSOL_sieglin(R0, a, Bpol_omp(eqt), Psol, ne_ped)
end

function widthSOL_sieglin(dd::IMAS.dd)
    return widthSOL_sieglin(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

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

function widthSOL_eich(eqt::IMAS.equilibrium__time_slice, Psol::Real)
    eq1d = eqt.profiles_1d
    R0 = (eq1d.r_outboard[end] .+ eq1d.r_inboard[end]) / 2.0
    a = (eq1d.r_outboard[end] .- eq1d.r_inboard[end])
    return widthSOL_eich(R0, a, Bpol_omp(eqt), Psol)
end

function widthSOL_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    Psol = power_sol(core_sources, cp1d)
    return widthSOL_eich(eqt, Psol)
end

function widthSOL_eich(dd::IMAS.dd)
    return widthSOL_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    q_pol_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)

Poloidal heat flux [W/m^2] at the outer midplane based on Eigh λ_q
"""
function q_pol_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eq1d = eqt.profiles_1d
    R_omp = eq1d.r_outboard[end]
    Psol = power_sol(core_sources, cp1d)
    channel_area = 2π * R_omp * widthSOL_eich(eqt, cp1d, core_sources)
    return Psol / channel_area
end

function q_pol_omp_eich(dd::IMAS.dd)
    return q_pol_omp_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

"""
    q_par_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)

Parallel heat flux [W/m^2] at the outer midplane based on Eigh λ_q
"""
function q_par_omp_eich(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, core_sources::IMAS.core_sources)
    eq1d = eqt.profiles_1d
    R0, B0 = vacuum_r0_b0(eqt)
    R_omp = eq1d.r_outboard[end]
    Bt_omp = B0 * R0 / R_omp
    return q_pol_omp_eich(eqt, cp1d, core_sources) / sin(atan(Bpol_omp(eqt) / Bt_omp))
end

function q_par_omp_eich(dd::IMAS.dd)
    return q_par_omp_eich(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[], dd.core_sources)
end

# ==== #
# Zhom #
# ==== #
"""
zohm_divertor_figure_of_merit(eqt::IMAS.equilibrium__time_slice)

Computes a figure of merit for the divertor (Zohm) PB/R/q/A [W T/m]
"""
function zohm_divertor_figure_of_merit(core_sources::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice, T::summary__global_quantities)
    R0 = eqt.boundary.geometric_axis.r
    a = eqt.boundary.minor_radius
    A = R0 / a
    q95 = eqt.global_quantities.q_95
    B0 = @ddtime(T.b0.value)
    Psol = power_sol(core_sources, cp1d)

    zohm = Psol * B0 / R0 / A / q95 # W T/m
    return zohm
end

function zohm_divertor_figure_of_merit(dd::IMAS.dd)
    return zohm_divertor_figure_of_merit(dd.core_sources, dd.core_profiles.profiles_1d[], dd.equilibrium.time_slice[], dd.summary.global_quantities)
end

# ============= #
# Strike points #
# ============= #
"""
    find_strike_points(wall_outline_r::T, wall_outline_z::T, pr::T, pz::T) where {T<:AbstractVector{<:Real}}

Finds strike points and angles of incidence between two paths
"""
function find_strike_points(wall_outline_r::T, wall_outline_z::T, pr::T, pz::T) where {T<:AbstractVector{<:Real}}
    indexes, crossings = intersection(wall_outline_r, wall_outline_z, pr, pz)
    pvx = [cr[1] for cr in crossings]
    pvy = [cr[2] for cr in crossings]
    angles = intersection_angles(wall_outline_r, wall_outline_z, pr, pz, indexes)
    return pvx, pvy, angles
end

"""
    find_strike_points(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}

Finds equilibrium strike points and angle of incidence between wall and strike leg
"""
function find_strike_points(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}
    Rx = Float64[]
    Zx = Float64[]
    θx = Float64[]

    private, _ = flux_surface(eqt, eqt.profiles_1d.psi[end], :open)
    for (pr, pz) in private
        pvx, pvy, angles = find_strike_points(wall_outline_r, wall_outline_z, pr, pz)
        append!(Rx, pvx)
        append!(Zx, pvy)
        append!(θx, angles)
    end

    return Rx, Zx, θx
end

"""
    find_strike_points!(eqt::IMAS.equilibrium__time_slice, dv::IMAS.divertors)

Adds strike points location to equilibrium IDS and the tilt_angle_pol in the divertors IDS
"""
function find_strike_points!(eqt::IMAS.equilibrium__time_slice, dv::IMAS.divertors)
    Rx = Float64[]
    Zx = Float64[]
    θx = Float64[]

    time = eqt.time

    for divertor in dv.divertor
        for target in divertor.target
            Rx0, Zx0, θx0 = find_strike_points(eqt, target.tile[1].surface_outline.r, target.tile[1].surface_outline.z)
            # allow for strike points to miss the divertors
            if isempty(Rx0)
                continue
            end
            push!(Rx, Rx0[1])
            push!(Zx, Zx0[1])
            push!(θx, θx0[1])
            set_time_array(target.tilt_angle_pol, :data, time, θx0[1])
        end
    end

    resize!(eqt.boundary_separatrix.strike_point, length(Rx))
    for (k, strike_point) in enumerate(eqt.boundary_separatrix.strike_point)
        strike_point.r = Rx[k]
        strike_point.z = Zx[k]
    end

    return Rx, Zx, θx
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall_outline_r::T, wall_outline_z::T) where {T<:AbstractVector{<:Real}}
    Rx, Zx, θx = find_strike_points(eqt, wall_outline_r, wall_outline_z)

    resize!(eqt.boundary_separatrix.strike_point, length(Rx))
    for (k, strike_point) in enumerate(eqt.boundary_separatrix.strike_point)
        strike_point.r = Rx[k]
        strike_point.z = Zx[k]
    end

    return Rx, Zx, θx
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, bd::IMAS.build)
    wall_outline = get_build_layer(bd.layer; type=_plasma_).outline
    return find_strike_points!(eqt, wall_outline.r, wall_outline.z)
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice, wall::IMAS.wall)
    wall_outline = first_wall(wall)
    if wall_outline !== missing
        return find_strike_points!(eqt, wall_outline.r, wall_outline.z)
    end
end

function find_strike_points!(eqt::IMAS.equilibrium__time_slice)
    dd = top_dd(eqt)
    wall_outline = first_wall(dd.wall)
    if wall_outline !== missing
        return find_strike_points!(eqt, wall_outline.r, wall_outline.z)
    elseif !isempty(dd.build.layer)
        return find_strike_points!(eqt, dd.build)
    end
end
