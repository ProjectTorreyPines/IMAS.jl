document[Symbol("Physics thermal loads")] = Symbol[]

"""
    r::Vector{Float64}                  # R coordinate of the wall mesh                     - m
    z::Vector{Float64}                  # Z coordinate of the wall mesh                     - m
    q_wall::Vector{Float64}             # total heat flux on the wall                       - W/m2
    q_part::Vector{Float64}             # heat flux on the wall due to particles            - W/m2
    q_core_rad::Vector{Float64}         # heat flux on the wall due to core radiation       - W/m2
    q_parallel::Vector{Float64}         # parallel heat flux due to particles at the wall   - W/m2
    s::Vector{Float64}                  # wall curvilinear abscissa
"""
Base.@kwdef mutable struct WallHeatFlux{T}
    r::Vector{Float64}=Float64[]
    z::Vector{Float64}=Float64[]
    q_wall::Vector{T}=T[]
    q_part::Vector{T}=T[]
    q_core_rad::Vector{T}=T[]
    q_parallel::Vector{T}=T[]
    s::Vector{Float64}=Float64[]
end

"""
    particle_heat_flux(
        eqt::IMAS.equilibrium__time_slice,
        SOL::OrderedCollections.OrderedDict{Symbol,Vector{OpenFieldLine}},
        wall_r::AbstractVector{<:Real},
        wall_z::AbstractVector{<:Real},
        r::Vector{<:Real},
        q::Vector{<:Real};
        merge_wall::Bool=true)

Computes the heat flux on the wall due to the influx of charged particles, using the magnetic equilibrium,
the Scrape Off-Layer, the wall, and an hypothesis of the decay of the parallel heat flux at the OMP
"""
function particle_heat_flux(
    eqt::IMAS.equilibrium__time_slice,
    SOL::OrderedCollections.OrderedDict{Symbol,Vector{OpenFieldLine}},
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    r::Vector{<:Real},
    q::Vector{<:Real};
    merge_wall::Bool=true)

    @assert length(r) == length(q)
    @assert all(q .>= 0) # q is all positive
    @assert all(r .>= 0) # r is all positive  
    @assert r[1] == 0

    Rwall = Float64[]
    Zwall = Float64[]
    Qwall = Float64[]
    Qpara = Float64[]
    indexes = Int64[]

    if isempty(wall_r) || isempty(wall_z)
        error("Impossible to map the heat flux onto the wall because dd.wall is empty")
    end

    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    psi_separatrix = find_psi_boundary(eqt, wall_r, wall_z; raise_error_on_not_open=true).first_open # psi at LCFS
    surface = flux_surface(eqt, psi_separatrix, :open, wall_r, wall_z)
    r_separatrix = Float64[]
    for (rr, zz) in surface
        if isempty(rr) || all(zz .> ZA) || all(zz .< ZA)
            continue
        end

        crossings = intersection(rr, zz, [RA, 10 * RA], [ZA, ZA]).crossings # find intersection with midplane

        if isempty(crossings)
            continue
        end

        push!(r_separatrix, crossings[1][1])
    end
    r_separatrix = r_separatrix[1]
    r = r .+ r_separatrix
    @assert maximum(r) > maximum(wall_r)

    q_interp = interp1d(r, q, :cubic)

    if eqt.boundary.x_point[1].z < ZA
        case = :lower
    else
        case = :upper
    end
    if isempty(SOL[:lfs_far])
        case = :double
    end

    # lower single null case
    if case == :lower || case == :upper
        #order clockwise starting from midplane
        for sol in reverse(SOL[:lfs_far])
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(Qwall, qmid / (sol.total_flux_expansion[end]) * sin(sol.grazing_angles[end]))
            push!(Qpara, qmid / (sol.total_flux_expansion[end]))
            push!(indexes, sol.wall_index[end])
        end
        for sol in reverse(SOL[:lfs])
            # Outer target
            qmid = q_interp(sol.r[sol.midplane_index])
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(Qwall, qmid / (sol.total_flux_expansion[end]) * sin(sol.grazing_angles[end]))
            push!(Qpara, qmid / (sol.total_flux_expansion[end]))
            push!(indexes, sol.wall_index[end])
        end
        for sol in SOL[:lfs]
            # Inner target
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(Qwall, qmid / (sol.total_flux_expansion[1]) * sin(sol.grazing_angles[1]))
            push!(Qpara, qmid / (sol.total_flux_expansion[1]))
            push!(indexes, sol.wall_index[1])
        end

        for sol in SOL[:lfs_far]
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(Qwall, qmid / (sol.total_flux_expansion[1]) * sin(sol.grazing_angles[1]))
            push!(Qpara, qmid / (sol.total_flux_expansion[1]))
            push!(indexes, sol.wall_index[1])
        end

        Rwall_hfs = Float64[]
        Zwall_hfs = Float64[]
        Qwall_hfs = Float64[]
        Qpara_hfs = Float64[]
        indexes_hfs = Int64[]

        if !isempty(SOL[:hfs])
            for sol in SOL[:hfs]
                #order clockwise starting from midplane
                push!(Rwall_hfs, sol.r[1]) # R of points after midplane
                push!(Zwall_hfs, sol.z[1]) # z of points after midplane
                push!(Qwall_hfs, 0.0)
                push!(Qpara_hfs, 0.0)
                push!(indexes_hfs, sol.wall_index[1])
            end
            for sol in reverse(SOL[:hfs])
                push!(Rwall_hfs, sol.r[end]) # R of points after midplane
                push!(Zwall_hfs, sol.z[end]) # z of points after midplane
                push!(Qwall_hfs, 0.0)
                push!(Qpara_hfs, 0.0)
                push!(indexes_hfs, sol.wall_index[end])
            end

            # insert hfs
            Rwall = vcat(Rwall[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                Rwall_hfs,
                Rwall[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
            Zwall = vcat(Zwall[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                Zwall_hfs,
                Zwall[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
            Qwall = vcat(Qwall[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                Qwall_hfs,
                Qwall[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
            Qpara = vcat(Qpara[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                Qpara_hfs,
                Qpara[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
            indexes = vcat(indexes[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                indexes_hfs,
                indexes[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
        end
    end

    # double null case - SOL[:lfs_far] is empty
    if case == :double
        #order clockwise starting from midplane
        for sol in reverse(SOL[:lfs])
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(Qwall, qmid / (sol.total_flux_expansion[end]) * sin(sol.grazing_angles[end]))
            push!(Qpara, qmid / (sol.total_flux_expansion[end]))
            push!(indexes, sol.wall_index[end])
        end
        for sol in SOL[:hfs]
            #order clockwise starting from midplane
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(Qwall, 0)
            push!(Qpara, 0)
            push!(indexes, sol.wall_index[1])
        end
        for sol in reverse(SOL[:hfs])
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(Qwall, 0)
            push!(Qpara, 0)
            push!(indexes, sol.wall_index[end])
        end
        for sol in SOL[:lfs]
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(Qwall, qmid / (sol.total_flux_expansion[1]) * sin(sol.grazing_angles[1]))
            push!(Qpara, qmid / (sol.total_flux_expansion[1]))
            push!(indexes, sol.wall_index[1])
        end
    end

    if merge_wall
        _, _, PSI_interpolant = ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)

        #! format: off
        if wall_z[1]>ZA # correction if first point is above midplane 
            indexes[indexes .== 1 .&& Zwall.>ZA] .= length(wall_r) # put it after index = end  
            indexes[indexes .== 1 .&& Zwall.>wall_z[2] .&& Zwall.<=ZA] .= 0   # put before index = 1
        end
        #! format: on

        crossings = intersection([(minimum(wall_r) + maximum(wall_r)) / 2, maximum(wall_r) * 1.05], [ZA, ZA], wall_r, wall_z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
        r_wall_midplane = crossings[1][1] # R coordinate of the wall at OMP
        psi_wall_midplane = PSI_interpolant(r_wall_midplane, ZA)
        _, psi_first_lfs_far, null_within_wall = find_psi_last_diverted(eqt, wall_r, wall_z, PSI_interpolant) # psi of grazing surface
        psi_wall = PSI_interpolant.(wall_r, wall_z)
        tollZ = max(abs(ZA), 2 * abs(wall_z[2] - wall_z[1]))
        if Zwall[1] >= -tollZ # check if the particle reach the wall around the midplane
            # if so, do not add any point
            add_omp = false
        else
            # add points
            add_omp = true
        end

        #! format: off
        add_indexes = collect((1:length(psi_wall)))
        add_indexes = add_indexes[(psi_wall .< psi_separatrix .|| psi_wall .> psi_wall_midplane) .|| # add private flux region  around first null (psi<psi_sep) + add everyhting above psi midplane
                                (psi_wall .>= psi_separatrix .&& psi_wall .<= psi_first_lfs_far .&&  # add also points inside the private region around second null (only if null is within wall)
                                sign(eqt.boundary.x_point[end].z).*wall_z.>abs(eqt.boundary.x_point[end].z)).&& null_within_wall .||
                                psi_wall .> psi_first_lfs_far .&& (wall_r.<eqt.boundary.x_point[end].r) .|| # add point in :hfs 
                                (psi_wall .< psi_wall_midplane .&& (wall_z.<=tollZ) .&& (wall_z .>= -tollZ) .&& (wall_r.>eqt.boundary.x_point[end].r) .&& add_omp) # add also points close to OMP but with psi lower than psi_wall midplane within tollZ from omp
                                ] 
        L = length(add_indexes)
        #! format: on

        average_step = sum(sqrt.(diff(wall_r) .^ 2 + diff(wall_z) .^ 2)) ./ (length(wall_r) - 1)
        #filter out point that are for some reason already inside Rwall,Zwall
        for (k, ind) in enumerate(reverse(add_indexes))
            #check if these points have been already saved in Rwall, Zwall
            dist = sqrt.((Rwall .- wall_r[ind]) .^ 2 + (Zwall .- wall_z[ind]) .^ 2)
            if minimum(dist) <= average_step # if a point is at less than average_step  from an already saved point, do not save it
                #point is already in and must be removed
                deleteat!(add_indexes, L + 1 - k)
            end
        end

        add_indexes = reverse!(add_indexes)

        for ind in add_indexes
            if ind != 1
                Rwall = append!(Rwall[indexes.<ind], [wall_r[ind]], Rwall[indexes.>=ind])
                Zwall = append!(Zwall[indexes.<ind], [wall_z[ind]], Zwall[indexes.>=ind])
                Qwall = append!(Qwall[indexes.<ind], [0.0], Qwall[indexes.>=ind])
                Qpara = append!(Qpara[indexes.<ind], [0.0], Qpara[indexes.>=ind])
                indexes = append!(indexes[indexes.<ind], [ind], indexes[indexes.>=ind])
            end
        end

        # include points shadowed by corners
        dist = (diff(Rwall)) .^ 2 + (diff(Zwall)) .^ 2
        avg = sum(dist) / length(dist) # mean distance between adjacent points in mesh
        std = sqrt(sum((dist .- avg) .^ 2) / (length(dist) - 1)) #standard deviation of distance adjacent points in mesh
        save_length = length(Rwall)
        counter = 0
        while sum(dist .> avg + 4 * std) > 0 && counter < 10
            # we have outliers, which are shadowed areas.
            # in the outliears dist > average + 4 standard deviations
            indexes = 1:length(dist)
            indexes = indexes[dist.>avg+4*std]

            for ind in reverse(indexes)
                # if vertical lines
                if abs.(diff([Rwall[ind], Rwall[ind+1]]))[1] > 0.01
                    dr = 0.0
                else
                    # dr = 1.0*maximum(sqrt.(diff(wall_r).^2 + diff(wall_z).^2))
                    dr = sqrt(2) * sum(sqrt.(diff(wall_r) .^ 2 + diff(wall_z) .^ 2)) / length(dist)
                end
                #if horizontal lines
                if abs.(diff([Zwall[ind], Zwall[ind+1]]))[1] > 0.01
                    dz = 0.0
                else
                    # dz = 1.0*maximum(sqrt.(diff(wall_r).^2 + diff(wall_z).^2))
                    dz = sqrt(2) * sum(sqrt.(diff(wall_r) .^ 2 + diff(wall_z) .^ 2)) / length(dist)
                end

                #search points in rectangle between two points
                #! format: off
                add_r = wall_r[wall_r.>(minimum([Rwall[ind],Rwall[ind+1]])-dr ).&& wall_r .< (maximum([Rwall[ind],Rwall[ind+1]])+dr ).&&
                               wall_z.>(minimum([Zwall[ind],Zwall[ind+1]])-dz) .&& wall_z .< (maximum([Zwall[ind],Zwall[ind+1]])+dz )  ]
                add_z = wall_z[wall_r.>(minimum([Rwall[ind],Rwall[ind+1]])-dr) .&& wall_r .< (maximum([Rwall[ind],Rwall[ind+1]])+dr ).&&
                               wall_z.>(minimum([Zwall[ind],Zwall[ind+1]])-dz) .&& wall_z .< (maximum([Zwall[ind],Zwall[ind+1]])+dz)   ]
                #! format: on
                Rwall = append!(Rwall[1:ind], add_r, Rwall[ind+1:end])
                Zwall = append!(Zwall[1:ind], add_z, Zwall[ind+1:end])
                Qwall = append!(Qwall[1:ind], add_r .* 0.0, Qwall[ind+1:end])
                Qpara = append!(Qpara[1:ind], add_r .* 0.0, Qpara[ind+1:end])
            end
            if length(Rwall) == save_length
                # I am not finding anymore points to add
                break
            end
            dist = (diff(Rwall)) .^ 2 + (diff(Zwall)) .^ 2
            avg = sum(dist) / length(dist) # mean distance between adjacent points in mesh
            std = sqrt(sum((dist .- avg) .^ 2) / (length(dist) - 1)) #standard deviation of distance adjacent points in mesh
            counter += 1
            save_length = length(Rwall)
        end
        # add OMP point
        # NOTE: this treatment of the OMP_wall point must be revised (physics missing)
        if (r_wall_midplane, ZA) in (Rwall, Zwall)
        else
            Rwall = vcat(r_wall_midplane, Rwall)
            Zwall = vcat(ZA, Zwall)
            Qwall = vcat(0.0, Qwall) # this must be modified according to physics
            Qpara = vcat(q_interp(r_wall_midplane), Qpara)
        end
    end

    # Compute the heat flux due to the core radiation - Qrad - W/m2
    # Make sure (Rwall,Zwall) is closed
    if !((Rwall[1], Zwall[1]) == (Rwall[end], Zwall[end]))
        # Mesh not closed!
        push!(Rwall, Rwall[1])
        push!(Zwall, Zwall[1])
        push!(Qwall, Qwall[1])
        push!(Qpara, Qpara[1])
    end

    return (Qwall=Qwall, Qpara=Qpara)
end

@compat public particle_heat_flux
push!(document[Symbol("Physics thermal loads")], :particle_heat_flux)

"""
    core_radiation_heat_flux(
        eqt::IMAS.equilibrium__time_slice, 
        psi::Vector{<:Real}, 
        source_1d::Vector{<:Real}, 
        N::Int,
        wall_r::AbstractVector{<:Real}, 
        wall_z::AbstractVector{<:Real},
        Prad_core::Float64)
"""
function core_radiation_heat_flux(
    eqt::IMAS.equilibrium__time_slice,
    psi::Vector{<:Real},
    source_1d::Vector{<:Real},
    N::Int,
    wall_r::AbstractVector{<:Real},
    wall_z::AbstractVector{<:Real},
    Prad_core::Float64)

    photons, W_per_trace, dr, dz = define_particles(eqt, psi, source_1d, N)

    qflux_r, qflux_z, wall_s = find_flux(photons, W_per_trace, wall_r, wall_z, dr, dz)

    qq = sqrt.(qflux_r .^ 2 + qflux_z .^ 2) # norm of the heat flux
    power = qq .* wall_s

    # normalization to match perfectly the power in core_sources
    norm = Prad_core / sum(power)
    qq *= norm

    # qq is defined in the midpoints of the grid (wall_r, wall_z), which is the center of the cells where the heat flux is computed
    # Interpolate the values on the nodes (wall_r, wall_z)
    Qrad = similar(wall_r)
    Qrad[1] = (qq[1] + qq[end]) / 2

    if ((wall_r[1], wall_z[1]) == (wall_r[end], wall_z[end]))
        # Wall is closed, therefore length(wall_s) = length(wall_r) - 1
        # length(q) = length(wall_s)

        # qq is defined in the midpoint between nodes in the wall mesh (ss vector)
        # Qrad shall be instead defined on the Rwall,Zwall mesh (s vector)

        s = similar(wall_r)
        ds = sqrt.(diff(wall_r) .^ 2 + diff(wall_z) .^ 2)
        s[1] = 0
        for i in 1:length(ds)
            s[i+1] = s[i] + ds[i]
        end

        ss = (s[2:end] + s[1:end-1]) ./ 2

        ss = vcat(0.0, ss)
        qq = vcat((qq[1] + qq[end]) / 2, qq)

        interp = interp1d(vcat(ss .- s[end], ss, ss .+ s[end]), vcat(qq, qq, qq), :cubic)
        Qrad = interp.(s)

    else
        # Wall is NOT closed, therefore length(wall_s) = length(wall_r)
        error("wall is not closed")
    end

    return Qrad
end

"""
    mesher_heat_flux(dd::IMAS.dd; 
        r::AbstractVector{T}=Float64[], 
        q::AbstractVector{T}=Float64[], 
        merge_wall::Bool = true, 
        levels::Union{Int,AbstractVector} = 20, 
        step::T = 0.1) where {T<:Real}

Computes the wall mesh for the heat flux deposited on the wall. Returns:

(Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
(rwall, zwall)        evenly spaced wall mesh - m
s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
SOL                   list of OpenFieldLines used to compute (Rwall, Zwall)
(r,q)                 Hypothesis of power density decay at omp for definition of SOL
"""
function mesher_heat_flux(dd::IMAS.dd;
    r::AbstractVector{T}=Float64[],
    q::AbstractVector{T}=Float64[],
    merge_wall::Bool=true,
    levels::Union{Int,AbstractVector}=20,
    step::T=0.1) where {T<:Real}

    eqt = dd.equilibrium.time_slice[]
    fw = first_wall(dd.wall)

    if isempty(fw.r) || isempty(fw.z)
        error("Impossible to map the heat flux onto the wall because dd.wall is empty")
    end

    R0 = eqt.global_quantities.magnetic_axis.r # R magentic axis 
    Z0 = eqt.global_quantities.magnetic_axis.z # Z magnetic axis
    psi_separatrix = find_psi_boundary(eqt, fw.r,fw.z; raise_error_on_not_open=true).first_open # psi at LCFS

    if isempty(r) || isempty(q)
        ##########################################################################
        ####### NOTE: do it better with functions lambdaq 1 e lambdaq 2 and P_elms
        ##########################################################################

        # define q(r), but could be arbitrary. Here double exponential
        a = dd.equilibrium.time_slice[].boundary.minor_radius      # minor radius
        l = widthSOL_eich(dd) # decay length of first exponential
        l2 = 0.1  # decay length of second exponential - NOTE put scaling second lambda q - Loarte 2008
        frac = 0.2 # fraction of power flowing in the second esponential 
        NN = 2 # imbalance factor (in/out)  1 < N < 2
        Bt_omp = dd.equilibrium.time_slice[].global_quantities.vacuum_toroidal_field.b0 * R0 / (R0 + a)

        crossings = intersection([R0, 2 * maximum(fw.r)], [Z0, Z0], fw.r, fw.z).crossings # (r,z) point of intersection btw outer midplane (OMP) with wall
        r_wall_omp = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
        r_wall_omp = r_wall_omp[1] # make it float

        surface = flux_surface(eqt, psi_separatrix, :encircling, fw.r, fw.z)
        (rr, zz) = surface[1]
        crossings = intersection([R0, 2 * maximum(fw.r)], [Z0, Z0], rr, zz).crossings
        r_sep = [cr[1] for cr in crossings] # R coordinate of points in SOL surface at MP (inner and outer)
        r_sep = r_sep[1]

        r = collect(LinRange(0, r_wall_omp - r_sep + 2 * l2, 10000)) # Ideally: from 0 to r max grid - r_separatrix_OMP
        # double exponential
        q =
            power_sol(dd) / 2 / NN / π / (R0 + a) / l / sin(atan(Bpol_omp(dd.equilibrium.time_slice[]) / abs(Bt_omp))) * exp.(-r ./ l) +
            power_sol(dd) * frac / 2 / NN / π / (R0 + a) / l2 / sin(atan(Bpol_omp(dd.equilibrium.time_slice[]) / abs(Bt_omp))) * exp.(-r ./ l2)
    end

    step = minimum([step, sum(sqrt.(diff(fw.r) .^ 2 + diff(fw.z) .^ 2)) / 250]) # ensure decent resolution of the wall
    # resample wall and make sure it's clockwise (for COCOS = 11)
    rwall, zwall = resample_2d_path(fw.r, fw.z; step, method=:linear, retain_original_xy=true)
    reorder_flux_surface!(rwall, zwall, R0, Z0; force_close=true)

    # Parameters for particle heat flux
    if typeof(levels) <: Int
        #levels is an Int, build a vector of psi_levels of that size
        eqt2d = findfirst(:rectangular, eqt.profiles_2d)
        _, _, PSI_interpolant = ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)
        psi_levels, _, _ = find_levels_from_P(eqt, rwall, zwall, PSI_interpolant, r, q, levels)
        add_psi = find_levels_from_wall(eqt, rwall, zwall, PSI_interpolant)
        psi_sign = sign(psi_levels[end] - psi_levels[1])
        psi_levels = unique!(sort!(vcat(psi_levels, add_psi)))


        if psi_sign == -1
            psi_levels = reverse!(psi_levels) # if psi is decreasing, sort in descending order
        end
    else
        #levels is a vector, use it as it is
        psi_levels = levels
    end

    psi_levels[1] = psi_separatrix
    # build SOL
    SOL = sol(eqt, rwall, zwall; levels=psi_levels, use_wall=true)

    Rwall = Float64[]
    Zwall = Float64[]
    indexes = Int64[]

    if eqt.boundary.x_point[1].z < Z0
        case = :lower
    else
        case = :upper
    end
    if isempty(SOL[:lfs_far])
        case = :double
    end

    # lower single null case
    if case == :lower || case == :upper
        #order clockwise starting from midplane
        for sol in reverse(SOL[:lfs_far])
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(indexes, sol.wall_index[end])
        end
        for sol in reverse(SOL[:lfs])
            # Outer target
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(indexes, sol.wall_index[end])
        end
        for sol in SOL[:lfs]
            # Inner target
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(indexes, sol.wall_index[1])
        end

        for sol in SOL[:lfs_far]
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(indexes, sol.wall_index[1])
        end

        Rwall_hfs = Float64[]
        Zwall_hfs = Float64[]
        indexes_hfs = Int64[]

        if !isempty(SOL[:hfs])
            for sol in SOL[:hfs]
                #order clockwise starting from midplane
                push!(Rwall_hfs, sol.r[1]) # R of points after midplane
                push!(Zwall_hfs, sol.z[1]) # z of points after midplane
                push!(indexes_hfs, sol.wall_index[1])
            end
            for sol in reverse(SOL[:hfs])
                push!(Rwall_hfs, sol.r[end]) # R of points after midplane
                push!(Zwall_hfs, sol.z[end]) # z of points after midplane
                push!(indexes_hfs, sol.wall_index[end])
            end

            # insert hfs
            Rwall = vcat(
                Rwall[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                Rwall_hfs,
                Rwall[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
            Zwall = vcat(
                Zwall[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                Zwall_hfs,
                Zwall[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
            indexes = vcat(
                indexes[1:argmin(abs.(indexes .- maximum(indexes_hfs)))-1],
                indexes_hfs,
                indexes[argmin(abs.(indexes .- maximum(indexes_hfs))):end])
        end
    end

    # double null case - SOL[:lfs_far] is empty
    if case == :double
        #order clockwise starting from midplane
        for sol in reverse(SOL[:lfs])
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(indexes, sol.wall_index[end])
        end
        for sol in SOL[:hfs]
            #order clockwise starting from midplane
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(indexes, sol.wall_index[1])
        end
        for sol in reverse(SOL[:hfs])
            push!(Rwall, sol.r[end]) # R of points after midplane
            push!(Zwall, sol.z[end]) # z of points after midplane
            push!(indexes, sol.wall_index[end])
        end
        for sol in SOL[:lfs]
            push!(Rwall, sol.r[1]) # R of points after midplane
            push!(Zwall, sol.z[1]) # z of points after midplane
            push!(indexes, sol.wall_index[1])
        end
    end

    if merge_wall
        _, _, PSI_interpolant = ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)

        #! format: off
        if zwall[1] > Z0 # correction if first point is above midplane 
            indexes[indexes .== 1 .&& Zwall.>Z0] .= length(rwall) # put it after index = end  
            indexes[indexes .== 1 .&& Zwall.>zwall[2] .&& Zwall.<=Z0] .= 0   # put before index = 1
        end
        #! format: on

        crossings = intersection([(minimum(rwall) + maximum(rwall)) / 2, maximum(rwall) * 1.05], [Z0, Z0], rwall, zwall)[2] # (r,z) point of intersection btw outer midplane (OMP) with wall
        r_wall_midplane = crossings[1][1] # R coordinate of the wall at OMP
        psi_wall_midplane = PSI_interpolant(r_wall_midplane, Z0)
        _, psi_first_lfs_far, null_within_wall = find_psi_last_diverted(eqt, rwall, zwall, PSI_interpolant) # psi of grazing surface
        psi_wall = PSI_interpolant.(rwall, zwall)
        tollZ = max(abs(Z0), 2 * abs(zwall[2] - zwall[1]))
        if Zwall[1] >= -tollZ # check if the particle reach the wall around the midplane
            # if so, do not add any point
            add_omp = false
        else
            # add points
            add_omp = true
        end

        add_indexes = collect((1:length(psi_wall)))
        #! format: off
        add_indexes = add_indexes[(psi_wall .< psi_separatrix .|| psi_wall .> psi_wall_midplane) .|| # add private flux region  around first null (psi<psi_sep) + add everyhting above psi midplane
                                (psi_wall .>= psi_separatrix .&& psi_wall .<= psi_first_lfs_far .&&  # add also points inside the private region around second null (only if null is within wall)
                                sign(eqt.boundary.x_point[end].z).*zwall.>abs(eqt.boundary.x_point[end].z)) .&& null_within_wall .||
                                psi_wall .> psi_first_lfs_far .&& (rwall.<eqt.boundary.x_point[end].r) .|| # add point in :hfs 
                                (psi_wall .< psi_wall_midplane .&& (zwall.<=tollZ) .&& (zwall .>= -tollZ) .&& (rwall.>eqt.boundary.x_point[end].r) .&& add_omp) # add also points close to OMP but with psi lower than psi_wall midplane within tollZ from omp
                                ]
        #! format: on

        L = length(add_indexes)
        average_step = sum(sqrt.(diff(rwall) .^ 2 + diff(zwall) .^ 2)) ./ (length(rwall) - 1)
        #filter out point that are for some reason already inside Rwall,Zwall
        for (k, ind) in enumerate(reverse(add_indexes))
            #check if these points have been already saved in Rwall, Zwall
            dist = sqrt.((Rwall .- rwall[ind]) .^ 2 + (Zwall .- zwall[ind]) .^ 2)
            if minimum(dist) <= average_step # if a point is at less than average_step  from an already saved point, do not save it
                #point is already in and must be removed
                deleteat!(add_indexes, L + 1 - k)
            end
        end

        add_indexes = reverse!(add_indexes)
        for ind in add_indexes
            if ind != 1
                Rwall = append!(Rwall[indexes.<ind], [rwall[ind]], Rwall[indexes.>=ind])
                Zwall = append!(Zwall[indexes.<ind], [zwall[ind]], Zwall[indexes.>=ind])
                indexes = append!(indexes[indexes.<ind], [ind], indexes[indexes.>=ind])
            end
        end

        # include points shadowed by corners
        dist = (diff(Rwall)) .^ 2 + (diff(Zwall)) .^ 2
        avg = sum(dist) / length(dist) # mean distance between adjacent points in mesh
        std = sqrt(sum((dist .- avg) .^ 2) / (length(dist) - 1)) #standard deviation of distance adjacent points in mesh
        save_length = length(Rwall)
        counter = 0
        while sum(dist .> avg + 4 * std) > 0 && counter < 10
            # we have outliers, which are shadowed areas.
            # in the outliears dist > average + 4 standard deviations
            indexes = 1:length(dist)
            indexes = indexes[dist.>avg+4*std]

            for ind in reverse(indexes)
                # if vertical lines
                if abs.(diff([Rwall[ind], Rwall[ind+1]]))[1] > 0.01
                    dr = 0.0
                else
                    # dr = 1.0*maximum(sqrt.(diff(rwall).^2 + diff(zwall).^2))
                    dr = sqrt(2) * sum(sqrt.(diff(rwall) .^ 2 + diff(zwall) .^ 2)) / length(dist)
                end
                #if horizontal lines
                if abs.(diff([Zwall[ind], Zwall[ind+1]]))[1] > 0.01
                    dz = 0.0
                else
                    # dz = 1.0*maximum(sqrt.(diff(rwall).^2 + diff(zwall).^2))
                    dz = sqrt(2) * sum(sqrt.(diff(rwall) .^ 2 + diff(zwall) .^ 2)) / length(dist)
                end

                #search points in rectangle between two points
                #! format: off
                add_r = rwall[rwall.>(minimum([Rwall[ind],Rwall[ind+1]])-dr ).&& rwall .< (maximum([Rwall[ind],Rwall[ind+1]])+dr ).&&
                              zwall.>(minimum([Zwall[ind],Zwall[ind+1]])-dz) .&& zwall .< (maximum([Zwall[ind],Zwall[ind+1]])+dz )  ]
                add_z = zwall[rwall.>(minimum([Rwall[ind],Rwall[ind+1]])-dr) .&& rwall .< (maximum([Rwall[ind],Rwall[ind+1]])+dr ).&&
                              zwall.>(minimum([Zwall[ind],Zwall[ind+1]])-dz) .&& zwall .< (maximum([Zwall[ind],Zwall[ind+1]])+dz)   ]
                #! format: on
                
                # check that (add_r,add_z) are unique points
                add_rz = unique!(collect(zip(add_r,add_z)))
                if length(add_r) > length(add_rz)
                    add_r  = [rz[1] for rz in add_rz]
                    add_z  = [rz[2] for rz in add_rz]
                end

                Rwall = append!(Rwall[1:ind], add_r, Rwall[ind+1:end])
                Zwall = append!(Zwall[1:ind], add_z, Zwall[ind+1:end])
            end
            if length(Rwall) == save_length
                # I am not finding anymore points to add
                break
            end
            dist = (diff(Rwall)) .^ 2 + (diff(Zwall)) .^ 2
            avg = sum(dist) / length(dist) # mean distance between adjacent points in mesh
            std = sqrt(sum((dist .- avg) .^ 2) / (length(dist) - 1)) #standard deviation of distance adjacent points in mesh
            counter += 1
            save_length = length(Rwall)
        end
        # add OMP point
        if (r_wall_midplane, Z0) in (Rwall, Zwall)
        else
            Rwall = vcat(r_wall_midplane, Rwall)
            Zwall = vcat(Z0, Zwall)
        end
    end

    s = similar(Rwall)
    ds = sqrt.(diff(Rwall) .^ 2 + diff(Zwall) .^ 2)
    s[1] = 0.0
    for i in 1:length(ds)
        s[i+1] = s[i] + ds[i]
    end

    # Make sure (Rwall,Zwall) is closed
    if !((Rwall[1], Zwall[1]) == (Rwall[end], Zwall[end]))
        # Mesh not closed!
        push!(Rwall, Rwall[1])
        push!(Zwall, Zwall[1])
        push!(s, s[end] + sqrt((Rwall[1] - Rwall[end])^2 + (Zwall[1] - Zwall[end])^2))
    end

    return (Rwall=Rwall, Zwall=Zwall, rwall=rwall, zwall=zwall, s=s, SOL=SOL, r=r, q=q)
end

#= ============= =#
#  plotting       #
#= ============= =#
"""
Recipe for plot of heat flux
  - which_plot = :twoD, :oneD
  - plot_type  = :path, :scatter (only for 2D)
  - q          =
    :all (for 1D),
    :wall,
    :parallel
    :particle
    :core_radiation
    :both (= :particle + :parallel)
"""
@recipe function plot_heat_flux(HF::WallHeatFlux; which_plot=:twoD, plot_type=:path, q=:wall)
    @assert which_plot in (:oneD, :twoD)
    @assert plot_type in (:path, :scatter)
    @assert q in (:wall, :parallel, :particle, :core_radiation, :both, :all)

    if q in (:both, :all)
        if q == :both
            qs = (:parallel, :particle)
        else
            qs = (:wall, :particle, :core_radiation)
        end
        if which_plot == :twoD
            layout := (1, length(qs))
        end
        for (k, q) in enumerate(qs)
            @series begin
                if which_plot == :twoD
                    subplot := k
                end
                which_plot := which_plot
                plot_type := plot_type
                q := q
                HF
            end
        end

    elseif which_plot == :oneD
        @series begin
            xlabel --> "Clockwise distance along wall [m]"
            ylabel --> "Wall flux [W/m²]"
            yscale --> :log10

            x = []
            y = []
            if q == :particle
                label --> "particles"
                if !isempty(HF.q_part)
                    x = HF.s
                    y = HF.q_part .+ 1.0
                end

            elseif q == :core_radiation
                label --> "core radiation"
                if !isempty(HF.q_core_rad)
                    x = HF.s
                    y = HF.q_core_rad .+ 1.0
                end

            elseif q == :wall
                label --> "wall"
                if !isempty(HF.q_wall)
                    x = HF.s
                    y = HF.q_wall .+ 1.0
                end

            elseif q == :parallel
                label --> "parallel"
                if !isempty(HF.q_parallel)
                    x = HF.s
                    y = HF.q_parallel .+ 1.0
                end
            end

            x, y
        end

    elseif which_plot == :twoD
        @series begin
            aspect_ratio := :equal
            legend --> false
            colorbar --> true
            xlim --> [0.95 * minimum(HF.r) - 0.05 * (maximum(HF.r)), 1.05 * maximum(HF.r) - 0.05 * (minimum(HF.r))]
            ylim --> [1.05 * minimum(HF.z) - 0.05 * (maximum(HF.z)), 1.05 * maximum(HF.z) - 0.05 * (minimum(HF.z))]
            if plot_type == :path
                seriestype --> :path
                linewidth --> 8
            elseif plot_type == :scatter
                seriestype --> :scatter
                markersize --> 2
                markerstrokewidth --> 0
            end

            if q == :wall && !isempty(HF.q_wall)
                colorbar_title := "log₁₀(q wall [W/m²])"
                if plot_type == :path
                    line_z := log10.(HF.q_wall .+ 1)
                elseif plot_type == :scatter
                    zcolor := log10.(HF.q_wall .+ 1)
                end

            elseif q == :core_radiation && !isempty(HF.q_core_rad)
                colorbar_title := "log₁₀(q core rad [W/m²])"
                if plot_type == :path
                    line_z := log10.(HF.q_core_rad .+ 1)
                elseif plot_type == :scatter
                    zcolor := log10.(HF.q_core_rad .+ 1)
                end

            elseif q == :particle && !isempty(HF.q_part)
                colorbar_title := "log₁₀(q particle [W/m²])"
                floor = minimum(HF.q_part[HF.q_part.>0.0]) / 100.0
                if plot_type == :path
                    line_z --> log10.(HF.q_part .+ floor)
                elseif plot_type == :scatter
                    zcolor --> log10.(HF.q_part .+ floor)
                end

            elseif q == :parallel && !isempty(HF.q_parallel)
                colorbar_title := "log₁₀(q parallel [W/m²])"
                if plot_type == :path
                    line_z --> log10.(HF.q_parallel .+ 1)
                end
                if plot_type == :scatter
                    zcolor --> log10.(HF.q_parallel .+ 1)
                end
            end

            HF.r, HF.z
        end
    end
end
