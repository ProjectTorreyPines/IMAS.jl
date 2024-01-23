mutable struct WallHeatFlux
    r::Vector{Float64}
    z::Vector{Float64}
    q_wall::Vector{Float64}
    q_parallel::Vector{Float64}
    s::Vector{Float64}
end

function WallHFMapper(eqt::IMAS.equilibrium__time_slice, 
                      SOL::OrderedCollections.OrderedDict{Symbol, Vector{IMAS.OpenFieldLine}}, 
                      wall_r::Vector{<:Real}, 
                      wall_z::Vector{<:Real}, 
                      r::Vector{<:Real}, 
                      q::Vector{<:Real}; 
                      merge_wall::Bool = true)

    @assert length(r) == length(q)
    @assert all(q .>=0) # q is all positive
    @assert all(r .>=0) # r is all positive  
    @assert r[1] == 0
                      
    Rwall =  Float64[]
    Zwall =  Float64[]
    Qwall =  Float64[]
    Qpara =  Float64[]
    indexes =  Int64[] 

    if isempty(wall_r) || isempty(wall_z)
        error("Impossible to map the heat flux onto the wall because dd.wall is empty")
    end

    ZA = eqt.global_quantities.magnetic_axis.z # Z of magnetic axis
    RA = eqt.global_quantities.magnetic_axis.r # R of magnetic axis

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)
    _, psi_separatrix   = IMAS.find_psi_boundary(eqt; raise_error_on_not_open=true) # psi at LCFS
    surface, _ =  IMAS.flux_surface(eqt, psi_separatrix, :open)
    r_separatrix =  Float64[]
    for (rr,zz) in surface
        if isempty(rr) || all(zz .> ZA) || all(zz .< ZA) 
            continue
        end

        _, crossings = IMAS.intersection(rr, zz, [1, 10]*RA, [1, 1]*ZA) # find intersection with midplane

        if isempty(crossings)
            continue
        end
        
        rsep = [cr[1] for cr in crossings]
        
        push!(r_separatrix,rsep[1])

    end
    r_separatrix =  r_separatrix[1]
    r = r .+ r_separatrix
    q_interp = IMAS.interp1d(r, q, :cubic) 

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
            push!(Rwall,sol.r[end]) # R of points after midplane
            push!(Zwall,sol.z[end]) # z of points after midplane
            push!(Qwall, qmid/(sol.total_flux_expansion[end])*sin(sol.grazing_angles[end]))
            push!(Qpara, qmid/(sol.total_flux_expansion[end]))
            push!(indexes,sol.wall_index[end])
        end
        for  sol in reverse(SOL[:lfs])
            # Outer target
            qmid = q_interp(sol.r[sol.midplane_index])
            push!(Rwall,sol.r[end]) # R of points after midplane
            push!(Zwall,sol.z[end]) # z of points after midplane
            push!(Qwall,qmid/(sol.total_flux_expansion[end])*sin(sol.grazing_angles[end]))
            push!(Qpara, qmid/(sol.total_flux_expansion[end]))
            push!(indexes,sol.wall_index[end])
        end
        for  sol in SOL[:lfs]
            # Inner target
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall,sol.r[1]) # R of points after midplane
            push!(Zwall,sol.z[1]) # z of points after midplane
            push!(Qwall, qmid/(sol.total_flux_expansion[1])*sin(sol.grazing_angles[1]))
            push!(Qpara, qmid/(sol.total_flux_expansion[1]))
            push!(indexes,sol.wall_index[1])
        end

        for sol in SOL[:lfs_far]
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall,sol.r[1]) # R of points after midplane
            push!(Zwall,sol.z[1]) # z of points after midplane
            push!(Qwall, qmid/(sol.total_flux_expansion[1])*sin(sol.grazing_angles[1]))
            push!(Qpara, qmid/(sol.total_flux_expansion[1]))
            push!(indexes,sol.wall_index[1])
        end

        Rwall_hfs =  Float64[]
        Zwall_hfs =  Float64[]
        Qwall_hfs =  Float64[]
        Qpara_hfs =  Float64[]
        indexes_hfs =  Int64[] 

        if !isempty(SOL[:hfs])
            for sol in SOL[:hfs]
                #order clockwise starting from midplane
                push!(Rwall_hfs,sol.r[1]) # R of points after midplane
                push!(Zwall_hfs,sol.z[1]) # z of points after midplane
                push!(Qwall_hfs,0.0)
                push!(Qpara_hfs,0.0)
                push!(indexes_hfs,sol.wall_index[1])
            end
            for sol in reverse(SOL[:hfs])
                push!(Rwall_hfs,sol.r[end]) # R of points after midplane
                push!(Zwall_hfs,sol.z[end]) # z of points after midplane
                push!(Qwall_hfs,0.0)
                push!(Qpara_hfs,0.0)
                push!(indexes_hfs,sol.wall_index[end])
            end

            
            # insert hfs
            Rwall   = vcat(Rwall[1:argmin(abs.(indexes.-maximum(indexes_hfs)))-1], 
                        Rwall_hfs, 
                        Rwall[argmin(abs.(indexes.-maximum(indexes_hfs))):end])
            Zwall   = vcat(Zwall[1:argmin(abs.(indexes.-maximum(indexes_hfs)))-1], 
                        Zwall_hfs, 
                        Zwall[argmin(abs.(indexes.-maximum(indexes_hfs))):end])
            Qwall   = vcat(Qwall[1:argmin(abs.(indexes.-maximum(indexes_hfs)))-1], 
                        Qwall_hfs, 
                        Qwall[argmin(abs.(indexes.-maximum(indexes_hfs))):end])
            Qpara   = vcat(Qpara[1:argmin(abs.(indexes.-maximum(indexes_hfs)))-1], 
                        Qpara_hfs, 
                        Qpara[argmin(abs.(indexes.-maximum(indexes_hfs))):end])
            indexes = vcat(indexes[1:argmin(abs.(indexes.-maximum(indexes_hfs)))-1], 
                        indexes_hfs, 
                        indexes[argmin(abs.(indexes.-maximum(indexes_hfs))):end])
        end
    end

    # double null case - SOL[:lfs_far] is empty
    if case == :double
        #order clockwise starting from midplane
        for sol in reverse(SOL[:lfs])
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall,sol.r[end]) # R of points after midplane
            push!(Zwall,sol.z[end]) # z of points after midplane
            push!(Qwall, qmid/(sol.total_flux_expansion[end])*sin(sol.grazing_angles[end]))
            push!(Qpara, qmid/(sol.total_flux_expansion[end]))
            push!(indexes,sol.wall_index[end])
        end
        for sol in SOL[:hfs]
            #order clockwise starting from midplane
            push!(Rwall,sol.r[1]) # R of points after midplane
            push!(Zwall,sol.z[1]) # z of points after midplane
            push!(Qwall,0)
            push!(Qpara,0)
            push!(indexes,sol.wall_index[1])
        end
        for sol in reverse(SOL[:hfs])
            push!(Rwall,sol.r[end]) # R of points after midplane
            push!(Zwall,sol.z[end]) # z of points after midplane
            push!(Qwall,0)
            push!(Qpara,0)
            push!(indexes,sol.wall_index[end])
        end
        for sol in SOL[:lfs]
            qmid = q_interp(sol.r[sol.midplane_index]) # compute parallell heat flux at omp
            push!(Rwall,sol.r[1]) # R of points after midplane
            push!(Zwall,sol.z[1]) # z of points after midplane
            push!(Qwall, qmid/(sol.total_flux_expansion[1])*sin(sol.grazing_angles[1]))
            push!(Qpara, qmid/(sol.total_flux_expansion[1]))
            push!(indexes,sol.wall_index[1])
        end
    end


    if merge_wall
        if wall_z[1]>ZA # correction if first point is above midplane (wall not respecting COCOS11)
            indexes[indexes .== 1 .&& Zwall.>ZA] .= length(wall_r) # put it after index = end  
            indexes[indexes .== 1 .&& Zwall.>wall_z[2] .&& Zwall.<=ZA] .= 0   # put before index = 1
        end

        crossings = IMAS.intersection([(minimum(wall_r)+maximum(wall_r))/2, maximum(wall_r)*1.05], [ZA, ZA], wall_r, wall_z)[2] # (r,z) point of intersection btw outer midplane (OMP) with wall
        r_wall_midplane = [cr[1] for cr in crossings] # R coordinate of the wall at OMP
        r_wall_midplane = r_wall_midplane[1];
        psi_wall_midplane = PSI_interpolant(r_wall_midplane,ZA); 
        _, psi_first_lfs_far, null_within_wall = IMAS.find_psi_last_diverted(eqt, wall_r, wall_z, PSI_interpolant) # psi of grazing surface
        psi_wall = PSI_interpolant.(wall_r,wall_z)
        add_indexes =  collect((1:length(psi_wall)))
        add_indexes = add_indexes[(psi_wall .< psi_separatrix .|| psi_wall .> psi_wall_midplane) .||
                                (psi_wall .>= psi_separatrix .&& psi_wall .<= psi_first_lfs_far .&& 
                                sign(eqt.boundary.x_point[end].z).*wall_z.>abs(eqt.boundary.x_point[end].z)).&& null_within_wall]
        add_indexes = reverse!(add_indexes)

        for ind in add_indexes
            if ind !=1
            Rwall = append!(Rwall[indexes.<ind],[wall_r[ind]], Rwall[indexes.>=ind])
            Zwall = append!(Zwall[indexes.<ind],[wall_z[ind]], Zwall[indexes.>=ind])
            Qwall = append!(Qwall[indexes.<ind],[0.0], Qwall[indexes.>=ind])
            Qpara = append!(Qpara[indexes.<ind],[0.0], Qpara[indexes.>=ind])
            indexes = append!(indexes[indexes.<ind],[ind], indexes[indexes.>=ind])
            end
        end

        # include points shadowed by corners
        dist = (diff(Rwall)).^2 + (diff(Zwall)).^2
        avg = sum(dist)/length(dist) # mean distance between adjacent points in mesh
        std = sqrt(sum((dist.-avg).^2)/(length(dist)-1)) #standard deviation of distance adjacent points in mesh
        save_length = length(Rwall)
        counter = 0
        while sum(dist.> avg+3*std)>0 && counter < 10
            # we have outliers, which are shadowed areas.
            # in the outliears dist > average + 3 standard deviations
            indexes = 1:length(dist)
            indexes = indexes[dist.> avg+3*std]

            for ind in reverse(indexes)
                # if vertical lines
                if abs.(diff([Rwall[ind],Rwall[ind+1]]))[1]>0.01
                    dr = 0
                else
                    # dr = 1.0*maximum(sqrt.(diff(wall_r).^2 + diff(wall_z).^2))
                    dr = sqrt(2)*sum(sqrt.(diff(wall_r).^2 + diff(wall_z).^2))/length(dist)
                end
                #if horizontal lines
                if abs.(diff([Zwall[ind],Zwall[ind+1]]))[1]>0.01
                    dz = 0
                else
                    # dz = 1.0*maximum(sqrt.(diff(wall_r).^2 + diff(wall_z).^2))
                    dz = sqrt(2)*sum(sqrt.(diff(wall_r).^2 + diff(wall_z).^2))/length(dist)
                end

                #search points in rectangle between two points 
                add_r = wall_r[wall_r.>(minimum([Rwall[ind],Rwall[ind+1]])-dr ).&& wall_r .< (maximum([Rwall[ind],Rwall[ind+1]])+dr ).&&
                               wall_z.>(minimum([Zwall[ind],Zwall[ind+1]])-dz) .&& wall_z .< (maximum([Zwall[ind],Zwall[ind+1]])+dz )  ]
                add_z = wall_z[wall_r.>(minimum([Rwall[ind],Rwall[ind+1]])-dr) .&& wall_r .< (maximum([Rwall[ind],Rwall[ind+1]])+dr ).&&
                               wall_z.>(minimum([Zwall[ind],Zwall[ind+1]])-dz) .&& wall_z .< (maximum([Zwall[ind],Zwall[ind+1]])+dz)   ]
                Rwall = append!(Rwall[1:ind],add_r, Rwall[ind+1:end])
                Zwall = append!(Zwall[1:ind],add_z, Zwall[ind+1:end])
                Qwall = append!(Qwall[1:ind],add_r.*0.0, Qwall[ind+1:end])
                Qpara = append!(Qpara[1:ind],add_r.*0.0, Qpara[ind+1:end])        
            end
            if length(Rwall) == save_length
                # I am not finding anymore points to add
                    break
            end
            dist = (diff(Rwall)).^2 + (diff(Zwall)).^2
            avg = sum(dist)/length(dist) # mean distance between adjacent points in mesh
            std = sqrt(sum((dist.-avg).^2)/(length(dist)-1)) #standard deviation of distance adjacent points in mesh
            counter += 1
            save_length = length(Rwall)
        end
    end
    

    s = similar(Rwall)
    ds= sqrt.(diff(Rwall).^2 + diff(Zwall).^2)
    s[1] = 0
    for i in 1:length(ds)
        s[i+1]= s[i]+ds[i]
    end

    HF = WallHeatFlux(zeros(2), zeros(2), zeros(2), zeros(2), zeros(2))
    HF.r = Rwall
    HF.z = Zwall
    HF.q_wall = Qwall
    HF.q_parallel = Qpara
    HF.s = s

    return HF
end

function WallHFMapper(dd::IMAS.dd, 
    r::Vector{<:Real}, 
    q::Vector{<:Real}; 
    merge_wall::Bool = true, levels::Int = 20)

    eqt = dd.equilibrium.time_slice[]
    wall_r = IMAS.first_wall(dd.wall).r
    wall_z = IMAS.first_wall(dd.wall).z
    if isempty(wall_r) || isempty(wall_z)
        error("Impossible to map the heat flux onto the wall because dd.wall is empty")
    end
    
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt2d)  #interpolation of PSI in equilirium at locations (r,z)
    psi_levels, _, _  =  IMAS.find_levels_from_P(eqt,wall_r,wall_z,PSI_interpolant,r,q,levels)
    add_psi =  IMAS.find_levels_from_wall(eqt,wall_r,wall_z,PSI_interpolant)
   
    psi_levels = unique!(sort!(vcat(psi_levels, add_psi)))
    if psi_sign == -1
        psi_levels = reverse!(psi_levels) # if psi is decreasing, sort in descending order
    end

    _, psi_separatrix   = IMAS.find_psi_boundary(eqt; raise_error_on_not_open=true) # psi at LCFS
    
    psi_levels[1] = psi_separatrix

    SOL = IMAS.sol(eqt, wall_r, wall_z, levels = psi_levels, use_wall = true)

    return WallHFMapper(eqt, SOL, wall_r, wall_z, r, q; merge_wall)

end




@recipe function plot_heat_flux(HF::WallHeatFlux; which_plot= :twoD, plot_type = :path, cat = false, q=:wall) 
    if cat
        cat_jet = cgrad(:jet, categorical = true)
    else
        cat_jet = :jet
    end

    if which_plot == :twoD
        if q == :both
            layout := (1,2)  
        end

        if q == :wall .|| q == :both
            @series begin
                aspect_ratio --> :equal
                legend --> false
                colorbar --> true
                xlabel --> "R [m]"
                ylabel --> "Z [m]"
                colorbar_title := "log₁₀(q wall [W/m^2])"
                xlim --> [0.95*minimum(HF.r)-0.05*(maximum(HF.r)), 1.05*maximum(HF.r)-0.05*(minimum(HF.r))]
                ylim --> [1.05*minimum(HF.z)-0.05*(maximum(HF.z)), 1.05*maximum(HF.z)-0.05*(minimum(HF.z))]
                clim --> (floor(log10(maximum(HF.q_wall)))-8, floor(log10(maximum(HF.q_wall)))+1)
                if plot_type == :path
                    seriestype --> :path
                    line_z --> log10.(HF.q_wall.+1)
                    linewidth --> 3
                    color --> cat_jet
                end
                if plot_type == :scatter
                    seriestype --> :scatter
                    zcolor --> log10.(HF.q_wall.+1)
                    markersize --> 2
                    markerstrokewidth --> 0
                    color --> cat_jet
                end
                HF.r, HF.z
            end
         end

        if q == :parallel .|| q == :both
            @series begin
                aspect_ratio --> :equal
                legend --> false
                colorbar --> true
                xlabel --> "R [m]"
                ylabel --> "Z [m]"
                colorbar_title := "log₁₀(q parallel [W/m^2])"
                xlim --> [0.95*minimum(HF.r)-0.05*(maximum(HF.r)), 1.05*maximum(HF.r)-0.05*(minimum(HF.r))]
                ylim --> [1.05*minimum(HF.z)-0.05*(maximum(HF.z)), 1.05*maximum(HF.z)-0.05*(minimum(HF.z))]
                clim --> (floor(log10(maximum(HF.q_parallel)))-8, floor(log10(maximum(HF.q_parallel)))+1)
                if plot_type == :path
                    seriestype --> :path
                    line_z --> log10.(HF.q_parallel.+1)
                    linewidth --> 3
                    color --> cat_jet
                end
                if plot_type == :scatter
                    seriestype --> :scatter
                    zcolor --> log10.(HF.q_wall.+1)
                    markersize --> 2
                    markerstrokewidth --> 0
                    color --> cat_jet
                end
                HF.r, HF.z
            end
        end    
    end

    if which_plot == :oneD 
        if q == :wall .|| q == :both
            @series begin
                if q == :wall
                    legend --> :none
                else
                    label --> ""
                end
                xlabel --> "s [m]"
                ylabel --> "log₁₀(q wall [W/m^2])"
                color --> :dodgerblue
                yticks --> vcat([-Inf, 0],collect(1:maximum(log10.(HF.q_wall.+1))+1))
                HF.s,log10.(HF.q_wall .+ 1)
            end
            
            @series begin
                if q == :wall
                    legend --> :none
                else
                    label --> "q wall"
                end
                xlabel --> "s [m]"
                ylabel --> "log₁₀(q wall [W/m^2])"
                yticks --> vcat([-Inf, 0],collect(1:maximum(log10.(HF.q_wall.+1))+1))
                seriestype --> :scatter
                markersize --> 2
                color --> :dodgerblue
                HF.s,log10.(HF.q_wall .+ 1)
            end
        end

        if q == :parallel .|| q == :both
            @series begin
                if q == :parallel
                    legend --> :none
                else
                    label --> ""
                end
                xlabel --> "s [m]"
                ylabel --> "log₁₀(q parallel [W/m^2])"
                color --> :darkorange3
                yticks --> vcat([-Inf, 0],collect(1:maximum(log10.(HF.q_parallel.+1))+1))
                HF.s,log10.(HF.q_parallel .+ 1)
            end

            @series begin
                if q == :parallel
                    legend --> :none
                    ylabel --> "log₁₀(q parallel [W/m^2])"
                else
                    label --> "q parallel"
                    ylabel --> "log₁₀(q [W/m^2])"
                end
                xlabel --> "s [m]"
                yticks --> vcat([-Inf, 0],collect(1:maximum(log10.(HF.q_wall.+1))+1))
                seriestype --> :scatter
                markersize --> 2
                color --> :darkorange3
                HF.s,log10.(HF.q_parallel .+ 1)
            end
        end

    end
end 
