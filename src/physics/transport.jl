document[Symbol("Physics transport")] = Symbol[]

"""
    profile_from_z_transport(
        profile_old::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        transport_grid::AbstractVector{<:Real},
        z_transport_grid::AbstractVector{<:Real},
        rho_ped::Real=0.0)

Updates profile_old with the scale lengths given by z_transport_grid

If `rho_ped > transport_grid[end]` then scale-length is linearly interpolated between `transport_grid[end]` and `rho_ped`
if `rho_ped < transport_grid[end]` then scale-length then boundary condition is at `transport_grid[end]`
"""
function profile_from_z_transport(
    profile_old::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    transport_grid::AbstractVector{<:Real},
    z_transport_grid::AbstractVector{<:Real},
    rho_ped::Real=0.0)

    transport_indices = [argmin(abs.(rho .- rho_x)) for rho_x in transport_grid]
    index_ped = argmin(abs.(rho .- rho_ped))
    index_last = transport_indices[end]
    if index_ped > index_last
        z_old = calc_z(rho, profile_old, :backward)
        transport_indices = vcat(1, transport_indices, index_ped)
        z_transport_grid = vcat(0.0, z_transport_grid, z_old[index_ped])
    else
        transport_indices = vcat(1, transport_indices)
        z_transport_grid = vcat(0.0, z_transport_grid)
    end

    z = interp1d(transport_indices, z_transport_grid).(1:index_last)

    profile_new = similar(profile_old)
    profile_new[index_last:end] = @views profile_old[index_last:end]
    profile_new[1:index_last] = @views integ_z(rho[1:index_last], -z, profile_new[index_last])

    return profile_new
end

@compat public profile_from_z_transport
push!(document[Symbol("Physics transport")], :profile_from_z_transport)

"""
    total_fluxes(
        core_transport::IMAS.core_transport{T},
        cp1d::IMAS.core_profiles__profiles_1d,
        rho_total_fluxes::AbstractVector{<:Real};
        time0::Float64) where {T<:Real}

Sums up all the fluxes and returns it as a core_transport.model IDS
"""
function total_fluxes(
    core_transport::IMAS.core_transport{T},
    cp1d::IMAS.core_profiles__profiles_1d,
    rho_total_fluxes::AbstractVector{<:Real};
    time0::Float64) where {T<:Real}

    total_flux = resize!(core_transport.model, :combined; wipe=false)
    total_flux1d = resize!(total_flux.profiles_1d, time0; wipe=false)
    total_flux1d.grid_flux.rho_tor_norm = rho_total_fluxes

    skip_flux_list = [:unknown, :unspecified, :combined]
    index_to_name = index_2_name(core_transport.model)

    # initialize ions
    for ion in cp1d.ion
        resize!(total_flux1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
    end
    for model in core_transport.model
        if !isempty(model.profiles_1d) && time0 >= model.profiles_1d[1].time
            model1d = model.profiles_1d[time0]
            for ion in model1d.ion
                resize!(total_flux1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
            end
        end
    end

    # defines paths to fill
    paths = []
    push!(paths, (:electrons, :energy))
    push!(paths, (:electrons, :particles))
    for k in eachindex(total_flux1d.ion)
        push!(paths, (:ion, k, :energy))
        push!(paths, (:ion, k, :particles))
    end
    push!(paths, (:momentum_tor,))
    push!(paths, (:total_ion_energy,))

    # zero out total_flux
    for path in paths
        ids1 = goto(total_flux1d, path)
        if hasdata(ids1, :flux)
            getproperty(ids1, :flux) .*= 0.0
        else
            setproperty!(ids1, :flux, zeros(T, size(rho_total_fluxes)))
        end
    end

    if !isempty(rho_total_fluxes)
        for model in core_transport.model
            # Make sure we don't double count a specific flux type
            if index_to_name[model.identifier.index] ∈ skip_flux_list
                @debug "skipped model.identifier.index = $(model.identifier.index) [$(index_to_name[model.identifier.index])]"
                continue
            end
            push!(skip_flux_list, index_to_name[model.identifier.index])

            m1d = model.profiles_1d[time0]
            x = m1d.grid_flux.rho_tor_norm
            x_1 = argmin(abs.(rho_total_fluxes .- x[1]))
            x_2 = argmin(abs.(rho_total_fluxes .- x[end]))

            for path in paths
                ids1 = try
                    goto(m1d, path)
                catch e
                    if isa(e, InterruptException)
                        retrhow(e)
                    end
                    continue
                end
                if !hasdata(ids1, :flux)
                    continue
                end
                ids2 = goto(total_flux1d, path)
                y = getproperty(ids1, :flux)
                if rho_total_fluxes == x
                    getproperty(ids2, :flux)[x_1:x_2] .+= y
                else
                    getproperty(ids2, :flux)[x_1:x_2] .+= interp1d(x, y, :linear).(rho_total_fluxes[x_1:x_2])
                end
            end
        end
    end

    return total_flux1d
end

"""
    total_fluxes(dd::IMAS.dd, rho_total_fluxes::AbstractVector{<:Real}=dd.core_profiles.profiles_1d[].grid.rho_tor_norm; time0::Float64=dd.global_time)
"""
function total_fluxes(dd::IMAS.dd, rho_total_fluxes::AbstractVector{<:Real}=dd.core_profiles.profiles_1d[].grid.rho_tor_norm; time0::Float64=dd.global_time)
    return total_fluxes(dd.core_transport, dd.core_profiles.profiles_1d[], rho_total_fluxes; time0)
end

@compat public total_fluxes
push!(document[Symbol("Physics transport")], :total_fluxes)