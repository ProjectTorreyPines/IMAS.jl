"""
    setup_transport_grid!(model::IMAS.core_transport__model, rho_gridpoints::Vector{<:Real})

Sets up the transport grid of core_transport.model and initializes the fluxes
"""
function setup_transport_grid!(m1d::IMAS.core_transport__model___profiles_1d, rho_gridpoints::Vector{<:Real}; setup_grid_for::Vector{Symbol}=[:electrons_energy, :electrons_particles, :total_ion_energy, :momentum_tor])
    m1d.grid_flux.rho_tor_norm = rho_gridpoints
    if :electrons_particles ∈ setup_grid_for
        m1d.electrons.particles.flux = zeros(length(rho_gridpoints))
    end
    if :electrons_energy ∈ setup_grid_for
        m1d.electrons.energy.flux = zeros(length(rho_gridpoints))
    end
    if :total_ion_energy ∈ setup_grid_for
        m1d.total_ion_energy.flux = zeros(length(rho_gridpoints))
    end
    if :momentum_tor ∈ setup_grid_for
        m1d.momentum_tor.flux = zeros(length(rho_gridpoints))
    end
end

"""
    profile_from_z_transport(
        profile_old::AbstractVector{<:Real},
        rho::Vector{<:Real},
        transport_grid::Vector{<:Real},
        z_transport_grid::Vector{<:Real})

Updates profile_old with the scale lengths given by z_transport_grid
"""
function profile_from_z_transport(
    profile_old::AbstractVector{<:Real},
    rho::Vector{<:Real},
    transport_grid::Vector{<:Real},
    z_transport_grid::Vector{<:Real})

    transport_idices = [argmin((rho_x .- rho) .^ 2) for rho_x in transport_grid]
    transport_idices = vcat(1, transport_idices)
    z_transport_grid = vcat(0.0, z_transport_grid)

    rho_transport_grid = rho[transport_idices]

    z_old = calc_z(rho, profile_old)
    z_new = similar(z_old)
    z_new[transport_idices[end]:end] = z_old[transport_idices[end]:end]

    z_new[1:transport_idices[end]] = IMAS.interp1d(rho_transport_grid, z_transport_grid).(rho[1:transport_idices[end]])
    profile_new = similar(profile_old)
    profile_new[transport_idices[end]:end] = profile_old[transport_idices[end]:end]
    for i in transport_idices[end]-1:-1:1
        profile_new[i] = profile_new[i+1] * exp(0.5 * (z_new[i] + z_new[i+1]) * (rho[i] - rho[i+1]))
    end
    return profile_new
end

function total_fluxes(dd::IMAS.dd)
    return total_fluxes(dd.core_transport)
end


"""
    total_fluxes!(ct::IMAS.core_transport,rho_total_fluxes::AbstractVector{<:Real} = collect(0.0:0.05:1.0))

Sums up all the fluxes and returns it as a core_transport.model IDS
"""
function total_fluxes(ct::IMAS.core_transport,rho_total_fluxes::AbstractVector{<:Real} = collect(0.0:0.05:1.0))
    total_fluxes = IMAS.core_transport__model___profiles_1d()
    total_fluxes.grid_flux.rho_tor_norm = rho_total_fluxes
    skip_flux_list = [0, 1, 25]
    for model in ct.model
        if model.identifier.index ∈ skip_flux_list
            if model.identifier.index ∈ [0, 25]
                @warn "skipped model.identifier.index = $(model.identifier.index), do not use this index"
            end
            continue
        end
        append!(skip_flux_list, model.identifier.index) # Make sure we don't double count a specific flux type

        m1d = model.profiles_1d[]
        for sub in [:electrons, :momentum_tor, :total_ion_energy]
            ids1 = m1d
            ids2 = total_fluxes

            if sub == :electrons
                ids1 = getproperty(ids1, sub)
                ids2 = getproperty(ids2, sub)
                iterator = deleteat!(keys(ids1), findall(x-> x ∈ [:electrons, :grid_flux, :time],keys(ids1)))
            else
                iterator = [sub]
            end
            for field in iterator
                y = getproperty(getproperty(ids1, field, missing),:flux, missing)
                if typeof(y) <: AbstractVector{<:Real}
                    old_value = getproperty(getproperty(ids2, field),:flux,zeros(length(rho_total_fluxes)))
                    x = m1d.grid_flux.rho_tor_norm
                    x_1 = argmin(abs.(rho_total_fluxes .- x[1]))
                    x_2 = argmin(abs.(rho_total_fluxes .- x[end]))
                    if !ismissing(getproperty(ids2,field,missing))
                        if length(x) == 1
                            old_value[x_1] += y[1]
                            setproperty!(getproperty(ids2,field), :flux, old_value)
                        else
                            setproperty!(getproperty(ids2,field), :flux, 
                                vcat(old_value[1:x_1-1],
                                    old_value[x_1:x_2] .+  IMAS.interp1d(x, y,:cubic).(rho_total_fluxes)[x_1:x_2],
                                    old_value[x_2+1:end]))
                        end
                    end
                end
            end
        end
    end
    return total_fluxes
end


