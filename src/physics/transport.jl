"""
    setup_transport_grid!(model::IMAS.core_transport__model, rho_gridpoints::Vector{<:Real})

Sets up the transport grid of core_transport.model and initializes the fluxes
"""
function setup_transport_grid!(model::IMAS.core_transport__model, rho_gridpoints::Vector{<:Real})
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = rho_gridpoints
    m1d.electrons.particles.flux = zeros(length(rho_gridpoints))
    m1d.electrons.energy.flux = zeros(length(rho_gridpoints))
    m1d.total_ion_energy.flux = zeros(length(rho_gridpoints))
    m1d.momentum_tor.flux = zeros(length(rho_gridpoints))
end

"""
update_profile_with_z_transport(profile_old::Vector{<:Real},
                                rho::Vector{<:Real},
                                transport_grid::Vector{<:Real},
                                z_transport_grid::Vector{<:Real})

Updates profile_old with the scale lengths given by z_transport_grid
"""
function update_profile_with_z_transport(profile_old::Vector{<:Real}, rho::Vector{<:Real}, transport_grid::Vector{<:Real}, z_transport_grid::Vector{<:Real})
    profile_transport_grid = profile_old[transport_grid]
    rho_transport_grid = rho[transport_grid]

    z_old = calcz(rho, profile_old)
    z_new = similar(z_old)
    z_new[transport_grid[end]:end] = z_old[transport_grid[end]:end]

    z_new[1:transport_grid[end]] = IMAS.interp1d(rho_transport_grid, z_transport_grid).(rho[1:transport_grid[end]])
    profile_transport_grid = profile_old[transport_grid]
    profile_new = similar(profile_old)
    profile_new[transport_grid[end]:end] = profile_old[transport_grid[end]:end]

    for i in transport_grid[end]-1:-1:1
        profile_new[i] = profile_new[i+1] * exp(0.5 * (z_new[i] + z_new[i+1]) * (rho[i] - rho[i+1]))
    end

    return profile_new
end

