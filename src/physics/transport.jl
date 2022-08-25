function setup_transport_grid!(model::IMAS.core_transport__model, rho_gridpoints::Vector{<:Real})
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = rho_gridpoints
    m1d.electrons.particles.flux = zeros(length(rho_gridpoints))
    m1d.electrons.energy.flux = zeros(length(rho_gridpoints))
    m1d.total_ion_energy.flux = zeros(length(rho_gridpoints))
    m1d.momentum_tor.flux = zeros(length(rho_gridpoints))
end
