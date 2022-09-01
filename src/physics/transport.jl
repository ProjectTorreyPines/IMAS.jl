"""
    setup_transport_grid!(model::IMAS.core_transport__model, rho_gridpoints::Vector{<:Real})

Sets up the transport grid of core_transport.model and initializes the fluxes
"""
function setup_transport_grid!(m1d::IMAS.core_transport__model___profiles_1d, rho_gridpoints::Vector{<:Real})
    m1d.grid_flux.rho_tor_norm = rho_gridpoints
    m1d.electrons.particles.flux = zeros(length(rho_gridpoints))
    m1d.electrons.energy.flux = zeros(length(rho_gridpoints))
    m1d.total_ion_energy.flux = zeros(length(rho_gridpoints))
    m1d.momentum_tor.flux = zeros(length(rho_gridpoints))
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

"""
    blend_core_pedestal(
            profile::Vector{<:Real},
            rho::Vector{<:Real},
            ped_height::Real,
            ped_width::Real,
            rho_bound::Real,
            ; do_plot=false)

Blends the core and pedestal for given profile to match ped_height, ped_width using rho_bound as blending boundary
"""
function blend_core_pedestal(
    profile::Vector{<:Real},
    rho::Vector{<:Real},
    ped_height::Real,
    ped_width::Real,
    rho_bound::Real,
    ; do_plot=false)

    ngrid = length(rho)
    profile_ped = Hmode_profiles(profile[end], ped_height, profile[1], ngrid, 2.0, 1.4, ped_width)

    rho_top = 1 - 1.5 * ped_width

    if rho_top < rho_bound
        error("Unable to blend the core-pedestal becuase the rho_bound $rho_bound > rho_pedestal top $rho_top")
    end

    irho_top = argmin((rho .- rho_top) .^ 2)
    irho_bound = argmin((rho .- rho_bound) .^ 2)

    z_profile = -IMAS.calc_z(rho, profile)
    z_bound = z_profile[irho_bound]
    z_top = -IMAS.calc_z(rho, profile_ped)[irho_top]

    drho_nml = rho[irho_top] - rho[irho_bound]

    profile_new = deepcopy(profile_ped)


    for i in irho_bound:irho_top
        profile_new[i] = profile_ped[irho_top] * exp(0.5 * z_top / drho_nml * (drho_nml^2 - (rho[i] - rho[irho_bound])^2) + 0.5 * z_bound / drho_nml * (rho[i] - rho[irho_top])^2)
    end

    for i in irho_bound-1:-1:1
        profile_new[i] = profile_new[i+1] * exp(0.5 * (z_profile[i] + z_profile[i+1]) * (-rho[i] + rho[i+1]))
    end

    if do_plot
        import Plots
        Plots.plot(rho, profile, label="old profile")
        Plots.plot!(rho[irho_bound:end], profile_ped[irho_bound:end], lw=2.5, linestyle=:dash, label="new pedestal to blend")
        display(Plots.plot!(rho, profile_new, label="blended profile"))
    end

    return profile_new
end
