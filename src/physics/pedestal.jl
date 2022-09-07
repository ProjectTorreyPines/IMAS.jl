"""
    blend_core_pedestal_Hmode(
            profile::Vector{<:Real},
            rho::Vector{<:Real},
            ped_height::Real,
            ped_width::Real,
            rho_bound::Real,
            ; do_plot=false)

Blends the core and pedestal for given profile to match ped_height, ped_width using rho_bound as blending boundary
"""
function blend_core_pedestal_Hmode(
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
    """ got to make this to recipes
    if do_plot
        plot(rho, profile, label="old profile")
        plot!(rho[irho_bound:end], profile_ped[irho_bound:end], lw=2.5, linestyle=:dash, label="new pedestal to blend")
        display(plot!(rho, profile_new, label="blended profile"))
    end
    """

    return profile_new
end

"""
    blend_core_pedestal_Hmode(dd::IMAS.dd)

Blends Te,Ti, ne, nis core_pedestal for Hmodeprofiles
"""
function blend_core_pedestal_Hmode(dd::IMAS.dd)
    dd_ped = dd.summary.local.pedestal
    cp1d = dd.core_profiles.profiles_1d[]
    rho = cp1d.grid.rho_tor_norm
    w_ped = 1 - @ddtime(dd_ped.position.rho_tor_norm)

    ion_fractions = zeros(Float64, length(cp1d.ion), length(cp1d.electrons.density))
    for (ii,ion) in enumerate(cp1d.ion)
        ion_fractions[ii,:] = ion.density_thermal ./ cp1d.electrons.density
    end
    cp1d.electrons.temperature = blend_core_pedestal_Hmode(cp1d.electrons.temperature, rho, @ddtime(dd_ped.t_e.value),w_ped, 0.8, do_plot=true)
    cp1d.electrons.density_thermal = blend_core_pedestal_Hmode(cp1d.electrons.density, rho, @ddtime(dd_ped.n_e.value),w_ped, 0.8, do_plot=true)
    ti_avg_new  = blend_core_pedestal_Hmode(cp1d.ion[1].temperature, rho, @ddtime(dd_ped.t_i_average.value),w_ped, 0.8, do_plot=true)

    for (ii,ion) in enumerate(cp1d.ion)
        ion.density_thermal = ion_fractions[ii,:] .*  cp1d.electrons.density
        ion.temperature = ti_avg_new
    end
end