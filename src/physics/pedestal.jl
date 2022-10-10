"""
    blend_core_edge_Hmode(
            profile::AbstractVector{<:Real},
            rho::AbstractVector{<:Real},
            ped_height::Real,
            ped_width::Real,
            nml_bound::Real)

Blends the core and pedestal for given profile to match ped_height, ped_width using nml_bound as blending boundary
"""
function blend_core_edge_Hmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    nml_bound::Real,
    ped_bound::Real=1.0 - 1.5 * ped_width;
    expin::Real,
    expout::Real)

    @assert nml_bound < ped_bound "Unable to blend the core-pedestal because the nml_bound $nml_bound > ped_bound top $ped_bound"
    iped = argmin(abs.(rho .- ped_bound))
    inml = argmin(abs.(rho .- nml_bound))

    z_profile = -calc_z(rho, profile)
    z_nml = z_profile[inml]

    # H-mode profile used for pedestal
    profile_ped = Hmode_profiles(profile[end], ped_height, -1.0, length(rho), expin, expout, ped_width)

    # linear z between nml and pedestal
    z_profile_ped = -calc_z(rho, profile_ped)
    z_ped = z_profile_ped[iped]
    z_profile[inml:iped] = (z_nml - z_ped) ./ (rho[inml] - rho[iped]) .* (rho[inml:iped] .- rho[inml]) .+ z_nml

    # integrate from pedestal inward
    profile_new = deepcopy(profile_ped)
    profile_new[1:iped] = integ_z(rho[1:iped], z_profile[1:iped], profile_ped[iped])

    return profile_new
end

"""
    blend_core_edge_Hmode(cp1d::IMAS.core_profiles__profiles_1d, dd_ped::IMAS.summary__local__pedestal, edge_bound::Real)

Blends Te, Ti, ne, and nis in core_profiles with H-mode like pedestal defined in summary
"""
function blend_core_edge_Hmode(cp1d::IMAS.core_profiles__profiles_1d, dd_ped::IMAS.summary__local__pedestal, edge_bound::Real)
    rho = cp1d.grid.rho_tor_norm
    w_ped = 1 - @ddtime(dd_ped.position.rho_tor_norm)

    ion_fractions = zeros(Float64, length(cp1d.ion), length(cp1d.electrons.density))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ cp1d.electrons.density
    end

    cp1d.electrons.temperature = blend_core_edge_Hmode(cp1d.electrons.temperature, rho, @ddtime(dd_ped.t_e.value), w_ped, edge_bound; expin=1.2, expout=1.4)
    cp1d.electrons.density_thermal = blend_core_edge_Hmode(cp1d.electrons.density, rho, @ddtime(dd_ped.n_e.value), w_ped, edge_bound; expin=1.1, expout=1.1)

    ti_avg_new = blend_core_edge_Hmode(cp1d.ion[1].temperature, rho, @ddtime(dd_ped.t_i_average.value), w_ped, edge_bound; expin=1.2, expout=1.4)
    for (ii, ion) in enumerate(cp1d.ion)
        ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density
        ion.temperature = ti_avg_new
    end
end

"""
    blend_core_edge_Lmode(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        edge_bound::Real)

Generates L-mode like profile by enforcing constant gradient between edge_bound and edge
"""
function blend_core_edge_Lmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    edge_bound::Real)

    inml = argmin(abs.(rho .- edge_bound))

    grad_profile = -gradient(rho, profile; method=:backward)
    grad_edge = grad_profile[inml]
    edge = profile[end] .+ (1.0 .- rho) .* grad_edge

    profile = profile .- profile[inml] .+ edge[inml]

    profile[inml:end] = edge[inml:end]
    return profile
end

"""
    blend_core_edge_Lmode(cp1d::IMAS.core_profiles__profiles_1d)

Blends Te, Ti, ne, and nis in core_profiles with L-mode like pedestal
"""
function blend_core_edge_Lmode(cp1d::IMAS.core_profiles__profiles_1d, edge_bound::Real)
    rho = cp1d.grid.rho_tor_norm

    ion_fractions = zeros(Float64, length(cp1d.ion), length(cp1d.electrons.density))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ cp1d.electrons.density
    end

    cp1d.electrons.temperature = blend_core_edge_Lmode(cp1d.electrons.temperature, rho, edge_bound)
    cp1d.electrons.density_thermal = blend_core_edge_Lmode(cp1d.electrons.density, rho, edge_bound)

    ti_avg_new = blend_core_edge_Lmode(cp1d.ion[1].temperature, rho, edge_bound)
    for (ii, ion) in enumerate(cp1d.ion)
        ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density
        ion.temperature = ti_avg_new
    end
end