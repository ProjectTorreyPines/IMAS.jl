"""
    blend_core_edge_Hmode(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        nml_bound::Real,
        ped_bound::Real;
        expin::Real,
        expout::Real)

Blends the core and pedestal for given profile to match ped_height, ped_width using nml_bound as blending boundary
"""
function blend_core_edge_Hmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    nml_bound::Real,
    ped_bound::Real,
    expin::Real,
    expout::Real)

    @assert nml_bound <= ped_bound "Unable to blend the core-pedestal because the nml_bound $nml_bound > ped_bound top $ped_bound"
    iped = argmin(abs.(rho .- ped_bound))
    inml = argmin(abs.(rho .- nml_bound))

    z_profile = -calc_z(rho, profile)
    z_nml = z_profile[inml]

    # H-mode profile used for pedestal  (DON"T CHANGE THE -1 to profile[1] as this causes the optimizizer optimize with a lag)
    profile_ped = Hmode_profiles(profile[end], ped_height, -1, length(rho), expin, expout, ped_width)

    # linear z between nml and pedestal
    if nml_bound < ped_bound
        z_profile_ped = -calc_z(rho, profile_ped)
        z_ped = z_profile_ped[iped]
        z_profile[inml:iped] = (z_nml - z_ped) ./ (rho[inml] - rho[iped]) .* (rho[inml:iped] .- rho[inml]) .+ z_nml
    end

    # integrate from pedestal inward
    profile_new = deepcopy(profile_ped)
    profile_new[1:iped] = integ_z(rho[1:iped], z_profile[1:iped], profile_ped[iped])

    return profile_new
end


function cost_find_exps(x::AbstractVector, ped_height::Real, ped_width::Real, rho::AbstractVector, profile::AbstractVector, z_targets::AbstractVector, rho_targets::AbstractVector)
    x = abs.(x)
    profile_ped = Hmode_profiles(profile[end], ped_height, profile[1], length(rho), x[1], x[2], ped_width)

    z_ped = -calc_z(rho, profile_ped)
    z_ped_values = interp1d(rho, z_ped).(rho_targets)
    return sum(abs.(z_targets .- z_ped_values))
end


"""
    blend_core_edge_Hmode(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        tr_bound0::Real,
        tr_bound1::Real)

Blends the core profiles to the pedestal for H-mode profiles, this version makes sure the Z's of tr_bound0 and tr_bound1 match exactly with the Z's from the original profile
"""
function blend_core_edge_Hmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    tr_bound0::Real,
    tr_bound1::Real)

    z_profile = -calc_z(rho, profile)
    z_targets = interp1d(rho, z_profile).([tr_bound0, tr_bound1])

    # figure out expin and expout such that the Z's of Hmode_profiles match the z_targets from transport
    x_guess = [1.0, 1.0]
    res = Optim.optimize(x -> cost_find_exps(x, ped_height, ped_width, rho, profile, z_targets,
            [tr_bound0, tr_bound1]), x_guess, Optim.NelderMead(), Optim.Options(; g_tol=1E-3))

    expin = abs(res.minimizer[1])
    expout = abs(res.minimizer[2])

    return blend_core_edge_Hmode(
        profile,
        rho,
        ped_height,
        ped_width,
        tr_bound0,
        tr_bound1,
        expin,
        expout)
end

"""
    blend_core_edge_Hmode(cp1d::IMAS.core_profiles__profiles_1d, dd_ped::IMAS.summary__local__pedestal, edge_bound::Real)

Blends Te, Ti, ne, and nis in core_profiles with H-mode like pedestal defined in summary
"""
function blend_core_edge_Hmode(cp1d::IMAS.core_profiles__profiles_1d, dd_ped::IMAS.summary__local__pedestal, rho_nml::Real, rho_ped::Real)
    rho = cp1d.grid.rho_tor_norm
    w_ped = 1 - @ddtime(dd_ped.position.rho_tor_norm)

    ion_fractions = zeros(Float64, length(cp1d.ion), length(cp1d.electrons.density))
    for (ii, ion) in enumerate(cp1d.ion)
        ion_fractions[ii, :] = ion.density_thermal ./ cp1d.electrons.density
    end

    cp1d.electrons.temperature = blend_core_edge_Hmode(cp1d.electrons.temperature, rho, @ddtime(dd_ped.t_e.value), w_ped, rho_nml, rho_ped)
    cp1d.electrons.density_thermal = blend_core_edge_Hmode(cp1d.electrons.density, rho, @ddtime(dd_ped.n_e.value), w_ped, rho_nml, rho_ped)

    ti_avg_new = blend_core_edge_Hmode(cp1d.ion[1].temperature, rho, @ddtime(dd_ped.t_i_average.value), w_ped, rho_nml, rho_ped)
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

"""
    pedestal_finder(profile::AbstractVector{<:Real}, psi_norm::AbstractVector{<:Real})

Finds the pedetal height and width using the EPED1 definition
returns ped_height, ped_width
"""
function pedestal_finder(profile::Vector{T}, psi_norm::Vector{T}) where {T<:Real}
    function cost_function(profile, weight_func, params)
        if any(x -> (x < 0.0), params)
            return 1e10
        end
        profile_fit = Hmode_profiles(profile[end], params[1], profile[1], length(profile), 2.0, 2.0, params[2])
        return sqrt(sum(((profile .- profile_fit) .^ 2 ./ profile[1]^2) .* weight_func))
    end
    ngrid = length(profile)
    half_grid = Int(floor(ngrid / 2))

    inversion_point = argmin(gradient(psi_norm[half_grid:end], profile[half_grid:end])) + half_grid
    inversion_point_margin = inversion_point - Int(floor(0.1 * ngrid))

    weight_func = zeros(ngrid)
    weight_func[inversion_point_margin:end] .+= 1.0

    width0 = 1.0 - psi_norm[inversion_point]
    guess = [interp1d(psi_norm, profile)(1.0 - 2.0 * width0), width0]
    res = Optim.optimize(params -> cost_function(profile, weight_func, params), guess, Optim.NelderMead(), Optim.Options(; g_tol=1E-5))

    return res.minimizer
end