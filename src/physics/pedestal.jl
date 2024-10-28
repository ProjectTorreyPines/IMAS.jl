"""
    blend_core_edge_EPED(
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
function blend_core_edge_EPED(
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

    z_profile = -calc_z(rho, profile, :backward)
    z_nml = z_profile[inml]

    # H-mode profile used for pedestal
    # NOTE: Note that we do not provide a core value as this causes the pedestal solution to depend on the core solution
    #       which breaks finding self-consistent core-pedestal solution through an optimizer
    profile_ped = Hmode_profiles(profile[end], ped_height, length(rho), expin, expout, ped_width)

    # linear z between nml and pedestal
    if nml_bound < ped_bound
        z_profile_ped = -calc_z(rho, profile_ped, :backward)
        z_ped = z_profile_ped[iped]
        z_profile[inml:iped] = (z_nml - z_ped) ./ (rho[inml] - rho[iped]) .* (rho[inml:iped] .- rho[inml]) .+ z_nml
    end

    # integrate from pedestal inward
    profile_new = deepcopy(profile_ped)
    profile_new[inml:iped] = integ_z(rho[inml:iped], z_profile[inml:iped], profile_ped[iped])

    # we avoid integ_z in the core region to avoid drift of profiles
    # when calling blend_core_edge_EPED multiple times
    profile_new[1:inml-1] = profile[1:inml-1] .- profile[inml] .+ profile_new[inml]

    return profile_new
end

"""
    blend_core_edge_Hmode(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        tr_bound0::Real,
        tr_bound1::Real)

Blends the core profiles to the pedestal for H-mode profiles, making sure the Z's at tr_bound0 and tr_bound1 match the Z's from the original profile
"""
function blend_core_edge_Hmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    tr_bound0::Real,
    tr_bound1::Real)

    function cost_find_EPED_exps(
        x::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        rho::AbstractVector{<:Real},
        profile::AbstractVector{<:Real},
        p_targets::AbstractVector{<:Real},
        z_targets::AbstractVector{<:Real},
        rho_targets::AbstractVector{<:Real})

        x = abs.(x)
        profile_ped = Hmode_profiles(profile[end], ped_height, length(rho), x[1], x[2], ped_width)
        z_ped = -calc_z(rho, profile_ped, :backward)
        z_ped_values = interp1d(rho, z_ped).(rho_targets)
        p_values = interp1d(rho, profile_ped).(rho_targets)

        # Z's can be matched both by varying the value as well as the gradients
        # Here we want to keep the values as fixed as possible, while varying the gradients
        return norm(z_targets .- z_ped_values) / sum(abs.(z_targets)) .+ norm(p_targets .- p_values) / sum(abs.(p_targets))
    end

    z_profile = -calc_z(rho, profile, :backward)
    rho_targets = [tr_bound0, tr_bound1]
    z_targets = interp1d(rho, z_profile).(rho_targets)
    p_targets = interp1d(rho, profile).(rho_targets)

    # figure out expin and expout such that the Z's of Hmode_profiles match the z_targets from transport
    x_guess = [1.0, 1.0]
    res = Optim.optimize(x -> cost_find_EPED_exps(x, ped_height, ped_width, rho, profile, p_targets, z_targets, rho_targets), x_guess, Optim.NelderMead())
    expin = abs(res.minimizer[1])
    expout = abs(res.minimizer[2])

    return blend_core_edge_EPED(
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
    blend_core_edge(mode::Symbol, cp1d::IMAS.core_profiles__profiles_1d, summary_ped::IMAS.summary__local__pedestal, rho_nml::Real, rho_ped::Real; what::Symbol=:all)

Blends Te, Ti, ne, and nis in core_profiles with :H_mode or :L_mode like pedestal defined in summary IDS
"""
function blend_core_edge(mode::Symbol, cp1d::IMAS.core_profiles__profiles_1d, summary_ped::IMAS.summary__local__pedestal, rho_nml::Real, rho_ped::Real; what::Symbol=:all)
    if mode == :L_mode
        blend_function = blend_core_edge_Lmode
    elseif mode == :H_mode
        blend_function = blend_core_edge_Hmode
    else
        @assert (mode ∈ (:L_mode, :H_mode)) "Mode can be either :L_mode or :H_mode"
    end
    rho = cp1d.grid.rho_tor_norm
    w_ped = 1.0 - @ddtime(summary_ped.position.rho_tor_norm)

    # NOTE! this does not take into account summary.local.pedestal.zeff.value
    if what ∈ (:all, :densities)
        ion_fractions = zeros(Float64, length(cp1d.ion), length(cp1d.electrons.density))
        for (ii, ion) in enumerate(cp1d.ion)
            ion_fractions[ii, :] = ion.density_thermal ./ cp1d.electrons.density
        end
        cp1d.electrons.density_thermal = blend_function(cp1d.electrons.density, rho, @ddtime(summary_ped.n_e.value), w_ped, rho_nml, rho_ped)
        for (ii, ion) in enumerate(cp1d.ion)
            ion.density_thermal = ion_fractions[ii, :] .* cp1d.electrons.density
        end
    end

    if what ∈ (:all, :temperatures)
        cp1d.electrons.temperature = blend_function(cp1d.electrons.temperature, rho, @ddtime(summary_ped.t_e.value), w_ped, rho_nml, rho_ped)
        ti_avg_new = blend_function(cp1d.t_i_average, rho, @ddtime(summary_ped.t_i_average.value), w_ped, rho_nml, rho_ped)
        for ion in cp1d.ion
            if !ismissing(ion, :temperature)
                ion.temperature = ti_avg_new
            end
        end
    end
end

function blend_core_edge_Lmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    value::Real,
    rho_bound::Real)

    res = Optim.optimize(α -> cost_WPED_α!(rho, profile, α, value, rho_bound), -500, 500, Optim.GoldenSection(); rel_tol=1E-3)
    cost_WPED_α!(rho, profile, res.minimizer, value, rho_bound)

    return profile
end

"""
    blend_core_edge_Lmode(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        tr_bound0::Real,
        tr_bound1::Real)

Blends the core profiles to the pedestal for L-mode profiles, making sure the Z's at tr_bound1 matches the Z's from the original profile

NOTE: ped_width, tr_bound0, tr_bound1 are not utilized
"""
function blend_core_edge_Lmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    tr_bound0::Real,
    tr_bound1::Real)
    return blend_core_edge_Lmode(profile, rho, ped_height, tr_bound1)
end

function cost_WPED_α!(rho::AbstractVector{<:Real}, profile::AbstractVector{<:Real}, α::Real, value_ped::Real, rho_ped::Real)
    rho_ped_idx = argmin(abs.(rho .- rho_ped))

    profile_ped = edge_profile(rho, rho_ped, value_ped, profile[end], α)
    z_profile_ped = calc_z(rho, profile_ped, :backward)

    profile .+= (-profile[rho_ped_idx] + value_ped)
    z_profile = calc_z(rho, profile, :backward)

    profile[rho_ped_idx+1:end] .= interp1d(rho, profile_ped).(rho[rho_ped_idx+1:end])

    cost = abs.((z_profile[rho_ped_idx] - z_profile_ped[rho_ped_idx]) / z_profile[rho_ped_idx])
    return cost
end

"""
    pedestal_finder(profile::AbstractVector{<:Real}, psi_norm::AbstractVector{<:Real})

Finds the pedetal height and width using the EPED1 definition

returns ped_height, ped_width
"""
function pedestal_finder(profile::Vector{T}, psi_norm::Vector{T}) where {T<:Real}
    psi_norm_fit = range(0, 1, length(profile))
    mask = psi_norm .> 0.5
    expin = 1.0
    expout = 1.0

    function cost_function(profile, params)
        width = abs(params[1]) + 0.01
        core = abs(params[2])
        height = interp1d(psi_norm, profile)(1.0 - width)

        tmp = Hmode_profiles(profile[end], height, core, length(profile), expin, expout, width)
        profile_fit = interp1d(psi_norm_fit, tmp).(psi_norm)

        cost = sqrt(trapz(psi_norm, mask .* (profile .- profile_fit) .^ 2)) / trapz(psi_norm, mask .* abs.(profile))
        return cost
    end

    res = Optim.optimize(params -> cost_function(profile, params), [0.04, profile[1]], Optim.NelderMead())

    width = abs(res.minimizer[1]) + 0.01
    core = abs(res.minimizer[2])
    height = interp1d(psi_norm, profile)(1.0 - width)

    # profile_fit = Hmode_profiles(profile[end], height, core, length(profile), expin, expout, width)
    # p = plot(psi_norm, profile; label="profile", marker=:circle, markersize=1)
    # plot!(p, psi_norm_fit, profile_fit; label="fit")
    # hline!(p, [height]; ls=:dash, primary=false)
    # vline!(p, [1.0 .- width]; ls=:dash, primary=false)
    # vline!(p, [1.0 - 2.0 * width]; ls=:dash, primary=false)
    # scatter!(p, [1.0 - width], [interp1d(psi_norm_fit, profile_fit)(1.0 - width)]; primary=false)
    # display(p)

    return height, width
end
