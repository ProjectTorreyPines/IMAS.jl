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
    profile_new[1:inml] = profile[1:inml] .- profile[inml+1] .+ profile_new[inml+1]

    return profile_new
end

function cost_find_EPED_exps(
    x::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    rho::AbstractVector{<:Real},
    profile::AbstractVector{<:Real},
    z_targets::AbstractVector{<:Real},
    rho_targets::AbstractVector{<:Real}
)
    x = abs.(x)
    profile_ped = Hmode_profiles(profile[end], ped_height, length(rho), x[1], x[2], ped_width)
    z_ped = -calc_z(rho, profile_ped, :backward)
    z_ped_values = interp1d(rho, z_ped).(rho_targets)

    return norm(z_targets .- z_ped_values)
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
    tr_bound1::Real
)
    z_profile = -calc_z(rho, profile, :backward)
    rho_targets = [tr_bound0, tr_bound1]
    z_targets = interp1d(rho, z_profile).(rho_targets)
    
    # figure out expin and expout such that the Z's of Hmode_profiles match the z_targets from transport
    x_guess = [1.0, 1.0]
    res = Optim.optimize(x -> cost_find_EPED_exps(x, ped_height, ped_width, rho, profile, z_targets, rho_targets), x_guess, Optim.NelderMead())
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
    rho_bound::Real
)
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
    tr_bound1::Real
)
    return blend_core_edge_Lmode(profile, rho, ped_height, tr_bound1)
end

function cost_WPED_α!(rho::AbstractVector{<:Real}, profile::AbstractVector{<:Real}, α::Real, value::Real, rho_ped::Real)
    rho_ped_idx = argmin(abs.(rho .- rho_ped))

    profile_ped = IMAS.edge_profile(rho, rho_ped, value, profile[end], α)
    z_profile_ped = IMAS.calc_z(rho, profile_ped, :backward)

    profile .+= (-profile[rho_ped_idx] + value)
    z_profile = IMAS.calc_z(rho, profile, :backward)

    profile[rho_ped_idx+1:end] .= IMAS.interp1d(rho, profile_ped).(rho[rho_ped_idx+1:end])

    cost = abs.((z_profile[rho_ped_idx] - z_profile_ped[rho_ped_idx]) / z_profile[rho_ped_idx])
    return cost
end

"""
    pedestal_finder(profile::AbstractVector{<:Real}, psi_norm::AbstractVector{<:Real})

Finds the pedetal height and width using the EPED1 definition

returns ped_height, ped_width
"""
function pedestal_finder(profile::Vector{T}, psi_norm::Vector{T}) where {T<:Real}
    function cost_function(profile, params)
        params = abs.(params)
        height = params[1]
        width = min(max(0.0, params[2]), 0.5)
        expin = params[3]
        expout = params[4]
        psi_norm_fit = range(0, 1, length(profile))

        tmp = Hmode_profiles(profile[end], height, profile[1], length(profile), expin, expout, width)
        profile_fit = interp1d(psi_norm_fit, tmp).(psi_norm)

        p1 = abs.(profile .- profile_fit) ./ maximum(profile)
        p2 = abs.(gradient(psi_norm, profile) .- gradient(psi_norm, profile_fit)) ./ maximum(abs.(gradient(psi_norm, profile)))

        return sum(((p1 .+ p2) .* weight) .^ 2.0)
    end

    weight = gradient(psi_norm)

    width0 = 0.05
    height0 = interp1d(psi_norm, profile)(1.0 - width0)
    expin0 = 1.0
    expout0 = 1.0
    guess = [height0, width0, expin0, expout0]
    res = Optim.optimize(params -> cost_function(profile, params), guess, Optim.NelderMead(), Optim.Options(; g_tol=1E-5))

    params = abs.(res.minimizer)
    height = params[1]
    width = min(max(0.0, params[2]), 0.5)
    expin = params[3]
    expout = params[4]

    # psi_norm_fit = range(0, 1, length(profile_fit))
    # profile_fit = Hmode_profiles(profile[end], height, profile[1], length(profile), expin, expout, width)
    # p = plot(psi_norm, profile)
    # plot!(p, psi_norm_fit, profile_fit)
    # hline!([height])
    # vline!([1.0 .- width])
    # scatter!(p, [1.0 - width], [interp1d(psi_norm_fit, profile_fit)(1.0 - width)])
    # display(p)

    return height, width
end
