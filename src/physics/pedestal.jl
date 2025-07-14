document[Symbol("Physics pedestal")] = Symbol[]

function diffnorm(A::AbstractVector{<:Real}, B::AbstractVector{<:Real})
    @assert length(A) == length(B)
    return @inbounds norm(A[k] - B[k] for k in eachindex(A))
end

"""
    blend_core_edge_Hmode(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        tr_bound0::Real,
        tr_bound1::Real;
        method::Symbol=:shift)

Blends the core profiles to the pedestal for H-mode profiles, making sure the Z's at tr_bound0 and tr_bound1 match the Z's from the original profile
"""
function blend_core_edge_Hmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    tr_bound0::Real,
    tr_bound1::Real;
    method::Symbol=:shift)

    function cost_find_EPED_exps(
        x::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        rho::AbstractVector{<:Real},
        profile::AbstractVector{<:Real},
        p_targets::AbstractVector{<:Real},
        z_targets::AbstractVector{<:Real},
        rho_targets::AbstractVector{<:Real})

        expin = abs(x[1])
        expout = abs(x[2])
        profile_ped = Hmode_profiles(profile[end], ped_height, length(rho), expin, expout, ped_width)
        z_ped = -calc_z(rho, profile_ped, :backward)
        z_ped_values = interp1d(rho, z_ped).(rho_targets)
        p_values = interp1d(rho, profile_ped).(rho_targets)

        # Z's can be matched both by varying the value as well as the gradients
        # Here we want to keep the values as fixed as possible, while varying the gradients
        return diffnorm(z_targets, z_ped_values) / sum(abs, z_targets) + diffnorm(p_targets, p_values) / sum(abs, p_targets)
    end

    @assert 0.0 < ped_height "invalid ped_height = $(ped_height)"
    @assert 0.0 < ped_width < 1.0 "invalid ped_width = ($ped_width)"
    @assert rho[end] == 1.0

    z_profile = -calc_z(rho, profile, :backward)
    rho_targets = StaticArrays.@SVector[tr_bound0, tr_bound1]
    z_targets = interp1d(rho, z_profile).(rho_targets)
    p_targets = interp1d(rho, profile).(rho_targets)

    # figure out expin and expout such that the Z's of Hmode_profiles match the z_targets from transport
    expin = 1.0
    expout = 1.0
    x_guess = StaticArrays.@MVector[expin, expout]
    res = Optim.optimize(
        x -> cost_find_EPED_exps(x, ped_height, ped_width, rho, profile, p_targets, z_targets, rho_targets),
        x_guess,
        Optim.NelderMead(),
        Optim.Options(; g_tol=1E-6)
    )
    expin = abs(res.minimizer[1])
    expout = abs(res.minimizer[2])

    return blend_core_edge_EPED(profile, rho, ped_height, ped_width, tr_bound0, tr_bound1, expin, expout; method)
end

@compat public blend_core_edge_Hmode
push!(document[Symbol("Physics pedestal")], :blend_core_edge_Hmode)

"""
    blend_core_edge_EPED(
        profile::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        ped_height::Real,
        ped_width::Real,
        nml_bound::Real,
        ped_bound::Real;
        expin::Real,
        expout::Real;
        method::Symbol=:shift)

Blends the core and pedestal for given profile to match ped_height, ped_width using nml_bound and ped_bound as blending boundaries
"""
function blend_core_edge_EPED(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    ped_height::Real,
    ped_width::Real,
    nml_bound::Real,
    ped_bound::Real,
    expin::Real,
    expout::Real;
    method::Symbol=:shift
)
    # H-mode profile used for pedestal
    # NOTE: Note that we do not provide a core value as this causes the pedestal solution to depend on the core solution
    #       which breaks finding self-consistent core-pedestal solution through an optimizer
    profile_ped = ped_height_at_09(rho, Hmode_profiles(profile[end], ped_height, length(rho), expin, expout, ped_width), ped_height)
    return blend_core_edge(profile, profile_ped, rho, nml_bound, ped_bound; method)
end

@compat public blend_core_edge_EPED
push!(document[Symbol("Physics pedestal")], :blend_core_edge_EPED)

"""
    blend_core_edge(
        profile::AbstractVector{<:Real},
        profile_ped::AbstractVector{<:Real},
        rho::AbstractVector{<:Real},
        nml_bound::Real,
        ped_bound::Real;
        method::Symbol=:shift)

Blends core and edge profiles via inverse-scale-lengths method using `nml_bound` and `ped_bound` as blending boundaries

Different methods for connecting core region are:

  - z: inverse scale length
  - shift: add/subract constant to core region
  - scale: multiply/divide by a constant the core region
"""
function blend_core_edge(
    profile::AbstractVector{<:Real},
    profile_ped::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    nml_bound::Real,
    ped_bound::Real;
    method::Symbol=:shift
)
    @assert rho[end] == 1.0
    @assert nml_bound <= ped_bound "Unable to blend the core-pedestal because the nml_bound $nml_bound > ped_bound top $ped_bound"
    @assert length(profile) == length(profile_ped) == length(rho)
    iped = argmin_abs(rho, ped_bound)
    inml = argmin_abs(rho, nml_bound)

    z_profile = -calc_z(rho, profile, :backward)
    z_nml = z_profile[inml]

    # linear z between nml and pedestal
    if nml_bound < ped_bound
        z_profile_ped = -calc_z(rho, profile_ped, :backward)
        z_ped = z_profile_ped[iped]
        z_profile[inml:iped] = (z_nml - z_ped) ./ (rho[inml] - rho[iped]) .* (rho[inml:iped] .- rho[inml]) .+ z_nml
    end

    # integrate from pedestal inward
    profile_new = deepcopy(profile_ped)
    profile_new[inml:iped] .= integ_z(rho[inml:iped], z_profile[inml:iped], profile_ped[iped])

    # different blending strategies for the core
    @assert method in (:shift, :scale, :z)
    if method == :shift
        profile_new[inml:iped] .= integ_z(rho[inml:iped], z_profile[inml:iped], profile_ped[iped])
        profile_new[1:inml-1] = profile[1:inml-1] .- profile[inml] .+ profile_new[inml]
    elseif method == :scale
        profile_new[inml:iped] .= integ_z(rho[inml:iped], z_profile[inml:iped], profile_ped[iped])
        profile_new[1:inml-1] .= profile[1:inml-1] ./ profile[inml] .* profile_new[inml]
    elseif method == :z
        profile_new[1:iped] .= integ_z(rho[1:iped], z_profile[1:iped], profile_ped[iped])
    end

    return profile_new
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
    @assert rho[end] == 1.0
    return blend_core_edge_Lmode(profile, rho, ped_height, tr_bound1)
end

function blend_core_edge_Lmode(
    profile::AbstractVector{<:Real},
    rho::AbstractVector{<:Real},
    value::Real,
    rho_bound::Real
)
    @assert rho[end] == 1.0

    res = Optim.optimize(α -> cost_WPED_α!(rho, profile, α, value, rho_bound), -500, 500, Optim.Brent(); rel_tol=1E-3)
    cost_WPED_α!(rho, profile, res.minimizer, value, rho_bound)

    return profile
end

@compat public blend_core_edge_Lmode
push!(document[Symbol("Physics pedestal")], :blend_core_edge_Lmode)

function cost_WPED_α!(rho::AbstractVector{<:Real}, profile::AbstractVector{<:Real}, α::Real, value_ped::Real, rho_ped::Real)
    @assert rho[end] == 1.0

    rho_ped_idx = argmin_abs(rho, rho_ped)

    profile_ped = exponential_profile(rho, rho_ped, value_ped, profile[end], α)
    z_profile_ped = calc_z(rho, profile_ped, :backward)

    profile .+= (-profile[rho_ped_idx] + value_ped)
    z_profile = calc_z(rho, profile, :backward)

    profile[rho_ped_idx+1:end] .= interp1d(rho, profile_ped).(rho[rho_ped_idx+1:end])

    cost = abs.((z_profile[rho_ped_idx] - z_profile_ped[rho_ped_idx]) / z_profile[rho_ped_idx])
    return cost
end

"""
    pedestal_finder(profile::Vector{T}, psi_norm::Vector{T}; do_plot::Bool=false) where {T<:Real}

Finds the pedetal height and width using the EPED1 definition.

NOTE: The width is limited to be between 0.01 and 0.1.
If the width is at the 0.1 boundary it is likely an indication that the profile is not a typical H-mode profile.
The height is the value of the profile evaluated at (1.0 - width)
"""
function pedestal_finder(profile::Vector{T}, psi_norm::Vector{T}; do_plot::Bool=false, guess=nothing) where {T<:Real}
    @assert psi_norm[end] == 1.0

    psi_norm0 = range(0, 1, length(profile))
    profile0 = interp1d(psi_norm, profile).(psi_norm0)

    mask = psi_norm0

    function cost_function(params)
        width0 = mirror_bound(params[1], 0.01, 0.1)
        height0 = interp1d(psi_norm0, profile0)(1.0 - width0)
        core0 = abs(params[2])
        expin0 = abs(params[3]) + 1.0
        expout0 = abs(params[4]) + 1.0
        offset0 = mirror_bound(params[5], 0.0, width0 * 0.4)

        profile_fit0 = Hmode_profiles(profile0[end], height0, core0, length(profile0), expin0, expout0, width0 - offset0; offset=offset0)

        return norm(mask .* (profile_fit0 .- profile0))
    end

    if guess != nothing
        width = guess.width
        height = interp1d(psi_norm0, profile0)(1.0 - width)
        core = profile[1]
        expin = guess.expin
        expout = guess.expout
        offset = guess.offset
    else
        width = (1.0 - psi_norm0[argmax(psi_norm0 .* abs.(gradient(psi_norm0, profile0)))]) * 2.0
        height = interp1d(psi_norm0, profile0)(1.0 - width)
        core = profile[1]
        expin = 1.0
        expout = 1.0
        offset = 0.0
    end

    res = Optim.optimize(params -> cost_function(params), [width, core, expin - 1.0, expout - 1.0, offset], Optim.NelderMead())

    width = mirror_bound(res.minimizer[1], 0.01, 0.1)
    height = interp1d(psi_norm0, profile0)(1.0 - width)
    core = abs(res.minimizer[2])
    expin = 1.0 + abs(res.minimizer[3])
    expout = 1.0 + abs(res.minimizer[4])
    offset = mirror_bound(res.minimizer[1], 0.0, width * 0.4)

    if do_plot
        profile_fit = Hmode_profiles(profile[end], height, core, length(profile), expin, expout, width - offset; offset)
        p = plot(psi_norm0, profile0; label="profile", marker=:circle, markersize=1)
        plot!(p, psi_norm0, profile_fit; label="fit")
        hline!(p, [height]; ls=:dash, primary=false)
        vline!(p, [1.0 .- width]; ls=:dash, primary=false)
        scatter!(p, [1.0 - width], [height]; primary=false)
        display(p)
    end

    return (height=height, width=width, expin=expin, expout=expout, offset=offset, core=core, edge=profile[end])
end

@compat public pedestal_finder
push!(document[Symbol("Physics pedestal")], :pedestal_finder)

"""
    ped_height_at_09(rho::AbstractVector{T}, profile::AbstractVector{T}, height09::T) where {T<:Real}

Scale density profile so that value at rho=0.9 is height09
"""
function ped_height_at_09(rho::AbstractVector{T}, profile::AbstractVector{T}, ped_height::T) where {T<:Real}
    return profile / interp1d(rho, profile).(0.9) * ped_height
end

@compat public ped_height_at_09
push!(document[Symbol("Physics pedestal")], :ped_height_at_09)

"""
    pedestal_tanh_width_half_maximum(rho::AbstractVector{T}, profile::AbstractVector{T}; rho_pedestal_full_height::Float64=0.9, tanh_width_to_09_factor::Float64=0.85)

Estimates the tanh-width at half-maximum of a pedestal profile, based on a hyperbolic tangent fit

  - `rho_pedestal_full_height`: The normalized flux coordinate at which the pedestal is considered to reach full height.
  - `tanh_width_to_09_factor`: Factor to convert the pedestal width at full height (ρ = 0.9) to a hyperbolic tangent profile width.
"""
function pedestal_tanh_width_half_maximum(rho::AbstractVector{T}, profile::AbstractVector{T}; rho_pedestal_full_height::Float64=0.9, tanh_width_to_09_factor::Float64=0.85) where {T<:Real}
    index = rho .>= rho_pedestal_full_height
    profile09 = (IMAS.interp1d(rho, profile)(rho_pedestal_full_height) - profile[end]) * tanh_width_to_09_factor + profile[end]
    profile10 = profile[end]
    profilex = (profile10 + profile09) / 2
    tanh_width = (1.0 - IMAS.intersection(rho[index], profile[index], [0.0, 1.0], [profilex, profilex]).crossings[end][1]) * 2
    return tanh_width
end

@compat public pedestal_tanh_width_half_maximum
push!(document[Symbol("Physics pedestal")], :pedestal_tanh_width_half_maximum)

"""
    h_mode_detector(rho::AbstractVector{T}, electrons_pressure::AbstractVector{T}; threshold::Float64=0.4) where {T<:Real}

Given a profile (works well with electron pressure) it identifies the presence of a pedestal.

This function works by comparing the inverse scalelength at the pedestal
(defined as where the inverse scalelength is maximum)
against the inverse scalelength at the top of the pedestal.
"""
function h_mode_detector(rho::AbstractVector{T}, electrons_pressure::AbstractVector{T}; threshold::Float64=0.4) where {T<:Real}
    p = electrons_pressure
    n = length(p)
    v = p .+ sum(p) / n * 0.1 # let the inverse scalelength be large because the gradients are high, not because of very small values
    z = -IMAS.calc_z(rho, v, :backward)

    imaxZ = argmax(z)
    maxz = z[imaxZ]
    rho0 = rho[imaxZ]

    rhoσ = 1.0 - rho0
    zwell = IMAS.interp1d(rho, z).(rho0 - 2 * rhoσ)

    if rho0 < 1.0 && zwell / maxz < threshold
        hmode = true
    else
        hmode = false
    end

    # plot(rho, p / maximum(p); label="", ylim=(0, 1), color=hmode ? :red : :blue)
    # plot!(rho, z / maxz; label="", ylim=(0, 1), color=:black)
    # return hline!([zwell]; label="")

    return hmode
end
