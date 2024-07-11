"""
    elongation_limit(R0_over_a::Real)

Returns elongation limit due to control limit from simple aspect ratio scaling
"""
function elongation_limit(R0_over_a::Real)
    return 2.43 + 65.0 * exp(-R0_over_a / 0.376)
end

"""
    elongation_limit(eqt::IMAS.equilibrium__time_slice)

Returns elongation limit due to control limit from simple aspect ratio scaling
"""
function elongation_limit(eqt::IMAS.equilibrium__time_slice)
    return elongation_limit(eqt.global_quantities.magnetic_axis.r / eqt.boundary.minor_radius)
end

"""
    optimal_kappa_delta(li::T1, βp::T1, ϵ::T1, γτw::T2, ∆o::T2) where {T1<:Real,T2<:Real}

An analytic scaling relation for the maximum tokamak elongation against n=0 MHD resistive wall modes
Jungpyo Lee, Jeffrey P. Freidberg, Antoine J. Cerfon, Martin Greenwald
https://doi.org/10.1088%2F1741-4326%2Faa6877

NOTE:

  - γτw is the feedback capability parameter and represents how fast a instability is controllable (𝛾 is the instability growth rate and τw is the wall diffusion time)

  - ∆o is the outer gap (NOTE: assumes ∆o = ∆i = 1/3 * ∆v) detemines the relation between κ and δ of the plasma boundary and the κw=(κ+3∆o)(1+∆o) and δw=δ(1+∆o) of the wall boundary
"""
function optimal_kappa_delta(li::T1, βp::T1, ϵ::T1, γτw::T2, ∆o::T2) where {T1<:Real,T2<:Real}
    δδ_ = [0.0, 0.33, 0.50, 0.70]

    k̂0_ = Float64[0.54, 0.54, 0.55, 0.63]
    ν2_ = Float64[-0.68, -0.47, -0.08, 1.20] # li
    ν1_ = Float64[0.62, 0.71, 0.82, 1.14] # γτw
    ν3_ = Float64[-3.52, -4.00, -4.74, -6.67] # (1.0 + ∆o)
    k0_ = Float64[]
    for (k̂0, ν2, ν1, ν3) in zip(k̂0_, ν2_, ν1_, ν3_)
        push!(k0_, 1.0 + k̂0 * li^ν2 * γτw^ν1 * (1.0 + ∆o)^ν3)
    end

    k̂1_ = Float64[0.04, 0.35, 0.41, 0.52]
    μ1_ = Float64[-6.98, -1.42, -1.21, -2.00] # li
    μ2_ = Float64[-2.67, -0.04, 0.06, 0.17] # βp
    μ3_ = Float64[-1.47, -0.27, -0.18, -0.50] # γτw
    μ4_ = Float64[1.84, 0.42, 0.68, 2.32] # (1.0 + ∆o)
    k1_ = Float64[]
    for (k̂1, μ1, μ2, μ3, μ4) in zip(k̂1_, μ1_, μ2_, μ3_, μ4_)
        push!(k1_, k̂1 * li^μ1 * βp^μ2 * γτw^μ3 * (1.0 + ∆o)^μ4)
    end

    δ_opt = 2.30 * li^1.27 * βp^-0.01 * ϵ^(1.21 − 0.76 * li - 1.22 * βp - 0.001 * γτw + 1.21 * (1.0 + ∆o))
    δ_opt = max(min(δ_opt, maximum(δδ_)), minimum(δδ_))

    k0 = interp1d(δδ_, k0_, :cubic).(δ_opt)
    k1 = interp1d(δδ_, k1_, :cubic).(δ_opt)
    k_max = k0 + k1 * ((2.0 * ϵ) / (1.0 + ϵ^2))^2

    return k_max, δ_opt
end

function optimal_kappa_delta(eqt::IMAS.equilibrium__time_slice, γτw::T, ∆o::T) where {T<:Real}
    li = eqt.global_quantities.li_3
    βp = eqt.global_quantities.beta_pol
    ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    return optimal_kappa_delta(li, βp, ϵ, γτw, ∆o)
end

"""
    A_effective(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

A_effective towards L to H scaling see G. Birkenmeier et al 2022 Nucl. Fusion 62 086005
"""
function A_effective(cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}
    numerator = T[]
    denominator = T[]
    for ion in cp1d.ion
        if ion.element[1].z_n == 1
            n_int = trapz(cp1d.grid.volume, ion.density)
            push!(numerator, n_int * ion.element[1].a)
            push!(denominator, n_int)
        end
    end
    return reduce(+, numerator) / reduce(+, denominator)
end

"""
    scaling_L_to_H_power(A_effective::Real, ne_volume::Real, B0::Real, surface_area::Real)

L to H transition power scaling for metal walls and isotope effect according to : G. Birkenmeier et al 2022 Nucl. Fusion 62 086005

inputs in SI and returns power in W
"""
function scaling_L_to_H_power(A_effective::Real, ne_volume::Real, B0::Real, surface_area::Real)
    return 1e6 * 0.8 * 2.0 / A_effective * 0.049 * (ne_volume / 1e20)^0.72 * abs(B0)^0.8 * surface_area^0.94
end

function scaling_L_to_H_power(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    B0 = B0_geo(eqt)
    return scaling_L_to_H_power(
        A_effective(cp1d),
        trapz(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end],
        B0,
        eqt.profiles_1d.surface[end]
    )
end

function scaling_L_to_H_power(dd::IMAS.dd)
    return scaling_L_to_H_power(dd.core_profiles.profiles_1d[], dd.equilibrium.time_slice[])
end

function L_H_threshold(dd::IMAS.dd)
    return L_H_threshold(dd.core_sources, dd.core_profiles.profiles_1d[], dd.equilibrium.time_slice[])
end

"""
    L_H_threshold(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Returns ratio of Psol to Plh
"""
function L_H_threshold(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    Psol = power_sol(cs, cp1d)
    Plh = scaling_L_to_H_power(cp1d, eqt)
    return Psol / Plh
end

"""
    satisfies_h_mode_conditions(dd::IMAS.dd)

Returns true if the plasma is diverted, has positive triangularity and Psol>Plh
"""
function satisfies_h_mode_conditions(dd::IMAS.dd)
    eqt = dd.equilibrium.time_slice[]
    diverted = length(eqt.boundary.x_point) > 0
    Psol_gt_Plh = IMAS.L_H_threshold(dd) > getproperty(dd.requirements, :lh_power_threshold_fraction, 1.0)
    positive_triangularity = eqt.boundary.triangularity > 0.0
    if Psol_gt_Plh && diverted && positive_triangularity
        return true
    else
        return false
    end
end
