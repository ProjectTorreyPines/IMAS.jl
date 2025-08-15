document[Symbol("Physics equilibrium")] = Symbol[]

"""
    ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)

Returns r, z, and ψ interpolant named tuple
"""
function ψ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d)
    dim1 = convert(Vector{Float64}, eqt2d.grid.dim1)
    dim2 = convert(Vector{Float64}, eqt2d.grid.dim2)
    r = range(dim1[1], dim1[end], length(dim1))
    z = range(dim2[1], dim2[end], length(dim2))
    return ψ_interpolant(r, z, eqt2d.psi)
end

"""
    ψ_interpolant(r::AbstractRange{T1}, z::AbstractRange{T1}, psi::Matrix{T2}) where {T1<:Real,T2<:Real}
"""
function ψ_interpolant(r::AbstractRange{T1}, z::AbstractRange{T1}, psi::Matrix{T2}) where {T1<:Real,T2<:Real}
    PSI_interpolant = Interpolations.cubic_spline_interpolation((r, z), psi; extrapolation_bc=Interpolations.Line())
    return (r=r, z=z, PSI_interpolant=PSI_interpolant)
end

"""
    ψ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
"""
function ψ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    eqt2d = findfirst(:rectangular, eqt2dv)
    return ψ_interpolant(eqt2d)
end

"""
    ψ_interpolant(eqt::IMAS.equilibrium__time_slice)
"""
function ψ_interpolant(eqt::IMAS.equilibrium__time_slice)
    return ψ_interpolant(eqt.profiles_2d)
end

@compat public ψ_interpolant
push!(document[Symbol("Physics equilibrium")], :ψ_interpolant)

"""
    ρ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d{T}, phi_norm::T) where {T<:Real}

Returns r, z, and ρ interpolant named tuple
"""
function ρ_interpolant(eqt2d::IMAS.equilibrium__time_slice___profiles_2d{T}, phi_norm::T) where {T<:Real}
    r = range(eqt2d.grid.dim1[1], eqt2d.grid.dim1[end], length(eqt2d.grid.dim1))
    z = range(eqt2d.grid.dim2[1], eqt2d.grid.dim2[end], length(eqt2d.grid.dim2))
    return ρ_interpolant(r, z, sqrt.(abs.(eqt2d.phi ./ phi_norm)))
end

"""
    ρ_interpolant(r::AbstractRange{T}, z::AbstractRange{T}, rho::Matrix{T}) where {T<:Real}
"""
function ρ_interpolant(r::AbstractRange{T}, z::AbstractRange{T}, rho::Matrix{T}) where {T<:Real}
    RHO_interpolant = Interpolations.cubic_spline_interpolation((r, z), rho; extrapolation_bc=Interpolations.Line())
    return (r=r, z=z, RHO_interpolant=RHO_interpolant)
end

"""
    ρ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d{T}}, phi_norm::T) where {T<:Real}
"""
function ρ_interpolant(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d{T}}, phi_norm::T) where {T<:Real}
    eqt2d = findfirst(:rectangular, eqt2dv)
    return ρ_interpolant(eqt2d, phi_norm)
end

"""
    ρ_interpolant(eqt::IMAS.equilibrium__time_slice)
"""
function ρ_interpolant(eqt::IMAS.equilibrium__time_slice)
    return ρ_interpolant(eqt.profiles_2d, eqt.profiles_1d.phi[end])
end

@compat public ρ_interpolant
push!(document[Symbol("Physics equilibrium")], :ρ_interpolant)

"""
    calc_pprime_ffprim_f(
        psi::T,
        R::T,
        one_R::T,
        one_R2::T,
        R0::Real,
        B0::Real;
        pressure::Union{Missing,T}=missing,
        dpressure_dpsi::Union{Missing,T}=missing,
        j_tor::Union{Missing,T}=missing,
        j_tor_over_R::Union{Missing,T}=missing,
        f::Union{Missing,T}=missing,
        f_df_dpsi::Union{Missing,T}=missing,
        ) where {T<:AbstractVector{<:Real}}

This method returns the `P'` and `FF'` given `P` or `P'` and `J` or `J/R` based on the current equilibrium fluxsurfaces geometry

Arguments:
* `psi`: poloidal flux
* `R`: <R>
* `one_R`: <1/R>
* `one_R2`: <1/R²>
* `R0`: R at which B0 is defined
* `B0`: vacuum magnetic field
* `pressure`: P
* `dpressure_dpsi`: P'
* `j_tor`: toroidal current
* `j_tor_over_R`: flux surface averaged toroidal current density over major radius
* `f`: F
* `f_df_dpsi`: FF'

returns (P', FF', F)
"""
function calc_pprime_ffprim_f(
    psi::T,
    R::T,
    one_R::T,
    one_R2::T,
    R0::Real,
    B0::Real;
    pressure::Union{Missing,T}=missing,
    dpressure_dpsi::Union{Missing,T}=missing,
    j_tor::Union{Missing,T}=missing,
    j_tor_over_R::Union{Missing,T}=missing,
    f::Union{Missing,T}=missing,
    f_df_dpsi::Union{Missing,T}=missing
) where {T<:AbstractVector{<:Real}}

    COCOS = cocos(11)

    if pressure !== missing && dpressure_dpsi === missing
        dpressure_dpsi = gradient(psi, pressure)
    end

    if f !== missing && f_df_dpsi === missing
        f_df_dpsi = gradient(psi, f) .* f
    end

    if j_tor !== missing && f_df_dpsi === missing
        f_df_dpsi = j_tor .* COCOS.sigma_Bp ./ (2π)^COCOS.exp_Bp .+ dpressure_dpsi .* R
        f_df_dpsi .*= -mks.μ_0 ./ one_R
    elseif j_tor_over_R !== missing && f_df_dpsi === missing
        f_df_dpsi = j_tor_over_R * COCOS.sigma_Bp / (2π)^COCOS.exp_Bp + dpressure_dpsi
        f_df_dpsi .*= -mks.μ_0 ./ one_R2
    end

    if f === missing
        # calculate F from FF' with boundary condition so that f[end]=R0*B0
        f = 2 * cumtrapz(psi, f_df_dpsi)
        if minimum(f) < f[1] && minimum(f) < f[end]
            f .-= maximum(f)
        else
            f .-= minimum(f)
        end
        C = (R0 * B0)^2 - f[end]
        f = sqrt.(abs.(f .+ C))
    end

    return (dpressure_dpsi=dpressure_dpsi, f_df_dpsi=f_df_dpsi, f=f)
end

@compat public calc_pprime_ffprim_f
push!(document[Symbol("Physics equilibrium")], :calc_pprime_ffprim_f)

"""
    p_jtor_2_pprime_ffprim_f(eqt1::IMAS.equilibrium__time_slice___profiles_1d, R0::Real, B0::Real)

Calculate `P'`, `FF'` and `F` from `pressure` and `j_tor` in equilibrium
"""
function p_jtor_2_pprime_ffprim_f(eqt1::IMAS.equilibrium__time_slice___profiles_1d, R0::Real, B0::Real)
    psi = eqt1.psi
    R = eqt1.gm8
    one_R = eqt1.gm9
    one_R2 = eqt1.gm1

    presssure = getproperty(eqt1, :pressure, missing)
    j_tor = getproperty(eqt1, :j_tor, missing)

    return calc_pprime_ffprim_f(psi, R, one_R, one_R2, R0, B0; presssure, j_tor)
end

@compat public p_jtor_2_pprime_ffprim_f
push!(document[Symbol("Physics equilibrium")], :p_jtor_2_pprime_ffprim_f)

"""
    symmetrize_equilibrium!(eqt::IMAS.equilibrium__time_slice)

Update equilibrium time slice in place to be symmetric with respect to its magnetic axis.

This is done by averaging the upper and lower parts of the equilibrium.

Flux surfaces should re-traced after this operation.

NOTE: Use with care! This operation will change the flux surfaces (LCFS included) and as such quantities may change
"""
function symmetrize_equilibrium!(eqt::IMAS.equilibrium__time_slice)
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = ψ_interpolant(eqt2d)

    psi_sign = sign(eqt.profiles_1d.psi[end] - eqt.profiles_1d.psi[1])
    RA, ZA = find_magnetic_axis(r, z, PSI_interpolant, psi_sign)

    Z1 = (maximum(z) + minimum(z)) / 2.0
    zz = z .- Z1 .+ ZA
    zz = range(max(minimum(z), minimum(zz)), min(maximum(z), maximum(zz)), length(z))

    psi = PSI_interpolant(r, zz)

    eqt2d.grid.dim2 = zz
    eqt2d.psi = (psi[1:end, end:-1:1] .+ psi) ./ 2.0

    return eqt
end

@compat public symmetrize_equilibrium!
push!(document[Symbol("Physics equilibrium")], :symmetrize_equilibrium!)

"""
    B0_geo(eqt::IMAS.equilibrium__time_slice{T})::T where{T<:Real}

Returns vacuum B0 at the plasma geometric center, which may or may not be `eqt.global_quantities.vacuum_toroidal_field.r0`
"""
function B0_geo(eqt::IMAS.equilibrium__time_slice{T})::T where {T<:Real}
    r0 = eqt.global_quantities.vacuum_toroidal_field.r0
    b0 = eqt.global_quantities.vacuum_toroidal_field.b0
    Rgeo = eqt.boundary.geometric_axis.r
    return r0 .* b0 ./ Rgeo
end

@compat public B0_geo
push!(document[Symbol("Physics equilibrium")], :B0_geo)

"""
    elongation_limit(R0_over_a::Real)

Returns elongation limit due to control limit from simple aspect ratio scaling
"""
function elongation_limit(R0_over_a::Real)
    return 2.43 + 65.0 * exp(-R0_over_a / 0.376)
end

"""
    elongation_limit(eqt::IMAS.equilibrium__time_slice)
"""
function elongation_limit(eqt::IMAS.equilibrium__time_slice)
    return elongation_limit(eqt.global_quantities.magnetic_axis.r / eqt.boundary.minor_radius)
end

@compat public elongation_limit
push!(document[Symbol("Physics equilibrium")], :elongation_limit)


"""
    optimal_kappa_delta(li::T1, βp::T1, ϵ::T1, γτw::T2, ∆o::T2) where {T1<:Real,T2<:Real}

An analytic scaling relation for the maximum tokamak elongation against n=0 MHD resistive wall modes
Jungpyo Lee, Jeffrey P. Freidberg, Antoine J. Cerfon, Martin Greenwald https://doi.org/10.1088%2F1741-4326%2Faa6877

NOTE:

* `γτw` is the feedback capability parameter and represents how fast a instability is controllable (`𝛾` is the instability growth rate and `τw` is the wall diffusion time)
   (typically `γτw < 10` is assumed for controllability, see `VacuumFields.normalized_growth_rate()`)

* `∆o` is the outer gap (NOTE: assumes `∆o = ∆i = 1/3 * ∆v`) detemines the relation between `κ` and `δ` of the plasma boundary and the `κw=(κ+3∆o)(1+∆o)` and `δw=δ(1+∆o)` of the wall boundary
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

    k0 = cubic_interp1d(δδ_, k0_).(δ_opt)
    k1 = cubic_interp1d(δδ_, k1_).(δ_opt)
    κ_opt = k0 + k1 * ((2.0 * ϵ) / (1.0 + ϵ^2))^2

    return (κ_opt=κ_opt, δ_opt=δ_opt)
end

"""
    optimal_kappa_delta(eqt::IMAS.equilibrium__time_slice, γτw::T, ∆o::T) where {T<:Real}
"""
function optimal_kappa_delta(eqt::IMAS.equilibrium__time_slice, γτw::T, ∆o::T) where {T<:Real}
    li = eqt.global_quantities.li_3
    βp = eqt.global_quantities.beta_pol
    ϵ = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    return optimal_kappa_delta(li, βp, ϵ, γτw, ∆o)
end

@compat public optimal_kappa_delta
push!(document[Symbol("Physics equilibrium")], :optimal_kappa_delta)