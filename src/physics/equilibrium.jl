"""
    calc_pprime_ffprim_f(
        psi::T,
        R::T,
        one_R::T,
        one_R2::T,
        R0::Real,
        B0::Real;
        press::Union{Missing, T}=missing,
        pprime::Union{Missing, T}=missing,
        jtor::Union{Missing, T}=missing,
        jtor_over_R::Union{Missing, T}=missing,
        fpol::Union{Missing, T}=missing) where {T<:AbstractVector{<:Real}}

This method returns the P' and FF' given P or P' and J or J/R based on the current equilibrium fluxsurfaces geometry

Arguments:
psi::T,

  - R: <R>
  - one_R: <1/R>
  - one_R2: <1/R²>
  - R0: R at which B0 is defined
  - B0: vacuum magnetic field
  - press: pressure
  - pprime: pressure * pressure'
  - jtor: toroidal current
  - jtor_over_R: flux surface averaged toroidal current density over major radius
  - fpol: F

returns (P', FF', F)
"""
function calc_pprime_ffprim_f(
    psi::T,
    R::T,
    one_R::T,
    one_R2::T,
    R0::Real,
    B0::Real;
    press::Union{Missing,T}=missing,
    pprime::Union{Missing,T}=missing,
    jtor::Union{Missing,T}=missing,
    jtor_over_R::Union{Missing,T}=missing,
    fpol::Union{Missing,T}=missing) where {T<:AbstractVector{<:Real}}

    COCOS = cocos(11)

    if press !== missing
        pprime = gradient(psi, press)
    end

    if fpol !== missing
        ffprim = gradient(psi, fpol) .* fpol
    end

    if jtor !== missing
        ffprim = jtor .* COCOS.sigma_Bp ./ (2π)^COCOS.exp_Bp .+ pprime .* R
        ffprim .*= -constants.μ_0 ./ one_R
    elseif jtor_over_R !== missing
        ffprim = jtor_over_R * COCOS.sigma_Bp / (2π)^COCOS.exp_Bp + pprime
        ffprim .*= -constants.μ_0 ./ one_R2
    end

    # calculate F from FF' with boundary condition so that f[end]=R0*B0
    f = 2 * cumtrapz(psi, ffprim)
    if minimum(f) < f[1] && minimum(f) < f[end]
        f .-= maximum(f)
    else
        f .-= minimum(f)
    end
    C = (R0 * B0)^2 - f[end]
    f = sqrt.(abs.(f .+ C))

    return pprime, ffprim, f
end

function p_jtor_2_pprime_ffprim_f(eqt1::IMAS.equilibrium__time_slice___profiles_1d, R0::Real, B0::Real)
    psi = eqt1.psi
    R = eqt1.gm8
    one_R = eqt1.gm9
    one_R2 = eqt1.gm1

    press = getproperty(eqt1, :pressure, missing)
    jtor = getproperty(eqt1, :j_tor, missing)

    return calc_pprime_ffprim_f(psi, R, one_R, one_R2, R0, B0; press, jtor)
end

function p_jtor_2_pprime_ffprim_f!(eqt1::IMAS.equilibrium__time_slice___profiles_1d, R0::Real, B0::Real)
    pprime, ffprim, f = p_jtor_2_pprime_ffprim_f(eqt1, R0, B0)
    eqt1.dpressure_dpsi = pprime
    eqt1.f_df_dpsi = ffprim
    eqt1.f = f
    return pprime, ffprim, f
end

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
    R0, Z0 = find_magnetic_axis(r, z, PSI_interpolant, psi_sign)

    Z1 = (maximum(z) + minimum(z)) / 2.0
    zz = z .- Z1 .+ Z0
    zz = range(max(minimum(z), minimum(zz)), min(maximum(z), maximum(zz)), length(z))

    psi = PSI_interpolant(r, zz)

    eqt2d.grid.dim2 = zz
    eqt2d.psi = (psi[1:end, end:-1:1] .+ psi) ./ 2.0

    return eqt
end

"""
    B0_geo(eqt::IMAS.equilibrium__time_slice{T})::T where{T<:Real}

Returns vacuum B0 at the plasma geometric center
"""
function B0_geo(eqt::IMAS.equilibrium__time_slice{T})::T where {T<:Real}
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0
    Rgeo = eqt.boundary.geometric_axis.r
    return R0 .* B0 ./ Rgeo
end
