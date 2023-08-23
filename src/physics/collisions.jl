"""
    lnΛ_ee(ne::Real, Te::Real)

Calculate Couloumb logarithm for thermal electron-electron collisions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:return lnΛ: Coloumb logarithm
"""
function lnΛ_ee(ne::Real, Te::Real)
    ne *= 1e-6 # cm^-3
    return 23.5 - log(sqrt(ne) * (Te^(-5 / 4))) - sqrt(1e-5 + ((log(Te) - 2)^2) / 16)
end

"""
    nΛ_ei(ne::S, Te::P, ni::AbstractVector{Q}, Ti::AbstractVector{R}, mi::AbstractVector{T}, Zi::AbstractVector{Int}) where {S<:Real,P<:Real,Q<:Real,R<:Real,T<:Real}

Calculate Couloumb logarithm for thermal electron-ion collisions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:return lnΛ: list of Coloumb logarithms for each provided ion
"""
function lnΛ_ei(ne::S, Te::P, ni::AbstractVector{Q}, Ti::AbstractVector{R}, mi::AbstractVector{T}, Zi::AbstractVector{Int}) where {S<:Real,P<:Real,Q<:Real,R<:Real,T<:Real}
    ne *= 1e-6  #cm^-3
    ni *= 1e-6 #cm^-3

    lnΛ = zeros(promote_type(S, P, Q, R, T), ni)
    m_e = constants.m_e / constants.m_u

    for i = 1:length(ni)
        mr = m_e / mi[i]
        if Ti[i] * mr < Te < 10Zi[i]^2
            lnΛ[i] = 23 - log(sqrt(ne) * Zi[i] * Te^(-3 / 2))
        elseif Ti[i] * mr < 10Zi[i]^22 < Te
            lnΛ[i] = 24 - log(sqrt(ne) / Te)
        elseif Te < Ti[i] * mr
            lnΛ[i] = 16 - log(mi[i] * sqrt(ni) * (Ti[i]^(-3 / 2)) * Zi[i]^2)
        else
            lnΛ[i] = 16.0
        end
    end

    return lnΛ
end

"""
    lnΛ_ei(ne::Real, Te::Real)

Calculate Couloumb logarithm for thermal electron-ion collisions where Ti*me/mi < 10Zi^2 eV < Te [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:return lnΛ: Coloumb logarithm
"""
function lnΛ_ei(ne::Real, Te::Real)
    ne *= 1e-6  #cm^-3
    return 24 - log(sqrt(ne) / Te)
end

"""
    lnΛ_ii(ne::Real, Te::Real, ni1::Real, Ti1::Real, mi1::Real, Zi1::Int, ni2::Real, Ti2::Real, mi2::Real, Zi2::Int; beta_D::Real=0.0)

Calculate Couloumb logarithm for mixed thermal ion1-ion2 collisions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperaturs [eV]

:param ni1: ion1 density [m^-3]

:param Ti1: ion1 temperature [eV]

:param mi1: ion1 mass [amu]

:param Zi1: ion1 charge

:param ni1: ion2 density [m^-3]

:param Ti1: ion2 temperature [eV]

:param mi1: ion2 mass [amu]

:param Zi1: ion2 charge

:param beta_D: relative drift velocities between ion species v_D = beta_D*c

:return lnΛ: Coloumb logarithm for ion1-ion2 collisions
"""
function lnΛ_ii(ne::Real, Te::Real, ni1::Real, Ti1::Real, mi1::Real, Zi1::Int, ni2::Real, Ti2::Real, mi2::Real, Zi2::Int; beta_D::Real=0.0)
    ni1 *= 1e-6 #cm^-3
    ni2 *= 1e-6 #cm^-3
    ne *= 1e-6 #cm^-3

    c = constants.c
    m_u = constants.m_u
    m_e = constants.m_e
    e = constants.e

    U = Te * e / m_e

    L1 = Ti1 * e / (mi1 * m_u)
    L2 = Ti2 * e / (mi2 * m_u)
    if max(L1, L2) < (beta_D * c)^2 < U
        lnΛ = 43 - log((Zi1 * Zi2 * (mi1 + mi2) / (mi1 * mi2 * beta_D^2)) * sqrt(ne / Te))
    else
        lnΛ = 23 - log(((Zi1 * Zi2 * (mi1 + mi2)) / (mi1 * Ti2 + mi2 * Ti1)) * sqrt((ni1 * Zi1^2) / Ti1 + (ni2 * Zi2^2) / Ti2))
    end


    return lnΛ
end

"""
    lnΛ_fi(ne::Real, Te::Real, n_i::Real, T_i::Real, m_i::Real, Z_i::Int, beta_f::Real, mf::Real, Zf::Int; verbose=true)

Calculate Couloumb logarithm for beam/fast ion in the presence of warm electrons/ions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperaturs [eV]

:param n_i: ion density [m^-3]

:param Ti: ion temperature [eV]

:param mi: ion mass [amu]

:param Zi: ion charge

:param beta_f: relative fast ion velocity v_f = beta_f*c

:param mf: mass of fast ion

:param Zf: charge of fast ion

:return lnΛ: Coloumb logarithm for thermal ion species
"""
function lnΛ_fi(ne::Real, Te::Real, n_i::Real, T_i::Real, m_i::Real, Z_i::Int, beta_f::Real, mf::Real, Zf::Int; verbose=true)
    ne *= 1e-6 #cm^-3

    c = constants.c
    m_u = constants.m_u
    m_e = constants.m_e
    e = constants.e

    U = Te * e / m_e
    L = T_i * e / (m_i * m_u)
    if L < (beta_f * c)^2 < U
        lnΛ = 43 - log((Zf * Z_i * (mf + m_i) / (mf * m_i * beta_f^2)) * sqrt(ne / Te))
    else
        verbose && @warn "Fast ion velocity outside of applicable range: $L < $((beta_f*c)^2) < $U. Assuming Ef=Ti"
        lnΛ = lnΛ_ii(ne * 1e6, Te, n_i, T_i, mf, Zf, n_i, T_i, m_i, Z_i)
    end
    return lnΛ
end