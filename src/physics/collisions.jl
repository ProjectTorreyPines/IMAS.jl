"""
    lnΛ_ee(ne, Te)

Calculate Couloumb logarithm for thermal electron-electron collisions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:return lnΛ: Coloumb logarithm
"""
function lnΛ_ee(ne::S, Te::T) where {S<:Real,T<:Real}
    ne *= 1e-6 # cm^-3
    return 23.5 - log(sqrt(ne) * (Te^(-5 / 4))) - sqrt(1e-5 + ((log(Te) - 2)^2) / 16)
end

"""
    lnΛ_ei(ne, Te, ni::Vector, Ti::Vector, mi::Vector, Zi::Vector{Int})

Calculate Couloumb logarithm for thermal electron-ion collisions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:return lnΛ: list of Coloumb logarithms for each provided ion
"""
function lnΛ_ei(ne::S, Te::P, ni::Vector{Q}, Ti::Vector{R}, mi::Vector{T}, Zi::Vector{Int}) where
{S<:Real,P<:Real,Q<:Real,R<:Real,T<:Real}
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
    lnΛ_ei(ne, T)

Calculate Couloumb logarithm for thermal electron-ion collisions where Ti*me/mi < 10Zi^2 eV < Te [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperature [eV]

:return lnΛ: Coloumb logarithm
"""
function lnΛ_ei(ne::S, Te::P) where {S<:Real,P<:Real}
    ne *= 1e-6  #cm^-3
    return 24 - log(sqrt(ne) / Te)
end

"""
    lnΛ_ii(ne, Te, ni::Vector, Ti::Vector, mi::Vector, Zi::Vector{Int}; beta_D=nothing)

Calculate Couloumb logarithm for mixed thermal ion-ion collisions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperaturs [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:param beta_D: matrix of relative drift velocities between ion species v_D = beta_D*c

:return lnΛ: matrix of Coloumb logarithms for each provided ion
"""
function lnΛ_ii(ne::S, Te::P, ni::Vector{Q}, Ti::Vector{R}, mi::Vector{T}, Zi::Vector{Int}; beta_D::Union{Nothing,Matrix{<:Real}}=nothing) where
{S<:Real,P<:Real,Q<:Real,R<:Real,T<:Real}
    ni *= 1e-6 #cm^-3
    ne *= 1e-6 #cm^-3

    c = constants.c
    m_u = constants.m_u
    m_e = constants.m_e
    e = constants.e

    U = Te * e / m_e

    N = length(ni)
    type = promote_type(S, P, Q, R, T)

    if beta_D == nothing
        beta_D = zeros(type, N, N)
    else
        if size(beta_D) != (N, N)
            @error "Incorrect size of relative drift velocity matrix: beta_D($N,$N) expected."
        end
    end

    lnΛ = zeros(type, N, N)
    for i = 1:N
        ni1 = ni[i]
        Z1 = Zi[i]
        mu1 = mi[i]
        Ti1 = Ti[i]
        L1 = Ti1 * e / (mu1 * m_u)
        for j = 1:N
            ni2 = ni[j]
            Z2 = Zi[j]
            mu2 = mi[j]
            Ti2 = Ti[j]
            L2 = Ti2 * e / (mu2 * m_u)
            if max(L1, L2) < (beta_D[i, j] * c)^2 < U
                lnΛ[i, j] = 43 - log((Z1 * Z2 * (mu1 + mu2) / (mu1 * mu2 * beta_D[i, j]^2)) * sqrt(ne / Te))
            else
                lnΛ[i, j] = 23 - log(((Z1 * Z2 * (mu1 + mu2)) / (mu1 * Ti2 + mu2 * Ti1)) * sqrt((ni1 * Z1^2) / Ti1 + (ni2 * Z2^2) / Ti2))
            end
        end
    end
    return lnΛ
end

"""
    lnΛ_fi(ne, Te, ni::Vector, Ti::Vector, mi::Vector, Zi::Vector{Int}, beta_f, mf, Zf)

Calculate Couloumb logarithm for beam/fast ion in the presence of warm electrons/ions [NRL Plasma Formulary]

:param ne: electron density [m^-3]

:param Te: electron temperaturs [eV]

:param ni: list of ion densities [m^-3]

:param Ti: list of ion temperaturs [eV]

:param mi: list of ion masses [amu]

:param Zi: list of ion charges

:param beta_f: relative fast ion velocity v_f = beta_f*c

:param mf: mass of fast ion

:param Zf: charge of fast ion

:return lnΛ: list of Coloumb logarithms for each provided thermal ion species
"""
function lnΛ_fi(ne::S, Te::P, ni::Vector{Q}, Ti::Vector{R}, mi::Vector{O}, Zi::Vector{Int}, beta_f::T, mf::V, Zf::Int; verbose=true) where
{S<:Real,P<:Real,Q<:Real,R<:Real,O<:Real,T<:Real,V<:Real}
    ne *= 1e-6 #cm^-3

    c = constants.c
    m_u = constants.m_u
    m_e = constants.m_e
    e = constants.e

    U = Te * e / m_e

    m_f = constants.m_u * mf

    N = length(ni)
    lnΛ = zeros(promote_type(S, P, Q, R, O, T, V), N)
    for i = 1:N
        Z_i = Zi[i]
        m_i = mi[i] * m_u
        T_i = Ti[i]
        L = T_i * e / m_i
        if L < (beta_f * c)^2 < U
            lnΛ[i] = 43 - log((Zf * Z_i * (mf + mi[i]) / (mf * mi[i] * beta_f^2)) * sqrt(ne / Te))
        else
            verbose && @warn "Fast ion velocity outside of applicable range: $L < $((beta_f*c)^2) < $U. Assuming Ef=Ti"
            lnΛ[i] = lnΛ_ii(ne * 1e6, Te, [ni[i], ni[i]], [Ti[i], Ti[i]], [mf, mi[i]], [Zf, Zi[i]])[1, 2]
        end
    end
    return lnΛ
end
