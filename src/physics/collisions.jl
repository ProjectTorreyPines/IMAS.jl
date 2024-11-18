document[Symbol("Physics collisions")] = Symbol[]

"""
    lnΛ_ee(ne::Real, Te::Real)

Calculate Couloumb logarithm (lnΛ) for thermal electron-electron collisions [NRL Plasma Formulary]

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]
"""
function lnΛ_ee(ne::Real, Te::Real)
    ne *= 1e-6 # cm^-3
    return 23.5 - log(sqrt(ne) * (Te^(-5 / 4))) - sqrt(1e-5 + ((log(Te) - 2)^2) / 16)
end

@compat public lnΛ_ee
push!(document[Symbol("Physics collisions")], :lnΛ_ee)

"""
    nΛ_ei(ne::S, Te::P, ni::AbstractVector{Q}, Ti::AbstractVector{R}, mi::AbstractVector{T}, Zi::AbstractVector{Int}) where {S<:Real,P<:Real,Q<:Real,R<:Real,T<:Real}

Calculate Couloumb logarithm (lnΛ) for thermal electron-ion collisions [NRL Plasma Formulary]

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `ni`: list of ion densities [m^-3]

* `Ti`: list of ion temperaturs [eV]

* `mi`: list of ion masses [amu]

* `Zi`: list of ion charges
"""
function lnΛ_ei(ne::S, Te::P, ni::AbstractVector{Q}, Ti::AbstractVector{R}, mi::AbstractVector{T}, Zi::AbstractVector{Int}) where {S<:Real,P<:Real,Q<:Real,R<:Real,T<:Real}
    ne *= 1e-6  #cm^-3
    ni *= 1e-6 #cm^-3

    lnΛ = zeros(promote_type(S, P, Q, R, T), ni)
    m_e = mks.m_e / mks.m_u

    for i in 1:length(ni)
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

Calculate Couloumb logarithm (lnΛ) for thermal electron-ion collisions where Ti*me/mi < 10Zi^2 eV < Te [NRL Plasma Formulary]

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]
"""
function lnΛ_ei(ne::Real, Te::Real)
    ne *= 1e-6  #cm^-3
    return 24 - log(sqrt(ne) / Te)
end

@compat public lnΛ_ei
push!(document[Symbol("Physics collisions")], :lnΛ_ei)

"""
    lnΛ_ii(ne::Real, Te::Real, ni1::Real, Ti1::Real, mi1::Real, Zi1::Int, ni2::Real, Ti2::Real, mi2::Real, Zi2::Int; beta_D::Real=0.0)

Calculate Couloumb logarithm (lnΛ) for mixed thermal ion1-ion2 collisions [NRL Plasma Formulary]

* `ne`: electron density [m^-3]

* `Te`: electron temperaturs [eV]

* `ni1`: ion1 density [m^-3]

* `Ti1`: ion1 temperature [eV]

* `mi1`: ion1 mass [amu]

* `Zi1`: ion1 charge

* `ni1`: ion2 density [m^-3]

* `Ti1`: ion2 temperature [eV]

* `mi1`: ion2 mass [amu]

* `Zi1`: ion2 charge

* `beta_D`: relative drift velocities between ion species `v_D = beta_D*c`
"""
function lnΛ_ii(ne::Real, Te::Real, ni1::Real, Ti1::Real, mi1::Real, Zi1::Int, ni2::Real, Ti2::Real, mi2::Real, Zi2::Int; beta_D::Real=0.0)
    ni1 *= 1e-6 #cm^-3
    ni2 *= 1e-6 #cm^-3
    ne *= 1e-6 #cm^-3

    c = mks.c
    m_u = mks.m_u
    m_e = mks.m_e
    e = mks.e

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

@compat public lnΛ_ii
push!(document[Symbol("Physics collisions")], :lnΛ_ii)

"""
    lnΛ_fi(ne::Real, Te::Real, n_i::Real, T_i::Real, m_i::Real, Z_i::Int, beta_f::Real, mf::Real, Zf::Int; verbose=true)

Calculate Couloumb logarithm (lnΛ) for beam/fast ion in the presence of warm electrons/ions [NRL Plasma Formulary]

* `ne`: electron density [m^-3]

* `Te`: electron temperaturs [eV]

* `n_i`: ion density [m^-3]

* `Ti`: ion temperature [eV]

* `mi`: ion mass [amu]

* `Zi`: ion charge

* `beta_f`: relative fast ion velocity `v_f = beta_f*c`

* `mf`: mass of fast ion [amu]

* `Zf`: charge of fast ion
"""
function lnΛ_fi(ne::Real, Te::Real, n_i::Real, T_i::Real, m_i::Real, Z_i::Int, beta_f::Real, mf::Real, Zf::Int; verbose=true)
    ne *= 1e-6 #cm^-3

    c = mks.c
    m_u = mks.m_u
    m_e = mks.m_e
    e = mks.e

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

@compat public lnΛ_fi
push!(document[Symbol("Physics collisions")], :lnΛ_fi)