"""
    ω_pe(ne::Real)

Returns electron plasma frequency [rad/s] given electron density in m⁻³
"""
function ω_pe(ne::Real)
    return sqrt(ne * mks.e^2 / (mks.ϵ_0 * mks.m_e))
end

"""
    ω_ce(B::Real)

Returns electron cyclotron frequency [rad/s] given magnetic field B in T
"""
function ω_ce(B::Real)
    return mks.e * abs(B) / mks.m_e
end

"""
    B_ω_ce(ω::Real)

Returns magnetic field B in T for a given electron cyclotron frequency [rad/s]
"""
function B_ω_ce(ω::Real)
    return ω / mks.e * mks.m_e
end

"""
    ω_ci(B::Real, Z::Real, A::Real)

Returns ion cyclotron frequency [rad/s] given magnetic field B in T and the ion charge and mass in amu
"""
function ω_ci(B::Real, Z::Real, A::Real)
    return mks.e * abs(B) * Z / A * mks.m_p
end

"""
    stix_P(ω::Real, ne::Real)

P (Plasma term) of the Stix dielectric tensor
"""
function stix_P(ω::Real, ne::Real)
    return 1.0 - (w_p(ne) / ω)^2
end

"""
    stix_S(ω::Real, ne::Real, B::Real)

S (Sum term) of the Stix dielectric tensor
"""
function stix_S(ω::Real, ne::Real, B::Real)
    return 1.0 - ω_p(ne)^2 / (ω^2 - ω_ce(B)^2)
end

"""
    stix_D(ω::Real, ne::Real, B::Real)

D (Difference term) of the Stix dielectric tensor
"""
function stix_D(ω::Real, ne::Real, B::Real)
    return 1.0 - ω_ce(B) * ω_p(ne)^2 / (ω * (ω^2 - ω_ce(B)^2))
end