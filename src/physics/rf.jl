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

"""
    frequency(beam::IMAS.ec_launchers__beam)

Returns operating frequency of a given Gyrotron in Hz
"""
function frequency(beam::IMAS.ec_launchers__beam)
    return frequency(beam.frequency)
end

function frequency(beam_freq::IMAS.ec_launchers__beam___frequency)
    @assert minimum(beam_freq.data) == maximum(beam_freq.data)
    return maximum(beam_freq.data)
end

"""
    ech_resonance_layer(eqt::IMAS.equilibrium__time_slice{T}, frequency::T; harmonics::AbstractVector{Int}=[1, 2, 3], only_vacuum_Bt::Bool=false) where {T<:Real}

Returns named tuple with (r, z, harmonic) of ECH resonance layer closest to the plasma axis
"""
function ech_resonance_layer(eqt::IMAS.equilibrium__time_slice{T}, frequency::T; harmonics::AbstractVector{Int}=[1, 2, 3], only_vacuum_Bt::Bool=false) where {T<:Real}
    fundamental_B_ec_resonance = B_ω_ce(frequency * 2π)

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    R = eqt2d.grid.dim1
    Z = eqt2d.grid.dim2
    if only_vacuum_Bt
        RR, ZZ = meshgrid(R, Z)
        vacuum_Bt = abs(eqt.global_quantities.vacuum_toroidal_field.b0 * eqt.global_quantities.vacuum_toroidal_field.r0) ./ RR'
        Bωce_over_Btot = fundamental_B_ec_resonance ./ vacuum_Bt
    else
        Br, Bz = Br_Bz(eqt2d)
        Bωce_over_Btot = fundamental_B_ec_resonance ./ sqrt.(eqt2d.b_field_tor .^ 2 .+ Br .^ 2 .+ Bz .^ 2)
    end

    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.r

    layers = Dict{T,Tuple{Vector{T},Vector{T},Int}}()
    for harmonic in harmonics
        for line in Contour.lines(Contour.contour(R, Z, Bωce_over_Btot, harmonic))
            pr, pz = Contour.coordinates(line)
            dl = minimum((pr .- RA) .^ 2 .+ (pz .- ZA) .^ 2)
            layers[dl] = (pr, pz, harmonic)
        end
    end

    if isempty(layers)
        return (r=T[], z=T[], harmonic=0)
    else
        core_layer = layers[minimum(keys(layers))]
        return (r=core_layer[1], z=core_layer[2], harmonic=core_layer[3])
    end
end

"""
    ech_resonance(eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}

Returns the ech_resonance desirable for a given equilibrium, that is the layer that resonates with the vacuum field at the magnetic axis
"""
function ech_resonance(eqt::IMAS.equilibrium__time_slice{T}) where {T<:Real}
    Bgeo = eqt.global_quantities.vacuum_toroidal_field.b0 .* eqt.global_quantities.vacuum_toroidal_field.r0 / eqt.global_quantities.magnetic_axis.r
    f_fundamental = ω_ce(Bgeo) / 2π
    if f_fundamental * 2 < 170E9
        return (frequency=f_fundamental * 2, harmonic=2, mode="X")
    else
        return (frequency=f_fundamental, harmonic=1, mode="O")
    end
end
