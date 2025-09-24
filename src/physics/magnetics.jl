document[Symbol("Physics magnetics")] = Symbol[]

"""
    magnetics!(mag::IMAS.magnetics, eq::IMAS.equilibrium)

Calculates synthetic magnetic probes and flux loops data at all time slices where equilibrium is available
"""
function magnetics!(mag::IMAS.magnetics, eq::IMAS.equilibrium)
    for leaf in IMASdd.AbstractTrees.Leaves(mag)
        if leaf.field ∈ (:data_σ, :data, :time)
            empty!(leaf.ids, leaf.field)
        end
    end
    for eqt in eq.time_slice
        magnetics!(mag, eqt)
    end
end

"""
    magnetics!(mag::IMAS.magnetics, eqt::IMAS.equilibrium__time_slice)

Calculates synthetic magnetic probes and flux loops data for a single time slice
"""
function magnetics!(mag::IMAS.magnetics, eqt::IMAS.equilibrium__time_slice)
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt.profiles_2d)
    field!(mag.b_field_pol_probe, PSI_interpolant, eqt.time)
    return flux!(mag.flux_loop, PSI_interpolant, eqt.time)
end

"""
    field!(probes::IMAS.IDSvector{b_field_pol_probe{T}}, PSI_interpolant::Interpolations.AbstractInterpolation, time0::Float64) where {T<:Real

Calculates synthetic magnetic probes data
"""
function field!(probes::IMAS.IDSvector{IMAS.magnetics__b_field_pol_probe{T}}, PSI_interpolant::Interpolations.AbstractInterpolation, time0::Float64) where {T<:Real}
    for probe in probes
        # local Br and Bz at probe location
        Br, Bz = IMAS.Br_Bz(PSI_interpolant, probe.position.r, probe.position.z)
        #component of polodial field measured by the probe (poloidal angle is clock-wise angle from the horizontal)
        B = Br .* cos.(probe.poloidal_angle) .- Bz .* sin.(probe.poloidal_angle)
        set_time_array(probe.field, :data, time0, B)
    end
end

"""
    flux!(loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}, PSI_interpolant::Interpolations.AbstractInterpolation, time0::Float64) where {T <: Real}

Calculates synthetic flux loops data
"""
function flux!(loops::IMAS.IDSvector{IMAS.magnetics__flux_loop{T}}, PSI_interpolant::Interpolations.AbstractInterpolation, time0::Float64) where {T<:Real}
    for loop in loops
        r = sum(pos.r for pos in loop.position) / length(loop.position)
        z = sum(pos.z for pos in loop.position) / length(loop.position)
        psi = PSI_interpolant.(r, z)
        set_time_array(loop.flux, :data, time0, psi)
    end
end
