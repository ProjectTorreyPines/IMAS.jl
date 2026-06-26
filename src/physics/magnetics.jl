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
    if PSI_interpolant !== nothing
        field!(mag.b_field_pol_probe, PSI_interpolant, eqt.time)
        return flux!(mag.flux_loop, PSI_interpolant, eqt.time)
    end
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

"""
    sortperm(probes::AbstractVector{<:IMAS.magnetics__b_field_pol_probe})

Sort field probes around the vessel
"""
function Base.sortperm(probes::AbstractVector{<:IMAS.magnetics__b_field_pol_probe})
    R0 = sum([probe.position.r for probe in probes]) / length(probes)
    Z0 = sum([probe.position.z for probe in probes]) / length(probes)
    r = [probe.position.r - R0 for probe in probes]
    z = [probe.position.z - Z0 for probe in probes]
    return sortperm(atan.(z, r))
end

"""
    sortperm(probes::AbstractVector{<:IMAS.magnetics__flux_loop})

Sort flux loops around the vessel
"""
function Base.sortperm(loops::AbstractVector{<:IMAS.magnetics__flux_loop})
    R0 = sum([loop.position[1].r for loop in loops]) / length(loops)
    Z0 = sum([loop.position[1].z for loop in loops]) / length(loops)
    r = [loop.position[1].r - R0 for loop in loops]
    z = [loop.position[1].z - Z0 for loop in loops]
    return sortperm(atan.(z, r))
end
