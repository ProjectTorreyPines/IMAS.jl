document[Symbol("Physics fast")] = Symbol[]

"""
    estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Int)

Returns solution to `i2` and `i4` integrals from  [Estrada et al.,  Phys of Plasm. 13, 112303 (2006)] Eq. 9 & 10

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `ni`: list of ion densities [m^-3]

* `Ti`: list of ion temperatures [eV]

* `mi`: list of ion masses [amu]

* `Zi`: list of ion charges

* `Ef`: fast ion energy [eV]

* `mf`: mass of fast ion [amu]

* `Zf`: fast ion charge
"""
function estrada_I_integrals(
    ne::Real,
    Te::Real,
    ni::AbstractVector{<:Real},
    Ti::AbstractVector{<:Real},
    mi::AbstractVector{<:Real},
    Zi::AbstractVector{Int},
    Ef::Real,
    mf::Real,
    Zf::Int
)
    Ec = critical_energy(ne, Te, ni, Ti, mi, Zi, mf, Zf)
    i2, i4 = estrada_I_integrals(Ec, Ef)
    return (i2=i2, i4=i4)
end

"""
    estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Int)

Returns solution to i2 and i4 integrals from  [Estrada et al.,  Phys of Plasm. 13, 112303 (2006)] Eq. 9 & 10

* `Ef`: fast ion energy [eV]

* `Ec`: critical energy [eV]
"""
function estrada_I_integrals(Ec::Real, Ef::Real)
    a = sqrt.(Ec ./ Ef)
    i2 = (1 / 3.0) .* log.((1 .+ a .^ 3) ./ (a .^ 3))
    i4 = 0.5 .- a .^ 2 .* ((1 / 6.0) .* log.((1 .- a .+ a .^ 2) ./ (1 .+ a) .^ 2) .+
                           1 ./ sqrt(3.0) .* (atan.((2 .- a) ./ (a .* sqrt(3.0))) .+ pi / 6))
    return (i2=i2, i4=i4)
end

@compat public estrada_I_integrals
push!(document[Symbol("Physics fast")], :estrada_I_integrals)

"""
    slowing_down_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int)

Returns the slowing down time `τ_s` in seconds [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `ni`: list of ion densities [m^-3]

* `Ti`: list of ion temperatures [eV]

* `mi`: list of ion masses [amu]

* `Zi`: list of ion charges

* `mf`: mass of fast ion [amu]

* `Zf`: fast ion charge
"""
function slowing_down_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int)
    lnΛ = lnΛ_ei(ne, Te, ni, Ti, mi, Zi)
    ne_cm3 = 1e-6 * ne
    τ_s = 6.27e8 * mf * (Te^1.5) ./ (ne_cm3 * lnΛ * Zf^2)
    return τ_s
end

"""
    slowing_down_time(ne::Real, Te::Real, mf::Real, Zf::Int)

Calculates the slowing down time `τ_s` in seconds for `Ti*me/mi < 10Zi^2 eV < Te` [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `mf`: mass of fast ion [amu]

* `Zf`: fast ion charge
"""
function slowing_down_time(ne::Real, Te::Real, mf::Real, Zf::Int)
    lnΛ = lnΛ_ei(ne, Te)
    ne_cm3 = 1e-6 * ne
    τ_s = 6.27e8 * mf * (Te^1.5) / (ne_cm3 * lnΛ * Zf^2)
    return τ_s
end

"""
    α_slowing_down_time(cp1d::IMAS.core_profiles__profiles_1d)

Returns the slowing down time in seconds of α particles evaluated on axis
"""
function α_slowing_down_time(cp1d::IMAS.core_profiles__profiles_1d)
    α = ion_properties(:α)
    return slowing_down_time(cp1d.electrons.density_thermal[1], cp1d.electrons.temperature[1], α.a, Int(α.z_n))
end

@compat public α_slowing_down_time
push!(document[Symbol("Physics fast")], :α_slowing_down_time)

"""
    _drag_coefficient(n::Real, Z::Int, mf::Real, Zf::Int, lnΛ::Real)

Returns drag coefficient `Γ` for a fast-ions interacting with a thermal species as defined Eq. 8 in [Gaffey, J. D. (1976). Energetic ion distribution resulting from neutral beam injection in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

* `n`: density of the thermal species [m^-3]

* `Z`: charge of the thermal species

* `mf`: fast-ion mass [amu]

* `Zf`: fast-ion charge

* `lnΛ`: Couloumb logarithm
"""
function _drag_coefficient(n::Real, Z::Int, mf::Real, Zf::Int, lnΛ::Real)
    n *= 1e-6 #cm^-3
    mf *= mks.m_u # kg
    Γ = (2 * pi * n * (mks.e)^4 * Z^2 * Zf^2 * lnΛ) / (mf^2)
    return Γ
end

"""
    _electron_ion_drag_difference(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real, mf::Real, Zf::Int)

Returns `ΔD` the difference of the electron and ion drag terms in the collision operator defined in Eq. 19 in [Gaffey, J.D (1976). Energetic ion distribution resulting from neutral beam injectioin in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `ni`: list of ion densities [m^-3]

* `Ti`: list of ion temperatures [eV]

* `mi`: list of ion masses [amu]

* `Zi`: list of ion charges

* `Ef`: fast ion energy [eV]

* `mf`: mass of fast ion [amu]

* `Zf`: fast ion charge
"""
function _electron_ion_drag_difference(
    ne::Real,
    Te::Real,
    ni::AbstractVector{<:Real},
    Ti::AbstractVector{<:Real},
    mi::AbstractVector{<:Real},
    Zi::AbstractVector{Int},
    Ef::Real,
    mf::Real,
    Zf::Int
)
    m_e = mks.m_e
    m_f = mf * mks.m_u

    v_f = sqrt(2 * Ef * mks.e / m_f)
    v_e = sqrt(2 * Te * mks.e / m_e)

    lnΛ_fe = lnΛ_ei(ne, Te)
    Γ_fe = _drag_coefficient(ne, -1, mf, Zf, lnΛ_fe)
    electron_drag = ((8 * Γ_fe * m_f) / (3 * sqrt(pi) * m_e * v_e^3)) * v_f^3

    lnΛ_fis = k -> lnΛ_fi(ne, Te, ni[k], Ti[k], mi[k], Zi[k], v_f / mks.c, mf, Zf; verbose=false)
    Γ_fi = k -> _drag_coefficient(ni[k], Zi[k], mf, Zf, lnΛ_fis(k))
    ion_drag = 2 * m_f * sum(Γ_fi(k) / mi[k] for k in eachindex(mi)) / mks.m_u

    return electron_drag - ion_drag
end

"""
    critical_energy(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Int; approximate::Bool=false)

Returns `Ec` the critical energy by finding the root of the difference between the electron and ion drag

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `ni`: list of ion densities [m^-3]

* `Ti`: list of ion temperatures [eV]

* `mi`: list of ion masses [amu]

* `Zi`: list of ion charges

* `mf`: mass of fast ion [amu]

* `Zf`: fast ion charge

* `approximate`: calculate critical energy assuming `lnΛ_fe == lnΛ_fi`. For DIII-D this results in a correction factor of (lnΛ_fi/lnΛ_fe)^(2/3) ≈ 1.2.
"""
function critical_energy(
    ne::Real,
    Te::Real,
    ni::AbstractVector{<:Real},
    Ti::AbstractVector{<:Real},
    mi::AbstractVector{<:Real},
    Zi::AbstractVector{Int},
    mf::Real,
    Zf::Int;
    approximate::Bool=false
)
    avg_cmr = sum(ni .* (Zi .^ 2) ./ mi) / ne
    Ec = 14.8 * mf * Te * avg_cmr^(2.0 / 3.0)
    if !approximate
        index = ni .> 0.0 .&& Ti .> 0.0
        Ec = Roots.find_zero(Ec0 -> _electron_ion_drag_difference(ne, Te, ni[index], Ti[index], mi[index], Zi[index], Ec0, mf, Zf), (0.5 * Ec, 2 * Ec))
    end
    return Ec
end

@compat public critical_energy
push!(document[Symbol("Physics fast")], :critical_energy)

"""
    thermalization_time(v_f, v_c, tau_s)

Calculate thermalization time in seconds

* `v_f`: fast ion velocity

* `v_c`: critical velocity

* `tau_s`: slowing down time
"""
function thermalization_time(v_f::Real, v_c::Real, tau_s::Real)
    vf3 = v_f^3
    vc3 = v_c^3
    return tau_s * log((vf3 + vc3) / vc3) / 3
end

"""
    thermalization_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real, mf::Real, Zf::Int)

Calculate thermalization time in seconds of a fast ion with energy Ef and `Ti*me/mi < 10Zi^2 eV < Te`

* `ne`: electron density [m^-3]

* `Te`: electron temperature [eV]

* `ni`: list of ion densities [m^-3]

* `Ti`: list of ion temperatures [eV]

* `mi`: list of ion masses [amu]

* `Zi`: list of ion charges

* `Ef`: fast ion energy [eV]

* `mf`: mass of fast ion [amu]

* `Zf`: fast ion charge
"""
function thermalization_time(
    ne::Real,
    Te::Real,
    ni::AbstractVector{<:Real},
    Ti::AbstractVector{<:Real},
    mi::AbstractVector{<:Real},
    Zi::AbstractVector{Int},
    Ef::Real,
    mf::Real,
    Zf::Int
)
    m_f = mf * mks.m_u

    tau_s = slowing_down_time(ne, Te, mf, Zf)
    Ec = critical_energy(ne, Te, ni, Ti, mi, Zi, mf, Zf)

    v_f = sqrt(2 * mks.e * Ef / m_f)
    v_c = sqrt(2 * mks.e * Ec / m_f)

    return thermalization_time(v_f, v_c, tau_s)
end

@compat public thermalization_time
push!(document[Symbol("Physics fast")], :thermalization_time)

"""
    α_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d)

Returns the thermalization time in seconds of α particles evaluated on axis
"""
function α_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d)
    α = ion_properties(:α)
    return thermalization_time(
        cp1d.electrons.density_thermal[1],
        cp1d.electrons.temperature[1],
        [ion.density_thermal[1] for ion in cp1d.ion],
        [ion.temperature[1] for ion in cp1d.ion],
        [ion.element[1].a for ion in cp1d.ion],
        [Int(ion.element[1].z_n) for ion in cp1d.ion],
        mks.E_α,
        α.a,
        Int(α.z_n))
end

@compat public α_thermalization_time
push!(document[Symbol("Physics fast")], :α_thermalization_time)

"""
    beam_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, ion, ion_energy::Real)

Returns the beam thermalization time in seconds evaluated on axis
"""
function beam_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, ion::IMAS.nbi__unit___species, ion_energy::Real)
    return thermalization_time(
        cp1d.electrons.density_thermal[1],
        cp1d.electrons.temperature[1],
        [ion.density_thermal[1] for ion in cp1d.ion],
        [ion.temperature[1] for ion in cp1d.ion],
        [ion.element[1].a for ion in cp1d.ion],
        [Int(ion.element[1].z_n) for ion in cp1d.ion],
        ion_energy,
        ion.a,
        Int(ion.z_n))
end

@compat public beam_thermalization_time
push!(document[Symbol("Physics fast")], :beam_thermalization_time)

"""
    fast_particles!(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; verbose::Bool=false)

Fills the core_profiles fast ion densities and pressures that result from fast ion sources (eg. fusion and nbi)

This calculation is done based on the `slowing_down_time` and `thermalization_time` of the fast ion species.
"""
function fast_particles!(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; verbose::Bool=false)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    Npsi = length(ne)
    Nions = length(cp1d.ion)

    # prepare inputs for slowing_down_time() and thermalization_time() functions
    ni = zeros(Nions, Npsi)
    Ti = zeros(Nions, Npsi)
    Zi = zeros(Int, Nions)
    mi = zeros(Nions)
    for (idx, ion) in enumerate(cp1d.ion)
        if !ismissing(ion, :temperature) # ion temperature may be missing for purely fast-ions species
            ni[idx, :] = ion.density_thermal
            Ti[idx, :] = ion.temperature
            Zi[idx] = Int(ion.element[1].z_n)
            mi[idx] = ion.element[1].a
        end
    end

    # empty cp1d pressures (expressions)
    empty!(cp1d, :pressure)
    empty!(cp1d, :pressure_parallel)
    empty!(cp1d, :pressure_perpendicular)
    empty!(cp1d, :pressure_ion_total)
    # empty all cp1d fast-ion related quantities (expressions)
    for ion in cp1d.ion
        empty!(ion, :pressure)
        empty!(ion, :density)
    end
    # zero out all cp1d fast-ion related quantities
    for ion in cp1d.ion
        ion.pressure_fast_parallel = zeros(Npsi)
        ion.pressure_fast_perpendicular = zeros(Npsi)
        ion.density_fast = zeros(Npsi)
    end
    # zero out cp1d.electrons.density_fast to "close" the expression
    if ismissing(cp1d.electrons, :density_fast)
        cp1d.electrons.density_fast = zeros(Npsi)
    end

    # go through sources and look for ones that have ion particles source at given energy
    taus = zeros(Npsi)
    taut = zeros(Npsi)
    i4 = zeros(Npsi)
    for source in cs.source
        source1d = source.profiles_1d[]
        for sion in source1d.ion
            if !ismissing(sion, :particles) && sum(sion.particles) > 0.0 && !ismissing(sion, :fast_particles_energy) && sion.fast_particles_energy > 0.0
                particle_mass = sion.element[1].a
                particle_charge = Int(sion.element[1].z_n)
                particle_energy = sion.fast_particles_energy

                # find the corresponding thermal ion in core_profiles (we use z_n and a since that's more reliable than labels)
                # NOTE: contribution of non-thermal components in core_profiles comes in through pressure_fast and density_fast
                cindex = findfirst(
                    cion ->
                        (Int(round(cion.element[1].z_n)) == Int(round(particle_charge)) && Int(round(cion.element[1].a)) == Int(round(particle_mass))) ||
                            (Int(round(particle_charge)) == 1 && Int(round(particle_mass)) == 2 && cion.label == "DT") ||
                            (Int(round(particle_charge)) == 1 && Int(round(particle_mass)) == 3 && cion.label == "DT"), cp1d.ion)

                if cindex === nothing
                    if verbose
                        println("$(sion.label) --> ?")
                    end
                else
                    cion = cp1d.ion[cindex]
                    if verbose
                        println("$(sion.label) --> $(cion.label) @ $(sion.fast_particles_energy)")
                    end

                    taus .*= 0.0
                    taut .*= 0.0
                    for i in 1:Npsi
                        taus[i] = slowing_down_time(ne[i], Te[i], particle_mass, particle_charge)
                        taut[i] = @views thermalization_time(ne[i], Te[i], ni[:, i], Ti[:, i], mi, Zi, particle_energy, particle_mass, particle_charge)
                        _, i4[i] = @views estrada_I_integrals(ne[i], Te[i], ni[:, i], Ti[:, i], mi, Zi, particle_energy, particle_mass, particle_charge)
                    end

                    sion_particles = interp1d(source1d.grid.rho_tor_norm, sion.particles).(cp1d.grid.rho_tor_norm)
                    pressa = i4 .* taus .* 2.0 ./ 3.0 .* (sion_particles .* particle_energy .* mks.e)
                    cion.pressure_fast_parallel += pressa ./ 3.0
                    cion.pressure_fast_perpendicular += pressa ./ 3.0
                    cion.density_fast += sion_particles .* taut
                end
            end
        end
    end
end

@compat public fast_particles!
push!(document[Symbol("Physics fast")], :fast_particles!)

"""
    sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)

Compute a low-accuracy but fast approximation of the ion to electron heating fraction for fast particles (like alpha particles and beam particles)
"""
function sivukhin_fraction(cp1d::IMAS.core_profiles__profiles_1d, particle_energy::Real, particle_mass::Real)
    Te = cp1d.electrons.temperature
    ne = cp1d.electrons.density_thermal
    rho = cp1d.grid.rho_tor_norm

    tp = typeof(promote(Te[1], ne[1], rho[1])[1])
    c_a = zeros(tp, length(rho))
    for ion in cp1d.ion
        if !ismissing(ion, :temperature) # ion temperature may be missing for purely fast-ions species
            ni = ion.density_thermal
            @assert all(ni .>= 0.0) "$ni"
            Zi = avgZ(ion.element[1].z_n, ion.temperature)
            mi = ion.element[1].a
            c_a .+= (ni ./ ne) .* Zi .^ 2 ./ (mi ./ particle_mass)
        end
    end

    W_crit = Te .* (4.0 .* sqrt.(mks.m_e / (mks.m_p * particle_mass)) ./ (3.0 * sqrt(pi) .* c_a)) .^ (-2.0 / 3.0)
    ion_to_electron_fraction = similar(W_crit)
    x = particle_energy ./ W_crit
    for (idx, x_i) in enumerate(x)
        if x_i > 4.0
            # Large-x asymptotic formula
            f = (2 * pi / 3) / sin(2 * pi / 3) - 2.0 / sqrt(x_i) + 0.5 / x_i^2
            ion_to_electron_fraction[idx] = f / x_i
        elseif x_i < 0.1
            # Small-x asymptotic series
            ion_to_electron_fraction[idx] = 1.0 - 0.4 * x_i^1.5
        else
            y = x_i .* range(0, 1, 12)
            f = trapz(y, 1.0 ./ (1.0 .+ y .^ 1.5))
            ion_to_electron_fraction[idx] = f / x_i
        end
    end

    return ion_to_electron_fraction
end

@compat public sivukhin_fraction
push!(document[Symbol("Physics fast")], :sivukhin_fraction)

"""
    smooth_beam_power(power::Vector{Float64}, time::AbstractVector{Float64}, taus::Float64) where {T<:Real}

Smooths out the beam power history based on a given thermalization constant `taus`

Such smoothing mimics the delayed contribution from the instantaneous source,
as well as the gradual decay of the previous fast ion population over time.
"""
function smooth_beam_power(time::AbstractVector{Float64}, power::AbstractVector{T}, taus::Float64) where {T<:Real}
    @assert taus >= 0.0
    if taus == 0.0
        return power
    end

    n = length(power)
    smoothed_power = similar(power)

    if !isempty(power)
        smoothed_power[1] = power[1]
        for i in 2:n
            dt = time[i] - time[i - 1]

            # Calculate the decay factor for the current time step
            decay_factor = exp(-dt / taus)

            # Calculate the source contribution
            source_term = power[i] * (1.0 - decay_factor)

            # Update the smoothed density
            smoothed_power[i] = smoothed_power[i - 1] * decay_factor + source_term
        end
    end

    return smoothed_power
end

"""
    smooth_beam_power(time::AbstractVector{Float64}, power::AbstractVector{T}, time0::Float64, taus::Float64) where {T<:Real}

return smoothed beam power at time0
"""
function smooth_beam_power(time::AbstractVector{Float64}, power::AbstractVector{T}, time0::Float64, taus::Float64) where {T<:Real}
    n = length(power)

    smoothed_power = power[1]
    if time0 > time[1]
        for i in 2:n
            dt = time[i] - time[i - 1]

            # Calculate the decay factor for the current time step
            decay_factor = exp(-dt / taus)

            # Calculate the source contribution
            source_term = power[i] * (1.0 - decay_factor)

            # Update the smoothed density
            smoothed_power = smoothed_power * decay_factor + source_term

            if time[i] >= time0
                break
            end
        end
    end

    return smoothed_power
end

@compat public smooth_beam_power
push!(document[Symbol("Physics fast")], :smooth_beam_power)
