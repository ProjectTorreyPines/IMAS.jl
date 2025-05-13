document[Symbol("Physics fast")] = Symbol[]

"""
    estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Real)

Returns solution to `i2` and `i4` integrals from  [Estrada et al.,  Phys of Plasm. 13, 112303 (2006)] Eq. 9 & 10

  - `ne`: electron density [m^-3]
  - `Te`: electron temperature [eV]
  - `ni`: list of ion densities [m^-3]
  - `Ti`: list of ion temperatures [eV]
  - `mi`: list of ion masses [amu]
  - `Zi`: list of ion charges
  - `Ef`: fast ion energy [eV]
  - `mf`: mass of fast ion [amu]
  - `Zf`: fast ion charge
"""
function estrada_I_integrals(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, Ef::Real, mf::Real, Zf::Real)
    Ec = critical_energy(cp1d, irho, mf, Zf)
    i2, i4 = estrada_I_integrals(Ec, Ef)
    return (i2=i2, i4=i4)
end

function estrada_I_integrals(cp1d::IMAS.core_profiles__profiles_1d, Ef::Real, mf::Real, Zf::Real)
    Ec = critical_energy(cp1d, mf, Zf)
    return estrada_I_integrals(Ec, Ef)
end

"""
    estrada_I_integrals(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, Ef::Real,  mf::Real, Zf::Real)

Returns solution to i2 and i4 integrals from  [Estrada et al.,  Phys of Plasm. 13, 112303 (2006)] Eq. 9 & 10

  - `Ef`: fast ion energy [eV]
  - `Ec`: critical energy [eV]
"""
function estrada_I_integrals(Ec::Union{Real,Vector{<:Real}}, Ef::Real)
    a = sqrt.(Ec ./ Ef)
    i2 = (1 / 3.0) .* log.((1 .+ a .^ 3) ./ (a .^ 3))
    i4 = 0.5 .- a .^ 2 .* ((1 / 6.0) .* log.((1 .- a .+ a .^ 2) ./ (1 .+ a) .^ 2) .+
                           1 ./ sqrt(3.0) .* (atan.((2 .- a) ./ (a .* sqrt(3.0))) .+ pi / 6))
    return (i2=i2, i4=i4)
end

@compat public estrada_I_integrals
push!(document[Symbol("Physics fast")], :estrada_I_integrals)

"""
    slowing_down_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Real)

Returns the slowing down time `τ_s` in seconds [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

  - `ne`: electron density [m^-3]
  - `Te`: electron temperature [eV]
  - `ni`: list of ion densities [m^-3]
  - `Ti`: list of ion temperatures [eV]
  - `mi`: list of ion masses [amu]
  - `Zi`: list of ion charges
  - `mf`: mass of fast ion [amu]
  - `Zf`: fast ion charge
"""
function slowing_down_time(ne::Real, Te::Real, ni::AbstractVector{<:Real}, Ti::AbstractVector{<:Real}, mi::AbstractVector{<:Real}, Zi::AbstractVector{Int}, mf::Real, Zf::Real)
    lnΛ = lnΛ_ei(ne, Te, ni, Ti, mi, Zi)
    ne_cm3 = 1e-6 * ne
    τ_s = 6.27e8 * mf * (Te^1.5) ./ (ne_cm3 * lnΛ * Zf^2)
    return τ_s
end

"""
    slowing_down_time(ne::Real, Te::Real, mf::Real, Zf::Real)

Calculates the slowing down time `τ_s` in seconds for `Ti*me/mi < 10Zi^2 eV < Te` [Stix, Plasma Phys. 14 (1972) 367] Eq. 16

  - `ne`: electron density [m^-3]
  - `Te`: electron temperature [eV]
  - `mf`: mass of fast ion [amu]
  - `Zf`: fast ion charge
"""
function slowing_down_time(ne::Real, Te::Real, mf::Real, Zf::Real)
    lnΛ = lnΛ_ei(ne, Te)
    ne_cm3 = 1e-6 * ne
    τ_s = 6.27e8 * mf * (Te^1.5) / (ne_cm3 * lnΛ * Zf^2)
    return τ_s
end

"""
    α_slowing_down_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int=1)

Returns the slowing down time in seconds of α particles evaluated at irho (by default on axis)
"""
function α_slowing_down_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int=1)
    α = ion_properties(:α)
    return slowing_down_time(cp1d.electrons.density_thermal[irho], cp1d.electrons.temperature[irho], α.a, α.z_n)
end

@compat public α_slowing_down_time
push!(document[Symbol("Physics fast")], :α_slowing_down_time)

"""
    drag_coefficient(n::Real, Z::Real, mf::Real, Zf::Real, lnΛ::Real)

Returns drag coefficient `Γ` for a fast-ions interacting with a thermal species as defined Eq. 8 in [Gaffey, J. D. (1976). Energetic ion distribution resulting from neutral beam injection in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

  - `n`: density of the thermal species [m^-3]
  - `Z`: charge of the thermal species
  - `mf`: fast-ion mass [amu]
  - `Zf`: fast-ion charge
  - `lnΛ`: Couloumb logarithm
"""
function drag_coefficient(n::Real, Z::Real, mf::Real, Zf::Real, lnΛ::Real)
    n *= 1e-6 #cm^-3
    mf *= mks.m_u # kg
    Γ = (2 * pi * n * (mks.e)^4 * Z^2 * Zf^2 * lnΛ) / (mf^2)
    return Γ
end

"""
    electron_ion_drag_difference(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, Ef::Real, mf::Real, Zf::Real)

Returns `ΔD` the difference of the electron and ion drag terms in the collision operator defined in Eq. 19 in [Gaffey, J.D (1976). Energetic ion distribution resulting from neutral beam injectioin in tokamaks. Journal of Plasma Physics, 16(02), 149. doi:10.1017/s0022377800020134]

  - `Ef`: fast ion energy [eV]
  - `mf`: mass of fast ion [amu]
  - `Zf`: fast ion charge
"""
function electron_ion_drag_difference(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, Ef::Real, mf::Real, Zf::Real)
    ne = cp1d.electrons.density_thermal[irho]
    Te = cp1d.electrons.temperature[irho]

    m_e = mks.m_e
    m_f = mf * mks.m_u

    v_f = sqrt(2 * Ef * mks.e / m_f)
    v_e = sqrt(2 * Te * mks.e / m_e)

    lnΛ_fe = lnΛ_ei(ne, Te)
    Γ_fe = drag_coefficient(ne, -1, mf, Zf, lnΛ_fe)
    electron_drag = ((8 * Γ_fe * m_f) / (3 * sqrt(pi) * m_e * v_e^3)) * v_f^3

    ion_drag = 0.0
    for ion in cp1d.ion
        ni = ion.density_thermal[irho]
        Ti = ion.temperature[irho]
        Zi = ion.element[1].z_n
        mi = ion.element[1].a
        if ni > 0.0 && Ti > 0.0
            lnΛ_fis = lnΛ_fi(ne, Te, ni, Ti, mi, Zi, v_f / mks.c, mf, Zf)
            Γ_fi = drag_coefficient(ni, Zi, mf, Zf, lnΛ_fis)
            ion_drag += Γ_fi / mi
        end
    end
    ion_drag = 2 * m_f * ion_drag / mks.m_u

    return electron_drag - ion_drag
end

"""
    critical_energy(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, mf::Real, Zf::Real; approximate::Bool=true)

Returns `Ec` the critical energy by finding the root of the difference between the electron and ion drag

  - `mf`: mass of fast ion [amu]
  - `Zf`: fast ion charge
  - `approximate`: calculate critical energy assuming `lnΛ_fe == lnΛ_fi`. For DIII-D this results in a correction factor of (lnΛ_fi/lnΛ_fe)^(2/3) ≈ 1.2.
"""
function critical_energy(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, mf::Real, Zf::Real; approximate::Bool=true)
    ne = cp1d.electrons.density_thermal[irho]
    Te = cp1d.electrons.temperature[irho]

    avg_cmr = sum(ion.density_thermal[irho] * (ion.element[1].z_n^2) / ion.element[1].a for ion in cp1d.ion) / ne
    Ec = 14.8 * mf * Te * avg_cmr^(2.0 / 3.0)
    if !approximate
        Ec = Roots.find_zero(Ec0 -> electron_ion_drag_difference(cp1d, irho, Ec0, mf, Zf), (0.5 * Ec, 2 * Ec))
    end
    return Ec
end

function critical_energy(cp1d::IMAS.core_profiles__profiles_1d, mf::Real, Zf::Real; approximate::Bool=true)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    avg_cmr = sum(ion.density_thermal .* (ion.element[1].z_n^2) / ion.element[1].a for ion in cp1d.ion) ./ ne
    Ec = 14.8 .* mf .* Te .* avg_cmr .^ 1.5
    if !approximate
        Ec1 = critical_energy(cp1d, 1, mf, Zf; approximate)
        Ec = Ec .* (Ec1 ./ Ec[1])
    end
    return Ec
end

@compat public critical_energy
push!(document[Symbol("Physics fast")], :critical_energy)

"""
    thermalization_time(v_f, v_c, tau_s)

Calculate thermalization time in seconds

  - `v_f`: fast ion velocity
  - `v_c`: critical velocity
  - `tau_s`: slowing down time
"""
function thermalization_time(v_f::Real, v_c::Real, tau_s::Real)
    vf3 = v_f^3
    vc3 = v_c^3
    return tau_s * log((vf3 + vc3) / vc3) / 3
end

"""
    thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, Ef::Real, mf::Real, Zf::Real)

Calculate thermalization time in seconds of a fast ion with energy Ef and `Ti*me/mi < 10Zi^2 eV < Te`

  - `Ef`: fast ion energy [eV]
  - `mf`: mass of fast ion [amu]
  - `Zf`: fast ion charge
"""
function thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, Ef::Real, mf::Real, Zf::Real)
    m_f = mf * mks.m_u

    tau_s = slowing_down_time(cp1d.electrons.density_thermal[irho], cp1d.electrons.temperature[irho], mf, Zf)
    Ec = critical_energy(cp1d, irho, mf, Zf)

    v_f = sqrt(2 * mks.e * Ef / m_f)
    v_c = sqrt(2 * mks.e * Ec / m_f)

    return thermalization_time(v_f, v_c, tau_s)
end

function thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, Ef::Real, mf::Real, Zf::Real)
    m_f = mf * mks.m_u

    tau_s = slowing_down_time.(cp1d.electrons.density_thermal, cp1d.electrons.temperature, mf, Zf)
    Ec = critical_energy(cp1d, mf, Zf)

    v_f = sqrt(2 * mks.e * Ef / m_f)
    v_c = sqrt.(2 .* mks.e .* Ec ./ m_f)

    return thermalization_time.(v_f, v_c, tau_s)
end

@compat public thermalization_time
push!(document[Symbol("Physics fast")], :thermalization_time)

"""
    α_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int)

Returns the thermalization time in seconds of α particles evaluated on axis
"""
function α_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int)
    α = ion_properties(:α)
    return thermalization_time(cp1d, irho, mks.E_α, α.a, Int(α.z_n))
end

function α_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d)
    α = ion_properties(:α)
    return thermalization_time(cp1d, mks.E_α, α.a, Int(α.z_n))
end

@compat public α_thermalization_time
push!(document[Symbol("Physics fast")], :α_thermalization_time)

"""
    fast_ion_thermalization_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, ion::Union{IMAS.nbi__unit___species,IMAS.core_profiles__profiles_1d___ion___element}, ion_energy::Real)

Returns the fast ion thermalization time in seconds evaluated on axis
"""
function fast_ion_thermalization_time(
    cp1d::IMAS.core_profiles__profiles_1d,
    irho::Int,
    ion::Union{IMAS.nbi__unit___species,IMAS.core_profiles__profiles_1d___ion___element},
    ion_energy::Real
)
    return thermalization_time(cp1d, irho, ion_energy, ion.a, ion.z_n)
end

function fast_ion_thermalization_time(
    cp1d::IMAS.core_profiles__profiles_1d,
    ion::Union{IMAS.nbi__unit___species,IMAS.core_profiles__profiles_1d___ion___element},
    ion_energy::Real
)
    return thermalization_time(cp1d, ion_energy, ion.a, ion.z_n)
end

@compat public fast_ion_thermalization_time
push!(document[Symbol("Physics fast")], :fast_ion_thermalization_time)

"""
    fast_particles_profiles!(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; verbose::Bool=false)

Fills the core_profiles fast ion densities and pressures that result from fast ion sources (eg. fusion and nbi)

This calculation is done based on the `slowing_down_time` and `thermalization_time` of the fast ion species.
"""
function fast_particles_profiles!(cs::IMAS.core_sources, cp1d::IMAS.core_profiles__profiles_1d; verbose::Bool=false)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    Npsi = length(ne)

    # empty cp1d pressures (expressions)
    IMAS.unfreeze!(cp1d, :pressure)
    IMAS.unfreeze!(cp1d, :pressure_parallel)
    IMAS.unfreeze!(cp1d, :pressure_perpendicular)
    IMAS.unfreeze!(cp1d, :pressure_ion_total)
    # empty all cp1d fast-ion related quantities (expressions)
    for ion in cp1d.ion
        IMAS.unfreeze!(ion, :pressure)
        IMAS.unfreeze!(ion, :density)
    end

    # zero out all cp1d fast-ion related quantities
    for ion in cp1d.ion
        ion.pressure_fast_parallel = zeros(Npsi)
        ion.pressure_fast_perpendicular = zeros(Npsi)
        ion.density_fast = zeros(Npsi)
    end

    # go through sources and look for ones that have ion particles source at given energy
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

                    taus = slowing_down_time.(ne, Te, particle_mass, particle_charge)
                    taut = thermalization_time(cp1d, particle_energy, particle_mass, particle_charge)
                    _, i4 = estrada_I_integrals(cp1d, particle_energy, particle_mass, particle_charge)

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

function fast_particles_profiles!(dd::IMAS.dd; verbose::Bool=false)
    return fast_particles_profiles!(dd.core_sources, dd.core_profiles.profiles_1d[]; verbose)
end

@compat public fast_particles_profiles!
push!(document[Symbol("Physics fast")], :fast_particles_profiles!)

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
            @assert all(ni .>= 0.0) "Ion `$(ion.label)` has negative densities\n$ni"
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
    smooth_beam_power(time::AbstractVector{Float64}, power::AbstractVector{T}, taus::Float64) where {T<:Real}

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
            dt = time[i] - time[i-1]

            # Calculate the decay factor for the current time step
            decay_factor = exp(-dt / taus)

            # Calculate the source contribution
            source_term = power[i] * (1.0 - decay_factor)

            # Update the smoothed density
            smoothed_power[i] = smoothed_power[i-1] * decay_factor + source_term
        end
    end

    return smoothed_power
end

"""
    smooth_beam_power(time::AbstractVector{Float64}, power::AbstractVector{T}, time0::Float64, taus::Float64) where {T<:Real}

Return smoothed beam power at time0
"""
function smooth_beam_power(time::AbstractVector{Float64}, power::AbstractVector{T}, time0::Float64, taus::Float64) where {T<:Real}
    n = length(power)

    smoothed_power = power[1]
    if time0 > time[1]
        for i in 2:n
            dt = time[i] - time[i-1]

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

"""
    banana_width(T::Real, Bt::Real, Z::Real, m::Real, epsilon::Real, q::Real)

Estimates the banana orbit width [m]

  - `T`: Temperature [eV]
  - `Bt`: Magnetic field [T]
  - `Z`: Charge
  - `m`: Mass [AMU]
  - `epsilon`: Inverse aspect ratio
  - `q`: safety factor
"""
function banana_width(T::Real, Bt::Real, Z::Real, m::Real, epsilon::Real, q::Real)
    r_gyro = gyroradius(T, Bt, Z, m)
    return 2.0 * max(0.0, epsilon)^(-0.5) * abs(q) * r_gyro
end

"""
    gyroradius(T::Real, Bt::Real, Z::Real, m::Real) 

Calculates plasma gyroradius [m]

  - `T`: Ion temperature [eV]
  - `Bt`: Magnetic field [T]
  - `Z`: charge
  - `m`: Mass [AMU]
"""
function gyroradius(T::Real, Bt::Real, Z::Real, m::Real)
    M = m * mks.m_p
    vt = sqrt(T / M * mks.e)
    return M * vt / abs(mks.e * Z * Bt)
end


"""
    imfp_charge_exchange(atw::Real, e::Real, zni::Real)

Calculates the inverse mean free path due to charge exchange based on the fitted results of Freeman and Jones (1974).

Routine pulled from freya_fsgxn.f90

  - `atw`: Beam mass [AMU]
  - `e`: Beam energy divided by beam mass  [eV/AMU]
  - `zni`: Plasma ion density [m^-3]
"""
function imfp_charge_exchange(atw::Real, e::Real, zni::Real)
    if atw > 3.01
        return 0.0
    else
        aloge = log10(e)
        sigcx = (0.6937e-14 * (1.0 - 0.155 * aloge)^2 / (1.0 + 0.1112e-14 * e^3.3))
    end
    return 1e2 * sigcx * (zni * 1e-6)
end

"""
    imfp_ion_collisions(atw::Real, eova::Real, zni::Real, zzi::Real)

Calculates inverse mean free path due to proton and impurity impact ionization.

Routine pulled from freya_fsgxn.f90

  - `atw`: beam mass [AMU]
  - `eova`: Beam energy divided by beam mass  [eV/AMU]
  - `zni`: Plasma ion density [m^-3]
  - `zzi`: Plasma ion charge [AMU]
"""
function imfp_ion_collisions(atw::Real, eova::Real, zni::Real, zzi::Real)
    cfionp = (-4.203309e+01, 3.557321, -1.045134, 0.3139238, -0.07454475, 0.008459113, -3.495444e-04)
    if atw <= 3.01
        aloge = log10(eova) * 2.302585093 - 6.907755279
        if aloge <= -2.30258
            sigi = 0.0
        else
            expo = (((((cfionp[7] * aloge + cfionp[6]) * aloge + cfionp[5]) * aloge
                      +
                      cfionp[4]) * aloge + cfionp[3]) * aloge + cfionp[2]) * aloge + cfionp[1]
            sigi = exp(expo)
        end
        return 1e2 * sigi * (zni * 1e-6)
    else
        ekev = 1.0e-3 * eova
        return 1e2 * 1.0e-17 * (zni * 1e-6) * 46.0 * zzi * (32.0 * zzi / ekev) * (1.0 - exp(-ekev / (32.0 * zzi)))
    end
end

"""
    imfp_electron_collisions(vb::Real, te::Real, zne::Real)

Evaluates local inverse mean free path for electron impact ionization.

Routine pulled from freya_fsgxn.f90

  - `vb`: Velocity of neutral beam [m/s]
  - `te`: Electron temperature [eV]
  - `zne`: Electron density [m^-3]
"""
function imfp_electron_collisions(vb::Real, te::Real, zne::Real)
    cfione = (-3.173850e+01, 1.143818e+01, -3.833998, 0.7046692, -0.07431486, 0.004153749, -9.486967e-05)

    alogt = te > 1.0 ? log(te) : 0.0
    alogt = te > 1.0e+05 ? 11.51 : alogt

    expo = (((((cfione[7] * alogt + cfione[6]) * alogt + cfione[5]) * alogt
              +
              cfione[4]) * alogt + cfione[3]) * alogt + cfione[2]) * alogt + cfione[1]

    return 1e2 * exp(expo) * (zne * 1e-6) / (vb * 1E2)
end

"""
    bkefun(y::Real, vcvo::Real, tstcx::Real, emzrat::Real)
"""
function bkefun(y::Real, vcvo::Real, tstcx::Real, emzrat::Real)
    if y <= 0
        return 0.0
    end
    y3 = y^3
    v3 = vcvo^3
    y3v3 = y3 + v3
    inv_y3v3 = 1 / y3v3
    arg = (1 + v3) * inv_y3v3
    alogarg = log(arg)
    pcxlog = -tstcx * alogarg / 3
    log_y = log(y)
    alog3y = 3 * log_y
    alogy3v3 = log(y3v3)
    blog = (alog3y + alogarg) * emzrat / 3
    bkeflog = alog3y + pcxlog + blog - alogy3v3
    return bkeflog < -30 ? 0.0 : exp(bkeflog)
end

"""
    ion_momentum_fraction(vpar::Real, tpar::Real, emzpar::Real; N=100)
"""
function ion_momentum_fraction(vpar::Real, tpar::Real, emzpar::Real; N=100)
    vcvo = vpar
    tstcx = tpar
    emzrat = emzpar
    out = 0.0
    @inbounds @simd for y1 in range(0.0, 1.0, N)
        out += bkefun(y1, vcvo, tstcx, emzrat)
    end
    return out / N
end

"""
    ion_momentum_slowingdown_time(cp1d::IMAS.core_profiles__profiles_1d, Ef::Real, mf::Real, Zf::Real)
"""
function ion_momentum_slowingdown_time(cp1d::IMAS.core_profiles__profiles_1d, Ef::Real, mf::Real, Zf::Real)
    E_c = critical_energy(cp1d, mf, Zf)
    taus = slowing_down_time.(cp1d.electrons.density_thermal, cp1d.electrons.temperature, mf, Zf)
    emzrat = cp1d.ion[1].element[1].a .* cp1d.zeff ./ (mf * Zf)
    bke = ion_momentum_fraction.(sqrt.(E_c / Ef), 0.0, emzrat)
    tau_mom = taus .*  bke
    return tau_mom
end

function ion_momentum_slowingdown_time(cp1d::IMAS.core_profiles__profiles_1d, irho::Int, Ef::Real, mf::Real, Zf::Real)
    E_c = critical_energy(cp1d, irho, mf, Zf)
    taus = slowing_down_time.(cp1d.electrons.density_thermal[irho], cp1d.electrons.temperature[irho], mf, Zf)
    emzrat = cp1d.ion[1].element[1].a .* cp1d.zeff[irho] ./ (mf * Zf)
    bke = ion_momentum_fraction(sqrt(E_c / Ef), 0.0, emzrat)
    tau_mom = taus * bke
    return tau_mom
end