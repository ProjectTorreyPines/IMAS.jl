"""
    beta_tor_thermal_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)

Normalised toroidal beta from thermal pressure only, defined as 100 * beta_tor_thermal * a[m] * B0 [T] / ip [MA]
"""
function beta_tor_thermal_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)
    return beta_tor(eq, cp1d; norm=true, thermal=true)
end

"""
    beta_tor_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)

Normalised toroidal beta from total pressure, defined as 100 * beta_tor * a[m] * B0 [T] / ip [MA]
"""
function beta_tor_norm(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d)
    return beta_tor(eq, cp1d; norm=true, thermal=false)
end

"""
    beta_tor(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d; norm::Bool, thermal::Bool)

Toroidal beta, defined as the volume-averaged total perpendicular pressure divided by (B0^2/(2*mu0)), i.e. beta_toroidal = 2 mu0 int(p dV) / V / B0^2
"""
function beta_tor(eq::IMAS.equilibrium, cp1d::IMAS.core_profiles__profiles_1d; norm::Bool, thermal::Bool)
    eqt = eq.time_slice[cp1d.time]
    eq1d = eqt.profiles_1d
    if thermal
        pressure = cp1d.pressure_thermal
    else
        pressure = cp1d.pressure
    end
    rho = cp1d.grid.rho_tor_norm
    B0 = get_time_array(eq.vacuum_toroidal_field, :b0, eqt.time, :constant)
    Ip = eqt.global_quantities.ip
    volume_cp = interp1d(eq1d.rho_tor_norm, eq1d.volume).(rho)
    pressure_avg = trapz(volume_cp, pressure) / volume_cp[end]
    beta_tor = 2.0 * constants.μ_0 * pressure_avg / B0^2
    if norm
        return beta_tor * eqt.boundary.minor_radius * abs(B0) / abs(Ip / 1e6) * 1.0e2
    else
        return beta_tor
    end
end

function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_z::Int;
    fast::Bool=false)
    return ion_element!(ion, elements[ion_z].symbol; fast)
end

function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_string::AbstractString;
    fast::Bool=false)
    return ion_element!(ion, Symbol(ion_string); fast)
end

"""
    ion_element!(
        ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
        ion_symbol::Symbol)

Fills the `ion.element` structure with the a and z_n information, also updates the `ion.label`
"""
function ion_element!(
    ion::Union{IMAS.core_profiles__profiles_1d___ion,IMAS.core_sources__source___profiles_1d___ion},
    ion_symbol::Symbol;
    fast::Bool=false)

    element = resize!(ion.element, 1)[1]

    # H isotopes
    if ion_symbol ∈ (:H, :H1)
        element.z_n = 1.0
        element.a = 1.00797
        label = "H"

    elseif ion_symbol ∈ (:D, :H2)
        element.z_n = 1.0
        element.a = 2.014
        label = "D"

    elseif ion_symbol ∈ (:T, :H3)
        element.z_n = 1.0
        element.a = 3.016
        label = "T"

    elseif ion_symbol == :DT
        element.z_n = 1.0
        element.a = (2.014 + 3.016) / 2.0
        label = "DT"

    else
        # all other ions
        ion_name, ion_a = match(r"(.*?)(\d*)$", string(ion_symbol))
        element_ion = elements[Symbol(ion_name)]
        element.z_n = float(element_ion.number)
        if isempty(ion_a)
            element.a = element_ion.atomic_mass.val
        else
            z = element_ion.number
            n = parse(Int, ion_a) - z
            element.a = atomic_mass(z, n)
        end
        label = "$(ion_name)$(Int(floor(element.a)))"
    end

    if fast
        ion.label = "$(label)_fast"
    else
        ion.label = label
    end

    return ion.element
end

"""
    binding_energy(Z::Int, N::Int)

Return the estimate binding energy in MeV
"""
function binding_energy(Z::Int, N::Int)
    # Constants (in MeV)
    a_v = 15.8
    a_s = 18.3
    a_c = 0.714
    a_a = 23.2

    # Total number of nucleons
    A = Z + N

    # Calculate the pairing term
    if A % 2 != 0
        δ = 0.0
    elseif Z % 2 == 0
        δ = 12.0 / sqrt(A)
    else
        δ = -12.0 / sqrt(A)
    end

    # Calculate the binding energy (in MeV)
    E = a_v * A - a_s * A^(2 / 3) - a_c * Z * (Z - 1) / A^(1 / 3) - a_a * (N - Z)^2 / A + δ

    return E
end

"""
    atomic_mass(Z::Int, N::Int)

Returns the estimated nucleus mass including the estimated effect of binding energy
"""
function atomic_mass(Z::Int, N::Int)
    # Mass of proton and neutron (in amu)
    mass_proton = IMAS.constants.m_p / IMAS.constants.m_u
    mass_neutron = IMAS.constants.m_n / IMAS.constants.m_u

    # Conversion factor from MeV/c^2 to atomic mass units (amu)
    mev_to_amu = IMAS.constants.m_u * IMAS.constants.c^2 / (IMAS.constants.e * 1E6)

    # Estimate binding energy
    E = binding_energy(Z, N)

    # Calculate mass of nucleus (in amu)
    mass = Z * mass_proton + N * mass_neutron - E / mev_to_amu

    return mass
end

function energy_thermal(cp1d::IMAS.core_profiles__profiles_1d)
    return 3.0 / 2.0 * trapz(cp1d.grid.volume, cp1d.pressure_thermal)
end

"""
    energy_thermal_ped(cp1d::IMAS.core_profiles__profiles_1d, su::IMAS.summary)

Calculates the pedestal contribution to the thermal stored energy by integrating over entire domain but with pedestal pressure in the core
"""
function energy_thermal_ped(cp1d::IMAS.core_profiles__profiles_1d, su::IMAS.summary)
    rho_index = argmin(abs.(cp1d.grid.rho_tor_norm .- su.local.pedestal.position.rho_tor_norm))
    pressure = cp1d.pressure_thermal
    pressure[1:rho_index] .= cp1d.pressure_thermal[rho_index]
    return 3.0 / 2.0 * trapz(cp1d.grid.volume, pressure)
end

function ne_vol_avg(cp1d::IMAS.core_profiles__profiles_1d)
    return trapz(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end]
end

"""
    tau_e_thermal(cp1d::IMAS.core_profiles__profiles_1d, sources::IMAS.core_sources)

Evaluate thermal energy confinement time

NOTE: power losses due to radiation are neglected, as done for tau_e_h98 scaling
"""
function tau_e_thermal(cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources)
    total_source = total_sources(cs, cp1d; fields=[:power_inside, :total_ion_power_inside])
    total_power_inside = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end]
    return energy_thermal(cp1d) / (total_power_inside - radiation_losses(cs))
end


function tau_e_h98(dd::IMAS.dd; time0::Float64=dd.global_time)
    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    cs = dd.core_sources
    return tau_e_h98(eqt, cp1d, cs)
end

"""
    tau_e_h98(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.plot_core_sources)

H98y2 ITER elmy H-mode confinement time scaling
"""
function tau_e_h98(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources)
    total_source = total_sources(cs, cp1d; fields=[:power_inside, :total_ion_power_inside])
    total_power_inside = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end] - radiation_losses(cs)
    isotope_factor =
        trapz(cp1d.grid.volume, sum(ion.density .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) /
        trapz(cp1d.grid.volume, sum(ion.density for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    tau98 = (
        0.0562 *
        abs(eqt.global_quantities.ip / 1e6)^0.93 *
        abs(B0)^0.15 *
        (total_power_inside / 1e6)^-0.69 *
        (ne_vol_avg(cp1d) / 1e19)^0.41 *
        isotope_factor^0.19 *
        R0^1.97 *
        (R0 / eqt.boundary.minor_radius)^-0.58 *
        eqt.boundary.elongation^0.78
    )
    return tau98
end

function tau_e_ds03(dd::IMAS.dd; time0::Float64=dd.global_time)
    eqt = dd.equilibrium.time_slice[time0]
    cp1d = dd.core_profiles.profiles_1d[time0]
    cs = dd.core_sources
    return tau_e_ds03(eqt, cp1d, cs)
end

"""
    tau_e_ds03(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources)

Petty's 2003 confinement time scaling
"""
function tau_e_ds03(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d, cs::IMAS.core_sources)
    total_source = total_sources(cs, cp1d; fields=Symbol[:power_inside, :total_ion_power_inside])
    total_power_inside = total_source.electrons.power_inside[end] + total_source.total_ion_power_inside[end] - radiation_losses(cs)
    isotope_factor =
        trapz(cp1d.grid.volume, sum(ion.density .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) /
        trapz(cp1d.grid.volume, sum(ion.density for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    tauds03 = (
        0.028 *
        abs(eqt.global_quantities.ip / 1e6)^0.83 *
        abs(B0)^0.07 *
        (total_power_inside / 1e6)^-0.55 *
        (ne_vol_avg(cp1d) / 1e19)^0.49 *
        isotope_factor^0.14 *
        R0^2.11 *
        (R0 / eqt.boundary.minor_radius)^-0.30 *
        eqt.boundary.elongation^0.75
    )

    return tauds03
end

"""
    bunit(eqt::IMAS.equilibrium__time_slice)

Calculate bunit from equilibrium
"""
function bunit(eqt::IMAS.equilibrium__time_slice)
    eq1d = eqt.profiles_1d
    rmin = 0.5 .* (eq1d.r_outboard .- eq1d.r_inboard)
    phi = eq1d.phi
    return gradient(2π * rmin, phi) ./ rmin
end

"""
    greenwald_density(eqt::IMAS.equilibrium__time_slice)

Simple greenwald line-averaged density limit
"""
function greenwald_density(eqt::IMAS.equilibrium__time_slice)
    return greenwald_density(eqt.global_quantities.ip, eqt.boundary.minor_radius)
end

function greenwald_density(ip::T, minor_radius::T) where {T<:Real}
    return (ip / 1e6) / (pi * minor_radius^2) * 1e20
end

function greenwald_density(dd::IMAS.dd)
    return greenwald_density(dd.equilibrium.time_slice[])
end

function greenwald_density(ps::IMAS.pulse_schedule; time0=global_time(ps))
    ip = IMAS.get_time_array(ps.flux_control.i_plasma, :reference, time0, :linear)
    minor_radius = IMAS.get_time_array(ps.position_control.minor_radius, :reference, time0, :linear)
    return greenwald_density(ip, minor_radius)
end

"""
    greenwald_fraction(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Greewald fraction
"""
function greenwald_fraction(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    nel = IMAS.geometric_midplane_line_averaged_density(eqt, cp1d)
    ngw = greenwald_density(eqt)
    return nel / ngw
end

function greenwald_fraction(dd::IMAS.dd)
    return greenwald_fraction(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])
end

function geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, ne_profile::AbstractVector{<:Real}, rho_ne::AbstractVector{<:Real})
    a_cp = interp1d(eqt.profiles_1d.rho_tor_norm, (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) / 2.0).(rho_ne)
    return trapz(a_cp, ne_profile) / a_cp[end]
end

"""
    geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)

Calculates the line averaged density from a midplane horizantal line
"""
function geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d)
    return geometric_midplane_line_averaged_density(eqt, cp1d.electrons.density, cp1d.grid.rho_tor_norm)
end

"""
    beta_tor(pressure_average::Real, Bt::Real)

Calculates Beta_tor from pressure and Bt
"""
function beta_tor(pressure_average::Real, Bt::Real)
    return pi * 8.0e-7 * pressure_average / Bt^2
end

"""
    beta_n(beta_tor::Real, minor_radius::Real, Bt::Real, Ip::Real)

Calculates BetaN from beta_tor
"""
function beta_n(beta_tor::Real, minor_radius::Real, Bt::Real, Ip::Real)
    return beta_tor * minor_radius * abs(Bt) / abs(Ip / 1e6) * 1.0e2 # [%]
end

"""
    pressure_avg_from_beta_n(beta_n::Real, minor_radius::Real, Bt::Real, Ip::Real)

Calculates average pressure from BetaN
"""
function pressure_avg_from_beta_n(beta_n::Real, minor_radius::Real, Bt::Real, Ip::Real)
    return beta_n * abs(Bt) * abs(Ip / 1e6) / (minor_radius * pi * 8.0e-7 * 1.0e2)
end

"""
    Hmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real)

Generate H-mode density and temperature profiles evenly spaced in your favorite radial coordinate

:param edge: separatrix height

:param ped: pedestal height

:param core: on-axis profile height

:param ngrid: number of radial grid points

:param expin: inner core exponent for H-mode pedestal profile

:param expout: outer core exponent for H-mode pedestal profile

:param width: width of pedestal
"""
function Hmode_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real)
    @assert core >= 0.0

    xpsi = range(0.0, 1.0, ngrid)

    w_E1 = 0.5 * widthp  # width as defined in eped
    xphalf = 1.0 - w_E1
    pconst = 1.0 - tanh((1.0 - xphalf) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)
    coretanh = 0.5 * a_t * (1.0 - tanh(-xphalf / w_E1) - pconst) + edge

    # edge tanh part
    val = @. 0.5 * a_t * (1.0 - tanh((xpsi - xphalf) / w_E1) - pconst) + edge + core * 0.0

    # core tanh+polynomial part
    xped = xphalf - w_E1
    xtoped = xpsi ./ xped
    for i in 1:ngrid
        if xtoped[i] < 1.0
            @inbounds val[i] += (core - coretanh) * (1.0 - xtoped[i]^expin)^expout
        end
    end

    return val
end

"""
    Hmode_profiles(edge::Real, ped::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real)

Generate H-mode density and temperature profiles evenly spaced in your favorite radial coordinate

NOTE: The core value is allowed to float

:param edge: separatrix height

:param ped: pedestal height

:param ngrid: number of radial grid points

:param expin: inner core exponent for H-mode pedestal profile

:param expout: outer core exponent for H-mode pedestal profile

:param width: width of pedestal
"""
function Hmode_profiles(edge::Real, ped::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real)

    @assert expin >= 0.0
    @assert expout >= 0.0

    xpsi = range(0.0, 1.0, ngrid)

    w_E1 = 0.5 * widthp  # width as defined in eped
    xphalf = 1.0 - w_E1
    pconst = 1.0 - tanh((1.0 - xphalf) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)

    # edge tanh part
    val = @. 0.5 * a_t * (1.0 - tanh((xpsi - xphalf) / w_E1) - pconst) + edge

    # core tanh+polynomial part
    xped = xphalf - w_E1
    xtoped = xpsi ./ xped
    integral = 0.0
    factor = Inf
    for i in ngrid:-1:1
        if xtoped[i] < 1.0
            @inbounds val[i] += integral
            if i > 1
                factor = min(factor, val[i])
                xi = (xtoped[i] + xtoped[i-1]) / 2.0
                dx = (xtoped[i] - xtoped[i-1])
                if expin == 0.0
                    v1 = 0.0
                else
                    v1 = expin * expout * xi^(expin - 1.0) * (1.0 - xi^expin)^(expout - 1.0)
                end
                integral += v1 * dx * factor
            end
        end
    end

    return val
end


"""
    ITB_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real, ITBr::Real, ITBw::Real, ITBh::Real)

Generate H-mode density and temperature profiles with Internal Transport Barrier (ITB).
This makes an H-mode profile and adds a tanh offset at the ITB radial location.
Note that core, ped, and ITBh are all in absolute units, so if core < ped + ITBh there will be negative gradients.

:param edge: separatrix height

:param ped: pedestal height

:param core: on-axis profile height

:param ngrid: number of radial grid points

:param expin: inner core exponent for H-mode pedestal profile

:param expout: outer core exponent for H-mode pedestal profile

:param widthp: width of pedestal

:param ITBr: radial location of ITB center

:param ITBw: width of ITB

:param ITBh: height of ITB

"""
function ITB_profiles(edge::Real, ped::Real, core::Real, ngrid::Int, expin::Real, expout::Real, widthp::Real, ITBr::Real, ITBw::Real, ITBh::Real)

    xpsi = range(0.0, 1.0, ngrid)

    # H mode part

    val = Hmode_profiles(edge, ped, core-ITBh, ngrid, expin, expout, widthp)

    # ITB part

    itb = @. 0.5 * ITBh * (1.0 - tanh((xpsi - ITBr) / ITBw))

    return val + itb
end


"""
    species(cp1d::IMAS.core_profiles__profiles_1d; only_electrons_ions::Symbol=:all, only_thermal_fast::Symbol=:all)

Returns species index and names (followed by "_fast" if density_fast is present), for example:

    (0, :electrons)
    (0, :electrons_fast)
    (1, :DT)
    (2, :Kr83)
    (3, :He4)
    (1, :DT_fast)
    (3, :He4_fast)
"""
function species(cp1d::IMAS.core_profiles__profiles_1d; only_electrons_ions::Symbol=:all, only_thermal_fast::Symbol=:all)
    @assert only_electrons_ions ∈ (:all, :electrons, :ions) "only_electrons_ions can be one of (:all, :electrons, :ions)"
    @assert only_thermal_fast ∈ (:all, :thermal, :fast) "only_thermal_fast can be one of (:all, :thermal, :fast)"
    out = []
    if only_electrons_ions ∈ (:all, :electrons)
        if only_thermal_fast ∈ (:all, :thermal) && sum(cp1d.electrons.density_thermal) > 0.0
            push!(out, (0, :electrons))
        end
        if only_thermal_fast ∈ (:all, :fast) && sum(cp1d.electrons.density_fast) > 0.0
            push!(out, (0, :electrons_fast))
        end
    end
    if only_electrons_ions ∈ (:all, :ions)
        if only_thermal_fast ∈ (:all, :thermal)
            dd_thermal = ((k, Symbol(ion.label)) for (k, ion) in enumerate(cp1d.ion) if sum(ion.density_thermal) > 0.0)
            for item in dd_thermal
                push!(out, item)
            end
        end
        if only_thermal_fast ∈ (:all, :fast)
            dd_fast = ((k, Symbol("$(ion.label)_fast")) for (k, ion) in enumerate(cp1d.ion) if sum(ion.density_fast) > 0.0)
            for item in dd_fast
                push!(out, item)
            end
        end
    end
    return out
end

"""
    is_quasi_neutral(dd::IMAS.dd)
"""
function is_quasi_neutral(dd::IMAS.dd; rtol::Float64=0.001)
    return is_quasi_neutral(dd.core_profiles.profiles_1d[]; rtol)
end

"""
    is_quasi_neutral(cp1d::IMAS.core_profiles__profiles_1d; rtol::Float64=0.001)

Returns true if quasi neutrality is satisfied within a relative tolerance
"""
function is_quasi_neutral(cp1d::IMAS.core_profiles__profiles_1d; rtol::Float64=0.001)
    Nis = sum(sum(ion.density .* ion.z_ion for ion in cp1d.ion))
    Ne = sum(cp1d.electrons.density)
    if (1.0 + rtol) >= (Ne / Nis) >= (1.0 - rtol)
        return true
    else
        return false
    end
end

"""
    enforce_quasi_neutrality!(dd::IMAS.dd, species::Symbol)
"""
function enforce_quasi_neutrality!(dd::IMAS.dd, species::Symbol)
    return enforce_quasi_neutrality!(dd.core_profiles.profiles_1d[], species)
end

"""
    enforce_quasi_neutrality!(cp1d::IMAS.core_profiles__profiles_1d, species::Symbol)

Evaluates the difference in number of charges needed to reach quasineutrality,
and assigns positive difference  to target ion density_thermal species
and negative difference to electrons density_thermal

Also, makes sure `density` is set to the original expression
"""
function enforce_quasi_neutrality!(cp1d::IMAS.core_profiles__profiles_1d, species::Symbol)
    # identify ion species
    species_indx = findfirst(Symbol(ion.label) == species for ion in cp1d.ion)
    @assert species_indx !== nothing
    ion0 = cp1d.ion[species_indx]

    # Make sure expressions are used for total densities
    empty!(cp1d.electrons, :density)
    for ion in cp1d.ion
        empty!(ion, :density)
    end

    # evaluate the difference in number of charges needed to reach quasineutrality
    q_density_difference = cp1d.electrons.density .- sum(ion.density .* ion.z_ion for ion in cp1d.ion if ion != ion0) .- ion0.density_fast .* ion0.z_ion

    # positive difference is assigned to target ion density_thermal
    index = q_density_difference .> 0.0
    ion0.density_thermal = zero(cp1d.electrons.density_thermal)
    ion0.density_thermal[index] .= q_density_difference[index] .* ion0.z_ion

    # negative difference is assigned to electrons density_thermal
    index = q_density_difference .< 0.0
    cp1d.electrons.density_thermal[index] .+= abs.(q_density_difference[index])

    return ion0.density_thermal
end

"""
    broken_lump_ions_as_bulk_and_impurity((ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion}, rho_tor_norm::Vector{<:Real})

Changes core_profiles.ion to 2 species, bulk specie (H, D, T) and combined impurity specie by weigthing masses and densities

NOTE: This version of lump_ions_as_bulk_and_impurity is not correct, but we keep it around for now because that's how some NNs have been trained.
DO NOT USE IF YOU DON'T KNOW WHAT YOU ARE DOING.
"""
function broken_lump_ions_as_bulk_and_impurity(
    ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion{T}},
    rho_tor_norm::Vector{<:T}
)::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion{T}} where {T<:Real}
    if length(ions) < 2
        error("lump_ions_as_bulk_and_impurity requires at least two ion species")
    end

    zs = [ion.element[1].z_n for ion in ions]
    as = [ion.element[1].a for ion in ions]

    bulk_index = findall(zs .== 1)
    impu_index = findall(zs .!= 1)

    ratios = zeros(length(ions[1].density_thermal), length(ions))
    for index in (bulk_index, impu_index)
        ntot = zeros(length(ions[1].density_thermal))
        for ix in index
            tmp = ions[ix].density_thermal * sum(avgZ(zs[ix], ions[ix].temperature)) / length(ions[ix].temperature)
            ratios[:, ix] = tmp
            ntot .+= tmp
        end
        for ix in index
            ratios[:, ix] ./= ntot
        end
    end

    ions2 = IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion{T}}()

    # bulk ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    bulk = ions2[end]
    resize!(bulk.element, 1)
    bulk.label = "bulk"

    # impurity ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    impu = ions2[end]
    resize!(impu.element, 1)
    impu.label = "impurity"

    # weight different ion quantities based on their density
    for (index, ion2) in ((bulk_index, bulk), (impu_index, impu))
        ion2.element[1].z_n = 0.0
        ion2.element[1].a = 0.0
        for ix in index # z_average is tricky since it's a single constant for the whole profile
            ion2.element[1].z_n += sum(zs[ix] .* ratios[:, ix]) / length(ratios[:, ix])
            ion2.element[1].a += sum(as[ix] .* ratios[:, ix]) / length(ratios[:, ix])
        end
        for item in (:density_thermal, :temperature, :rotation_frequency_tor)
            value = rho_tor_norm .* 0.0
            IMAS.setraw!(ion2, item, value)
            for ix in index
                tp = Vector{typeof(ions[ix]).parameters[1]}
                tmp = getproperty(ions[ix], item, tp())::tp
                if !isempty(tmp)
                    value .+= tmp .* ratios[:, ix]
                end
            end
        end
    end

    return ions2
end

"""
    lump_ions_as_bulk_and_impurity(ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion}, rho_tor_norm::Vector{<:Real})

Changes core_profiles.ion to 2 species, one bulk species (H, D, T) and one combined impurity species
"""
function lump_ions_as_bulk_and_impurity(cp1d::IMAS.core_profiles__profiles_1d{T})::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion{T}} where {T<:Real}

    ions = cp1d.ion

    if length(ions) < 2
        error("lump_ions_as_bulk_and_impurity requires at least two ion species")
    end

    zs = [ion.element[1].z_n for ion in ions]
    as = [ion.element[1].a for ion in ions]
    bulk_index = findall(zs .== 1)
    impu_index = findall(zs .!= 1)

    rho_tor_norm = cp1d.grid.rho_tor_norm
    ratios = zeros(length(rho_tor_norm), length(ions))
    for index in (bulk_index, impu_index)
        ntot = zeros(length(rho_tor_norm))
        for ix in index
            tmp = ions[ix].density_thermal
            ratios[:, ix] = tmp
            ntot .+= tmp
        end
        for ix in index
            ratios[:, ix] ./= ntot
        end
    end

    ne = cp1d.electrons.density_thermal
    n1 = zero(ne)
    for ion in cp1d.ion
        if ion.element[1].z_n == 1
            n1 .+= ion.density_thermal
        end
    end
    Zi = (n1 .- ne .* cp1d.zeff) ./ (n1 .- ne)
    Zi = sum(Zi) / length(Zi) # make Zi constant as function of radius
    ni = (ne .- n1) ./ Zi

    ions2 = IMAS.IDSvector{IMAS.core_profiles__profiles_1d___ion{T}}()

    # bulk ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    bulk = ions2[end]
    resize!(bulk.element, 1)
    bulk.label = "bulk"
    bulk.element[1].z_n = 1.0
    IMAS.setraw!(bulk, :density_thermal, n1)

    # impurity ions
    push!(ions2, IMAS.core_profiles__profiles_1d___ion{T}())
    impu = ions2[end]
    resize!(impu.element, 1)
    impu.label = "impurity"
    impu.element[1].z_n = Zi
    IMAS.setraw!(impu, :density_thermal, ni)

    # weight different ion quantities based on their density
    for (index, ion2) in ((bulk_index, bulk), (impu_index, impu))
        ion2.element[1].a = 0.0
        for ix in index
            ion2.element[1].a += as[ix] * sum(ratios[:, ix]) / length(ratios[:, ix])
        end
        for item in (:temperature, :rotation_frequency_tor)
            value = rho_tor_norm .* 0.0
            IMAS.setraw!(ion2, item, value)
            for ix in index
                tmp = getproperty(ions[ix], item, T[])
                if !isempty(tmp)
                    value .+= tmp .* ratios[:, ix]
                end
            end
        end
    end

    return ions2
end

"""
    zeff(cp1d::IMAS.core_profiles__profiles_1d; temperature_dependent_ionization_state::Bool=true)

Returns plasma effective charge

`temperature_dependent_ionization_state` evaluates Zeff with average ionization state of an ion at a given temperature
"""
function zeff(cp1d::IMAS.core_profiles__profiles_1d; temperature_dependent_ionization_state::Bool=true)
    num = zero(cp1d.grid.rho_tor_norm)
    den = zero(cp1d.grid.rho_tor_norm)
    for ion in cp1d.ion
        if temperature_dependent_ionization_state
            Zi = avgZ(ion.element[1].z_n, ion.temperature)
        else
            Zi = ion.element[1].z_n
        end
        num .+= ion.density .* Zi .^ 2
        den .+= ion.density .* Zi
    end
    return num ./ cp1d.electrons.density
end

"""
    avgZ(Z::Float64,Ti::T)::T

Returns average ionization state of an ion at a given temperature
"""
function avgZ(Z::Float64, Ti::T)::T where {T}
    func = avgZinterpolator(joinpath(dirname(dirname(pathof(@__MODULE__))), "data", "Zavg_z_t.dat"))
    return 10.0 .^ (func.(log10.(Ti ./ 1E3), Z)) .- 1.0
end

Memoize.@memoize function avgZinterpolator(filename::String)
    txt = open(filename, "r") do io
        return strip(read(io, String))
    end

    iion = Vector{Int}()
    Ti = Vector{Float64}()
    for (k, line) in enumerate(split(txt, "\n"))
        if k == 1
            continue
        elseif k == 2
            continue
        elseif k == 3
            append!(iion, collect(map(x -> parse(Int, x), split(line))))
        elseif k == 4
            append!(Ti, collect(map(x -> parse(Float64, x), split(line))))
        end
    end

    data = zeros(length(Ti), length(iion))
    for (k, line) in enumerate(split(txt, "\n"))
        if k > 4
            data[k-4, :] = map(x -> parse(Float64, x), split(line))
        end
    end

    return Interpolations.extrapolate(Interpolations.interpolate((log10.(Ti), iion), log10.(data .+ 1.0), Interpolations.Gridded(Interpolations.Linear())), Interpolations.Flat())

end

"""
    t_i_average(cp1d::IMAS.core_profiles__profiles_1d)::Vector{<:Real}

Returns the average ion temperature weighted by each species density over the total number of ion particles
"""
function t_i_average(cp1d::IMAS.core_profiles__profiles_1d)::Vector{<:Real}
    n = length(cp1d.grid.rho_tor_norm)
    t_i_a = zeros(n)
    ntot = zeros(n)
    for ion in cp1d.ion
        t_i_a += ion.density_thermal .* ion.temperature
        ntot += ion.density_thermal
    end
    return t_i_a ./= ntot
end