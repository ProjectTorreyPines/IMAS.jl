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
    pressure_avg = integrate(volume_cp, pressure) / volume_cp[end]
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
    return 3.0 / 2.0 * integrate(cp1d.grid.volume, cp1d.pressure_thermal)
end

function ne_vol_avg(cp1d::IMAS.core_profiles__profiles_1d)
    return integrate(cp1d.grid.volume, cp1d.electrons.density) / cp1d.grid.volume[end]
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
        integrate(cp1d.grid.volume, sum(ion.density .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) / integrate(cp1d.grid.volume, sum(ion.density for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = vacuum_r0_b0(eqt)

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
        integrate(cp1d.grid.volume, sum(ion.density .* ion.element[1].a for ion in cp1d.ion if ion.element[1].z_n == 1.0)) / integrate(cp1d.grid.volume, sum(ion.density for ion in cp1d.ion if ion.element[1].z_n == 1.0))

    R0, B0 = vacuum_r0_b0(eqt)

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
    return (eqt.global_quantities.ip / 1e6) / (pi * eqt.boundary.minor_radius^2) * 1e20
end

function greenwald_fraction(dd::IMAS.dd)
    return greenwald_fraction(dd.equilibrium.time_slice[], dd.core_profiles.profiles_1d[])
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

function geometric_midplane_line_averaged_density(eqt::IMAS.equilibrium__time_slice, ne_profile::AbstractVector{<:Real}, rho_ne::AbstractVector{<:Real})
    a_cp = interp1d(eqt.profiles_1d.rho_tor_norm, (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) / 2.0).(rho_ne)
    return integrate(a_cp, ne_profile) / a_cp[end]
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
    xpsi = LinRange(0.0, 1.0, ngrid)

    w_E1 = 0.5 * widthp  # width as defined in eped
    xphalf = 1.0 - w_E1
    pconst = 1.0 - tanh((1.0 - xphalf) / w_E1)
    a_t = 2.0 * (ped - edge) / (1.0 + tanh(1.0) - pconst)
    coretanh = 0.5 * a_t * (1.0 - tanh(-xphalf / w_E1) - pconst) + edge

    # edge tanh part
    val = @. 0.5 * a_t * (1.0 - tanh((xpsi - xphalf) / w_E1) - pconst) + edge + core * 0.0

    # core tanh+polynomial part
    if core >= 0.0
        xped = xphalf - w_E1
        xtoped = xpsi ./ xped
        for i in 1:ngrid
            if xtoped[i] < 1.0
                @inbounds val[i] += (core - coretanh) * (1.0 - xtoped[i]^expin)^expout
            end
        end
    end

    return val
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
    if (1 + rtol) > (Ne / Nis) > (1 - rtol)
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

Enforce quasi neutrality by using density_thermal of species and makes sure density is set to the original expression
"""
function enforce_quasi_neutrality!(cp1d::IMAS.core_profiles__profiles_1d, species::Symbol)
    # Make sure expression is used for density
    empty!(cp1d.electrons, :density)
    for ion in cp1d.ion
        empty!(ion, :density)
    end
    species_indx = findfirst(Symbol(ion.label) == species for ion in cp1d.ion)
    @assert species_indx !== nothing
    cp1d.ion[species_indx].density_thermal = (cp1d.electrons.density .+ cp1d.ion[species_indx].density .* cp1d.ion[species_indx].z_ion .- sum(ion.density .* ion.z_ion for ion in cp1d.ion)) ./ cp1d.ion[species_indx].z_ion
end

"""
    lump_ions_as_bulk_and_impurity!(ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion})

Changes core_profiles.ion to 2 species, bulk specie (H, D, T) and combined impurity specie by weigthing masses and densities 
"""
function lump_ions_as_bulk_and_impurity!(ions::IMAS.IDSvector{<:IMAS.core_profiles__profiles_1d___ion}, rho_tor_norm::Vector{<:Real})
    if length(ions) < 2
        error("TAUENN requires two ion species to run")
        # elseif any(!ismissing(ion, :density_fast) for ion in ions)
        #     error("lump_ions_as_bulk_and_impurity! is not setup for handling fast ions")
    elseif length(ions) == 2
        return ions
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

    # bulk ions
    push!(ions, IMAS.core_profiles__profiles_1d___ion())
    bulk = ions[end]
    resize!(bulk.element, 1)
    bulk.label = "bulk"

    # impurity ions
    push!(ions, IMAS.core_profiles__profiles_1d___ion())
    impu = ions[end]
    resize!(impu.element, 1)
    impu.label = "impurity"

    # weight different ion quantities based on their density
    for (index, ion) in ((bulk_index, bulk), (impu_index, impu))
        ion.element[1].z_n = 0.0
        ion.element[1].a = 0.0
        for ix in index # z_average is tricky since it's a single constant for the whole profile
            ion.element[1].z_n += sum(zs[ix] .* ratios[:, ix]) / length(ratios[:, ix])
            ion.element[1].a += sum(as[ix] .* ratios[:, ix]) / length(ratios[:, ix])
        end
        for item in (:density_thermal, :temperature, :rotation_frequency_tor)
            value = rho_tor_norm .* 0.0
            IMAS.setraw!(ion, item, value)
            for ix in index
                if !ismissing(ions[ix], item)
                    value .+= getproperty(ions[ix], item) .* ratios[:, ix]
                end
            end
        end
    end

    for k in reverse(1:length(ions)-2)
        deleteat!(ions, k)
    end

    return ions
end

"""
    avgZ(Z::Float64,Ti::T)::T

Returns average ionization state of an ion at a given temperature
"""
function avgZ(Z::Float64, Ti::T,)::T where {T}
    return 10.0 .^ (avgZinterpolator(joinpath(dirname(dirname(pathof(@__MODULE__))), "data", "Zavg_z_t.dat")).(log10.(Ti ./ 1E3), Z)) .- 1.0
end

Memoize.@memoize function avgZinterpolator(filename::String)
    txt = open(filename, "r") do io
        strip(read(io, String))
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

    return Interpolations.extrapolate(Interpolations.interpolate((log10.(Ti), iion), log10.(data .+ 1), Interpolations.Gridded(Interpolations.Linear())), Interpolations.Flat())

end

"""
    t_i_average(cp1d::IMAS.core_profiles__profiles_1d)::Vector{<:Real}

Returns the average ion temperature weighted by each species density over the total number of ion particles
"""
function t_i_average(cp1d::IMAS.core_profiles__profiles_1d)::Vector{<:Real}
    t_i_a = zeros(length(cp1d.ion[1].density_thermal))
    ntot = zeros(length(cp1d.ion[1].density_thermal))
    for ion in cp1d.ion
        t_i_a += ion.density_thermal .* ion.temperature
        ntot += ion.density_thermal
    end
    return t_i_a ./= ntot
end