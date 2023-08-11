"""
    radiation_losses(sources::IMAS.core_sources)

Evaluate total plasma radiation losses [W] due to bremsstrahlung, synchrotron, and line radiation
"""
function radiation_losses(sources::IMAS.core_sources)
    n2i = name_2_index(sources.source)
    radiation_indices = [n2i[name] for name in (:bremsstrahlung, :synchrotron_radiation, :line_radiation)]
    radiation_energy = 0.0
    for source in sources.source
        if source.identifier.index ∈ radiation_indices
            radiation_energy += source.profiles_1d[].electrons.power_inside[end]
        end
    end
    return radiation_energy
end

"""
    bremsstrahlung_source!(dd::IMAS.dd)

Calculates approximate NRL Bremsstrahlung radiation source and modifies dd.core_sources
"""
function bremsstrahlung_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    powerDensityBrem = -1.690e-38 .* ne .^ 2 .* cp1d.zeff .* sqrt.(Te)

    source = resize!(dd.core_sources.source, :bremsstrahlung; wipe=false, allow_multiple_matches=true)
    new_source(source, source.identifier.index, "brem", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=powerDensityBrem)
    return dd
end

"""
    rad_sync(ϵ::T, a::T, B0::T, ne::T, Te::T; wall_reflection_coefficient=0.8) where {T<:Real}

Synchrotron radiation from Trubnikov, JETP Lett. 16 (1972) 25.0
Transpiled from gacode/tgyro/src/tgyro_rad.f90
See also: Study of heat and synchrotron radiation transport in fusion tokamak plasmas (C. Villar 1997)
"""
function rad_sync(ϵ::T, a::T, B0::T, ne::T, Te::T; wall_reflection_coefficient=0.8) where {T<:Real}
    #---------------------------------------------------
    # MKS to CGS
    aspect_ratio = 1 / ϵ
    r_min = a * gacode_units.m_to_cm # [cm]
    b_ref = B0 * gacode_units.T_to_Gauss # [G]
    ne = ne * gacode_units.m³_to_cm³ # [1/cm^3]
    e = gacode_units.e # [statcoul]
    k = gacode_units.k # [erg / eV]
    m_e = gacode_units.me # [g]
    c = gacode_units.c # [cm / s]
    #---------------------------------------------------
    wpe = sqrt(4.0 * pi * ne * e^2 / m_e)
    wce = e * abs(b_ref) / (m_e * c)
    g = k * Te / (m_e * c^2)
    phi = 60.0 * g^1.5 * sqrt((1.0 - wall_reflection_coefficient) * (1.0 + 1.0 / aspect_ratio / sqrt(g)) / (r_min * wpe^2 / c / wce))
    qsync = m_e / (3.0 * pi * c) * g * (wpe * wce)^2 * phi # [erg/cm^3/s]
    return -qsync * 1E-7 * 1E6 #[W/m^3]
end

"""
    synchrotron_source!(dd::IMAS.dd)

Calculates Synchrotron radiation source and modifies dd.core_sources
"""
function synchrotron_source!(dd::IMAS.dd; wall_reflection_coefficient=0.8)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    a = eqt.boundary.minor_radius
    R = eqt.global_quantities.magnetic_axis.r
    B0 = abs(@ddtime(eq.vacuum_toroidal_field.b0)) * eq.vacuum_toroidal_field.r0 / R
    ϵ = a ./ R

    # Synchrotron radiation
    powerDensitySync = rad_sync.(ϵ, a, B0, ne, Te; wall_reflection_coefficient)

    source = resize!(dd.core_sources.source, :synchrotron_radiation; wipe=false, allow_multiple_matches=true)
    new_source(source, source.identifier.index, "synch", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=powerDensitySync)
    return source
end

"""
    line_radiation_source!(dd::IMAS.dd)

Calculates line radiation source and modifies dd.core_sources
"""
function line_radiation_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density
    Te = cp1d.electrons.temperature

    linerad = zero(Te)
    for ion in cp1d.ion
        ni = ion.density
        zi = ion.z_ion
        namei = string(elements[Int(floor(ion.z_ion))].symbol)
        linerad .+= rad_ion_adas(Te, ne, ni, zi, namei)
    end

    source = resize!(dd.core_sources.source, :line_radiation; wipe=false, allow_multiple_matches=true)
    new_source(source, source.identifier.index, "line", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=linerad)
    return source
end

function rad_ion_adas(Te, ne, ni, zi, namei)
    ne = ne ./ 1E6 # [1 / cm^3]
    ni = ni ./ 1E6 # [1 / cm^3]

    # * Approximate * NRL Bremsstrahlung radiation [erg / cm^3 / s]
    qbremi = @. 1e7 * 1.69e-32 * ne * ni * zi^2 * sqrt(Te)

    # Chebyshev interpolation of ADAS data
    Lz = adas21(Te / 1E3, namei)

    # ADAS returns combined qcool = qbrem + qline
    qcool = @. ne * ni * Lz

    # Estimate (qline) by sutracting approximate brem
    qline = @. qcool - qbremi

    qline .*= -1e-7 .* 1e6 # in W/m^3
    qcool .*= -1e-7 .* 1e6 # in W/m^3

    return qline
end


"""
    adas21(Te, name)

NOTE: Te in [keV] and output is in [erg cm^3 / s]

Chebyshev polynomial fits to ADAS data
- Transpiled from gacode/tgyro/src/tgyro_rad.f90
- Lz = Lz_line + Lz_continuum 
- Aurora follows the radiation nomenclature of ADAS (as described here), separating
  "line" and "continuum" radiation. Line radiation basically comes from ADF11 PLT
  files and continuum radiation comes from ADF11 PRB files. Bremsstrahlung is
  included in the continuum term.
- Supports ["W","Xe","Mo","Kr","Ni","Fe","Ca","Ar","Si","Al","Ne","F","N","O","C","Be","He","H","T","D","DT"]
"""
function adas21(Te, name)
    # Min and max values of Te
    t0 = 0.05 # keV
    t1 = 50.0 # keV

    # Chebyshev expansion coefficients
    if name == "W"
        coefficients = [-4.093426327035e+01, -8.887660631564e-01, -3.780990284830e-01, -1.950023337795e-01, +3.138290691843e-01, +4.782989513315e-02, -9.942946187466e-02, +8.845089763161e-03, +9.069526573697e-02, -5.245048352825e-02, -1.487683353273e-02, +1.917578018825e-02]
    elseif name == "Xe"
        coefficients = [-4.126366679797e+01, -1.789569183388e+00, -2.380331458294e-01, +2.916911530426e-01, -6.217313390606e-02, +1.177929596352e-01, +3.114580325620e-02, -3.551020007260e-02, -4.850122964780e-03, +1.132323304719e-02, -5.275312157892e-02, -9.051568201374e-03]
    elseif name == "Mo"
        coefficients = [-4.178151951275e+01, -1.977018529373e+00, +5.339155696054e-02, +1.164267551804e-01, +3.697881990263e-01, -9.594816048640e-02, -1.392054581553e-01, +1.272648056277e-01, -1.336366483240e-01, +3.666060293888e-02, +9.586025795242e-02, -7.210209944439e-02]
    elseif name == "Kr"
        coefficients = [-4.235332287815e+01, -1.508707679199e+00, -3.300772886398e-01, +6.166385849657e-01, +1.752687990068e-02, -1.004626261246e-01, +5.175682671490e-03, -1.275380183939e-01, +1.087790584052e-01, +6.846942959545e-02, -5.558980841419e-02, -6.669294912560e-02]
    elseif name == "Ni"
        coefficients = [-4.269403899818e+01, -2.138567547684e+00, +4.165648766103e-01, +2.507972619622e-01, -1.454986877598e-01, +4.044612562765e-02, -1.231313167536e-01, +1.307076922327e-01, +1.176971646853e-01, -1.997449027896e-01, -8.027057678386e-03, +1.583614529900e-01]
    elseif name == "Fe"
        coefficients = [-4.277490044241e+01, -2.232798257858e+00, +2.871183684045e-01, +2.903760139426e-01, -4.662374777924e-02, -4.436273974526e-02, -1.004882554335e-01, +1.794710746088e-01, +3.168699330882e-02, -1.813266337535e-01, +5.762415716395e-02, +6.379542965373e-02]
    elseif name == "Ca"
        coefficients = [-4.390083075521e+01, -1.692920511934e+00, +1.896825846094e-01, +2.333977195162e-01, +5.307786998918e-02, -2.559420140904e-01, +4.733492400000e-01, -3.788430571182e-01, +3.375702537147e-02, +1.030183684347e-01, +1.523656115806e-02, -7.482021324342e-02]
    elseif name == "Ar"
        coefficients = [-4.412345259739e+01, -1.788450950589e+00, +1.322515262175e-01, +4.876947389538e-01, -2.869002749245e-01, +1.699452914498e-01, +9.950501421570e-02, -2.674585184275e-01, +7.451345261250e-02, +1.495713760953e-01, -1.089524173155e-01, -4.191575231760e-02]
    elseif name == "Si"
        coefficients = [-4.459983387390e+01, -2.279998599897e+00, +7.703525425589e-01, +1.494919348709e-01, -1.136851457700e-01, +2.767894295326e-01, -3.577491771736e-01, +7.013841334798e-02, +2.151919651291e-01, -2.052895326141e-01, +2.210085804088e-02, +9.270982150548e-02]
    elseif name == "Al"
        coefficients = [-4.475065090279e+01, -2.455868594007e+00, +9.468903008039e-01, +6.944445017599e-02, -4.550919134508e-02, +1.804382546971e-01, -3.573462505157e-01, +2.075274089736e-01, +1.024482383310e-01, -2.254367207993e-01, +1.150695613575e-01, +3.414328980459e-02]
    elseif name == "Ne"
        coefficients = [-4.599844680574e+01, -1.684860164232e+00, +9.039325377493e-01, -7.544604235334e-02, +2.849631706915e-01, -4.827471944126e-01, +3.138177972060e-01, +2.876874062690e-03, -1.809607030192e-01, +1.510609882754e-01, -2.475867654255e-02, -6.269602018004e-02]
    elseif name == "F"
        coefficients = [-4.595870691474e+01, -2.176917325041e+00, +1.176783264877e+00, -7.712313240060e-02, +1.847534287214e-01, -4.297192280031e-01, +3.374503944631e-01, -5.862051731844e-02, -1.363051725174e-01, +1.580531615737e-01, -7.677594113938e-02, -5.498186771891e-03]
    elseif name == "N"
        coefficients = [-4.719917668483e+01, -1.128938430123e+00, +5.686617156868e-01, +5.565647850806e-01, -6.103105546858e-01, +2.559496676285e-01, +3.204394187397e-02, -1.347036917773e-01, +1.166192946931e-01, -6.001774708924e-02, +1.078186024405e-02, +1.336864982060e-02]
    elseif name == "O"
        coefficients = [-4.688092238361e+01, -1.045540847894e+00, +3.574644442831e-01, +6.007860794100e-01, -3.812470436912e-01, -9.944716626912e-02, +3.141455586422e-01, -2.520592337580e-01, +9.745206757309e-02, +1.606664371633e-02, -5.269687016804e-02, +3.726780755484e-02]
    elseif name == "C"
        coefficients = [-4.752370087442e+01, -1.370806613078e+00, +1.119762977201e+00, +6.244262441360e-02, -4.172077577493e-01, +3.237504483005e-01, -1.421660253114e-01, +2.526893756273e-02, +2.320010310338e-02, -3.487271688767e-02, +2.758311539699e-02, -1.063180164276e-02]
    elseif name == "Be"
        coefficients = [-4.883447566291e+01, -8.543314577695e-01, +1.305444973614e+00, -4.830394934711e-01, +1.005512839480e-01, +1.392590190604e-02, -1.980609625444e-02, +5.342857189984e-03, +2.324970825974e-03, -2.466382923947e-03, +1.073116177574e-03, -9.834117466066e-04]
    elseif name == "He"
        coefficients = [-5.128490291648e+01, +7.743125302555e-01, +4.674917416545e-01, -2.087203609904e-01, +7.996303682551e-02, -2.450841492530e-02, +4.177032799848e-03, +1.109529527611e-03, -1.080271138220e-03, +1.914061606095e-04, +2.501544833223e-04, -3.856698155759e-04]
    elseif name in ("H", "D", "T", "DT")
        # Hydrogen - like ions (H, D, T, DT)
        coefficients = [-5.307012989032e+01, +1.382271913121e+00, +1.111772196884e-01, -3.989144654893e-02, +1.043427394534e-02, -3.038480967797e-03, +5.851591993347e-04, +3.472228652286e-04, -8.418918897927e-05, +3.973067124523e-05, -3.853620366361e-05, +2.005063821667e-04]
    else
        error("No line radiation for $name")
    end

    # Convert Te to Chebyshev grid x = [-1, 1]
    x = @. -1.0 + 2.0 * log(Te / t0) / log(t1 / t0)

    # If outside domain use boundary value
    x[x.<=-1.0] .= -1.0
    x[x.>=1.0] .= 1.0

    # Sum the Chebyshev series where T_n(x) = cos[n * arccos(x)]
    s = zero(x)
    for i in eachindex(coefficients)
        s = s .+ (coefficients[i] .* cos.((i - 1) .* acos.(x)))
    end

    # Lz (cooling rate) in [erg cm^3 / s]
    Lz = exp.(s)

    return Lz
end
