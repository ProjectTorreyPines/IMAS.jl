document[Symbol("Physics radiation")] = Symbol[]

"""
    radiation_losses(sources::IMAS.core_sources)

Evaluate total plasma radiation losses [W] due to bremsstrahlung, synchrotron, and line radiation
"""
function radiation_losses(sources::IMAS.core_sources)
    n2i = name_2_index(sources.source)
    radiation_indices =
        [n2i[name] for name in (:bremsstrahlung, :synchrotron_radiation, :line_radiation, :radiation, :cyclotron_radiation, :cyclotron_synchrotron_radiation, :impurity_radiation)]
    radiation_energy = 0.0
    for source in sources.source
        if source.identifier.index ∈ radiation_indices
            radiation_energy += source.profiles_1d[].electrons.power_inside[end]
        end
    end
    return radiation_energy
end

@compat public radiation_losses
push!(document[Symbol("Physics radiation")], :radiation_losses)

"""
    bremsstrahlung_source!(dd::IMAS.dd)

Calculates approximate NRL Bremsstrahlung radiation source and modifies dd.core_sources
"""
function bremsstrahlung_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    powerDensityBrem = -1.690e-38 .* ne .^ 2 .* cp1d.zeff .* sqrt.(Te)

    source = resize!(dd.core_sources.source, :bremsstrahlung; wipe=false)
    new_source(source, source.identifier.index, "brem", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=powerDensityBrem)
    return dd
end

@compat public bremsstrahlung_source!
push!(document[Symbol("Physics radiation")], :bremsstrahlung_source!)

"""
    rad_sync(ϵ::T, a::T, B0::T, ne::T, Te::T; wall_reflection_coefficient) where {T<:Real}

Synchrotron radiation from Trubnikov, JETP Lett. 16 (1972) 25.0

Transpiled from gacode/tgyro/src/tgyro_rad.f90

See also: Study of heat and synchrotron radiation transport in fusion tokamak plasmas (C. Villar 1997)
"""
function rad_sync(ϵ::T, a::T, B0::T, ne::T, Te::T; wall_reflection_coefficient) where {T<:Real}
    #---------------------------------------------------
    # MKS to CGS
    aspect_ratio = 1 / ϵ
    r_min = a * cgs.m_to_cm # [cm]
    b_ref = B0 * cgs.T_to_Gauss # [G]
    ne = ne / cgs.m³_to_cm³ # [1/cm^3]
    e = cgs.e # [statcoul]
    k = cgs.k # [erg / eV]
    m_e = cgs.me # [g]
    c = cgs.c # [cm / s]
    #---------------------------------------------------
    wpe = sqrt(4.0 * pi * ne * e^2 / m_e)
    wce = e * abs(b_ref) / (m_e * c)
    g = k * Te / (m_e * c^2)
    phi = 60.0 * g^1.5 * sqrt((1.0 - wall_reflection_coefficient) * (1.0 + 1.0 / aspect_ratio / sqrt(g)) / (r_min * wpe^2 / c / wce))
    qsync = m_e / (3.0 * pi * c) * g * (wpe * wce)^2 * phi # [erg/cm^3/s]
    return -qsync * 1E-7 * 1E6 #[W/m^3]
end

@compat public rad_sync
push!(document[Symbol("Physics radiation")], :rad_sync)

"""
    synchrotron_source!(dd::IMAS.dd; wall_reflection_coefficient=0.0)

Calculates synchrotron radiation source and modifies dd.core_sources
"""
function synchrotron_source!(dd::IMAS.dd; wall_reflection_coefficient=0.8)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    eq = dd.equilibrium
    eqt = eq.time_slice[]
    a = eqt.boundary.minor_radius
    R = eqt.global_quantities.magnetic_axis.r
    B0 = abs(@ddtime(eq.vacuum_toroidal_field.b0)) * eq.vacuum_toroidal_field.r0 / R
    ϵ = a ./ R

    # Synchrotron radiation
    powerDensitySync = rad_sync.(ϵ, a, B0, ne, Te; wall_reflection_coefficient)

    source = resize!(dd.core_sources.source, :synchrotron_radiation; wipe=false)
    new_source(source, source.identifier.index, "synch", cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=powerDensitySync)
    return source
end

@compat public synchrotron_source!
push!(document[Symbol("Physics radiation")], :synchrotron_source!)

"""
    line_radiation_source!(dd::IMAS.dd)

Calculates line radiation sources and modifies dd.core_sources
"""
function line_radiation_source!(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d[]
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature

    sources = [for ion in cp1d.ion
        ni = ion.density_thermal
        zi = ion.z_ion
        namei = string(elements[Int(floor(ion.z_ion))].symbol)
        linerad = rad_ion_adas(Te, ne, ni, zi, namei)

        name = "line $namei"
        source = resize!(dd.core_sources.source, :line_radiation, "identifier.name" => name; wipe=false)
        new_source(source, source.identifier.index, name, cp1d.grid.rho_tor_norm, cp1d.grid.volume, cp1d.grid.area; electrons_energy=linerad)
        source
    end]

    return sources
end

@compat public line_radiation_source!
push!(document[Symbol("Physics radiation")], :line_radiation_source!)

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

NOTE: Te in [keV] and output is in `[erg / cm^3 / s]`

12-term Chebyshev polynomial fits to ADAS data

    ln[ Lz(x) ] = sum c_n T_n(x)

where 

    Lz = cooling rate in erg/cm^3/s

    T_n(x) = cos[n*arccos(x)] (Chebyshev polynomials)

    T = electron temperature
 
                 ln(T/T_min)
    x = -1 + 2 --------------- 
               ln(T_max/T_min)

    c_n = tabulated polynomial coefficients for each ion

Acknowledgements:
 - F. Sciortino for providing access to ADAS data via Aurora 
 - T. Pütterich for up-to-data ADAS data
 - T. Odstrčil for help/checking of 2025 updates

References:
 - Open ADAS: https://open.adas.ac.uk
 - T. Pütterich et al 2019 Nucl. Fusion 59 056013

 Notes:
- Lz = Lz_line + Lz_continuum 
- Aurora follows the radiation nomenclature of ADAS (as described here), separating
  "line" and "continuum" radiation. Line radiation basically comes from ADF11 PLT
  files and continuum radiation comes from ADF11 PRB files. Bremsstrahlung is
  included in the continuum term.
- For generation of fit coefficients, see tgyro/tools/radiation
"""
function adas21(Te, name)
    # Min and max values of Te
    t0 = 0.05 # keV
    t1 = 50.0 # keV

    # Chebyshev expansion coefficients
    if name == "W"
        coefficients = [
            -4.098179995282e+01,
            -7.936786403815e-01,
            -4.619044787359e-01,
            -1.456232613687e-01,
            +2.937985817197e-01,
            +5.063695264214e-02,
            -1.062685541668e-01,
            +2.661063713322e-02,
            +8.700125448336e-02,
            -6.988848399826e-02,
            +4.169797142892e-03,
            +4.169797142892e-03
        ]
    elseif name == "B"
        coefficients = [
            -4.824500817201e+01,
            -9.583639579115e-01,
            +1.121730883301e+00,
            -1.860489734455e-01,
            -1.336561575251e-01,
            +1.485164070054e-01,
            -9.072467667403e-02,
            +4.806968067742e-02,
            -2.497690268099e-02,
            +1.137544241417e-02,
            -2.413605670955e-03,
            -2.413605670955e-03
        ]
    elseif name == "He"
        coefficients = [
            -5.134700413297e+01,
            +8.218593504938e-01,
            +4.572258948229e-01,
            -1.864060260034e-01,
            +6.186336465363e-02,
            -1.425022648657e-02,
            -3.228599626925e-04,
            +2.723193135411e-03,
            -1.770156311492e-03,
            +5.783337251866e-04,
            +2.308975845735e-05,
            +2.308975845735e-05
        ]
    elseif name == "Be"
        coefficients = [
            -4.901561102918e+01,
            -5.746662236007e-01,
            +1.102477316209e+00,
            -3.421385730883e-01,
            +3.464880456112e-02,
            +3.087936176685e-02,
            -1.824676177229e-02,
            +2.267469723172e-03,
            +2.551848339714e-03,
            -1.268897504553e-03,
            +2.514022413007e-04,
            +2.514022413007e-04
        ]
    elseif name == "C"
        coefficients = [
            -4.783829635152e+01,
            -8.537270849071e-01,
            +7.710123846565e-01,
            +2.354004692170e-01,
            -4.440438456161e-01,
            +2.927614297369e-01,
            -1.156983192437e-01,
            +1.733993408112e-02,
            +2.170627066061e-02,
            -3.195887243365e-02,
            +2.481800173124e-02,
            +2.481800173124e-02
        ]
    elseif name == "O"
        coefficients = [
            -4.709052058267e+01,
            -7.634895015203e-01,
            +2.965651964690e-01,
            +4.801687007883e-01,
            -1.986647439460e-01,
            -2.151688086625e-01,
            +3.313680950846e-01,
            -2.113807929523e-01,
            +4.535099018906e-02,
            +5.682622389415e-02,
            -7.724691172118e-02,
            -7.724691172118e-02
        ]
    elseif name == "N"
        coefficients = [
            -4.748406078178e+01,
            -6.895926848662e-01,
            +3.403024330483e-01,
            +5.832113228037e-01,
            -5.085976351090e-01,
            +1.385618342393e-01,
            +1.094185690421e-01,
            -1.741954140729e-01,
            +1.367404249744e-01,
            -7.076042601889e-02,
            +1.377707920484e-02,
            +1.377707920484e-02
        ]
    elseif name == "F"
        coefficients = [
            -4.608334126395e+01,
            -2.032569469519e+00,
            +1.174720067392e+00,
            -1.601584693147e-01,
            +2.522457754350e-01,
            -4.105666049867e-01,
            +2.639488083059e-01,
            -1.240057037307e-04,
            -1.483369633643e-01,
            +1.362062267100e-01,
            -4.689893298094e-02,
            -4.689893298094e-02
        ]
    elseif name == "Ne"
        coefficients = [
            -4.616289844396e+01,
            -1.476875032140e+00,
            +8.751481882578e-01,
            -1.554907576477e-01,
            +3.453234274569e-01,
            -4.448210670608e-01,
            +2.235707626168e-01,
            +6.309160185099e-02,
            -1.804901948539e-01,
            +1.138691735006e-01,
            +1.238414608039e-02,
            +1.238414608039e-02
        ]
    elseif name == "Al"
        coefficients = [
            -4.486524782614e+01,
            -2.322475619179e+00,
            +9.328274092632e-01,
            +4.744804675756e-02,
            -6.858095854371e-02,
            +2.345581654349e-01,
            -3.716843250349e-01,
            +1.795673265946e-01,
            +1.194090655176e-01,
            -2.087103454071e-01,
            +8.992619526224e-02,
            +8.992619526224e-02
        ]
    elseif name == "Si"
        coefficients = [
            -4.473392580102e+01,
            -2.118859872224e+00,
            +7.535819784959e-01,
            +1.199722562015e-01,
            -1.300605991542e-01,
            +3.223304318098e-01,
            -3.674984537405e-01,
            +5.351426696879e-02,
            +2.108242428513e-01,
            -1.792043287525e-01,
            +9.913893641685e-03,
            +9.913893641685e-03
        ]
    elseif name == "Ar"
        coefficients = [
            -4.425505195350e+01,
            -1.616269823924e+00,
            +8.194549090984e-02,
            +4.830273331631e-01,
            -2.675106132694e-01,
            +1.174875472786e-01,
            +1.676830849068e-01,
            -2.891215850568e-01,
            +4.842130798317e-02,
            +1.715110727683e-01,
            -1.077751857577e-01,
            -1.077751857577e-01
        ]
    elseif name == "Ca"
        coefficients = [
            -4.410322851491e+01,
            -1.396392282770e+00,
            +6.228605606294e-02,
            +2.603343415878e-01,
            +6.868991709327e-02,
            -3.101682260197e-01,
            +5.258665807121e-01,
            -3.702755461767e-01,
            -1.771137529772e-02,
            +1.469773755999e-01,
            -1.213876050533e-02,
            -1.213876050533e-02
        ]
    elseif name == "Fe"
        coefficients = [
            -4.289975750516e+01,
            -2.074724550183e+00,
            +2.592709540821e-01,
            +2.670554581981e-01,
            -3.245193180198e-02,
            -3.447451297355e-02,
            -1.365795512615e-01,
            +2.145851896399e-01,
            +2.140105004146e-02,
            -1.785057925001e-01,
            +5.091860663164e-02,
            +5.091860663164e-02
        ]
    elseif name == "Ni"
        coefficients = [
            -4.281682927312e+01,
            -1.989503354777e+00,
            +4.039460832605e-01,
            +2.224196803176e-01,
            -1.451323786293e-01,
            +7.644152992866e-02,
            -1.717626198140e-01,
            +1.431632992224e-01,
            +1.513611521059e-01,
            -2.308718720128e-01,
            -9.023390843523e-03,
            -9.023390843523e-03
        ]
    elseif name == "Kr"
        coefficients = [
            -4.248982974133e+01,
            -1.299186536746e+00,
            -4.512694457002e-01,
            +6.804566946335e-01,
            -6.310168823628e-03,
            -1.033700187711e-01,
            +2.624047807761e-02,
            -1.531016488962e-01,
            +1.229092533413e-01,
            +5.721604110683e-02,
            -3.858573093198e-02,
            -3.858573093198e-02
        ]
    elseif name == "Mo"
        coefficients = [
            -4.192185652052e+01,
            -1.792215448077e+00,
            +2.732017007262e-03,
            +1.000290040334e-01,
            +4.133083599928e-01,
            -1.495307361845e-01,
            -1.026532480504e-01,
            +1.359298888695e-01,
            -1.749938253055e-01,
            +7.007588723552e-02,
            +8.557768114101e-02,
            +8.557768114101e-02
        ]
    elseif name == "Xe"
        coefficients = [
            -4.136087494072e+01,
            -1.651105892275e+00,
            -2.928373643947e-01,
            +2.928170785526e-01,
            -5.910859658247e-02,
            +1.267994598292e-01,
            +2.256208498121e-02,
            -4.369518838603e-02,
            +1.738552444143e-02,
            +3.245918863202e-03,
            -6.535620998184e-02,
            -6.535620998184e-02
        ]
    elseif name == "Li"
        coefficients = [
            -5.007779179318e+01,
            +1.379180885582e-01,
            +7.886186847887e-01,
            -2.784155208154e-01,
            +5.393512924579e-02,
            +8.568540178773e-03,
            -8.339965761965e-03,
            -3.235331808761e-03,
            +7.058675207572e-03,
            -3.451108596954e-03,
            -2.663068572132e-03,
            -2.663068572132e-03
        ]
    elseif name == "H"
        coefficients = [
            -5.313111260147e+01,
            +1.418007193691e+00,
            +1.261973276523e-01,
            -4.362701135592e-02,
            +1.071044131666e-02,
            -3.046777668538e-03,
            +7.410284109890e-04,
            +5.926084113166e-05,
            +4.909076663983e-05,
            -5.820357621765e-05,
            +8.674684934388e-05,
            +8.674684934388e-05
        ]
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

    # Lz (cooling rate) in [erg / cm^3 / s]
    Lz = exp.(s)

    return Lz
end

@compat public adas21
push!(document[Symbol("Physics radiation")], :adas21)
