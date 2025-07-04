document[Symbol("Physics nuclear")] = Symbol[]

# D+T→He4   (3.518 MeV)  + n (14.072 MeV) + 17.59MeV
#
# D+D→He3   (0.8175 MeV) + n (2.4525 MeV) + 3.27MeV
# D+D→T     (1.0075 MeV) + H (3.0225 MeV) + 4.03MeV
#
# D+He3→He4 (3.66 MeV)   + H (14.64 MeV) + 18.3MeV

"""
    reactivity(Ti::AbstractVector{<:Real}, model::String; polarized_fuel_fraction::Real=0.0)

Fusion reactivity coming from H.-S. Bosch and G.M. Hale, Nucl. Fusion 32 (1992) 611.
Model can be ["D+T→He4", "D+He3→He4", "D+D→T", "D+D→He3"]")
"""
function reactivity(Ti::AbstractVector{<:Real}, model::String; polarized_fuel_fraction::Real=0.0)
    @assert 0.0 <= polarized_fuel_fraction <= 1.0 "Polarized fuel fraction should be between 0.0 and 1.0"
    spf = 1.0 # default value for non spin polarized fuel
    if model == "D+T→He4"
        # Table VII
        c1 = 1.17302e-9
        c2 = 1.51361e-2
        c3 = 7.51886e-2
        c4 = 4.60643e-3
        c5 = 1.3500e-2
        c6 = -1.06750e-4
        c7 = 1.36600e-5
        bg = 34.3827
        er = 1.124656e6
        if polarized_fuel_fraction > 0.0
            spf = 1.5 #spin polarization factor - 1.5 chosen according to GACP 20010393
        end
    elseif model == "D+He3→He4"
        bg = 68.7508
        mc2 = 1124572.0
        c1 = 5.51036e-10
        c2 = 6.41918e-3
        c3 = -2.02896e-3
        c4 = -1.91080e-5
        c5 = 1.35776e-4
        c6 = 0.0
        c7 = 0.0
        er = 18.3e6
        if polarized_fuel_fraction > 0.0
            spf = 1.5 #also 1.5 for D+He3→He4 according to Kulsrud (1982), PRL 49(17), 1248-1251
        end
    elseif model == "D+D→T"
        bg = 31.3970
        mc2 = 937814.0
        c1 = 5.65718e-12
        c2 = 3.41267e-3
        c3 = 1.99167e-3
        c4 = 0.0
        c5 = 1.05060e-5
        c6 = 0.0
        c7 = 0.0
        er = 4.03e6
        if polarized_fuel_fraction > 0.0
            error("Sorry, spin polarized fuel option is not available for $(model)")
        end
    elseif model == "D+D→He3"
        bg = 31.3970
        mc2 = 937814.0
        c1 = 5.43360e-12
        c2 = 5.85778e-3
        c3 = 7.68222e-3
        c4 = 0.0
        c5 = -2.96400e-6
        c6 = 0.0
        c7 = 0.0
        er = 0.82e6
        if polarized_fuel_fraction > 0.0
            error("Sorry, spin polarized fuel option is not available for $(model)")
        end
    else
        error("Reactivity model can be either [\"D+T→He4\", \"D+He3→He4\", \"D+D→T\", \"D+D→He3\"]")
    end

    Ti = Ti ./ 1e3  # from eV to keV
    r0 = Ti .* (c2 .+ Ti .* (c4 .+ Ti .* c6)) ./ (1.0 .+ Ti .* (c3 .+ Ti .* (c5 .+ Ti .* c7)))
    if model == "D+D→He3"
        # r0 > 1.0 is undefined, cutoff at r0 of 0.95
        r0 .= ifelse.(r0 .> 0.95, 0.95, r0)
        if any(r0 .== 0.95)
            @warn "r0 in D+D→He3 is incorrect. This is likely a symptom of very high Ti [keV]: $Ti"
        end
    end
    theta = Ti ./ (1.0 .- r0)

    xi = (bg .^ 2 ./ (4.0 .* theta)) .^ (1.0 ./ 3.0)
    sigv = c1 .* theta .* sqrt.(xi ./ (er .* Ti .^ 3)) .* exp.(-3.0 .* xi)

    return ((1.0 .- polarized_fuel_fraction) .+ spf * polarized_fuel_fraction) .* sigv / 1e6  # m^3/s
end

@compat public reactivity
push!(document[Symbol("Physics nuclear")], :reactivity)

#===========#
# REACTIONS #
#===========#
"""
    D_T_to_He4_reactions(cp1d::IMAS.core_profiles__profiles_1d)

Calculates the number of D-T thermal fusion reactions to He4 in [reactions/m³/s]
"""
function D_T_to_He4_reactions(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)
    ion_list = (ion.label for ion in cp1d.ion)
    result = zero(cp1d.electrons.density_thermal)

    if "D" in ion_list && "T" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density_thermal
        T_index = findfirst(ion -> isequal(ion, "T"), ion_list)
        n_tritium = cp1d.ion[T_index].density_thermal
        Ti = (cp1d.ion[D_index].temperature .+ cp1d.ion[T_index].temperature) ./ 2.0
        sigv = reactivity(Ti, "D+T→He4"; polarized_fuel_fraction)
        result .= n_deuterium .* n_tritium .* sigv  #  reactions/m³/s

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = n_tritium = cp1d.ion[DT_index].density_thermal ./ 2.0
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D+T→He4"; polarized_fuel_fraction)
        result .= n_deuterium .* n_tritium .* sigv  #  reactions/m³/s
    end

    return result
end

@compat public D_T_to_He4_reactions
push!(document[Symbol("Physics nuclear")], :D_T_to_He4_reactions)

"""
    D_D_to_He3_reactions(dd::IMAS.dd)

Calculates the number of D-D thermal fusion reactions to He3 in [reactions/m³/s]
"""
function D_D_to_He3_reactions(cp1d::IMAS.core_profiles__profiles_1d)
    ion_list = (ion.label for ion in cp1d.ion)
    result = zero(cp1d.electrons.density_thermal)

    if "D" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density_thermal
        Ti = cp1d.ion[D_index].temperature
        sigv = reactivity(Ti, "D+D→He3")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = cp1d.ion[DT_index].density_thermal ./ 2.0
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D+D→He3")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s
    end

    return result
end

@compat public D_D_to_He3_reactions
push!(document[Symbol("Physics nuclear")], :D_D_to_He3_reactions)

"""
    D_D_to_T_reactions(dd::IMAS.dd)

Calculates the number of D-D thermal fusion reactions to T in [reactions/m³/s]
"""
function D_D_to_T_reactions(cp1d::IMAS.core_profiles__profiles_1d)
    ion_list = (ion.label for ion in cp1d.ion)
    result = zero(cp1d.electrons.density_thermal)

    if "D" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density_thermal
        Ti = cp1d.ion[D_index].temperature
        sigv = reactivity(Ti, "D+D→T")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = cp1d.ion[DT_index].density_thermal ./ 2.0
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D+D→T")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s
    end

    return result
end

@compat public D_D_to_T_reactions
push!(document[Symbol("Physics nuclear")], :D_D_to_T_reactions)

#=========#
# SOURCES #
#=========#
"""
    fusion_reaction_source(
        s1d::IMAS.core_sources__source___profiles_1d,
        reactivity::Vector{<:Real},
        in1::Symbol,
        in2::Symbol,
        out::Symbol,
        eV::Float64
    )

Add a fusion reaction source for two isotopes coming in and one coming out
"""
function fusion_reaction_source(
    s1d::IMAS.core_sources__source___profiles_1d,
    reactivity::Vector{<:Real},
    in1::Symbol,
    in2::Symbol,
    out::Symbol,
    eV::Float64
)

    k = 0

    if in1 == in2
        k += 1
        ion = resize!(s1d.ion, k)[k]
        ion_element!(ion, in1)
        ion.particles = -2.0 * reactivity
    else
        k += 1
        ion = resize!(s1d.ion, k)[k]
        ion_element!(ion, in1)
        ion.particles = -reactivity

        k += 1
        ion = resize!(s1d.ion, k)[k]
        ion_element!(ion, in2)
        ion.particles = -reactivity
    end

    k += 1
    ion = resize!(s1d.ion, k)[k]
    ion_element!(ion, out; fast=true)
    ion.particles = reactivity
    ion.fast_particles_energy = eV

    return nothing
end

@compat public fusion_reaction_source
push!(document[Symbol("Physics nuclear")], :fusion_reaction_source)

"""
    D_T_to_He4_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles; combine_DT::Bool)

Calculates DT fusion heating with an estimation of the alpha slowing down to the ions and electrons, modifies dd.core_sources
"""
function D_T_to_He4_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles; combine_DT::Bool)
    cp1d = cp.profiles_1d[]
    polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)

    name = "D+T→He4"
    eV1 = mks.E_α # 3.518e6
    eV2 = mks.E_n # 14.072e6
    reactivity = D_T_to_He4_reactions(cp1d; polarized_fuel_fraction)
    ion_to_electron_fraction = sivukhin_fraction(cp1d, eV1, 4.0)
    energy = reactivity .* eV1 * mks.e # J/m^3/s = W/m^3
    source = resize!(cs.source, :fusion, "identifier.name" => name; wipe=false)
    new_source(
        source,
        source.identifier.index,
        name,
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume,
        cp1d.grid.area;
        electrons_energy=energy .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=energy .* ion_to_electron_fraction
    )

    if combine_DT
        fusion_reaction_source(source.profiles_1d[], reactivity, :DT, :DT, :He4, eV1)
    else
        fusion_reaction_source(source.profiles_1d[], reactivity, :D, :T, :He4, eV1)
    end

    return source
end

@compat public D_T_to_He4_source!
push!(document[Symbol("Physics nuclear")], :D_T_to_He4_source!)

"""
    D_D_to_He3_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)

Calculates the He-3 heating source from D-D fusion reactions, estimates energy transfer to ions and electrons, modifies dd.core_sources
"""
function D_D_to_He3_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)
    cp1d = cp.profiles_1d[]

    name = "D+D→He3"
    eV1 = 0.8175e6
    eV2 = 2.4525e6
    reactivity = D_D_to_He3_reactions(cp1d)
    ion_to_electron_fraction = sivukhin_fraction(cp1d, eV1, 3.0)
    energy = reactivity .* eV1 * mks.e # J/m^3/s = W/m^3
    source = resize!(cs.source, :fusion, "identifier.name" => name; wipe=false)
    new_source(
        source,
        source.identifier.index,
        name,
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume,
        cp1d.grid.area;
        electrons_energy=energy .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=energy .* ion_to_electron_fraction
    )

    fusion_reaction_source(source.profiles_1d[], reactivity, :D, :D, :He3, eV1)

    return source
end

@compat public D_D_to_He3_source!
push!(document[Symbol("Physics nuclear")], :D_D_to_He3_source!)

"""
    D_D_to_T_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)

Calculates the T heating source from D-D fusion reactions, estimates energy transfer to ions and electrons, modifies dd.core_sources
"""
function D_D_to_T_source!(cs::IMAS.core_sources, cp::IMAS.core_profiles)
    cp1d = cp.profiles_1d[]

    name = "D+D→T"
    eV1 = 1.0075e6
    reactivity = D_D_to_T_reactions(cp1d)
    ion_to_electron_fraction = sivukhin_fraction(cp1d, eV1, 3.0)
    energy = reactivity .* eV1 * mks.e # J/m^3/s = W/m^3
    source = resize!(cs.source, :fusion, "identifier.name" => name; wipe=false)
    new_source(
        source,
        source.identifier.index,
        name,
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume,
        cp1d.grid.area;
        electrons_energy=energy .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=energy .* ion_to_electron_fraction
    )

    fusion_reaction_source(source.profiles_1d[], reactivity / 2.0, :D, :D, :T, eV1)

    name = "D+D→H"
    eV2 = 3.0225e6
    ion_to_electron_fraction = sivukhin_fraction(cp1d, eV2, 1.0)
    energy = reactivity .* eV2 * mks.e # J/m^3/s = W/m^3
    source = resize!(cs.source, :fusion, "identifier.name" => name; wipe=false)
    new_source(
        source,
        source.identifier.index,
        name,
        cp1d.grid.rho_tor_norm,
        cp1d.grid.volume,
        cp1d.grid.area;
        electrons_energy=energy .* (1.0 .- ion_to_electron_fraction),
        total_ion_energy=energy .* ion_to_electron_fraction
    )

    fusion_reaction_source(source.profiles_1d[], reactivity / 2.0, :D, :D, :H, eV2)

    return source
end

@compat public D_D_to_T_source!
push!(document[Symbol("Physics nuclear")], :D_D_to_T_source!)

#========#
# TOTALS #
#========#
"""
    D_T_to_He4_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

Volumetric heating source of He4 particles coming from D-T reactions [W m⁻³]
"""
function D_T_to_He4_heating(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)
    energy = 3.518e6 * mks.e  # Joules
    return D_T_to_He4_reactions(cp1d; polarized_fuel_fraction) .* energy
end

@compat public D_T_to_He4_heating
push!(document[Symbol("Physics nuclear")], :D_T_to_He4_heating)

"""
    D_D_to_He3_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

Volumetric heating source of He3 particles coming from D-D reactions [W m⁻³]
"""
function D_D_to_He3_heating(cp1d::IMAS.core_profiles__profiles_1d)
    energy = 0.8175e6 * mks.e  # Joules
    return D_D_to_He3_reactions(cp1d) .* energy
end

@compat public D_D_to_He3_heating
push!(document[Symbol("Physics nuclear")], :D_D_to_He3_heating)

"""
    D_D_to_T_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

Volumetric heating source of T and H particles coming from D-D reactions [W m⁻³]
"""
function D_D_to_T_heating(cp1d::IMAS.core_profiles__profiles_1d)
    energy = (1.0075e6 + 3.0225e6) * mks.e  # Joules
    return D_D_to_T_reactions(cp1d) .* energy
end

@compat public D_D_to_T_heating
push!(document[Symbol("Physics nuclear")], :D_D_to_T_heating)

"""
    D_T_to_He4_plasma_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)

Total power in He4 from D-T reaction [W]
"""
function D_T_to_He4_plasma_power(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)
    return trapz(cp1d.grid.volume, D_T_to_He4_heating(cp1d; polarized_fuel_fraction))
end

@compat public D_T_to_He4_plasma_power
push!(document[Symbol("Physics nuclear")], :D_T_to_He4_plasma_power)

"""
    D_D_to_He3_plasma_power(cp1d::IMAS.core_profiles__profiles_1d)

Total power in He3 from D-D reaction [W]
"""
function D_D_to_He3_plasma_power(cp1d::IMAS.core_profiles__profiles_1d)
    return trapz(cp1d.grid.volume, D_D_to_He3_heating(cp1d))
end

@compat public D_D_to_He3_plasma_power
push!(document[Symbol("Physics nuclear")], :D_D_to_He3_plasma_power)

"""
    D_D_to_T_plasma_power(cp1d::IMAS.core_profiles__profiles_1d)

Total power in T from D-D reaction [W]
"""
function D_D_to_T_plasma_power(cp1d::IMAS.core_profiles__profiles_1d)
    return trapz(cp1d.grid.volume, D_D_to_T_heating(cp1d))
end

@compat public D_D_to_T_plasma_power
push!(document[Symbol("Physics nuclear")], :D_D_to_T_plasma_power)

"""
    fusion_plasma_power(cp1d::IMAS.core_profiles__profiles_1d)

Calculates the total fusion power in the plasma in [W]

If D+T plasma, then D+D is neglected
"""
function fusion_plasma_power(cp1d::IMAS.core_profiles__profiles_1d)
    cp = parent(parent(cp1d))
    ion_list = (ion.label for ion in cp1d.ion)
    if "T" in ion_list || "DT" in ion_list
        polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)
        tot_pow = D_T_to_He4_plasma_power(cp1d; polarized_fuel_fraction)
    else
        tot_pow = D_D_to_He3_plasma_power(cp1d)
        tot_pow += D_D_to_T_plasma_power(cp1d)
    end
    return tot_pow
end

"""
    fusion_plasma_power(dd::IMAS.dd)
"""
function fusion_plasma_power(dd::IMAS.dd)
    return fusion_plasma_power(dd.core_profiles.profiles_1d[])
end

@compat public fusion_plasma_power
push!(document[Symbol("Physics nuclear")], :fusion_plasma_power)

"""
    fusion_power(cp1d::IMAS.core_profiles__profiles_1d{T}; DD_fusion::Bool=false) where {T<:Real}

Calculates the total fusion power in [W]

If D+T plasma, then D+D is neglected
"""
function fusion_power(cp1d::IMAS.core_profiles__profiles_1d{T}; DD_fusion::Bool=false) where {T<:Real}
    cp = parent(parent(cp1d))
    ion_list = (ion.label for ion in cp1d.ion)
    if "T" in ion_list || "DT" in ion_list
        polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)
        tot_pow = D_T_to_He4_plasma_power(cp1d; polarized_fuel_fraction) * 5.0
    elseif DD_fusion
        tot_pow = D_D_to_He3_plasma_power(cp1d) * 4.0
        tot_pow += D_D_to_T_plasma_power(cp1d)
    else
        tot_pow = zero(T)
    end
    return tot_pow
end

"""
    fusion_power(dd::IMAS.dd; DD_fusion::Bool=false)
"""
function fusion_power(dd::IMAS.dd; DD_fusion::Bool=false)
    return fusion_power(dd.core_profiles.profiles_1d[]; DD_fusion)
end

@compat public fusion_power
push!(document[Symbol("Physics nuclear")], :fusion_power)

"""
    fusion_neutron_power(cp1d::IMAS.core_profiles__profiles_1d)

Calculates the total fusion power in the neutrons [W]

If D+T plasma, then D+D is neglected
"""
function fusion_neutron_power(cp1d::IMAS.core_profiles__profiles_1d)
    cp = parent(parent(cp1d))
    ion_list = (ion.label for ion in cp1d.ion)
    if "T" in ion_list || "DT" in ion_list
        polarized_fuel_fraction = getproperty(cp.global_quantities, :polarized_fuel_fraction, 0.0)
        tot_pow = D_T_to_He4_plasma_power(cp1d; polarized_fuel_fraction) * 4.0
    else
        tot_pow = D_D_to_He3_plasma_power(cp1d) * 3.0
    end
    return tot_pow
end

"""
    fusion_neutron_power(dd::IMAS.dd)
"""
function fusion_neutron_power(dd::IMAS.dd)
    return fusion_neutron_power(dd.core_profiles.profiles_1d[])
end

@compat public fusion_neutron_power
push!(document[Symbol("Physics nuclear")], :fusion_neutron_power)