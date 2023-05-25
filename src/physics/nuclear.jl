#### REACTIONS ####
# 
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
    spf = 1 #default value for non spin polarized fuel
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
    theta = Ti ./ (1.0 .- r0)
    xi = (bg .^ 2 ./ (4.0 .* theta)) .^ (1.0 ./ 3.0)
    sigv = c1 .* theta .* sqrt.(xi ./ (er .* Ti .^ 3)) .* exp.(-3.0 .* xi)

    @assert 0.0 <= polarized_fuel_fraction <= 1.0 "Polarized fuel fraction should be between 0.0 and 1.0"
    return ((1.0 .- polarized_fuel_fraction) .+ spf * polarized_fuel_fraction) .* sigv / 1e6  # m^3/s
end

"""
    D_T_to_He4_reactions(cp1d::IMAS.core_profiles__profiles_1d)

Calculates the number of D-T thermal fusion reactions to He4 in [reactions/m³/s]
"""
function D_T_to_He4_reactions(cp1d::IMAS.core_profiles__profiles_1d; polarized_fuel_fraction::Real=0.0)
    ion_list = (ion.label for ion in cp1d.ion)
    result = zero(cp1d.electrons.density)

    if "D" in ion_list && "T" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density
        T_index = findfirst(ion -> isequal(ion, "T"), ion_list)
        n_tritium = cp1d.ion[T_index].density
        Ti = (cp1d.ion[D_index].temperature .+ cp1d.ion[T_index].temperature) ./ 2.0
        sigv = reactivity(Ti, "D+T→He4"; polarized_fuel_fraction)
        result .= n_deuterium .* n_tritium .* sigv  #  reactions/m³/s

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = n_tritium = cp1d.ion[DT_index].density ./ 2.0
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D+T→He4"; polarized_fuel_fraction)
        result .= n_deuterium .* n_tritium .* sigv  #  reactions/m³/s
    end

    return result
end

"""
    D_D_to_He3_reactions(dd::IMAS.dd)

Calculates the number of D-D thermal fusion reactions to He3 in [reactions/m³/s]
"""
function D_D_to_He3_reactions(cp1d::IMAS.core_profiles__profiles_1d)
    ion_list = (ion.label for ion in cp1d.ion)
    result = zero(cp1d.electrons.density)

    if "D" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density
        Ti = cp1d.ion[D_index].temperature
        sigv = reactivity(Ti, "D+D→He3")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = cp1d.ion[DT_index].density ./ 2.0
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D+D→He3")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s
    end

    return result
end

"""
    D_D_to_T_reactions(dd::IMAS.dd)

Calculates the number of D-D thermal fusion reactions to T in [reactions/m³/s]
"""
function D_D_to_T_reactions(cp1d::IMAS.core_profiles__profiles_1d)
    ion_list = (ion.label for ion in cp1d.ion)
    result = zero(cp1d.electrons.density)

    if "D" in ion_list
        D_index = findfirst(ion -> isequal(ion, "D"), ion_list)
        n_deuterium = cp1d.ion[D_index].density
        Ti = cp1d.ion[D_index].temperature
        sigv = reactivity(Ti, "D+D→T")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s

    elseif "DT" in ion_list
        DT_index = findfirst(ion -> isequal(ion, "DT"), ion_list)
        n_deuterium = cp1d.ion[DT_index].density ./ 2.0
        Ti = cp1d.ion[DT_index].temperature
        sigv = reactivity(Ti, "D+D→T")
        result .= n_deuterium .^ 2 .* sigv  #  reactions/m³/s
    end

    return result
end
