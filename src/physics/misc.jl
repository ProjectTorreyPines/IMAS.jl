"""
    area(coil::IMAS.pf_active__coil)

Returns cross sectional area of PF coils
"""
function area(coil::IMAS.pf_active__coil)
    return coil.element[1].geometry.rectangle.width * coil.element[1].geometry.rectangle.height
end

"""
    volume(coil::IMAS.pf_active__coil)

Returns volume of PF coils
"""
function volume(coil::IMAS.pf_active__coil)
    return area(coil) * 2π * coil.element[1].geometry.rectangle.r
end

"""
    elongation_limit(R0_over_a::Real)

Returns elongation limit due to control limit from simple aspect ratio scaling
"""
function elongation_limit(R0_over_a::Real)
    return 2.43 + 65.0 * exp(-R0_over_a / 0.376)
end

"""
    elongation_limit(eqt::IMAS.equilibrium__time_slice)

Returns elongation limit due to control limit from simple aspect ratio scaling
"""
function elongation_limit(eqt::IMAS.equilibrium__time_slice)
    return elongation_limit(eqt.global_quantities.magnetic_axis.r / eqt.boundary.minor_radius)
end

"""
    A_effective(cp1d::IMAS.core_profiles__profiles_1d)

A_effective towards L to H scaling see G. Birkenmeier et al 2022 Nucl. Fusion 62 086005
"""
function A_effective(cp1d::IMAS.core_profiles__profiles_1d)
    numerator = []
    denominator = []
    for ion in cp1d.ion
        if ion.element[1].z_n == 1
            n_int = integrate(cp1d.grid.volume, ion.density)
            push!(numerator, n_int * ion.element[1].a)
            push!(denominator, n_int)
        end
    end
    return reduce(+, numerator) / reduce(+, denominator)
end


"""
    scaling_L_to_H_power(dd)

L to H transition power scaling for metal walls and isotope effect according to : G. Birkenmeier et al 2022 Nucl. Fusion 62 086005
returns power in W
"""
function scaling_L_to_H_power(dd)
    return 1e6 * 0.8 * 2 / A_effective(dd.core_profiles.profiles_1d[]) *
           (0.049 * (dd.summary.volume_average.n_e.value[end] / 1e20)^0.72 *
            dd.equilibrium.vacuum_toroidal_field.b0[end]^0.8 *
            dd.equilibrium.time_slice[].profiles_1d.surface[end]^0.94)
end


