"""
    area(coil::IMAS.pf_active__coil)

returns cross sectional area of PF coils
"""
function area(coil::IMAS.pf_active__coil)
    return coil.element[1].geometry.rectangle.width * coil.element[1].geometry.rectangle.height
end

"""
    volume(coil::IMAS.pf_active__coil)

returns volume of PF coils
"""
function volume(coil::IMAS.pf_active__coil)
    return area(coil) * 2pi * coil.element[1].geometry.rectangle.r
end

"""
    elongation_limit(A::Real)

returns elongation limit due to control limit from simple aspect ratio scaling
"""
function elongation_limit(A::Real)
     return 2.43 + 65.0 * exp(-A / 0.376)
end

function elongation_limit(eqt::IMAS.equilibrium__time_slice)
    return elongation_limit(eqt.global_quantities.magnetic_axis.r/eqt.boundary.minor_radius)
end
