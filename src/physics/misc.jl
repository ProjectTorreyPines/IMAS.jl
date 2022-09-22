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
