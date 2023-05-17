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
    labels = [ion.label for ion in cp1d.ion]
    volume = cp1d.grid.volume

    if "DTH" ∈ labels
        error("you can't have DTH in ion for this function")
    elseif "DT" ∈ labels && "D" ∈ labels || "DT" ∈ labels && "T" ∈ labels
        error("combine your DT and D or DT and T species")
    end

    if "H" ∉ labels
        nh = 0.0
    else
        nh = integrate(volume, cp1d.ion[findfirst(==("H"), labels)].density)
    end

    if "D" ∈ labels
        nd = integrate(volume, cp1d.ion[findfirst(==("D"), labels)].density)
    else
        nd = 0.0
    end

    if "T" ∈ labels
        nt = integrate(volume, cp1d.ion[findfirst(==("T"), labels)].density)
    else
        nt = 0.0
    end

    if "DT" ∈ labels
        nd = nt = integrate(volume, cp1d.ion[findfirst(==("DT"), labels)].density) / 2
    end

    return (nh + 2nd + 3nt) / (nh + nd + nt)
end
