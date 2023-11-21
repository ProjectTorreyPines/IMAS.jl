
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
    is_ohmic_coil(coil::IMAS.pf_active__coil)

Returns true/false if coil is part of the OH
"""
function is_ohmic_coil(coil::IMAS.pf_active__coil)
    if isempty(coil.function)
        error("`$(location(coil)).function` is not set. You may need to run `IMAS.set_coils_function(dd.pf_active)`.")
    end
    return findfirst(:flux, coil.function) !== nothing
end

"""
    set_coils_function(coils::IDSvector{<:IMAS.pf_active__coil})

Setup the pf_active.coil[:].function
"""
function set_coils_function(coils::IDSvector{<:IMAS.pf_active__coil})
    # find innermost coil
    min_radius = Inf
    for coil in coils
        min_radius = min(min_radius, coil.element[1].geometry.rectangle.r)
    end
    # define as OH coils all the coils that share that minimum radial coordinate
    for coil in coils
        func = resize!(coil.function, :b_field_shaping; wipe=false)
        func.description = "PF"
        if coil.element[1].geometry.rectangle.r == min_radius
            func = resize!(coil.function, :flux; wipe=false)
            func.description = "OH"
        else
            deleteat!(coil.function, :flux)
        end
    end
    return coils
end
