document[Symbol("Physics pf_active")] = Symbol[]

"""
    area(coil::IMAS.pf_active__coil)

Returns cross sectional area of PF coils
"""
function area(coil::IMAS.pf_active__coil)
    A = 0.0
    for element in coil.element
        oute = outline(element)
        A += area(oute.r, oute.z)
    end
    return A
end

@compat public area
push!(document[Symbol("Physics pf_active")], :area)

"""
    volume(coil::IMAS.pf_active__coil)

Returns volume of PF coils
"""
function volume(coil::IMAS.pf_active__coil)
    V = 0.0
    for element in coil.element
        oute = outline(element)
        A = area(oute.r, oute.z)
        RC, ZC = centroid(oute.r, oute.z)
        V += A * 2π * RC
    end
    return V
end

@compat public volume
push!(document[Symbol("Physics pf_active")], :volume)

"""
    is_ohmic_coil(coil::IMAS.pf_active__coil)

Returns true/false if coil is part of the OH
"""
function is_ohmic_coil(coil::IMAS.pf_active__coil)
    if isempty(coil.function)
        error("`$(location(coil)).function` is not set. You may need to run `IMAS.set_coils_function(dd.pf_active.coil, R0)`.")
    end
    return findfirst(:flux, coil.function) !== nothing
end

@compat public is_ohmic_coil
push!(document[Symbol("Physics pf_active")], :is_ohmic_coil)

"""
    set_coils_function(coils::IDSvector{<:IMAS.pf_active__coil}, R0::Float64; force::Bool=false)

Setup the pf_active.coil[:].function
"""
function set_coils_function(coils::IDSvector{<:IMAS.pf_active__coil}, R0::Float64; force::Bool=false)
    # set coil elements geometry_type attribute
    geometry_types = name_2_index(coils[1].element[1].geometry)
    for coil in coils
        for element in coil.element
            for geometry_type in (:outline, :rectangle, :oblique)
                if geometry_type != :geometry_type && !ismissing(getproperty(element.geometry, geometry_type), :r)
                    element.geometry.geometry_type = geometry_types[geometry_type]
                    break
                end
            end
        end
    end

    # find innermost coil... that should be an OH!
    oh_min_radius = Inf
    for coil in coils
        for element in coil.element
            oute = outline(element)
            oh_min_radius = min(oh_min_radius, minimum(oute.r))
        end
    end

    # define as OH coils all the coils that share that minimum radial coordinate
    oh_does_shaping = true
    for coil in coils
        if force
            empty!(coil.function)
        end
        if isempty(coil.function)
            if (!ismissing(coil, :name) && startswith(uppercase(coil.name), "V")) || (!ismissing(coil, :identifier) && startswith(uppercase(coil.identifier), "V"))
                func = resize!(coil.function, :vertical; wipe=false)
                func.description = "Vertical stability"
                continue
            else
                func = resize!(coil.function, :shaping; wipe=false)
                func.description = "PF"
            end

            min_radius = Inf
            max_radius = 0.0
            for element in coil.element
                oute = outline(element)
                min_radius = min(min_radius, minimum(oute.r))
                max_radius = max(max_radius, maximum(oute.r))
            end

            if (!ismissing(coil, :name) && (contains(uppercase(coil.name), "OH") || contains(uppercase(coil.name), "CS"))) ||
               (!ismissing(coil, :identifier) && (contains(uppercase(coil.identifier), "OH") || contains(uppercase(coil.identifier), "CS"))) || min_radius == oh_min_radius
                func = resize!(coil.function, :flux; wipe=false)
                func.description = "OH"

                # here we assume that if some coils are labeled as OH and are not stacked vertically in the central solenoid
                # then the configuration is like in DIII-D where there are OH compensation coils designed to decouple the
                # flux-swing and the plasma shaping functionalities. This is really a DIII-D thing...
                if max_radius > R0
                    oh_does_shaping = false
                end
            end
        end

        for element in coil.element
            if isempty(element, :turns_with_sign)
                element.turns_with_sign = 1.0
            end
        end
    end

    if oh_does_shaping === false
        for coil in coils
            if !isempty(findall(:flux, coil.function))
                deleteat!(coil.function, :shaping)
            end
        end
    end

    return coils
end

@compat public set_coils_function
push!(document[Symbol("Physics pf_active")], :set_coils_function)

"""
    outline(element::Union{IMAS.pf_active__coil___element{T},IMAS.pf_passive__loop___element{T}}) where {T<:Real}

Returns named tuple with r and z arrays outline, independent of geometry_type used to describe the element
"""
function outline(element::Union{IMAS.pf_active__coil___element{T},IMAS.pf_passive__loop___element{T}}) where {T<:Real}
    geometry_type = index_2_name(element.geometry)[element.geometry.geometry_type]

    if geometry_type == :outline
        oute = element.geometry.outline
        r = StaticArrays.SVector{length(oute.r)}(oute.r)
        z = StaticArrays.SVector{length(oute.z)}(oute.z)

    elseif geometry_type == :rectangle
        rect = element.geometry.rectangle
        rr = rect.r
        zz = rect.z
        Δr = 0.5 * rect.width
        Δz = 0.5 * rect.height
        r = StaticArrays.SVector(rr - Δr, rr + Δr, rr + Δr, rr - Δr)
        z = StaticArrays.SVector(zz - Δz, zz - Δz, zz + Δz, zz + Δz)

    elseif geometry_type == :annulus
        # approximate annulus as square coil that has the same conducting area of the original coil
        ann = element.geometry.annulus
        Aout = π * ann.radius_outer^2
        Ain = π * ann.radius_inner^2
        A = Aout - Ain
        Δr = Δz = sqrt(A) / 2.0
        ar = ann.r
        az = ann.z
        r = StaticArrays.SVector(ar - Δr, ar + Δr, ar + Δr, ar - Δr)
        z = StaticArrays.SVector(az - Δz, az - Δz, az + Δz, az + Δz)

    else
        error("pf_active geometry type `$geometry_type` is not yet supported")
    end

    return (r=r, z=z)
end

@compat public outline
push!(document[Symbol("Physics pf_active")], :outline)