document[Symbol("Physics build")] = Symbol[]

# BuildLayerType
@enum BuildLayerType::Int _plasma_ = -1 _gap_ _oh_ _tf_ _shield_ _blanket_ _wall_ _vessel_ _cryostat_ _divertor_ _port_
@compat public BuildLayerType
Base.Docs.@doc """
    Enum BuildLayerType

Used for `dd.build.layer[:].type`

* $(join(["`$m` -> $(Int(m))" for m in instances(IMAS.BuildLayerType)],"\n* "))
""" BuildLayerType
push!(document[Symbol("Physics build")], :BuildLayerType)

# BuildLayerSide
@enum BuildLayerSide::Int _lfs_ = -1 _lhfs_ _hfs_ _in_ _out_
@compat public BuildLayerSide
Base.Docs.@doc """
    Enum BuildLayerSide

Used for `dd.build.layer[:].side`

* $(join(["`$m` -> $(Int(m))" for m in instances(IMAS.BuildLayerSide)],"\n* "))
""" BuildLayerSide
push!(document[Symbol("Physics build")], :BuildLayerSide)

# BuildLayerShape
@enum BuildLayerShape::Int _offset_ _negative_offset_ _convex_hull_ _mirror_princeton_D_exact_ _princeton_D_ _mirror_princeton_D_ _princeton_D_scaled_ _mirror_princeton_D_scaled_ _rectangle_ _double_ellipse_ _mirror_double_ellipse_ _rectangle_ellipse_ _mirror_rectangle_ellipse_ _circle_ellipse_ _mirror_circle_ellipse_ _triple_arc_ _mirror_triple_arc_ _miller_ _silo_ _racetrack_ _undefined_
@compat public BuildLayerShape
Base.Docs.@doc """
    Enum BuildLayerShape

Used for `dd.build.layer[:].shape`

* $(join(["`$m` -> $(Int(m))" for m in instances(IMAS.BuildLayerShape)],"\n* "))
""" BuildLayerShape
push!(document[Symbol("Physics build")], :BuildLayerShape)

"""
    build_radii(layers::IMAS.IDSvector{<:IMAS.build__layer{T}}) where {T<:Real}

Convert thicknesses to absolute radii in the build layers
"""
function build_radii(layers::IMAS.IDSvector{<:IMAS.build__layer{T}}) where {T<:Real}
    layers_radii = zeros(T, length(layers) + 1)
    for (k, l) in enumerate(layers)
        @inbounds layers_radii[k+1] = layers_radii[k] + l.thickness
    end
    return layers_radii
end

@compat public build_radii
push!(document[Symbol("Physics build")], :build_radii)

"""
    get_build_indexes(
        layers::IMAS.IDSvector{<:IMAS.build__layer};
        type::Union{Nothing,BuildLayerType}=nothing,
        name::Union{Nothing,String}=nothing,
        identifier::Union{Nothing,Integer}=nothing,
        fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

Returns indexes of layer(s) in build based on a series of selection criteria
"""
function get_build_indexes(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

    valid_layers_indexes = Int[]
    for (k, l) in enumerate(layers)
        if (name === nothing || l.name == name) &&
           (type === nothing || l.type == Int(type)) &&
           (identifier === nothing || l.identifier == identifier) &&
           (fs === nothing || (typeof(fs) <: AbstractVector{BuildLayerSide} && l.side in map(Int, fs)) || (typeof(fs) <: BuildLayerSide && l.side == Int(fs)))
            push!(valid_layers_indexes, k)
        end
    end

    return valid_layers_indexes
end

@compat public get_build_indexes
push!(document[Symbol("Physics build")], :get_build_indexes)

"""
    get_build_index(
        layers::IMAS.IDSvector{<:IMAS.build__layer};
        type::Union{Nothing,BuildLayerType}=nothing,
        name::Union{Nothing,String}=nothing,
        identifier::Union{Nothing,Integer}=nothing,
        fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

Returns index of layer in build based on a series of selection criteria

It raises an error if none or more than one layer matches.
"""
function get_build_index(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

    valid_layers_indexes = get_build_indexes(layers; type, name, identifier, fs)

    if isempty(valid_layers_indexes)
        error("Did not find build.layer: name=$(repr(name)) type=$type identifier=$identifier fs=$fs")
    elseif length(valid_layers_indexes) > 1
        error("Found multiple layers that satisfy name:$name type:$type identifier:$identifier fs:$fs")
    end

    return valid_layers_indexes[1]
end

@compat public get_build_index
push!(document[Symbol("Physics build")], :get_build_index)

"""
    get_build_layers(
        layers::IMAS.IDSvector{<:IMAS.build__layer};
        type::Union{Nothing,BuildLayerType}=nothing,
        name::Union{Nothing,String}=nothing,
        identifier::Union{Nothing,Integer}=nothing,
        fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

Select layer(s) in build based on a series of selection criteria
"""
function get_build_layers(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

    valid_layers_indexes = get_build_indexes(layers; type, name, identifier, fs)

    return layers[valid_layers_indexes]
end

@compat public get_build_layers
push!(document[Symbol("Physics build")], :get_build_layers)

"""
    get_build_layer(
        layers::IMAS.IDSvector{<:IMAS.build__layer};
        type::Union{Nothing,BuildLayerType}=nothing,
        name::Union{Nothing,String}=nothing,
        identifier::Union{Nothing,Integer}=nothing,
        fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

Select layer in build based on a series of selection criteria

It raises an error if none or more than one layer matches.
"""
function get_build_layer(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

    valid_layer_index = get_build_index(layers; type, name, identifier, fs)

    return layers[valid_layer_index]
end

@compat public get_build_layer
push!(document[Symbol("Physics build")], :get_build_layer)

"""
    opposite_side_layer(layer::IMAS.build__layer)

Returns corresponding layer on the high/low field side
"""
function opposite_side_layer(layer::IMAS.build__layer)
    if layer.side in (Int(_hfs_), Int(_lfs_))
        return get_build_layer(parent(layer); identifier=layer.identifier, fs=(layer.side == Int(_lfs_)) ? _hfs_ : _lfs_)
    else
        return layer
    end
end

@compat public opposite_side_layer
push!(document[Symbol("Physics build")], :opposite_side_layer)

"""
    structures_mask(bd::IMAS.build; ngrid::Int = 257, border_fraction::Real = 0.1)

return rmask, zmask, mask of structures that are not vacuum
"""
function structures_mask(bd::IMAS.build; ngrid::Int=257, border_fraction::Real=0.1, layer_check::Function=layer -> layer.material == "vacuum")
    border = maximum(bd.layer[end].outline.r) * border_fraction
    xlim = [0.0, maximum(bd.layer[end].outline.r) + border]
    ylim = [minimum(bd.layer[end].outline.z) - border, maximum(bd.layer[end].outline.z) + border]
    rmask = range(xlim[1], xlim[2], ngrid)
    zmask = range(ylim[1], ylim[2], ngrid * round(Int, (ylim[2] - ylim[1]) / (xlim[2] - xlim[1])))
    mask = ones(length(rmask), length(zmask))

    # start from the first vacuum that goes to zero outside of the TF
    start_from = -1
    for k in get_build_indexes(bd.layer; fs=_out_)
        if layer_check(bd.layer[k]) && minimum(bd.layer[k].outline.r) < bd.layer[1].end_radius
            start_from = k
            break
        end
    end

    # assign boolean mask going through the layers
    valid = true
    for layer in vcat(bd.layer[start_from], bd.layer[1:start_from])
        if layer.type == Int(_plasma_)
            valid = false
        end
        if valid && !ismissing(layer.outline, :r)
            outline = collect(zip(layer.outline.r, layer.outline.z))
            if layer_check(layer) && (layer.side != Int(_in_))
                for (kr, rr) in enumerate(rmask)
                    for (kz, zz) in enumerate(zmask)
                        if PolygonOps.inpolygon((rr, zz), outline) != 0
                            mask[kr, kz] = 0.0
                        end
                    end
                end
            else
                for (kr, rr) in enumerate(rmask)
                    for (kz, zz) in enumerate(zmask)
                        if PolygonOps.inpolygon((rr, zz), outline) == 1
                            mask[kr, kz] = 1.0
                        end
                    end
                end
            end
        end
    end
    rlim_oh = get_build_layer(bd.layer; type=_oh_).start_radius
    for (kr, rr) in enumerate(rmask)
        for (kz, zz) in enumerate(zmask)
            if rr < rlim_oh
                mask[kr, kz] = 1.0
            end
        end
    end
    return rmask, zmask, mask
end

"""
    func_nested_layers(layer::IMAS.build__layer{D}, func::Function)::D where {D<:Real}

Apply function `func` to a layer, then subtract `func` applied to the layer inside of it.

This is used to caclulate `area`, `volume`, etc.. of each layer.
"""
function func_nested_layers(layer::IMAS.build__layer{D}, func::Function)::D where {D<:Real}
    i = index(layer)
    layers = parent(layer)
    # _in_ layers or plasma
    if layer.side ∈ (Int(_in_), Int(_lhfs_))
        return func(layer)
        # anular layers
    elseif layer.side ∈ (Int(_hfs_), Int(_lfs_))
        i = index(layer)
        if layer.side == Int(_hfs_)
            layer_in = layers[i+1]
        else
            layer_in = layers[i-1]
        end
        return func(layer) - func(layer_in)
        # _out_ layers
    elseif layer.side == Int(_out_)
        layer_in = layers[i-1]
        if layer_in.side == Int(_out_)
            return func(layer) - func(layer_in)
        elseif layer_in.side == Int(_lfs_)
            func_in = 0.0
            for l in layers
                if l.side == Int(_in_)
                    func_in += func(l)
                end
            end
            return func(layer) - func(layer_in) - func_in
        end
    end
end

@compat public func_nested_layers
push!(document[Symbol("Physics build")], :func_nested_layers)


"""
    area(layer::IMAS.build__layer)

Calculate cross-sectional area of a build layer
"""
function area(layer::IMAS.build__layer)
    return func_nested_layers(layer, l -> area(l.outline.r, l.outline.z))
end

@compat public area
push!(document[Symbol("Physics build")], :area)

"""
    volume(layer::IMAS.build__layer)

Calculate volume of a build layer revolved around z axis
"""
function volume(layer::IMAS.build__layer)
    if layer.type == Int(_tf_)
        build = parent(parent(layer))
        return func_nested_layers(layer, l -> area(l.outline.r, l.outline.z)) * build.tf.wedge_thickness * build.tf.coils_n
    else
        return func_nested_layers(layer, l -> revolution_volume(l.outline.r, l.outline.z))
    end
end

"""
    volume(layer::IMAS.build__structure)

Calculate volume of a build structure outline revolved around z axis
"""
function volume(structure::IMAS.build__structure)
    if structure.toroidal_extent == 2 * pi
        toroidal_angles = [0.0]
    else
        toroidal_angles = structure.toroidal_angles
    end
    return area(structure.outline.r, structure.outline.z) * structure.toroidal_extent * length(toroidal_angles)
end

@compat public volume
push!(document[Symbol("Physics build")], :volume)

Memoize.@memoize function memoized_inscribed_polygon(pr::AbstractVector{Float64}, pz::AbstractVector{Float64}, target_area_fraction::Float64)
    return inscribed_polygon(pr, pz; target_area_fraction)
end

"""
    first_wall(wall::IMAS.wall{T}; simplify_to_inscribed_fractional_area::Float64=1 - 1E-6) where {T<:Real}

Returns named tuple with outline of the official contiguous first wall limiter contour, or an empty outline if not present

NOTE: in IMAS `wall.description_2d[].limiter.type.index == 0` indicates an official contiguous limiter contour

NOTE: By default the first wall simplifies points on straight lines
"""
function first_wall(wall::IMAS.wall{T}; simplify_to_inscribed_fractional_area::Float64=1 - 1E-6) where {T<:Real}
    for d2d in wall.description_2d
        # d2d.limiter.type.index != 0 indicates a disjoint limiter
        if !ismissing(d2d.limiter.type, :index) && d2d.limiter.type.index != 0
            continue
        end
        for unit in d2d.limiter.unit
            @assert length(d2d.limiter.unit) == 1 # there should be just one official contiguous limiter contour
            oute = closed_polygon(unit.outline.r, unit.outline.z)
            if simplify_to_inscribed_fractional_area < 1.0
                simplified_r, simplified_z = memoized_inscribed_polygon(oute.r, oute.z, simplify_to_inscribed_fractional_area)
                return (r=simplified_r, z=simplified_z)
            else
                return (r=oute.r, z=oute.z)
            end
        end
    end
    return (r=Float64[], z=Float64[])
end

"""
    first_wall(pf_active::IMAS.pf_active{T}) where {T<:Real}

Returns named tuple with outline of the first wall, defined by the pf_active.coils
"""
function first_wall(pf_active::IMAS.pf_active{T}) where {T<:Real}
    if isempty(pf_active.coil)
        return (r=T[], z=T[])
    end

    # find center of the coils
    min_r = Inf
    max_r = -Inf
    min_z = Inf
    max_z = -Inf
    for coil in pf_active.coil
        for element in coil.element
            r, z = outline(element)
            min_r = min(min_r, minimum(r))
            max_r = max(max_r, maximum(r))
            min_z = min(min_z, minimum(z))
            max_z = max(max_z, maximum(z))
        end
    end
    R0 = (min_r + max_r) / 2.0
    Z0 = (min_z + max_z) / 2.0

    # pick points closer to this center
    rce = T[]
    zce = T[]
    for coil in pf_active.coil
        re = T[]
        ze = T[]
        for element in coil.element
            r, z = outline(element)
            append!(re, r)
            append!(ze, z)
        end
        index = argmin(sqrt.((re .- R0) .^ 2 .+ (ze .- Z0) .^ 2))
        push!(rce, re[index])
        push!(zce, ze[index])
    end
    index = sortperm(atan.(zce .- Z0, rce .- R0))
    r, z = closed_polygon(rce[index], zce[index]).rz

    r = T[]
    z = T[]
    for k in 1:length(rce)-1
        r0 = rce[k]
        z0 = zce[k]
        r1 = rce[k+1]
        z1 = zce[k+1]
        push!(r, r0)
        push!(z, z0)
        if r1 >= r0 && z1 <= z0
            push!(r, r1)
            push!(z, z0)
        elseif r1 <= r0 && z1 <= z0
            push!(r, r0)
            push!(z, z1)
        elseif r1 <= r0 && z1 >= z0
            push!(r, r1)
            push!(z, z0)
        elseif r1 >= r0 && z1 >= z0
            push!(r, r0)
            push!(z, z1)
        end
    end
    push!(r, rce[1])
    push!(z, zce[1])

    return (r=r, z=z)
end

"""
    first_wall(eqt::IMAS.equilibrium__time_slice; precision::Float=1E-3) where {T<:Real}

Returns named tuple with outline of the first wall, defined by equilibrium computation domain
"""
function first_wall(eqt::IMAS.equilibrium__time_slice{T}; precision::Float64=1E-3) where {T<:Real}
    # equilibrium compute box
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    if eqt2d === nothing
        return (r=T[], z=T[])
    else
        r_start = eqt2d.grid.dim1[1] + precision
        r_end = eqt2d.grid.dim1[end] - precision
        z_low = eqt2d.grid.dim2[1] + precision
        z_high = eqt2d.grid.dim2[end] - precision
        r = T[r_start, r_start, r_end, r_end, r_start]
        z = T[z_low, z_high, z_high, z_low, z_low]
        return (r=r, z=z)
    end
end

"""
    first_wall(eqt::IMAS.equilibrium__time_slice, pf_active::IMAS.pf_active{T}) where {T<:Real}

Returns named tuple with outline of the first wall, defined as the equilibrium computational domain with cutouts for the pf_active.coils that fall in it
"""
function first_wall(eqt::IMAS.equilibrium__time_slice, pf_active::IMAS.pf_active{T}) where {T<:Real}
    r_eq, z_eq = first_wall(eqt)
    eq_domain = collect(zip(r_eq, z_eq))

    R0 = sum(extrema(r_eq)) / 2.0
    Z0 = sum(extrema(z_eq)) / 2.0

    # pick point closes to this center
    rce = T[]
    zce = T[]
    for coil in pf_active.coil
        re = T[]
        ze = T[]
        for element in coil.element
            r, z = outline(element)
            append!(re, r)
            append!(ze, z)
        end
        index = argmin(sqrt.((re .- R0) .^ 2 .+ (ze .- Z0) .^ 2))
        if PolygonOps.inpolygon((re[index], ze[index]), eq_domain) == 1
            push!(rce, re[index])
            push!(zce, ze[index])
        end
    end

    r_eq, z_eq = resample_2d_path(r_eq, z_eq; method=:linear, n_points=100)
    append!(rce, r_eq)
    append!(zce, z_eq)

    # reorder
    index = sortperm(atan.(zce .- Z0, rce .- R0))

    return (r=rce[[index; index[1]]], z=zce[[index; index[1]]])
end

@compat public first_wall
push!(document[Symbol("Physics build")], :first_wall)

"""
    first_wall!(wall::IMAS.wall{T}, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}

Set `wall.description_2d[?].limiter.unit[1].outline` from input `r` and `z`

Returns the limiter.unit with the new outline
"""
function first_wall!(wall::IMAS.wall{T}, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
    oute = open_polygon(r, z)
    d2d = resize!(wall.description_2d, "limiter.type.index" => 0) # there should be only one limiter with type.index=0
    unit = resize!(d2d.limiter.unit, 1)[1]
    unit.closed = 1
    unit.name = "Contiguous first wall"
    unit.outline.r = oute.r
    unit.outline.z = oute.z
    return unit
end

@compat public first_wall!
push!(document[Symbol("Physics build")], :first_wall!)

"""
    contiguous_limiter_from_open_limiters!(wall::IMAS.wall{T}) where {T<:Real}

Combines multiple open limiter units into a single contiguous closed limiter by connecting 
their outlines in sequence. The function chooses the connection order to minimize gaps 
between adjacent limiter segments.

Modifies the wall structure in-place by replacing existing limiters with one closed 
contiguous limiter named "contiguous first wall" and limiter.type.index=0

Returns the limiter.unit with the contiguous outline.
"""
function contiguous_limiter_from_open_limiters!(wall::IMAS.wall{T}) where {T<:Real}
    r = T[]
    z = T[]
    for d2d in wall.description_2d
        for unit in d2d.limiter.unit
            if unit.closed == 0
                if isempty(r) || sqrt((r[end] - unit.outline.r[1])^2 + (z[end] - unit.outline.z[1])^2) < sqrt((r[end] - unit.outline.r[end])^2 + (z[end] - unit.outline.z[end])^2)
                    append!(r, unit.outline.r)
                    append!(z, unit.outline.z)
                else
                    append!(r, reverse(unit.outline.r))
                    append!(z, reverse(unit.outline.z))
                end
            end
        end
    end
    return first_wall!(wall, r, z)
end

@compat public contiguous_limiter_from_open_limiters!
push!(document[Symbol("Physics build")], :contiguous_limiter_from_open_limiters!)

"""
    build_max_R0_B0(bd::IMAS.build)

Returns the plasma geometric center (r0) and the maximum vacuum toroidal magnetic field (b0) evaluated at (r0) that the TF build allows
"""
function build_max_R0_B0(bd::IMAS.build)
    TFhfs = get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    plasma = get_build_layer(bd.layer; type=_plasma_)
    r0 = (plasma.start_radius + plasma.end_radius) / 2.0
    b0 = bd.tf.max_b_field / TFhfs.end_radius * r0
    return (r0=r0, b0=b0)
end

@compat public build_max_R0_B0
push!(document[Symbol("Physics build")], :build_max_R0_B0)

"""
    vertical_maintenance(bd::IMAS.build; tor_modularity::Int=2, pol_modularity::Int=1, gap_VV_BL::Float64=0.1)

Returns the radial dimensions of the vertical vacuum port for blanket maintenance

gap_VV_BL is the margin between the blanket module and the wall (10cm seems reasonable)
"""
function vertical_maintenance(bd::IMAS.build; tor_modularity::Int=2, pol_modularity::Int=1, gap_VV_BL::Float64=0.1)
    TFhfs = get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    VV = get_build_layers(bd.layer; type=_vessel_, fs=_lfs_)[1]
    BLhfs = get_build_layer(bd.layer; type=_blanket_, fs=_hfs_)
    BLlfs = get_build_layer(bd.layer; type=_blanket_, fs=_lfs_)
    iVVhfs = get_build_indexes(bd.layer; type=_vessel_, fs=_hfs_)[end]
    gap_VV_TF = bd.layer[iVVhfs-1]

    n_TF = bd.tf.coils_n
    n_modules = n_TF * tor_modularity
    phi_TF = 2π / n_TF
    phi_module = 2π / n_modules

    tor_thickness_TF = 2 * TFhfs.end_radius * tan(phi_TF / 2)
    if pol_modularity == 1
        rBL_ib = BLhfs.start_radius
        rBL_ob = BLlfs.end_radius
    elseif pol_modularity == 2
        rBL_ib = (BLlfs.end_radius + BLhfs.start_radius) / 2
        rBL_ob = BLlfs.end_radius
    else
        error("pol_modularity is an unallowed value - must be either 1 or 2")
    end

    w = tor_thickness_TF / 2 + gap_VV_TF.thickness + VV.thickness + gap_VV_BL
    r = (w - rBL_ib * sin(0.5 * (phi_TF - phi_module))) / sin(phi_TF / 2)

    rVP_hfs_ob = rBL_ib + r - gap_VV_BL
    rVP_hfs_ib = rVP_hfs_ob - VV.thickness
    rVP_lfs_ib = rVP_hfs_ob + 2 * gap_VV_BL + (rBL_ob - rBL_ib)
    rVP_lfs_ob = rVP_lfs_ib + VV.thickness

    return (rVP_hfs_ib=rVP_hfs_ib, rVP_hfs_ob=rVP_hfs_ob, rVP_lfs_ib=rVP_lfs_ib, rVP_lfs_ob=rVP_lfs_ob)
end

@compat public vertical_maintenance
push!(document[Symbol("Physics build")], :vertical_maintenance)

"""
    outline(layer::Union{IMAS.build__layer, IMAS.build__structure})

    outline(out::Union{IMAS.build__layer___outline,IMAS.build__structure___outline})

Returns a closed polygon as a named tuple with (r,z) of a `dd.build.layer` or `dd.build.structure`
"""
function outline(layer::Union{IMAS.build__layer,IMAS.build__structure})
    return outline(layer.outline)
end

function outline(out::Union{IMAS.build__layer___outline,IMAS.build__structure___outline})
    tmp = closed_polygon(out.r, out.z)
    return (r=tmp.R, z=tmp.Z)
end

@compat public outline
push!(document[Symbol("Physics build")], :outline)
