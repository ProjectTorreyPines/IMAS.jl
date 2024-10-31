"""
    build_radii(bd::IMAS.build)

Return list of radii in the build
"""
function build_radii(bd::IMAS.build)
    return _build_radii(bd.layer)
end

function _build_radii(layer)
    layers_radii = zeros(typeof(layer[1].thickness), length(layer) + 1)
    for (k, l) in enumerate(layer)
        @inbounds layers_radii[k+1] = layers_radii[k] + l.thickness
    end
    return layers_radii
end

function get_build(bd::IMAS.build; kw...)
    return get_build(bd.layer; kw...)
end

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

"""
    opposite_side_layer(layer::IMAS.build__layer)

Returns corresponding layer on the high/low field side 
"""
function opposite_side_layer(layer::IMAS.build__layer)
    @assert layer.side in (Int(_hfs_), Int(_lfs_)) "Asking for opposite_side_layer of $(layer.name), which does not make sense"
    return get_build_layer(parent(layer); identifier=layer.identifier, fs=(layer.side == Int(_lfs_)) ? _hfs_ : _lfs_)
end

"""
    structures_mask(bd::IMAS.build; ngrid::Int = 257, border_fraction::Real = 0.1)

return rmask, zmask, mask of structures that are not vacuum
"""
function structures_mask(bd::IMAS.build; ngrid::Int=257, border_fraction::Real=0.1, layer_check::Function=layer -> layer.material == "vacuum")
    border = maximum(bd.layer[end].outline.r) * border_fraction
    xlim = [0.0, maximum(bd.layer[end].outline.r) + border]
    ylim = [minimum(bd.layer[end].outline.z) - border, maximum(bd.layer[end].outline.z) + border]
    rmask = range(xlim[1], xlim[2], ngrid)
    zmask = range(ylim[1], ylim[2], ngrid * Int(round((ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))
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

"""
    area(layer::IMAS.build__layer)

Calculate area of a build layer outline
"""
function area(layer::IMAS.build__layer)
    return func_nested_layers(layer, l -> area(l.outline.r, l.outline.z))
end

"""
    volume(layer::IMAS.build__layer)

Calculate volume of a build layer outline revolved around z axis
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

"""
    first_wall(wall::IMAS.wall)

returns named tuple with outline of first wall, or an empty outline if not present
"""
function first_wall(wall::IMAS.wall{T}) where {T<:Real}
    if (!ismissing(wall.description_2d, ["1", "limiter", "unit", "1", "outline", "r"])) && (length(wall.description_2d[1].limiter.unit[1].outline.r) > 4)
        outline = wall.description_2d[1].limiter.unit[1].outline
        tmp = closed_polygon(outline.r, outline.z)
        return (r=tmp.r, z=tmp.z)
    else
        return (r=Float64[], z=Float64[])
    end
end

"""
    build_max_R0_B0(bd::IMAS.build)

Returns the plasma geometric center and the maximum vacuum toroidal magnetic field at the plasma geometric center that the TF build allows
"""
function build_max_R0_B0(bd::IMAS.build)
    TFhfs = get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    plasma = get_build_layer(bd.layer; type=_plasma_)
    R0 = (plasma.start_radius + plasma.end_radius) / 2.0
    B0 = bd.tf.max_b_field / TFhfs.end_radius * R0
    return R0, B0
end

"""
    vertical_maintenance(bd::IMAS.build; tor_modularity::Int=2, pol_modularity::Int=1, gap_VV_BL::Float64=0.1)

Returns the radial dimensions of the vertical vacuum port for blanket maintenance

gap_VV_BL is the margin between the blanket module and the wall (10cm seems reasonable)
"""
function vertical_maintenance(bd::IMAS.build; tor_modularity::Int=2, pol_modularity::Int=1, gap_VV_BL::Float64=0.1)
    TFhfs = get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    VV = get_build_layer(bd.layer; type=_vessel_, fs=_lfs_)
    BLhfs = get_build_layer(bd.layer; type=_blanket_, fs=_hfs_)
    BLlfs = get_build_layer(bd.layer; type=_blanket_, fs=_lfs_)
    iVVhfs = get_build_index(bd.layer; type=_vessel_, fs=_hfs_)
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

"""
outline(layer::Union{IMAS.build__layer, IMAS.build__structure})

Returns outline as named tuple with (r,z)

NOTE: returns a polygon that always closes
"""
function outline(layer::Union{IMAS.build__layer,IMAS.build__structure})
    return outline(layer.outline)
end

function outline(out::Union{IMAS.build__layer___outline,IMAS.build__structure___outline})
    tmp = closed_polygon(out.r, out.z)
    return (r=tmp.R, z=tmp.Z)
end
