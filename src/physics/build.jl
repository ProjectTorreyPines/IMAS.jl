"""
    build_radii(bd::IMAS.build)

Return list of radii in the build
"""
function build_radii(bd::IMAS.build)
    layers_radii = typeof(bd.layer[1].thickness)[]
    layer_start = 0.0
    for l in bd.layer
        push!(layers_radii, layer_start)
        layer_start = layer_start + l.thickness
    end
    push!(layers_radii, layer_start)
    return layers_radii
end

function get_build(bd::IMAS.build; kw...)
    get_build(bd.layer; kw...)
end

"""
    get_build(
        layers::IMAS.IDSvector{<:IMAS.build__layer};
        type::Union{Nothing,BuildLayerType}=nothing,
        name::Union{Nothing,String}=nothing,
        identifier::Union{Nothing,Integer}=nothing,
        fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing,
        return_only_one::Bool=true,
        return_index::Bool=false,
        raise_error_on_missing::Bool=true)

Select layer(s) in build based on a series of selection criteria
With `raise_error_on_missing=false` will returns `missing` if layer is missing
"""
function get_build(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing,
    return_only_one::Bool=true,
    return_index::Bool=false,
    raise_error_on_missing::Bool=true)
    error("`IMAS.get_build()` is obsolete. Use `get_build_layer()` and `get_build_index()` or `get_build_layers()` and `get_build_indexes()` instead")
end

function get_build_layers(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

    valid_layers = eltype(layers)[]
    for (k, l) in enumerate(layers)
        if (name === nothing || l.name == name) &&
           (type === nothing || l.type == Int(type)) &&
           (identifier === nothing || l.identifier == identifier) &&
           (fs === nothing || (typeof(fs) <: AbstractVector{BuildLayerSide} && l.fs in map(Int, fs)) || (typeof(fs) <: BuildLayerSide && l.fs == Int(fs)))
            push!(valid_layers, l)
        end
    end

    return valid_layers
end

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
           (fs === nothing || (typeof(fs) <: AbstractVector{BuildLayerSide} && l.fs in map(Int, fs)) || (typeof(fs) <: BuildLayerSide && l.fs == Int(fs)))
            push!(valid_layers_indexes, k)
        end
    end

    return valid_layers_indexes
end

function get_build_layer(
    layers::IMAS.IDSvector{<:IMAS.build__layer};
    type::Union{Nothing,BuildLayerType}=nothing,
    name::Union{Nothing,String}=nothing,
    identifier::Union{Nothing,Integer}=nothing,
    fs::Union{Nothing,BuildLayerSide,AbstractVector{BuildLayerSide}}=nothing)

    valid_layers = get_build_layers(layers; type, name, identifier, fs)

    if isempty(valid_layers)
        error("Did not find build.layer: name=$(repr(name)) type=$type identifier=$identifier fs=$fs")
    elseif length(valid_layers) > 1
        error("Found multiple layers that satisfy name:$name type:$type identifier:$identifier fs:$fs")
    end

    return valid_layers[1]
end

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
    structures_mask(bd::IMAS.build; ngrid::Int = 257, border_fraction::Real = 0.1, one_is_for_vacuum::Bool = false)

return rmask, zmask, mask of structures that are not vacuum
"""
function structures_mask(bd::IMAS.build; ngrid::Int=257, border_fraction::Real=0.1, one_is_for_vacuum::Bool=false)
    border = maximum(bd.layer[end].outline.r) * border_fraction
    xlim = [0.0, maximum(bd.layer[end].outline.r) + border]
    ylim = [minimum(bd.layer[end].outline.z) - border, maximum(bd.layer[end].outline.z) + border]
    rmask = range(xlim[1], xlim[2], length=ngrid)
    zmask = range(ylim[1], ylim[2], length=ngrid * Int(round((ylim[2] - ylim[1]) / (xlim[2] - xlim[1]))))
    mask = ones(length(rmask), length(zmask))

    # start from the first vacuum that goes to zero outside of the TF
    start_from = -1
    for k in get_build_indexes(bd.layer, fs=_out_)
        if bd.layer[k].material == "Vacuum" && minimum(bd.layer[k].outline.r) < bd.layer[1].end_radius
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
            if (layer.material == "Vacuum") && (layer.fs != Int(_in_))
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
    rlim_oh = get_build_layer(bd.layer, type=_oh_).start_radius
    for (kr, rr) in enumerate(rmask)
        for (kz, zz) in enumerate(zmask)
            if rr < rlim_oh
                mask[kr, kz] = 1.0
            end
        end
    end
    if one_is_for_vacuum
        return rmask, zmask, 1.0 .- mask
    else
        return rmask, zmask, mask
    end
end

function func_nested_layers(layer::IMAS.build__layer, func::Function)
    i = index(layer)
    layers = parent(layer)
    # _in_ layers or plasma
    if layer.fs ∈ (Int(_in_), Int(_lhfs_))
        return func(layer)
        # anular layers
    elseif layer.fs ∈ (Int(_hfs_), Int(_lfs_))
        i = index(layer)
        if layer.fs == Int(_hfs_)
            layer_in = layers[i+1]
        else
            layer_in = layers[i-1]
        end
        return func(layer) - func(layer_in)
        # _out_ layers
    elseif layer.fs == Int(_out_)
        layer_in = layers[i-1]
        if layer_in.fs == Int(_out_)
            return func(layer) - func(layer_in)
        elseif layer_in.fs == Int(_lfs_)
            func_in = 0.0
            for l in layers
                if l.fs == Int(_in_)
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
    func_nested_layers(layer, l -> area(l.outline.r, l.outline.z))
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
    tf_ripple(r, R_tf::Real, N_tf::Integer)

Evaluate fraction of toroidal magnetic field ripple at `r` [m]
generated from `N_tf` toroidal field coils with outer leg at `R_tf` [m]
"""
function tf_ripple(r, R_tf::Real, N_tf::Integer)
    eta = (r ./ R_tf) .^ N_tf
    return eta ./ (1.0 .- eta)
end

"""
    R_tf_ripple(r, ripple::Real, N_tf::Integer)

Evaluate location of toroidal field coils outer leg `R_tf`` [m] at which `N_tf`
toroidal field coils generate a given fraction of toroidal magnetic field ripple at `r` [m]
"""
function R_tf_ripple(r, ripple::Real, N_tf::Integer)
    return r .* (ripple ./ (ripple .+ 1.0)) .^ (-1 / N_tf)
end

"""
    first_wall(wall::IMAS.wall)

return outline of first wall
"""
function first_wall(wall::IMAS.wall)
    if (!ismissing(wall.description_2d, ["1", "limiter", "unit", "1", "outline", "r"])) && (length(wall.description_2d[1].limiter.unit[1].outline.r) > 5)
        return wall.description_2d[1].limiter.unit[1].outline
    else
        return missing
    end
end
