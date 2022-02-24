import JSON

include("dd.jl")

include("utils.jl")

include("f2.jl")

"""
    dict2imas(dct, ids::IDS=dd() ;verbose::Bool=false, path::Vector{String}=String[])::IDS

Populate IMAS data structure `ids` based on data contained in Julia dictionary `dct`.

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
- `skip_non_coordinates::Bool=false`: only assign coordinates to the data structure
"""
function dict2imas(dct, ids::T; verbose::Bool = false, path::Vector{String} = String[], skip_non_coordinates::Bool = false) where {T<:IDS}
    # recursively traverse `dtc` structure
    level = length(path)
    for (k, v) in dct
        if !hasfield(typeof(ids), Symbol(k))
            if !skip_non_coordinates
                @warn("$(f2i(ids)).$(k) was skipped in IMAS.jl data dictionary", maxlog = 1)
            end
            continue
        end
        # Struct
        if typeof(v) <: Dict
            if verbose
                println(("｜"^level) * string(k))
            end
            ff = getfield(ids, Symbol(k))
            dict2imas(v, ff; path = vcat(path, [string(k)]), verbose = verbose, skip_non_coordinates = skip_non_coordinates)

            # Array of struct
        elseif (typeof(v) <: Array) && (length(v) > 0) && (typeof(v[1]) <: Dict)
            ff = getfield(ids, Symbol(k))
            if verbose
                println(("｜"^level) * string(k))
            end
            if length(ff) < length(v)
                resize!(ff, length(v))
            end
            for i = 1:length(v)
                if verbose
                    println(("｜"^(level + 1)) * string(i))
                end
                dict2imas(v[i], ff[i]; path = vcat(path, [string(k), "[$i]"]), verbose = verbose, skip_non_coordinates = skip_non_coordinates)
            end

            # Leaf
        else
            if verbose
                print(("｜"^level) * string(k) * " → ")
            end
            target_type = Core.Compiler.typesubtract(struct_field_type(typeof(ids), Symbol(k)), Union{Missing,Function}, 1)
            if target_type <: AbstractArray
                if tp_ndims(target_type) == 2
                    v = transpose(reduce(hcat, v))
                end
                if tp_eltype(target_type) <: Real
                    v = convert(Array{Float64,tp_ndims(target_type)}, v)
                else
                    v = convert(Array{tp_eltype(target_type),tp_ndims(target_type)}, v)
                end
            end
            setproperty!(ids, Symbol(k), v; skip_non_coordinates = skip_non_coordinates)
            if verbose
                println(typeof(v))
            end
        end
    end

    return ids
end

tp_ndims(::Type{AbstractArray{T,N} where {T,N}}) = N
tp_ndims(v::UnionAll) = ndims(v.body)
tp_eltype(::Type{AbstractArray{T,N} where {T,N}}) = T
tp_eltype(v::UnionAll) = v.var.ub

"""
    json2imas(filename::String; verbose::Bool=false)::IDS

Load from a file with give `filename` the IMAS data structure saved in JSON format 

# Arguments
- `verbose::Bool=false`: print structure hierarchy as it is filled
"""
function json2imas(filename::String; verbose::Bool = false)::IDS
    ids_data = dd()
    json_data = JSON.parsefile(filename)
    dict2imas(json_data, ids_data; verbose = verbose, skip_non_coordinates = true)
    dict2imas(json_data, ids_data; verbose = verbose, skip_non_coordinates = false)
    return ids_data
end

"""
    imas2dict(ids::Union{IDS,IDSvector})

Populate Julia structure of dictionaries and vectors with data from IMAS data structure `ids`
"""
function imas2dict(ids::Union{IDS,IDSvector})
    if typeof(ids) <: IDSvector
        dct = Any[]
    else
        dct = Dict()
    end
    return imas2dict(ids, dct)
end

function imas2dict(ids::Union{IDS,IDSvector}, dct::Union{Dict,Vector})
    items = sort(keys(ids))
    for item in items
        if typeof(ids) <: IDSvector # arrays of structures
            push!(dct, Dict())
            imas2dict(ids[item], dct[item])
        else
            value = missing
            try
                value = getproperty(ids, item)
            catch e
                if typeof(e) <: IMASexpressionError
                    #@debug(sprint(showerror, e))
                    continue
                else
                    throw
                end
            end
            if typeof(value) <: Union{IDS,IDSvector} # structures
                if typeof(value) <: IDS
                    dct[item] = Dict()
                else
                    dct[item] = Any[]
                end
                imas2dict(value, dct[item])
            else # field
                dct[item] = value
            end
        end
    end
    return dct
end

"""
    imas2json(ids::Union{IDS,IDSvector}, filename::String; kw...)

Save the IMAS data structure to a JSON file with give `filename`
"""
function imas2json(ids::Union{IDS,IDSvector}, filename::String; kw...)
    open(filename, "w") do io
        JSON.print(io, imas2dict(ids); kw...)
    end
end


"""
    top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return top-level IDS in the DD hierarchy.
Considers IDS as maximum top level if IDS_is_absolute_top=true
"""
function top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool = true)
    if IDS_is_absolute_top & (typeof(ids) <: dd)
        error("Cannot call top(x::IMAS.dd,IDS_is_absolute_top=true). Use `IDS_is_absolute_top=false`.")
    elseif !(typeof(ids._parent.value) <: Union{IDS,IDSvector})
        return ids
    elseif IDS_is_absolute_top & (typeof(ids._parent.value) <: dd)
        return ids
    else
        return top(ids._parent.value; IDS_is_absolute_top = IDS_is_absolute_top)
    end
end

"""
    top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return top-level IDS in the DD hierarchy.
Considers IDS as maximum top level if IDS_is_absolute_top=true
"""

function top_ids(ids::Union{IDS,IDSvector})
    ids = top(ids::Union{IDS,IDSvector}; IDS_is_absolute_top = true)
    if length(f2p(ids)) == 1
        return ids
    else
        return missing
    end
end

function top_dd(ids::Union{IDS,IDSvector})
    ids = top(ids; IDS_is_absolute_top = false)
    if typeof(ids) <: dd
        return ids
    else
        return missing
    end
end

"""
    parent(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool=true)

Return parent IDS/IDSvector in the hierarchy
If IDS_is_absolute_top then returns `missing` instead of IMAS.dd()
"""
function parent(ids::Union{IDS,IDSvector}; IDS_is_absolute_top::Bool = true)
    if !(typeof(ids._parent.value) <: Union{IDS,IDSvector})
        return missing
    elseif IDS_is_absolute_top & (typeof(ids._parent.value) <: dd)
        return missing
    else
        return ids._parent.value
    end
end

"""
    children(ids::Union{IDS,IDSvector})::Vector{Symbol}

Return children of a IDS/IDSvector
"""
function children(ids::Union{IDS,IDSvector})::Vector{Symbol}
    return [k for k in fieldnames(typeof(ids)) if k != :_parent]
end


"""
    assign_expressions(ids::Union{IDS,IDSvector})

Assign expressions to a IDS/IDSvector
NOTE: This is done not recursively
"""
function assign_expressions(ids::Union{IDS,IDSvector})
    struct_name = f2u(ids)
    for item in children(ids)
        if typeof(getfield(ids, item)) <: Union{IDS,IDSvector}
            continue
        elseif "$(struct_name).$(item)" in keys(expressions)
            setproperty!(ids, item, expressions["$(struct_name).$(item)"])
        end
    end
    return ids
end

#= ===================== =#
#  IDS related functions  #
#= ===================== =#

"""
    goto(ids::IDS, location::String)

Reach location in a given IDS

# Arguments
- `f2::Function=f2i`: function used to process the IDS path to be compared to `location`
"""
function goto(ids::IDS, location::String; f2::Function = f2i)
    # find common ancestor
    cs, s1, s2 = IMAS.common_base_string(f2(ids), location)
    cs0 = replace(cs, r"\.$" => "")
    # go upstream until common acestor
    h = ids
    while f2(h) != cs0
        h = h._parent.value
        if h === missing
            error("Could not reach `$(location)` from `$(f2(ids))`")
        end
    end
    # then dive into the location branch
    for k in i2p(s2)
        h = getfield(h, Symbol(k))
    end
    return h
end

"""
    coordinates(ids::IDS, field::Symbol)

Return two lists, one of coordinate names and the other with their values in the data structure
Coordinate value is `nothing` when the data does not have a coordinate
Coordinate value is `missing` if the coordinate is missing in the data structure
"""
function coordinates(ids::IDS, field::Symbol)
    info = imas_info(ids, field)
    # handle scalar quantities (which do not have coordinate)
    if !("coordinates" in keys(info))
        return Dict(:names => [], :values => [])
    end
    coord_names = deepcopy(info["coordinates"])
    coord_values = []
    for coord in coord_names
        if occursin("...", coord)
            push!(coord_values, nothing)
        else
            h = goto(ids, u2fs(p2i(i2p(coord)[1:end-1])); f2 = f2fs)
            coord_leaf = Symbol(i2p(coord)[end])
            h = getfield(h, coord_leaf)
            # add value to the coord_values
            push!(coord_values, h)
        end
    end
    return Dict(:names => coord_names, :values => coord_values)
end

"""
    units(ids::IDS, field::Symbol)

Return string with units
"""
function units(location::String)
    info = imas_info(location)
    return get(info, "units", "")
end

"""
    units(ids::IDS, field::Symbol)

Return string with units
"""
function units(ids::IDS, field::Symbol)
    location = "$(f2u(ids)).$(field)"
    return units(location)
end

"""
    coords ids.path.to.array.y => interpolating_x

Macro for interpolating data
"""
macro coords(ex)
    return _coords(ex)
end

function _coords(ex)
    quote
        local expr = $(esc(Meta.QuoteNode(ex)))
        if (expr.head !== :call) || (expr.args[1] != :(=>))
            error("@coords must use `dd.ids.field => xx` syntax")
        end
        local ids = $(esc(ex.args[2].args[1]))
        local field = $(esc(ex.args[2].args[2]))
        local xx = $(esc(ex.args[3]))
        local yy = interp(ids, field)(xx)
        yy
    end
end

"""
    Base.keys(ids::IDS)

Returns list of fields with data in a IDS
"""
function Base.keys(ids::IDS)
    kkk = Symbol[]
    for k in fieldnames(typeof(ids))
        # hide the _parent field
        if k === :_parent
            continue
        end
        v = getfield(ids, k)
        # empty entries
        if v === missing
            continue
            # empty structures/arrays of structures (recursive)
        elseif typeof(v) <: Union{IDS,IDSvector}
            if length(keys(v)) > 0
                push!(kkk, k)
            end
            # entries with data
        else
            push!(kkk, k)
        end
    end
    return kkk
end

function Base.keys(ids::IDSvector)
    return collect(1:length(ids))
end

function Base.show(io::IO, ids::Union{IDS,IDSvector}, depth::Int)
    items = keys(ids)
    for (k, item) in enumerate(sort(items))
        printstyled(io, "$('｜'^depth)"; color = :yellow)
        # arrays of structurs
        if typeof(ids) <: IDSvector
            printstyled(io, "[$(item)]\n"; bold = true, color = :green)
            show(io, ids[item], depth + 1)
        else
            value = getfield(ids, item)
            if typeof(value) <: Union{IDS,IDSvector}
                # structures
                if (typeof(ids) <: dd)
                    printstyled(io, "$(uppercase(string(item)))\n"; bold = true)
                else
                    printstyled(io, "$(string(item))\n"; bold = true)
                end
                show(io, value, depth + 1)
            else
                # field
                printstyled(io, "$(item)")
                printstyled(io, " ➡ "; color = :red)
                if typeof(value) <: Function
                    printstyled(io, "Function\n"; color = :blue)
                elseif typeof(value) <: String
                    printstyled(io, "\"$(value)\"\n"; color = :magenta)
                elseif typeof(value) <: Integer
                    printstyled(io, "$(value)\n"; color = :yellow)
                elseif typeof(value) <: AbstractFloat
                    printstyled(io, @sprintf("%g\n", value); color = :red)
                elseif typeof(value) <: AbstractArray
                    if length(value) < 5
                        if eltype(value) <: AbstractFloat
                            printstyled(io, "[$(join([@sprintf("%g",v) for v in value],","))]\n"; color = :green)
                        else
                            printstyled(io, "$value\n"; color = :green)
                        end
                    else
                        printstyled(io, "$(Base.summary(value))\n"; color = :green)
                    end
                else
                    printstyled(io, "$(Base.summary(value))\n"; color = :blue)
                end
            end
        end
        if (typeof(ids) <: dd) & (k < length(items))
            println(io, "")
        end
    end
end

# show function for the Jupyter notebook
function Base.show(io::IO, ::MIME"text/plain", ids::Union{IDS,IDSvector})
    return show(io, ids, 0)
end

# show function for inline prints
function Base.show(io::IO, ids::IDS)
    fnames = []
    for item in keys(ids)
        push!(fnames, item)
    end
    return println(io, "$(f2i(ids)){$(join(collect(map(x->"$x",fnames)),", "))}")
end

# show function for inline prints
function Base.show(io::IO, ids::IDSvector)
    fnames = []
    for item in keys(ids)
        push!(fnames, item)
    end
    if length(ids) < 2
        return println(io, "$(p2i(f2p(ids)[1:end-1]))[$(join(collect(map(x->"$x",fnames)),", "))]")
    else
        return println(io, "$(p2i(f2p(ids)[1:end-1]))[1...$(length(ids))]")
    end
end