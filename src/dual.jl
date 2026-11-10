document[:Dual] = Symbol[]

# ======================================================== #
# Data dictionaries whose element type is ForwardDiff.Dual #
# ======================================================== #
# A `dd{D}` with `D<:ForwardDiff.Dual` carries, for every real-valued leaf, the derivatives with
# respect to a small set of seeded inputs. `time` arrays, `global_time`, integers and strings stay
# as they are. The functions here convert between element types, seed, and split value and
# derivative parts; they are the only place where element types of a whole IDS are changed.

"""
    InputsTag

Tag of the dual numbers that carry derivatives with respect to the inputs of a whole simulation.
It is ordered before every other ForwardDiff tag, so that the derivative calculations that some
physics routines do internally (which create their own tags) nest inside it without ambiguity.
Use `InputsDual{N}` as the element type of a `dd` that carries `N` partial derivatives.
"""
struct InputsTag end
const InputsDual{N} = Dual{InputsTag,Float64,N}
ForwardDiff.:≺(::Type{InputsTag}, ::Type{<:ForwardDiff.Tag}) = true
ForwardDiff.:≺(::Type{<:ForwardDiff.Tag}, ::Type{InputsTag}) = false
ForwardDiff.:≺(::Type{InputsTag}, ::Type{InputsTag}) = false

@compat public InputsTag, InputsDual
push!(document[:Dual], :InputsTag)
push!(document[:Dual], :InputsDual)

"""
    retype(::Type{T}, ids::Union{IDS,IDSvector}) where {T<:Real}

Copy of `ids` with element type `T`. For a `dd` the `global_time` and the auxiliary storage are carried over.

  - `Float64 -> Dual`: values are kept, partials are zero
  - `Dual -> Float64`: values are kept, partials are dropped
  - `Dual -> Dual` of a different type: values are kept, partials are zero

Only filled fields are copied, so expressions remain expressions in the new element type.
"""
function retype(::Type{T}, @nospecialize(ids::IDS)) where {T<:Real}
    ids_new = Base.typename(typeof(ids)).wrapper{T}(; frozen=getfield(ids, :_frozen))
    fill!(ids_new, ids)
    if ids isa DD
        setfield!(ids_new, :global_time, getfield(ids, :global_time))
        aux_new = getfield(ids_new, :_aux)
        for (k, v) in getfield(ids, :_aux)
            aux_new[k] = (v isa Union{IDS,IDSvector}) ? retype(T, v) : deepcopy(v)
        end
    end
    return ids_new
end

function retype(::Type{T}, @nospecialize(idsv::IDSvector)) where {T<:Real}
    idsv_new = Base.typename(typeof(idsv)).wrapper{Base.typename(eltype(idsv)).wrapper{T}}()
    return fill!(idsv_new, idsv)
end

@compat public retype
push!(document[:Dual], :retype)

# `convert(D, ids)` is defined in IMASdd for any `Type{<:Real}` and in ForwardDiff for any argument,
# so for a Dual element type the call is ambiguous unless it is resolved here
Base.convert(::Type{D}, @nospecialize(ids::IDS)) where {D<:Dual} = retype(D, ids)
Base.convert(::Type{D}, @nospecialize(ids::IDSvector)) where {D<:Dual} = retype(D, ids)

"""
    value_dd(ids::Union{IDS,IDSvector})

Float64 copy of a Dual-typed IDS holding the value parts
"""
value_dd(@nospecialize(ids::Union{IDS{<:Dual},IDSvector{<:IDS{<:Dual}}})) = retype(Float64, ids)
value_dd(@nospecialize(ids::Union{IDS{Float64},IDSvector{<:IDS{Float64}}})) = deepcopy(ids)

@compat public value_dd
push!(document[:Dual], :value_dd)

# element-wise helpers that leave non-Dual data untouched
_dual_value(x::Dual) = ForwardDiff.value(x)
_dual_value(x) = x
_dual_value(x::AbstractArray{<:Dual}) = ForwardDiff.value.(x)
_dual_partial(x::Dual, k::Int) = ForwardDiff.partials(x, k)
_dual_partial(x::AbstractArray{<:Dual}, k::Int) = ForwardDiff.partials.(x, k)

"""
    fill!(ids_new::IDS{Float64}, ids::IDS{<:Dual}, field::Symbol)

Maps fields of an `IDS{Dual}` to an `IDS{Float64}`, keeping the value parts
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Float64,T2<:Dual}
    value = getraw(ids, field)
    if field != :time && eltype(value) <: T2
        _setproperty!(ids_new, field, _dual_value(value); from_cocos=internal_cocos)
    else
        _setproperty!(ids_new, field, deepcopy(value); from_cocos=internal_cocos)
    end
    return nothing
end

"""
    fill!(ids_new::IDS{<:Dual}, ids::IDS{<:Dual}, field::Symbol)

Maps fields between Dual-typed IDSs. Same type: plain copy. Different Dual types: values are kept and partials are zero.
"""
function Base.fill!(@nospecialize(ids_new::IDS{<:T1}), @nospecialize(ids::IDS{<:T2}), field::Symbol) where {T1<:Dual,T2<:Dual}
    value = getraw(ids, field)
    D1 = eltype(ids_new)
    if D1 === eltype(ids) || field == :time || !(eltype(value) <: T2)
        _setproperty!(ids_new, field, deepcopy(value); from_cocos=internal_cocos)
    else
        _setproperty!(ids_new, field, D1.(_dual_value(value)); from_cocos=internal_cocos)
    end
    return nothing
end

# visit every filled data leaf of `ids`, calling f(ids_node, other_node, field); `other` must have the same structure
function _foreach_leaf(f::F, @nospecialize(ids::Union{IDS,IDSvector}), @nospecialize(other::Union{IDS,IDSvector})) where {F}
    stack = Tuple{Any,Any}[]
    if ids isa IDSvector
        @assert length(ids) == length(other) "IDS structures do not match"
        for i in eachindex(ids)
            push!(stack, (ids[i], other[i]))
        end
    else
        push!(stack, (ids, other))
    end
    while !isempty(stack)
        a, b = pop!(stack)
        filled = getfield(a, :_filled)
        T = typeof(a)
        for i in 1:fieldcount(T)
            field = fieldname(T, i)
            if hasfield(typeof(filled), field) && getfield(filled, field)
                field_type = fieldtype_typeof(a, field)
                if field_type <: IDS
                    push!(stack, (getfield(a, field), getfield(b, field)))
                elseif field_type <: IDSvector
                    va = getfield(a, field)
                    vb = getfield(b, field)
                    @assert length(va) == length(vb) "IDS structures do not match at $(location(a, field))"
                    for j in eachindex(va)
                        push!(stack, (va[j], vb[j]))
                    end
                else
                    f(a, b, field)
                end
            end
        end
    end
    return ids
end

"""
    partials_dd(ids::Union{IDS,IDSvector}, k::Int)

Float64 copy of a Dual-typed IDS holding the k-th partial derivative of every real-valued leaf.
Leaves that are not Dual-typed (times, integers, strings) are copied unchanged.
"""
function partials_dd(@nospecialize(ids::Union{IDS{<:Dual},IDSvector{<:IDS{<:Dual}}}), k::Int)
    out = retype(Float64, ids)
    _foreach_leaf(ids, out) do a, b, field
        value = getraw(a, field)
        if field != :time && eltype(value) <: Dual
            setfield!(b, field, _dual_partial(value, k))
        end
    end
    return out
end

@compat public partials_dd
push!(document[:Dual], :partials_dd)

"""
    replace_values!(ids::Union{IDS,IDSvector}, ids_values::Union{IDS,IDSvector})

Overwrite the value part of every real-valued leaf of the Dual-typed `ids` with the corresponding
entry of the Float64 `ids_values`, keeping the partials of `ids`. Both must have the same structure
and array sizes. Used to iterate derivative parts at a fixed solution.
"""
function replace_values!(@nospecialize(ids::Union{IDS{<:Dual},IDSvector{<:IDS{<:Dual}}}), @nospecialize(ids_values::Union{IDS{Float64},IDSvector{<:IDS{Float64}}}))
    _foreach_leaf(ids, ids_values) do a, b, field
        value = getraw(a, field)
        if field != :time && eltype(value) <: Dual
            newval = getraw(b, field)
            if value isa AbstractArray
                @assert size(value) == size(newval) "array sizes do not match at $(location(a, field))"
                setfield!(a, field, map((d, v) -> typeof(d)(v, ForwardDiff.partials(d)), value, newval))
            else
                setfield!(a, field, typeof(value)(newval, ForwardDiff.partials(value)))
            end
        end
    end
    return ids
end

@compat public replace_values!
push!(document[:Dual], :replace_values!)

"""
    npartials(ids::Union{IDS,IDSvector})

Number of partial derivatives carried by the element type of `ids` (0 for a non-Dual IDS)
"""
npartials(@nospecialize(ids::Union{IDS,IDSvector})) = _npartials(ids isa IDSvector ? eltype(eltype(ids)) : eltype(ids))
_npartials(::Type{<:Dual{T,V,N}}) where {T,V,N} = N
_npartials(::Type) = 0

@compat public npartials
push!(document[:Dual], :npartials)

"""
    seed(::Type{D}, k::Int) where {D<:Dual}

Dual number of value 1 whose k-th partial is 1 and all others 0. Multiplying an input by `seed(D, k)`
makes it the k-th differentiation variable in relative (logarithmic) terms: the k-th partial of any
output `y` is then ∂y/∂ln(input).
"""
function seed(::Type{D}, k::Int) where {T,V,N,D<:Dual{T,V,N}}
    @assert 1 <= k <= N "partial index $k is outside 1:$N"
    return D(one(V), ForwardDiff.Partials{N,V}(ntuple(i -> i == k ? one(V) : zero(V), N)))
end

@compat public seed
push!(document[:Dual], :seed)

# ---- #
# I/O  #
# ---- #
# JSON and HDF5 files hold plain numbers: a Dual-typed IDS is written as one file with the values
# and one file per partial derivative, each of which reads back as a regular Float64 IDS.

"""
    dual_filenames(filename::AbstractString, n::Int)

Names of the files written for a Dual-typed IDS: `filename` for the values and `<base>_d<k><ext>` for the k-th partial
"""
function dual_filenames(filename::AbstractString, n::Int)
    base, ext = splitext(filename)
    return (values=filename, partials=["$(base)_d$(k)$(ext)" for k in 1:n])
end

function IMASdd.imas2json(@nospecialize(ids::Union{IDS{<:Dual},IDSvector{<:IDS{<:Dual}}}), filename::AbstractString; partials::Bool=true, kw...)
    names = dual_filenames(filename, npartials(ids))
    out = IMASdd.imas2json(retype(Float64, ids), names.values; kw...)
    if partials
        for (k, name) in enumerate(names.partials)
            IMASdd.imas2json(partials_dd(ids, k), name; kw...)
        end
    end
    return out
end

function IMASdd.imas2hdf(@nospecialize(ids::Union{IDS{<:Dual},IDSvector{<:IDS{<:Dual}}}), filename::AbstractString; partials::Bool=true, kw...)
    names = dual_filenames(filename, npartials(ids))
    out = IMASdd.imas2hdf(retype(Float64, ids), names.values; kw...)
    if partials
        for (k, name) in enumerate(names.partials)
            IMASdd.imas2hdf(partials_dd(ids, k), name; kw...)
        end
    end
    return out
end
