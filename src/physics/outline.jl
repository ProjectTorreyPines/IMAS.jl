#NOTE: outlines in IMAS ontology are always open

struct OutlineClosedVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    data_is_closed::Bool
    """
        OutlineClosedVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

    Outlines in IMAS ontology are always open.

    Use OutlineClosedVector to have a vector that behaves like a it's a closed outline

    When OutlineClosedVector is put into the dd, it will automatically save it there as an open outline
    """
    function OutlineClosedVector(data::A, data_is_closed::Bool) where {T,A<:AbstractVector{T}}
        if data_is_closed
            @assert data[1] == data[end]
        end
        return new{T,A}(data, data_is_closed)
    end
end

# Define the length method, which adds 1 if the polygon is not closed
Base.size(o::OutlineClosedVector) = (length(o.data) + (o.data_is_closed ? 0 : 1),)

# Define indexing, where the last point is the first point if the polygon is not closed
function Base.getindex(o::OutlineClosedVector, i::Int)
    if i < 1 || i > length(o)
        throw(BoundsError(o, i))
    end
    if i == length(o) && !o.data_is_closed
        # Return the first point if this is the virtual closing point
        return o.data[1]
    else
        return o.data[i]
    end
end

function Base.setindex!(o::OutlineClosedVector, v::Any, i::Int)
    if i < 1 || i > length(o)
        throw(BoundsError(o, i))
    end
    if i == length(o) && !o.data_is_closed
        # Return the first point if this is the virtual closing point
        return o.data[1] = v
    else
        return o.data[i] = v
    end
end

# always save the OutlineClosedVector to the dd as a the vector of a open polygon
function Base.setproperty!(@nospecialize(ids::IMASdd.IDS), field::Symbol, o::IMAS.OutlineClosedVector)
    if o.data_is_closed
        return setproperty!(ids, field, o.data[1:end-1])
    else
        return setproperty!(ids, field, o.data)
    end
end

# ==================

struct OutlineOpenVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    data_is_open::Bool
    """
        OutlineOpenVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

    Outlines in IMAS ontology are always open.

    Use OutlineOpenVector to have a vector that behaves like a it's a open outline

    When OutlineOpenVector is put into the dd, it will automatically save it there as an open outline
    """
    function OutlineOpenVector(data::A, data_is_open::Bool) where {T,A<:AbstractVector{T}}
        if !data_is_open
            @assert data[1] == data[end]
        end
        return new{T,A}(data, data_is_open)
    end
end

# Define the length method, which subtracts 1 if the polygon is not open
Base.size(o::OutlineOpenVector) = (length(o.data) + (o.data_is_open ? 0 : -1),)

# Define indexing, where the last point is omitted if the polygon is not open
function Base.getindex(o::OutlineOpenVector, i::Int)
    if i < 1 || i > length(o)
        throw(BoundsError(o, i))
    end
    return o.data[i]
end

function Base.setindex!(o::OutlineOpenVector, v::Any, i::Int)
    if i < 1 || i > length(o)
        throw(BoundsError(o, i))
    end
    o.data[i] = v
    if i == 1
        o.data[end] = v
    end
    return v
end

# always save the OutlineOpenVector to the dd as a the vector of a open polygon
function Base.setproperty!(@nospecialize(ids::IMASdd.IDS), field::Symbol, o::IMAS.OutlineOpenVector)
    if o.data_is_open
        return setproperty!(ids, field, o.data)
    else
        return setproperty!(ids, field, o.data[1:end-1])
    end
end

# ===================

struct CircularVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    offset::Int
    """
        CircularVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

    Vector with circular indexes
    """
    function CircularVector(data::A, offset::Int=0) where {T,A<:AbstractVector{T}}
        return new{T,A}(data, offset)
    end
end

function Base.checkbounds(o::CircularVector, I...)
    nothing
end

# Define the length method, which subtracts 1 if the polygon is not open
Base.size(o::CircularVector) = size(o.data)

# Define indexing, where the last point is omitted if the polygon is not open
function Base.getindex(o::CircularVector, i::Int)
    i = i + o.offset
    i = mod(i - 1, length(o)) + 1
    return o.data[i]
end

function Base.setindex!(o::CircularVector, v::Any, i::Int)
    i = i + o.offset
    i = mod(i - 1, length(o)) + 1
    o.data[i] = v
    return v
end

# always save the CircularVector to the dd as a the vector of a open polygon
function Base.setproperty!(@nospecialize(ids::IMASdd.IDS), field::Symbol, o::IMAS.CircularVector)
    v = setproperty!(ids, field, o.data)
    circshift!(v, o.offset)
end