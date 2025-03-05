document[Symbol("Outlines")] = Symbol[]

"""
    OutlineClosedVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

Outlines in IMAS v3 ontology are always open (NOTE: in v4 they will always be closed!)

Use OutlineClosedVector to have a vector that behaves like a it's a closed outline, without creating a new array (allocation free)

When OutlineClosedVector is put into the `dd`, it will automatically save it there as an open outline
"""
struct OutlineClosedVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    data_is_closed::Bool
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

# Custom materialize! for OutlineClosedVector
function Base.materialize!(dest::OutlineClosedVector, bc::Base.Broadcast.Broadcasted)
    if dest.data_is_closed
        # Avoid reassigning the last point for closed outlines
        for i in 1:length(dest.data)-1
            dest.data[i] = bc[i]
        end
        dest.data[end] = dest.data[1]  # Ensures closure without reapplication
    else
        # Direct materialize! for open outlines
        for i in 1:length(dest.data)
            dest.data[i] = bc[i]
        end
    end
    return dest
end

@compat public OutlineClosedVector
push!(document[Symbol("Outlines")], :OutlineClosedVector)

# ==================

"""
    OutlineOpenVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

Outlines in IMAS v3 ontology are always open (NOTE: in v4 they will always be closed!)

Use OutlineOpenVector to have a vector that behaves like a it's a open outline, without creating a new array (allocation free)

When OutlineOpenVector is put into the dd, it will automatically save it there as an open outline
"""
struct OutlineOpenVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    data_is_open::Bool
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

@compat public OutlineOpenVector
push!(document[Symbol("Outlines")], :OutlineOpenVector)

# ===================

"""
    CircularVector{T,A<:AbstractVector{T}} <: AbstractVector{T}

Vector with circular indexes, without creating a new array (allocation free)

NOTE: It can be combined with `OutlineClosedVector` and/or `OutlineOpenVector`, like this: `CircularVector(OutlineClosedVector(x, x_is_closed))`
"""
struct CircularVector{T,A<:AbstractVector{T}} <: AbstractVector{T}
    data::A
    offset::Int
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

@compat public CircularVector
push!(document[Symbol("Outlines")], :CircularVector)

# ===================
"""
    is_open_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}

Determine if a polygon, defined by separate vectors for R and Z coordinates, is open.
"""
function is_open_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}
    return R[1] != R[end] || Z[1] != Z[end]
end

"""
    is_open_polygon(vertices::AbstractVector)::Bool
"""
function is_open_polygon(vertices::AbstractVector)::Bool
    r1 = vertices[1][1]
    z1 = vertices[1][2]
    rend = vertices[end][1]
    zend = vertices[end][2]
    return r1 != rend || z1 != zend
end

@compat public is_open_polygon
push!(document[Symbol("Outlines")], :is_open_polygon)

"""
    is_closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}

Determine if a polygon, defined by separate vectors for R and Z coordinates, is closed.
"""
function is_closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T})::Bool where {T<:Real}
    return !is_open_polygon(R, Z)
end

"""
    is_closed_polygon(vertices::AbstractVector)::Bool
"""
function is_closed_polygon(vertices::AbstractVector)::Bool
    return !is_open_polygon(vertices)
end

@compat public is_closed_polygon
push!(document[Symbol("Outlines")], :is_closed_polygon)

"""
    open_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, args...) where {T<:Real}

Returns a view of the vectors R and Z such that they are a open polygon

Returns a named tuple containing the status of the polygon (was_closed, was_open) and the views of the R and Z vectors.
"""
function open_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, args...) where {T<:Real}
    was_open = is_open_polygon(R, Z)
    R = OutlineOpenVector(R, was_open)
    Z = OutlineOpenVector(Z, was_open)
    args = collect(map(x -> OutlineOpenVector(x, was_open), args))
    return (was_closed=!was_open, was_open=was_open, R=R, Z=Z, r=R, z=Z, rz=(R, Z), args=args)
end

@compat public open_polygon
push!(document[Symbol("Outlines")], :open_polygon)

"""
    closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, args...) where {T<:Real}

Returns a view of the vectors R and Z such that they are a closed polygon

Returns a named tuple containing the status of the polygon (was_closed, was_open) and the views of the R and Z vectors.
"""
function closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, args...) where {T<:Real}
    was_closed = is_closed_polygon(R, Z)
    R = OutlineClosedVector(R, was_closed)
    Z = OutlineClosedVector(Z, was_closed)
    args = collect(map(x -> OutlineClosedVector(x, was_closed), args))
    return (was_closed=was_closed, was_open=!was_closed, R=R, Z=Z, r=R, z=Z, rz=(R, Z), args=args)
end

"""
    closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, closed::Bool, args...) where {T<:Real}

Returns a closed polygon depending on `closed` otherwise returns an open polygon
"""
function closed_polygon(R::AbstractVector{T}, Z::AbstractVector{T}, closed::Bool, args...) where {T<:Real}
    was_closed = is_closed_polygon(R, Z)
    was_open = !was_closed
    if closed
        R = OutlineClosedVector(R, was_closed)
        Z = OutlineClosedVector(Z, was_closed)
        args = collect(map(x -> OutlineClosedVector(x, was_closed), args))
    else
        R = OutlineOpenVector(R, was_open)
        Z = OutlineOpenVector(Z, was_open)
        args = collect(map(x -> OutlineOpenVector(x, was_open), args))
    end
    return (was_closed=was_closed, was_open=was_open, R=R, Z=Z, r=R, z=Z, rz=(R, Z), args=args)
end

@compat public closed_polygon
push!(document[Symbol("Outlines")], :closed_polygon)

"""
    is_clockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}

Returns true/false if polygon is defined clockwise
"""
function is_clockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}
    # Check if the vectors are of the same length
    if length(r) != length(z)
        error("Vectors must be of the same length")
    end

    area = zero(T)
    n = length(r)

    for i in 1:n-1
        area += (r[i] * z[i+1] - r[i+1] * z[i])
    end
    area += (r[n] * z[1] - r[1] * z[n])

    return area < 0
end

@compat public is_clockwise
push!(document[Symbol("Outlines")], :is_clockwise)

"""
    is_counterclockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}

Returns true/false if polygon is defined counterclockwise
"""
function is_counterclockwise(r::AbstractVector{T}, z::AbstractVector{T})::Bool where {T<:Real}
    return !is_clockwise(r, z)
end

@compat public is_counterclockwise
push!(document[Symbol("Outlines")], :is_counterclockwise)
