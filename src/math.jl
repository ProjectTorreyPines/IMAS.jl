document[Symbol("Math")] = Symbol[]

"""
    norm01(x::T)::T where {T<:AbstractVector{<:Real}}

Normalize a vector so that the first item in the array is 0 and the last one is 1
This is handy where psi_norm should be used (and IMAS does not define a psi_norm array)
"""
function norm01(x::T)::T where {T<:AbstractVector{<:Real}}
    return (x .- x[1]) ./ (x[end] .- x[1])
end

@compat public norm01
push!(document[Symbol("Math")], :norm01)

"""
    to_range(vector::AbstractVector)

Turn a vector into a range (if possible)
"""
function to_range(vector::AbstractVector{<:Real})
    N = length(vector)
    dv = vector[2] - vector[1]

    if maximum(abs(vector[k] - vector[k-1]) - dv for k in 2:N) > dv / 1000.0
        error("to_range requires vector data to be equally spaced: $(vector)")
    end
    return range(vector[1], vector[end], N)
end

@compat public to_range
push!(document[Symbol("Math")], :to_range)

"""
    meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}

Return coordinate matrices from coordinate vectors
"""
function meshgrid(x::AbstractVector{T}, y::AbstractVector{T}) where {T}
    return last.(Iterators.product(y, x)), first.(Iterators.product(y, x))
end

@compat public meshgrid
push!(document[Symbol("Math")], :meshgrid)

"""
    calc_z(x::AbstractVector{<:Real}, f::AbstractVector{<:Real}, method::Symbol)

Returns the gradient scale lengths of vector f on x

The finite difference `method` of the gradient can be one of [:backward, :central, :forward, :second_order, :third_order]

NOTE: the inverse scale length is NEGATIVE for typical density/temperature profiles
"""
function calc_z(x::AbstractVector{<:Real}, f::AbstractVector{<:Real}, method::Symbol)
    g = gradient(x, f; method)
    return g ./ f
end

@compat public calc_z
push!(document[Symbol("Math")], :calc_z)

"""
    integ_z(rho::AbstractVector{<:Real}, z_profile::AbstractVector{<:Real}, bc::Real)

Backward integration of inverse scale length vector with given edge boundary condition
"""
function integ_z(rho::AbstractVector{<:Real}, z_profile::AbstractVector{<:Real}, bc::Real)
    profile_new = similar(rho)
    profile_new[end] = bc
    for i in length(rho)-1:-1:1
        a = rho[i]
        fa = z_profile[i]
        b = rho[i+1]
        fb = z_profile[i+1]
        trapz_integral = (b - a) * (fa + fb) / 2.0
        profile_new[i] = profile_new[i+1] * exp(trapz_integral)
    end
    return profile_new
end

@compat public integ_z
push!(document[Symbol("Math")], :integ_z)

"""
    pack_grid_gradients(x::AbstractVector{T}, y::AbstractVector{T}; n_points::Int=length(x), l::Float64=1E-2) where {T<:Float64}

Returns grid between `minimum(x)` and `maximum(x)` with `n_points` points positioned to
sample `y(x)` in such a way to pack more points where gradients are greates.

`l` controls how much the adaptive gradiant sampling should approach linear sampling.
"""
function pack_grid_gradients(x::AbstractVector{T}, y::AbstractVector{T}; n_points::Int=length(x), l::Float64=1E-2) where {T<:Float64}
    tmp = abs.(gradient(x, y))
    m, M = extrema(tmp)
    tmp .+= (M - m) * l
    cumsum!(tmp, tmp)
    tmp .-= tmp[1]
    tmp ./= tmp[end]
    return interp1d(tmp, x).(range(0.0, 1.0, n_points))
end

@compat public pack_grid_gradients
push!(document[Symbol("Math")], :pack_grid_gradients)

"""
    unique_indices(vec::AbstractVector)::Vector{Int}

Return the indices of the first occurrence of each unique element in the input vector `vec`
"""
function unique_indices(vec::AbstractVector)::Vector{Int}
    uniq_elements = unique(vec)
    return [findfirst(==(elem), vec) for elem in uniq_elements]
end

@compat public unique_indices
push!(document[Symbol("Math")], :unique_indices)

"""
    index_circular(N::Int, idx::Int)

If `idx` is beyond N or less than 1, it wraps around in a circular manner.
"""
function index_circular(N::Int, idx::Int)
    return mod(idx - 1, N) + 1
end

@compat public index_circular
push!(document[Symbol("Math")], :index_circular)

"""
    getindex_circular(vec::AbstractVector{T}, idx::Int)::T where {T}

Return the element of the vector `vec` at the position `idx`.

If `idx` is beyond the length of `vec` or less than 1, it wraps around in a circular manner.
"""
function getindex_circular(vec::AbstractVector{T}, idx::Int)::T where {T}
    cidx = index_circular(length(vec), idx)
    return vec[cidx]
end

@compat public getindex_circular
push!(document[Symbol("Math")], :getindex_circular)

"""
    chunk_indices(dims::Tuple{Vararg{Int}}, N::Int)

Split the indices of an array with dimensions `dims` into `N` chunks of similar size.
Each chunk is a generator of `CartesianIndex` objects.

# Arguments

  - `dims::Tuple{Vararg{Int}}`: A tuple specifying the dimensions of the array. For a 2D array, this would be `(rows, cols)`.
  - `N::Int`: The number of chunks to split the indices into.

# Returns

  - Vector with chunks of cartesian indices. Each chunk can be iterated over to get the individual `CartesianIndex` objects.
"""
function chunk_indices(dims::Tuple{Vararg{Int}}, N::Int)
    # Total number of elements
    total_elements = prod(dims)

    # Calculate chunk size
    chunk_size, remainder = divrem(total_elements, N)

    # Split indices into N chunks
    chunks = []
    start_idx = 1
    for i in 1:N
        end_idx = start_idx + chunk_size - 1
        # Distribute the remainder among the first few chunks
        if i <= remainder
            end_idx += 1
        end
        chunk_1d = start_idx:end_idx
        chunk_multi = (CartesianIndices(dims)[i] for i in chunk_1d)
        push!(chunks, chunk_multi)
        start_idx = end_idx + 1
    end

    return chunks
end

@compat public chunk_indices
push!(document[Symbol("Math")], :chunk_indices)