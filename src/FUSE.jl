__precompile__()

module FUSE


include("dd.jl")

"""
Resize StructArray to contain n elements.
If n is smaller than the current collection length, the first n elements will be retained.
If n is larger, the new elements are not guaranteed to be initialized.
"""
function resize!(collection::StructArray, n::Int)
    if n > length(collection)
        for k in length(collection):n - 1
            push!(collection, eltype(collection)())
            # println("push $(length(collection))")
        end
    elseif n < length(collection)
        for k in n:length(collection) - 1
            pop!(collection)
            # println("pop $(length(collection))")
        end
    end
    return collection
end

end # module
