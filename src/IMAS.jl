__precompile__()

module IMAS

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#

include("data.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export f2p, f2u, p2i, i2p, u2f, f2f, f2fs, u2fs, f2i
export coordinates, set_field_time_array, get_time_index, is_missing
export IDS, IDSvector, IDSvectorElement
export top, parent, children, expressions

end # module
