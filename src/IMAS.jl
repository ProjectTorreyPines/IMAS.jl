__precompile__()

module IMAS

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#

include("data.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export f2p, f2u, p2i, i2p, u2f, f2f, f2fs, u2fs
export coordinates
export IDS, IDSvector, IDSvectorElement
export dd
export top, parent, children, expressions

end # module
