__precompile__()

module IMAS

struct GlobalTime end
const τ = GlobalTime()

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#

include("data.jl")

#= ================= =#
#= PHYSICS FUNCTIONS =#
#= ================= =#

include("time.jl")

include("physics.jl")

#= ================== =#
#= PLOTTING FUNCTIONS =#
#= ================== =#

include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export f2p, f2u, p2i, i2p, u2f, f2f, f2fs, u2fs, f2i
export coordinates, set_timedep_value!
export IDS, IDSvector, IDSvectorElement
export top, top_dd, top_ids, parent, children, expressions

export is_missing, timedep

end # module
