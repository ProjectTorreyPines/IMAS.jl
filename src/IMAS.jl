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
export coordinates, set_timedep_value!, insert_time_index, is_missing, insert_time_slice!
export IDS, IDSvector, IDSvectorElement
export top, top_dd, top_ids, parent, children, expressions

end # module
