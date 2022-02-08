__precompile__()

module IMAS

using Printf

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#

include("data.jl")

#= ================= =#
#= PHYSICS FUNCTIONS =#
#= ================= =#

include("time.jl")
include("constants.jl")
include("physics.jl")

#= ================== =#
#= PLOTTING FUNCTIONS =#
#= ================== =#

include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export is_missing, @ddtime, @coords, constants

end # module
