__precompile__()

module IMAS

using Printf

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#
include("real.jl")
include("data.jl")

#= ================= =#
#= PHYSICS FUNCTIONS =#
#= ================= =#
include("constants.jl")
include("time.jl")
include("math.jl")
include("physics.jl")
include("expressions.jl")

#= ================== =#
#= PLOTTING FUNCTIONS =#
#= ================== =#
include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export is_missing, @ddtime, @coords, constants, ±, force_float

end # module
