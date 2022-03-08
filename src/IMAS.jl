__precompile__()

module IMAS

using Printf

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#
include("data.jl")
include("time.jl")

#= ================= =#
#= PHYSICS FUNCTIONS =#
#= ================= =#
include("real.jl")
include("constants.jl")
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
