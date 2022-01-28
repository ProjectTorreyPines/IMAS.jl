__precompile__()

module IMAS

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
export is_missing, @ddtime

end # module
