__precompile__()

module FUSE

#= ============== =#
#= DATA STRUCTURE =#
#= ============== =#

include("data.jl")

#= ========= =#
#= BOOTSTRAP =#
#= ========= =#

include("bootstrap.jl")

#= ========== =#
#= PARAMETERS =#
#= ========== =#

include("parameters.jl")

end # module
