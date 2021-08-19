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

#= ====== =#
#= EXPORT =#
#= ====== =#
export f2p,f2u,p2i,i2p,u2f,f2f,f2fs,u2fs
export coordinates
export FDS,FDSvector
export dd
export NumericalFDVector,AnalyticalFDVector,top

end # module
