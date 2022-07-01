__precompile__()

module IMAS

using Printf
import BSON

#= ======= =#
#= IMAS DD =#
#= ======= =#
import IMASDD
# import all IMASDD.jl as if it was defined in IMAS.jl
for n in names(IMASDD; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(IMASDD), :eval, :include)
        @eval import IMASDD: $n
    end
end
import IMASDD: @ddtime, interp1d

#= ========= =#
#= SAVE/LOAD =#
#= ========= =#
"""
    save(@nospecialize(ids::IDS), filename::AbstractString)

Save IDS data to BSON file
"""
function save(@nospecialize(ids::IDS), filename::AbstractString)
    BSON.bson(filename, Dict(:data=>ids))
end

"""
    load(filename::AbstractString)

Load IDS data from BSON file
"""
function load(filename::AbstractString)
    BSON.load(filename, @__MODULE__)[:data]
end

#= ===== =#
#= UTILS =#
#= ===== =#
include("real.jl")
include("constants.jl")
include("math.jl")

#= ======= =#
#= PHYSICS =#
#= ======= =#
include("physics.jl")
include("expressions.jl")

#= ======== =#
#= PLOTTING =#
#= ======== =#
include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export @ddtime, constants, ±, force_float

end # module
