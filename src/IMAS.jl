__precompile__()

module IMAS

using Printf
import JLD2

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

Save IDS data to file (JLD2 format)
"""
function save(@nospecialize(ids::IDS), filename::AbstractString)
    JLD2.jldsave(filename; ids)
end

"""
    load(filename::AbstractString)

Load IDS data from file (JLD2 format)
"""
function load(filename::AbstractString)
    JLD2.jldopen(filename, "r") do file
        file["ids"]
    end
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

#= ====== =#
#= DIGEST =#
#= ====== =#
include("digest.jl")

#= ======== =#
#= PLOTTING =#
#= ======== =#
include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export @ddtime, constants, ±, force_float

end # module
