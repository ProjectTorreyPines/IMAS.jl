module IMAS

using Printf

#= ====== =#
#= IMASDD =#
#= ====== =#
import IMASDD
# import all IMASDD.jl as if it was defined in IMAS.jl
for n in names(IMASDD; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(IMASDD), :eval, :include, :document)
        @eval import IMASDD: $n
    end
end
import IMASDD: @ddtime

#= ===== =#
#= UTILS =#
#= ===== =#
include("real.jl")
include("constants.jl")
include("math.jl")
include("constraints.jl")
include("objectives.jl")
include("extract.jl")
include("get_from.jl")
include("fxp.jl")

#= ======= =#
#= PHYSICS =#
#= ======= =#
include("physics.jl")

#= =========== =#
#= EXPRESSIONS =#
#= =========== =#
include(joinpath(["expressions", "onetime.jl"]))
include(joinpath(["expressions", "dynamic.jl"]))

#= ======== =#
#= PLOTTING =#
#= ======== =#
include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export @ddtime, constants, ±, force_float, extract

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end # module
