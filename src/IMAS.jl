module IMAS

using Printf

#= ====== =#
#= IMASdd =#
#= ====== =#
import IMASdd
# import all IMASdd.jl as if it was defined in IMAS.jl
for n in names(IMASdd; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(IMASdd), :eval, :include, :document)
        @eval import IMASdd: $n
    end
end
import IMASdd: @ddtime

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
export @ddtime, constants, ±, extract, help_plot, help_plot!

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end # module
