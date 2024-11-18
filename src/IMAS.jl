module IMAS

using Printf
using Compat:@compat
import OrderedCollections
const document = OrderedCollections.OrderedDict()

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
import IMASdd: @ddtime, @findall

#= ===== =#
#= UTILS =#
#= ===== =#
document[:Real] = Symbol[]
include("real.jl")
include("math.jl")
document[Symbol("get from")] = Symbol[]
include("get_from.jl")
include("fxp.jl")

#= ======= =#
#= EXTRACT =#
#= ======= =#
document[Symbol("Functions library")] = Symbol[]
document[:Extract] = Symbol[]
include("constraints.jl")
include("objectives.jl")
include("extract.jl")

#= ======= =#
#= PHYSICS =#
#= ======= =#
document[:Physics] = Symbol[]
include("constants.jl")
include("physics.jl")

#= =========== =#
#= EXPRESSIONS =#
#= =========== =#
include(joinpath(["expressions", "onetime.jl"]))
include(joinpath(["expressions", "dynamic.jl"]))

#= ======== =#
#= PLOTTING =#
#= ======== =#
document[:Plot] = Symbol[]
include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export @ddtime, @findall

end # module
