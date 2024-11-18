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
document[Symbol("get from")] = Symbol[]
include("get_from.jl")

#= ==== =#
#= MATH =#
#= ==== =#
include(joinpath("math", "geometry.jl"))
include(joinpath("math", "math.jl"))

#= ======= =#
#= EXTRACT =#
#= ======= =#
document[Symbol("Functions library")] = Symbol[]
document[:Extract] = Symbol[]
include(joinpath("extract", "constraints.jl"))
include(joinpath("extract", "objectives.jl"))
include(joinpath("extract", "extract.jl"))

#= ======= =#
#= CONTROL =#
#= ======= =#
include(joinpath("control", "control.jl"))
include(joinpath("control", "fxp.jl"))

#= ======= =#
#= PHYSICS =#
#= ======= =#
document[:Physics] = Symbol[]
include("physics.jl")

#= =========== =#
#= EXPRESSIONS =#
#= =========== =#
document[:Expressions] = Symbol[]
include(joinpath("expressions", "onetime.jl"))
include(joinpath("expressions", "dynamic.jl"))

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
