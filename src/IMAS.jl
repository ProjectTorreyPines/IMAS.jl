module IMAS

using Printf
using Compat:@compat
import OrderedCollections
import MacroTools
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
include("real.jl")
include("get_from.jl")
include("geometry.jl")
include("math.jl")

#= ======= =#
#= EXTRACT =#
#= ======= =#
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
include("physics.jl")

#= =========== =#
#= EXPRESSIONS =#
#= =========== =#
include(joinpath("expressions", "onetime.jl"))
include(joinpath("expressions", "dynamic.jl"))

#= ======== =#
#= PLOTTING =#
#= ======== =#
include("plot.jl")

#= ====== =#
#= EXPORT =#
#= ====== =#
export @ddtime, @findall, @explain

end # module
