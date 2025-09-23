module IMAS

using Printf
using Compat:@compat
import OrderedCollections
const document = OrderedCollections.OrderedDict()

macro import_all(mod)
    syms = names(eval(mod); all=true)
    syms = filter(n -> Base.isidentifier(n) && n ∉ (mod, :eval, :include, :document), syms)
    imports = [:(import $mod: $(s)) for s in syms]
    return Expr(:toplevel, imports...)
end

function __init__()
    # register IMAS expressions with IMASdd
    IMASdd.set_dynamic_expressions(dynamic_expressions)
    IMASdd.set_onetime_expressions(onetime_expressions)
end

#= ====== =#
#= IMASdd =#
#= ====== =#
import IMASdd
@import_all IMASdd
import IMASdd: @ddtime, @findall
import AbstractTrees: print_tree

#= ===== =#
#= UTILS =#
#= ===== =#
include("real.jl")
include("get_from.jl")
include("geometry.jl")
include("math.jl")
include("benchmark.jl")
include("signal.jl")

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
include("experiments.jl")

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
export @ddtime, @findall, help

import HelpPlots: help_plot, help_plot!
export help_plot, help_plot!
export print_tree

end # module
