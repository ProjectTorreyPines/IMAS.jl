module IMAS

using Printf

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

#= ===== =#
#= UTILS =#
#= ===== =#
include("real.jl")
include("constants.jl")
include("math.jl")
include("extract.jl")

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
export @ddtime, constants, ±, force_float, extract

end # module
