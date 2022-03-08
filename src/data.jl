import JSON

include("dd.jl")

include("f2.jl")

"""
    assign_expressions(ids::Union{IDS,IDSvector})

Assign expressions to a IDS/IDSvector
NOTE: This is done not recursively
"""
function assign_expressions(ids::Union{IDS,IDSvector})
    struct_name = f2u(ids)
    for item in children(ids)
        if typeof(getfield(ids, item)) <: Union{IDS,IDSvector}
            continue
        elseif "$(struct_name).$(item)" in keys(expressions)
            setproperty!(ids, item, expressions["$(struct_name).$(item)"])
        end
    end
    return ids
end

#= ===================== =#
#  IDS related functions  #
#= ===================== =#

"""
    coords ids.path.to.array.y => interpolating_x

Macro for interpolating data
"""
macro coords(ex)
    return _coords(ex)
end

function _coords(ex)
    quote
        local expr = $(esc(Meta.QuoteNode(ex)))
        if (expr.head !== :call) || (expr.args[1] != :(=>))
            error("@coords must use `dd.ids.field => xx` syntax")
        end
        local ids = $(esc(ex.args[2].args[1]))
        local field = $(esc(ex.args[2].args[2]))
        local xx = $(esc(ex.args[3]))
        local yy = interp1d(ids, field).(xx)
        yy
    end
end
