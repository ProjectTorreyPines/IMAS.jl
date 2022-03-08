include("dd.jl")

include("f2.jl")

expressions = Dict{String,Function}()

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
