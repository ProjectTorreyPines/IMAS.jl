function IMASdd.get_expr_info_dict(::Type{Val{:onetime_and_dynamic}})
    return expr_info_dict
end

struct ExprInfo
    name::String
    args::String
    body::String
end

const expr_info_dict = Dict{String,ExprInfo}()

function Base.show(io::IO, ::MIME"text/plain", expr_info::ExprInfo)
    printstyled(io, "-"^length(expr_info.name) * "\n"; color=:red)
    printstyled(io, "$(expr_info.name)\n"; color=:red, bold=true)
    printstyled(io, "-"^length(expr_info.name) * "\n"; color=:red)
    printstyled(io, "[Args]\n"; color=:blue)
    print(io, "$(expr_info.args)\n")
    printstyled(io, "[Body]\n"; color=:blue)
    return print(io, "$(expr_info.body)\n")
end

"""
    @store_expr expr

Takes an assignment expression of the form `dict["key"] = (args...) -> body` and:

 1. Executes the original assignment.
 2. Extracts the argument list and body from the given anonymous function.
 3. Stores this information as an `ExprInfo` object in `expr_info_dict` under the specified key.
"""
macro store_expr(expr)
    @assert expr.head == :(=)
    lhs = expr.args[1]
    rhs = expr.args[2]

    @assert lhs.head == :ref && length(lhs.args) == 2
    key_expr = lhs.args[2]

    expr_args = sprint(show, MacroTools.prettify(rhs.args[1]))
    expr_body = sprint(show, MacroTools.prettify(rhs.args[2]))

    expr_info_dict = :expr_info_dict

    return quote
        # execute the given original expression
        $(esc(expr))
        # store the ExprInfo into expr_info_dict
        $expr_info_dict[$(key_expr)] = ExprInfo($key_expr, $expr_args, $expr_body)
    end
end

"""
    @explain expr

Temporarily sets a flag to explain the execution of `expr`. When `expr` runs under this macro,
`IMASdd.FLAG_EXPLAIN_EXPRESSION` is enabled, allowing additional diagnostic or explanatory
information to be displayed. After `expr` completes (successfully or not), the flag is reset.
"""
macro explain(expr)
    return quote
        IMASdd.FLAG_EXPLAIN_EXPRESSION[] = true
        try
            $(esc(expr))
        finally
            IMASdd.FLAG_EXPLAIN_EXPRESSION[] = false
        end
    end
end