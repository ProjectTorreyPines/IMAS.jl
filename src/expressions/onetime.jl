document[:Expressions] = Symbol[]

function IMASdd.get_expressions(::Type{Val{:onetime}})
    return onetime_expressions
end

const onetime_expressions = Dict{String,Function}()
otexp = onetime_expressions

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
# These expressions are frozen the first time they are accessed.
# This is necessary to ensure that core_profiles, core_sources, and core_transport grids do not change after changing the equilibrium.
# The idea is that we want to freeze in the core_profiles, core_sources, and core_transport grids the rho, psi, volume, area, ... info that were used when those IDSs were filled.
# While this is generally ok, this is not desirable when iterating the equilibrium solver with other actors.
# In this case, at each iteration we want core_profiles, core_sources, and core_transport to take the grids from the latest equilibrium.
#
# NOTE: make sure that expressions accept as argument (not keyword argument)
# the coordinates of the quantitiy you are writing the expression of
#
# For example, this will FAIL:
#    otexp["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
#         (; electrons, _...) -> electrons.temperature .* electrons.density_thermal * 1.60218e-19
#
# This is GOOD:
#    otexp["core_profiles.profiles_1d[:].electrons.pressure_thermal"] =
#         (rho_tor_norm; electrons, _...) -> electrons.temperature .* electrons.density_thermal * 1.60218e-19

#= =========== =#
# core_profiles #
#= =========== =#
@store_expr otexp["core_profiles.profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

@store_expr otexp["core_profiles.profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(volume), :cubic).(rho_tor_norm) .^ 2
    end

@store_expr otexp["core_profiles.profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(area), :cubic).(rho_tor_norm) .^ 2
    end

@store_expr otexp["core_profiles.profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end


@store_expr otexp["core_profiles.profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        sign_psi = sign(psi[end] - psi[1])
        return sign_psi .* (interp1d(eqt.profiles_1d.rho_tor_norm, sqrt.(abs.(psi .- psi[1])), :cubic).(rho_tor_norm) .^ 2) .+ psi[1]
    end

#= ============ =#
# core_transport #
#= ============ =#
@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_flux.rho_tor_norm"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.grid_d.rho_tor_norm

@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_d.rho_tor_norm"] =
    (rho_tor_norm; profiles_1d, _...) -> profiles_1d.grid_flux.rho_tor_norm

@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_flux.psi_norm"] =
    (rho_tor_norm; grid_flux, _...) -> norm01(grid_flux.psi)

@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_flux.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_flux.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_flux.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

@store_expr otexp["core_transport.model[:].profiles_1d[:].grid_flux.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

#= ============ =#
#  core_sources  #
#= ============ =#
@store_expr otexp["core_sources.source[:].profiles_1d[:].grid.psi_norm"] =
    (rho_tor_norm; grid, _...) -> norm01(grid.psi)

@store_expr otexp["core_sources.source[:].profiles_1d[:].grid.volume"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        volume = eqt.profiles_1d.volume
        return interp1d(eqt.profiles_1d.rho_tor_norm, volume, :cubic).(rho_tor_norm)
    end

@store_expr otexp["core_sources.source[:].profiles_1d[:].grid.area"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        area = eqt.profiles_1d.area
        return interp1d(eqt.profiles_1d.rho_tor_norm, area, :cubic).(rho_tor_norm)
    end

@store_expr otexp["core_sources.source[:].profiles_1d[:].grid.surface"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        surface = eqt.profiles_1d.surface
        return interp1d(eqt.profiles_1d.rho_tor_norm, surface, :cubic).(rho_tor_norm)
    end

@store_expr otexp["core_sources.source[:].profiles_1d[:].grid.psi"] =
    (rho_tor_norm; dd, profiles_1d, _...) -> begin
        eqt = dd.equilibrium.time_slice[Float64(profiles_1d.time)]
        psi = eqt.profiles_1d.psi
        return interp1d(eqt.profiles_1d.rho_tor_norm, psi, :cubic).(rho_tor_norm)
    end

# ============ #

Base.Docs.@doc """
    onetime_expressions = Dict{String,Function}()

Expressions that are frozen after first evaluation
* `$(join(sort!(collect(keys(onetime_expressions))),"`\n* `"))`
""" onetime_expressions

push!(document[:Expressions], :onetime_expressions)