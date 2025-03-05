"""
    name::Symbol
    units::String
    func::Function
    target::Float64
"""
mutable struct ObjectiveFunction
    name::Symbol
    units::String
    func::Function
    target::Float64
    # inner constructor to register ObjectiveFunction in ObjectiveFunctionsLibrary
    ObjectiveFunction(name::Symbol, units::String, func::Function, target::Float64) = begin
        objf = new(name, units, func, target)
        ObjectiveFunctionsLibrary[objf.name] = objf
        return objf
    end
end

@compat public ObjectiveFunction
push!(document[Symbol("Functions library")], :ObjectiveFunction)

"""
    (objf::ObjectiveFunction)(dd::IMAS.dd)

From real domain to objective domain (Metaheuristics will always minimize)
"""
function (objf::ObjectiveFunction)(dd::IMAS.dd)
    if isinf(objf.target)
        if objf.target < 0
            return objf.func(dd)
        else
            return -objf.func(dd)
        end
    elseif objf.target == 0.0
        return abs(objf.func(dd))
    else
        return abs(objf.func(dd) - objf.target) / objf.target
    end
end

"""
    (objf::ObjectiveFunction)(x::Float64)

From objective domain to real domain
"""
function (objf::ObjectiveFunction)(x::Float64)
    if isinf(objf.target)
        if objf.target < 0
            return x
        else
            return -x
        end
    elseif objf.target == 0.0
        return x
    else
        return x * objf.target + objf.target
    end
end

function Base.show(io::IO, ::MIME"text/plain", objf::ObjectiveFunction)
    printstyled(io, objf.name; bold=true, color=:blue)
    print(io, " →")
    print(io, " $(objf.target)")
    return print(io, " [$(objf.units)]")
end

function Base.show(io::IO, x::MIME"text/plain", objfs::AbstractDict{Symbol,ObjectiveFunction})
    for objf in objfs
        show(io, x, objf)
        println(io, "")
    end
end

function Base.string(objf::ObjectiveFunction)
    return "$(objf.name) → $(objf.target) [$(objf.units)]"
end

# ========================= #
# ObjectiveFunctionsLibrary #
# ========================= #
const ObjectiveFunctionsLibrary = OrderedCollections.OrderedDict{Symbol,ObjectiveFunction}()

function update_ObjectiveFunctionsLibrary!()
    empty!(ObjectiveFunctionsLibrary)
    #! format: off
    ObjectiveFunction(:min_levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE, -Inf)
    ObjectiveFunction(:min_log10_levelized_CoE, "log₁₀(\$/kW)", dd -> log10(dd.costing.levelized_CoE), -Inf)
    ObjectiveFunction(:min_capital_cost, "\$B", dd -> dd.costing.cost_direct_capital.cost / 1E3, -Inf)
    ObjectiveFunction(:max_fusion, "MW", dd -> fusion_power(dd.core_profiles.profiles_1d[]) / 1E6, Inf)
    ObjectiveFunction(:max_power_electric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6, Inf)
    ObjectiveFunction(:req_power_electric_net, "ΔMW", dd -> abs(@ddtime(dd.balance_of_plant.power_electric_net) - dd.requirements.power_electric_net) / 1E6, 0.0)
    ObjectiveFunction(:max_flattop, "hours", dd -> dd.build.oh.flattop_duration / 3600.0, Inf)
    ObjectiveFunction(:max_q95, "", dd -> dd.equilibrium.time_slice[].global_quantities.q_95, Inf)
    ObjectiveFunction(:req_flattop, "Δhours", dd -> abs(dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / 3600.0, 0.0)
    ObjectiveFunction(:max_log10_flattop, "log₁₀(hours)", dd -> log10(dd.build.oh.flattop_duration / 3600.0), Inf)
    ObjectiveFunction(:min_βn, "", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal, -Inf)
    ObjectiveFunction(:min_R0, "m", dd -> dd.equilibrium.time_slice[].boundary.geometric_axis.r, -Inf)
    ObjectiveFunction(:max_zeff, "", dd -> @ddtime(dd.summary.volume_average.zeff.value), Inf)
    #! format: on
    return ObjectiveFunctionsLibrary
end
update_ObjectiveFunctionsLibrary!()

Base.Docs.@doc """
    ObjectiveFunctionsLibrary::OrderedCollections.OrderedDict{Symbol,ObjectiveFunction}

Collection of ObjectiveFunction
* `$(join(collect(map(string,values(ObjectiveFunctionsLibrary))),"`\n* `"))`
""" ObjectiveFunctionsLibrary

@compat public ObjectiveFunctionsLibrary
push!(document[Symbol("Functions library")], :ObjectiveFunctionsLibrary)