import OrderedCollections

# ==================== #
# extract data from dd #
# ==================== #
mutable struct ExtractFunction
    group::Symbol
    name::Symbol
    units::String
    func::Function
    value::Any
end

function ExtractFunction(group::Symbol, name::Symbol, units::String, func::Function)
    return ExtractFunction(group, name, units, func, NaN)
end

function ExtractLibFunction(group::Symbol, name::Symbol, units::String, func::Function)
    xfun = ExtractFunction(group, name, units, func)
    push!(ExtractFunctionsLibrary, xfun)
end

const ExtractFunctionsLibrary = ExtractFunction[]
function update_ExtractFunctionsLibrary!()
    empty!(ExtractFunctionsLibrary)
    ExtractLibFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
    ExtractLibFunction(:equilibrium, :δ, "-", dd -> dd.equilibrium.time_slice[].boundary.triangularity)
    ExtractLibFunction(:equilibrium, :ζ, "-", dd -> dd.equilibrium.time_slice[].boundary.squareness)
    ExtractLibFunction(:equilibrium, :B0, "T", dd -> @ddtime(dd.summary.global_quantities.b0.value))
    ExtractLibFunction(:equilibrium, :ip, "MA", dd -> @ddtime(dd.summary.global_quantities.ip.value) / 1e6)
    ExtractLibFunction(:equilibrium, :R0, "m", dd -> dd.summary.global_quantities.r0.value)
    ExtractLibFunction(:equilibrium, :βn, "-", dd -> @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
    ExtractLibFunction(:profiles, :Qfusion, "-", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / @ddtime(dd.summary.heating_current_drive.power_launched_total.value))
    ExtractLibFunction(:profiles, :Pfusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6)
    ExtractLibFunction(:profiles, :zeff, "-", dd -> @ddtime(dd.summary.volume_average.zeff.value))
    ExtractLibFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ExtractLibFunction(:profiles, :Ti0, "keV", dd -> dd.core_profiles.profiles_1d[].ion[1].temperature[1] / 1E3)
    ExtractLibFunction(:balance_of_plant, :Pelectric_net, "MWe", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6)
    ExtractLibFunction(:balance_of_plant, :Qplant, "-", dd -> @ddtime(dd.balance_of_plant.Q_plant))
    ExtractLibFunction(:heating_current_drive, :Pec, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ec.value))
    ExtractLibFunction(:heating_current_drive, :Pnbi, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_nbi.value))
    ExtractLibFunction(:heating_current_drive, :Pic, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ic.value))
    ExtractLibFunction(:heating_current_drive, :Plh, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_lh.value))
    ExtractLibFunction(:heating_current_drive, :Paux_total, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_total.value))
    ExtractLibFunction(:costing, :levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE)
    ExtractLibFunction(:costing, :capital_cost, "\$M", dd -> dd.costing.cost_direct_capital.cost)
    ExtractLibFunction(:build, :flattop, "Hours", dd -> dd.build.oh.flattop_duration / 3600.0)
end
update_ExtractFunctionsLibrary!()

"""
    (xfun::ExtractFunction)(dd::IMAS.dd)

Run the extract function and store its result in xfun.value

NOTE: NaN is assigned on error
"""
function (xfun::ExtractFunction)(dd::IMAS.dd)
    try
        xfun.value = xfun.func(dd)
    catch
        xfun.value = NaN
    end
    return xfun.value
end

function Base.show(io::IO, f::ExtractFunction)
    printstyled(io, f.group; bold=true)
    printstyled(io, "."; bold=true)
    printstyled(io, f.name; bold=true, color=:blue)
    printstyled(io, " → ")
    if any(isnan.(f.value))
        printstyled(io, f.value; color=:red)
    else
        printstyled(io, @sprintf("%.3g", f.value))
    end
    printstyled(io, " [$(f.units)]")
end

function Base.show(io::IO, ::MIME"text/plain", xtract::Vector{ExtractFunction})
    for xfun in xtract
        if !isnan(xfun.value)
            show(io, xfun)
            println(io, "")
        end
    end
end

"""
    extract(dd::IMAS.dd, xtract::Vector{ExtractFunction}=ExtractFunctionsLibrary)::Vector{ExtractFunction}

Extract data from `dd`. Each of the `ExtractFunction` should accept `dd` as input, like this:

    xtract = IMAS.ExtractFunction[
        :κ => ExtractFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
        :Te0 => ExtractFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ]

By default, the `ExtractFunctionsLibrary` is used.
"""
function extract(dd::IMAS.dd, xtract::Vector{ExtractFunction}=ExtractFunctionsLibrary)::Vector{ExtractFunction}
    results = deepcopy(xtract)
    for xfun in xtract
        xfun(dd)
    end
    return results
end
