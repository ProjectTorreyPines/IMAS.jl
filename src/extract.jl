import OrderedCollections

# ==================== #
# extract data from dd #
# ==================== #
mutable struct ExtractFunction
    group::Symbol
    name::Symbol
    units::String
    func::Function
    # inner constructor to register ExtractFunction in ExtractFunctionsLibrary
    ExtractFunction(group::Symbol, name::Symbol, units::String, func::Function) = begin
        objf = new(group, name, units, func)
        ExtractFunctionsLibrary[objf.name] = objf
        return objf
    end
end

const ExtractFunctionsLibrary = OrderedCollections.OrderedDict{Symbol,ExtractFunction}()
function update_ExtractFunctionsLibrary!()
    empty!(ExtractFunctionsLibrary)
    ExtractFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
    ExtractFunction(:equilibrium, :δ, "-", dd -> dd.equilibrium.time_slice[].boundary.triangularity)
    ExtractFunction(:equilibrium, :ζ, "-", dd -> dd.equilibrium.time_slice[].boundary.squareness)
    ExtractFunction(:equilibrium, :B0, "T", dd -> @ddtime(dd.summary.global_quantities.b0.value))
    ExtractFunction(:equilibrium, :ip, "MA", dd -> @ddtime(dd.summary.global_quantities.ip.value) / 1e6)
    ExtractFunction(:equilibrium, :R0, "m", dd -> dd.summary.global_quantities.r0.value)
    ExtractFunction(:equilibrium, :βn, "-", dd -> @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
    ExtractFunction(:profiles, :Qfusion, "-", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / @ddtime(dd.summary.heating_current_drive.power_launched_total.value))
    ExtractFunction(:profiles, :Pfusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6)
    ExtractFunction(:profiles, :zeff, "-", dd -> @ddtime(dd.summary.volume_average.zeff.value))
    ExtractFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ExtractFunction(:profiles, :Ti0, "keV", dd -> dd.core_profiles.profiles_1d[].ion[1].temperature[1] / 1E3)
    ExtractFunction(:balance_of_plant, :Pelectric_net, "MWe", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6)
    ExtractFunction(:balance_of_plant, :Qplant, "-", dd -> @ddtime(dd.balance_of_plant.Q_plant))
    ExtractFunction(:heating_current_drive, :Pec, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ec.value))
    ExtractFunction(:heating_current_drive, :Pnbi, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_nbi.value))
    ExtractFunction(:heating_current_drive, :Pic, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ic.value))
    ExtractFunction(:heating_current_drive, :Plh, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_lh.value))
    ExtractFunction(:heating_current_drive, :Paux_total, "W", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_total.value))
    ExtractFunction(:costing, :levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE)
    ExtractFunction(:costing, :capital_cost, "\$M", dd -> dd.costing.cost_direct_capital.cost)
    ExtractFunction(:build, :flattop, "Hours", dd -> dd.build.oh.flattop_duration / 3600.0)
end
update_ExtractFunctionsLibrary!()

"""
    (ef::ExtractFunction)(dd::IMAS.dd)

Run the extract function
"""
function (ef::ExtractFunction)(dd::IMAS.dd)
    return ef.func(dd)
end

function Base.show(io::IO, f::ExtractFunction)
    printstyled(io, f.group; bold=true)
    printstyled(io, "."; bold=true)
    printstyled(io, f.name; bold=true, color=:blue)
    print(io, " →")
end

"""
    extract(dd::IMAS.dd, xtract::AbstractDict{Symbol,T}=ExtractFunctionsLibrary)::Dict{Symbol,Any} where {T<:Union{Function,ExtractFunction}}

Extract data from `dd``.
By default, the `ExtractFunctionsLibrary` is used.
Each of the `xtract` functions should accept `dd` as input, like this:

    xtract = Dict(
            :beta_normal => dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal,
            :time => dd -> @ddtime(dd.equilibrium.time)
        )
"""
function extract(dd::IMAS.dd, xtract::AbstractDict{Symbol,T}=ExtractFunctionsLibrary)::Dict{Symbol,Any} where {T<:Union{Function,ExtractFunction}}
    results = Dict{Symbol,Any}()
    for key in keys(xtract)
        if dd === missing
            results[key] = NaN
            continue
        end
        try
            res[key] = xtract[key](dd)
        catch
            res[key] = NaN
        end
    end
    return results
end
