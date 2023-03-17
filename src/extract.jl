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
    ExtractFunctionsLibrary[Symbol("$(group)_$name")] = xfun
    return xfun
end

const ExtractFunctionsLibrary = EFL = OrderedCollections.OrderedDict{Symbol,ExtractFunction}()

function update_ExtractFunctionsLibrary!()
    empty!(EFL)
    ExtractLibFunction(:geometry, :R0, "m", dd -> dd.equilibrium.time_slice[].boundary.geometric_axis.r)
    ExtractLibFunction(:geometry, :a, "m", dd -> dd.equilibrium.time_slice[].boundary.minor_radius)
    ExtractLibFunction(:geometry, Symbol("1/ϵ"), "m", dd -> EFL[:geometry_R0](dd) / EFL[:geometry_a](dd))
    ExtractLibFunction(:geometry, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
    ExtractLibFunction(:geometry, :δ, "-", dd -> dd.equilibrium.time_slice[].boundary.triangularity)
    ExtractLibFunction(:geometry, :ζ, "-", dd -> dd.equilibrium.time_slice[].boundary.squareness)

    ExtractLibFunction(:equilibrium, :B0, "T", dd -> @ddtime(dd.summary.global_quantities.b0.value))
    ExtractLibFunction(:equilibrium, :ip, "MA", dd -> @ddtime(dd.summary.global_quantities.ip.value) / 1e6)
    ExtractLibFunction(:equilibrium, :q95, "-", dd -> dd.equilibrium.time_slice[].global_quantities.q_95)
    ExtractLibFunction(:equilibrium, :βpol, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_pol)
    ExtractLibFunction(:equilibrium, :βtor, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_tor)
    ExtractLibFunction(:equilibrium, :βn, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal)

    ExtractLibFunction(:profiles, :Te0, "keV", dd -> @ddtime(dd.summary.local.magnetic_axis.t_e.value) / 1E3)
    ExtractLibFunction(:profiles, :Ti0, "keV", dd -> @ddtime(dd.summary.local.magnetic_axis.t_i_average.value) / 1E3)
    ExtractLibFunction(:profiles, :ne0, "m⁻³", dd -> @ddtime(dd.summary.local.magnetic_axis.n_e.value))
    ExtractLibFunction(:profiles, :P0, "MPa", dd -> dd.core_profiles.profiles_1d[].pressure[1] / 1E6)

    ExtractLibFunction(:profiles, Symbol("<Te>"), "keV", dd -> @ddtime(dd.summary.volume_average.t_e.value) / 1E3)
    ExtractLibFunction(:profiles, Symbol("<Ti>"), "keV", dd -> @ddtime(dd.summary.volume_average.t_i_average.value) / 1E3)
    ExtractLibFunction(:profiles, Symbol("<ne>"), "m⁻³", dd -> @ddtime(dd.summary.volume_average.n_e.value))
    ExtractLibFunction(:profiles, Symbol("<P>"), "MPa", dd -> begin
        cp1d = dd.core_profiles.profiles_1d[]
        integrate(cp1d.grid.volume, cp1d.pressure) / cp1d.grid.volume[end] / 1E6
    end)

    ExtractLibFunction(:profiles, Symbol("Te0/<Te>"), "-", dd -> EFL[:profiles_Te0](dd) / EFL[Symbol("profiles_<Te>")](dd))
    ExtractLibFunction(:profiles, Symbol("Ti0/<Ti>"), "-", dd -> EFL[:profiles_Ti0](dd) / EFL[Symbol("profiles_<Ti>")](dd))
    ExtractLibFunction(:profiles, Symbol("ne0/<ne>"), "-", dd -> EFL[:profiles_ne0](dd) / EFL[Symbol("profiles_<ne>")](dd))
    ExtractLibFunction(:profiles, Symbol("P0/<P>"), "-", dd -> EFL[:profiles_P0](dd) / EFL[Symbol("profiles_<P>")](dd))
    ExtractLibFunction(:profiles, :zeff, "-", dd -> @ddtime(dd.summary.volume_average.zeff.value))

    ExtractLibFunction(:transport, :τe, "s", dd -> @ddtime(dd.summary.global_quantities.tau_energy.value))
    ExtractLibFunction(:transport, :τe98, "s", dd -> @ddtime(dd.summary.global_quantities.tau_energy_98.value))
    ExtractLibFunction(:transport, :H98y2, "-", dd -> EFL[:transport_τe](dd) / EFL[:transport_τe98](dd))

    ExtractLibFunction(:hcd, :Pec, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ec.value) / 1E6)
    ExtractLibFunction(:hcd, :Pnbi, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_nbi.value) / 1E6)
    ExtractLibFunction(:hcd, :Pic, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ic.value) / 1E6)
    ExtractLibFunction(:hcd, :Plh, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_lh.value) / 1E6)
    ExtractLibFunction(:hcd, :Paux_tot, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_total.value) / 1E6)
    ExtractLibFunction(:hcd, :ip, "MA", dd -> (EFL[:hcd_ip_ni](dd) + EFL[:hcd_ip_aux](dd) + EFL[:hcd_ip_ohm](dd)))
    ExtractLibFunction(:hcd, :ip_ni, "MA", dd -> @ddtime(dd.summary.global_quantities.current_non_inductive.value) / 1E6)
    ExtractLibFunction(:hcd, :ip_bs, "MA", dd -> @ddtime(dd.summary.global_quantities.current_bootstrap.value) / 1E6)
    ExtractLibFunction(:hcd, :ip_aux, "MA", dd -> EFL[:hcd_ip_ni](dd) - EFL[:hcd_ip_bs](dd))
    ExtractLibFunction(:hcd, :ip_ohm, "MA", dd -> @ddtime(dd.summary.global_quantities.current_ohm.value) / 1E6)
    ExtractLibFunction(:hcd, :flattop, "Hours", dd -> dd.build.oh.flattop_duration / 3600.0)

    ExtractLibFunction(:bop, :Pfusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6)
    ExtractLibFunction(:bop, :Qfusion, "-", dd -> EFL[:bop_Pfusion](dd) / EFL[:hcd_Paux_tot](dd))
    ExtractLibFunction(:bop, :Pelectric_net, "MWe", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6)
    ExtractLibFunction(:bop, :Qplant, "-", dd -> @ddtime(dd.balance_of_plant.Q_plant))

    ExtractLibFunction(:costing, :levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE)
    ExtractLibFunction(:costing, :capital_cost, "\$B", dd -> dd.costing.cost_direct_capital.cost / 1E3)
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

"""
    extract(dd::IMAS.dd, xtract::Vector{ExtractFunction}=ExtractFunctionsLibrary)::Vector{ExtractFunction}

Extract data from `dd`. Each of the `ExtractFunction` should accept `dd` as input, like this:

    xtract = IMAS.ExtractFunction[
        :κ => ExtractFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
        :Te0 => ExtractFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ]

By default, the `ExtractFunctionsLibrary` is used.
"""
function extract(dd::IMAS.dd, xtract::AbstractVector{ExtractFunction}=values(ExtractFunctionsLibrary))::Vector{ExtractFunction}
    results = deepcopy(xtract)
    for xfun in xtract
        xfun(dd)
    end
    return results
end

function extract(dd::IMAS.dd, xtract::AbstractDict{Symbol,<:ExtractFunction}=ExtractFunctionsLibrary)::Vector{ExtractFunction}
    return extract(dd, collect(values(xtract)))
end

# ================= #
# show extract data #
# ================= #
function Base.show(io::IO, xfun::ExtractFunction; group::Bool=true, indent::Integer=0)
    printstyled(io, " "^indent; bold=true)
    if group
        printstyled(io, " "^indent * "$(xfun.group)."; bold=true)
    end
    printstyled(io, xfun.name; bold=true, color=:blue)
    printstyled(io, " → ")
    if any(isnan.(xfun.value))
        printstyled(io, xfun.value; color=:red)
    else
        printstyled(io, @sprintf("%.3g", xfun.value))
    end
    if xfun.units != "-"
        printstyled(io, " [$(xfun.units)]")
    end
end

function Base.show(io::IO, ::MIME"text/plain", xtract::Vector{ExtractFunction})
    last_group = ""
    for xfun in xtract
        if !isnan(xfun.value)
            if last_group != xfun.group
                if last_group != ""
                    printstyled(io, "\n"; bold=true)
                end
                printstyled(io, "$(xfun.group)\n"; bold=true)
                last_group = xfun.group
            end
            show(io, xfun; group=false, indent=4)
            println(io, "")
        end
    end
end

function Base.show(io::IO, x::MIME"text/plain", xtract::AbstractDict{Symbol,ExtractFunction})
    return Base.show(io, x, collect(values(xtract)))
end

function print_tiled(xtract::Vector{ExtractFunction}, terminal_width::Int=160)
    lists = OrderedCollections.OrderedDict{Symbol,Vector}()
    for xfun in xtract
        #if !isnan(xfun.value)
        group = xfun.group
        if group ∉ keys(lists)
            lists[group] = ExtractFunction[]
        end
        push!(lists[group], xfun)
        #end
    end

    if isempty(lists)
        return
    end

    function length_(xfun::ExtractFunction)
        buffer = IOBuffer()
        show(buffer, xfun; group=false)
        return length(String(take!(buffer)))
    end

    max_title_width = maximum([length(string(title)) for title in keys(lists)])
    max_item_width = maximum([maximum([length_(item) for item in list]) for list in values(lists)])
    max_width = max(max_title_width, max_item_width) + 4  # Add some padding
    max_height = maximum([length(list) for list in values(lists)]) + 2  # Add space for the title and separator

    ncols = max(1, floor(Int, terminal_width / max_width))
    nrows = ceil(Int, length(lists) / ncols)

    idx = 1
    for row in 1:nrows
        idxs = idx:min(idx + ncols - 1, length(lists))
        for title_row in idxs
            title = collect(keys(lists))[title_row]
            printstyled(title, " "^(max_width - length(string(title))); bold=true)
        end
        println()
        for list_row in 1:max_height
            for col in idxs
                list = collect(values(lists))[col]
                if list_row == 1
                    print(("─"^max_item_width * " "^(max_width - max_item_width)))
                elseif list_row - 1 <= length(list)
                    item = list[list_row-1]
                    show(stdout, item; group=false)
                    print(" "^(max_width - length_(item)))
                else
                    print(" "^max_width)
                end
            end
            println()
        end
        idx += ncols
    end
end