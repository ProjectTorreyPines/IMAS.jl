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
    ExtractFunctionsLibrary[Symbol(name)] = xfun
    return xfun
end

const ExtractFunctionsLibrary = EFL = OrderedCollections.OrderedDict{Symbol,ExtractFunction}()

function update_ExtractFunctionsLibrary!()
    empty!(EFL)
    ExtractLibFunction(:geometry, :R0, "m", dd -> dd.equilibrium.time_slice[].boundary.geometric_axis.r)
    ExtractLibFunction(:geometry, :a, "m", dd -> dd.equilibrium.time_slice[].boundary.minor_radius)
    ExtractLibFunction(:geometry, Symbol("1/ϵ"), "m", dd -> EFL[:R0](dd) / EFL[:a](dd))
    ExtractLibFunction(:geometry, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
    ExtractLibFunction(:geometry, :δ, "-", dd -> dd.equilibrium.time_slice[].boundary.triangularity)
    ExtractLibFunction(:geometry, :ζ, "-", dd -> dd.equilibrium.time_slice[].boundary.squareness)

    ExtractLibFunction(:equilibrium, :B0, "T", dd -> @ddtime(dd.summary.global_quantities.b0.value))
    ExtractLibFunction(:equilibrium, :ip, "MA", dd -> @ddtime(dd.summary.global_quantities.ip.value) / 1e6)
    ExtractLibFunction(:equilibrium, :q95, "-", dd -> dd.equilibrium.time_slice[].global_quantities.q_95)
    ExtractLibFunction(:equilibrium, :βpol, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_pol)
    ExtractLibFunction(:equilibrium, :βtor, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_tor)
    ExtractLibFunction(:equilibrium, :βn, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal)

    ExtractLibFunction(:temperatures, :Te0, "keV", dd -> @ddtime(dd.summary.local.magnetic_axis.t_e.value) / 1E3)
    ExtractLibFunction(:temperatures, :Ti0, "keV", dd -> @ddtime(dd.summary.local.magnetic_axis.t_i_average.value) / 1E3)
    ExtractLibFunction(:temperatures, Symbol("<Te>"), "keV", dd -> @ddtime(dd.summary.volume_average.t_e.value) / 1E3)
    ExtractLibFunction(:temperatures, Symbol("<Ti>"), "keV", dd -> @ddtime(dd.summary.volume_average.t_i_average.value) / 1E3)
    ExtractLibFunction(:temperatures, Symbol("Te0/<Te>"), "-", dd -> EFL[:Te0](dd) / EFL[Symbol("<Te>")](dd))
    ExtractLibFunction(:temperatures, Symbol("Ti0/<Ti>"), "-", dd -> EFL[:Ti0](dd) / EFL[Symbol("<Ti>")](dd))

    ExtractLibFunction(:densities, :ne0, "m⁻³", dd -> @ddtime(dd.summary.local.magnetic_axis.n_e.value))
    ExtractLibFunction(:densities, Symbol("<ne>"), "m⁻³", dd -> @ddtime(dd.summary.volume_average.n_e.value))
    ExtractLibFunction(:densities, Symbol("ne0/<ne>"), "-", dd -> EFL[:ne0](dd) / EFL[Symbol("<ne>")](dd))
    ExtractLibFunction(:densities, :zeff, "-", dd -> @ddtime(dd.summary.volume_average.zeff.value))

    ExtractLibFunction(:pressures, :P0, "MPa", dd -> dd.core_profiles.profiles_1d[].pressure[1] / 1E6)
    ExtractLibFunction(:pressures, Symbol("<P>"), "MPa", dd -> begin
        cp1d = dd.core_profiles.profiles_1d[]
        integrate(cp1d.grid.volume, cp1d.pressure) / cp1d.grid.volume[end] / 1E6
    end)
    ExtractLibFunction(:pressures, Symbol("P0/<P>"), "-", dd -> EFL[:P0](dd) / EFL[Symbol("<P>")](dd))

    ExtractLibFunction(:transport, :τe, "s", dd -> @ddtime(dd.summary.global_quantities.tau_energy.value))
    ExtractLibFunction(:transport, :τe98, "s", dd -> @ddtime(dd.summary.global_quantities.tau_energy_98.value))
    ExtractLibFunction(:transport, :H98y2, "-", dd -> EFL[:τe](dd) / EFL[:τe98](dd))

    ExtractLibFunction(:sources, :Pec, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ec.value) / 1E6)
    ExtractLibFunction(:sources, :Pnbi, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_nbi.value) / 1E6)
    ExtractLibFunction(:sources, :Pic, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ic.value) / 1E6)
    ExtractLibFunction(:sources, :Plh, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_lh.value) / 1E6)
    ExtractLibFunction(:sources, :Paux_tot, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_total.value) / 1E6)

    ExtractLibFunction(:currents, :ip_bs_aux_ohm, "MA", dd -> (EFL[:ip_bs](dd) + EFL[:ip_aux](dd) + EFL[:ip_ohm](dd)))
    ExtractLibFunction(:currents, :ip_ni, "MA", dd -> @ddtime(dd.summary.global_quantities.current_non_inductive.value) / 1E6)
    ExtractLibFunction(:currents, :ip_bs, "MA", dd -> @ddtime(dd.summary.global_quantities.current_bootstrap.value) / 1E6)
    ExtractLibFunction(:currents, :ip_aux, "MA", dd -> EFL[:ip_ni](dd) - EFL[:ip_bs](dd))
    ExtractLibFunction(:currents, :ip_ohm, "MA", dd -> @ddtime(dd.summary.global_quantities.current_ohm.value) / 1E6)
    ExtractLibFunction(:currents, :flattop, "Hours", dd -> dd.build.oh.flattop_duration / 3600.0)

    ExtractLibFunction(:bop, :Pfusion, "MW", dd -> IMAS.fusion_power(dd.core_profiles) / 1E6)
    ExtractLibFunction(:bop, :Qfusion, "-", dd -> EFL[:Pfusion](dd) / EFL[:Paux_tot](dd))
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
    catch e
        #display(e)
        xfun.value = NaN
    end
    return xfun.value
end

"""
    extract(dd::IMAS.dd, xtract::AbstractDict{Symbol,<:ExtractFunction}=ExtractFunctionsLibrary)::Vector{ExtractFunction}

Extract data from `dd`. Each of the `ExtractFunction` should accept `dd` as input, like this:

    xtract = IMAS.ExtractFunction[
        :κ => ExtractFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
        :Te0 => ExtractFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ]

By default, the `ExtractFunctionsLibrary` is used.
"""
function extract(dd::IMAS.dd, xtract::T=ExtractFunctionsLibrary)::T where {T<:AbstractDict{Symbol,<:ExtractFunction}}
    xtract = deepcopy(xtract)
    for xfun in values(xtract)
        xfun(dd)
    end
    return xtract
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

function Base.show(io::IO, x::MIME"text/plain", xtract::AbstractDict{Symbol,ExtractFunction})
    return print_tiled(io, xtract)
end

function print_vertical(io, xtract::AbstractDict{Symbol,ExtractFunction})
    last_group = ""
    for xfun in values(xtract)
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

function print_tiled(xtract::AbstractDict{Symbol,ExtractFunction}, terminal_width::Int=160)
    return print_tiled(stdout, xtract, terminal_width)
end

function print_tiled(io::IO, xtract::AbstractDict{Symbol,ExtractFunction}, terminal_width::Int=160)
    lists = OrderedCollections.OrderedDict{Symbol,Vector}()
    for xfun in values(xtract)
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

    ncols = max(1, floor(Int, terminal_width / max_width))
    nrows = ceil(Int, length(lists) / ncols)

    idx = 1
    max_heights = []
    for row in 1:nrows
        idxs = idx:min(idx + ncols - 1, length(lists))
        max_height = 0
        for col in idxs
            list = collect(values(lists))[col]
            max_height = max(max_height, length(list))
        end
        push!(max_heights, max_height)
        idx += ncols
    end

    idx = 1
    for row in 1:nrows
        idxs = idx:min(idx + ncols - 1, length(lists))
        for title_row in idxs
            title = collect(keys(lists))[title_row]
            printstyled(io,title, " "^(max_width - length(string(title))); bold=true)
        end
        println(io)
        for list_row in 1:max_heights[row]+2
            for col in idxs
                list = collect(values(lists))[col]
                if list_row == 1
                    print(io,("─"^max_item_width * " "^(max_width - max_item_width)))
                elseif list_row - 1 <= length(list)
                    item = list[list_row-1]
                    show(io, item; group=false)
                    print(io," "^(max_width - length_(item)))
                else
                    print(io, " "^max_width)
                end
            end
            println(io)
        end
        idx += ncols
    end
end