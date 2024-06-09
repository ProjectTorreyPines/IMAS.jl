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
    ExtractFunctionsLibrary[name] = xfun
    return xfun
end

const ExtractFunctionsLibrary = EFL = OrderedCollections.OrderedDict{Symbol,ExtractFunction}()

function update_ExtractFunctionsLibrary!()
    empty!(EFL)
    CFL = ConstraintFunctionsLibrary
    #! format: off
    ExtractLibFunction(:geometry, :R0, "m", dd -> dd.equilibrium.time_slice[].boundary.geometric_axis.r)
    ExtractLibFunction(:geometry, :a, "m", dd -> dd.equilibrium.time_slice[].boundary.minor_radius)
    ExtractLibFunction(:geometry, Symbol("1/ϵ"), "-", dd -> EFL[:R0](dd) / EFL[:a](dd))
    ExtractLibFunction(:geometry, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
    ExtractLibFunction(:geometry, :δ, "-", dd -> dd.equilibrium.time_slice[].boundary.triangularity)
    ExtractLibFunction(:geometry, :ζ, "-", dd -> dd.equilibrium.time_slice[].boundary.squareness)
    ExtractLibFunction(:geometry, :Volume, "m³", dd -> dd.equilibrium.time_slice[].profiles_1d.volume[end])
    ExtractLibFunction(:geometry, :Surface, "m²", dd -> dd.equilibrium.time_slice[].profiles_1d.surface[end])

    ExtractLibFunction(:equilibrium, :B0, "T", dd -> @ddtime(dd.equilibrium.vacuum_toroidal_field.b0))
    ExtractLibFunction(:equilibrium, :ip, "MA", dd -> @ddtime(dd.summary.global_quantities.ip.value) / 1e6)
    ExtractLibFunction(:equilibrium, :q95, "-", dd -> dd.equilibrium.time_slice[].global_quantities.q_95)
    ExtractLibFunction(:equilibrium, Symbol("<Bpol>"), "T", dd -> Bpol(EFL[:a](dd), EFL[:κ](dd), EFL[:ip](dd) * 1e6))
    ExtractLibFunction(:equilibrium, :βpol_MHD, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_pol)
    ExtractLibFunction(:equilibrium, :βtor_MHD, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_tor)
    ExtractLibFunction(:equilibrium, :βn_MHD, "-", dd -> dd.equilibrium.time_slice[].global_quantities.beta_normal)

    ExtractLibFunction(:temperatures, :Te0, "keV", dd -> @ddtime(dd.summary.local.magnetic_axis.t_e.value) / 1E3)
    ExtractLibFunction(:temperatures, :Ti0, "keV", dd -> @ddtime(dd.summary.local.magnetic_axis.t_i_average.value) / 1E3)
    ExtractLibFunction(:temperatures, Symbol("<Te>"), "keV", dd -> @ddtime(dd.summary.volume_average.t_e.value) / 1E3)
    ExtractLibFunction(:temperatures, Symbol("<Ti>"), "keV", dd -> @ddtime(dd.summary.volume_average.t_i_average.value) / 1E3)
    ExtractLibFunction(:temperatures, Symbol("Te0/<Te>"), "-", dd -> EFL[:Te0](dd) / EFL[Symbol("<Te>")](dd))
    ExtractLibFunction(:temperatures, Symbol("Ti0/<Ti>"), "-", dd -> EFL[:Ti0](dd) / EFL[Symbol("<Ti>")](dd))

    ExtractLibFunction(:densities, :ne0, "m⁻³", dd -> @ddtime(dd.summary.local.magnetic_axis.n_e.value))
    ExtractLibFunction(:densities, :ne_ped, "m⁻³", dd -> @ddtime(dd.summary.local.pedestal.n_e.value))
    ExtractLibFunction(:densities, Symbol("<ne>"), "m⁻³", dd -> @ddtime(dd.summary.volume_average.n_e.value))
    ExtractLibFunction(:densities, Symbol("ne0/<ne>"), "-", dd -> EFL[:ne0](dd) / EFL[Symbol("<ne>")](dd))
    ExtractLibFunction(:densities, Symbol("fGW"), "-", dd -> greenwald_fraction(dd))
    ExtractLibFunction(:densities, :zeff_ped, "-", dd -> @ddtime(dd.summary.local.pedestal.zeff.value))
    ExtractLibFunction(:densities, Symbol("<zeff>"), "-", dd -> @ddtime(dd.summary.volume_average.zeff.value))

    ExtractLibFunction(:pressures, :P0, "MPa", dd -> dd.core_profiles.profiles_1d[].pressure[1] / 1E6)
    ExtractLibFunction(:pressures, Symbol("<P>"), "MPa", dd -> begin
        cp1d = dd.core_profiles.profiles_1d[]
        trapz(cp1d.grid.volume, cp1d.pressure) / cp1d.grid.volume[end] / 1E6
    end)
    ExtractLibFunction(:pressures, Symbol("P0/<P>"), "-", dd -> EFL[:P0](dd) / EFL[Symbol("<P>")](dd))
    ExtractLibFunction(:pressures, :βn, "-", dd -> @ddtime(dd.summary.global_quantities.beta_tor_norm.value))
    ExtractLibFunction(:pressures, :βn_th, "-", dd -> @ddtime(dd.summary.global_quantities.beta_tor_thermal_norm.value))

    ExtractLibFunction(:transport, :τe, "s", dd -> @ddtime(dd.summary.global_quantities.tau_energy.value))
    ExtractLibFunction(:transport, :H98y2, "-", dd -> EFL[:τe](dd) / tau_e_h98(dd))
    ExtractLibFunction(:transport, :Hds03, "-", dd -> EFL[:τe](dd) / tau_e_ds03(dd))

    ExtractLibFunction(:sources, :Pec, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ec.value) / 1E6)
    ExtractLibFunction(:sources, :rho0_ec, "MW", dd -> findfirst(:ec, dd.core_sources.source).profiles_1d[].grid.rho_tor_norm[argmax(findfirst(:ec, dd.core_sources.source).profiles_1d[].electrons.energy)])
    ExtractLibFunction(:sources, :Pnbi, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_nbi.value) / 1E6)
    ExtractLibFunction(:sources, :Enbi1, "MeV", dd -> @ddtime(dd.nbi.unit[1].energy.data)/1e6)

    ExtractLibFunction(:sources, :Pic, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_ic.value) / 1E6)
    ExtractLibFunction(:sources, :Plh, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_lh.value) / 1E6)
    ExtractLibFunction(:sources, :Paux_tot, "MW", dd -> @ddtime(dd.summary.heating_current_drive.power_launched_total.value) / 1E6)
    ExtractLibFunction(:sources, :Pα, "MW", dd -> fusion_plasma_power(dd) / 1E6)
    ExtractLibFunction(:sources, :Pohm, "MW", dd -> total_power_source(ohmic_source!(dd).profiles_1d[])/ 1E6)
    ExtractLibFunction(:sources, :Pheat, "MW", dd ->  EFL[:Paux_tot](dd) + EFL[:Pα](dd) + EFL[:Pohm](dd))
    ExtractLibFunction(:sources, :Prad_tot, "MW", dd -> radiation_losses(dd.core_sources) / 1E6)

    ExtractLibFunction(:neutronics, :Nw_peak, "MW/m²", dd -> maximum(@. sqrt(dd.neutronics.time_slice[].wall_loading.flux_r^2+dd.neutronics.time_slice[].wall_loading.flux_z^2)) / 1E6)

    ExtractLibFunction(:exhaust, :Psol, "MW", dd -> power_sol(dd) / 1E6)
    ExtractLibFunction(:exhaust, :PLH, "MW", dd -> scaling_L_to_H_power(dd) / 1E6)
    ExtractLibFunction(:exhaust, :Bpol_omp, "T", dd -> Bpol_omp(dd.equilibrium.time_slice[]))
    ExtractLibFunction(:exhaust, :λq, "mm", dd -> widthSOL_eich(dd) * 1E3)
    ExtractLibFunction(:exhaust, :qpol, "MW/m²", dd -> q_pol_omp_eich(dd) / 1E6)
    ExtractLibFunction(:exhaust, :qpar, "MW/m²", dd -> q_par_omp_eich(dd) / 1E6)
    ExtractLibFunction(:exhaust, Symbol("P/R0"), "MW/m", dd -> EFL[:Psol](dd) / EFL[:R0](dd))
    ExtractLibFunction(:exhaust, Symbol("PB/R0"), "MW T/m", dd -> EFL[:Psol](dd) * EFL[:B0](dd) / EFL[:R0](dd))
    ExtractLibFunction(:exhaust, Symbol("PBp/R0"), "MW T/m", dd -> EFL[:Psol](dd) * EFL[Symbol("<Bpol>")](dd) / EFL[:R0](dd))
    ExtractLibFunction(:exhaust, Symbol("PBϵ/R0q95"), "MW T/m", dd -> zohm_divertor_figure_of_merit(dd)/1E6)

    ExtractLibFunction(:currents, :ip_bs_aux_ohm, "MA", dd -> (EFL[:ip_bs](dd) + EFL[:ip_aux](dd) + EFL[:ip_ohm](dd)))
    ExtractLibFunction(:currents, :ip_ni, "MA", dd -> @ddtime(dd.summary.global_quantities.current_non_inductive.value) / 1E6)
    ExtractLibFunction(:currents, :ip_bs, "MA", dd -> @ddtime(dd.summary.global_quantities.current_bootstrap.value) / 1E6)
    ExtractLibFunction(:currents, :ip_aux, "MA", dd -> EFL[:ip_ni](dd) - EFL[:ip_bs](dd))
    ExtractLibFunction(:currents, :ip_ohm, "MA", dd -> @ddtime(dd.summary.global_quantities.current_ohm.value) / 1E6)
    ExtractLibFunction(:currents, :ejima, "-", dd -> @ddtime(dd.core_profiles.global_quantities.ejima))
    ExtractLibFunction(:currents, :flattop, "Hours", dd -> max(0.0, dd.build.oh.flattop_duration / 3600.0))

    ExtractLibFunction(:bop, :Pfusion, "MW", dd -> fusion_power(dd.core_profiles.profiles_1d[]) / 1E6)
    ExtractLibFunction(:bop, :Qfusion, "-", dd -> EFL[:Pfusion](dd) / EFL[:Paux_tot](dd))
    ExtractLibFunction(:bop, :thermal_cycle_type, "-", dd -> dd.balance_of_plant.power_plant.power_cycle_type)
    ExtractLibFunction(:bop, :thermal_efficiency_plant, "MW", dd -> @ddtime(dd.balance_of_plant.thermal_efficiency_plant) )
    ExtractLibFunction(:bop, :thermal_efficiency_cycle, "MW", dd -> @ddtime(dd.balance_of_plant.thermal_efficiency_cycle) )
    ExtractLibFunction(:bop, :power_electric_generated, "MW", dd -> @ddtime(dd.balance_of_plant.power_plant.power_electric_generated) / 1E6)
    ExtractLibFunction(:bop, :Pelectric_net, "MW", dd -> @ddtime(dd.balance_of_plant.power_electric_net) / 1E6)
    ExtractLibFunction(:bop, :Qplant, "-", dd -> @ddtime(dd.balance_of_plant.Q_plant))
    ExtractLibFunction(:bop, :TBR, "-", dd -> @ddtime(dd.blanket.tritium_breeding_ratio))

    ExtractLibFunction(:build, :PF_material, "-", dd -> dd.build.pf_active.technology.material)
    ExtractLibFunction(:build, :TF_material, "-", dd -> dd.build.tf.technology.material)
    ExtractLibFunction(:build, :OH_material, "-", dd -> dd.build.oh.technology.material)
    ExtractLibFunction(:build, :TF_max_b, "T", dd -> dd.build.tf.max_b_field)
    ExtractLibFunction(:build, :OH_max_b, "T", dd -> dd.build.oh.max_b_field)
    ExtractLibFunction(:build, :TF_j_margin, "-", dd -> dd.build.tf.critical_j/dd.build.tf.max_j)
    ExtractLibFunction(:build, :OH_j_margin, "-", dd -> dd.build.oh.critical_j/dd.build.oh.max_j)
    ExtractLibFunction(:build, :TF_stress_margin, "-", dd -> dd.solid_mechanics.center_stack.properties.yield_strength.tf/maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf))
    ExtractLibFunction(:build, :OH_stress_margin, "-", dd -> dd.solid_mechanics.center_stack.properties.yield_strength.oh/maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh))

    ExtractLibFunction(:costing, :levelized_CoE, "\$/kWh", dd -> dd.costing.levelized_CoE)
    ExtractLibFunction(:costing, :TF_of_total, "%", dd -> 100 * select_direct_captial_cost(dd,"TF") / dd.costing.cost_direct_capital.cost)
    ExtractLibFunction(:costing, :BOP_of_total, "%", dd -> 100 * select_direct_captial_cost(dd,"balance of plant equipment") / dd.costing.cost_direct_capital.cost)
    ExtractLibFunction(:costing, :blanket_of_total, "%", dd -> 100 * select_direct_captial_cost(dd,"blanket") / dd.costing.cost_direct_capital.cost)
    ExtractLibFunction(:costing, :cryostat_of_total, "%", dd -> 100 * select_direct_captial_cost(dd,"cryostat") / dd.costing.cost_direct_capital.cost)

    ExtractLibFunction(:costing, :capital_cost, "\$B", dd -> dd.costing.cost_direct_capital.cost / 1E3)
    ExtractLibFunction(:constraint, :min_required_power_electric_net, "-", dd -> CFL[:min_required_power_electric_net](dd))
    ExtractLibFunction(:constraint, :required_power_electric_net, "-", dd -> CFL[:required_power_electric_net](dd))
    ExtractLibFunction(:constraint, :min_q95, "-", dd -> CFL[:min_q95](dd))
    ExtractLibFunction(:constraint, :max_tf_j, "-", dd -> CFL[:max_tf_j](dd))
    ExtractLibFunction(:constraint, :max_oh_j, "-", dd -> CFL[:max_oh_j](dd))
    ExtractLibFunction(:constraint, :max_pl_stress, "-", dd -> CFL[:max_pl_stress](dd))
    ExtractLibFunction(:constraint, :max_tf_stress, "-", dd -> CFL[:max_tf_stress](dd))
    ExtractLibFunction(:constraint, :max_oh_stress, "-", dd -> CFL[:max_oh_stress](dd))

    ExtractLibFunction(:constraint, :min_required_power_electric_net, "-", dd -> CFL[:min_required_power_electric_net](dd))
    ExtractLibFunction(:constraint, :required_power_electric_net, "-", dd -> CFL[:required_power_electric_net](dd))
    ExtractLibFunction(:constraint, :min_q95, "-", dd -> CFL[:min_q95](dd))
    ExtractLibFunction(:constraint, :max_tf_j, "-", dd -> CFL[:max_tf_j](dd))
    ExtractLibFunction(:constraint, :max_oh_j, "-", dd -> CFL[:max_oh_j](dd))
    ExtractLibFunction(:constraint, :max_pl_stress, "-", dd -> CFL[:max_pl_stress](dd))
    ExtractLibFunction(:constraint, :max_tf_stress, "-", dd -> CFL[:max_tf_stress](dd))
    ExtractLibFunction(:constraint, :max_oh_stress, "-", dd -> CFL[:max_oh_stress](dd))
    #! format: on

    return EFL
end
update_ExtractFunctionsLibrary!()

"""
    (xfun::ExtractFunction)(dd::IMAS.dd)

Run the extract function and store its result in xfun.value

NOTE: NaN is assigned on error
"""
function (xfun::ExtractFunction)(dd::IMAS.dd)
    try
        return xfun.func(dd)
    catch e
        return NaN
    end
end

"""
    extract(dd::IMAS.dd, xtract::AbstractDict{Symbol,<:ExtractFunction}=ExtractFunctionsLibrary)

Extract data from `dd`. Each of the `ExtractFunction` should accept `dd` as input, like this:

    xtract = IMAS.ExtractFunction[
        :κ => ExtractFunction(:equilibrium, :κ, "-", dd -> dd.equilibrium.time_slice[].boundary.elongation)
        :Te0 => ExtractFunction(:profiles, :Te0, "keV", dd -> dd.core_profiles.profiles_1d[].electrons.temperature[1] / 1E3)
    ]

By default, the `ExtractFunctionsLibrary` is used.
"""
function extract(dd::IMAS.dd, xtract::AbstractDict{Symbol,<:ExtractFunction}=ExtractFunctionsLibrary)
    xtract_out = OrderedCollections.OrderedDict{Symbol,ExtractFunction}()
    for xfun in values(xtract)
        xtract_out[xfun.name] = deepcopy(xfun)
        xtract_out[xfun.name].value = xfun(dd)
    end
    return xtract_out
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
    if !(typeof(xfun.value) <: Real)
        printstyled(io, xfun.value)
    elseif any(isnan.(xfun.value))
        printstyled(io, xfun.value; color=:red)
    else
        printstyled(io, @sprintf("%.3g", xfun.value))
    end
    if xfun.units != "-"
        printstyled(io, " [$(xfun.units)]")
    end
end

function Base.show(io::IO, x::MIME"text/plain", xtract::AbstractDict{Symbol,ExtractFunction}; terminal_width::Int=136)
    return print_tiled(io, xtract; terminal_width)
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

function print_tiled(xtract::AbstractDict{Symbol,ExtractFunction}; terminal_width::Int, line_char::Char='─')
    return print_tiled(stdout, xtract; terminal_width, line_char)
end

function print_tiled(io::IO, xtract::AbstractDict{Symbol,ExtractFunction}; terminal_width::Int, line_char::Char='─')
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
            title = uppercase(string(collect(keys(lists))[title_row]))
            printstyled(io, rpad(title, max_width); bold=true)
        end
        println(io)
        for list_row in 1:max_heights[row]+2
            for col in idxs
                list = collect(values(lists))[col]
                if list_row == 1
                    print(io, (line_char^max_item_width * " "^(max_width - max_item_width)))
                elseif list_row - 1 <= length(list)
                    item = list[list_row-1]
                    show(io, item; group=false)
                    print(io, " "^(max_width - length_(item)))
                else
                    print(io, " "^max_width)
                end
            end
            println(io)
        end
        idx += ncols
    end
end

function select_direct_captial_cost(dd::IMAS.dd, what::String)
    for sys in dd.costing.cost_direct_capital.system
        idx = findfirst(x-> x.name ==what, sys.subsystem)
        if !isnothing(idx)
            return sys.subsystem[idx].cost
        end
    end
end
