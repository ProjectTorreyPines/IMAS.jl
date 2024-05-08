mutable struct ConstraintFunction
    name::Symbol
    units::String
    func::Function
    operation::Function
    limit::Float64
    tolerance::Float64
    # inner constructor to register ConstraintFunction in ConstraintsFunctionsLibrary
    ConstraintFunction(name::Symbol, units::String, func::Function, operation::Function, limit::Float64, tolerance::Float64) = begin
        @assert ===(operation, ==) "tolerance specification only used for == constraint"
        cnst = new(name, units, func, operation, limit, tolerance)
        ConstraintFunctionsLibrary[cnst.name] = cnst
        return cnst
    end
    ConstraintFunction(name::Symbol, units::String, func::Function, operation::Function, limit::Float64) = begin
        @assert !==(operation, ==) "Must specify tolerance of == constraint"
        cnst = new(name, units, func, operation, limit, 0.0)
        ConstraintFunctionsLibrary[cnst.name] = cnst
        return cnst
    end
end



const ConstraintFunctionsLibrary = Dict{Symbol,ConstraintFunction}() #s
function update_ConstraintFunctionsLibrary!()
    empty!(ConstraintFunctionsLibrary)
    #! format: off
    # logic for creating constraint functions with margin
    # we want: value/max_value + margin < 1
    ConstraintFunction(:required_power_electric_net, "%", dd -> abs(@ddtime(dd.balance_of_plant.power_electric_net) - dd.requirements.power_electric_net) / dd.requirements.power_electric_net, ==, 0.0, 0.2) # relative tolerance
    ConstraintFunction(:min_required_power_electric_net, "%", dd -> (@ddtime(dd.balance_of_plant.power_electric_net) - dd.requirements.power_electric_net) / dd.requirements.power_electric_net, >, 0.0)
    ConstraintFunction(:required_flattop, "%", dd -> abs(dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / dd.requirements.flattop_duration, ==, 0.0, 0.01) # relative tolerance
    ConstraintFunction(:min_required_flattop, "%", dd -> (dd.build.oh.flattop_duration - dd.requirements.flattop_duration) / dd.requirements.flattop_duration, >, 0.0)
    ConstraintFunction(:min_required_B0, "%", dd -> (abs(prod(IMAS.build_max_R0_B0(dd.build))) - dd.equilibrium.vacuum_toroidal_field.r0 * maximum(abs, dd.equilibrium.vacuum_toroidal_field.b0)) / (dd.equilibrium.vacuum_toroidal_field.r0 * maximum(abs, dd.equilibrium.vacuum_toroidal_field.b0)), >, 0.0)
    ConstraintFunction(:zero_ohmic, "MA", dd -> abs(sum(integrate(dd.core_profiles.profiles_1d[].grid.area, dd.core_profiles.profiles_1d[].j_ohmic))) / 1E6, ==, 0.0, 0.1) # absolute tolerance
    ConstraintFunction(:max_ne_peaking, "%", dd -> ((@ddtime(dd.summary.local.magnetic_axis.n_e.value) / @ddtime(dd.summary.volume_average.n_e.value)) - dd.requirements.ne_peaking) / dd.requirements.ne_peaking, <, 0.0)
    ConstraintFunction(:min_lh_power_threshold, "%", dd -> (IMAS.power_sol(dd) / dd.requirements.lh_power_threshold_fraction - IMAS.scaling_L_to_H_power(dd)) / IMAS.scaling_L_to_H_power(dd), >, 0.0)
    ConstraintFunction(:max_ωpe_ωce, "%", dd -> IMAS.ω_pe(@ddtime(dd.summary.local.magnetic_axis.n_e.value)) / IMAS.ω_ce(@ddtime(dd.equilibrium.vacuum_toroidal_field.b0)), <, 1.0)
    ConstraintFunction(:max_qpol_omp, "%", dd -> (IMAS.q_pol_omp_eich(dd) - dd.requirements.q_pol_omp) / dd.requirements.q_pol_omp, <, 0.0)
    ConstraintFunction(:max_tf_j, "%", dd -> dd.build.tf.max_j / dd.build.tf.critical_j + dd.requirements.coil_j_margin, <, 1.0)
    ConstraintFunction(:max_oh_j, "%", dd -> dd.build.oh.max_j / dd.build.oh.critical_j + dd.requirements.coil_j_margin, <, 1.0)
    ConstraintFunction(:max_pl_stress, "%", dd -> maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl) / (ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl) ? 0.0 : dd.solid_mechanics.center_stack.properties.yield_strength.pl) + dd.requirements.coil_stress_margin, <, 1.0)
    ConstraintFunction(:max_tf_stress, "%", dd -> maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) / dd.solid_mechanics.center_stack.properties.yield_strength.tf + dd.requirements.coil_stress_margin, <, 1.0)
    ConstraintFunction(:max_oh_stress, "%", dd -> maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) / dd.solid_mechanics.center_stack.properties.yield_strength.oh + dd.requirements.coil_stress_margin, <, 1.0)
    ConstraintFunction(:max_hds03, "%", dd -> (@ddtime(dd.summary.global_quantities.tau_energy.value)/IMAS.tau_e_ds03(dd) - dd.requirements.hds03)/dd.requirements.hds03, <, 0.0)
    ConstraintFunction(:min_q95, "%", dd -> (dd.equilibrium.time_slice[].global_quantities.q_95 - dd.requirements.q95) / dd.requirements.q95, >, 0.0)
    #! format: on
    return ConstraintFunctionsLibrary
end
update_ConstraintFunctionsLibrary!()

function (cnst::ConstraintFunction)(dd::IMAS.dd)
    return constraint_cost_transform(cnst.func(dd), cnst.operation, cnst.limit, cnst.tolerance)
end

function constraint_cost_transform(value::Float64, operation::Function, limit::Float64, tolerance::Float64)
    if ===(operation, ==)
        out = abs(value - limit) - tolerance
    elseif operation(1.0, 0.0) # > or >=
        out = limit - value
    else # < or <=
        out = value - limit
    end

    if out < 0.0
        out = 0.0
    end

    return out
end

function Base.show(io::IO, cnst::ConstraintFunction)
    printstyled(io, cnst.name; bold=true, color=:blue)
    print(io, " $(cnst.operation)")
    print(io, " $(cnst.limit)")
    if ===(cnst.operation, ==)
        print(io, " ± $(cnst.tolerance)")
    end
    return print(io, " [$(cnst.units)]")
end

function Base.show(io::IO, x::MIME"text/plain", cnsts::AbstractDict{Symbol,ConstraintFunction})
    for cnst in cnsts
        show(io, x, cnst)
        println(io, "")
    end
end

update_ConstraintFunctionsLibrary!()