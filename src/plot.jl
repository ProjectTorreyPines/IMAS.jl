import PlotUtils
using RecipesBase
using LaTeXStrings
import Measures
import Graphs
using GraphRecipes

# ========= #
# plot_help #
# ========= #
mutable struct PlotHelpParameter
    dispatch::String
    argument::Symbol
    type::Any
    description::String
    value::Any
end

const _help_plot = Vector{PlotHelpParameter}()

function plot_help_id(args...)
    dispatches = String[]
    for arg in args
        if typeof(arg) <: IDS
            push!(dispatches, ulocation(arg))
        else
            push!(dispatches, string(typeof(arg)))
        end
    end
    return join(dispatches, ", ")
end

function assert_type_and_record_argument(dispatch, type::Type, description::String; kw...)
    @assert length(kw) == 1
    argument = collect(keys(kw))[1]
    value = collect(values(kw))[1]

    this_plotpar = PlotHelpParameter(dispatch, argument, type, description, value)

    only_match = true
    all_same_values = all(plotpar.value == this_plotpar.value for plotpar in _help_plot if "$(plotpar.dispatch), $(plotpar.argument)" == "$(dispatch), $(argument)")
    for plotpar in _help_plot
        if "$(plotpar.dispatch), $(plotpar.argument)" == "$(this_plotpar.dispatch), $(this_plotpar.argument)"
            if !all_same_values
                plotpar.value = :__MIXED__
            end
            only_match = false
        end
    end

    if only_match
        push!(_help_plot, this_plotpar)
    end

    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", plotpar::PlotHelpParameter)
    printstyled(io, plotpar.argument; bold=true)
    printstyled(io, "::$(plotpar.type)"; color=:blue)
    printstyled(io, " = ")
    if plotpar.value == :__MIXED__
        printstyled(io, "MIXED"; color=:magenta, bold=true)
    else
        printstyled(io, "$(repr(plotpar.value))"; color=:red, bold=true)
    end
    return printstyled(io, "  # $(plotpar.description)"; color=248)
end

function Base.show(io::IO, x::MIME"text/plain", plotpars::AbstractVector{<:PlotHelpParameter})
    old_dispatch = nothing
    for (k, plotpar) in enumerate(plotpars)
        if plotpar.dispatch != old_dispatch
            if k > 1
                print(io, "\n")
            end
            printstyled(io, "$(plotpar.dispatch)"; color=:green)
            print(io, "\n")
        end
        old_dispatch = plotpar.dispatch
        print(io, "     ")
        show(io, x, plotpar)
        if k < length(plotpars)
            print(io, "\n")
        end
    end
end

"""
    help_plot(args...; kw...)

Prints plotting arguments for object that are part of IMAS.

Call `help_plot(...)` just like you would `plot(...)`

Returns the plot.
"""
function help_plot(args...; kw...)
    empty!(_help_plot)
    p = plot(args...; kw...)
    display(_help_plot)
    return p
end

export help_plot
push!(document[:Plot], :help_plot)

"""
    help_plot!(args...; kw...)

Prints plotting arguments for object that are part of IMAS.

Call `help_plot!(...)` just like you would `plot!(...)`

Returns the plot.
"""
function help_plot!(args...; kw...)
    empty!(_help_plot)
    p = plot!(args...; kw...)
    display(_help_plot)
    return p
end

export help_plot!
push!(document[:Plot], :help_plot!)

# ========= #
# pf_active #
# ========= #
@recipe function plot_pf_active_cx(pfa::pf_active, what::Symbol=:cx; time0=global_time(pfa), cname=:vik)
    id = plot_help_id(pfa, what)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Symbol, "Colormap name"; cname)

    if ismissing(pfa.coil[1].current, :time) || isempty(pfa.coil[1].current.time) || time0 < pfa.coil[1].current.time[1]
        currents = [0.0 for c in pfa.coil]
        c_unit = "A"
    else
        if time0 == -Inf && pfa.coil[1].current.time[1] == -Inf
            index = 1
        elseif pfa.coil[1].current.time[1] == -Inf
            index = 2:length(pfa.coil[1].current.time)
        else
            index = 1:length(pfa.coil[1].current.time)
        end

        currents = [get_time_array(c.current, :data, time0) * getproperty(c.element[1], :turns_with_sign, 1.0) for c in pfa.coil]
        CURRENT = maximum((maximum(abs, c.current.data[index] * getproperty(c.element[1], :turns_with_sign, 1.0)) for c in pfa.coil))
        if maximum(currents) > 1e6
            currents = currents ./ 1e6
            CURRENT = CURRENT ./ 1e6
            c_unit = "MA"
        elseif maximum(currents) > 1e3
            currents = currents ./ 1e3
            CURRENT = CURRENT ./ 1e3
            c_unit = "kA"
        else
            c_unit = "A"
        end
    end

    if what ∈ (:cx, :coils_flux)
        label --> ""
        aspect_ratio := :equal

        # dummy markers to get the colorbar right
        if any(currents .!= 0.0)
            @series begin
                colorbar_title := "PF currents [$c_unit]"
                seriestype --> :scatter
                color --> cname
                clim --> (-CURRENT, CURRENT)
                marker_z --> [-CURRENT, CURRENT]
                [(NaN, NaN), (NaN, NaN)]
            end
        end

        # plot individual coils
        for (k, c) in enumerate(pfa.coil)
            @series begin
                colorbar_entry := false
                if all(currents .== 0.0)
                    color --> hash_to_color(sort!([func.index for func in c.function]); seed=4)
                else
                    current_color_index = (currents[k] + CURRENT) / (2 * CURRENT)
                    color --> PlotUtils.cgrad(cname)[current_color_index]
                end
                c
            end
        end

    elseif what == :currents
        label --> "$(time0) s"
        @series begin
            linestyle --> :dash
            marker --> :circle
            ylabel := "[$c_unit]"
            ["$k" for k in eachindex(currents)], currents
        end

        Imax = []
        for c in pfa.coil
            if !ismissing(c.b_field_max_timed, :data)
                b_max = get_time_array(c.b_field_max_timed, :data, time0)
                # issue: IMAS does not have a way to store the current pf coil temperature
                #temperature = c.temperature[1]
                #Icrit = Interpolations.cubic_spline_interpolation((to_range(c.b_field_max), to_range(c.temperature)), c.current_limit_max * c.element[1].turns_with_sign)(b_max, temperature)
                Icrit = interp1d(c.b_field_max, c.current_limit_max[:, 1] * getproperty(c.element[1], :turns_with_sign, 1.0))(b_max)
                push!(Imax, Icrit / 1E6)
            else
                push!(Imax, NaN)
            end
        end

        if !all(isnan.(Imax))
            @series begin
                marker --> :cross
                label := "Max current"
                ylabel := "[$c_unit]"
                ["$k" for k in eachindex(currents)], Imax
            end
            @series begin
                marker --> :cross
                primary := false
                ylabel := "[$c_unit]"
                ["$k" for k in eachindex(currents)], -Imax
            end
        end

    else
        error("IMAS.pf_active `what` to plot can only be :cx or :currents")
    end

end

@recipe function plot_coil(coil::pf_active__coil{T}; coil_names=false, coil_identifiers=false, coil_numbers=false) where {T<:Real}
    id = plot_help_id(coil)
    assert_type_and_record_argument(id, Bool, "Show coil names"; coil_names)
    assert_type_and_record_argument(id, Bool, "Show coil identifiers"; coil_identifiers)
    assert_type_and_record_argument(id, Bool, "Show coil numbers"; coil_numbers)

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    r = T[]
    z = T[]
    for (k, element) in enumerate(coil.element)
        oute = outline(element)
        @series begin
            primary := k == 1
            seriestype --> :shape
            linewidth := 0.25 * base_linewidth
            colorbar --> :right
            label --> ""
            oute.r, oute.z
        end
        append!(r, oute.r)
        append!(z, oute.z)
    end

    if coil_names || coil_identifiers || coil_numbers
        r_avg = sum(r) / length(r)
        z_avg = sum(z) / length(z)
        @series begin
            label := ""
            if coil_names
                series_annotations := [(coil.name, :center, :middle, :red, 6)]
            elseif coil_numbers
                series_annotations := [(index(coil), :center, :middle, :red, 6)]
            else
                series_annotations := [(coil.identifier, :center, :middle, :red, 6)]
            end
            [r_avg], [z_avg]
        end
    end
end

@recipe function plot_coil(oute::pf_active__coil___element___geometry__outline{T}) where {T<:Real}
    @series begin
        seriestype --> :shape
        oute.r, oute.z
    end
end

@recipe function plot_pf_active_rail(rail::IMAS.build__pf_active__rail)
    if !ismissing(rail.outline, :r)
        @series begin
            seriestype --> :scatter
            color --> :gray
            marker --> :circle
            markerstrokewidth --> 0
            rail.outline.r, rail.outline.z
        end
    end
end

@recipe function plot_pf_active_rail(rails::IDSvector{<:IMAS.build__pf_active__rail})
    for (krail, rail) in enumerate(rails)
        if !ismissing(rail.outline, :r)
            @series begin
                label --> "Coil opt. rail"
                primary --> krail == 1 ? true : false
                rail
            end
        end
    end
end

@recipe function plot_circuit(circuit::IMAS.pf_active__circuit)
    pf_active = top_ids(circuit)
    supplies_labels = [supply.identifier for supply in pf_active.supply]
    coils_labels = [coil.identifier for coil in pf_active.coil]
    @series begin
        (circuit, supplies_labels, coils_labels)
    end
end

@recipe function plot_circuit(circuit::IMAS.pf_active__circuit, supplies_labels::Vector{String}, coils_labels::Vector{String})
    title_str = circuit.name
    conn_matrix = circuit.connections
    num_supplies = length(supplies_labels)
    num_coils = length(coils_labels)

    conn_matrix = conn_matrix[:, 1:2:end] .+ conn_matrix[:, 2:2:end]
    num_nodes, num_elements = size(conn_matrix)
    @assert num_elements == num_supplies + num_coils "num_elements=$(num_elements) while (num_supplies + num_coils) * 2 = $((num_supplies + num_coils) * 2)"

    membership = [fill(2, num_nodes); fill(1, num_supplies); fill(3, num_coils)]

    names = []
    for k in 1:num_nodes
        push!(names, " $k ")
    end
    if isempty(supplies_labels)
        for k in 1:num_supplies
            push!(names, "S$k")
        end
    else
        append!(names, supplies_labels)
    end
    if isempty(coils_labels)
        for k in 1:num_coils
            push!(names, "C$k")
        end
    else
        append!(names, coils_labels)
    end
    names = collect(map(x -> replace(x, " " => "\n"), names))

    # hide unused supplies and coils
    index = (1:size(conn_matrix)[2])[dropdims(sum(conn_matrix; dims=1); dims=1).>0]
    conn_matrix = conn_matrix[:, index]
    membership = [membership[1:num_nodes]; membership[num_nodes.+index]]
    names = [names[1:num_nodes]; names[num_nodes.+index]]
    num_supplies = sum(membership .== 1)
    num_coils = sum(membership .== 3)

    # Pastel versions of red, green, and blue
    markercolor = [PlotUtils.Colors.colorant"#FFAAAA", PlotUtils.Colors.colorant"#AAFFAA", PlotUtils.Colors.colorant"#AAAAFF"]
    markercolor = markercolor[membership]
    node_weights = [1.0, 0.25, 1.0]
    node_weights = node_weights[membership]
    nodeshape = [:rect, :circle, :rect]
    nodeshape = nodeshape[membership]

    # Create an empty graph with nodes
    g = Graphs.SimpleDiGraph(num_nodes + num_supplies + num_coils)

    # Iterate through the matrix to add edges
    for i in 1:num_nodes
        connected_elements = findall(x -> x == 1, conn_matrix[i, :])
        for j in connected_elements
            Graphs.add_edge!(g, i, j + num_nodes)
        end
    end

    @series begin
        method --> :sfdp
        markercolor --> markercolor
        curves --> false
        arrows --> false
        names --> names
        nodeshape --> nodeshape
        node_weights --> node_weights
        titlefont --> font(12, "Arial", "bold")
        title --> title_str
        input_type := :sourcedestiny
        markerstrokewidth --> 0.0
        GraphRecipes.GraphPlot((g,))
    end
end

@recipe function plot_circuits(circuits::IMAS.IDSvector{<:IMAS.pf_active__circuit})
    layout := RecipesBase.@layout length(circuits)
    for (kkk, circuit) in enumerate(circuits)
        @series begin
            subplot := kkk
            circuit
        end
    end
end

# ========== #
# pf_passive #
# ========== #
@recipe function plot_pf_passive(pf_passive::IMAS.pf_passive)
    @series begin
        pf_passive.loop
    end
end

@recipe function plot_pf_passive(loops::IMAS.IDSvector{<:IMAS.pf_passive__loop})
    for loop in loops
        @series begin
            loop
        end
    end
end

@recipe function plot_loop(loop::pf_passive__loop{T}; loop_names=false) where {T<:Real}
    id = plot_help_id(loop)
    assert_type_and_record_argument(id, Bool, "Show loop names"; loop_names)

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    aspect_ratio := :equal
    r = T[]
    z = T[]
    for (k, element) in enumerate(loop.element)
        oute = outline(element)
        @series begin
            primary := k == 1
            seriestype --> :shape
            linewidth := 0.25 * base_linewidth
            colorbar --> :right
            label --> ""
            oute.r, oute.z
        end
        append!(r, oute.r)
        append!(z, oute.z)
    end

    if loop_names
        r_avg = sum(r) / length(r)
        z_avg = sum(z) / length(z)
        @series begin
            label := ""
            series_annotations := [(loop.name, :center, :middle, :red, 6)]
            [r_avg], [z_avg]
        end
    end
end

# ======= #
# costing #
# ======= #
@recipe function plot_costing(cstdc::IMAS.costing__cost_direct_capital)
    cstdcM = convert(Measurement{eltype(cstdc)}, cstdc)
    @series begin
        cstdcM
    end
end

@recipe function plot_costing(cstdc::IMAS.costing__cost_direct_capital{T}) where {T<:Measurement}
    costing = parent(cstdc)
    cst = cstdc.system

    names = ["$(sys_cst.name)" for sys_cst in cst]
    uncertain_costs = [sys_cst.cost for sys_cst in cst]
    costs = Float64[Measurements.value(sys_cst.cost) for sys_cst in cst]
    perc = ["$(round(sys_cst/sum(costs)*100))%" for sys_cst in costs]

    name_series = [["$(sub_cst.name)" for sub_cst in reverse(sys_cst.subsystem)] for sys_cst in cst]
    uncertain_cost_series = [[sub_cst.cost for sub_cst in reverse(sys_cst.subsystem)] for sys_cst in cst]
    cost_series = [Float64[Measurements.value(sub_cst.cost) for sub_cst in reverse(sys_cst.subsystem)] for sys_cst in cst]
    perc_series = [["$(round(sub_cst/sum(sys_cst)*1000)/10)%" for sub_cst in sys_cst] for sys_cst in cost_series]

    size --> (1000, 300)
    cols = 1 + length(filter!(!isempty, cost_series))
    layout := RecipesBase.@layout (1, cols)
    margin --> 5 * Measures.mm

    for s in (-1, 1)
        @series begin
            subplot := 1
            seriestype := :bar
            orientation := :horizontal
            alpha := 0.25
            linecolor := :match
            label := ""
            color := [PlotUtils.palette(:tab10)[c] for c in length(costs):-1:1]
            reverse(names), reverse([Measurements.value(c) + s * Measurements.uncertainty(c) for c in uncertain_costs])
        end
    end
    @series begin
        subplot := 1
        seriestype := :bar
        orientation := :horizontal
        title := "\n" * @sprintf("Direct Capital Cost in %d [%.3g \$\$B]", costing.construction_start_year, sum(costs) / 1E3)
        titlefontsize := 10
        ylim := (0, length(costs))
        label := ""
        annotation := [(0.0, kk - 0.5, ("   $x  $(titlecase(n,strict=false))", :left, 12)) for (kk, (c, x, n)) in enumerate(reverse!(collect(zip(costs, perc, names))))]
        annotationvalign := :center
        label := ""
        xticks := 0:1000:1000E3
        xlabel := "[\$M]"
        showaxis := :x
        yaxis := nothing
        alpha := 0.25
        primary := false
        linecolor := :match
        color := [PlotUtils.palette(:tab10)[c] for c in length(costs):-1:1]
        reverse(names), reverse(costs)
    end

    for (k, (sub_names, sub_perc, sub_costs, uncertain_sub_costs)) in enumerate(zip(name_series, perc_series, cost_series, uncertain_cost_series))
        if !isempty(sub_costs)
            for s in (-1, 1)
                @series begin
                    subplot := 1 + k
                    seriestype := :bar
                    orientation := :horizontal
                    alpha := 0.25
                    linecolor := :match
                    label := ""
                    color := PlotUtils.palette(:tab10)[k]
                    sub_names, [Measurements.value(c) + s * Measurements.uncertainty(c) for c in uncertain_sub_costs]
                end
            end
            @series begin
                subplot := 1 + k
                seriestype := :bar
                orientation := :horizontal
                title := titlecase("\n" * word_wrap(string(names[k]), 30) * " " * @sprintf("[%.3g \$\$B]", sum(sub_costs) / 1E3))
                titlefontsize := 10
                annotation := [(0.0, kk - 0.5, ("   $x  $(titlecase(n,strict=false))", :left, 8)) for (kk, (c, x, n)) in enumerate(zip(sub_costs, sub_perc, sub_names))]
                annotationvalign := :center
                label := ""
                yaxis := nothing
                if maximum(sub_costs) > 1E3
                    xticks := 0:500:1000E3
                else
                    xticks := 0:100:1000E3
                end
                xlabel := "[\$M]"
                showaxis := :x
                primary := false
                ylim := (0, length(sub_costs))
                alpha := 0.25
                linecolor := :match
                color := PlotUtils.palette(:tab10)[k]
                sub_names, sub_costs
            end
        end
    end
end

@recipe function plot_costing(costing::IMAS.costing)
    @series begin
        return costing.cost_direct_capital
    end
end

# =========== #
# equilibrium #
# =========== #
@recipe function plot_eq(eq::IMAS.equilibrium; time0=global_time(eq))
    id = plot_help_id(eq)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    @series begin
        return eq.time_slice[time0]
    end
end

@recipe function plot_eqt(eqt::IMAS.equilibrium__time_slice; cx=false, coordinate=:psi_norm, core_profiles_overlay=false)
    id = plot_help_id(eqt)
    assert_type_and_record_argument(id, Bool, "Display 2D cross section"; cx)
    assert_type_and_record_argument(id, Symbol, "X coordinate for 1D profiles"; coordinate)
    assert_type_and_record_argument(id, Bool, "Overlay core_profiles data"; core_profiles_overlay)

    if !cx
        layout := RecipesBase.@layout [a{0.35w} [a b; c d]]
        if Plots.backend_name() == :unicodeplots
            layout := 5
        end
        size --> (800, 500)
    end

    @series begin
        if !cx
            subplot := 1
        end
        eqt.profiles_2d
    end

    if !cx
        coordinate := coordinate

        if core_profiles_overlay
            cp = top_dd(eqt).core_profiles
            cp1d = top_dd(eqt).core_profiles.profiles_1d[argmin(abs.(cp.time .- eqt.time))]
        end

        # pressure
        if !ismissing(eqt.profiles_1d, :pressure)
            @series begin
                label := ""
                xlabel := ""
                subplot := 2
                normalization := 1E-6
                ylabel := ""
                title := latex_support() ? L"P~~\mathrm{[MPa]}" : "P [MPa]"
                eqt.profiles_1d, :pressure
            end
        end
        if core_profiles_overlay
            if !ismissing(cp1d, :pressure)
                @series begin
                    primary := false
                    ls := :dash
                    lw := 2
                    label := ""
                    xlabel := ""
                    subplot := 2
                    normalization := 1E-6
                    ylabel := ""
                    title := latex_support() ? L"P~~\mathrm{[MPa]}" : "P [MPa]"
                    cp1d, :pressure
                end
            end
        end

        # j_tor
        if !ismissing(eqt.profiles_1d, :j_tor)
            @series begin
                label := ""
                xlabel := ""
                subplot := 3
                normalization := 1E-6
                ylabel := ""
                title := latex_support() ? L"J_\mathrm{tor}~[\mathrm{MA/m^2}]" : "Jtor [MA/m²]"
                eqt.profiles_1d, :j_tor
            end
        end
        if core_profiles_overlay
            if !ismissing(cp1d, :pressure)
                @series begin
                    primary := false
                    ls := :dash
                    lw := 2
                    label := ""
                    xlabel := ""
                    subplot := 3
                    normalization := 1E-6
                    ylabel := ""
                    title := latex_support() ? L"J_\mathrm{tor}~[\mathrm{MA/m^2}]" : "Jtor [MA/m²]"
                    cp1d, :j_tor
                end
            end
        end

        # psi or rho_tor_norm
        if !ismissing(eqt.profiles_1d, :rho_tor_norm)
            @series begin
                label := ""
                subplot := 4
                ylabel := ""
                normalization := 1.0
                if contains(string(coordinate), "psi")
                    title := latex_support() ? L"\rho" : "ρ"
                    eqt.profiles_1d, :rho_tor_norm
                else
                    title := latex_support() ? L"\psi~~[\mathrm{Wb}]" : "Ψ [Wb]"
                    eqt.profiles_1d, :psi
                end
            end
        end

        # q
        if !ismissing(eqt.profiles_1d, :q)
            @series begin
                label := ""
                primary := false
                subplot := 5
                seriestype --> :hline
                color := :gray
                linestyle := :dash
                [sum(eqt.profiles_1d.q) > 0 ? 1.0 : -1.0]
            end
            @series begin
                label := ""
                subplot := 5
                ylabel := ""
                normalization := 1.0
                title := latex_support() ? L"q" : "q"
                eqt.profiles_1d, :q
            end
        end
    end
end

@recipe function plot_eqt2dv(eqt2dv::IDSvector{<:IMAS.equilibrium__time_slice___profiles_2d})
    if !isempty(eqt2dv)
        if ismissing(eqt2dv[1], :psi)
            @series begin
                eqt2dv[1].boundary.outline
            end
        else
            return eqt2dv[1]
        end
    end
end

@recipe function plot_eqt2(
    eqt2d::IMAS.equilibrium__time_slice___profiles_2d;
    psi_levels_in=11,
    psi_levels_out=11,
    show_secondary_separatrix=false,
    show_x_points=true,
    show_strike_points=true,
    show_magnetic_axis=true)

    id = plot_help_id(eqt2d)
    assert_type_and_record_argument(id, Union{Int,AbstractVector{<:Real}}, "Psi levels inside LCFS"; psi_levels_in)
    assert_type_and_record_argument(id, Union{Int,AbstractVector{<:Real}}, "Psi levels outside LCFS"; psi_levels_out)
    assert_type_and_record_argument(id, Bool, "Plot secondary separatrix"; show_secondary_separatrix)
    assert_type_and_record_argument(id, Bool, "Show X points"; show_x_points)
    assert_type_and_record_argument(id, Bool, "Show strike points"; show_strike_points)
    assert_type_and_record_argument(id, Bool, "Show magnetic axis"; show_magnetic_axis)

    label --> ""
    aspect_ratio := :equal
    @series begin
        [], []
    end
    primary --> false

    eqt = parent(parent(eqt2d))

    # plot cx
    # handle psi levels
    psi__boundary_level = eqt.profiles_1d.psi[end]
    tmp = eqt.profiles_1d.psi[end]
    if tmp !== nothing
        psi__boundary_level = tmp
    end
    if typeof(psi_levels_in) <: Int
        if psi_levels_in == 1
            psi_levels_in = [psi__boundary_level]
        elseif psi_levels_in > 0
            psi_levels_in = range(eqt.profiles_1d.psi[1], psi__boundary_level, psi_levels_in)
        else
            psi_levels_in = []
        end
    end
    delta_psi = (psi__boundary_level - eqt.profiles_1d.psi[1])
    if typeof(psi_levels_out) <: Int
        if psi_levels_out > 0
            psi_levels_out = delta_psi / length(psi_levels_in) .* collect(0:psi_levels_out) .+ psi__boundary_level
        else
            psi_levels_out = []
        end
    end
    psi_levels = unique(vcat(psi_levels_in, psi_levels_out))

    fw = first_wall(top_dd(eqt).wall)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    for psi_level in psi_levels
        for (pr, pz) in flux_surface(eqt, psi_level, :any, fw.r, fw.z)
            @series begin
                seriestype --> :path
                if psi_level == psi__boundary_level
                    linewidth := base_linewidth * 1.5
                else
                    linewidth := base_linewidth * 0.5
                    if psi_level in psi_levels_in &&
                       (is_closed_polygon(pr, pz) && (PolygonOps.inpolygon((RA, ZA), collect(zip(pr, pz))) == 1) && !IMAS.intersects(pr, pz, fw.r, fw.z))
                        linestyle --> :solid
                    else
                        linestyle --> :dash
                    end
                end
                pr, pz
            end
        end
    end

    if show_secondary_separatrix
        @series begin
            primary --> false
            linewidth := base_linewidth * 1.0
            eqt.boundary_secondary_separatrix.outline.r, eqt.boundary_secondary_separatrix.outline.z
        end
    end

    if show_magnetic_axis
        @series begin
            primary --> false
            eqt.global_quantities.magnetic_axis
        end
    end

    if show_x_points
        @series begin
            primary --> false
            eqt.boundary.x_point
        end
    end

    if show_strike_points
        @series begin
            primary --> false
            eqt.boundary.strike_point
        end
    end

end

@recipe function plot_eqtb(eqtb::IMAS.equilibrium__time_slice___boundary)
    @series begin
        eqtb.outline
    end
    @series begin
        primary --> false
        eqtb.x_point
    end
    @series begin
        primary --> false
        eqtb.strike_point
    end
end

@recipe function plot_eqtb(eqtbo::IMAS.equilibrium__time_slice___boundary__outline)
    label --> ""
    aspect_ratio := :equal
    @series begin
        eqtbo.r, eqtbo.z
    end
end

@recipe function plot_x_points(x_points::IDSvector{<:IMAS.equilibrium__time_slice___boundary__x_point})
    for x_point in x_points
        @series begin
            x_point
        end
    end
end

@recipe function plot_x_point(x_point::IMAS.equilibrium__time_slice___boundary__x_point)
    @series begin
        seriestype := :scatter
        marker --> :circle
        markerstrokewidth --> 0
        label --> ""
        aspect_ratio := :equal
        [(x_point.r, x_point.z)]
    end
end

@recipe function plot_strike_points(s_points::IDSvector{<:IMAS.equilibrium__time_slice___boundary__strike_point})
    for s_point in s_points
        @series begin
            s_point
        end
    end
end

@recipe function plot_strike_point(s_point::IMAS.equilibrium__time_slice___boundary__strike_point)
    @series begin
        seriestype := :scatter
        marker --> :cross
        markerstrokewidth --> 0
        label --> ""
        aspect_ratio := :equal
        [(s_point.r, s_point.z)]
    end
end

@recipe function plot_mag_axis(mag_axis::IMAS.equilibrium__time_slice___global_quantities__magnetic_axis)
    @series begin
        seriestype --> :scatter
        markershape --> :cross
        label --> ""
        aspect_ratio := :equal
        [(mag_axis.r, mag_axis.z)]
    end
end

# ========= #
# divertors #
# ========= #
@recipe function plot_divertors(divertors::IMAS.divertors)
    @series begin
        divertors.divertor
    end
end

@recipe function plot_divertors(divertors::IDSvector{<:IMAS.divertors__divertor})
    for divertor in divertors
        @series begin
            divertor
        end
    end
end

@recipe function plot_divertors(divertor::IMAS.divertors__divertor)
    @series begin
        divertor.target
    end
end

@recipe function plot_divertors(targets::IDSvector{<:IMAS.divertors__divertor___target})
    for target in targets
        @series begin
            target
        end
    end
end

@recipe function plot_divertors(target::IMAS.divertors__divertor___target)
    @series begin
        label --> target.name
        target.tile
    end
end

@recipe function plot_divertors(tiles::IDSvector{<:IMAS.divertors__divertor___target___tile})
    for tile in tiles
        @series begin
            tile
        end
    end
end

@recipe function plot_divertors(tile::IMAS.divertors__divertor___target___tile)
    @series begin
        label --> tile.name
        tile.surface_outline
    end
end

@recipe function plot_divertors(surface_outline::IMAS.divertors__divertor___target___tile___surface_outline)
    base_linewidth = get(plotattributes, :linewidth, 1.0)
    @series begin
        linewidth := base_linewidth * 2
        aspect_ratio := :equal
        surface_outline.r, surface_outline.z
    end
end

# ===== #
# build #
# ===== #
function join_outlines(r1::AbstractVector{T}, z1::AbstractVector{T}, r2::AbstractVector{T}, z2::AbstractVector{T}) where {T<:Real}
    i1, i2 = minimum_distance_polygons_vertices(r1, z1, r2, z2; return_index=true)
    r = vcat(reverse(r1[i1:end]), r2[i2:end], r2[1:i2], reverse(r1[1:i1]))
    z = vcat(reverse(z1[i1:end]), z2[i2:end], z2[1:i2], reverse(z1[1:i1]))
    return r, z
end

function layer_attrs(l::IMAS.build__layer)
    name = l.name
    color = :gray
    if l.material == "water"
        name = ""
        color = :lightblue
    elseif l.type == Int(_gap_)
        name = ""
        color = :white
    elseif l.type == Int(_oh_)
        if l.side == _in_
            color = :gray
        else
            color = :white
        end
    elseif l.type == Int(_tf_)
        color = :green
    elseif l.type == Int(_shield_)
        color = :red
    elseif l.type == Int(_blanket_)
        color = :orange
    elseif l.type == Int(_wall_)
        color = :yellow
    elseif l.type == Int(_vessel_)
        color = :yellow
    elseif l.type == Int(_cryostat_)
        color = :lightgray
    end
    for nm in ("inner", "outer", "vacuum", "hfs", "lfs", "gap")
        name = replace(name, Regex("$(nm) ", "i") => "")
    end
    return name, color
end



@recipe function plot_build_cx(bd::IMAS.build; cx=true, wireframe=false, equilibrium=true, pf_active=true, pf_passive=true, only=Symbol[], exclude_layers=Symbol[])
    id = plot_help_id(bd)
    assert_type_and_record_argument(id, Bool, "Plot cross section"; cx)
    assert_type_and_record_argument(id, Bool, "Use wireframe"; wireframe)
    assert_type_and_record_argument(id, Bool, "Include plot of equilibrium"; equilibrium)
    assert_type_and_record_argument(id, Bool, "Include plot of pf_active"; pf_active)
    assert_type_and_record_argument(id, Bool, "Include plot of pf_passive"; pf_passive)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "Only include certain layers"; only)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "Exclude certain layers"; exclude_layers)

    dd = top_dd(bd)

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    legend_position --> :outerbottomright
    aspect_ratio := :equal
    grid --> :none

    # cx
    if cx

        if dd !== nothing && equilibrium
            @series begin
                cx := true
                dd.equilibrium
            end
        end

        rmax = maximum(bd.layer[end].outline.r)

        # everything after first vacuum in _out_
        if (isempty(only) || (:cryostat in only)) && :cryostat ∉ exclude_layers
            for k in get_build_indexes(bd.layer; fs=_out_)[2:end]
                if !wireframe
                    @series begin
                        seriestype --> :shape
                        linewidth := 0.0
                        color --> :lightgray
                        label --> (!wireframe ? "Cryostat" : "")
                        xlim --> [0, rmax]
                        join_outlines(
                            bd.layer[k].outline.r,
                            bd.layer[k].outline.z,
                            bd.layer[k-1].outline.r,
                            bd.layer[k-1].outline.z
                        )
                    end
                end
                @series begin
                    seriestype --> :path
                    linewidth := base_linewidth
                    color --> :black
                    label --> ""
                    xlim --> [0, rmax]
                    bd.layer[k].outline.r[1:end-1], bd.layer[k].outline.z[1:end-1]
                end
            end
        end

        # first vacuum in _out_
        if !wireframe
            k = get_build_indexes(bd.layer; fs=_out_)[1]
            @series begin
                seriestype --> :shape
                linewidth := 0.0
                color --> :white
                linecolor --> :white
                label --> ""
                xlim --> [0, rmax]
                join_outlines(
                    bd.layer[k].outline.r,
                    bd.layer[k].outline.z,
                    get_build_layer(bd.layer; type=_plasma_).outline.r,
                    get_build_layer(bd.layer; type=_plasma_).outline.z
                )
            end
            if (isempty(only) || (:cryostat in only)) && :cryostat ∉ exclude_layers
                @series begin
                    seriestype --> :path
                    linewidth := base_linewidth
                    color --> :black
                    label --> ""
                    xlim --> [0, rmax]
                    bd.layer[k].outline.r[1:end-1], bd.layer[k].outline.z[1:end-1]
                end
            end
        end

        # axis of symmetry
        @series begin
            seriestype --> :path
            linewidth := base_linewidth
            label --> ""
            linestyle --> :dash
            color --> :black
            [0.0, 0.0], [minimum(bd.layer[end].outline.z), maximum(bd.layer[end].outline.z)]
        end

        # all layers inside of the TF
        if (isempty(only) || (:oh in only)) && (!(:oh in exclude_layers))
            for k in get_build_indexes(bd.layer; fs=_in_)
                layer = bd.layer[k]
                if getproperty(layer, :material, "not_vacuum") != "vacuum"
                    if !wireframe
                        @series begin
                            seriestype --> :shape
                            linewidth := 0.0
                            color --> :gray
                            label --> (!wireframe ? layer.name : "")
                            xlim --> [0, rmax]
                            layer.outline.r, layer.outline.z
                        end
                    end
                    @series begin
                        seriestype --> :path
                        linewidth := 0.5 * base_linewidth
                        color --> :black
                        label --> ""
                        xlim --> [0, rmax]
                        layer.outline.r, layer.outline.z
                    end
                end
            end
        end

        # all layers between the OH and the plasma
        for k in get_build_indexes(bd.layer; fs=_hfs_)
            l = bd.layer[k]
            l1 = bd.layer[k+1]
            poly = join_outlines(l.outline.r, l.outline.z, l1.outline.r, l1.outline.z)

            # setup labels and colors
            name, color = layer_attrs(l)

            if (isempty(only) || (Symbol(name) in only)) && (!(Symbol(name) in exclude_layers))
                if !wireframe
                    @series begin
                        seriestype --> :shape
                        linewidth := 0.0
                        color --> color
                        label --> uppercasefirst(name)
                        xlim --> [0, rmax]
                        poly[1], poly[2]
                    end
                end
                @series begin
                    seriestype --> :path
                    linewidth := 0.5 * base_linewidth
                    color --> :black
                    label --> ""
                    xlim --> [0, rmax]
                    l.outline.r, l.outline.z
                end
            end
        end

        # plasma
        if (isempty(only) || (:plasma in only)) && (!(:plasma in exclude_layers))
            plasma_outline = outline(get_build_layer(bd.layer; type=_plasma_))
            @series begin
                seriestype --> :path
                linewidth := base_linewidth
                color --> :black
                label --> ""
                xlim --> [0, rmax]
                plasma_outline.r, plasma_outline.z
            end
        end

        if any(structure.type == Int(_divertor_) for structure in bd.structure)
            if (isempty(only) || (:divertor in only)) && (!(:divertor in exclude_layers))
                for (k, index) in enumerate(findall(x -> x.type == Int(_divertor_), bd.structure))
                    structure_outline = outline(bd.structure[index])
                    if !wireframe
                        @series begin
                            seriestype --> :shape
                            linewidth := 0.0
                            color --> :mediumpurple1
                            label --> (!wireframe && k == 1 ? "Divertor" : "")
                            xlim --> [0, rmax]
                            structure_outline.r, structure_outline.z
                        end
                    end
                    @series begin
                        seriestype --> :path
                        linewidth := 0.5 * base_linewidth
                        color --> :black
                        label --> ""
                        xlim --> [0, rmax]
                        structure_outline.r, structure_outline.z
                    end
                end
            end
        end

        # maintenance port (if any)
        if any(structure.type == Int(_port_) for structure in bd.structure)
            if (isempty(only) || (:port in only)) && (!(:port in exclude_layers))
                for (k, index) in enumerate(findall(x -> x.type == Int(_port_), bd.structure))
                    if !wireframe
                        @series begin
                            seriestype --> :path
                            linewidth := 2.0 * base_linewidth
                            color --> :lightblue
                            label --> (!wireframe && k == 1 ? "Vessel Port" : "")
                            xlim --> [0, maximum(bd.structure[index].outline.r)]
                            bd.structure[index].outline.r, bd.structure[index].outline.z
                        end
                    end
                    @series begin
                        seriestype --> :path
                        linewidth := 0.5 * base_linewidth
                        color --> :black
                        label --> ""
                        xlim --> [0, maximum(bd.structure[index].outline.r)]
                        bd.structure[index].outline.r, bd.structure[index].outline.z
                    end
                end
            end
        end

        if dd !== nothing && pf_active
            @series begin
                colorbar --> :false
                xlim --> [0, rmax]
                dd.pf_active
            end
        end

        if dd !== nothing && pf_passive
            @series begin
                colorbar --> :false
                xlim --> [0, rmax]
                color := :yellow
                dd.pf_passive
            end
        end

    else  # not-cx

        @series begin
            seriestype --> :vline
            linewidth := 2 * base_linewidth
            label --> ""
            linestyle --> :dash
            color --> :black
            [0]
        end

        at = 0
        for l in bd.layer
            name, color = layer_attrs(l)
            @series begin
                seriestype --> :vspan
                label --> name
                color --> color
                alpha --> 0.2
                xlim --> [0, at]
                [at, at + l.thickness]
            end
            at += l.thickness
            @series begin
                seriestype --> :vline
                linewidth := 0.5 * base_linewidth
                label --> ""
                color --> :black
                xlim --> [0, at]
                [at]
            end
        end
    end
end

@recipe function plot_build_layer_outline(layer::IMAS.build__layer)
    @series begin
        aspect_ratio := :equal
        label --> layer.name
        layer.outline.r, layer.outline.z
    end
end

@recipe function plot_build_structure_outline(structure::IMAS.build__structure)
    @series begin
        aspect_ratio := :equal
        label --> structure.name
        structure.outline.r, structure.outline.z
    end
end

# ======== #
# build tf #
# ======== #
@recipe function plot_build_tf(tf::IMAS.build__tf)
    layers = parent(tf).layer
    TF = get_build_layers(layers; type=IMAS._tf_)

    for n in 1:tf.coils_n
        x, y = top_outline(tf, n)
        @series begin
            linewidth := 2
            label := ""
            seriestype := :shape
            x, y
        end
    end

    ϕ = range(0, 2pi, 101)
    @series begin
        color := :black
        label := ""
        TF[1].start_radius .* cos.(ϕ), TF[1].start_radius .* sin.(ϕ)
    end
    @series begin
        primary := false
        color := :black
        label := ""
        TF[1].end_radius .* cos.(ϕ), TF[1].end_radius .* sin.(ϕ)
    end
    @series begin
        primary := false
        color := :black
        label := ""
        TF[2].start_radius .* cos.(ϕ), TF[2].start_radius .* sin.(ϕ)
    end
    @series begin
        primary := false
        color := :black
        label := ""
        aspect_ratio := :equal
        TF[2].end_radius .* cos.(ϕ), TF[2].end_radius .* sin.(ϕ)
    end
end

# ========= #
# transport #
# ========= #
@recipe function plot_ct1d(ct1d__electrons__energy::IMAS.core_transport__model___profiles_1d___electrons__energy, ::Val{:flux})
    @series begin
        markershape --> :none
        title := "Electron energy flux"
        label --> ""
        normalization := 1E-6
        ylabel := "[MW/m²]"
        ct1d__electrons__energy, :flux
    end
end

@recipe function plot_ct1d(ct1d__total_ion_energy::IMAS.core_transport__model___profiles_1d___total_ion_energy, ::Val{:flux})
    @series begin
        markershape --> :none
        title := "Total ion energy flux"
        label --> ""
        normalization := 1E-6
        ylabel := "[MW/m²]"
        ct1d__total_ion_energy, :flux
    end
end

@recipe function plot_ct1d(ct1d__electrons__particles::IMAS.core_transport__model___profiles_1d___electrons__particles, ::Val{:flux})
    @series begin
        markershape --> :none
        title := "Electron particle flux"
        label --> ""
        ylabel := "[s⁻¹/m²]"
        ct1d__electrons__particles, :flux
    end
end

@recipe function plot_ct1d(ct1d__ion___particles::IMAS.core_transport__model___profiles_1d___ion___particles, ::Val{:flux})
    ion = parent(ct1d__ion___particles)
    @series begin
        markershape --> :none
        title := "$(ion.label) particle flux"
        label --> ""
        ylabel := "[s⁻¹/m²]"
        ct1d__ion___particles, :flux
    end
end

@recipe function plot_ct1d(ct1d__momentum_tor::IMAS.core_transport__model___profiles_1d___momentum_tor, ::Val{:flux})
    @series begin
        markershape --> :none
        title := "Momentum flux"
        label --> ""
        ylabel := "[Nm/m²]"
        ct1d__momentum_tor, :flux
    end
end

function transport_channel_paths(ions::AbstractVector{<:IDS}, ions_list::AbstractVector{Symbol})
    paths = Any[
        [:electrons, :energy],
        [:total_ion_energy],
        [:electrons, :particles],
        [:momentum_tor]]
    available_ions = [Symbol(ion.label) for ion in ions]
    for (kion, ion) in enumerate(ions_list)
        if ion in available_ions
            push!(paths, [:ion, kion, :particles])
        else
            push!(paths, [nothing])
        end
    end
    return paths
end

@recipe function plot_ct1d(ct1d::IMAS.core_transport__model___profiles_1d; only=nothing, ions=Symbol[])
    id = plot_help_id(ct1d)
    assert_type_and_record_argument(id, Union{Nothing,Int}, "Plot only this subplot number"; only)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)

    paths = transport_channel_paths(ct1d.ion, ions)

    if only === nothing
        layout := 4 + length(ions)
        Nr = Int(ceil(sqrt(4 + length(ions))))
        size --> (1100, Nr * 290)
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)
    end

    for (k, path) in enumerate(paths)
        if k == only || only === nothing
            if path[end] !== nothing && !ismissing(ct1d, [path; :flux])
                @series begin
                    if only === nothing
                        subplot := k
                    end
                    goto(ct1d, path), Val(:flux)
                end
            end
        end
    end
end

@recipe function plot_core_transport(ct::IMAS.core_transport{D}; ions=Symbol[:my_ions], time0=global_time(ct)) where {D<:Real}
    id = plot_help_id(ct)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    dd = top_dd(ct)
    cp1d = dd.core_profiles.profiles_1d[]

    if ions == [:my_ions]
        ions = list_ions(ct, cp1d)
    end

    model_type = name_2_index(ct.model)
    rhos = D[]
    for model in ct.model
        if model.identifier.index ∈ (model_type[k] for k in (:combined, :unspecified, :transport_solver, :unknown))
            continue
        end
        ct1d = model.profiles_1d[]
        append!(rhos, ct1d.grid_flux.rho_tor_norm)
    end
    rhos = unique(rhos)

    if dd !== nothing
        @series begin
            linewidth := 2
            color := :blue
            name := "Total source"
            flux := true
            show_zeros := true
            ions := ions
            total_sources(dd; time0)
        end
    end

    @series begin
        linewidth := 2
        color := :red
        label := "Total transport"
        ions := ions
        total_fluxes(ct, cp1d, rhos; time0)
    end

    for model in ct.model
        if model.identifier.index ∈ (model_type[k] for k in (:combined, :unspecified, :transport_solver, :unknown))
            continue
        end
        @series begin
            ions := ions
            label := model.identifier.name
            if model.identifier.index == model_type[:anomalous]
                markershape := :diamond
                markerstrokewidth := 0.5
                linewidth := 0
                color := :orange
            elseif model.identifier.index == model_type[:neoclassical]
                markershape := :cross
                linewidth := 0
                color := :purple
            end
            model.profiles_1d[time0]
        end
    end
end

# ============ #
# core_sources #
# ============ #
@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:electrons__energy})
    return cs1d.electrons, Val(:energy)
end

@recipe function plot_source1d(cs1de::IMAS.core_sources__source___profiles_1d___electrons, v::Val{:energy};
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    min_power=1e3,
    only_positive_negative=0,
    show_source_number=false)

    id = plot_help_id(cs1de, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Float64, "Minimum power threshold"; min_power)
    assert_type_and_record_argument(id, Int, "Show only positive or negative values (0 for all)"; only_positive_negative)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    cs1d = parent(cs1de)

    name, identifier, idx = source_name_identifier(parent(parent(cs1d); error_parent_of_nothing=false), name, show_source_number)

    tot = 0.0
    if !ismissing(cs1de, :energy)
        tot = trapz(cs1d.grid.volume, cs1de.energy)
    end
    show_condition =
        show_zeros || identifier in [:total, :collisional_equipartition, :time_derivative] ||
        (abs(tot) > min_power && (only_positive_negative == 0 || sign(tot) == sign(only_positive_negative)))
    @series begin
        if identifier == :collisional_equipartition
            linestyle --> :dash
        end
        color := idx
        title --> "Electron Energy"
        if show_condition
            label := "$name " * @sprintf("[%.3g MW]", tot / 1E6) * label
            if identifier in [:ec, :ic, :lh, :nbi, :pellet]
                fill0 --> true
            end
            if !ismissing(cs1de, :power_inside) && flux
                ylabel := "[MW/m²]"
                normalization = 1E-6 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1de, :power_inside
            elseif !integrated && !ismissing(cs1de, :energy)
                cs1de, :energy
            elseif integrated && !ismissing(cs1de, :power_inside)
                cs1de, :power_inside
            else
                label := ""
                [NaN], [NaN]
            end
        else
            label := ""
            [NaN], [NaN]
        end
    end
end

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:total_ion_energy};
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    min_power=1e3,
    only_positive_negative=0,
    show_source_number=false)

    id = plot_help_id(cs1d, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Float64, "Minimum power threshold"; min_power)
    assert_type_and_record_argument(id, Int, "Show only positive or negative values (0 for all)"; only_positive_negative)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)
    name, identifier, idx = source_name_identifier(parent(parent(cs1d); error_parent_of_nothing=false), name, show_source_number)

    tot = 0.0
    if !ismissing(cs1d, :total_ion_energy)
        tot = trapz(cs1d.grid.volume, cs1d.total_ion_energy)
    end
    show_condition =
        show_zeros || identifier in [:total, :collisional_equipartition, :time_derivative] ||
        (abs(tot) > min_power && (only_positive_negative == 0 || sign(tot) == sign(only_positive_negative)))
    @series begin
        if identifier == :collisional_equipartition
            linestyle --> :dash
        end
        color := idx
        title --> "Total Ion Energy"
        if show_condition
            label := "$name " * @sprintf("[%.3g MW]", tot / 1E6) * label
            if identifier in [:ec, :ic, :lh, :nbi, :pellet]
                fill0 --> true
            end
            if !ismissing(cs1d, :total_ion_power_inside) && flux
                ylabel := "[MW/m²]"
                normalization = 1E-6 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1d, :total_ion_power_inside
            elseif !integrated && !ismissing(cs1d, :total_ion_energy)
                cs1d, :total_ion_energy
            elseif integrated && !ismissing(cs1d, :total_ion_power_inside)
                cs1d, :total_ion_power_inside
            else
                label := ""
                [NaN], [NaN]
            end
        else
            label := ""
            [NaN], [NaN]
        end
    end
end

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:electrons__particles})
    return cs1d.electrons, Val(:particles)
end

@recipe function plot_source1d(cs1de::IMAS.core_sources__source___profiles_1d___electrons, v::Val{:particles};
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    show_source_number=false)

    id = plot_help_id(cs1de, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    cs1d = parent(cs1de)

    name, identifier, idx = source_name_identifier(parent(parent(cs1d); error_parent_of_nothing=false), name, show_source_number)

    tot = 0.0
    if !ismissing(cs1de, :particles)
        tot = trapz(cs1d.grid.volume, cs1de.particles)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        color := idx
        title --> "Electron Particles"
        if show_condition
            label := "$name " * @sprintf("[%.3g s⁻¹]", tot) * label
            if identifier in [:ec, :ic, :lh, :nbi, :pellet]
                fill0 --> true
            end
            if !ismissing(cs1de, :particles_inside) && flux
                ylabel := "[s⁻¹/m²]"
                normalization = 1.0 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1de, :particles_inside
            elseif !integrated && !ismissing(cs1de, :particles)
                cs1de, :particles
            elseif integrated && !ismissing(cs1de, :particles_inside)
                cs1de, :particles_inside
            else
                label := ""
                [NaN], [NaN]
            end
        else
            label := ""
            [NaN], [NaN]
        end
    end
end

@recipe function plot_source1d(cs1di::IMAS.core_sources__source___profiles_1d___ion, v::Val{:particles};
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    show_source_number=false)

    id = plot_help_id(cs1di, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    cs1d = parent(parent(cs1di))
    name, identifier, idx = source_name_identifier(parent(parent(cs1d); error_parent_of_nothing=false), name, show_source_number)

    tot = 0.0
    if !ismissing(cs1di, :particles)
        tot = trapz(cs1d.grid.volume, cs1di.particles)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        color := idx
        title --> "$(cs1di.label) Particles"
        if show_condition
            label := "$name $(cs1di.label) " * @sprintf("[%.3g s⁻¹]", tot) * label
            if identifier in [:ec, :ic, :lh, :nbi, :pellet]
                fill0 --> true
            end
            if !ismissing(cs1di, :particles_inside) && flux
                ylabel := "[s⁻¹/m²]"
                normalization = 1.0 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1di, :particles_inside
            elseif !integrated && !ismissing(cs1di, :particles)
                cs1di, :particles
            elseif integrated && !ismissing(cs1di, :particles_inside)
                cs1di, :particles_inside
            else
                label := ""
                [NaN], [NaN]
            end
        else
            label := ""
            [NaN], [NaN]
        end
    end
end

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:momentum_tor};
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    show_source_number=false)

    id = plot_help_id(cs1d, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    name, identifier, idx = source_name_identifier(parent(parent(cs1d); error_parent_of_nothing=false), name, show_source_number)

    tot = 0.0
    if !ismissing(cs1d, :momentum_tor)
        tot = trapz(cs1d.grid.volume, cs1d.momentum_tor)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        color := idx
        title --> "Momentum"
        if show_condition
            label := "$name " * @sprintf("[%.3g N m]", tot / 1E6) * label
            if identifier in [:ec, :ic, :lh, :nbi, :pellet]
                fill0 --> true
            end
            if !ismissing(cs1d, :torque_tor_inside) && flux
                ylabel := "[Nm/m²]"
                normalization = 1.0 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1d, :torque_tor_inside
            elseif !integrated && !ismissing(cs1d, :momentum_tor)
                cs1d, :momentum_tor
            elseif integrated && !ismissing(cs1d, :torque_tor_inside)
                cs1d, :torque_tor_inside
            else
                label := ""
                [NaN], [NaN]
            end
        else
            label := ""
            [NaN], [NaN]
        end
    end
end

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:j_parallel};
    name="",
    label="",
    integrated=false,
    show_zeros=false,
    show_source_number=false)

    id = plot_help_id(cs1d, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    name, identifier, idx = source_name_identifier(parent(parent(cs1d); error_parent_of_nothing=false), name, show_source_number)

    tot = 0.0
    if !ismissing(cs1d, :j_parallel)
        tot = trapz(cs1d.grid.area, cs1d.j_parallel)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        color := idx
        title --> "Parallel Current"
        if show_condition
            label := "$name " * @sprintf("[%.3g MA]", tot / 1E6) * label
            if !integrated && !ismissing(cs1d, :j_parallel)
                if identifier in [:ec, :ic, :lh, :nbi, :pellet]
                    fill0 --> true
                end
                cs1d, :j_parallel
            elseif integrated && !ismissing(cs1d, :current_parallel_inside)
                cs1d, :current_parallel_inside
            else
                label := ""
                [NaN], [NaN]
            end
        else
            label := ""
            [NaN], [NaN]
        end
    end
end

function source_name_identifier(source::Union{IMAS.core_sources__source,Nothing}, name::AbstractString, show_source_number::Bool)
    if parent(source; error_parent_of_nothing=false) !== nothing
        idx = index(source)
        if show_source_number
            name = "[$idx] $name"
        end
    else
        source = nothing
        idx = 1
    end
    if source !== nothing
        identifier = identifier_name(source)
    else
        identifier = :undefined
    end
    return name, identifier, idx
end

@recipe function plot_source1d(
    cs1d::IMAS.core_sources__source___profiles_1d;
    name="",
    label="",
    integrated=false,
    flux=false,
    only=nothing,
    show_zeros=false,
    min_power=1e3,
    only_positive_negative=0,
    show_source_number=false,
    ions=Symbol[])

    id = plot_help_id(cs1d)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Union{Nothing,Int}, "Plot only this subplot number"; only)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Float64, "Minimum power threshold"; min_power)
    assert_type_and_record_argument(id, Int, "Show only positive or negative values (0 for all)"; only_positive_negative)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)

    paths = transport_channel_paths(cs1d.ion, ions)
    if !flux
        insert!(paths, 5, [:j_parallel])
    end
    if only === nothing
        if flux
            N = 4 + length(ions)
        else
            N = 5 + length(ions)
        end
        layout := N
        Nr = Int(ceil(sqrt(N)))
        size --> (1100, Nr * 290)
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)
    end

    for (k, path) in enumerate(paths)
        if k == only || only === nothing
            @series begin
                if only === nothing
                    subplot := k
                end
                if path[end] !== nothing && !ismissing(cs1d, path)
                    name := name
                    show_source_number := show_source_number
                    label := label
                    integrated := integrated
                    flux := flux
                    show_zeros := show_zeros
                    min_power := min_power
                    only_positive_negative := only_positive_negative
                    if path[1] == :ion
                        goto(cs1d, path[1:end-1]), Val(path[end])
                    else
                        cs1d, Val(Symbol(join(filter(x -> !(typeof(x) <: Int), path), "__")))
                    end
                else
                    [NaN], [NaN]
                end
            end
        end
    end
end

@recipe function plot_core_sources(cs::IMAS.core_sources{T}; ions=[:my_ions], time0=global_time(cs), aggregate_radiation=false) where {T<:Real}
    id = plot_help_id(cs)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Bool, "Aggregate radiation sources"; aggregate_radiation)

    dd = top_dd(cs)

    if ions == [:my_ions]
        if dd !== nothing && !isempty(dd.core_profiles.profiles_1d)
            cp1d = dd.core_profiles.profiles_1d[]
            ions = list_ions(cs, cp1d)
        else
            ions = list_ions(cs, nothing)
        end
    end

    all_indexes = [source.identifier.index for source in cs.source]
    for source in cs.source
        if !retain_source(source, all_indexes, Int[], Int[])
            continue
        end
        @series begin
            if aggregate_radiation
                only_positive_negative := 1
            end
            nozeros := true
            time0 := time0
            ions := ions
            source
        end
    end

    dd = top_dd(cs)
    if dd !== nothing

        if aggregate_radiation
            @series begin
                rad_source = IMAS.core_sources__source{T}()
                resize!(rad_source.profiles_1d, 1)
                merge!(rad_source.profiles_1d[1], total_radiation_sources(dd; time0))
                rad_source.identifier.index = 200
                rad_source.identifier.name = "radiation"
                ions := ions
                rad_source
            end
        end

        @series begin
            name := "total"
            linewidth := 2
            color := :black
            min_power := 0.0
            ions := ions
            total_sources(dd; time0)
        end
    end
end

@recipe function plot_source(source::IMAS.core_sources__source; time0=global_time(source))
    id = plot_help_id(source)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    if !isempty(source.profiles_1d) && source.profiles_1d[1].time <= time0
        @series begin
            name := source.identifier.name
            source.profiles_1d[time0]
        end
    end
end

# ============= #
# core_profiles #
# ============= #
@recipe function plot_core_profiles(cp::IMAS.core_profiles; time0=global_time(cp))
    id = plot_help_id(cp)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    @series begin
        return cp.profiles_1d[time0]
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d; label=nothing, only=nothing, greenwald=false)
    id = plot_help_id(cpt)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Union{Nothing,Int}, "Plot only this subplot number"; only)
    assert_type_and_record_argument(id, Bool, "Include Greenwald density"; greenwald)

    if label === nothing
        label = ""
    end

    if only === nothing
        layout := (1, 3)
        size --> (1100, 290)
        margin --> 5 * Measures.mm
    end

    # temperatures
    if only === nothing || only == 1
        @series begin
            if only === nothing
                subplot := 1
            end
            title := "Temperatures"
            label := "e" * label
            ylim --> (0, Inf)
            cpt.electrons, :temperature
        end

        same_temps = false
        if length(cpt.ion) > 1
            same_temps = !any(x -> x == false, [iion.temperature == cpt.ion[1].temperature for iion in cpt.ion[2:end] if !ismissing(iion, :temperature)])
            if same_temps
                @series begin
                    if only === nothing
                        subplot := 1
                    end
                    title := "Temperatures"
                    label := "Ions" * label
                    linestyle --> :dash
                    ylim --> (0, Inf)
                    cpt.ion[1], :temperature
                end
            end
        end

        if !same_temps
            for ion in cpt.ion
                if !ismissing(ion, :temperature)
                    @series begin
                        if only === nothing
                            subplot := 1
                        end
                        title := "Temperatures"
                        label := ion.label * label
                        linestyle --> :dash
                        ylim --> (0, Inf)
                        ion, :temperature
                    end
                end
            end
        end
    end

    # densities
    if only === nothing || only == 2
        @series begin
            if only === nothing
                subplot := 2
            end
            title := "Densities"
            label := "e" * label
            ylim --> (0.0, Inf)
            cpt.electrons, :density
        end
        if greenwald
            @series begin
                if only === nothing
                    subplot := 2
                end
                seriestype := :hline
                primary := false
                style := :dashdotdot
                [IMAS.greenwald_density(IMAS.top_dd(cpt).equilibrium.time_slice[cpt.time])]
            end
        end
        for ion in cpt.ion
            @series begin
                Z = ion.element[1].z_n
                if only === nothing
                    subplot := 2
                end
                title := "Densities"
                if Z == 1.0
                    label := ion.label * label
                else
                    label := "$(ion.label) × " * @sprintf("%.3g", Z) * label
                end
                linestyle --> :dash
                ylim --> (0.0, Inf)
                normalization --> Z
                ion, :density
            end
        end
    end

    # rotation
    if only === nothing || only == 3
        @series begin
            if only === nothing
                subplot := 3
            end
            title := "Rotation"
            label := "" * label
            if !ismissing(cpt, :rotation_frequency_tor_sonic)
                cpt, :rotation_frequency_tor_sonic
            else
                [NaN], [NaN]
            end
        end
    end
end

# =============== #
# solid_mechanics #
# =============== #
@recipe function plot_solid_mechanics(
    stress::Union{IMAS.solid_mechanics__center_stack__stress__hoop,IMAS.solid_mechanics__center_stack__stress__radial,IMAS.solid_mechanics__center_stack__stress__vonmises}
)
    smcs = parent(parent(stress))

    r_oh = smcs.grid.r_oh
    r_tf = smcs.grid.r_tf
    r_pl = getproperty(smcs.grid, :r_pl, missing)

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    @series begin
        r_oh, stress.oh ./ 1E6
    end
    @series begin
        primary := false
        r_tf, stress.tf ./ 1E6
    end
    if r_pl !== missing
        @series begin
            primary := false
            r_pl, stress.pl ./ 1E6
        end
    end

    if typeof(stress) <: IMAS.solid_mechanics__center_stack__stress__vonmises
        @series begin
            linestyle := :dash
            linewidth := 0.75 * base_linewidth
            label := "Yield strength"
            r_oh, r_oh .* 0.0 .+ smcs.properties.yield_strength.oh / 1E6
        end
        @series begin
            primary := false
            linestyle := :dash
            linewidth := 0.75 * base_linewidth
            r_tf, r_tf .* 0.0 .+ smcs.properties.yield_strength.tf / 1E6
        end
        if r_pl !== missing
            @series begin
                primary := false
                linestyle := :dash
                linewidth := 0.75 * base_linewidth
                r_pl, r_pl .* 0.0 .+ smcs.properties.yield_strength.pl / 1E6
            end
        end
    end
end

@recipe function plot_solid_mechanics(stress::IMAS.solid_mechanics__center_stack__stress)
    base_linewidth = get(plotattributes, :linewidth, 1.0)

    legend_position --> :bottomleft
    ylabel --> "Stresses [MPa]"
    xlabel --> "Radius [m]"
    title --> "Center Stack stresses"

    center_stack = parent(stress)

    config = []
    if center_stack.bucked == 1
        push!(config, "bucked")
    end
    if center_stack.noslip == 1
        push!(config, "noslip")
    end
    if center_stack.plug == 1
        push!(config, "plug")
    end
    if length(config) > 0
        config = " ($(join(config,",")))"
    else
        config = ""
    end

    @series begin
        label := "Von Mises" * config
        linewidth := base_linewidth * 3
        linestyle := :solid
        stress.vonmises
    end
    @series begin
        label := "Hoop" * config
        linewidth := base_linewidth * 2
        linestyle := :dash
        stress.hoop
    end
    @series begin
        label := "Radial" * config
        linewidth := base_linewidth * 2
        linestyle := :dashdot
        stress.radial
    end

    smcs = parent(stress)
    if smcs !== nothing
        for radius in (smcs.grid.r_oh[1], smcs.grid.r_oh[end], smcs.grid.r_tf[1], smcs.grid.r_tf[end])
            @series begin
                seriestype := :vline
                label := ""
                linestyle := :dash
                color := :black
                [radius]
            end
        end
        dd = top_dd(smcs)
        if dd!== nothing
            r_nose = smcs.grid.r_tf[1] + (smcs.grid.r_tf[end] - smcs.grid.r_tf[1]) * dd.build.tf.nose_hfs_fraction
            @series begin
                seriestype := :vline
                label := ""
                linestyle := :dash
                color := :black
                linewidth := base_linewidth / 2.0
                [r_nose]
            end
        end
    end
end

# ================ #
# balance_of_plant #
# ================ #
@recipe function plot_balance_of_plant(bop::IMAS.balance_of_plant)
    base_linewidth = plotattributes[:linewidth]

    size --> (800, 600)
    legend_position --> :outertopright
    ylabel --> "Electricity [Watts Electric]"
    xlabel --> "Time [s]"

    @series begin
        label := "Net electric"
        linewidth := base_linewidth * 3
        color := "Black"
        bop, :power_electric_net
    end

    @series begin
        label := "Electricity generated"
        linewidth := base_linewidth * 2
        linestyle --> :dash
        color := "Black"
        bop.power_plant, :power_electric_generated
    end

    for sys in bop.power_electric_plant_operation.system
        @series begin
            label := string(sys.name)
            sys, :power
        end
    end
end

# ========== #
# neutronics #
# ========== #
@recipe function plot_neutron_wall_loading_cx(nwl::IMAS.neutronics__time_slice___wall_loading; cx=true, component=:norm)
    id = plot_help_id(nwl)
    assert_type_and_record_argument(id, Symbol, "Component to plot (:norm, :r, :z, :power)"; component)
    assert_type_and_record_argument(id, Bool, "Plot cross section"; cx)

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    neutronics = top_ids(nwl)
    title = "Wall flux"
    units = "[W/m²]"
    if component == :norm
        data = sqrt.(nwl.flux_r .^ 2.0 .+ nwl.flux_z .^ 2.0)
    elseif component == :r
        data = nwl.flux_r
    elseif component == :z
        data = nwl.flux_z
    elseif component == :power
        data = nwl.power
        title = "Wall neutron power"
        units = "[W]"
    end
    data = data .+ 1.0

    wall_r = (neutronics.first_wall.r[1:end-1] .+ neutronics.first_wall.r[2:end]) ./ 2.0
    wall_z = (neutronics.first_wall.z[1:end-1] .+ neutronics.first_wall.z[2:end]) ./ 2.0

    if cx
        seriestype --> :path
        line_z --> log10.(data)
        aspect_ratio := :equal
        linewidth := 8 * base_linewidth
        label --> ""
        if component in (:r, :z)
            linecolor --> :seismic
        end
        colorbar_title --> "log₁₀(q neutron $units)"
        HF = neutronics.first_wall
        xlim --> [0.95 * minimum(HF.r) - 0.05 * (maximum(HF.r)), 1.05 * maximum(HF.r) - 0.05 * (minimum(HF.r))]
        ylim --> [1.05 * minimum(HF.z) - 0.05 * (maximum(HF.z)), 1.05 * maximum(HF.z) - 0.05 * (minimum(HF.z))]
        wall_r, wall_z

    else
        # sort data so that we start from the outer midplane
        index = vcat(argmax(wall_r)+1:length(wall_r), 1:argmax(wall_r))
        wall_r = wall_r[index]
        wall_z = wall_z[index]

        d = sqrt.(IMAS.gradient(neutronics.first_wall.r) .^ 2.0 .+ IMAS.gradient(neutronics.first_wall.z) .^ 2.0)
        d = sqrt.(diff(neutronics.first_wall.r) .^ 2.0 .+ diff(neutronics.first_wall.z) .^ 2.0)
        l = cumsum(d)
        l .-= l[1]

        xlabel --> "Clockwise distance along wall [m]"
        ylabel --> "$title $units"
        @series begin
            label --> "neutrons"
            ylabel --> "Wall flux [W/m²]"
            yscale --> :log10
            l, data[index]
        end
    end
end

# ==== #
# wall #
# ==== #
@recipe function plot_wall(wall::IMAS.wall)
    aspect_ratio := :equal
    @series begin
        return wall.description_2d
    end
end

@recipe function plot_wd2d_list(wd2d_list::IDSvector{<:IMAS.wall__description_2d})
    for wd2d in wd2d_list
        @series begin
            return wd2d
        end
    end
end

@recipe function plot_wd2d(wd2d::IMAS.wall__description_2d; show_limiter=true, show_vessel=true)
    id = plot_help_id(wd2d)
    assert_type_and_record_argument(id, Bool, "Show limiter units"; show_limiter)
    assert_type_and_record_argument(id, Bool, "Show vessel units"; show_vessel)
    if show_limiter
        @series begin
            return wd2d.limiter
        end
    end
    if show_vessel
        @series begin
            return wd2d.vessel
        end
    end
end

@recipe function plot_wd2dl(wd2dl::IMAS.wall__description_2d___limiter)
    @series begin
        lw := 2
        return wd2dl.unit
    end
end

@recipe function plot_wd2dv(wd2dv::IMAS.wall__description_2d___vessel)
    @series begin
        return wd2dv.unit
    end
end

# --- Vessel

@recipe function plot_wd2dvu_list(wd2dvu_list::IDSvector{<:IMAS.wall__description_2d___vessel__unit})
    for wd2dvu in wd2dvu_list
        @series begin
            return wd2dvu
        end
    end
end

@recipe function plot_wd2dvu(wd2dvu::IMAS.wall__description_2d___vessel__unit)
    if !isempty(wd2dvu.annular)
        @series begin
            label --> getproperty(wd2dvu, :name, "")
            wd2dvu.annular
        end
    end
    if !isempty(wd2dvu.element)
        @series begin
            label --> getproperty(wd2dvu, :name, "")
            wd2dvu.element
        end
    end
end

@recipe function plot_wd2dvua(annular::IMAS.wall__description_2d___vessel__unit___annular)
    if !ismissing(annular, :thickness)
        tmp = closed_polygon(annular.centreline.r, annular.centreline.z, Bool(annular.centreline.closed), annular.thickness)
        r, z = tmp.rz
        thickness = tmp.args[1]

        for i in 1:length(r)-1
            pts = thick_line_polygon(r[i], z[i], r[i+1], z[i+1], thickness[i], thickness[i+1])
            # Extract x and y coordinates for plotting
            xs, ys = map(collect, zip(pts...))
            @series begin
                primary := i == 1
                aspect_ratio := :equal
                seriestype := :shape
                linewidth := 0.0
                xs, ys
            end
        end
        @series begin
            primary := false
            color := :black
            lw := 0.2
            linestyle := :dash
            r, z
        end

    else
        if Bool(annular.outline_inner.closed)
            poly = join_outlines(annular.outline_inner.r, annular.outline_inner.z, annular.outline_outer.r, annular.outline_outer.z)
        else
            poly = ([reverse(annular.outline_inner.r); annular.outline_outer.r], [reverse(annular.outline_inner.z); annular.outline_outer.z])
        end
        poly = closed_polygon(poly...).rz
        @series begin
            aspect_ratio := :equal
            seriestype := :shape
            linewidth := 0.0
            poly[1], poly[2]
        end
        @series begin
            primary := false
            color := :black
            lw := 0.2
            poly
        end
    end
end

@recipe function plot_wd2dvue(elements::IDSvector{<:IMAS.wall__description_2d___vessel__unit___element})
    for (k, element) in enumerate(elements)
        @series begin
            primary := k == 1
            element
        end
    end
end

@recipe function plot_wd2dvue(element::IMAS.wall__description_2d___vessel__unit___element)
    aspect_ratio := :equal
    seriestype := :shape
    @series begin
        element.outline.r, element.outline.z
    end
end

# --- Limiter

@recipe function plot_wd2dlu_list(wd2dlu_list::IDSvector{<:IMAS.wall__description_2d___limiter__unit})
    for wd2dlu in wd2dlu_list
        @series begin
            return wd2dlu
        end
    end
end

@recipe function plot_wd2dlu(wd2dlu::IMAS.wall__description_2d___limiter__unit)
    @series begin
        label --> getproperty(wd2dlu, :name, "")
        wd2dlu.outline
    end
end

@recipe function plot_limiter_unit_outline(outline::IMAS.wall__description_2d___limiter__unit___outline)
    poly = closed_polygon(outline.r, outline.z)
    aspect_ratio := :equal
    @series begin
        poly.r, poly.z
    end
end

#= ============== =#
#  pulse_schedule  #
#= ============== =#
@recipe function plot_ps(ps::IMAS.pulse_schedule; time0=global_time(ps), simulation_start=nothing)
    id = plot_help_id(ps)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Union{Nothing,Float64}, "Simulation start time"; simulation_start)
    plots = []
    for (loc, (ids, field)) in filled_ids_fields(ps; eval_expr=true)
        if field != :reference
            continue
        end
        path = collect(f2p(ids))
        if "boundary_outline" in path
            continue
        end

        time_value = time_array_parent(ids)
        data_value = getproperty(ids, :reference)

        if length(collect(filter(x -> !isinf(x), time_value))) == 1
            continue
        end

        plt = Dict()
        plt[:ids] = ids
        plt[:x] = time_value
        plt[:y] = data_value
        remove = ("antenna", "flux_control", "beam", "unit", "density_control", "position_control")
        substitute =
            ("deposition_" => "",
                "_launched" => "",
                "rho_tor_norm_width" => "width",
                "rho_tor_norm" => "ρ",
                "pellet.launcher" => "pellet",
                "." => " ",
                "[" => " ",
                "]" => " ")
        plt[:label] = replace(p2i(filter(x -> x ∉ remove, path[2:end])), substitute...)
        push!(plots, plt)
    end

    layout := RecipesBase.@layout [length(plots) + 1]
    size --> (1000, 1000)

    @series begin
        subplot := 1
        label := "$(time0) [s]"
        aspect_ratio := :equal
        time0 := time0
        ps.position_control
    end

    eqt = try
        top_dd(ps).equilibrium.time_slice[time0]
    catch
        nothing
    end
    if eqt !== nothing
        @series begin
            subplot := 1
            legend := false
            cx := true
            alpha := 0.2
            color := :gray
            eqt
        end
    end

    wall = try
        top_dd(ps).wall
    catch
        nothing
    end
    if wall !== nothing
        @series begin
            show_limiter := true
            show_vessel := false
            subplot := 1
            legend := false
            wall
        end
    end

    # plotting at infinity does not work
    tmax = -Inf
    tmin = Inf
    for plt in plots
        xx = filter(x -> !isinf(x), plt[:x])
        tmax = max(tmax, maximum(xx))
        tmin = min(tmax, minimum(xx))
    end

    for (k, plt) in enumerate(plots)
        x = plt[:x]
        x[x.==Inf] .= tmax
        x[x.==-Inf] .= tmin
        y = plt[:y]
        n = min(length(x), length(y))
        @series begin
            subplot := k + 1
            label := ""
            x[1:n], y[1:n]
        end
        if simulation_start !== nothing
            @series begin
                subplot := k + 1
                primary := false
                seriestype := :vline
                linestyle := :dash
                [simulation_start]
            end
        end
        @series begin
            subplot := k + 1
            seriestype := :scatter
            primary := false
            marker := :circle
            markerstrokewidth := 0.0
            titlefontsize := 10
            title := nice_field(plt[:label])
            xlim := (tmin, tmax)
            [time0], interp1d(x[1:n], y[1:n]).([time0])
        end
    end
end

@recipe function plot_pc_time(pc::IMAS.pulse_schedule__position_control; time0=global_time(pc))
    id = plot_help_id(pc)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    aspect_ratio := :equal
    bnd = boundary(pc; time0)
    @series begin
        label := ""
        bnd.r, bnd.z
    end
    @series begin
        seriestype := :scatter
        marker --> :star
        markerstrokewidth --> 0
        primary := false
        Xs = x_points(pc.x_point; time0)
        [x[1] for x in Xs if x[1] != 0.0], [x[2] for x in Xs if x[1] != 0.0]
    end
    @series begin
        seriestype := :scatter
        marker --> :circle
        markerstrokewidth --> 0
        primary := false
        Ss = strike_points(pc.strike_point; time0)
        [x[1] for x in Ss if x[1] != 0.0], [x[2] for x in Ss if x[1] != 0.0]
    end
    if !ismissing(pc.magnetic_axis.r, :reference) && !ismissing(pc.magnetic_axis.z, :reference)
        @series begin
            seriestype := :scatter
            marker --> :cross
            markerstrokewidth --> 0
            primary := false
            RA = interp1d(pc.time, pc.magnetic_axis.r.reference).(time0)
            ZA = interp1d(pc.time, pc.magnetic_axis.z.reference).(time0)
            [RA], [ZA]
        end
    end
end

#= =========== =#
#  controllers  #
#= =========== =#
@recipe function plot_controllers(controller_outputs::controllers__linear_controller___outputs)
    controller = IMAS.parent(controller_outputs)

    data = controller.inputs.data[1, :]
    integral = cumsum((data[k+1] + data[k]) / 2.0 * (controller.inputs.time[k+1] - controller.inputs.time[k]) for k in eachindex(controller.inputs.time) - 1)
    derivative = diff(data) ./ diff(controller.inputs.time)

    @series begin
        label := "PID"
        color := :black
        lw := 2
        controller.outputs.time, controller.outputs.data[1, :]
    end

    @series begin
        label := "P=$(controller.pid.p.data[1])"
        controller.outputs.time, data .* controller.pid.p.data[1]
    end

    @series begin
        label := "I=$(controller.pid.i.data[1])"
        controller.outputs.time[2:end], integral .* controller.pid.i.data[1]
    end

    @series begin
        label := "D=$(controller.pid.d.data[1])"
        controller.outputs.time[2:end], derivative .* controller.pid.d.data[1]
    end

    @series begin
        label := "PI"
        controller.outputs.time[2:end], data[2:end] .* controller.pid.p.data[1] .+ integral .* controller.pid.i.data[1]
    end

    @series begin
        label := "PID"
        controller.outputs.time[2:end], data[2:end] .* controller.pid.p.data[1] .+ integral .* controller.pid.i.data[1] .+ derivative .* controller.pid.d.data[1]
    end

end

@recipe function plot_controllers(controller_inputs::controllers__linear_controller___inputs)
    controller = IMAS.parent(controller_inputs)
    label = "P=$(controller.pid.p.data[1]), I=$(controller.pid.i.data[1]), D=$(controller.pid.d.data[1])"
    @series begin
        label := "$(controller.input_names[1]) ($label)"
        controller.inputs.time, controller.inputs.data[1, :]
    end
end

@recipe function plot_controllers(controller::controllers__linear_controller)
    layout := (2, 1)
    size --> (800, 600)
    margin --> 5 * Measures.mm
    @series begin
        subplot := 1
        controller.inputs
    end
    @series begin
        subplot := 2
        controller.outputs
    end
end

#= ======= =#
#  summary  #
#= ======= =#
@recipe function plot(summary::IMASdd.summary)
    valid_leaves = [leaf for leaf in IMASdd.AbstractTrees.Leaves(summary) if typeof(leaf.value) <: Vector && leaf.field != :time]
    layout := length(valid_leaves)
    N = Int(ceil(sqrt(length(valid_leaves))))
    size --> (200 * N, 200 * N)
    for (k, leaf) in enumerate(valid_leaves)
        loc = replace(location(leaf.ids, leaf.field), ".value" => "")
        _, title, label = rsplit(loc, "."; limit=3)
        @series begin
            guidefontvalign := :top
            titlefont --> font(8, "Arial", "bold")
            subplot := k
            title := replace("$title.$label", "global_quantities." => "")
            label := ""
            xlabel := ""
            leaf.ids, leaf.field
        end
    end
end

#= ================ =#
#  generic plotting  #
#= ================ =#
@recipe function plot_field(ids::IMAS.IDS, field::Symbol; normalization=1.0, coordinate=nothing, weighted=nothing, fill0=false)
    id = plot_help_id(ids, field)
    @assert hasfield(typeof(ids), field) "$(location(ids)) does not have field `$field`. Did you mean: $(keys(ids))"
    assert_type_and_record_argument(id, Union{Real,AbstractVector{<:Real}}, "Normalization factor"; normalization)
    assert_type_and_record_argument(id, Union{Nothing,Symbol}, "Coordinate for x-axis"; coordinate)
    assert_type_and_record_argument(id, Union{Nothing,Symbol}, "Weighting field"; weighted)
    assert_type_and_record_argument(id, Bool, "Fill area under curve"; fill0)

    coords = coordinates(ids, field; coord_leaves=[coordinate])
    coordinate_name = coords.names[1]
    coordinate_value = coords.values[1]

    # If the field is the reference coordinate of the given IDS,
    # set the coordinate_value as its index
    if coordinate_name == "1...N"
        coordinate_value = 1:length(getproperty(ids, field))
    end

    xvalue = coordinate_value
    yvalue = getproperty(ids, field) .* normalization

    @series begin
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)

        if endswith(coordinate_name, "_norm")
            xlim --> (0.0, 1.0)
        end

        # multiply y by things like `:area` or `:volume`
        if weighted !== nothing
            weight = coordinates(ids, field; coord_leaves=[weighted])
            yvalue .*= weight.values[1]
            ylabel = nice_units(units(ids, field) * "*" * units(weight.names[1]))
            label = nice_field("$field*$weighted")
        else
            ylabel = nice_units(units(ids, field))
            label = nice_field(field)
        end

        if coordinate_name == "1...N"
            xlabel --> "index"
        else
            xlabel --> nice_field(i2p(coordinate_name)[end]) * nice_units(units(coordinate_name))
        end
        ylabel --> ylabel
        label --> label

        # plot 1D Measurements with ribbon
        if (eltype(yvalue) <: Measurement) && !((eltype(xvalue) <: Measurement))
            ribbon := Measurements.uncertainty.(yvalue)
            yvalue = Measurements.value.(yvalue)
        end

        xvalue, yvalue
    end
    if fill0
        @series begin
            label := ""
            primary := false
            alpha := 0.25
            fillrange := yvalue
            xvalue, yvalue .* 0.0
        end
    end
end

@recipe function plot(x::AbstractVector{<:Real}, y::AbstractVector{<:Measurement}, err::Symbol=:ribbon)
    if err == :ribbon
        ribbon := Measurements.uncertainty.(y)
    elseif err == :bar
        yerror := Measurements.uncertainty.(y)
    end
    return x, Measurements.value.(y)
end

#= ============= =#
#  time plotting  #
#= ============= =#
# This plot recipe handles time dependent quantities that are shorter than their time coordinate
@recipe function f(::Type{Val{:time}}, plt::AbstractPlot)
    x, y = plotattributes[:x], plotattributes[:y]

    # Determine the length of the shorter vector
    min_length = min(length(x), length(y))

    # Truncate x and y to the length of the shorter vector
    x = x[1:min_length]
    y = y[1:min_length]

    # plot
    @series begin
        if sum(.!isnan.(y .* x)) <= 1
            markerstrokewidth := 0.0
            seriestype := :scatter
        else
            seriestype := :path
        end
        x, y
    end
end

function latex_support()
    return Plots.backend_name() in [:gr, :pgfplotsx, :pythonplot, :inspectdr]
end

#= ================== =#
#  handling of labels  #
#= ================== =#
nice_field_symbols = Dict()
nice_field_symbols["rho_tor_norm"] = () -> latex_support() ? L"\rho" : "ρ"
nice_field_symbols["psi"] = () -> latex_support() ? L"\psi" : "ψ"
nice_field_symbols["psi_norm"] = () -> latex_support() ? L"\psi_\mathrm{N}" : "ψₙ"
nice_field_symbols["rotation_frequency_tor_sonic"] = "Rotation"
nice_field_symbols["i_plasma"] = "Plasma current"
nice_field_symbols["b_field_tor_vacuum_r"] = "B₀×R₀"
nice_field_symbols["rho_tor_norm_width"] = "w₀"
nice_field_symbols["geometric_axis.r"] = "Rgeo"
nice_field_symbols["geometric_axis.z"] = "Zgeo"

function nice_field(field::AbstractString)
    if field in keys(nice_field_symbols)
        field = nice_field_symbols[field]
        return field isa Function ? field() : field
    else
        field = replace(field,
            r"n_e" => "nₑ",
            r"n_i" => "nᵢ",
            r"T_e" => "Tₑ",
            r"T_i" => "Tᵢ",
            r"_tor" => " toroidal",
            r"_greenwald" => " GW",
            r"_pol" => " poloidal",
            r"\bip\b" => "Ip",
            r"\bec\b" => "EC",
            r"\bic\b" => "IC",
            r"\blh\b" => "LH",
            r"\bnb\b" => "NB",
            r"\bnbi\b" => "NBI",
            "_" => " ")
        if length(split(field, " ")[1]) > 2
            field = uppercasefirst(field)
        end
    end
    return field
end

function nice_field(field::Symbol)
    return nice_field(String(field))
end

function nice_units(units::String)
    if units == "-"
        units = ""
    end
    if length(units) > 0
        units = replace(units, r"\^([-+]?[0-9]+)" => s"^{\1}")
        units = replace(units, "." => s"\\,")
        units = latex_support() ? L"[\mathrm{%$units}]" : "[$units]"
        units = " " * units
    end
    return units
end

"""
    hash_to_color(input::Any; seed::Int=0)

Generate a unique RGB color based on the hash of an input, with an optional seed for color adjustment.

This function computes the hash of the given `input` and converts it to an RGB color in the [0,1] range.
The optional `seed` parameter shifts the hash to allow different color mappings for the same input.
Using a different `seed` value will produce a unique color set for the same input.
"""
function hash_to_color(input::Any; seed::Int=0)
    # Generate a hash for the input
    h = hash(hash(input) + seed)

    # Convert the hash to a range for RGB values (0 to 1)
    r = ((h >> 16) & 0xFF) / 255.0
    g = ((h >> 8) & 0xFF) / 255.0
    b = (h & 0xFF) / 255.0

    # Return an RGB color using Plots' RGB type
    return RGB(r, g, b)
end

# To fix a vscode's bug w.r.t UnicodePlots
# (see https://github.com/JuliaPlots/Plots.jl/issues/4956  for more details)
Base.showable(::MIME"image/png", ::Plots.Plot{Plots.UnicodePlotsBackend}) = applicable(UnicodePlots.save_image, devnull)

