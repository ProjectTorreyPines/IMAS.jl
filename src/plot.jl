document[:Plot] = Symbol[]

import PlotUtils
using RecipesBase
using LaTeXStrings
import Measures
import Graphs
using GraphRecipes
import HelpPlots
using HelpPlots

function HelpPlots.recipe_dispatch(arg::IDS)
    return ulocation(arg)
end

#= =========== =#
#  PlotExtrema  #
#= =========== =#
mutable struct PlotExtrema
    yext::Vector{Float64} # the extrema of y that are observed
    ylim::Vector{Float64} # the y limits of the plot
    active::Bool # should the y limits be touched
end

function PlotExtrema()
    return PlotExtrema([Inf, -Inf], [-Inf, Inf], false)
end

"""
    update_limits!(plot_extrema::PlotExtrema, y::Vector{T}, x::Vector{T}) where {T<:Real}

Update plot limits based on y value
"""
function update_limits!(plot_extrema::PlotExtrema, y::Vector{T}, x::Vector{T}) where {T<:Real}
    irho = argmin_abs(x, 0.85)

    if maximum(abs.(@views y[irho:end])) / (1 + maximum(abs.(@views y[1:irho]))) > 3
        plot_extrema.active = true
    end

    return update_limits!(plot_extrema, @views y[1:irho])
end

function update_limits!(plot_extrema::PlotExtrema, y::AbstractVector{T}) where {T<:Real}
    plot_extrema.yext[1] = min(plot_extrema.yext[1], minimum(y))
    plot_extrema.yext[2] = max(plot_extrema.yext[2], maximum(y))

    if plot_extrema.yext[1] < Inf
        plot_extrema.ylim[1] = plot_extrema.yext[1]
    end
    if plot_extrema.yext[2] > -Inf
        plot_extrema.ylim[2] = plot_extrema.yext[2]
    end

    if any(map(isinf, plot_extrema.ylim))
        plot_extrema.active = false
    else
        delta = plot_extrema.ylim[2] - plot_extrema.ylim[1]
        plot_extrema.ylim[1] -= delta / 30
        plot_extrema.ylim[2] += delta / 30
    end

    return nothing
end

#= ========= =#
#  pf_active  #
#= ========= =#
@recipe function plot_pf_active_cx(pfa::pf_active; what=:cx, time0=global_time(pfa), cname=:vik)
    id = recipe_dispatch(pfa)
    assert_type_and_record_argument(id, Symbol, "What to plot [:currents, :cx, :coils_flux]"; what)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Symbol, "Colormap name"; cname)

    currents = [ismissing(c.current, :data) ? 0.0 : get_time_array(c.current, :data, time0) * getproperty(c.element[1], :turns_with_sign, 1.0) for c in pfa.coil]
    if isempty(pfa.coil) || ismissing(pfa.coil[1].current, :time) || isempty(pfa.coil[1].current.time) || time0 < pfa.coil[1].current.time[1]
        c_unit = "A"
    else
        CURRENT = 0.0
        for c in pfa.coil
            if !ismissing(c.current, :data)
                if time0 == -Inf && c.current.time[1] == -Inf
                    CURRENT = max(CURRENT, maximum(abs, c.current.data[1] * getproperty(c.element[1], :turns_with_sign, 1.0)))
                else
                    CURRENT = max(CURRENT, maximum(abs, c.current.data * getproperty(c.element[1], :turns_with_sign, 1.0)))
                end
            end
        end
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
    id = recipe_dispatch(coil)
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

@recipe function plot_pf_active_rail(rails::AbstractVector{<:IMAS.build__pf_active__rail})
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

@recipe function plot_circuits(circuits::AbstractVector{<:IMAS.pf_active__circuit})
    layout := RecipesBase.@layout length(circuits)
    for (kkk, circuit) in enumerate(circuits)
        @series begin
            subplot := kkk
            circuit
        end
    end
end

#= ========== =#
#  pf_passive  #
#= ========== =#
@recipe function plot_pf_passive(pf_passive::IMAS.pf_passive)
    @series begin
        pf_passive.loop
    end
end

@recipe function plot_pf_passive(loops::AbstractVector{<:IMAS.pf_passive__loop})
    for loop in loops
        @series begin
            loop
        end
    end
end

@recipe function plot_loop(loop::pf_passive__loop{T}; loop_names=false) where {T<:Real}
    id = recipe_dispatch(loop)
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

#= ======= =#
#  costing  #
#= ======= =#
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

#= =========== =#
#  equilibrium  #
#= =========== =#
@recipe function plot_eq(eq::IMAS.equilibrium; time0=global_time(eq))
    id = recipe_dispatch(eq)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    @series begin
        return eq.time_slice[time0]
    end
end

@recipe function plot_eqt(eqt::IMAS.equilibrium__time_slice; cx=false, coordinate=:psi_norm, core_profiles_overlay=false)
    id = recipe_dispatch(eqt)
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

    coordinate := coordinate

    @series begin
        if !cx
            subplot := 1
        end
        eqt.profiles_2d
    end

    if !cx
        if core_profiles_overlay
            cp = top_dd(eqt).core_profiles
            cp1d = top_dd(eqt).core_profiles.profiles_1d[argmin_abs(cp.time, eqt.time)]
        end

        # pressure
        if !ismissing(eqt.profiles_1d, :pressure)
            @series begin
                label := ""
                xlabel := ""
                subplot := 2
                normalization := 1E-6
                ylabel := ""
                title := "P [MPa]"
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
                    title := "P [MPa]"
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
                title := "Jtor [MA/m²]"
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
                    title := "Jtor [MA/m²]"
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
                    title := "ρ"
                    eqt.profiles_1d, :rho_tor_norm
                else
                    title := "Ψ [Wb]"
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
                title := "q"
                if eqt.profiles_1d.q[end] > 0.0
                    ylim := (0.0, 5)
                else
                    ylim := (-5.0, 0.0)
                end
                eqt.profiles_1d, :q
            end
        end
    end
end

@recipe function plot_eqt2dv(eqt2dv::AbstractVector{<:IMAS.equilibrium__time_slice___profiles_2d})
    if !isempty(eqt2dv)
        @series begin
            return eqt2dv[1]
        end
    end
end

function _circle(a)
    return [a * cos(t) for t in range(0, 2π, 101)], [a * sin(t) for t in range(0, 2π, 101)]
end

@recipe function plot_eqt2(
    eqt2d::IMAS.equilibrium__time_slice___profiles_2d;
    levels_in=11,
    levels_out=11,
    show_secondary_separatrix=false,
    show_x_points=true,
    show_strike_points=true,
    show_magnetic_axis=true,
    coordinate=:psi,
    top=false)

    id = recipe_dispatch(eqt2d)
    assert_type_and_record_argument(id, Union{Int,AbstractVector{<:Real}}, "Levels inside LCFS"; levels_in)
    assert_type_and_record_argument(id, Union{Int,AbstractVector{<:Real}}, "Levels outside LCFS"; levels_out)
    assert_type_and_record_argument(id, Symbol, "X coordinate for 2D contours"; coordinate)
    assert_type_and_record_argument(id, Bool, "Plot secondary separatrix"; show_secondary_separatrix)
    assert_type_and_record_argument(id, Bool, "Show X points"; show_x_points)
    assert_type_and_record_argument(id, Bool, "Show strike points"; show_strike_points)
    assert_type_and_record_argument(id, Bool, "Show magnetic axis"; show_magnetic_axis)
    assert_type_and_record_argument(id, Bool, "Top view"; top)

    label --> ""
    aspect_ratio := :equal
    @series begin
        [], []
    end
    primary --> false

    eqt = parent(parent(eqt2d))

    # handle levels
    x_coord = getproperty(eqt.profiles_1d, coordinate)
    boundary_level = x_coord[end]
    axis_level = x_coord[1]
    if typeof(levels_in) <: Int
        if levels_in == 1
            levels_in = [boundary_level]
        elseif levels_in > 0
            levels_in = range(axis_level, boundary_level, levels_in)
        else
            levels_in = []
        end
    end
    delta_psi = (boundary_level - axis_level)
    if typeof(levels_out) <: Int
        if levels_out > 0
            levels_out = delta_psi / length(levels_in) .* collect(0:levels_out) .+ boundary_level
        else
            levels_out = []
        end
    end
    levels = unique(vcat(levels_in, levels_out))

    if coordinate == :psi
        psi_levels_in = levels_in
        psi_levels_out = levels_out
        psi_levels = levels
        psi__boundary_level = boundary_level
        psi__axis_level = axis_level
    else
        psi_interp = interp1d(x_coord, eqt.profiles_1d.psi)
        psi_levels_in = psi_interp.(levels_in)
        psi_levels_out = psi_interp.(levels_out)
        psi_levels = psi_interp.(levels)
        psi__boundary_level = psi_interp.(boundary_level)
        psi__axis_level = psi_interp.(axis_level)
    end

    fw = first_wall(top_dd(eqt).wall)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    if top
        xlabel --> "X [m]"
        ylabel --> "Y [m]"
        @series begin
            ls := :dash
            _circle(RA)
        end
        for psi_level in psi_levels_in
            for (pr, pz) in flux_surface(eqt, psi_level, :closed, fw.r, fw.z)
                for f in (minimum, maximum)
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
                        _circle(f(pr))
                    end
                end
            end
        end
    else
        xlabel --> "R [m]"
        ylabel --> "Z [m]"
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

@recipe function plot_x_points(x_points::AbstractVector{<:IMAS.equilibrium__time_slice___boundary__x_point})
    for (k, x_point) in enumerate(x_points)
        @series begin
            markersize --> ((length(x_points) + 1 - k) / length(x_points)) * 4
            x_point
        end
    end
end

@recipe function plot_x_point(x_point::IMAS.equilibrium__time_slice___boundary__x_point)
    @series begin
        seriestype := :scatter
        marker --> :star
        markerstrokewidth --> 0
        label --> ""
        aspect_ratio := :equal
        [(x_point.r, x_point.z)]
    end
end

@recipe function plot_strike_points(s_points::AbstractVector{<:IMAS.equilibrium__time_slice___boundary__strike_point})
    for s_point in s_points
        @series begin
            s_point
        end
    end
end

@recipe function plot_strike_point(s_point::IMAS.equilibrium__time_slice___boundary__strike_point)
    @series begin
        seriestype := :scatter
        marker --> :circle
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

#= ========= =#
#  divertors  #
#= ========= =#
@recipe function plot_divertors(divertors::IMAS.divertors)
    @series begin
        divertors.divertor
    end
end

@recipe function plot_divertors(divertors::AbstractVector{<:IMAS.divertors__divertor})
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

@recipe function plot_divertors(targets::AbstractVector{<:IMAS.divertors__divertor___target})
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

@recipe function plot_divertors(tiles::AbstractVector{<:IMAS.divertors__divertor___target___tile})
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

#= ===== =#
#  build  #
#= ===== =#
function join_outlines(r1::AbstractVector{T}, z1::AbstractVector{T}, r2::AbstractVector{T}, z2::AbstractVector{T}) where {T<:Real}
    distance, i1, i2 = minimum_distance_polygons_vertices(r1, z1, r2, z2)
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

@recipe function plot_build(bd::IMAS.build; equilibrium=true, pf_active=true, pf_passive=true, wall=true)
    id = recipe_dispatch(bd)
    assert_type_and_record_argument(id, Bool, "Include plot of dd.equilibrium"; equilibrium)
    assert_type_and_record_argument(id, Bool, "Include plot of dd.pf_active"; pf_active)
    assert_type_and_record_argument(id, Bool, "Include plot of dd.pf_passive"; pf_passive)
    assert_type_and_record_argument(id, Bool, "Include plot of dd.wall"; wall)

    dd = top_dd(bd)

    legend_position --> :outerbottomright
    aspect_ratio := :equal
    grid --> :none

    if !isempty(bd.layer)
        rmax = maximum(bd.layer[end].outline.r)
        xlim --> [0, rmax]
    else
        xlim --> [0, Inf]
    end

    if dd !== nothing && equilibrium
        @series begin
            cx := true
            dd.equilibrium
        end
    end

    if isempty(bd.layer)
        if dd !== nothing && wall
            fw = IMAS.first_wall(dd.wall)
            @series begin
                label := ""
                color := :black
                fw.r, fw.z
            end
        end
    else
        @series begin
            bd.layer
        end
    end

    if dd !== nothing && pf_active
        @series begin
            colorbar --> :false
            dd.pf_active
        end
    end

    if dd !== nothing && pf_passive
        @series begin
            colorbar --> :false
            color := :yellow
            dd.pf_passive
        end
    end

end

@recipe function plot_build_layer(layers::AbstractVector{<:IMAS.build__layer}; cx=true, wireframe=false, only=Symbol[], exclude_layers=Symbol[])
    id = recipe_dispatch(layers)
    assert_type_and_record_argument(id, Bool, "Plot cross section"; cx)
    assert_type_and_record_argument(id, Bool, "Use wireframe"; wireframe)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "Only include certain layers"; only)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "Exclude certain layers"; exclude_layers)

    bd = parent(layers)

    base_linewidth = get(plotattributes, :linewidth, 1.0)

    legend_position --> :outerbottomright
    aspect_ratio := :equal
    grid --> :none

    # cx
    if cx && !isempty(layers)
        rmax = maximum(layers[end].outline.r)

        # everything after first vacuum in _out_
        if (isempty(only) || (:cryostat in only)) && :cryostat ∉ exclude_layers
            for k in get_build_indexes(layers; fs=_out_)[2:end]
                if !wireframe
                    @series begin
                        seriestype --> :shape
                        linewidth := 0.0
                        color --> :lightgray
                        label --> (!wireframe ? "Cryostat" : "")
                        xlim --> [0, rmax]
                        join_outlines(
                            layers[k].outline.r,
                            layers[k].outline.z,
                            layers[k-1].outline.r,
                            layers[k-1].outline.z
                        )
                    end
                end
                @series begin
                    seriestype --> :path
                    linewidth := base_linewidth
                    color --> :black
                    label --> ""
                    xlim --> [0, rmax]
                    layers[k].outline.r[1:end-1], layers[k].outline.z[1:end-1]
                end
            end
        end

        # first vacuum in _out_
        if !wireframe
            k = get_build_indexes(layers; fs=_out_)[1]
            @series begin
                seriestype --> :shape
                linewidth := 0.0
                color --> :white
                linecolor --> :white
                label --> ""
                xlim --> [0, rmax]
                join_outlines(
                    layers[k].outline.r,
                    layers[k].outline.z,
                    get_build_layer(layers; type=_plasma_).outline.r,
                    get_build_layer(layers; type=_plasma_).outline.z
                )
            end
            if (isempty(only) || (:cryostat in only)) && :cryostat ∉ exclude_layers
                @series begin
                    seriestype --> :path
                    linewidth := base_linewidth
                    color --> :black
                    label --> ""
                    xlim --> [0, rmax]
                    layers[k].outline.r[1:end-1], layers[k].outline.z[1:end-1]
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
            [0.0, 0.0], [minimum(layers[end].outline.z), maximum(layers[end].outline.z)]
        end

        # all layers inside of the TF
        if (isempty(only) || (:oh in only)) && (!(:oh in exclude_layers))
            for k in get_build_indexes(layers; fs=_in_)
                layer = layers[k]
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
        for k in get_build_indexes(layers; fs=_hfs_)
            l = layers[k]
            l1 = layers[k+1]
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
            plasma_outline = outline(get_build_layer(layers; type=_plasma_))
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
        for l in layers
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

@recipe function plot_build_layer(layer::IMAS.build__layer)
    @series begin
        label --> layer.name
        layer.outline
    end
end

@recipe function plot_build_layer_outline(outline::IMAS.build__layer___outline)
    @series begin
        aspect_ratio := :equal
        outline.r, outline.z
    end
end

@recipe function plot_build_structure(structure::IMAS.build__structure)
    @series begin
        aspect_ratio := :equal
        structure.outline
    end
end

@recipe function plot_build_structure_outline(outline::IMAS.build__structure___outline)
    @series begin
        aspect_ratio := :equal
        outline.r, outline.z
    end
end

#= ======== =#
#  build tf  #
#= ======== =#
@recipe function plot_build_tf(tf::IMAS.build__tf)
    layers = parent(tf).layer
    TF = get_build_layers(layers; type=IMAS._tf_)

    for n in 1:tf.coils_n
        x, y = top_view_outline(tf, n)
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

#= ========= =#
#  transport  #
#= ========= =#
@recipe function plot_ct1d(ct1d__electrons__energy::IMAS.core_transport__model___profiles_1d___electrons__energy, ::Val{:flux}, plot_extrema::PlotExtrema=PlotExtrema())
    @series begin
        markershape --> :none
        title := "Electron energy flux"
        label --> ""
        normalization := 1E-6
        ylabel --> "[MW/m²]"
        ct1d__electrons__energy, :flux
    end
end

@recipe function plot_ct1d(ct1d__total_ion_energy::IMAS.core_transport__model___profiles_1d___total_ion_energy, ::Val{:flux}, plot_extrema::PlotExtrema=PlotExtrema())
    @series begin
        markershape --> :none
        title := "Total ion energy flux"
        label --> ""
        normalization := 1E-6
        ylabel --> "[MW/m²]"
        ct1d__total_ion_energy, :flux
    end
end

@recipe function plot_ct1d(ct1d__electrons__particles::IMAS.core_transport__model___profiles_1d___electrons__particles, ::Val{:flux}, plot_extrema::PlotExtrema=PlotExtrema())
    @series begin
        update_limits!(plot_extrema, ct1d__electrons__particles.flux)
        if plot_extrema.active
            ylim --> plot_extrema.ylim
        end
        markershape --> :none
        title := "Electron particle flux"
        label --> ""
        ylabel --> "[s⁻¹/m²]"
        ct1d__electrons__particles, :flux
    end
end

@recipe function plot_ct1d(ct1d__ion___particles::IMAS.core_transport__model___profiles_1d___ion___particles, ::Val{:flux}, plot_extrema::PlotExtrema=PlotExtrema())
    ion = parent(ct1d__ion___particles)
    @series begin
        update_limits!(plot_extrema, ct1d__ion___particles.flux)
        if plot_extrema.active
            ylim --> plot_extrema.ylim
        end
        markershape --> :none
        title := "$(ion.label) particle flux"
        label --> ""
        ylabel --> "[s⁻¹/m²]"
        ct1d__ion___particles, :flux
    end
end

@recipe function plot_ct1d(ct1d__momentum_tor::IMAS.core_transport__model___profiles_1d___momentum_tor, ::Val{:flux}, plot_extrema::PlotExtrema=PlotExtrema())
    @series begin
        markershape --> :none
        title := "Momentum flux"
        label --> ""
        ylabel --> "[Nm/m²]"
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
    for ion in ions_list
        if ion in available_ions
            kion = findfirst(x -> x == ion, available_ions)
            push!(paths, [:ion, kion, :particles])
        else
            push!(paths, [nothing])
        end
    end
    return paths
end

@recipe function plot_ct1d(ct1d::IMAS.core_transport__model___profiles_1d, plots_extrema::Vector{PlotExtrema}=PlotExtrema[]; only=nothing, ions=Symbol[])
    id = recipe_dispatch(ct1d)
    assert_type_and_record_argument(id, Union{Nothing,Int}, "Plot only this subplot number"; only)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)

    paths = transport_channel_paths(ct1d.ion, ions)

    if only === nothing
        layout := 4 + length(ions)
        Nr = round(Int, sqrt(4 + length(ions)), RoundUp)
        size --> (1100, Nr * 290)
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)
    end

    for (k, path) in enumerate(paths)
        if length(plots_extrema) < k
            resize!(plots_extrema, k)
            plots_extrema[k] = PlotExtrema()
        end
        if k == only || only === nothing
            if path[end] !== nothing && !ismissing(ct1d, [path; :flux]) && sum(abs.(getproperty(goto(ct1d, path), :flux))) > 0.0
                @series begin
                    if only === nothing
                        subplot := k
                    end
                    goto(ct1d, path), Val(:flux), plots_extrema[k]
                end
            end
        end
    end
end

@recipe function plot_core_transport_model(
    model::IMAS.core_transport__model{T},
    plots_extrema::Vector{PlotExtrema}=PlotExtrema[];
    ions=Symbol[:my_ions],
    time0=global_time(model)
) where {T<:Real}
    if nearest_causal_time(time_array_from_parent_ids(model.profiles_1d, Val(:get)), time0; bounds_error=false).index > 0
        model_type = name_2_index(model)
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
            model.profiles_1d[time0], plots_extrema
        end
    end
end

@recipe function plot_core_transport(ct::IMAS.core_transport{T}; ions=Symbol[:my_ions], time0=global_time(ct)) where {T<:Real}
    id = recipe_dispatch(ct)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    dd = top_dd(ct)
    cp1d = dd.core_profiles.profiles_1d[time0]

    if ions == [:my_ions]
        ions = list_ions(ct, dd.core_profiles; time0)
    end

    model_type = name_2_index(ct.model)
    rhos = T[]
    for model in ct.model
        if model.identifier.index ∈ (model_type[k] for k in (:combined, :unspecified, :transport_solver, :unknown))
            continue
        end
        if !isempty(model.profiles_1d) && time0 >= model.profiles_1d[1].time
            ct1d = model.profiles_1d[time0]
            append!(rhos, ct1d.grid_flux.rho_tor_norm)
        end
    end
    rhos = unique(rhos)

    plots_extrema = PlotExtrema[]

    if dd !== nothing
        tot_source = IMAS.core_sources__source{T}()
        resize!(tot_source.profiles_1d, 1)
        merge!(tot_source.profiles_1d[1], total_sources(dd; time0))
        tot_source.identifier.index = 1
        tot_source.identifier.name = "Total source"
        @series begin
            linewidth := 2
            color := :blue
            flux := true
            show_zeros := true
            ions := ions
            tot_source, plots_extrema
        end
    end

    @series begin
        linewidth := 2
        color := :red
        label := "Total transport"
        ions := ions
        total_fluxes(ct, cp1d, rhos; time0), plots_extrema
    end

    for model in ct.model
        if model.identifier.index ∈ (model_type[k] for k in (:combined, :unspecified, :transport_solver, :unknown))
            continue
        end
        if !isempty(model.profiles_1d) && time0 >= model.profiles_1d[1].time
            @series begin
                time0 := time0
                ions := ions
                model, plots_extrema
            end
        end
    end
end

#= ============ =#
#  core_sources  #
#= ============ =#
@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:electrons__energy}, plot_extrema::PlotExtrema=PlotExtrema())
    return cs1d.electrons, Val(:energy), plot_extrema
end

@recipe function plot_source1d(cs1de::IMAS.core_sources__source___profiles_1d___electrons, v::Val{:energy}, plot_extrema::PlotExtrema=PlotExtrema();
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    min_power=0.0,
    only_positive_negative=0,
    show_source_number=false)

    id = recipe_dispatch(cs1de, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Float64, "Minimum power threshold"; min_power)
    assert_type_and_record_argument(id, Int, "Show only positive or negative values (0 for all)"; only_positive_negative)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    cs1d = parent(cs1de; error_parent_of_nothing=true)
    source = parent(parent(cs1d; error_parent_of_nothing=false); error_parent_of_nothing=false)
    name, identifier, idx = source_name_identifier(source, name, show_source_number)

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
        elseif identifier == :time_derivative
            linestyle --> :dashdot
        end
        color := idx
        title --> "Electron Energy"
        if show_condition
            label := "$name " * @sprintf("[%.3g MW]", tot / 1E6) * label
            if identifier in (:ec, :ic, :lh, :nbi, :pellet)
                fill0 --> true
            end
            if !ismissing(cs1de, :power_inside) && flux
                ylabel --> "[MW/m²]"
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

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:total_ion_energy}, plot_extrema::PlotExtrema=PlotExtrema();
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    min_power=0.0,
    only_positive_negative=0,
    show_source_number=false)

    id = recipe_dispatch(cs1d, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Float64, "Minimum power threshold"; min_power)
    assert_type_and_record_argument(id, Int, "Show only positive or negative values (0 for all)"; only_positive_negative)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    source = parent(parent(cs1d; error_parent_of_nothing=false); error_parent_of_nothing=false)
    name, identifier, idx = source_name_identifier(source, name, show_source_number)

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
        elseif identifier == :time_derivative
            linestyle --> :dashdot
        end
        color := idx
        title --> "Total Ion Energy"
        if show_condition
            label := "$name " * @sprintf("[%.3g MW]", tot / 1E6) * label
            if identifier in (:ec, :ic, :lh, :nbi, :pellet)
                fill0 --> true
            end
            if !ismissing(cs1d, :total_ion_power_inside) && flux
                ylabel --> "[MW/m²]"
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

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:electrons__particles}, plot_extrema::PlotExtrema=PlotExtrema())
    return cs1d.electrons, Val(:particles), plot_extrema
end

@recipe function plot_source1d(cs1de::IMAS.core_sources__source___profiles_1d___electrons, v::Val{:particles}, plot_extrema::PlotExtrema=PlotExtrema();
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    show_source_number=false)

    id = recipe_dispatch(cs1de, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    cs1d = parent(cs1de; error_parent_of_nothing=true)
    source = parent(parent(cs1d; error_parent_of_nothing=false); error_parent_of_nothing=false)
    name, identifier, idx = source_name_identifier(source, name, show_source_number)

    tot = 0.0
    if !ismissing(cs1de, :particles)
        tot = trapz(cs1d.grid.volume, cs1de.particles)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        if identifier == :time_derivative
            linestyle --> :dashdot
        end
        color := idx
        title --> "Electron Particles"
        if show_condition
            label := "$name " * @sprintf("[%.3g s⁻¹]", tot) * label
            if identifier in (:ec, :ic, :lh, :nbi, :pellet)
                fill0 --> true
            end
            if !ismissing(cs1de, :particles_inside) && flux
                update_limits!(plot_extrema, cs1de.particles_inside ./ cs1d.grid.surface, cs1d.grid.rho_tor_norm)
                if plot_extrema.active
                    ylim --> plot_extrema.ylim
                end
                ylabel --> "[s⁻¹/m²]"
                normalization = 1.0 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1de, :particles_inside
            elseif !integrated && !ismissing(cs1de, :particles)
                update_limits!(plot_extrema, cs1de.particles, cs1d.grid.rho_tor_norm)
                if plot_extrema.active
                    ylim --> plot_extrema.ylim
                end
                cs1de, :particles
            elseif integrated && !ismissing(cs1de, :particles_inside)
                update_limits!(plot_extrema, cs1de.particles_inside, cs1d.grid.rho_tor_norm)
                if plot_extrema.active
                    ylim --> plot_extrema.ylim
                end
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

@recipe function plot_source1d(cs1di::IMAS.core_sources__source___profiles_1d___ion, v::Val{:particles}, plot_extrema::PlotExtrema=PlotExtrema();
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    show_source_number=false)

    id = recipe_dispatch(cs1di, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    cs1d = parent(parent(cs1di; error_parent_of_nothing=false); error_parent_of_nothing=true)
    source = parent(parent(cs1d; error_parent_of_nothing=false); error_parent_of_nothing=false)
    name, identifier, idx = source_name_identifier(source, name, show_source_number)

    tot = 0.0
    if !ismissing(cs1di, :particles)
        tot = trapz(cs1d.grid.volume, cs1di.particles)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        if identifier == :time_derivative
            linestyle --> :dashdot
        end
        color := idx
        title --> "$(cs1di.label) Particles"
        if show_condition
            label := "$name $(cs1di.label) " * @sprintf("[%.3g s⁻¹]", tot) * label
            if identifier in (:ec, :ic, :lh, :nbi, :pellet)
                fill0 --> true
            end
            if !ismissing(cs1di, :particles_inside) && flux
                update_limits!(plot_extrema, cs1di.particles_inside ./ cs1d.grid.surface, cs1d.grid.rho_tor_norm)
                if plot_extrema.active
                    ylim --> plot_extrema.ylim
                end
                ylabel --> "[s⁻¹/m²]"
                normalization = 1.0 ./ cs1d.grid.surface
                normalization[1] = NaN
                normalization := normalization
                cs1di, :particles_inside
            elseif !integrated && !ismissing(cs1di, :particles)
                update_limits!(plot_extrema, cs1di.particles, cs1d.grid.rho_tor_norm)
                if plot_extrema.active
                    ylim --> plot_extrema.ylim
                end
                cs1di, :particles
            elseif integrated && !ismissing(cs1di, :particles_inside)
                update_limits!(plot_extrema, cs1di.particles_inside, cs1d.grid.rho_tor_norm)
                if plot_extrema.active
                    ylim --> plot_extrema.ylim
                end
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

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:momentum_tor}, plot_extrema::PlotExtrema=PlotExtrema();
    name="",
    label="",
    integrated=false,
    flux=false,
    show_zeros=false,
    show_source_number=false)

    id = recipe_dispatch(cs1d, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Plot flux"; flux)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    source = parent(parent(cs1d; error_parent_of_nothing=false); error_parent_of_nothing=false)
    name, identifier, idx = source_name_identifier(source, name, show_source_number)

    tot = 0.0
    if !ismissing(cs1d, :momentum_tor)
        tot = trapz(cs1d.grid.volume, cs1d.momentum_tor)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        if identifier == :time_derivative
            linestyle --> :dashdot
        end
        color := idx
        title --> "Momentum"
        if show_condition
            label := "$name " * @sprintf("[%.3g N m]", tot) * label
            if identifier in (:ec, :ic, :lh, :nbi, :pellet)
                fill0 --> true
            end
            if !ismissing(cs1d, :torque_tor_inside) && flux
                ylabel --> "[Nm/m²]"
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

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, v::Val{:j_parallel}, plot_extrema::PlotExtrema=PlotExtrema();
    name="",
    label="",
    integrated=false,
    show_zeros=false,
    show_source_number=false)

    id = recipe_dispatch(cs1d, v)
    assert_type_and_record_argument(id, AbstractString, "Name of the source"; name)
    assert_type_and_record_argument(id, Union{Nothing,AbstractString}, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Plot integrated values"; integrated)
    assert_type_and_record_argument(id, Bool, "Show zeros"; show_zeros)
    assert_type_and_record_argument(id, Bool, "Show source number"; show_source_number)

    source = parent(parent(cs1d; error_parent_of_nothing=false); error_parent_of_nothing=false)
    name, identifier, idx = source_name_identifier(source, name, show_source_number)

    tot = 0.0
    if !ismissing(cs1d, :j_parallel)
        tot = trapz(cs1d.grid.area, cs1d.j_parallel)
    end
    show_condition = show_zeros || identifier in [:total, :time_derivative] || abs(tot) > 0.0
    @series begin
        if identifier == :time_derivative
            linestyle --> :dashdot
        end
        color := idx
        title --> "Parallel Current"
        if show_condition
            label := "$name " * @sprintf("[%.3g MA]", tot / 1E6) * label
            if !integrated && !ismissing(cs1d, :j_parallel)
                if identifier in (:ec, :ic, :lh, :nbi, :pellet)
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
        idx = 1
    end
    if source !== nothing
        identifier = identifier_name(source)
    else
        identifier = :undefined
    end
    return name, identifier, idx
end

@recipe function plot_source1d(cs1d::IMAS.core_sources__source___profiles_1d, plots_extrema::Vector{PlotExtrema}=PlotExtrema[];
    name="",
    label="",
    integrated=false,
    flux=false,
    only=nothing,
    show_zeros=false,
    min_power=0.0,
    only_positive_negative=0,
    show_source_number=false,
    ions=Symbol[])

    id = recipe_dispatch(cs1d)
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
        Nr = round(Int, sqrt(N), RoundUp)
        size --> (1100, Nr * 290)
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)
    end

    for (k, path) in enumerate(paths)
        if length(plots_extrema) < k
            resize!(plots_extrema, k)
            plots_extrema[k] = PlotExtrema()
        end
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
                        goto(cs1d, path[1:end-1]), Val(path[end]), plots_extrema[k]
                    else
                        cs1d, Val(Symbol(join(filter(x -> !(typeof(x) <: Int), path), "__"))), plots_extrema[k]
                    end
                else
                    [NaN], [NaN]
                end
            end
        end
    end
end

@recipe function plot_source(source::IMAS.core_sources__source, plots_extrema::Vector{PlotExtrema}=PlotExtrema[]; time0=global_time(source))
    id = recipe_dispatch(source)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    if !isempty(source.profiles_1d) && source.profiles_1d[1].time <= time0
        @series begin
            name := source.identifier.name
            source.profiles_1d[time0], plots_extrema
        end
    end
end

@recipe function plot_core_sources(cs::IMAS.core_sources{T}; ions=[:my_ions], time0=global_time(cs), aggregate_radiation=false, aggregate_hcd=false) where {T<:Real}
    id = recipe_dispatch(cs)
    assert_type_and_record_argument(id, AbstractVector{Symbol}, "List of ions"; ions)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Bool, "Aggregate radiation sources"; aggregate_radiation)
    assert_type_and_record_argument(id, Bool, "Aggregate heating and current drive sources"; aggregate_hcd)

    dd = top_dd(cs)

    if ions == [:my_ions]
        if dd !== nothing
            ions = list_ions(cs, dd.core_profiles; time0)
        else
            ions = list_ions(cs; time0)
        end
    end

    plots_extrema = PlotExtrema[]
    all_indexes = [source.identifier.index for source in cs.source]
    exclude_indexes = Int[]
    if aggregate_radiation
        append!(exclude_indexes, index_radiation_sources)
    end
    if aggregate_hcd
        append!(exclude_indexes, index_hcd_sources)
    end
    for source in cs.source
        if !retain_source(source, all_indexes, Int[], exclude_indexes)
            continue
        end
        @series begin
            nozeros := true
            time0 := time0
            ions := ions
            source, plots_extrema
        end
    end

    if dd !== nothing
        if aggregate_radiation
            rad_source = IMAS.core_sources__source{T}()
            resize!(rad_source.profiles_1d, 1)
            merge!(rad_source.profiles_1d[1], total_radiation_sources(dd; time0))
            rad_source.identifier.index = 200
            rad_source.identifier.name = "radiation"
            @series begin
                color := :orange
                ions := ions
                rad_source, plots_extrema
            end
        end

        if aggregate_hcd
            for (source_type, color) in ((:ec, :blue), (:ic, :green), (:lh, :magenta), (:nbi, :red), (:pellet, :purple))
                hcd_source = IMAS.core_sources__source{T}()
                idx = name_2_index(hcd_source)[source_type]
                resize!(hcd_source.profiles_1d, 1)
                merge!(hcd_source.profiles_1d[1], total_sources(dd; time0, include_indexes=[idx]))
                hcd_source.identifier.index = idx
                hcd_source.identifier.name = "$source_type"
                @series begin
                    color := color
                    ions := ions
                    hcd_source, plots_extrema
                end
            end
        end

        tot_source = IMAS.core_sources__source{T}()
        resize!(tot_source.profiles_1d, 1)
        merge!(tot_source.profiles_1d[1], total_sources(dd; time0))
        tot_source.identifier.index = 1
        tot_source.identifier.name = "total"
        @series begin
            linewidth := 2
            color := :black
            min_power := 0.0
            ions := ions
            tot_source, plots_extrema
        end
    end

end

#= ============= =#
#  core_profiles  #
#= ============= =#
@recipe function plot_core_profiles(cp::IMAS.core_profiles; time0=global_time(cp))
    id = recipe_dispatch(cp)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    @series begin
        return cp.profiles_1d[time0]
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d; only=nothing)
    id = recipe_dispatch(cpt)
    assert_type_and_record_argument(id, Union{Nothing,Int}, "Plot only this subplot number"; only)

    if only === nothing
        layout := (1, 3)
        size --> (1100, 290)
        margin --> 5 * Measures.mm
    end

    if only === nothing || only == 1
        @series begin
            if only === nothing
                subplot := 1
            end
            cpt, Val(:temperature)
        end
    end

    if only === nothing || only == 2
        @series begin
            if only === nothing
                subplot := 2
            end
            cpt, Val(:density)
        end
    end

    if only === nothing || only == 3
        @series begin
            if only === nothing
                subplot := 3
            end
            cpt, Val(:sonic_rotation)
        end
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:electrons__temperature}; label="", thomson_scattering=false)
    id = recipe_dispatch(cpt, v)
    assert_type_and_record_argument(id, AbstractString, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Overlay Thomson scattering data"; thomson_scattering)

    @series begin
        title --> "Temperatures"
        label := "e" * label
        ylim --> (0, Inf)
        cpt.electrons, :temperature
    end

    if thomson_scattering
        try
            dd = IMAS.top_dd(cpt)
            if !isempty(dd.thomson_scattering)
                @series begin
                    title --> "Temperatures"
                    time0 := cpt.time
                    dd.thomson_scattering, :t_e
                end
            end
        catch
        end
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:ions__temperature}; label="", charge_exchange=false)
    id = recipe_dispatch(cpt, v)
    assert_type_and_record_argument(id, AbstractString, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Overlay charge exchange data"; charge_exchange)

    same_temps = false
    if length(cpt.ion) > 1
        same_temps = !any(x -> x == false, [iion.temperature == cpt.ion[1].temperature for iion in cpt.ion[2:end] if !ismissing(iion, :temperature)])
        if same_temps
            @series begin
                title --> "Temperatures"
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
                    title --> "Temperatures"
                    label := ion.label * label
                    linestyle --> :dash
                    ylim --> (0, Inf)
                    ion, :temperature
                end
            end
        end
    end

    if charge_exchange
        try
            dd = IMAS.top_dd(cpt)
            if !isempty(dd.charge_exchange)
                @series begin
                    title --> "Temperatures"
                    markershape := :diamond
                    time0 := cpt.time
                    dd.charge_exchange, :t_i
                end
            end
        catch
        end
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:temperature}; label="")
    @series begin
        label := label
        cpt, Val(:electrons__temperature)
    end

    @series begin
        label := label
        cpt, Val(:ions__temperature)
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:electrons__density}; label="", what_density=:density, thomson_scattering=false)
    id = recipe_dispatch(cpt, v)
    assert_type_and_record_argument(id, AbstractString, "Label for the plot"; label)
    assert_type_and_record_argument(id, Symbol, "What density to plot: [:density, :density_thermal, :density_fast]"; what_density)
    assert_type_and_record_argument(id, Bool, "Overlay Thomson scattering data"; thomson_scattering)

    @series begin
        title --> "Densities"
        label := "e" * label
        ylim --> (0.0, Inf)
        cpt.electrons, what_density
    end

    if thomson_scattering
        try
            dd = IMAS.top_dd(cpt)
            if !isempty(dd.thomson_scattering)
                @series begin
                    title --> "Densities"
                    time0 := cpt.time
                    dd.thomson_scattering, :n_e
                end
            end
        catch
        end
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:ions__density}; label="", what_density=:density)
    id = recipe_dispatch(cpt, v)
    assert_type_and_record_argument(id, AbstractString, "Label for the plot"; label)
    assert_type_and_record_argument(id, Symbol, "What density to plot: [:density, :density_thermal, :density_fast]"; what_density)

    for ion in cpt.ion
        @series begin
            Z = ion.element[1].z_n
            title --> "Densities"
            if Z == 1.0
                label := ion.label * label
            else
                label := "$(ion.label) × " * @sprintf("%.3g", Z) * label
            end
            linestyle --> :dash
            ylim --> (0.0, Inf)
            normalization --> Z
            ion, what_density
        end
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:density}; label="", greenwald=false, what_density=:density)
    @series begin
        label := label
        what_density := what_density
        cpt, Val(:electrons__density)
    end

    if greenwald
        @series begin
            seriestype := :hline
            primary := false
            style := :dashdotdot
            [IMAS.greenwald_density(IMAS.top_dd(cpt).equilibrium.time_slice[cpt.time])]
        end
    end

    @series begin
        label := label
        what_density := what_density
        cpt, Val(:ions__density)
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:sonic_rotation}; label="")
    id = recipe_dispatch(cpt, v)
    assert_type_and_record_argument(id, AbstractString, "Label for the plot"; label)

    @series begin
        title --> "Sonic rotation"
        label := label
        if !ismissing(cpt, :rotation_frequency_tor_sonic)
            cpt, :rotation_frequency_tor_sonic
        else
            [NaN], [NaN]
        end
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d, v::Val{:ions__toroidal_rotation}; label="", charge_exchange=false)
    id = recipe_dispatch(cpt, v)
    assert_type_and_record_argument(id, AbstractString, "Label for the plot"; label)
    assert_type_and_record_argument(id, Bool, "Overlay charge exchange data"; charge_exchange)

    @series begin
        title --> "Toroidal rotation"
        label := "Ions" * label
        linestyle --> :dash
        ylim --> (0, Inf)
        cpt.ion[1], :rotation_frequency_tor
    end

    if charge_exchange
        try
            dd = IMAS.top_dd(cpt)
            if !isempty(dd.charge_exchange)
                @series begin
                    title --> "Toroidal rotation"
                    markershape := :circle
                    time0 := cpt.time
                    dd.charge_exchange, :ω_tor
                end
            end
        catch
        end
    end
end

#= ============ =#
#  ec_launchers  #
#= ============ =#
function label(beam::IMAS.ec_launchers__beam)
    name = beam.name
    freq = frequency(beam.frequency)
    if ismissing(beam, :mode)
        mode = "?"
    else
        mode = beam.mode == 1 ? "O" : "X"
    end
    if ismissing(beam, :steering_angle_pol)
        angle_pol = "?"
    else
        angle_pol = "$(@sprintf("%.1f", @ddtime(beam.steering_angle_pol) * 180 / pi))"
    end
    if ismissing(beam, :steering_angle_tor)
        angle_tor = "?"
    else
        angle_tor = "$(@sprintf("%.1f", @ddtime(beam.steering_angle_tor) * 180 / pi))"
    end
    return "$name @ $(@sprintf("%.0f", freq/1E9)) GHz $mode ∠ₜ=$(angle_tor)° ∠ₚ=$(angle_pol)°"
end

@recipe function plot_ec_beam(beam::IMAS.ec_launchers__beam; time0=global_time(beam))
    id = recipe_dispatch(beam)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    dd = top_dd(beam)
    if dd !== nothing && !isempty(dd.waves.time) && time0 >= dd.waves.time[1]
        beam_index = index(beam)
        @series begin
            label --> label(beam)
            time0 := time0
            dd.waves.coherent_wave[beam_index]
        end
    end
end

@recipe function plot_ec_beam_frequencies(beam_freqs::AbstractVector{<:IMAS.ec_launchers__beam___frequency{T}}) where {T<:Real}
    freqs_dict = Dict{Float64,ec_launchers__beam___frequency{T}}()
    for beam_freq in beam_freqs
        freqs_dict[frequency(beam_freq)] = beam_freq
    end
    for beam_freq in values(freqs_dict)
        @series begin
            beam_freq
        end
    end
end

@recipe function plot_ec_beam_frequency(beam_freq::IMAS.ec_launchers__beam___frequency; time0::global_time(beam_freq), show_vacuum_Bt=false)
    id = recipe_dispatch(beam_freq)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Bool, "Show EC resonance assuming only vacuum Bt"; show_vacuum_Bt)

    dd = top_dd(beam_freq)
    eqt = dd.equilibrium.time_slice[time0]

    freq = frequency(beam_freq)
    resonance_layer = ech_resonance_layer(eqt, freq)
    @series begin
        label --> "$(resonance_layer.harmonic) @ harmonic $(@sprintf("%.0f", freq/1E9)) GHz"
        resonance_layer.r, resonance_layer.z
    end
    if show_vacuum_Bt
        resonance_layer = ech_resonance_layer(eqt, freq; only_vacuum_Bt=true)
        @series begin
            primary := false
            linestyle := :dash
            resonance_layer.r, resonance_layer.z
        end
    end
end

@recipe function plot_ec(ec::IMAS.ec_launchers{T}) where {T<:Real}
    @series begin
        [beam.power_launched for beam in ec.beam]
    end
end

@recipe function plot_ec(beams::AbstractVector{<:IMAS.ec_launchers__beam{T}}; time0::global_time(beams)) where {T<:Real}
    id = recipe_dispatch(beams)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    dd = top_dd(beams)
    eqt = dd.equilibrium.time_slice[time0]

    @series begin
        cx := true
        color := :gray
        coordinate := :rho_tor_norm
        eqt
    end

    @series begin
        color := :black
        dd.wall
    end

    @series begin
        [beam.frequency for beam in beams]
    end

    for beam in beams
        @series begin
            beam
        end
    end
end

#= === =#
#  nbi  #
#= === =#
function label(unit::IMAS.nbi__unit)
    if !isempty(unit.beamlets_group) && !ismissing(unit.beamlets_group[1], :angle)
        angle_tor = unit.beamlets_group[1].direction * asin(unit.beamlets_group[1].tangency_radius / unit.beamlets_group[1].position.r) * 180 / pi
        angle_pol = unit.beamlets_group[1].angle
        return "$(unit.name) (∠ₜ=$(@sprintf("%.0f", angle_tor))°, ∠ₚ=$(@sprintf("%.0f", angle_pol))°)"
    else
        return unit.name
    end
end

@recipe function plot_nb(nb::IMAS.nbi{T}) where {T<:Real}
    @series begin
        [unit.power_launched for unit in nb.unit if sum(unit.power_launched.data) > 0.0]
    end
end

@recipe function plot_power_lauched(
    pl::Union{
        IMAS.ec_launchers__beam___power_launched{T},
        IMAS.nbi__unit___power_launched{T},
        IMAS.ic_antennas__antenna___power_launched{T},
        IMAS.lh_antennas__antenna___power_launched{T}
    };
    smooth_tau=0.0
) where {T<:Real}
    id = recipe_dispatch(pl)
    assert_type_and_record_argument(id, Float64, "Smooth instantaneous power assuming τ decorrelation time [s]"; smooth_tau)
    hcd = parent(pl)
    data1 = pl.data
    if smooth_tau > 0.0
        data1 = smooth_beam_power(pl.time, pl.data, smooth_tau)
    end
    @series begin
        label --> label(hcd)
        ylabel := "Power Launched [MW]"
        xlabel := "Time [s]"
        legend_position --> :topleft
        pl.time, data1 ./ 1E6
    end
end

@recipe function plot_hcd_power_lauched(
    pls::AbstractVector{
        <:Union{
            IMAS.ec_launchers__beam___power_launched{T},
            IMAS.nbi__unit___power_launched{T},
            IMAS.ic_antennas__antenna___power_launched{T},
            IMAS.lh_antennas__antenna___power_launched{T}
        }
    };
    show_total=true,
    smooth_tau=0.0
) where {T<:Real}
    id = recipe_dispatch(pls)
    assert_type_and_record_argument(id, Bool, "Show total power launched"; show_total)
    assert_type_and_record_argument(id, Float64, "Smooth instantaneous power"; smooth_tau)

    background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)

    if show_total
        # time basis for the total
        time_range = T[]
        for pl in pls
            append!(time_range, pl.time)
            unique!(time_range)
        end

        # accumulate total
        total = time_range .* 0.0
        for pl in pls
            if time_range == pl.time
                total += pl.data
            else
                total += interp1d(pl.time, pl.data).(time_range)
            end
        end
        if smooth_tau > 0.0
            total = smooth_beam_power(time_range, total, smooth_tau)
        end

        # plot total
        @series begin
            color := :black
            lw := 2
            label := "Total"
            time_range, total / 1E6
        end
    end

    for pl in pls
        @series begin
            smooth_tau := smooth_tau
            pl
        end
    end
end

#= ===== =#
#  waves  #
#= ===== =#
@recipe function plot_waves_profiles_1d(wvc1d::IMAS.waves__coherent_wave___profiles_1d)
    layout := RecipesBase.@layout [1, 1]
    if !ismissing(parent(parent(wvc1d)).identifier, :antenna_name)
        name = parent(parent(wvc1d)).identifier.antenna_name
    else
        name = "Antenna $(index(parent(parent(wvc1d))))"
    end

    @series begin
        subplot := 1
        title := "Power density"
        label := name
        wvc1d, :power_density
    end
    @series begin
        subplot := 2
        title := "Parallel current density"
        label := name
        wvc1d, :current_parallel_density
    end
end

@recipe function plot_waves_profiles_1d(wvc1ds::AbstractVector{<:IMAS.waves__coherent_wave___profiles_1d})
    for wvc1d in wvc1ds
        @series begin
            wvc1d
        end
    end
end

@recipe function plot_beam(beam::IMAS.waves__coherent_wave___beam_tracing___beam; top=false, min_power=0.0)
    id = recipe_dispatch(beam)
    assert_type_and_record_argument(id, Bool, "Top view"; top)
    assert_type_and_record_argument(id, Float64, "Minimum power above which to show the beam"; min_power)

    active = min_power < 0.0 || ismissing(beam, :power_initial) || beam.power_initial > min_power
    name = parent(parent(parent(parent(beam)))).identifier.antenna_name
    color = hash_to_color(name; seed=4)
    if top
        # top view
        if active
            color --> color
            @series begin
                seriestype --> :scatter
                markerstrokewidth := 0.0
                [beam.position.r[1] * cos(beam.position.phi[1])], [beam.position.r[1] * sin(beam.position.phi[1])]
            end
            @series begin
                primary := false
                beam.position.r .* cos.(beam.position.phi), beam.position.r .* sin.(beam.position.phi)
            end
            @series begin
                primary := false
                linewidth := 0.5
                color := :black
                beam.position.r .* cos.(beam.position.phi), beam.position.r .* sin.(beam.position.phi)
            end
        else
            @series begin
                color := :gray
                seriestype --> :scatter
                markerstrokewidth := 0.0
                [beam.position.r[1] * cos(beam.position.phi[1])], [beam.position.r[1] * sin(beam.position.phi[1])]
            end
        end
    else
        # cross-sectional view
        if active
            color --> color
            @series begin
                seriestype --> :scatter
                markerstrokewidth := 0.0
                [beam.position.r[1]], [beam.position.z[1]]
            end
            @series begin
                primary := false
                beam.position.r, beam.position.z
            end
            @series begin
                primary := false
                linewidth := 0.5
                color := :black
                beam.position.r, beam.position.z
            end
        else
            @series begin
                color := :gray
                seriestype --> :scatter
                markerstrokewidth := 0.0
                [beam.position.r[1]], [beam.position.z[1]]
            end
        end
    end
end

@recipe function beam_tracing(bt::IMAS.waves__coherent_wave___beam_tracing; max_beam=length(bt.beam))
    id = recipe_dispatch(bt)
    assert_type_and_record_argument(id, Int, "Maximum number of beams to show"; max_beam)

    aspect_ratio := :equal
    legend_position --> :outerbottomright
    for (ibeam, beam) in enumerate(bt.beam[1:max_beam])
        @series begin
            if ibeam == 1
                primary := true
                label --> parent(parent(bt)).identifier.antenna_name
                lw := 3.0
            else
                primary := false
                label := ""
                lw := 1.0
            end
            beam
        end
    end
end

@recipe function plot_wave(bts::IDSvector{<:IMAS.waves__coherent_wave___beam_tracing}; time0=global_time(bts))
    id = recipe_dispatch(bts)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    if !isempty(bts) && bts[1].time >= time0
        @series begin
            bts[time0]
        end
    end
end

@recipe function plot_wave(bts::AbstractVector{<:IMAS.waves__coherent_wave___beam_tracing})
    for bt in bts
        @series begin
            bt
        end
    end
end

@recipe function plot_wave(wv::IMAS.waves; max_beam=length(wv.coherent_wave), time0=global_time(wv))
    id = recipe_dispatch(wv)
    assert_type_and_record_argument(id, Int, "Maximum number of beams to show"; max_beam)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)

    dd = top_dd(wv)
    @series begin
        cx := true
        dd.equilibrium
    end
    @series begin
        alpha --> 0.75
        [beam.beam_tracing[time0] for beam in wv.coherent_wave[1:max_beam]]
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
        if dd !== nothing
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
    base_linewidth = get(plotattributes, :linewidth, 1.0)

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
    id = recipe_dispatch(nwl)
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

        l = arc_length(neutronics.first_wall.r, neutronics.first_wall.z; include_zero=false)

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
@recipe function plot_wall(wall::IMAS.wall; top=false)
    id = recipe_dispatch(wall)
    assert_type_and_record_argument(id, Bool, "Top view"; top)

    aspect_ratio := :equal
    if top
        fw = first_wall(wall)
        z0 = (maximum(fw.z) + minimum(fw.z)) / 2
        _, crossings = intersection(fw.r, fw.z, [0.0, 1000.0], [z0, z0])
        color := :black
        label := ""
        r1 = min(crossings[1][1], crossings[2][1])
        r2 = max(crossings[1][1], crossings[2][1])
        @series begin
            _circle(r1)
        end
        @series begin
            xlim --> (-r2 * 1.1, r2 * 1.1)
            ylim --> (-r2 * 1.1, r2 * 1.1)
            primary := false
            _circle(r2)
        end
    else
        @series begin
            wall.description_2d
        end
    end
end

@recipe function plot_wd2d_list(wd2d_list::AbstractVector{<:IMAS.wall__description_2d})
    for wd2d in wd2d_list
        @series begin
            return wd2d
        end
    end
end

@recipe function plot_wd2d(wd2d::IMAS.wall__description_2d; show_limiter=true, show_vessel=true)
    id = recipe_dispatch(wd2d)
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

@recipe function plot_wd2dvu_list(wd2dvu_list::AbstractVector{<:IMAS.wall__description_2d___vessel__unit})
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

@recipe function plot_wd2dvue(elements::AbstractVector{<:IMAS.wall__description_2d___vessel__unit___element})
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

@recipe function plot_wd2dlu_list(wd2dlu_list::AbstractVector{<:IMAS.wall__description_2d___limiter__unit})
    for wd2dlu in wd2dlu_list
        @series begin
            return wd2dlu
        end
    end
end

@recipe function plot_wd2dlu(wd2dlu::IMAS.wall__description_2d___limiter__unit)
    if ismissing(wd2dlu, :closed)
        poly = closed_polygon(wd2dlu.outline.r, wd2dlu.outline.z)
    else
        poly = closed_polygon(wd2dlu.outline.r, wd2dlu.outline.z, Bool(wd2dlu.closed))
    end
    @series begin
        label --> getproperty(wd2dlu, :name, "")
        poly.r, poly.z
    end
end

#= ============== =#
#  pulse_schedule  #
#= ============== =#
const UnionPulseScheduleSubIDS = Union{IMAS.pulse_schedule,(tp for tp in fieldtypes(IMAS.pulse_schedule) if tp <: IDS)...}
@recipe function plot_ps(ps::UnionPulseScheduleSubIDS; time0=global_time(ps), simulation_start=nothing)
    id = recipe_dispatch(ps)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Union{Nothing,Float64}, "Simulation start time"; simulation_start)
    plots = []
    for (loc, (ids, field)) in filled_ids_fields(ps; eval_expr=true)
        if field != :reference
            continue
        end
        path = f2p(ids)
        if "boundary_outline" in path
            continue
        end

        time_value = getproperty(coordinates(ids, :reference)[1])
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
        h = ids
        while h !== nothing
            if hasfield(typeof(h), :name) && hasdata(h, :name)
                plt[:label] = h.name
                break
            end
            h = parent(h)
        end
        push!(plots, plt)
    end

    size --> (1200, 1200)

    if typeof(ps) <: IMAS.pulse_schedule && !isempty(ps.position_control) || typeof(ps) <: IMAS.pulse_schedule__position_control
        plot_offset = 1
        layout := RecipesBase.@layout [length(plots) + plot_offset]
        @series begin
            subplot := 1
            label := "$(time0) [s]"
            aspect_ratio := :equal
            time0 := time0
            if typeof(ps) <: IMAS.pulse_schedule
                ps.position_control, nothing
            else
                ps, nothing
            end
        end
        eqt = try
            top_dd(pc).equilibrium.time_slice[time0]
        catch
            nothing
        end
        if eqt !== nothing
            @series begin
                legend := false
                cx := true
                alpha := 0.2
                color := :gray
                eqt
            end
        end

        wall = try
            top_dd(pc).wall
        catch
            nothing
        end
        if wall !== nothing
            @series begin
                show_limiter := true
                show_vessel := false
                legend := false
                wall
            end
        end
    else
        plot_offset = 0
        layout := RecipesBase.@layout [length(plots)]
    end

    # plotting at infinity does not show
    tmax = -Inf
    tmin = Inf
    for plt in plots
        xx = filter(x -> !isinf(x), plt[:x])
        tmax = max(tmax, maximum(xx))
        tmin = min(tmax, minimum(xx))
    end

    for (k, plt) in enumerate(plots)
        x = deepcopy(plt[:x])
        x[x.==Inf] .= tmax
        x[x.==-Inf] .= tmin
        y = plt[:y]
        n = min(length(x), length(y))
        @series begin
            subplot := k + plot_offset
            label := ""
            x[1:n], y[1:n]
        end
        if simulation_start !== nothing
            @series begin
                subplot := k + plot_offset
                primary := false
                seriestype := :vline
                linestyle := :dash
                [simulation_start]
            end
        end
        @series begin
            subplot := k + plot_offset
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

@recipe function plot_pc_time(pc::IMAS.pulse_schedule__position_control, ::Nothing; time0=global_time(pc))
    id = recipe_dispatch(pc)
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
        [x[1] for x in Xs], [x[2] for x in Xs]
    end
    @series begin
        seriestype := :scatter
        marker --> :circle
        markerstrokewidth --> 0
        primary := false
        Ss = strike_points(pc.strike_point; time0)
        [x[1] for x in Ss], [x[2] for x in Ss]
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
@recipe function plot_controllers(controller_outputs::IMAS.controllers__linear_controller___outputs, k::Int)
    controller = IMAS.parent(controller_outputs)

    data = controller.inputs.data[1, :]
    integral = cumtrapz(controller.inputs.time, data)
    derivative = [0.0, diff(data) ./ diff(controller.inputs.time)]

    @series begin
        label := controller.output_names[k]
        color := :black
        lw := 2
        controller.outputs.time, controller.outputs.data[1, :]
    end

    @series begin
        label := "P=$(join(map(x->@sprintf("%.3e", x), @ddtime(controller.pid.p.data)),", "))"
        controller.outputs.time, data .* controller.pid.p.data[k]
    end

    @series begin
        label := "I=$(join(map(x->@sprintf("%.3e", x), @ddtime(controller.pid.i.data)),", "))"
        controller.outputs.time, integral .* controller.pid.i.data[k]
    end

    if any(controller.pid.d.data[k] .!= 0.0)
        @series begin
            label := "D=$(join(map(x->@sprintf("%.3e", x), @ddtime(controller.pid.d.data)),", "))"
            controller.outputs.time, derivative .* controller.pid.d.data[k]
        end
    end

    if any(controller.pid.i.data[k] .!= 0.0)
        @series begin
            label := "PI"
            controller.outputs.time, data .* controller.pid.p.data[k] .+ integral .* controller.pid.i.data[k]
        end
    end

    if any(controller.pid.d.data[k] .!= 0.0) && any(controller.pid.i.data[k] .!= 0.0)
        @series begin
            label := "PID"
            controller.outputs.time, data .* controller.pid.p.data[k] .+ integral .* controller.pid.i.data[k] .+ derivative .* controller.pid.d.data[k]
        end
    end

end

@recipe function plot_controllers(controller_inputs::IMAS.controllers__linear_controller___inputs, k::Int)
    controller = IMAS.parent(controller_inputs)
    @series begin
        label := controller.input_names[k]
        controller_inputs.time, controller_inputs.data[k, :]
    end
end

@recipe function plot_controllers(controller::IMAS.controllers__linear_controller)
    n = size(controller.inputs.data)[1] + size(controller.outputs.data)[1]
    layout := (n, 1)
    size --> (800, n * 300)
    margin --> 5 * Measures.mm
    for k in 1:size(controller.inputs.data)[1]
        @series begin
            subplot := 1
            controller.inputs, k
        end
    end
    for k in 1:size(controller.outputs.data)[1]
        @series begin
            subplot := 2
            controller.outputs, k
        end
    end
end

#= ======= =#
#  getdata  #
#= ======= =#
@recipe function plot_getdata(ids::IDS, what::Val; time0=global_time(ids), time_averaging=0.05, normalization=1.0)
    id = recipe_dispatch(ids)
    assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    assert_type_and_record_argument(id, Float64, "Time averaging window"; time_averaging)
    assert_type_and_record_argument(id, Float64, "Normalization factor"; normalization)

    time, rho, data, weights, units = getdata(what, top_dd(ids), time0, time_averaging)
    if normalization != 1.0
        data = data .* normalization
        units = "$normalization $units"
    end
    units = nice_units(units)

    if eltype(data) <: Measurements.Measurement
        data_σ = [d.err for d in data]
        data = [d.val for d in data]
        if all(data_σ .== 0.0)
            data_σ = []
        end
    else
        data_σ = []
    end

    if isempty(data_σ)
        @series begin
            ylabel := "[$units]"
            label --> string(what)
            color := :transparent
            seriestype := :scatter
            seriesalpha := weights
            rho, data
        end
    else
        k = 1
        for w in unique(weights)
            index = findall(==(w), weights)
            primary := (k == 1)
            k += 1
            @series begin
                ylabel := units
                label --> string(what)
                color := :transparent
                seriestype := :scatter
                alpha := w
                yerror := isempty(data_σ) ? nothing : data_σ[index]
                rho[index], data[index]
            end
        end
    end
end

#= ================== =#
#  thomson_scattering  #
#= ================== =#
@recipe function plot_thomson(ts::IMAS.thomson_scattering, what::Symbol; time0=global_time(ts), time_averaging=0.05, normalization=1.0)
    @assert what in (:n_e, :t_e)

    @series begin
        label := "Thomson scattering"
        seriestype := :scatter
        markersize := 3
        color := :black
        [NaN], [NaN]
    end

    @series begin
        time0 := time0
        label := ""
        time_averaging := time_averaging
        normalization := normalization
        ts, Val(what)
    end
end

#= =============== =#
#  charge_exchange  #
#= =============== =#
@recipe function plot_charge_exchange(cer::IMAS.charge_exchange, what::Symbol; time0=global_time(cer), time_averaging=0.05, normalization=1.0)
    @assert what in (:t_i, :n_i_over_n_e, :zeff, :n_imp, :ω_tor)

    @series begin
        label := "Charge exchange"
        seriestype := :scatter
        markersize := 3
        color := :black
        [NaN], [NaN]
    end

    @series begin
        time0 := time0
        label := ""
        time_averaging := time_averaging
        normalization := normalization
        cer, Val(what)
    end
end

#= ============== =#
#  interferometer  #
#= ============== =#
@recipe function plot_interferometer_line_of_sight(icls::IMAS.interferometer__channel___line_of_sight; top=false)
    id = recipe_dispatch(icls)
    assert_type_and_record_argument(id, Bool, "Top view"; top)

    label --> parent(icls).name

    if isempty(icls.third_point)
        if top
            a1 = icls.first_point.r * cos(icls.first_point.phi)
            a2 = icls.second_point.r * cos(icls.second_point.phi)
            b1 = icls.first_point.r * sin(icls.first_point.phi)
            b2 = icls.second_point.r * sin(icls.second_point.phi)
        else
            a1 = icls.first_point.r
            a2 = icls.second_point.r
            b1 = icls.first_point.z
            b2 = icls.second_point.z
        end

        if a1 == a2 && b1 == b2
            @series begin
                seriestype := :scatter
                [a1, a2], [b1, b2]
            end
        else
            @series begin
                [a1, a2], [b1, b2]
            end
        end

    else
        if top
            a1 = icls.first_point.r * cos(icls.first_point.phi)
            a2 = icls.second_point.r * cos(icls.second_point.phi)
            a3 = icls.third_point.r * cos(icls.third_point.phi)
            b1 = icls.first_point.r * sin(icls.first_point.phi)
            b2 = icls.second_point.r * sin(icls.second_point.phi)
            b3 = icls.third_point.r * sin(icls.third_point.phi)
        else
            a1 = icls.first_point.r
            a2 = icls.second_point.r
            a3 = icls.third_point.r
            b1 = icls.first_point.z
            b2 = icls.second_point.z
            b3 = icls.third_point.z
        end

        if a1 == a2 == a3 && b1 == b2 == b3
            @series begin
                seriestype := :scatter
                [a1], [b1]
            end
        else
            @series begin
                [a1, a2, a3], [b1, b2, b3]
            end
            @series begin
                primary := false
                seriestype := :scatter
                marker := :square
                markerstrokewidth := 0.0
                [a2], [b2]
            end
        end
    end
end

@recipe function plot_interferometer_line_of_sights(iclss::AbstractVector{<:IMAS.interferometer__channel___line_of_sight})
    for icls in iclss
        @series begin
            icls
        end
    end
end

@recipe function plot_interferometer(ifrmtr::IMAS.interferometer)
    @series begin
        [ch.line_of_sight for ch in ifrmtr.channel]
    end
end

#= ======= =#
#  summary  #
#= ======= =#
@recipe function plot(summary::IMAS.summary)
    valid_leaves = [leaf for leaf in IMASdd.AbstractTrees.Leaves(summary) if typeof(leaf.value) <: Vector && leaf.field != :time]
    layout := length(valid_leaves)
    N = round(Int, sqrt(length(valid_leaves)), RoundUp)
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

@recipe function plot(hcd::IMAS.summary__heating_current_drive)
    @series begin
        color := :black
        linewidth := 2
        label := "total"
        normalization := 1E-6
        ylabel --> "[MW]"
        getproperty(hcd, :power_launched_total), :value
    end
    for (field, label) in ((:power_launched_ec, "ec"), (:power_launched_ic, "ic"), (:power_launched_lh, "lh"), (:power_launched_nbi, "nbi"))
        if !ismissing(getproperty(hcd, field), :value)
            @series begin
                label := label
                normalization := 1E-6
                ylabel --> "[MW]"
                getproperty(hcd, field), :value
            end
        end
    end
end

#= ================ =#
#  generic plotting  #
#= ================ =#
@recipe function plot_multiple_fields(ids::Union{IDS,IDSvector}, target_fields::Union{AbstractArray{Symbol},Regex})
    @series begin
        # calls "plot_IFF_list" recipe
        IFF_list = findall(ids, target_fields)
    end
end

@recipe function plot_IFF_list(IFF_list::AbstractArray{IDS_Field_Finder}; nrows=:auto, ncols=:auto, each_size=(500, 400))
    id = recipe_dispatch(IFF_list)
    assert_type_and_record_argument(id, Tuple{Integer,Integer}, "Size of each subplot. (Default=(500, 400))"; each_size)
    assert_type_and_record_argument(id, Union{Integer,Symbol}, "Number of rows for subplots' layout (Default = :auto)"; nrows)
    assert_type_and_record_argument(id, Union{Integer,Symbol}, "Number of columns for subplots' layout (Default = :auto)"; ncols)

    # Keep only non-empty Vector or Matrix type data
    IFF_list = filter(IFF -> IFF.field_type <: Union{AbstractVector{<:Real},AbstractMatrix{<:Real}}, IFF_list)
    IFF_list = filter(IFF -> length(IFF.value) > 0, IFF_list)

    # At least one valid filed name is required to proceed
    @assert length(IFF_list) > 0 "All field names are complex type, or invalid, or missing"

    my_layout = get(plotattributes, :layout) do
        # Fallback: Compute my_layout using nrows or ncols if provided
        layout_spec = length(IFF_list)

        if nrows !== :auto && ncols == :auto
            layout_spec = (layout_spec, (nrows, :))
        elseif ncols !== :auto && nrows == :auto
            layout_spec = (layout_spec, (:, ncols))
        elseif nrows !== :auto || ncols !== :auto
            layout_spec = (layout_spec, (nrows, ncols))
        end
        return Plots.layout_args(layout_spec...)
    end

    # Define layout and compute nrows & ncols
    if my_layout isa Tuple{Plots.GridLayout,Int}
        if nrows == :auto && ncols == :auto
            layout --> my_layout[2]
        else
            layout --> my_layout[1]
        end
        nrows = size(my_layout[1].grid)[1]
        ncols = size(my_layout[1].grid)[2]
    elseif my_layout isa Plots.GridLayout
        layout --> my_layout
        nrows = size(my_layout.grid)[1]
        ncols = size(my_layout.grid)[2]
    elseif my_layout isa Int
        layout --> my_layout
        tmp_grid = Plots.layout_args(plotattributes, my_layout)[1].grid
        nrows = size(tmp_grid)[1]
        ncols = size(tmp_grid)[2]
    end

    size := (ncols * each_size[1], nrows * each_size[2])

    legend_position --> :best
    left_margin --> [10 * Measures.mm 10 * Measures.mm]
    bottom_margin --> 10 * Measures.mm

    # Create subplots for the current group
    for IFF in IFF_list
        @series begin
            IFF # (calls "plot_IFF" recipe)
        end
    end
end

const default_abbreviations = Dict(
    "dd." => "",
    "equilibrium" => "eq",
    "description" => "desc",
    "time_slice" => "time",
    "profiles" => "prof",
    "source" => "src",
    "parallel" => "para"
)

# Function to shorten names directly in the original text
function shorten_ids_name(full_name::String, abbreviations::Dict=default_abbreviations)
    # Apply each abbreviation to the full name
    for (key, value) in abbreviations
        value = value isa Function ? value() : value
        full_name = replace(full_name, key => value)
    end
    return full_name
end

@recipe function plot_IFF(IFF::IDS_Field_Finder; abbreviations=default_abbreviations, seriestype_1d=:path, seriestype_2d=:contourf, nicer_title=true)
    if Plots.backend_name() == :unicodeplots && seriestype_2d == :contourf
        # unicodeplots cannot render :contourf, use :contour instead
        seriestype_2d = :contour
    end

    id = recipe_dispatch(IFF)
    assert_type_and_record_argument(id, Dict, "Abbreviations to shorten titles of subplots"; abbreviations)
    assert_type_and_record_argument(id, Symbol, "Seriestype for 1D data [:path (default), :scatter, :bar ...]"; seriestype_1d)
    assert_type_and_record_argument(id, Symbol, "Seriestype for 2D data [:contourf (default), :contour, :surface, :heatmap]"; seriestype_2d)
    assert_type_and_record_argument(id, Bool, "Flag to use nicer title"; nicer_title)

    field_name = shorten_ids_name(IFF.field_path, abbreviations)

    if nicer_title
        field_name = shorten_ids_name(field_name, nice_field_symbols)
        field_name = nice_field(field_name)
    end

    if IFF.field_type <: AbstractVector{<:Real} && length(IFF.value) > 0
        if length(IFF.value) == 1
            seriestype --> :scatter
        else
            seriestype --> seriestype_1d
        end

        title --> field_name
        IFF.parent_ids, IFF.field, Val(:plt_1d) # (calls "plot_field_1d" recipe)

    elseif IFF.field_type <: AbstractMatrix{<:Real} && length(IFF.value) > 0
        seriestype --> seriestype_2d

        if nicer_title
            title --> field_name * " " * nice_units(units(IFF.parent_ids, IFF.field))
        else
            title --> field_name * " [" * units(IFF.parent_ids, IFF.field) * "]"
        end

        IFF.parent_ids, IFF.field, Val(:plt_2d) # calls "plot_field_2d" recipe
    end
end

@recipe function plot_field(ids::IDS, field::Symbol)
    @assert hasfield(typeof(ids), field) "$(location(ids)) does not have field `$field`. Did you mean: $(keys(ids))"

    fType = fieldtype_typeof(ids, field)

    if fType <: Vector{<:Real}
        ids, field, Val(:plt_1d) # calls plot_field_1d recipe
    elseif fType <: Matrix{<:Real}
        ids, field, Val(:plt_2d) # calls plot_field_2d recipe
    end
end

@recipe function plot_field_1d(ids::IDS, field::Symbol, ::Val{:plt_1d}; normalization=1.0, coordinate=nothing, weighted=nothing, fill0=false)
    @assert hasfield(typeof(ids), field) "$(location(ids)) does not have field `$field`. Did you mean: $(keys(ids))"

    id = recipe_dispatch(ids, field)
    assert_type_and_record_argument(id, Union{Real,AbstractVector{<:Real}}, "Normalization factor"; normalization)
    assert_type_and_record_argument(id, Union{Nothing,Symbol}, "Coordinate for x-axis"; coordinate)
    assert_type_and_record_argument(id, Union{Nothing,Symbol}, "Weighting field"; weighted)
    assert_type_and_record_argument(id, Bool, "Fill area under curve"; fill0)

    coords = coordinates(ids, field; override_coord_leaves=[coordinate])
    coordinate_name = string(coords[1].field)
    coordinate_value = getproperty(coords[1]; return_missing_time=true)

    @assert coordinate_value !== missing "Missing coordinate for plotting $(coords[1])"

    # If the field is the reference coordinate of the given IDS,
    # set the coordinate_value as its index
    if coordinate_name == "1...N"
        coordinate_value = 1:length(getproperty(ids, field))
    end

    xvalue = coordinate_value
    yvalue = getproperty(ids, field)

    if hasdata(ids, Symbol("$(field)_σ"))
        yvalue = Measurements.measurement.(yvalue, getproperty(ids, Symbol("$(field)_σ")))
    end

    yvalue = yvalue .* normalization

    # figure out a good label
    lbl = join((name(ids), " ", field))
    h = ids
    while h !== nothing
        if hasfield(typeof(h), :name) && hasdata(h, :name)
            lbl = h.name
            break
        end
        h = parent(h)
    end

    @series begin
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)

        if endswith(coordinate_name, "_norm")
            xlim --> (0.0, 1.0)
        end

        # multiply y by things like `:area` or `:volume`
        if weighted !== nothing
            weight = coordinates(ids, field; override_coord_leaves=[weighted])[1]
            yvalue .*= getproperty(weight)
            ylabel = nice_units(units(ids, field) * "*" * units(weight.field))
            label = nice_field("$lbl*$weighted")
        else
            ylabel = nice_units(units(ids, field))
            label = nice_field(lbl)
        end

        if coordinate_name == "1...N"
            xlabel --> "index"
        else
            xlabel --> nice_field(coordinate_name) * nice_units(units(coords[1].ids, coords[1].field))
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

@recipe function plot_field_1d_manyDDs(ids::IDS, field::Symbol, DDs::AbstractVector{<:IMAS.dd})
    @series begin
        [IMAS.goto(dd, IMAS.location(ids)) for dd in DDs], field
    end
end

@recipe function plot_field_1d_manyIDSs(IDSs::Vector{<:IDS}, field::Symbol; alpha_of_individual_lines=0.1, normalization=1.0, coordinate=nothing)
    id = recipe_dispatch(IDSs, field)
    assert_type_and_record_argument(id, Real, "Alpha  individual lines"; alpha_of_individual_lines)
    assert_type_and_record_argument(id, Union{Real,AbstractVector{<:Real}}, "Normalization factor"; normalization)
    assert_type_and_record_argument(id, Union{Nothing,Symbol}, "Coordinate for x-axis"; coordinate)

    ids = IDSs[1]
    @assert hasfield(typeof(ids), field) "$(location(ids)) does not have field `$field`. Did you mean: $(keys(ids))"
    coords = coordinates(ids, field; override_coord_leaves=[coordinate])

    coords = hcat((getproperty(coordinates(ids, field; override_coord_leaves=[coordinate])[1]) for ids in IDSs)...)

    values = hcat((getproperty(ids, field) for ids in IDSs)...) .* normalization
    mv = Statistics.mean(values; dims=2)[:, 1]
    σv = Statistics.std(values; dims=2)[:, 1]
    mc = Statistics.mean(coords; dims=2)[:, 1]
    σc = Statistics.std(coords; dims=2)[:, 1]
    @series begin
        lw := 2
        ribbon := σv
        mc, mv
    end
    if alpha_of_individual_lines > 0.0
        @series begin
            primary := false
            alpha := alpha_of_individual_lines
            coords, values
        end
    end
end

@recipe function plot_field_2d(ids::IDS, field::Symbol, ::Val{:plt_2d}; normalization=1.0, seriestype=:contourf)
    @assert hasfield(typeof(ids), field) "$(location(ids)) does not have field `$field`. Did you mean: $(keys(ids))"

    id = recipe_dispatch(ids, field)
    assert_type_and_record_argument(id, Union{Real,AbstractVector{<:Real}}, "Normalization factor"; normalization)
    assert_type_and_record_argument(id, Symbol, "Seriestype for 2D data [:contourf (default), :contour, :surface, :heatmap]"; seriestype)

    zvalue = getproperty(ids, field) .* normalization

    @series begin
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)

        coord = coordinates(ids, field)
        dim1 = getproperty(coord[1])
        dim2 = getproperty(coord[2])
        if dim1 === nothing || dim1 === missing || dim2 === nothing || dim2 === missing
            # Plot zvalue as "matrix"
            xlabel --> "column"
            ylabel --> "row"

            xlim --> (1, size(zvalue, 2))
            ylim --> (1, size(zvalue, 1))

            yflip --> true # To make 'row' counting starts from the top
            zvalue
        else
            # Plot zvalue with the coordinates
            xvalue = dim1
            yvalue = dim2
            xlabel --> coord[1].field
            ylabel --> coord[2].field

            xlim --> (minimum(xvalue), maximum(xvalue))
            ylim --> (minimum(yvalue), maximum(yvalue))
            aspect_ratio --> :equal
            xvalue, yvalue, zvalue' # (calls Plots' default recipe for a given seriestype)
        end
    end
end

"""
    ylim(extrema::Dict{Int,Float64})

Set y-limits for n-th subplot in a layout

Negative `n`'s are interpreted as minimum value and positive `n`'s as maximum value

For example:

    ylim(Dict{Int,Float64}(
        3=>3.0, 4=>1E20,
        -6=>-.5, 6=>1.5, -7=>-.5, 7=>1.5, -8=>-2E20, 8=>2.5E20,
        -10=>0.0, 10=>0.101, -11=>0.0, 11=>0.101, -12=>-1.2E19, 12=>1.2E19))
"""
function ylim(extrema::Dict{Int,Float64})
    extrema = deepcopy(extrema)
    for n in collect(keys(extrema))
        if n > 0 && -n ∉ keys(extrema)
            extrema[-n] = -Inf
        elseif n < 0 && -n ∉ keys(extrema)
            extrema[-n] = Inf
        end
    end
    for n in unique!(abs.(collect(keys(extrema))))
        plot!(; ylim=(extrema[-n], extrema[n]), subplot=n)
    end
    return plot!()
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
const nice_field_symbols = Dict()
nice_field_symbols["rho_tor_norm"] = () -> latex_support() ? L"\rho" : "ρ"
nice_field_symbols["psi"] = () -> latex_support() ? L"\psi" : "ψ"
nice_field_symbols["psi_norm"] = () -> latex_support() ? L"\psi_\mathrm{N}" : "ψₙ"
nice_field_symbols["rotation_frequency_tor_sonic"] = "Sonic rotation"
nice_field_symbols["rotation_frequency_tor"] = "Toroidal rotation"
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
        units = replace(units, " " => s"~")
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
