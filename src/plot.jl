import PlotUtils
using RecipesBase
using LaTeXStrings
import Measures

# ========= #
# pf_active #
# ========= #
"""
    plot_pf_active_cx(pfa::pf_active)

Plots pf active cross-section

NOTE: Current plots are for the total current flowing in the coil (ie. it is multiplied by turns_with_sign)
"""
@recipe function plot_pf_active_cx(pfa::pf_active, what::Symbol=:cx; time0=global_time(pfa), cname=:vik)
    @assert typeof(time0) <: Float64
    @assert typeof(cname) <: Symbol

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
                    color --> :black
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

"""
    plot_coil(coil::pf_active__coil{T}; coil_names=false) where {T<:Real}

Plots cross-section of individual coils
"""
@recipe function plot_coil(coil::pf_active__coil{T}; coil_names=false) where {T<:Real}
    @assert typeof(coil_names) <: Bool

    r = T[]
    z = T[]
    for (k, element) in enumerate(coil.element)
        oute = outline(element)
        @series begin
            primary := k == 1
            seriestype --> :shape
            linewidth --> 0.25
            colorbar --> :right
            label --> ""
            oute.r, oute.z
        end
        append!(r, oute.r)
        append!(z, oute.z)
    end

    if coil_names
        r_avg = sum(r) / length(r)
        z_avg = sum(z) / length(z)
        @series begin
            label := ""
            series_annotations := [(coil.name, :center, :middle, :red, 6)]
            [r_avg], [z_avg]
        end
    end
end

@recipe function plot_pf_active_rail(rail::IMAS.build__pf_active__rail)
    if !ismissing(rail.outline, :r)
        @series begin
            color --> :gray
            linestyle --> :dash
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

"""
    plot_loop(loop::pf_passive__loop{T}; loop_names=false) where {T<:Real}

Plots cross-section of individual loops
"""
@recipe function plot_loop(loop::pf_passive__loop{T}; loop_names=false) where {T<:Real}
    @assert typeof(loop_names) <: Bool

    r = T[]
    z = T[]
    for (k, element) in enumerate(loop.element)
        oute = outline(element)
        @series begin
            primary := k == 1
            seriestype --> :shape
            aspect_ratio := :equal
            linewidth --> 0.25
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
@recipe function plot_costing(cst::IMAS.IDSvector{<:IMAS.costing__cost_direct_capital__system})
    costing = top_ids(cst)
    names = ["$(sys_cst.name)" for sys_cst in cst]
    costs = [sys_cst.cost for sys_cst in cst]
    perc = ["$(round(sys_cst.cost/sum(costs)*100))%" for sys_cst in cst]

    name_series = [["$(sub_cst.name)" for sub_cst in reverse(sys_cst.subsystem)] for sys_cst in cst]
    cost_series = [[sub_cst.cost for sub_cst in reverse(sys_cst.subsystem)] for sys_cst in cst]
    perc_series = [["$(round(sub_cst.cost/sys_cst.cost*1000)/10)%" for sub_cst in reverse(sys_cst.subsystem)] for sys_cst in cst]

    size --> (1000, 300)
    cols = 1 + length(filter!(!isempty, cost_series))
    layout := RecipesBase.@layout (1, cols)
    margin --> 5 * Measures.mm

    @series begin
        subplot := 1
        seriestype := :bar
        orientation := :horizontal
        title := "\n" * "Direct Capital Cost in $(costing.construction_start_year) " * @sprintf("[%.3g \$\$B]", sum(costs) / 1E3)
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
        alpha := 0.5
        linecolor := :match
        color := [PlotUtils.palette(:tab10)[c] for c in length(costs):-1:1]
        reverse(names), reverse(costs)
    end

    for (k, (sub_names, sub_perc, sub_costs)) in enumerate(zip(name_series, perc_series, cost_series))
        if !isempty(sub_costs)
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
                ylim := (0, length(sub_costs))
                alpha := 0.5
                linecolor := :match
                color := PlotUtils.palette(:tab10)[k]
                sub_names, sub_costs
            end
        end
    end

end

@recipe function plot_costing(cst::IMAS.costing__cost_direct_capital)
    @series begin
        return cst.system
    end
end

@recipe function plot_costing(cst::IMAS.costing)
    @series begin
        return cst.cost_direct_capital.system
    end
end

# =========== #
# equilibrium #
# =========== #
@recipe function plot_eq(eq::IMAS.equilibrium)
    @series begin
        return eq.time_slice[]
    end
end

@recipe function plot_eqt(eqt::IMAS.equilibrium__time_slice; cx=false, coordinate=:psi_norm, core_profiles_overlay=false)
    @assert typeof(cx) <: Bool
    @assert typeof(coordinate) <: Symbol

    if !cx
        layout := RecipesBase.@layout [a{0.35w} [a b; c d]]
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
                title := L"P~~[MPa]"
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
                    title := L"P~~[MPa]"
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
                title := L"J_{tor}~[MA/m^2]"
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
                    title := L"J_{tor}~[MA/m^2]"
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
                if contains(string(coordinate), "psi")
                    title := L"\rho"
                    eqt.profiles_1d, :rho_tor_norm
                else
                    title := L"\psi~~[Wb]"
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
                title := L"q"
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
    psi_levels_in=nothing,
    psi_levels_out=nothing,
    lcfs=false,
    secondary_separatrix=false,
    show_x_points=false,
    magnetic_axis=true)

    @assert typeof(psi_levels_in) <: Union{Nothing,Int,AbstractVector{<:Real}}
    @assert typeof(psi_levels_out) <: Union{Nothing,Int,AbstractVector{<:Real}}
    @assert typeof(lcfs) <: Bool
    @assert typeof(secondary_separatrix) <: Bool
    @assert typeof(show_x_points) <: Bool
    @assert typeof(magnetic_axis) <: Bool

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
    if lcfs
        psi_levels_in = [psi__boundary_level, psi__boundary_level]
        psi_levels_out = []
    else
        npsi = 11
        if psi_levels_in === nothing
            psi_levels_in = range(eqt.profiles_1d.psi[1], psi__boundary_level, npsi)
        elseif isa(psi_levels_in, Int)
            if psi_levels_in > 1
                npsi = psi_levels_in
                psi_levels_in = range(eqt.profiles_1d.psi[1], psi__boundary_level, psi_levels_in)
            else
                psi_levels_in = []
            end
        end
        delta_psi = (psi__boundary_level - eqt.profiles_1d.psi[1])
        if psi_levels_out === nothing
            psi_levels_out = delta_psi .* range(0.0, 1.0, npsi) .+ psi__boundary_level
        elseif isa(psi_levels_out, Int)
            if psi_levels_out > 1
                psi_levels_out = delta_psi / npsi .* collect(0:psi_levels_out) .+ psi__boundary_level
            else
                psi_levels_out = []
            end
        end
    end
    psi_levels = unique(vcat(psi_levels_in, psi_levels_out))

    for psi_level in psi_levels
        for (pr, pz) in flux_surface(eqt, psi_level, :any)[1]
            @series begin
                seriestype --> :path
                if psi_level == psi__boundary_level
                    linewidth --> 2.0
                elseif psi_level in psi_levels_in
                    linewidth --> 1.0
                else
                    linewidth --> 0.5
                end
                pr, pz
            end
        end
    end

    if secondary_separatrix
        @series begin
            primary --> false
            linewidth --> 1.5
            eqt.boundary_secondary_separatrix.outline.r, eqt.boundary_secondary_separatrix.outline.z
        end
    end

    if magnetic_axis
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

end

@recipe function plot_eqtb(eqtb::IMAS.equilibrium__time_slice___boundary)
    @series begin
        eqtb.outline
    end
    @series begin
        primary --> false
        eqtb.x_point
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
    @series begin
        linewidth --> 2
        aspect_ratio := :equal
        surface_outline.r, surface_outline.z
    end
end

# ===== #
# build #
# ===== #
function join_outlines(r1::AbstractVector{T}, z1::AbstractVector{T}, r2::AbstractVector{T}, z2::AbstractVector{T}) where {T<:Real}
    i1, i2 = minimum_distance_two_shapes(r1, z1, r2, z2; return_index=true)
    r = vcat(reverse(r1[i1:end]), r2[i2:end], r2[1:i2], reverse(r1[1:i1]))
    z = vcat(reverse(z1[i1:end]), z2[i2:end], z2[1:i2], reverse(z1[1:i1]))
    return r, z
end

"""
    plot_build_cx(bd::IMAS.build; cx=true, wireframe=false, only=Symbol[], exclude_layers=Symbol[])

Plot build cross-section
"""
@recipe function plot_build_cx(bd::IMAS.build; cx=true, wireframe=false, only=Symbol[], exclude_layers=Symbol[])

    @assert typeof(cx) <: Bool
    @assert typeof(wireframe) <: Bool
    @assert typeof(only) <: AbstractVector{Symbol}
    @assert typeof(exclude_layers) <: AbstractVector{Symbol}

    legend_position --> :outerbottomright
    aspect_ratio := :equal
    grid --> :none

    # cx
    if cx
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
                    linewidth --> 1
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
                    linewidth --> 1
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
            linewidth --> 1
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
                        linewidth --> 0.5
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
            name = l.name
            color = :gray
            if l.type == Int(_gap_)
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
                color = :lightblue
            elseif l.type == Int(_cryostat_)
                color = :lightgray
            end
            for nm in ("inner", "outer", "vacuum", "hfs", "lfs", "gap")
                name = replace(name, Regex("$(nm) ", "i") => "")
            end

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
                    linewidth --> 0.5
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
                linewidth --> 1.0
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
                        linewidth --> 0.5
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
                            linewidth := 2.0
                            color --> :lightblue
                            label --> (!wireframe && k == 1 ? "Vessel Port" : "")
                            xlim --> [0, maximum(bd.structure[index].outline.r)]
                            bd.structure[index].outline.r, bd.structure[index].outline.z
                        end
                    end
                    @series begin
                        seriestype --> :path
                        linewidth --> 0.5
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
            linewidth --> 2
            label --> ""
            linestyle --> :dash
            color --> :black
            [0]
        end

        at = 0
        for l in bd.layer
            @series begin
                if l.type == Int(_plasma_)
                    color --> :pink
                elseif l.type == Int(_gap_)
                    color --> :white
                elseif l.type == Int(_oh_)
                    if l.side == _in_
                        color --> :gray
                    else
                        color --> :white
                    end
                elseif l.type == Int(_tf_)
                    color --> :green
                elseif l.type == Int(_shield_)
                    color --> :red
                elseif l.type == Int(_blanket_)
                    color --> :orange
                elseif l.type == Int(_wall_)
                    color --> :yellow
                elseif l.type == Int(_vessel_)
                    color --> :lightblue
                elseif l.type == Int(_cryostat_)
                    color --> :lightgray
                end
                seriestype --> :vspan
                if contains(l.name, "gap ")
                    label --> ""
                else
                    label --> l.name
                end
                alpha --> 0.2
                xlim --> [0, at]
                [at, at + l.thickness]
            end
            at += l.thickness
            @series begin
                seriestype --> :vline
                linewidth --> 0.5
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
@recipe function plot_build_tf(tf::IMAS.build__tf; cutouts=false)
    @assert typeof(cutouts) <: Bool
    layers = parent(tf).layer
    TF = get_build_layers(layers; type=IMAS._tf_)

    for n in 1:tf.coils_n
        x, y = top_outline(tf, n; cutouts)
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
@recipe function plot_core_transport(ct::IMAS.core_transport{D}; time0=global_time(ct)) where {D<:Real}
    model_type = name_2_index(ct.model)

    rhos = D[]
    for model in ct.model
        if model.identifier.index ∈ (model_type[k] for k in (:combined, :unspecified, :transport_solver, :unknown))
            continue
        end
        append!(rhos, model.profiles_1d[].grid_flux.rho_tor_norm)
    end
    rhos = unique(rhos)

    dd = top_dd(ct)
    cp1d = dd.core_profiles.profiles_1d[]

    if dd !== nothing
        @series begin
            linewidth := 2
            color := :blue
            label := "Total source"
            flux := true
            total_sources(dd; time0)
        end
    end

    @series begin
        linewidth := 2
        color := :red
        label := "Total transport"
        total_fluxes(ct, cp1d, rhos; time0)
    end

    for model in ct.model
        if model.identifier.index ∈ (model_type[k] for k in (:combined, :unspecified, :transport_solver, :unknown))
            continue
        end
        @series begin
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

@recipe function plot_ct1d(ct1d::IMAS.core_transport__model___profiles_1d; label="", markershape=:none, color=:green, only=nothing)
    if only === nothing
        layout := (2, 2)
        size --> (800, 600)
        margin --> 5 * Measures.mm
    end

    if only === nothing || only == 1
        if !ismissing(ct1d.electrons.energy, :flux)
            @series begin
                if only === nothing
                    subplot := 1
                end
                color := color
                markershape := markershape
                title := "Electron energy flux"
                label := :none

                ct1d.electrons.energy, :flux
            end
        end
    end

    if only === nothing || only == 2
        if !ismissing(ct1d.total_ion_energy, :flux)
            @series begin
                if only === nothing
                    subplot := 2
                end
                color := color
                markershape := markershape
                title := "Ion energy flux"
                label := label

                ct1d.total_ion_energy, :flux
            end
        end
    end

    if only === nothing || only == 3
        if !ismissing(ct1d.electrons.particles, :flux)
            @series begin
                if only === nothing
                    subplot := 3
                end
                color := color
                markershape := markershape
                title := "Electron particle flux"
                label := :none

                ct1d.electrons.particles, :flux
            end
        end
    end

    if only === nothing || only == 4
        if !ismissing(ct1d.momentum_tor, :flux)
            @series begin
                if only === nothing
                    subplot := 4
                end
                color := color
                markershape := markershape
                title := "Toroidal momentum flux"
                label := :none

                ct1d.momentum_tor, :flux
            end
        end
    end
end

# ============ #
# core_sources #
# ============ #
@recipe function plot_core_sources(cs::IMAS.core_sources{T}; time0=global_time(cs), aggregate_radiation=false) where {T<:Real}
    @assert typeof(time0) <: Float64
    @assert typeof(aggregate_radiation) <: Bool

    for source in cs.source
        @series begin
            if aggregate_radiation
                only_positive_negative := 1
            end
            nozeros := true
            time0 := time0
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
                rad_source
            end
        end

        @series begin
            name := "total"
            linewidth := 2
            color := :black
            min_power := 0.0
            total_sources(dd; time0)
        end
    end
end

@recipe function plot_source(source::IMAS.core_sources__source; time0=global_time(source))
    if !isempty(source.profiles_1d) && source.profiles_1d[1].time <= time0
        @series begin
            name := source.identifier.name
            source.profiles_1d[time0]
        end
    end
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
    show_source_number=false
)
    @assert typeof(name) <: AbstractString
    @assert typeof(integrated) <: Bool
    @assert typeof(flux) <: Bool
    @assert typeof(label) <: Union{Nothing,AbstractString}
    @assert typeof(show_zeros) <: Bool
    @assert typeof(min_power) <: Float64
    @assert typeof(only_positive_negative) <: Int
    @assert typeof(show_source_number) <: Bool

    if label === nothing
        label = ""
    end

    if only === nothing
        if flux
            layout := (2, 2)
            size --> (800, 600)
        else
            layout := (1, 4)
            size --> (1100, 290)
            margin --> 5 * Measures.mm
        end
    end

    if parent(cs1d) !== nothing && parent(parent(cs1d)) !== nothing
        source = parent(parent(cs1d))
        idx = index(source)
        if show_source_number
            name = "[$idx] $name"
        end
    else
        source = nothing
        idx = 1
    end
    if source !== nothing
        source_name = identifier_name(source)
    else
        source_name = :undefined
    end

    # electron energy
    if only === nothing || only == 1
        tot = 0.0
        if !ismissing(cs1d.electrons, :energy) && !flux
            tot = trapz(cs1d.grid.volume, cs1d.electrons.energy)
        end
        show_condition =
            flux || show_zeros || source_name in [:collisional_equipartition, :time_derivative] ||
            (abs(tot) > min_power && (only_positive_negative == 0 || sign(tot) == sign(only_positive_negative)))
        @series begin
            if only === nothing
                subplot := 1
            end
            if source_name == :collisional_equipartition
                linestyle --> :dash
            end
            color := idx
            title --> "Electron Energy"
            if show_condition
                label := "$name " * @sprintf("[%.3g MW]", tot / 1E6) * label
                if !ismissing(cs1d.electrons, :power_inside) && flux
                    label := :none
                    cs1d.grid.rho_tor_norm[2:end], (cs1d.electrons.power_inside./cs1d.grid.surface)[2:end]
                elseif !integrated && !ismissing(cs1d.electrons, :energy)
                    if source_name in [:ec, :ic, :lh, :nbi, :pellet]
                        fill0 --> true
                    end
                    cs1d.electrons, :energy
                elseif integrated && !ismissing(cs1d.electrons, :power_inside)
                    cs1d.electrons, :power_inside
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

    # ion energy
    if only === nothing || only == 2
        tot = 0.0
        if !ismissing(cs1d, :total_ion_energy)
            tot = trapz(cs1d.grid.volume, cs1d.total_ion_energy)
        end
        show_condition = flux || show_zeros || source_name in [:collisional_equipartition, :time_derivative] || abs(tot) > min_power
        @series begin
            if only === nothing
                subplot := 2
            end
            if source_name == :collisional_equipartition
                linestyle --> :dash
            end
            color := idx
            title --> "Ion Energy"
            if show_condition
                label := "$name " * @sprintf("[%.3g MW]", tot / 1E6) * label
                if !ismissing(cs1d, :total_ion_power_inside) && flux
                    cs1d.grid.rho_tor_norm[2:end], (cs1d.total_ion_power_inside./cs1d.grid.surface)[2:end]
                elseif !integrated && !ismissing(cs1d, :total_ion_energy)
                    if source_name in [:ec, :ic, :lh, :nbi, :pellet]
                        fill0 --> true
                    end
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

    # particles
    if only === nothing || only == 3
        tot = 0.0
        if !ismissing(cs1d.electrons, :particles)
            tot = trapz(cs1d.grid.volume, cs1d.electrons.particles)
        end
        show_condition = flux || show_zeros || source_name == :time_derivative || abs(tot) > 0.0
        @series begin
            if only === nothing
                subplot := 3
            end
            color := idx
            title --> "Electron Particle"
            if show_condition
                label := "$name " * @sprintf("[%.3g s⁻¹]", tot) * label
                if !ismissing(cs1d.electrons, :particles_inside) && flux
                    label := :none
                    cs1d.grid.rho_tor_norm[2:end], (cs1d.electrons.particles_inside./cs1d.grid.surface)[2:end]
                elseif !integrated && !ismissing(cs1d.electrons, :particles)
                    if source_name in [:ec, :ic, :lh, :nbi, :pellet]
                        fill0 --> true
                    end
                    cs1d.electrons, :particles
                elseif integrated && !ismissing(cs1d.electrons, :particles_inside)
                    cs1d.electrons, :particles_inside
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

    # current (or momentum, if plotting flux)
    if only === nothing || only == 4
        if flux
            if only === nothing
                subplot := 4
            end
            color := idx
            title --> "Momentum Tor"
            if !ismissing(cs1d, :torque_tor_inside)
                label := :none
                cs1d.grid.rho_tor_norm[2:end], (cs1d.torque_tor_inside./cs1d.grid.surface)[2:end]
            else
                label := ""
                [NaN], [NaN]
            end
        else
            tot = 0.0
            if !ismissing(cs1d, :j_parallel)
                tot = trapz(cs1d.grid.area, cs1d.j_parallel)
            end
            show_condition = flux || show_zeros || abs(tot) > 0.0
            @series begin
                if only === nothing
                    subplot := 4
                end
                color := idx
                title --> "Parallel Current"
                if show_condition
                    label := "$name " * @sprintf("[%.3g MA]", tot / 1E6) * label
                    if !integrated && !ismissing(cs1d, :j_parallel)
                        if source_name in [:ec, :ic, :lh, :nbi, :pellet]
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
    end
end

# ============= #
# core_profiles #
# ============= #
@recipe function plot_core_profiles(cp::IMAS.core_profiles)
    @series begin
        return cp.profiles_1d[]
    end
end

@recipe function plot_core_profiles(cpt::IMAS.core_profiles__profiles_1d; label=nothing, only=nothing, greenwald=false)

    @assert typeof(label) <: Union{Nothing,AbstractString}
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
            linewidth := 0.75
            label := "Yield strength"
            r_oh, r_oh .* 0.0 .+ smcs.properties.yield_strength.oh / 1E6
        end
        @series begin
            primary := false
            linestyle := :dash
            linewidth := 0.75
            r_tf, r_tf .* 0.0 .+ smcs.properties.yield_strength.tf / 1E6
        end
        if r_pl !== missing
            @series begin
                primary := false
                linestyle := :dash
                linewidth := 0.75
                r_pl, r_pl .* 0.0 .+ smcs.properties.yield_strength.pl / 1E6
            end
        end
    end
end

@recipe function plot_solid_mechanics(stress::IMAS.solid_mechanics__center_stack__stress; linewidth=1)
    @assert typeof(linewidth) <: Real

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
        linewidth := linewidth + 2
        linestyle := :solid
        stress.vonmises
    end
    @series begin
        label := "Hoop" * config
        linewidth := linewidth + 1
        linestyle := :dash
        stress.hoop
    end
    @series begin
        label := "Radial" * config
        linewidth := linewidth + 1
        linestyle := :dashdot
        stress.radial
    end

    smcs = parent(stress)
    for radius in (smcs.grid.r_oh[1], smcs.grid.r_oh[end], smcs.grid.r_tf[1], smcs.grid.r_tf[end])
        @series begin
            seriestype := :vline
            linewidth := linewidth
            label := ""
            linestyle := :dash
            color := :black
            [radius]
        end
    end
end

# ================ #
# balance_of_plant #
# ================ #
@recipe function plot_balance_of_plant(bop::IMAS.balance_of_plant; linewidth=2)
    @assert typeof(linewidth) <: Real

    size --> (800, 600)
    legend_position --> :outertopright
    ylabel --> "Electricity [Watts Electric]"
    xlabel --> "Time [s]"

    @series begin
        label := "Net electric"
        linewidth := linewidth + 2
        color := "Black"
        bop, :power_electric_net
    end

    @series begin
        label := "Electricity generated"
        linewidth := linewidth + 1
        linestyle --> :dash
        color := "Black"
        bop.power_plant, :power_electric_generated
    end

    for sys in bop.power_electric_plant_operation.system
        @series begin
            label := string(sys.name)
            linewidth := linewidth
            sys, :power
        end
    end
end

# ========== #
# neutronics #
# ========== #
@recipe function plot_neutron_wall_loading_cx(nwl::IMAS.neutronics__time_slice___wall_loading, component::Symbol=:norm; cx=true)

    @assert typeof(cx) <: Bool

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
        linewidth --> 8
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
    @series begin
        color := :gray
        return wall.description_2d
    end
end

@recipe function plot_wd2d_list(wd2d_list::IDSvector{<:IMAS.wall__description_2d})
    label := ""
    for wd2d in wd2d_list
        @series begin
            return wd2d
        end
    end
end

@recipe function plot_wd2d(wd2d::IMAS.wall__description_2d; show_limiter=true, show_vessel=true)
    @assert typeof(show_limiter) <: Bool
    @assert typeof(show_vessel) <: Bool
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

"""
    thick_line_polygon(r1, z1, r2, z2, thickness1, thickness2)

Casts thick line as a polygon. Returns points of the quadrilateral polygon
"""
function thick_line_polygon(r1::Float64, z1::Float64, r2::Float64, z2::Float64, thickness1::Float64, thickness2::Float64)
    direction = normalize([z2 - z1, -(r2 - r1)]) # Perpendicular direction
    offset1 = direction .* thickness1 / 2
    offset2 = direction .* thickness2 / 2
    p1 = [r1, z1] + offset1
    p2 = [r2, z2] + offset2
    p3 = [r2, z2] - offset2
    p4 = [r1, z1] - offset1
    return [p1, p2, p3, p4, p1]
end

@recipe function plot_wd2dvu_list(wd2dvu_list::IDSvector{<:IMAS.wall__description_2d___vessel__unit})
    for wd2dvu in wd2dvu_list
        @series begin
            return wd2dvu
        end
    end
end

@recipe function plot_wd2dvu(wd2dvu::IMAS.wall__description_2d___vessel__unit)
    centreline = closed_polygon(wd2dvu.annular.centreline.r, wd2dvu.annular.centreline.z, Bool(wd2dvu.annular.centreline.closed))

    if !ismissing(wd2dvu.annular, :thickness)
        thickness = wd2dvu.annular.thickness
        if Bool(wd2dvu.annular.centreline.closed)
            thickness = [thickness; thickness[1]]
        end

        for i in 1:length(centreline.R)-1
            pts = thick_line_polygon(centreline.R[i], centreline.Z[i], centreline.R[i+1], centreline.Z[i+1], thickness[i], thickness[i+1])
            # Extract x and y coordinates for plotting
            xs, ys = map(collect, zip(pts...))
            @series begin
                primary := i == 1
                linewidth := 0.1
                aspect_ratio := :equal
                fill := (0, 0.5, :gray)
                xlim := (0, Inf)
                xs, ys
            end
        end
    else
        @series begin
            label --> getproperty(wd2dvu, :name, "")
            aspect_ratio := :equal
            centreline.R, centreline.Z
        end
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
    @series begin
        aspect_ratio := :equal
        poly.R, poly.Z
    end
end

#= ============== =#
#  pulse_schedule  #
#= ============== =#
@recipe function plot_ps(ps::IMAS.pulse_schedule; time0=global_time(ps))
    @assert typeof(time0) <: Float64
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
        substitute = ("deposition_" => "", "_launched" => "")
        plt[:label] = replace(p2i(filter(x -> x ∉ remove, path[2:end])), substitute...)
        push!(plots, plt)
    end

    layout := RecipesBase.@layout [length(plots) + 1]
    size --> (1000, 1000)
    margin --> 5 * Measures.mm

    @series begin
        subplot := 1
        label := "$(time) [s]"
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

    for (k, plt) in enumerate(plots)
        x = plt[:x]
        y = plt[:y]
        n = min(length(x), length(y))
        @series begin
            subplot := k + 1
            label := ""
            x[1:n], y[1:n]
        end
        @series begin
            subplot := k + 1
            seriestype := :scatter
            primary := false
            marker := :circle
            markerstrokewidth := 0.0
            title := nice_field(plt[:label])
            [time0], interp1d(x[1:n], y[1:n]).([time0])
        end
    end
end

@recipe function plot_pc_time(pc::IMAS.pulse_schedule__position_control; time0=global_time(pc))
    @assert typeof(time0) <: Float64
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

#= ================ =#
#  generic plotting  #
#= ================ =#
@recipe function plot_field(ids::IMAS.IDS, field::Symbol; normalization=1.0, coordinate=nothing, weighted=:none, fill0=false)
    @assert hasfield(typeof(ids), field) "$(location(ids)) does not have field `$field`. Did you mean: $(keys(ids))"
    @assert typeof(normalization) <: Real
    @assert typeof(coordinate) <: Union{Nothing,Symbol}
    @assert typeof(weighted) <: Symbol
    @assert typeof(fill0) <: Bool

    coords = coordinates(ids, field; coord_leaves=[coordinate])
    coordinate_name = coords.names[1]
    coordinate_value = coords.values[1]

    xvalue = coordinate_value
    yvalue = getproperty(ids, field) .* normalization

    @series begin
        background_color_legend := PlotUtils.Colors.RGBA(1.0, 1.0, 1.0, 0.6)

        if endswith(coordinate_name, "_norm")
            xlim --> (0.0, 1.0)
        end

        # multiply y by things like `:area` or `:volume`
        if weighted != :none
            weight = coordinates(ids, field; coord_leaves=[weighted])
            yvalue .*= weight.values[1]
            ylabel = nice_units(units(ids, field) * "*" * units(weight.names[1]))
            label = nice_field("$field*$weighted")
        else
            ylabel = nice_units(units(ids, field))
            label = nice_field(field)
        end

        xlabel --> nice_field(i2p(coordinate_name)[end]) * nice_units(units(coordinate_name))
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

@recipe function plot(x::AbstractVector{<:Real}, y::AbstractVector{<:Measurements.Measurement}, err::Symbol=:ribbon)
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

#= ================== =#
#  handling of labels  #
#= ================== =#
nice_field_symbols = Dict()
nice_field_symbols["rho_tor_norm"] = L"\rho"
nice_field_symbols["psi"] = L"\psi"
nice_field_symbols["psi_norm"] = L"\psi_N"
nice_field_symbols["rotation_frequency_tor_sonic"] = "Rotation"
nice_field_symbols["i_plasma"] = "Plasma current"
nice_field_symbols["b_field_tor_vacuum_r"] = "B₀×R₀"
nice_field_symbols["rho_tor_norm_width"] = "w₀"
nice_field_symbols["geometric_axis.r"] = "Rgeo"
nice_field_symbols["geometric_axis.z"] = "Zgeo"

function nice_field(field::AbstractString)
    if field in keys(nice_field_symbols)
        field = nice_field_symbols[field]
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
        units = L"[%$units]"
        units = " " * units
    end
    return units
end

#= ============= =#
#  Thermal loads #
#= ============= =#
"""
Recipe for plot of heat flux

  - which_plot = :twoD, :oneD
  - plot_type  = :path, :scatter (only for 2D)
  - q          =
    :all (for 1D),
    :wall,
    :parallel
    :particle
    :core_radiation
    :both (= :particle + :parallel)
"""
@recipe function plot_heat_flux(HF::WallHeatFlux; which_plot=:twoD, plot_type=:path, q=:wall)
    @assert which_plot in (:oneD, :twoD)
    @assert plot_type in (:path, :scatter)
    @assert q in (:wall, :parallel, :particle, :core_radiation, :both, :all)

    if q in (:both, :all)
        if q == :both
            qs = (:parallel, :particle)
        else
            qs = (:wall, :particle, :core_radiation)
        end
        if which_plot == :twoD
            layout := (1, length(qs))
        end
        for (k, q) in enumerate(qs)
            @series begin
                if which_plot == :twoD
                    subplot := k
                end
                which_plot := which_plot
                plot_type := plot_type
                q := q
                HF
            end
        end

    elseif which_plot == :oneD
        @series begin
            xlabel --> "Clockwise distance along wall [m]"
            ylabel --> "Wall flux [W/m²]"
            yscale --> :log10

            x = []
            y = []
            if q == :particle
                label --> "particles"
                if !isempty(HF.q_part)
                    x = HF.s
                    y = HF.q_part .+ 1.0
                end

            elseif q == :core_radiation
                label --> "core radiation"
                if !isempty(HF.q_core_rad)
                    x = HF.s
                    y = HF.q_core_rad .+ 1.0
                end

            elseif q == :wall
                label --> "wall"
                if !isempty(HF.q_wall)
                    x = HF.s
                    y = HF.q_wall .+ 1.0
                end

            elseif q == :parallel
                label --> "parallel"
                if !isempty(HF.q_parallel)
                    x = HF.s
                    y = HF.q_parallel .+ 1.0
                end
            end

            x, y
        end

    elseif which_plot == :twoD
        @series begin
            aspect_ratio := :equal
            legend --> false
            colorbar --> true
            xlim --> [0.95 * minimum(HF.r) - 0.05 * (maximum(HF.r)), 1.05 * maximum(HF.r) - 0.05 * (minimum(HF.r))]
            ylim --> [1.05 * minimum(HF.z) - 0.05 * (maximum(HF.z)), 1.05 * maximum(HF.z) - 0.05 * (minimum(HF.z))]
            if plot_type == :path
                seriestype --> :path
                linewidth --> 8
            elseif plot_type == :scatter
                seriestype --> :scatter
                markersize --> 2
                markerstrokewidth --> 0
            end

            if q == :wall && !isempty(HF.q_wall)
                colorbar_title := "log₁₀(q wall [W/m²])"
                if plot_type == :path
                    line_z := log10.(HF.q_wall .+ 1)
                elseif plot_type == :scatter
                    zcolor := log10.(HF.q_wall .+ 1)
                end

            elseif q == :core_radiation && !isempty(HF.q_core_rad)
                colorbar_title := "log₁₀(q core rad [W/m²])"
                if plot_type == :path
                    line_z := log10.(HF.q_core_rad .+ 1)
                elseif plot_type == :scatter
                    zcolor := log10.(HF.q_core_rad .+ 1)
                end

            elseif q == :particle && !isempty(HF.q_part)
                colorbar_title := "log₁₀(q particle [W/m²])"
                floor = minimum(HF.q_part[HF.q_part.>0.0]) / 100.0
                if plot_type == :path
                    line_z --> log10.(HF.q_part .+ floor)
                elseif plot_type == :scatter
                    zcolor --> log10.(HF.q_part .+ floor)
                end

            elseif q == :parallel && !isempty(HF.q_parallel)
                colorbar_title := "log₁₀(q parallel [W/m²])"
                if plot_type == :path
                    line_z --> log10.(HF.q_parallel .+ 1)
                end
                if plot_type == :scatter
                    zcolor --> log10.(HF.q_parallel .+ 1)
                end
            end

            HF.r, HF.z
        end
    end
end