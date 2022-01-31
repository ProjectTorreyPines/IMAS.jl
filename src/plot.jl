using Plots
using LaTeXStrings

"""
    plot_pf_active_cx(pfa::pf_active)

Plots pf active cross-section
"""
@recipe function plot_pf_active_cx(pfa::pf_active, what::Symbol = :cx; cname = :roma, time_index = 1)

    if what in [:cx, :coils_flux]
        label --> ""
        aspect --> :equal
        colorbar_title --> "PF currents [A]"

        currents = [c.current.data[time_index] for c in pfa.coil]
        CURRENT = maximum(abs.(currents))

        # dummy markers to get the colorbar right
        @series begin
            seriestype --> :scatter
            color --> cname
            clim --> (-CURRENT, CURRENT)
            marker_z --> [-CURRENT, CURRENT]
            [(NaN, NaN), (NaN, NaN)]
        end

        # plot individual coils
        for c in pfa.coil
            current_color_index = (c.current.data[time_index] + CURRENT) / (2 * CURRENT)
            @series begin
                color --> cgrad(cname)[current_color_index]
                c
            end
        end

    elseif what == :currents
        label --> "$(pfa.coil[1].current.time[time_index]) s"
        currents = [c.current.data[time_index] for c in pfa.coil]
        @series begin
            linestyle --> :dash
            marker --> :circle
            ["$k" for k = 1:length(currents)], currents
        end

    else
        error("IMAS.pf_active `what` to plot can only be :cx or :currents")
    end

end

"""
    plot_pf_active__coil_cx(coil::pf_active__coil; color::Symbol=:gray)

Plots cross-section of individual coils
"""
@recipe function plot_pf_active__coil_cx(coil::pf_active__coil; color = :gray)
    if (coil.element[1].geometry.rectangle.width == 0.0) || (coil.element[1].geometry.rectangle.height == 0.0)
        @series begin
            color --> color
            linewidth --> 0.0
            seriestype --> :scatter
            markershape --> :star
            colorbar --> :right
            [(coil.element[1].geometry.rectangle.r, coil.element[1].geometry.rectangle.z)]
        end
    else
        r = coil.element[1].geometry.rectangle.r
        z = coil.element[1].geometry.rectangle.z
        Δr = coil.element[1].geometry.rectangle.width / 2.0
        Δz = coil.element[1].geometry.rectangle.height / 2.0
        @series begin
            seriestype --> :shape
            linewidth --> 0.5
            colorbar --> :right
            color --> color
            label --> ""
            [-Δr, Δr, Δr, -Δr, -Δr] .+ r, [-Δz, -Δz, Δz, Δz, -Δz] .+ z
        end
    end
end

"""
    plot_eqcx(eqt::IMAS.equilibrium__time_slice; psi_levels=nothing, psi_levels_out=nothing, lcfs=false)

Plots equilibrium cross-section
"""
@recipe function plot_eqcx(eq::IMAS.equilibrium; psi_levels = nothing, psi_levels_out = nothing, lcfs = false)
    @series begin
        psi_levels --> psi_levels
        psi_levels_out --> psi_levels_out
        lcfs --> lcfs
        return eq.time_slice[]
    end
end

@recipe function plot_eqtcx(eqt::IMAS.equilibrium__time_slice; psi_levels = nothing, psi_levels_out = nothing, lcfs = false)

    label --> ""
    aspect_ratio --> :equal
    primary --> false

    # if there is no psi map then plot the boundary
    if (length(eqt.profiles_2d) == 0) || is_missing(eqt.profiles_2d[1], :psi)
        @series begin
            eqt.boundary.outline.r, eqt.boundary.outline.z
        end

        # plot cx
    else
        # handle psi levels
        psi__boundary_level = eqt.profiles_1d.psi[end]
        tmp = find_psi_boundary(eqt, raise_error_on_not_open = false) # do not trust eqt.profiles_1d.psi[end], and find boundary level that is closest to lcfs
        if tmp !== nothing
            psi__boundary_level = tmp
        end
        if lcfs
            psi_levels = [psi__boundary_level, psi__boundary_level]
            psi_levels_out = []
        else
            if psi_levels === nothing
                psi_levels = range(eqt.profiles_1d.psi[1], psi__boundary_level, length = 11)
            elseif isa(psi_levels, Int)
                if psi_levels > 1
                    psi_levels = range(eqt.profiles_1d.psi[1], psi__boundary_level, length = psi_levels)
                else
                    psi_levels = []
                end
            end
            if psi_levels_out === nothing
                psi_levels_out = (psi__boundary_level - eqt.profiles_1d.psi[1]) .* collect(range(0, 1, length = 11)) .+ psi__boundary_level
            elseif isa(psi_levels_out, Int)
                if psi_levels_out > 1
                    psi_levels_out = (psi__boundary_level - eqt.profiles_1d.psi[1]) .* collect(range(0, 1, length = psi_levels_out)) .+ psi__boundary_level
                else
                    psi_levels_out = []
                end
            end
        end
        psi_levels = unique(vcat(psi_levels, psi_levels_out))

        xlims --> eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end]
        ylims --> eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end]

        # @series begin
        #     seriestype --> :contour
        #     levels --> psi_levels
        #     eqt.profiles_2d[1].grid.dim1, eqt.profiles_2d[1].grid.dim2, transpose(eqt.profiles_2d[1].psi)
        # end
        for psi_level in psi_levels
            for (pr, pz) in IMAS.flux_surface(eqt, psi_level, nothing)
                @series begin
                    seriestype --> :path
                    if psi_level == psi__boundary_level
                        linewidth --> 2
                    else
                        linewidth --> 1
                    end
                    pr, pz
                end
            end
        end

        @series begin
            seriestype --> :scatter
            markershape --> :cross
            [(eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z)]
        end
    end
end

function join_outlines(r1, z1, r2, z2)
    i1 = 1
    i2 = argmin((r2 .- r1[i1]) .^ 2 + (z2 .- z1[i1]) .^ 2)
    p(x1, x2) = vcat(x1, x2[i2:end], x2[1:i2], x1[i1])
    return p(r1, r2), p(z1, z2)
end

"""
    plot_build_cx(bd::IMAS.build)

Plots build cross-section
"""
@recipe function plot_build_cx(bd::IMAS.build; cx = true, outlines = false, only_layers = nothing, exclude_layers = Symbol[])
    aspect_ratio --> :equal
    grid --> :none

    # cx
    if cx
        rmax = maximum(bd.layer[end].outline.r) * 2.0

        # Cryostat
        if ((only_layers === nothing) || (:cryostat in only_layers)) && (!(:cryostat in exclude_layers))
            if !outlines
                @series begin
                    seriestype --> :shape
                    linewidth --> 0.0
                    color --> :white
                    label --> ""
                    xlim --> [0, rmax]
                    join_outlines(
                        bd.layer[end].outline.r,
                        bd.layer[end].outline.z,
                        IMAS.get_build(bd, type = -1).outline.r,
                        IMAS.get_build(bd, type = -1).outline.z,
                    )
                end
            end

            @series begin
                seriestype --> :path
                linewidth --> 1
                color --> :black
                label --> (!outlines ? "Cryostat" : "")
                xlim --> [0, rmax]
                bd.layer[end].outline.r[1:end-1], bd.layer[end].outline.z[1:end-1]
            end

            @series begin
                seriestype --> :path
                linewidth --> 1
                label --> ""
                linestyle --> :dash
                color --> :black
                [0.0, 0.0], [minimum(bd.layer[end].outline.z), maximum(bd.layer[end].outline.z)]
            end
        end

        # OH
        if ((only_layers === nothing) || (:oh in only_layers)) && (!(:oh in exclude_layers))
            if !outlines
                @series begin
                    seriestype --> :shape
                    linewidth --> 0.0
                    color --> :gray
                    label --> (!outlines ? IMAS.get_build(bd, type = 1).name : "")
                    xlim --> [0, rmax]
                    IMAS.get_build(bd, type = 1).outline.r, IMAS.get_build(bd, type = 1).outline.z
                end
            end
            @series begin
                seriestype --> :path
                linewidth --> 0.5
                color --> :black
                label --> ""
                xlim --> [0, rmax]
                IMAS.get_build(bd, type = 1).outline.r, IMAS.get_build(bd, type = 1).outline.z
            end
        end

        # all layers between the OH and the vessel
        valid = false
        for (k, l) in enumerate(bd.layer[1:end-1])
            if (l.type == 2) && (l.hfs == 1)
                valid = true
            end
            if l.type == -1
                valid = false
            end
            if IMAS.is_missing(l.outline, :r) || !valid
                continue
            end
            l1 = bd.layer[k+1]
            poly = join_outlines(l.outline.r, l.outline.z, l1.outline.r, l1.outline.z)

            # setup labels and colors
            name = l.name
            color = :gray
            if !is_missing(l, :material) && l.material == "vacuum"
                name = ""
                color = :white
            elseif occursin("TF", l.name)
                color = :green
            elseif occursin("shield", l.name)
                color = :red
            elseif occursin("blanket", l.name)
                color = :orange
            elseif occursin("wall", l.name)
                color = :yellow
            end
            for nm in ["inner", "outer", "vacuum", "hfs", "lfs"]
                name = replace(name, "$nm " => "")
            end

            if ((only_layers === nothing) || (Symbol(name) in only_layers)) && (!(Symbol(name) in exclude_layers))
                if !outlines
                    @series begin
                        seriestype --> :shape
                        linewidth --> 0.0
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

        if ((only_layers === nothing) || (:vessel in only_layers)) && (!(:vessel in exclude_layers))
            @series begin
                seriestype --> :path
                linewidth --> 2
                color --> :black
                label --> (!outlines ? "Vessel" : "")
                xlim --> [0, rmax]
                IMAS.get_build(bd, type = -1).outline.r, IMAS.get_build(bd, type = -1).outline.z
            end
        end

        # not-cx
    else

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
                if !is_missing(l, :material) && l.material == "vacuum"
                    color --> :white
                elseif occursin("OH", l.name)
                    color --> :gray
                elseif occursin("TF", l.name)
                    color --> :green
                elseif occursin("shield", l.name)
                    color --> :red
                elseif occursin("blanket", l.name)
                    color --> :orange
                elseif occursin("wall", l.name)
                    color --> :yellow
                end
                seriestype --> :vspan
                label --> l.name
                alpha --> 0.2
                xlim --> [0, at * 2.0]
                [at, at + l.thickness]
            end
            at += l.thickness
            @series begin
                seriestype --> :vline
                linewidth --> 0.5
                label --> ""
                color --> :black
                xlim --> [0, at * 2.0]
                [at]
            end
        end
    end
end

@recipe function plot_cp(cp::IMAS.core_profiles)
    @series begin
        return cp.profiles_1d[]
    end
end

@recipe function plot_cp(cpt::IMAS.core_profiles__profiles_1d)
    layout := (1,3)
    size := (1100, 290)

    @series begin
        subplot := 1
        label --> "Electrons"
        cpt.electrons, :temperature
    end

    @series begin
        subplot := 1
        label --> "Ions"
        linestyle --> :dash
        cpt.electrons, :temperature
    end

    @series begin
        subplot := 2
        label := ""
        cpt, :rotation_frequency_tor_sonic
    end

    @series begin
        subplot := 3
        label := "e"
        cpt.electrons, :density
    end

    for ion in cpt.ion
        @series begin
            subplot := 3
            label --> ion.label
            ion, :density
        end
    end

end

#= ================ =#
#  generic plotting  #
#= ================ =#
@recipe function plot_field(ids::IMAS.IDS, field::Symbol)
    coords = coordinates(ids, field)
    @series begin
        xlabel --> nice_field(i2p(coords[:names][1])[end])*nice_units(units(coords[:names][1]))
        ylabel --> nice_units(units(ids, field))
        title --> nice_field(field)
        coords[:values][1], getproperty(ids, field)
    end
end

#= ================== =#
#  handling of labels  #
#= ================== =#
nice_field_symbols = Dict()
nice_field_symbols["rho_tor_norm"] = L"\rho"
nice_field_symbols["psi"] = L"\psi"
nice_field_symbols["rotation_frequency_tor_sonic"] = "Rotation"

function nice_field(field::String)
    if field in keys(nice_field_symbols)
        field = nice_field_symbols[field]
    else
        field = uppercasefirst(replace(field, "_" => " "))
    end
    return field
end

function nice_field(field::Symbol)
    nice_field(String(field))
end

function nice_units(units::String)
    if units == "-"
        units = ""
    end
    if length(units) > 0
        units=replace(units, r"\^([-+]?[0-9]+)" => s"^{\1}")
        units = L"[%$units]"
        units =" "*units
    end
    return units
end