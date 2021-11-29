using Plots

"""
    plot_pf_active_cx(pfa::pf_active)

Plots pf active cross-section
"""
@recipe function plot_pf_active_cx(pfa::pf_active)
    markershape --> :rect
    label --> ""
    aspect --> :equal
    seriestype --> :scatter
    r = [c.element[1].geometry.rectangle.r for c in pfa.coil]
    z = [c.element[1].geometry.rectangle.z for c in pfa.coil]
    if ! any([is_missing(c.current,:data) for c in pfa.coil])
        currents = [sum(c.current.data)/length(c.current.data) for c in pfa.coil]
        color --> :roma
        clim --> (-maximum(abs.(currents)),maximum(abs.(currents)))
        @series begin
            marker_z --> -currents
            [(r, z)]
        end
        @series begin
            marker_z --> currents
            [(r, z)]
        end
    else
        @series begin
            [(r, z)]
        end
    end
end


"""
    plot_eqcx(eqt::IMAS.equilibrium__time_slice; psi_levels=nothing, psi_levels_out=nothing, lcfs=false)

Plots equilibrium cross-section
"""
@recipe function plot_eqcx(eqt::IMAS.equilibrium__time_slice; psi_levels=nothing, psi_levels_out=nothing, lcfs=false)

    # do not trust eqt.profiles_1d.psi[end], and find boundary level that is closest to lcfs
    psi__boundary_level = eqt.profiles_1d.psi[end]
    tmp = find_psi_boundary(eqt, raise_error_on_not_open=false)
    if tmp !== nothing
        psi__boundary_level = tmp
    end

    if lcfs
        psi_levels = [psi__boundary_level, psi__boundary_level]
        psi_levels_out = []
    else
        if psi_levels === nothing
            psi_levels = range(eqt.profiles_1d.psi[1], psi__boundary_level, length=11)
        elseif isa(psi_levels, Int)
            psi_levels = range(eqt.profiles_1d.psi[1], psi__boundary_level, length=psi_levels)
        end
        if psi_levels_out === nothing
            psi_levels_out = (psi_levels[end] - psi_levels[1]) .* collect(range(0, 1, length=11)) .+ psi_levels[end]
        elseif isa(psi_levels_out, Int)
            psi_levels_out = (psi_levels[end] - psi_levels[1]) .* collect(range(0, 1, length=psi_levels_out)) .+ psi_levels[end]
        end
    end
    psi_levels = unique(vcat(psi_levels, psi_levels_out))

    label --> ""
    aspect_ratio --> :equal
    primary --> false

    xlims --> eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end]
    ylims --> eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end]

    # @series begin
    #     seriestype --> :contour
    #     levels --> psi_levels
    #     eqt.profiles_2d[1].grid.dim1, eqt.profiles_2d[1].grid.dim2, transpose(eqt.profiles_2d[1].psi)
    # end
    for psi_level in psi_levels
        for (pr, pz) in IMAS.flux_surface(eqt, psi_level, false)
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

function join_outlines(r1, z1, r2, z2)
        i1 = 1
    i2 = argmin((r2 .- r1[i1]).^2 + (z2 .- z1[i1]).^2)
    p(x1, x2) = vcat(x1, x2[i2:end], x2[1:i2], x1[i1])
    return p(r1, r2), p(z1, z2)
end

"""
    plot_radial_build_cx(rb::IMAS.radial_build)

Plots radial build cross-section
"""
@recipe function plot_radial_build_cx(rb::IMAS.radial_build; outline=true)
    aspect_ratio --> :equal
    grid --> :none

    @series begin
        seriestype --> :vline
        linewidth --> 2
        label --> "center"
        linestyle --> :dash
        color --> :black
        [0]
    end

    # outline
    if outline
        rmax = maximum(rb.layer[end].outline.r) * 2.0

        # Cryostat
        @series begin
            seriestype --> :shape
        linewidth --> 0.0
            color --> :white
            label --> ""
            xlim --> [0,rmax]
            join_outlines(rb.layer[end].outline.r, rb.layer[end].outline.z, IMAS.get_radial_build(rb, type=-1).outline.r, IMAS.get_radial_build(rb, type=-1).outline.z)
        end

        @series begin
            seriestype --> :path
            linewidth --> 1
            color --> :black
            label --> "Cryostat"
            xlim --> [0,rmax]
            rb.layer[end].outline.r[1:end - 1], rb.layer[end].outline.z[1:end - 1]
        end

        # OH
        @series begin
            seriestype --> :shape
            linewidth --> 0.5
            color --> :blue
            label --> IMAS.get_radial_build(rb, type=1).name
            xlim --> [0,rmax]
            IMAS.get_radial_build(rb, type=1).outline.r, IMAS.get_radial_build(rb, type=1).outline.z
        end

        # all layers between the OH and the vessel
        valid = false
        for (k, l) in enumerate(rb.layer[1:end - 1])
            if (l.type == 2) && (l.hfs == 1)
                valid = true
            end
            if l.type == -1
                valid = false
            end
            if IMAS.is_missing(l.outline, :r) || ! valid
                continue
            end
            l1 = rb.layer[k + 1]
            poly = join_outlines(l.outline.r, l.outline.z, l1.outline.r, l1.outline.z)

            # setup labels and colors
            name = l.name
            color = :gray
            if ! is_missing(l, :material) && l.material == "vacuum"
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

            @series begin
                seriestype --> :shape
                linewidth --> 0.0
                color --> color
                label --> uppercasefirst(name)
                xlim --> [0,rmax]
                poly[1], poly[2]
            end
            @series begin
                seriestype --> :path
                linewidth --> 0.5
                color --> :black
                label --> ""
                xlim --> [0,rmax]
                l.outline.r, l.outline.z
            end
        end

        @series begin
            seriestype --> :path
            linewidth --> 2
            color --> :black
            label --> "Vessel"
            xlim --> [0, rmax]
            IMAS.get_radial_build(rb, type=-1).outline.r, IMAS.get_radial_build(rb, type=-1).outline.z
        end

    # not-outlines
    else

        at = 0
        for l in rb.layer
            @series begin
                if ! is_missing(l, :material) && l.material == "vacuum"
                    color --> :white
                elseif occursin("OH", l.name)
                    color --> :blue
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
                xlim --> [0,at * 2.0]
                [at,at + l.thickness]
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