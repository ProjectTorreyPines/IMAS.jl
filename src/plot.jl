using Plots

"""
    function plot_pf_active_cx(pfa::pf_active)

Plots pf active cross-section
"""
@recipe function plot_pf_active_cx(pfa::pf_active)
    markershape --> :rect
    label --> ""
    aspect --> :equal
    seriestype --> :scatter
    r = [c.element[1].geometry.rectangle.r for c in pfa.coil]
    z = [c.element[1].geometry.rectangle.z for c in pfa.coil]
    [(r, z)]
end


"""
    function plot_eqcx(eqt::IMAS.equilibrium__time_slice; lcfs=false)

Plots equilibrium cross-section
"""
@recipe function plot_eqcx(eqt::IMAS.equilibrium__time_slice; lcfs=false)
    if lcfs
        psi_levels = [eqt.global_quantities.psi_boundary, eqt.global_quantities.psi_boundary]
        psi_levels_out = []
        
    else
        psi_levels = range(eqt.global_quantities.psi_axis, eqt.global_quantities.psi_boundary, length=11)
        psi_levels_out = (psi_levels[end] - psi_levels[1]) .* collect(range(0, 1, length=11)) .+ psi_levels[end]
    end

    seriestype --> :path
    label --> ""
    aspect_ratio --> :equal
    primary --> false

    xlims --> eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end]
    ylims --> eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end]

    for psi_level in vcat(psi_levels[2:end], psi_levels_out)
        for (pr, pz) in IMAS.flux_surface(eqt, psi_level, false)
            @series begin
                if psi_level == eqt.global_quantities.psi_boundary
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

"""
    function plot_radial_build_cx(rb::IMAS.radial_build)

Plots radial build cross-section
"""
@recipe function plot_radial_build_cx(rb::IMAS.radial_build)
    seriestype --> :vline
    aspect_ratio --> :equal
    linewidth --> 2
    at = 0
    for l in rb.center_stack
        @series begin
            label --> l.name
            at += l.thickness
            [at]
        end
    end
end