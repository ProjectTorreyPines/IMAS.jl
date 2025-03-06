using IMAS
using Test
using IMAS.Plots
using HelpPlots

@testset "plot_recipes" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; show_warnings = false)

    plot(dd.equilibrium)

    help_plot(dd.equilibrium)
    help_plot!(dd.equilibrium)
end