using IMAS
using Test

do_plot = false
if do_plot
    using Plots
end

@testset "flux_surfaces" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; show_warnings = false)

    dd_orig = deepcopy(dd)

    fw = IMAS.first_wall(dd.wall)
    IMAS.flux_surfaces(dd.equilibrium, fw.r, fw.z)

    ids_orig = dd_orig.equilibrium
    ids = dd.equilibrium

    diff(ids_orig, ids; tol = 1E-1)
end