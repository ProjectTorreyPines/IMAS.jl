using IMAS
using Test

do_plot = false
if do_plot
    using Plots
end

@testset "flux_surfaces" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; verbose = false)

    dd_orig = deepcopy(dd)

    IMAS.flux_surfaces(dd.equilibrium)

    ids_orig = dd_orig.equilibrium
    ids = dd.equilibrium

    diff(ids_orig, ids; tol = 1E-1)
end

@testset "eq_calcs" begin
    # Prep sample
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; verbose = false)
    eqt = dd.equilibrium.time_slice[1]
    if !ismissing(eqt.profiles_1d, :r_outboard)
        r_outboard_backup = eqt.profiles_1d.r_outboard
        eqt.profiles_1d.r_outboard .= NaN
    else
        r_outboard_backup = missing
    end
    if !ismissing(eqt.profiles_1d, :r_inboard)
        r_inboard_backup = eqt.profiles_1d.r_inboard
        eqt.profiles_1d.r_inboard .= NaN
    else
        r_inboard_backup = missing
    end

    IMAS.find_midplane_outboard_inboard!(eqt)
    @test !ismissing(eqt.profiles_1d, :r_outboard)
    @test !ismissing(eqt.profiles_1d, :r_inboard)
    @test length(eqt.profiles_1d.r_outboard) == length(eqt.profiles_1d.psi)
    @test length(eqt.profiles_1d.r_inboard) == length(eqt.profiles_1d.psi)
    @test minimum(eqt.profiles_1d.r_outboard) >= eqt.global_quantities.magnetic_axis.r
    @test maximum(eqt.profiles_1d.r_inboard) <= eqt.global_quantities.magnetic_axis.r
end
