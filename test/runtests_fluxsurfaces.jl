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

@testset "fluxsurface_extrema" begin
    # explicit polyline; min_z ties at idx 1 and 5 -> first index wins (matches findmin)
    pr = [1.0, 2.0, 1.5, 0.5, 1.0]
    pz = [0.0, 0.5, 1.0, 0.3, 0.0]
    #     imaxr iminr imaxz iminz  r@maxz max_z r@minz min_z  z@maxr max_r z@minr min_r
    @test IMAS.fluxsurface_extrema(pr, pz) ==
          (2, 4, 3, 1, 1.5, 1.0, 1.0, 0.0, 0.5, 2.0, 0.3, 0.5)

    # equivalence with findmax/findmin on a closed D-shaped polyline
    θ = range(0, 2π; length = 257)
    pr2 = collect(1.7 .+ 0.6 .* cos.(θ .+ 0.3 .* sin.(θ)))
    pz2 = collect(1.8 .* sin.(θ))
    pr2[end], pz2[end] = pr2[1], pz2[1]
    g = IMAS.fluxsurface_extrema(pr2, pz2)
    @test (g[1], g[2], g[3], g[4]) == (findmax(pr2)[2], findmin(pr2)[2], findmax(pz2)[2], findmin(pz2)[2])
    @test (g[6], g[8], g[10], g[12]) == (maximum(pz2), minimum(pz2), maximum(pr2), minimum(pr2))

    @test_throws DimensionMismatch IMAS.fluxsurface_extrema([1.0, 2.0], [0.0])
end