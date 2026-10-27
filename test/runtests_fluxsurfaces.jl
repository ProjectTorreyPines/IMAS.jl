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
# Luce superellipse boundary with prescribed per-quadrant squareness
# T.C. Luce, Plasma Phys. Control. Fusion 55 (2013) 095009
function luce_boundary(; R0=1.7, a=0.6, Zrm=0.0, κu=1.8, κl=1.9, δu=0.4, δl=0.5,
    ζ=(ou=0.10, iu=-0.20, il=0.25, ol=-0.05), n=401)
    nexp(ζ) = -log(2) / log(1 / sqrt(2) + ζ * (1 - 1 / sqrt(2)))
    pr = Float64[]
    pz = Float64[]
    for θ in range(0, 2π; length=n)[1:end-1]
        if θ < π / 2          # upper outer
            t = θ; A = a * (1 + δu); B = κu * a; nn = nexp(ζ.ou); Rc = R0 - a * δu; sx = +1; sz = +1
        elseif θ < π          # upper inner
            t = π - θ; A = a * (1 - δu); B = κu * a; nn = nexp(ζ.iu); Rc = R0 - a * δu; sx = -1; sz = +1
        elseif θ < 3π / 2     # lower inner
            t = θ - π; A = a * (1 - δl); B = κl * a; nn = nexp(ζ.il); Rc = R0 - a * δl; sx = -1; sz = -1
        else                  # lower outer
            t = 2π - θ; A = a * (1 + δl); B = κl * a; nn = nexp(ζ.ol); Rc = R0 - a * δl; sx = +1; sz = -1
        end
        x = A * abs(cos(t))^(2 / nn)
        y = B * abs(sin(t))^(2 / nn)
        push!(pr, Rc + sx * x)
        push!(pz, Zrm + sz * y)
    end
    return pr, pz
end

@testset "luce_squareness quadrant labels" begin
    # four DISTINCT squareness values so that a label swap between quadrants is detected
    ζ = (ou=0.10, iu=-0.20, il=0.25, ol=-0.05)
    pr, pz = luce_boundary(; ζ)
    g = IMAS.miller_R_a_κ_δ_ζ(pr, pz)
    @test isapprox(g.ζou, ζ.ou; atol=2e-3)
    @test isapprox(g.ζol, ζ.ol; atol=2e-3)
    @test isapprox(g.ζil, ζ.il; atol=2e-3)
    @test isapprox(g.ζiu, ζ.iu; atol=2e-3)
    @test isapprox(g.δu, 0.4; atol=1e-6)
    @test isapprox(g.δl, 0.5; atol=1e-6)
end
