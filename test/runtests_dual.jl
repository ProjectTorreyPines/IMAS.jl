using IMAS
using Test
import IMAS.ForwardDiff
import IMAS.ForwardDiff: value, partials

# Derivatives carried by dual numbers are compared with central finite differences of the Float64 code

const D2 = IMAS.InputsDual{2}

central(f, h) = (f(+h) .- f(-h)) ./ (2h)

@testset "dual-typed dd: conversion, seeding, value and partials" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd = IMAS.json2imas(filename; show_warnings=false)

    ddd = IMAS.retype(D2, dd)
    @test eltype(ddd) === D2
    @test ddd.global_time == dd.global_time
    @test IMAS.npartials(ddd) == 2
    @test convert(D2, dd.equilibrium) isa IMAS.equilibrium{D2} # no ambiguity with ForwardDiff's convert

    eqt = ddd.equilibrium.time_slice[]
    f0 = dd.equilibrium.time_slice[].profiles_1d.f
    eqt.profiles_1d.f = eqt.profiles_1d.f .* IMAS.seed(D2, 2)
    @test value.(eqt.profiles_1d.f) == f0
    @test partials.(eqt.profiles_1d.f, 2) == f0
    @test all(iszero, partials.(eqt.profiles_1d.f, 1))
    @test eltype(eqt.profiles_1d.psi) === D2
    @test eqt.time isa Float64 # times stay plain floats

    back = IMAS.value_dd(ddd)
    @test eltype(back) === Float64
    @test back.equilibrium.time_slice[].profiles_1d.f == f0
    @test back.equilibrium.time_slice[].profiles_2d[1].psi == dd.equilibrium.time_slice[].profiles_2d[1].psi

    d2 = IMAS.partials_dd(ddd, 2)
    @test d2.equilibrium.time_slice[].profiles_1d.f == f0
    @test all(iszero, d2.equilibrium.time_slice[].profiles_1d.psi)

    # values are replaced, partials are kept
    target = deepcopy(dd)
    target.equilibrium.time_slice[].profiles_1d.f = 2 .* f0
    @test IMAS.replace_values!(ddd, target) == 0
    @test value.(eqt.profiles_1d.f) == 2 .* f0
    @test partials.(eqt.profiles_1d.f, 2) == f0

    # the tag of the inputs nests outside any other ForwardDiff tag
    g = ForwardDiff.derivative(x -> x^2 * IMAS.seed(D2, 1), 3.0)
    @test value(g) == 6.0 && partials(g, 1) == 6.0

    # files hold plain numbers: one for the values and one per partial
    dir = mktempdir()
    IMAS.imas2json(ddd.equilibrium, joinpath(dir, "eq.json"))
    @test isfile(joinpath(dir, "eq.json")) && isfile(joinpath(dir, "eq_d1.json")) && isfile(joinpath(dir, "eq_d2.json"))
end

@testset "derivatives of iterative searches" begin
    a = 2.0 * IMAS.seed(D2, 1)
    b = 3.0 * IMAS.seed(D2, 2)

    # root of G(u) = [u1^2 - a, u1*u2 - b]: u1 = sqrt(a), u2 = b/sqrt(a)
    G = u -> [u[1]^2 - a, u[1] * u[2] - b]
    u = IMAS.implicit_solution(G, [sqrt(2.0), 3.0 / sqrt(2.0)], D2)
    @test partials(u[1], 1) ≈ 0.5 * sqrt(2.0) rtol = 1e-6 # d/dln a of sqrt(a)
    @test partials(u[2], 1) ≈ -0.5 * 3.0 / sqrt(2.0) rtol = 1e-6
    @test partials(u[2], 2) ≈ 3.0 / sqrt(2.0) rtol = 1e-6

    # a condition that is not met leaves its unknown fixed
    G2 = u -> [u[1]^2 - a, 1.0 + 0.0 * u[2] + 0.0 * b]
    u = IMAS.implicit_solution(G2, [sqrt(2.0), 7.0], D2; atol=1e-2)
    @test partials(u[1], 1) ≈ 0.5 * sqrt(2.0) rtol = 1e-6
    @test all(iszero, partials(u[2]))

    # minimum of (x1 - a)^2 + (x2 - a*b)^2 + x1*x2*0
    cost = x -> (x[1] - a)^2 + (x[2] - a * b)^2
    x = IMAS.implicit_argmin(cost, [2.0, 6.0], D2)
    @test partials(x[1], 1) ≈ 2.0
    @test partials(x[1], 2) ≈ 0.0 atol = 1e-12
    @test partials(x[2], 1) ≈ 6.0
    @test partials(x[2], 2) ≈ 6.0

    # switched off: values only
    x = Base.ScopedValues.with(IMAS.IMPLICIT_DERIVATIVES => false) do
        IMAS.implicit_argmin(cost, [2.0, 6.0], D2)
    end
    @test all(iszero, partials(x[1])) && value(x[1]) == 2.0

    # finite derivative parts where the plain function has an infinite derivative
    @test all(isfinite, partials(IMAS.sqrt0(0.0 * IMAS.seed(D2, 1))))
    @test all(isfinite, partials(IMAS._chebyshev_T(3, D2(1.0))))
end

@testset "pedestal blending with duals" begin
    rho = collect(range(0.0, 1.0, 101))
    profile(s) = s .* (200.0 .+ 2800.0 .* (1.0 .- rho .^ 2) .^ 1.5)
    tb, ts = IMAS.pedestal_search_tolerance.brent[], IMAS.pedestal_search_tolerance.simplex[]
    IMAS.pedestal_search_tolerance.brent[] = 1e-12
    IMAS.pedestal_search_tolerance.simplex[] = 1e-14
    try
        # H-mode blend: derivatives with respect to the pedestal height and to a scale of the profile
        blend(ped, s) = IMAS.blend_core_edge_Hmode(profile(s), rho .* one(ped), ped, 0.06 * one(ped), 0.775, 0.8)
        ad = blend(600.0 + (IMAS.seed(D2, 1) - 1), IMAS.seed(D2, 2))
        @test value.(ad) ≈ blend(600.0, 1.0) rtol = 1e-10
        fd1 = central(h -> blend(600.0 + h, 1.0), 0.1) # a step well above the noise of the simplex search
        fd2 = central(h -> blend(600.0, 1.0 + h), 1e-4)
        @test partials.(ad, 1) ≈ fd1 rtol = 2e-2
        @test partials.(ad, 2) ≈ fd2 rtol = 2e-2

        # L-mode blend
        lblend(ped, s) = IMAS.blend_core_edge_Lmode(profile(s), rho .* one(ped), ped, 0.9)
        ad = lblend(300.0 + (IMAS.seed(D2, 1) - 1), IMAS.seed(D2, 2))
        @test value.(ad) ≈ lblend(300.0, 1.0) rtol = 1e-10
        @test partials.(ad, 1) ≈ central(h -> lblend(300.0 + h, 1.0), 1e-3) rtol = 2e-2 atol = 1e-6
        @test partials.(ad, 2) ≈ central(h -> lblend(300.0, 1.0 + h), 1e-5) rtol = 2e-2 atol = 1e-3
    finally
        IMAS.pedestal_search_tolerance.brent[] = tb
        IMAS.pedestal_search_tolerance.simplex[] = ts
    end
end

@testset "flux surfaces with duals" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
    dd0 = IMAS.json2imas(filename; show_warnings=false)
    fw = IMAS.first_wall(dd0.wall)

    # perturbations: (1) toroidal field scale, (2) a change of shape of ψ(R,Z)
    function perturb!(dd, s)
        eqt = dd.equilibrium.time_slice[]
        e1 = eqt.profiles_1d
        e2 = eqt.profiles_2d[1]
        R = Float64[IMAS.force_float64(x) for x in e2.grid.dim1]
        Z = Float64[IMAS.force_float64(x) for x in e2.grid.dim2]
        bump = [((r - 1.7) / 0.6)^2 - 0.5 * (z / 0.6)^2 for r in R, z in Z]
        dpsi = IMAS.force_float64(e1.psi[end] - e1.psi[1])
        e2.psi = e2.psi .+ (s[2] - 1) .* dpsi .* bump
        e1.f_df_dpsi = e1.f_df_dpsi .* s[1]^2
        e1.f = e1.f .* s[1]
        eqt.global_quantities.vacuum_toroidal_field.b0 = eqt.global_quantities.vacuum_toroidal_field.b0 * s[1]
        return dd
    end
    function outputs(dd)
        eqt = dd.equilibrium.time_slice[]
        e1 = eqt.profiles_1d
        n = length(e1.psi)
        idx = [5, n ÷ 3, n ÷ 2, 3n ÷ 4, n - 3]
        g = eqt.global_quantities
        return vcat(e1.q[idx], e1.volume[idx], e1.gm1[idx], e1.gm9[idx], e1.elongation[idx], e1.r_outboard[idx], e1.rho_tor[idx], e1.psi[idx],
            [g.magnetic_axis.r, g.magnetic_axis.z, g.ip, g.beta_normal, eqt.boundary.psi])
    end
    function run(s)
        dd = s isa Vector{Float64} ? deepcopy(dd0) : IMAS.retype(D2, dd0)
        perturb!(dd, s)
        w = IMAS.first_wall(dd.wall)
        IMAS.flux_surfaces(dd.equilibrium.time_slice[], w.r, w.z)
        return outputs(dd)
    end

    y0 = run([1.0, 1.0])
    yd = run([IMAS.seed(D2, 1), IMAS.seed(D2, 2)])
    @test value.(yd) ≈ y0 rtol = 1e-6

    # The last closed surface is found by a search on the contouring grid, which makes every output a
    # piecewise function of ψ(R,Z) with small jumps. The finite-difference slope is therefore a fit over
    # several steps that span many pieces, not a two-point difference inside one of them.
    function slope(k, H)
        hs = H .* vcat(-4:-1, 1:4) ./ 4
        e = zeros(2)
        e[k] = 1.0
        Y = [run([1.0, 1.0] .+ h .* e) for h in hs]
        return [sum(hs[j] * (Y[j][i] - y0[i]) for j in eachindex(hs)) / sum(abs2, hs) for i in eachindex(y0)]
    end
    for (k, H) in ((1, 1e-4), (2, 2e-2))
        fd = slope(k, H)
        ad = partials.(yd, k)
        rel = abs.(ad .- fd) ./ max.(abs.(ad), abs.(fd), abs.(y0)) # error per unit change of the input, relative to the output
        @test count(<(0.02), rel) >= 0.95 * length(rel)
        @test sort(rel)[end÷2] < 5e-3
        @test maximum(rel) < 0.1
    end
end
