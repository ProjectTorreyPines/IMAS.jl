using IMAS
using Test

# Characterization tests for the 1-D interpolation wrappers cubic_interp1d /
# linear_interp1d (math.jl), currently DataInterpolations-backed and about to be
# migrated to FastInterpolations. Each test drives a wrapper through a real
# workflow (currents, interp_rmid_at_psi, neoclassical) plus a direct
# extrapolation-semantics check, so a broken/changed interpolation surfaces here.
#
# RTOL is the "essentially unchanged" band. Note: linear interpolation is the same
# math in both backends, so the neoclassical goldens should stay tight; the cubic
# spline construction differs (boundary conditions), so cubic goldens may drift
# slightly — inspect before re-baselining if a value moves past RTOL.
const RTOL = 1e-2

@testset "interp1d coverage" begin
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "D3D_eq_ods.json"); show_warnings=false)
    fw = IMAS.first_wall(dd.wall)
    IMAS.flux_surfaces(dd.equilibrium, fw.r, fw.z)
    eqt = dd.equilibrium.time_slice[1]

    # cubic_interp1d — currents: gm-profiles cubic-interpolated onto a 51-pt grid
    # (off the ~native eqt grid), then the Jtor↔Jpar transform.
    @testset "cubic_interp1d: currents" begin
        rho_out = collect(range(0.0, 1.0, 51))
        Jtor = 1.0e6 .* (1.0 .- rho_out .^ 2)
        Jpar = IMAS.Jtor_2_Jpar(rho_out, Jtor, false, eqt)

        # invariant: round-trip Jpar→Jtor recovers the input (backend-independent)
        Jtor_rt = IMAS.Jpar_2_Jtor(rho_out, Jpar, false, eqt)
        @test maximum(abs.(Jtor_rt .- Jtor)) < 10.0     # ~1e-10 in practice

        # regression goldens
        @test Jpar[1]  ≈ 998944.9011757972 rtol = RTOL
        @test Jpar[26] ≈ 777436.8712898433 rtol = RTOL
        @test abs(Jpar[end]) < 1.0                      # Jtor[end]=0 ⇒ Jpar[end]≈0
    end

    # cubic_interp1d — interp_rmid_at_psi (fluxsurfaces): R(ψ) on the outboard midplane
    @testset "cubic_interp1d: interp_rmid_at_psi" begin
        _, _, PSI = IMAS.ψ_interpolant(eqt.profiles_2d)
        ra, za = eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z
        R = collect(range(ra, maximum(fw.r) - 0.05, 40))
        rmid = IMAS.interp_rmid_at_psi(PSI, R, za)

        # invariant: at ψ_axis the outboard midplane radius is the magnetic axis
        @test rmid(eqt.profiles_1d.psi[1]) ≈ ra rtol = 1e-3

        # regression goldens
        @test rmid(eqt.profiles_1d.psi[1])   ≈ 1.7316542632361691 rtol = RTOL
        @test rmid(eqt.profiles_1d.psi[end]) ≈ 2.2665269416113203 rtol = RTOL
    end

    # linear_interp1d — neoclassical: eqt geometry (R, a, q, ...) mapped onto the
    # cp1d grid, which has a different length ⇒ a genuine off-node linear interp.
    @testset "linear_interp1d: neoclassical" begin
        dd2 = IMAS.json2imas(joinpath(samp, "omas_sample.json"); show_warnings=false)
        fw2 = IMAS.first_wall(dd2.wall)
        IMAS.flux_surfaces(dd2.equilibrium, fw2.r, fw2.z)
        eqt2 = dd2.equilibrium.time_slice[1]
        cp1d = dd2.core_profiles.profiles_1d[1]
        @test length(cp1d.grid.rho_tor_norm) != length(eqt2.profiles_1d.rho_tor_norm)

        nue = IMAS.nuestar(eqt2, cp1d)
        @test all(nue .> 0)                              # collisionality > 0
        @test nue[1]                 ≈ 97.56723198386155 rtol = RTOL
        @test nue[end]               ≈ 8.07637773177939  rtol = RTOL
        @test sum(nue) / length(nue) ≈ 9.697133802559067 rtol = RTOL

        cond = IMAS.neo_conductivity(eqt2, cp1d)
        @test all(cond .> 0)
        @test cond[1]   ≈ 1.7799469179575473e8 rtol = RTOL
        @test cond[end] ≈ 788970.6202650126    rtol = RTOL
    end

    # extrapolation semantics — DataInterpolations `Extension` continues the edge
    # polynomial. The FI replacement must keep linear continuation outside the grid
    # or the linear checks break (cf. the avgZ Flat() check in runtests_interpolations).
    @testset "interp1d extrapolation (Extension)" begin
        litp = IMAS.linear_interp1d([0.0, 1.0, 2.0], [0.0, 10.0, 20.0])
        @test litp(1.5) ≈ 15.0 rtol = 1e-9               # exact on a line
        @test litp(3.0) ≈ 30.0 rtol = 1e-9               # Extension = linear continuation

        citp = IMAS.cubic_interp1d([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 8.0, 27.0])  # = x^3
        @test citp(1.0) ≈ 1.0                            # interpolates through the nodes
        @test citp(2.0) ≈ 8.0
        # A cubic spline should reproduce a cubic exactly. DataInterpolations' default
        # (natural) BC does NOT for few points — citp(1.5)=3.15 — while FastInterpolations'
        # CubicFit BC does (3.375 exactly). So this is broken on the DataInterp baseline and
        # will start passing once cubic_interp1d moves to FI: flip @test_broken → @test then.
        @test_broken citp(1.5) ≈ 1.5^3 rtol = 1e-6
        @test citp(3.5) > 27.0                           # extrapolation continues upward
    end
end
