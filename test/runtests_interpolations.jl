using IMAS
using Test

# Characterization tests for the interpolation-backed equilibrium/profile workflows
# (flux_surfaces, ρ_interpolant, magnetics!, sol, line_average, avgZ/zeff). Each test
# drives a real interpolant through the workflow that consumes it, so a broken
# interpolation path makes the workflow error and a changed result trips the goldens.
#
# RTOL is the "essentially unchanged" band: interpolant results can drift slightly
# between spline backends, so small interior drift is acceptable — a value moving past
# RTOL is worth inspecting before re-baselining.
const RTOL = 1e-2

@testset "interpolation coverage" begin
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "D3D_eq_ods.json"); show_warnings=false)
    fw = IMAS.first_wall(dd.wall)
    IMAS.flux_surfaces(dd.equilibrium, fw.r, fw.z)
    eqt = dd.equilibrium.time_slice[1]
    gq = eqt.global_quantities
    b = eqt.boundary
    p1 = eqt.profiles_1d
    r0 = gq.magnetic_axis.r
    z0 = gq.magnetic_axis.z

    # ψ_interpolant, find_magnetic_axis, trace_simple_surfaces, interp_rmid_at_psi, Br_Bz, q_95
    @testset "flux_surfaces" begin
        # invariants (backend-independent)
        @test 1.0 < r0 < 2.5 && -0.5 < z0 < 0.5                          # axis inside grid
        @test p1.psi[1] < p1.psi[end]                                    # psi increases axis->edge
        @test abs(gq.q_95) > abs(gq.q_axis) > 1.0                        # q rises outward
        @test gq.ip > 0
        @test b.minor_radius > 0 && b.elongation > 1.0
        @test p1.volume[end] > 0 && p1.area[end] > 0

        # regression goldens
        @test r0                     ≈ 1.7316542632364058  rtol = RTOL
        @test z0                     ≈ -0.02512307340140489 atol = 1e-3
        @test gq.q_axis              ≈ -1.4006582729539356 rtol = RTOL
        @test gq.q_95                ≈ -4.358007499181327  rtol = RTOL
        @test gq.ip                  ≈ 1.0836612803480667e6 rtol = RTOL
        @test p1.psi[1]              ≈ -2.793742510699974  rtol = RTOL
        @test p1.psi[end]            ≈ -1.3391050935759865 rtol = RTOL
        @test b.minor_radius         ≈ 0.5992480020403466  rtol = RTOL
        @test b.elongation           ≈ 1.7824186727819513  rtol = RTOL
        @test b.triangularity_upper  ≈ 0.6026408848657151  rtol = RTOL
        @test b.triangularity_lower  ≈ 0.3272667542197012  rtol = RTOL
        @test b.geometric_axis.r     ≈ 1.6670702530108759  rtol = RTOL
        @test p1.volume[end]         ≈ 19.12590518063392   rtol = RTOL
        @test p1.area[end]           ≈ 1.8817059348806984  rtol = RTOL
    end

    # ρ_interpolant — 2-D cubic spline; consumed off-grid by interferometer/getdata
    @testset "ρ_interpolant" begin
        ri = IMAS.ρ_interpolant(eqt)
        rho_axis = ri.RHO_interpolant(r0, z0)
        rho_01 = ri.RHO_interpolant(r0 + 0.1, z0)
        rho_03 = ri.RHO_interpolant(r0 + 0.3, z0)
        rho_off = ri.RHO_interpolant(2.6, 0.0)   # off-grid (extrapolation)

        # invariants
        @test rho_axis < 0.05                  # ρ ≈ 0 on axis
        @test rho_axis < rho_01 < rho_03       # ρ increases outward
        @test isfinite(rho_off) && rho_off > 1.0   # outside LCFS, extrapolation finite

        # regression goldens (interior)
        @test rho_axis ≈ 0.01975067410323628 rtol = RTOL
        @test rho_01   ≈ 0.15193360172059409 rtol = RTOL
        @test rho_03   ≈ 0.48257505484846797 rtol = RTOL
        # off-grid extrapolation is the value most likely to drift; pinned loosely.
        @test rho_off  ≈ 1.7284309375336697 rtol = 0.15
    end

    # Br_Bz / Bp — poloidal field from ∇ψ
    @testset "Br_Bz / Bp" begin
        _, _, PSI = IMAS.ψ_interpolant(eqt.profiles_2d)
        Br, Bz = IMAS.Br_Bz(PSI, r0 + 0.2, z0)
        Bp = IMAS.Bp(PSI, r0 + 0.2, z0)

        # invariants
        @test Bp ≈ hypot(Br, Bz)                # definition
        @test Bp > 0
        @test PSI(r0, z0) ≈ p1.psi[1] rtol = 1e-6   # interpolant reproduces axis ψ

        # regression goldens
        @test Br ≈ 0.0016167584953939427 atol = 1e-4
        @test Bz ≈ -0.18795527297458944  rtol = RTOL
        @test Bp ≈ 0.18796222638334772   rtol = RTOL
    end

    # magnetics! → field! / flux! on a synthetic probe + loop
    @testset "magnetics! (field! / flux!)" begin
        pr, pz = r0 + 0.3, 0.2
        resize!(dd.magnetics.b_field_pol_probe, 1)
        probe = dd.magnetics.b_field_pol_probe[1]
        probe.position.r = pr
        probe.position.z = pz
        probe.poloidal_angle = 0.0
        resize!(dd.magnetics.flux_loop, 1)
        loop = dd.magnetics.flux_loop[1]
        resize!(loop.position, 1)
        loop.position[1].r = pr
        loop.position[1].z = pz

        IMAS.magnetics!(dd.magnetics, dd.equilibrium)
        _, _, PSI = IMAS.ψ_interpolant(eqt.profiles_2d)

        # invariants: probe (poloidal_angle=0) measures Br; loop measures ψ
        Br, _ = IMAS.Br_Bz(PSI, pr, pz)
        @test probe.field.data[1] ≈ Br rtol = 1e-6
        @test loop.flux.data[1]   ≈ PSI(pr, pz) rtol = 1e-6

        # regression goldens
        @test probe.field.data[1] ≈ 0.09799427799636301  rtol = RTOL
        @test loop.flux.data[1]   ≈ -2.1356278924110654  rtol = RTOL
    end

    # sol → find_levels_from_P / find_levels_from_wall / find_psi_* / _deduplicate_knots!
    @testset "sol (SOL field lines)" begin
        OFL = IMAS.sol(dd)

        # invariants
        @test Set(keys(OFL)) == Set([:hfs, :lfs, :lfs_far])
        @test !isempty(OFL[:lfs])
        n_total = sum(length(OFL[k]) for k in keys(OFL))
        @test 20 <= n_total <= 28
        @test OFL[:lfs][1].psi ≈ p1.psi[end] rtol = 1e-3   # first SOL surface ≈ separatrix

        # regression goldens (counts are structural; a small shift is worth inspecting)
        @test length(OFL[:hfs])     == 2
        @test length(OFL[:lfs])     == 17
        @test length(OFL[:lfs_far]) == 5
        @test OFL[:lfs][1].psi   ≈ -1.3391045383366853 rtol = RTOL
        @test OFL[:lfs][end].psi ≈ -1.3145960750468815 rtol = RTOL
        @test OFL[:lfs][1].Bp[OFL[:lfs][1].midplane_index] ≈ 0.3044471249238311 rtol = RTOL
    end

    # _deduplicate_knots! — SOL knot dedup. Contract (backend-independent): output strictly
    # increasing with minimal (nextfloat) perturbation, so any correct impl passes.
    @testset "_deduplicate_knots!" begin
        # trailing duplicate bumped by one ulp; earlier values untouched
        out = IMAS._deduplicate_knots!([-8.0, 0.0, 20.0, 20.0])
        @test out[1:3] == [-8.0, 0.0, 20.0]
        @test out[4] == nextfloat(20.0)
        @test issorted(out) && allunique(out)

        # realistic SOL case: a flat spot (P[i] == P[i+1]) in a cumulative quantity
        flat = IMAS._deduplicate_knots!([0.0, 1.0, 1.0, 2.0])
        @test flat == [0.0, 1.0, nextfloat(1.0), 2.0]
        @test issorted(flat) && allunique(flat)

        # runs of 3+ duplicates resolve to strictly increasing (not just a no-op)
        out3 = IMAS._deduplicate_knots!([1.0, 1.0, 1.0])
        @test out3[1] == 1.0
        @test issorted(out3) && allunique(out3)

        # already strictly increasing -> returned unchanged
        inc = [0.0, 1.0, 2.5, 9.0]
        @test IMAS._deduplicate_knots!(copy(inc)) == inc

        # mutates in place and returns the same object
        v = [0.0, 0.0]
        @test IMAS._deduplicate_knots!(v) === v
        @test v[2] == nextfloat(0.0)
    end

    # line_average — a real ρ_interpolant consumer (interferometer-style LOS)
    @testset "line_average (ρ_interpolant consumer)" begin
        ri = IMAS.ρ_interpolant(eqt)
        rho = p1.rho_tor_norm
        ne = collect(range(8.0e19, 2.0e19, length(rho)))   # synthetic peaked ne
        r1 = minimum(fw.r) + 0.05
        r2 = maximum(fw.r) - 0.05
        res = IMAS.line_average(ri, eqt.boundary.outline.r, eqt.boundary.outline.z,
            ne, rho, r1, z0, 0.0, r2, z0, 0.0)

        # invariants: chord-averaged ne must lie within the profile range
        @test 2.0e19 <= res.line_average <= 8.0e19
        @test res.path_length > 0

        # regression goldens
        @test res.line_average  ≈ 5.378698461960459e19 rtol = RTOL
        @test res.line_integral ≈ 6.863219237461546e19 rtol = RTOL
    end

    # avgZ / zeff — average ionization-state lookup (2-D gridded linear, log-log space)
    @testset "avgZ / zeff" begin
        # THEORETICAL anchors (backend-independent physics)
        @test IMAS.avgZ(1.0, 1000.0) == 1.0                  # H has no bound electrons
        @test IMAS.avgZ(6.0, 1000.0) ≈ 6.0 rtol = 1e-6       # carbon fully stripped @1keV
        for (Z, Ti) in ((6.0, 50.0), (18.0, 200.0), (18.0, 2000.0), (74.0, 1000.0))
            @test 1.0 <= IMAS.avgZ(Z, Ti) <= Z               # 0 < <Z> <= nuclear charge
        end
        # monotone in Ti (hotter => more ionized)
        zlo = IMAS.avgZ(18.0, 100.0)
        zhi = IMAS.avgZ(18.0, 5000.0)
        @test zlo < zhi

        # regression goldens
        @test IMAS.avgZ(18.0, 2000.0) ≈ 16.646891208683993 rtol = RTOL
        @test IMAS.avgZ(6.0, [100.0, 1000.0, 5000.0]) ≈ [5.280000000000001, 6.0, 6.0] rtol = RTOL

        # zeff over a real core_profiles slice (drives avgZ for every ion)
        dd2 = IMAS.json2imas(joinpath(samp, "omas_sample.json"); show_warnings=false)
        cp1d = dd2.core_profiles.profiles_1d[1]
        zf = IMAS.zeff(cp1d)
        @test all(zf .>= 1.0)                                # Zeff >= 1 always
        @test zf[1]                  ≈ 2.369490401395846 rtol = RTOL
        @test zf[end]                ≈ 3.612027923377595 rtol = RTOL
        @test sum(zf) / length(zf)   ≈ 2.6828619948019363 rtol = RTOL
    end

    # avgZ extrapolation — outside the table Ti range the lookup must stay CONSTANT (flat);
    # pins the extrapolation semantics, not just a number. Table Ti ≈ 0.5 eV .. 100 keV.
    @testset "avgZ extrapolation (flat beyond table)" begin
        # high-T: above the 100 keV table max -> flat -> equals the boundary value
        @test IMAS.avgZ(18.0, 500_000.0) == IMAS.avgZ(18.0, 100_000.0)
        @test IMAS.avgZ(18.0, 100_000.0) ≈ 18.0 rtol = 1e-6     # argon fully stripped
        @test IMAS.avgZ(6.0, 1_000_000.0) ≈ 6.0 rtol = 1e-6     # carbon fully stripped

        # low-T: below the 0.5 eV table min -> flat -> equals the boundary value
        @test IMAS.avgZ(18.0, 0.2) == IMAS.avgZ(18.0, 0.5)
        @test 0.0 <= IMAS.avgZ(18.0, 0.5) < 1.0                 # argon ~neutral when cold
        @test IMAS.avgZ(18.0, 0.5) ≈ 5.349999998571775e-8 rtol = RTOL
    end
end
