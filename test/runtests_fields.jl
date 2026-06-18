using IMAS
using Test

# Characterization tests that PIN the current behavior of field-line tracing
# (src/physics/fields.jl) before the implicit-midpoint root solver is swapped
# from NLsolve to SimpleNonlinearSolve. Two cases match closed-form theory
# (solver-independent); one records the trajectory on the bundled D3D
# equilibrium as a regression snapshot.

_vecnorm(v) = sqrt(sum(abs2, v))

# Unit rotation field has an exact solution: speed 1 ⇒ after arc length s=(k-1)*h
# the k-th point lies on the circle at angle s/r0. These measure TRUE global
# error (vs the exact circle) and the drift of the conserved radius invariant.
_circle_pos_error(traj, r0, h) = maximum(
    hypot(traj[k][1] - r0 * cos((k - 1) * h / r0), traj[k][2] - r0 * sin((k - 1) * h / r0))
    for k in eachindex(traj))
_radius_drift(traj, r0) = maximum(abs(hypot(p[1], p[2]) - r0) for p in traj)

@testset "field line tracing" begin

    # --- Case A: constant unit field => exact straight line --------------------
    # obj(next) = next - current - h*vf(midpoint); with vf ≡ const the implicit
    # midpoint solution is exactly next = current + h*dir, every step.
    @testset "constant field is a straight line" begin
        dir = [1.0, 0.0, 0.0]
        x0 = [1.0, 0.5, 0.0]
        h = 0.1
        # Pass x0 directly (no copy): neither integrator may mutate the caller's
        # start_point, so x0 doubles as the untouched comparison baseline below.
        traj = IMAS.trace_field_line(p -> dir, x0, h, obj -> obj.count >= 5)
        @test length(traj) == 6
        for k in eachindex(traj)
            @test traj[k] ≈ x0 .+ (k - 1) * h .* dir atol = 1e-7
        end
        @test x0 == [1.0, 0.5, 0.0]   # start_point must not be mutated in place
    end

    # --- Case B: rotation field => circle, radius is a conserved invariant -----
    # The implicit midpoint rule exactly conserves quadratic invariants, so the
    # radius ‖(x,y)‖ is preserved to solver tolerance regardless of the solver.
    @testset "rotation field preserves radius (circle)" begin
        vf = p -> (v = [-p[2], p[1], 0.0]; v ./ _vecnorm(v))
        x0 = [1.0, 0.0, 0.0]
        h = 0.05
        traj = IMAS.trace_field_line(vf, copy(x0), h, obj -> obj.count >= 20)
        @test length(traj) == 21
        radii = [hypot(p[1], p[2]) for p in traj]
        @test maximum(abs.(radii .- 1.0)) < 1e-6          # radius conserved
        @test all(p -> abs(p[3]) < 1e-12, traj)           # stays in z = 0 plane
        steps = [_vecnorm(traj[k + 1] .- traj[k]) for k in 1:length(traj) - 1]
        @test all(s -> isapprox(s, h; atol = 1e-6), steps)  # unit field => |step| = h
    end

    # --- Case C: real D3D equilibrium => recorded regression snapshot ----------
    # Golden values recorded from the current (NLsolve trust-region) code.
    # Per-step points are pinned at fixed indices; the total step count depends
    # on a `Δϕ ≥ 2π` stop test, so it is bounded rather than fixed exactly.
    @testset "D3D equilibrium field line (regression snapshot)" begin
        filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
        dd = IMAS.json2imas(filename; show_warnings = false)
        eqt = dd.equilibrium.time_slice[1]
        ar = eqt.global_quantities.magnetic_axis.r
        az = eqt.global_quantities.magnetic_axis.z
        @test ar ≈ 1.73165427 atol = 1e-6
        @test az ≈ -0.025123087 atol = 1e-6

        res = IMAS.trace_field_line(eqt, ar + 0.10, az; step_size = 0.05, max_turns = 1)
        n = length(res.x)

        # deterministic start point
        @test res.x[1] ≈ ar + 0.10 atol = 1e-8
        @test res.y[1] ≈ 0.0 atol = 1e-12
        @test res.z[1] ≈ az atol = 1e-8

        # ~one toroidal turn near the axis (allow ±2 from the Δϕ stop test)
        @test abs(n - 216) <= 2

        # per-step behavior pinned at fixed indices (insensitive to ±1 step count)
        golden = Dict(
            25  => (1.4419151476599015, -1.1086311551344152, -0.09653635125986251),
            50  => (0.39313192673181796, -1.7375050652548147, -0.15113318361646508),
            100 => (-1.579889165027515, -0.577249019215715, -0.15484682301378697),
            150 => (-0.6445336434296377, 1.4960557330994009, -0.03813688860109741),
        )
        for (k, (gx, gy, gz)) in golden
            @test res.x[k] ≈ gx atol = 1e-4
            @test res.y[k] ≈ gy atol = 1e-4
            @test res.z[k] ≈ gz atol = 1e-4
        end

        # coarse geometric extent of the trajectory
        @test minimum(res.r) ≈ 1.6284816238060678 atol = 1e-3
        @test maximum(res.r) ≈ 1.8316579372123754 atol = 1e-3
        @test minimum(res.z) ≈ -0.171835522057876 atol = 1e-3
        @test maximum(res.z) ≈ 0.11310876800067095 atol = 1e-3
    end

    # --- RK4 path (method=:rk4): explicit 4th-order, no nonlinear solve --------
    @testset "RK4 method" begin
        # constant field is still an exact straight line
        dir = [1.0, 0.0, 0.0]
        x0 = [1.0, 0.5, 0.0]
        h = 0.1
        traj = IMAS.trace_field_line(p -> dir, copy(x0), h, obj -> obj.count >= 5; method=:rk4)
        @test length(traj) == 6
        for k in eachindex(traj)
            @test traj[k] ≈ x0 .+ (k - 1) * h .* dir atol = 1e-7
        end

        # rotation field: radius preserved well (RK4 ~1e-9, vs implicit's machine ε)
        vf = p -> (v = [-p[2], p[1], 0.0]; v ./ _vecnorm(v))
        tc = IMAS.trace_field_line(vf, [1.0, 0.0, 0.0], 0.05, obj -> obj.count >= 20; method=:rk4)
        radii = [hypot(p[1], p[2]) for p in tc]
        @test maximum(abs.(radii .- 1.0)) < 1e-6
        @test all(p -> abs(p[3]) < 1e-12, tc)
        steplen = [_vecnorm(tc[k+1] .- tc[k]) for k in 1:length(tc)-1]
        @test all(s -> abs(s - 0.05) < 1e-4, steplen)   # RK4 chord ≈ h (not exact like implicit)

        # D3D snapshot with its own RK4 golden (differs from implicit by ~2e-4)
        filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "D3D_eq_ods.json")
        dd = IMAS.json2imas(filename; show_warnings = false)
        eqt = dd.equilibrium.time_slice[1]
        ar = eqt.global_quantities.magnetic_axis.r
        az = eqt.global_quantities.magnetic_axis.z
        res = IMAS.trace_field_line(eqt, ar + 0.10, az; step_size = 0.05, max_turns = 1, method = :rk4)
        @test res.x[1] ≈ ar + 0.10 atol = 1e-8
        @test abs(length(res.x) - 216) <= 2
        rk4_golden = Dict(
            25 => (1.4419207446530165, -1.1085824850895571, -0.09666252619176363),
            100 => (-1.579674064946065, -0.577035443181165, -0.155075593017569),
        )
        for (k, (gx, gy, gz)) in rk4_golden
            @test res.x[k] ≈ gx atol = 1e-4
            @test res.y[k] ≈ gy atol = 1e-4
            @test res.z[k] ≈ gz atol = 1e-4
        end
    end

    # --- smoke: both methods integrate the same ODE, so they must agree --------
    @testset "implicit vs RK4 agree (smoke)" begin
        vf = p -> (v = [-p[2], p[1], 0.0]; v ./ _vecnorm(v))   # unit rotation → circle
        stop = o -> o.count >= 10
        imp = IMAS.trace_field_line(vf, [1.0, 0.0, 0.0], 0.05, stop; method=:implicit_midpoint)
        rk4 = IMAS.trace_field_line(vf, [1.0, 0.0, 0.0], 0.05, stop; method=:rk4)
        @test length(imp) == length(rk4)
        @test maximum(_vecnorm(imp[k] .- rk4[k]) for k in eachindex(imp)) < 1e-3
    end

    # --- accuracy & invariants vs the EXACT circle ----------------------------
    # Measures true global error (not self-consistency) to pin the method
    # trade-off: RK4 is far more accurate in POSITION (4th- vs 2nd-order), while
    # implicit-midpoint (symplectic) conserves the RADIUS to ~machine precision.
    # A coefficient bug that demoted RK4 below 4th order would break the
    # position-accuracy and convergence-order assertions here.
    @testset "accuracy & invariants (analytic circle)" begin
        r0 = 1.0
        vf = p -> (v = [-p[2], p[1], 0.0]; v ./ _vecnorm(v))   # unit rotation → circle
        circle_traj(h, nsteps; method) =
            IMAS.trace_field_line(vf, [r0, 0.0, 0.0], h, o -> o.count >= nsteps; method)

        # five turns at h = 0.05
        h = 0.05
        nsteps = round(Int, 5 * 2pi * r0 / h)
        imp = circle_traj(h, nsteps; method=:implicit_midpoint)
        rk4 = circle_traj(h, nsteps; method=:rk4)
        imp_pos, rk4_pos = _circle_pos_error(imp, r0, h), _circle_pos_error(rk4, r0, h)
        imp_rad, rk4_rad = _radius_drift(imp, r0), _radius_drift(rk4, r0)

        # RK4 is the more accurate integrator (position vs the true circle)
        @test rk4_pos < 1e-5
        @test 100 * rk4_pos < imp_pos          # RK4 ≳ 100x better in position

        # implicit-midpoint is the better invariant-preserver (radius)
        @test imp_rad < rk4_rad                 # symplectic conserves radius better
        @test imp_rad < 1e-9                    # ...to ~machine precision

        # convergence order: halving h ⇒ RK4 ~4th-order, implicit ~2nd-order
        n1, n2 = round(Int, 2pi * r0 / 0.05), round(Int, 2pi * r0 / 0.025)
        rk4_e1 = _circle_pos_error(circle_traj(0.05, n1; method=:rk4), r0, 0.05)
        rk4_e2 = _circle_pos_error(circle_traj(0.025, n2; method=:rk4), r0, 0.025)
        imp_e1 = _circle_pos_error(circle_traj(0.05, n1; method=:implicit_midpoint), r0, 0.05)
        imp_e2 = _circle_pos_error(circle_traj(0.025, n2; method=:implicit_midpoint), r0, 0.025)
        @test log2(rk4_e1 / rk4_e2) > 3.0        # ~4th order (fails if RK4 demoted)
        @test 1.7 < log2(imp_e1 / imp_e2) < 2.3  # ~2nd order
    end

    # unknown method is rejected
    @test_throws ErrorException IMAS.trace_field_line(
        p -> [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, obj -> obj.count >= 1; method=:bogus)
end
