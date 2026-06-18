using IMAS
using Test

# Characterization tests that PIN the current behavior of field-line tracing
# (src/physics/fields.jl) before the implicit-midpoint root solver is swapped
# from NLsolve to SimpleNonlinearSolve. Two cases match closed-form theory
# (solver-independent); one records the trajectory on the bundled D3D
# equilibrium as a regression snapshot.

_vecnorm(v) = sqrt(sum(abs2, v))

@testset "field line tracing" begin

    # --- Case A: constant unit field => exact straight line --------------------
    # obj(next) = next - current - h*vf(midpoint); with vf ≡ const the implicit
    # midpoint solution is exactly next = current + h*dir, every step.
    @testset "constant field is a straight line" begin
        dir = [1.0, 0.0, 0.0]
        x0 = [1.0, 0.5, 0.0]
        h = 0.1
        # NOTE: trace_field_line mutates its start_point in place, so pass a copy
        # and keep x0 as the untouched comparison baseline.
        traj = IMAS.trace_field_line(p -> dir, copy(x0), h, obj -> obj.count >= 5)
        @test length(traj) == 6
        for k in eachindex(traj)
            @test traj[k] ≈ x0 .+ (k - 1) * h .* dir atol = 1e-7
        end
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

    # unknown method is rejected
    @test_throws ErrorException IMAS.trace_field_line(
        p -> [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.1, obj -> obj.count >= 1; method=:bogus)
end
