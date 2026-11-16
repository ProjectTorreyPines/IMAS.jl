using IMAS
using Test

include("helpers.jl")

# ---------------------------------------------------------------------------------------------------------------------
# Smoke tests: each function runs on the sample data and returns finite values of the right shape.
# ---------------------------------------------------------------------------------------------------------------------

finite(x) = all(isfinite, x)

@testset "smoke" begin
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "omas_sample.json"); show_warnings=false)
    fw = IMAS.first_wall(dd.wall)
    IMAS.flux_surfaces(dd.equilibrium, fw.r, fw.z)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    n = length(cp1d.grid.rho_tor_norm)

    @testset "equilibrium.jl" begin
        @test isfinite(IMAS.B0_geo(eqt))
        @test isfinite(IMAS.b_field_poloidal_average(eqt))
        @test 1.0 < IMAS.elongation_limit(eqt) < 3.0
        e2 = detached_copy(eqt)
        IMAS.symmetrize_equilibrium!(e2)
        @test finite(findfirst(:rectangular, e2.profiles_2d).psi)
    end

    @testset "currents.jl" begin
        @test isfinite(IMAS.Ip(eqt))
        @test isfinite(IMAS.Ip(cp1d, eqt))
        @test isfinite(IMAS.Ip_bootstrap(cp1d, eqt))
        @test isfinite(IMAS.Ip_ohmic(cp1d, eqt))
        @test isfinite(IMAS.Ip_non_inductive(cp1d, eqt))
        @test isfinite(IMAS.vloop(cp1d, eqt))
        @test isfinite(IMAS.plasma_lumped_resistance(dd))
        # the parallel/toroidal current conversions are inverses of each other
        rho = cp1d.grid.rho_tor_norm
        jtor = 1e6 .* (1.0 .- rho .^ 2)
        for includes_bootstrap in (false, true)
            jpar = IMAS.Jtor_2_Jpar(rho, jtor, includes_bootstrap, eqt)
            @test finite(jpar)
            @test isapprox(IMAS.Jpar_2_Jtor(rho, jpar, includes_bootstrap, eqt), jtor; rtol=1e-8)
            jparB = IMAS.JtoR_2_JparB(rho, jtor, includes_bootstrap, eqt)
            @test finite(jparB)
            @test isapprox(IMAS.JparB_2_JtoR(rho, jparB, includes_bootstrap, eqt), jtor; rtol=1e-8)
        end
    end

    @testset "neoclassical.jl" begin
        cf = IMAS.collision_frequencies(cp1d)
        @test all(finite(getproperty(cf, k)) for k in keys(cf))
        @test finite(IMAS.nuistar(eqt, cp1d))
        @test finite(IMAS.neo_conductivity(eqt, cp1d))
        @test isfinite(IMAS.collisionless_bootstrap_coefficient(eqt, cp1d))
        for neo_2021 in (false, true)
            jb = IMAS.Sauter_neo2021_bootstrap(eqt, cp1d; neo_2021)
            @test length(jb) == n && finite(jb)
        end
    end

    @testset "radiation.jl" begin
        dd2 = deepcopy(dd)
        IMAS.bremsstrahlung_source!(dd2)
        IMAS.line_radiation_source!(dd2)
        IMAS.synchrotron_source!(dd2)
        for index in (8, 9, 10) # bremsstrahlung, synchrotron, line radiation
            srcs = [src for src in dd2.core_sources.source if src.identifier.index == index]
            @test !isempty(srcs)
            for src in srcs
                qe = src.profiles_1d[].electrons.energy
                @test length(qe) == n && finite(qe) && all(qe .<= 0.0)
            end
        end
    end

    @testset "transport.jl" begin
        # (total_fluxes needs core_transport ion elements, which the sample does not carry)
        D = IMAS.calculate_diffusivities(dd)
        @test all(finite(getproperty(D, k)) for k in keys(D))
    end

    @testset "pedestal.jl" begin
        @test isfinite(IMAS.pedestal_poloidal_beta(eqt, 1e4))
    end

    @testset "rf.jl" begin
        res = IMAS.ech_resonance(eqt)
        @test isfinite(res.frequency) && res.frequency > 1e9
        layer = IMAS.ech_resonance_layer(eqt, res.frequency)
        @test !isempty(layer)
    end

    @testset "particles.jl" begin
        psi_norm = collect(range(0.0, 1.0; length=21))
        source_1d = exp.(-psi_norm ./ 0.3)
        parts = IMAS.define_particles(eqt, psi_norm, source_1d, 100)
        @test length(parts.particles) == 100
        @test isfinite(parts.I_per_trace) && isfinite(parts.dr) && isfinite(parts.dz)
    end

    @testset "tf.jl / pf_active.jl" begin
        r = collect(range(1.0, 3.0; length=11))
        rip = IMAS.tf_ripple(r, 4.0, 18)
        @test length(rip) == 11 && finite(rip) && issorted(rip)
        @test IMAS.R_tf_ripple(3.0, rip[end], 18) ≈ 4.0
        for coil in dd.pf_active.coil
            @test IMAS.area(coil) > 0.0
            @test IMAS.volume(coil) > 0.0
        end
    end

    @testset "sol.jl" begin
        # Eich regression #15 (NF 53 093031) with the ITER parameters of its table 1: the paper quotes λq = 0.73 mm
        @test IMAS.widthSOL_eich(6.2, 2.0, 1.185, 100E6) ≈ 0.73E-3 rtol = 0.05

        # the aspect ratio of the scaling is the one of the minor radius
        R0 = (eqt.profiles_1d.r_outboard[end] + eqt.profiles_1d.r_inboard[end]) / 2.0
        a = (eqt.profiles_1d.r_outboard[end] - eqt.profiles_1d.r_inboard[end]) / 2.0
        Bpol_omp = IMAS.Bpol_omp(eqt)
        Psol = 2E6
        @test a ≈ eqt.boundary.minor_radius
        @test IMAS.widthSOL_eich(eqt, Psol) ≈ IMAS.widthSOL_eich(R0, a, Bpol_omp, Psol)

        Psol = IMAS.power_sol(dd.core_sources, cp1d)
        ne_ped = IMAS.interp1d(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal).(0.95)
        @test IMAS.widthSOL_eich(dd) ≈ IMAS.widthSOL_eich(R0, a, Bpol_omp, Psol)
        @test IMAS.widthSOL_sieglin(dd) ≈ IMAS.widthSOL_sieglin(R0, a, Bpol_omp, Psol, ne_ped)
    end

    @testset "thermal_loads.jl" begin
        ddD, fwD = fresh_sample()
        IMAS.flux_surfaces(ddD.equilibrium, fwD.r, fwD.z)
        r = collect(range(0.0, 0.4; length=200))
        q = 1e7 .* exp.(-r ./ 0.01)
        mesh = IMAS.mesher_heat_flux(ddD; r, q)
        @test finite(mesh.Rwall) && finite(mesh.Zwall)
        @test length(mesh.Rwall) == length(mesh.Zwall) > length(fwD.r)
        hf = IMAS.particle_heat_flux(mesh.SOL, mesh.Rwall, mesh.Zwall, mesh.r, mesh.q)
        @test finite(hf[1]) && maximum(hf[1]) > 0.0

        # the mesh is closed, and the abscissa of its closing point includes the closing segment
        @test (mesh.Rwall[1], mesh.Zwall[1]) == (mesh.Rwall[end], mesh.Zwall[end])
        @test mesh.s[end] - mesh.s[end-1] ≈ hypot(mesh.Rwall[end] - mesh.Rwall[end-1], mesh.Zwall[end] - mesh.Zwall[end-1])
        # the outer midplane point (the first one) is in the mesh once
        @test count((mesh.Rwall[1:end-1] .== mesh.Rwall[1]) .& (mesh.Zwall[1:end-1] .== mesh.Zwall[1])) == 1

        # psi levels given as a vector: with the levels of the mesh above it is the same mesh, with and without
        # the wall, and the vector is not changed
        levels = [ofl.psi for ofl in vcat(mesh.SOL[:lfs], mesh.SOL[:lfs_far])]
        levels[1] += 1e-3 * (levels[2] - levels[1])  # the first level is always the separatrix
        levels0 = copy(levels)
        mesh_levels = IMAS.mesher_heat_flux(ddD; r, q, levels)
        @test levels == levels0
        @test mesh_levels.Rwall == mesh.Rwall && mesh_levels.Zwall == mesh.Zwall
        mesh_sol = IMAS.mesher_heat_flux(ddD; r, q, merge_wall=false)
        mesh_sol_levels = IMAS.mesher_heat_flux(ddD; r, q, levels, merge_wall=false)
        @test levels == levels0
        @test mesh_sol_levels.Rwall == mesh_sol.Rwall && mesh_sol_levels.Zwall == mesh_sol.Zwall
        @test length(mesh_sol.Rwall) < length(mesh.Rwall)
    end
end
