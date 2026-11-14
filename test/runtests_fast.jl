using IMAS
using Test

# Guards the ONLY Roots.find_zero call site in the package: critical_energy(...;
# approximate=false) in src/physics/fast.jl. Nothing else exercises Roots, so
# this test pins the result and keeps the widened `Roots = "2, 3"` compat honest
# (the value below was verified to match under both Roots 2.x and 3.x).
@testset "fast particles" begin
    filename = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample", "omas_sample.json")
    dd = IMAS.json2imas(filename; show_warnings=false)
    cp1d = dd.core_profiles.profiles_1d[1]

    @testset "critical_energy Roots path (approximate=false)" begin
        ec_approx = IMAS.critical_energy(cp1d, 1, 2.0, 1.0; approximate=true)
        ec_roots = IMAS.critical_energy(cp1d, 1, 2.0, 1.0; approximate=false)  # Roots.find_zero

        # the root-find must run and return a sensible, distinct value
        @test isfinite(ec_roots) && ec_roots > 0
        @test !isapprox(ec_roots, ec_approx; rtol=1e-3)   # Roots path actually changed the answer

        # regression golden (verified to match under both Roots 2.x and 3.x)
        @test ec_roots ≈ 93662.596 rtol = 1e-5
    end
end

@testset "fast_density_cap!" begin
    L = [1.0, 1.0, 1.0, 1.0, 0.0]
    nf = [0.5, 0.9, 1.5, 10.0, 0.3]
    @test IMAS.fast_density_cap!(copy(nf), L, 1.0) == min.(nf, L)
    soft = IMAS.fast_density_cap!(copy(nf), L, 0.7)
    @test soft[1] == 0.5
    @test all(0.7 .< soft[2:4] .<= 1.0) && issorted(soft[2:4])
    @test soft[5] == 0.0
    # no slope break at the knee
    y = IMAS.fast_density_cap!(collect(range(0.0, 3.0, 3001)), ones(3001), 0.7)
    @test maximum(abs, diff(diff(y))) < 1e-5
end

@testset "fast_particles_profiles! called repeatedly" begin
    # a deuterium + carbon plasma with one deuterium fast-ion source of amplitude S0
    function plasma_with_fast_source(S0)
        dd = IMAS.dd()
        rho = collect(range(0.0, 1.0, 21))
        n = 5e19 .* (1 .- 0.8 .* rho .^ 2)
        T = 3e3 .* (1 .- 0.9 .* rho .^ 2)
        cp1d = resize!(dd.core_profiles.profiles_1d, 1)[1]
        cp1d.time = 0.0
        dd.core_profiles.time = [0.0]
        dd.global_time = 0.0
        cp1d.grid.rho_tor_norm = rho
        cp1d.electrons.density_thermal = copy(n)
        cp1d.electrons.temperature = copy(T)
        resize!(cp1d.ion, 2)
        for (k, (label, z, a, fraction)) in enumerate((("D", 1.0, 2.0, 0.88), ("C", 6.0, 12.0, 0.02)))
            ion = cp1d.ion[k]
            ion.label = label
            resize!(ion.element, 1)
            ion.element[1].z_n = z
            ion.element[1].a = a
            ion.density_thermal = fraction .* n
            ion.temperature = copy(T)
        end
        source1d = resize!(resize!(dd.core_sources.source, 1)[1].profiles_1d, 1)[1]
        source1d.time = 0.0
        dd.core_sources.time = [0.0]
        source1d.grid.rho_tor_norm = rho
        sion = resize!(source1d.ion, 1)[1]
        sion.label = "D"
        resize!(sion.element, 1)
        sion.element[1].z_n = 1.0
        sion.element[1].a = 2.0
        sion.particles = S0 .* (1 .- rho .^ 2) .+ 1e15
        sion.fast_particles_energy = 80e3
        return dd
    end

    # S0 = 2e19: fast ions well below the cap; S0 = 2e21: the cap sets the fast-ion density
    for S0 in (2e19, 2e21)
        # electron density from quasi-neutrality: the thermal ions are left alone, however often this is called
        dd = plasma_with_fast_source(S0)
        cp1d = dd.core_profiles.profiles_1d[1]
        thermal0 = [copy(ion.density_thermal) for ion in cp1d.ion]
        fast = Vector{Float64}[]
        for _ in 1:8
            IMAS.fast_particles_profiles!(dd.core_sources, cp1d; modify_electron_density=true)
            D = cp1d.ion[1]
            @test all(ion.density_thermal == thermal0[k] for (k, ion) in enumerate(cp1d.ion))
            @test all(D.density_fast .<= 0.8 .* thermal0[1] .* (1 + 1e-12))
            # quasi-neutrality with the fast ions (carbon is not quite fully stripped at the edge)
            @test cp1d.electrons.density ≈ D.density_thermal .+ D.density_fast .+ 6.0 .* cp1d.ion[2].density_thermal rtol = 1e-3
            push!(fast, copy(D.density_fast))
        end
        @test fast[end] ≈ fast[end-1] rtol = 1e-3
        S0 == 2e21 && @test fast[end][1] ≈ 0.8 * thermal0[1][1]

        # default: the electron and total ion densities are left alone
        dd = plasma_with_fast_source(S0)
        cp1d = dd.core_profiles.profiles_1d[1]
        ne0 = copy(cp1d.electrons.density_thermal)
        total0 = copy(cp1d.ion[1].density_thermal)
        for _ in 1:3
            IMAS.fast_particles_profiles!(dd.core_sources, cp1d)
            @test cp1d.electrons.density ≈ ne0
            @test cp1d.ion[1].density ≈ total0
        end
    end
end
