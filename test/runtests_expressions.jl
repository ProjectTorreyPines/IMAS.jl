using IMAS
using Test

# ---------------------------------------------------------------------------------------------------------------------
# Dynamic expressions: evaluated on the sample data, compared with the functions they call.
# ---------------------------------------------------------------------------------------------------------------------

@testset "expressions" begin
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "omas_sample.json"); show_warnings=false)
    fw = IMAS.first_wall(dd.wall)
    IMAS.flux_surfaces(dd.equilibrium, fw.r, fw.z)
    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    @testset "core_profiles.global_quantities.v_loop" begin
        empty!(dd.core_profiles.global_quantities, :v_loop)
        v_loop = dd.core_profiles.global_quantities.v_loop
        @test length(v_loop) == length(dd.core_profiles.profiles_1d)
        @test all(isfinite, v_loop)
        k = findfirst(x -> x === cp1d, dd.core_profiles.profiles_1d)
        @test v_loop[k] == IMAS.vloop(cp1d, eqt)
    end

    @testset "pulse_schedule.ec.power_launched.reference" begin
        ec = dd.pulse_schedule.ec
        empty!(ec)
        ec.time = [0.0, 1.0]
        resize!(ec.beam, 2)
        ec.beam[1].power_launched.reference = [1.0e6, 2.0e6]
        ec.beam[2].power_launched.reference = [3.0e6, 5.0e6]
        @test ec.power_launched.reference == [4.0e6, 7.0e6]
    end

    @testset "summary.global_quantities.beta_tor.value" begin
        empty!(dd.summary)
        dd.summary.time = [eqt.time]
        beta_tor = dd.summary.global_quantities.beta_tor.value
        @test beta_tor == [IMAS.beta_tor(eqt, cp1d; norm=false, thermal=false)]
        @test all(isfinite, beta_tor)
        # the normalized value is the same pressure average, scaled by a * B0 / ip
        beta_tor_norm = dd.summary.global_quantities.beta_tor_norm.value
        B0 = eqt.global_quantities.vacuum_toroidal_field.b0
        ip = eqt.global_quantities.ip
        @test beta_tor_norm[1] ≈ beta_tor[1] * eqt.boundary.minor_radius * abs(B0) / abs(ip / 1e6) * 1.0e2
    end
end
