import IMAS
import IMAS: @timedep
using Test

@testset "time_ids" begin
    dd = IMAS.dd()

    dd.global_time = 1010.0
    eqt = resize!(dd.equilibrium.time_slice)[end]
    eqt.global_quantities.ip = 1010.0
    @test dd.equilibrium.time_slice[1010.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[2].global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip

    dd.global_time = 2020.0
    eqt = resize!(dd.equilibrium.time_slice)[end]
    dd.equilibrium.time_slice[].global_quantities.ip = 2020.0
    @test dd.equilibrium.time_slice[2020.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[3].global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip

    dd.global_time = 3030.0
    eqt = resize!(dd.equilibrium.time_slice, IMAS.τ)[end]
    dd.equilibrium.time_slice[].global_quantities.ip = 3030.0
    @test dd.equilibrium.time_slice[3030.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[3000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[3].global_quantities.ip === eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[4].global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip

    dd.global_time = 2020.0
    eqt = dd.equilibrium.time_slice[2]
    @test dd.equilibrium.time_slice[2020.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip !== eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[IMAS.τ].global_quantities.ip === eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip === eqt.global_quantities.ip
end

@testset "time_array" begin
    dd = IMAS.dd()

    dd.global_time = 1010.0
    @test @timedep(dd.equilibrium.vacuum_toroidal_field.b0 = 1.0) == 1.0
    @test dd.equilibrium.time == [1010.0]
    @test @timedep(dd.equilibrium.vacuum_toroidal_field.b0) == 1.0
    @test dd.equilibrium.vacuum_toroidal_field.b0 == [1.0]

    dd.global_time = 2010.0
    @test @timedep(dd.equilibrium.vacuum_toroidal_field.b0 = 2.0) == 2.0
    @test dd.equilibrium.time == [1010.0, 2010.0]
    @test @timedep(dd.equilibrium.vacuum_toroidal_field.b0) == 2.0
    @test dd.equilibrium.vacuum_toroidal_field.b0 == [1.0, 2.0]
end
