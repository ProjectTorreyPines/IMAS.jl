import IMAS
using Test

@testset "time" begin
    dd=IMAS.dd()

    dd.global_time = 1010.0
    eqt=resize!(dd.equilibrium.time_slice)[end]
    eqt.global_quantities.ip=1010.0
    @test dd.equilibrium.time_slice[1010.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1].global_quantities.ip===eqt.global_quantities.ip
    @test_throws Exception dd.equilibrium.time_slice[2].global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip===eqt.global_quantities.ip

    dd.global_time = 2020.0
    eqt=resize!(dd.equilibrium.time_slice)[end]
    dd.equilibrium.time_slice[].global_quantities.ip=2020.0
    @test dd.equilibrium.time_slice[2020.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip!==eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[2].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip===eqt.global_quantities.ip

    dd.global_time = 3030.0
    eqt=resize!(dd.equilibrium.time_slice, IMAS.τ)[end]
    dd.equilibrium.time_slice[].global_quantities.ip=3030.0
    @test dd.equilibrium.time_slice[3030.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[3000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[1000.0].global_quantities.ip!==eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[10000.0].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[3].global_quantities.ip===eqt.global_quantities.ip
    @test dd.equilibrium.time_slice[].global_quantities.ip===eqt.global_quantities.ip
end
