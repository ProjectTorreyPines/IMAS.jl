using IMAS
using Test

function rotation_sample()
    samp = joinpath(dirname(dirname(pathof(IMAS.IMASdd))), "sample")
    dd = IMAS.json2imas(joinpath(samp, "omas_sample.json"); show_warnings=false)
    dd.global_time = 0.0
    return dd
end

@testset "rotation: flux normalization" begin
    dd = rotation_sample()
    cp1d = dd.core_profiles.profiles_1d[]
    eqt1d = dd.equilibrium.time_slice[].profiles_1d
    # COCOS 11 stores the total poloidal flux: q = dΦ_tor/dψ_Wb, i.e. q = (1/2π) dΦ_tor/dψ_rad
    @test IMAS.gradient(eqt1d.psi, eqt1d.phi)[end÷2] ≈ eqt1d.q[end÷2] rtol = 1E-2
    @test IMAS.poloidal_flux_per_radian(cp1d) ≈ cp1d.grid.psi ./ (2π)
end

@testset "rotation: diamagnetic frequency" begin
    # synthetic profiles with p linear in ψ: ω_* = (1 / (Z e n)) dp/dψ_rad exactly
    dd = IMAS.dd()
    resize!(dd.core_profiles.profiles_1d, 1)
    cp1d = dd.core_profiles.profiles_1d[1]
    cp1d.time = 0.0
    dd.global_time = 0.0
    rho = collect(range(0.0, 1.0, 21))
    cp1d.grid.rho_tor_norm = rho
    cp1d.grid.psi = 2π .* (1.0 .+ 0.3 .* rho)  # ψ_rad = 1 + 0.3 ρ
    resize!(cp1d.ion, 1)
    ion = cp1d.ion[1]
    IMAS.ion_element!(ion, :C)
    ion.z_ion = 6.0
    n = 1E18
    ion.density_thermal = fill(n, length(rho))
    T = 2000.0 .- 1500.0 .* rho  # eV
    ion.temperature = T
    dp_dpsi = IMAS.mks.e * n * (-1500.0 / 0.3)
    @test IMAS.diamagnetic_rotation_frequency(cp1d, ion; rho_correction=0.0) ≈ fill(dp_dpsi / (6.0 * IMAS.mks.e * n), length(rho))
    @test IMAS.thermal_rotation_frequency(cp1d, ion; rho_correction=0.0) ≈ fill(-1500.0 / 0.3 / 6.0, length(rho))
    # the taper only acts inside rho_correction
    ωt = IMAS.diamagnetic_rotation_frequency(cp1d, ion; rho_correction=0.2)
    @test ωt[1] == 0.0
    @test ωt[end] ≈ dp_dpsi / (6.0 * IMAS.mks.e * n)
end

@testset "rotation: poloidal flow coefficient" begin
    @test IMAS.poloidal_flow_coefficient(0.0, 0.3) ≈ 1.17
    @test IMAS.poloidal_flow_coefficient(1E4, 0.3) ≈ -2.1 rtol = 1E-3
end

@testset "rotation: species round trip" begin
    dd = rotation_sample()
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    ωE = cp1d.rotation_frequency_tor_sonic
    @test any(!iszero, ωE)
    for poloidal_flow in (:none, :analytic)
        Ωs = [IMAS.sonic2ωtor(cp1d, ion; poloidal_flow, eqt) for ion in cp1d.ion]
        # main ion and carbon rotate differently with :none (different diamagnetic frequencies)
        poloidal_flow == :none && @test !(Ωs[1] ≈ Ωs[2])
        for (k, ion) in enumerate(cp1d.ion)
            ion.rotation_frequency_tor = Ωs[k]
        end
        # any species can be used as the reference to recover ω_E
        for k in eachindex(cp1d.ion)
            @test IMAS.ωtor2sonic(cp1d; ind=k, poloidal_flow, eqt) ≈ ωE
        end
        for ion in cp1d.ion
            empty!(ion, :rotation_frequency_tor)
        end
    end
    # :analytic: the main ion gets the k_i f ω_T term, impurities are locked to its parallel flow
    im = IMAS.main_ion_index(cp1d)
    @test cp1d.ion[im].label == "D"
    Δi = IMAS.rotation_frequency_offset(cp1d, cp1d.ion[im]; poloidal_flow=:analytic, eqt)
    Δz = IMAS.rotation_frequency_offset(cp1d, cp1d.ion[2]; poloidal_flow=:analytic, eqt)
    @test maximum(abs.(Δz .- Δi)) < 0.1 * maximum(abs.(Δi))
    @test_throws AssertionError IMAS.rotation_frequency_offset(cp1d, cp1d.ion[1]; poloidal_flow=:analytic)
    @test_throws ErrorException IMAS.rotation_frequency_offset(cp1d, cp1d.ion[1]; poloidal_flow=:bogus)
end

@testset "rotation: angular momentum and ∂/∂t source" begin
    dd = rotation_sample()
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    L = IMAS.toroidal_angular_momentum_density(cp1d, eqt)
    # rigid rotation part: Σ m n ⟨R²⟩ ω_E, corrected by the species diamagnetic offsets
    eqt1d = eqt.profiles_1d
    R2 = IMAS.interp1d(eqt1d.rho_tor_norm, eqt1d.gm8 .^ 2).(cp1d.grid.rho_tor_norm)
    L_rigid = sum(ion.element[1].a * IMAS.mks.m_p .* ion.density_thermal for ion in cp1d.ion) .* R2 .* cp1d.rotation_frequency_tor_sonic
    L_dia = sum(ion.element[1].a * IMAS.mks.m_p .* ion.density_thermal .* IMAS.rotation_frequency_offset(cp1d, ion) for ion in cp1d.ion) .* R2
    @test L ≈ L_rigid .+ L_dia

    # spin-up between time levels is a sink in the momentum equation, with units of a torque density
    Δt = 0.1
    cp1d_old = deepcopy(cp1d)
    cp1d_old.rotation_frequency_tor_sonic = 0.5 .* cp1d.rotation_frequency_tor_sonic
    IMAS.time_derivative_source!(dd, cp1d_old, Δt)
    source = dd.core_sources.source[findfirst(s -> s.identifier.name == "∂/∂t", dd.core_sources.source)]
    expected = -(L .- IMAS.toroidal_angular_momentum_density(cp1d_old, eqt; psi=IMAS.poloidal_flux_per_radian(cp1d))) ./ Δt
    @test source.profiles_1d[].momentum_tor ≈ expected
end
