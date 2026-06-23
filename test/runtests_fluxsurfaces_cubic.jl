using IMAS
using Test

# Analytic field: nested ellipses ψ = ((R-R0)/a0)^2 + (Z/b0)^2.
# Level-L contour is an ellipse, semi-axes (a0√L in R, b0√L in Z),
# enclosed area = π·a0·b0·L (closed-form check for tracing tests).
function ellipse_itp(R0, a0, b0; n=65)
    r = range(R0 - 1.4a0, R0 + 1.4a0, length=n)
    z = range(-1.4b0, 1.4b0, length=n)
    PSI = [((ri - R0) / a0)^2 + (zj / b0)^2 for ri in r, zj in z]
    return IMAS.ψ_interpolant(r, z, PSI).PSI_interpolant, r, z
end

@testset "fluxsurfaces_cubic" begin
    R0, a0, b0 = 1.7, 0.5, 1.0
    itp, _, _ = ellipse_itp(R0, a0, b0)

    @testset "_project_to_level snaps onto ψ=c" begin
        c = 0.36
        # a point off the c-contour (interior), should land on ψ=c
        seed = (R0 + 0.3a0, 0.1b0)
        (R, Z), ok = IMAS._project_to_level(itp, c, seed)
        @test ok
        val, _ = IMAS.FI.value_gradient(itp, (R, Z))
        @test isapprox(val, c; atol=1e-9)
    end

    @testset "_project_to_level reports failure at a critical point" begin
        # center of the ellipse: ∇ψ = 0 -> cannot project
        (_, _), ok = IMAS._project_to_level(itp, 0.36, (R0, 0.0))
        @test ok == false
    end

    @testset "_contour_tangent is unit and ⊥ ∇ψ" begin
        x = (R0 + a0 * sqrt(0.36), 0.0)         # OMP of the ψ=0.36 ellipse
        t = IMAS._contour_tangent(itp, x, 1)
        g = IMAS.FI.gradient(itp, x)
        @test isapprox(hypot(t[1], t[2]), 1.0; atol=1e-10)
        @test isapprox(t[1]*g[1] + t[2]*g[2], 0.0; atol=1e-8)   # tangent ⟂ gradient
    end

    @testset "_contour_curvature matches a circle (κ = 1/radius)" begin
        citp, _, _ = ellipse_itp(R0, 0.5, 0.5)   # a0=b0 -> circles
        c = 0.25                                   # radius = 0.5*sqrt(0.25) = 0.25
        x = (R0 + 0.5 * sqrt(c), 0.0)
        κ = IMAS._contour_curvature(citp, x)
        @test isapprox(abs(κ), 1 / 0.25; rtol=1e-3)
    end
end
