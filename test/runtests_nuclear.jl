using IMAS
using Test

@testset "nuclear" begin
    @testset "reactivity reproduces Bosch-Hale Table VIII" begin
        # Table VIII of Bosch & Hale, Nucl. Fusion 32 (1992) 611, as reprinted in the erratum
        # Nucl. Fusion 33 (1993) 1919: <σv> [cm³/s] at Ti [keV] for D(t,n)α, ³He(d,p)α, D(d,p)T, D(d,n)³He
        published = (
            1.0 => (6.857e-21, 3.057e-26, 1.017e-22, 9.933e-23),
            2.0 => (2.977e-19, 1.399e-23, 3.150e-21, 3.110e-21),
            5.0 => (1.366e-17, 6.377e-21, 9.024e-20, 9.128e-20),
            10.0 => (1.136e-16, 2.126e-19, 5.781e-19, 6.023e-19),
            20.0 => (4.330e-16, 3.482e-18, 2.399e-18, 2.603e-18),
            50.0 => (8.649e-16, 5.554e-17, 9.838e-18, 1.133e-17)
        )
        models = ("D+T→He4", "D+He3→He4", "D+D→T", "D+D→He3")
        for (T_keV, values) in published, (model, ref) in zip(models, values)
            @test IMAS.reactivity([1e3 * T_keV], model)[1] * 1e6 ≈ ref rtol = 1e-3
        end
    end
end
