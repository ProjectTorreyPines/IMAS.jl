"""
Material properties
"""
Base.@kwdef struct MaterialProperties
    yield_strength::Float64 = NaN
    young_modulus::Float64 = NaN
    poisson_ratio::Float64 = NaN
end

const stainless_steel = MaterialProperties(;
    yield_strength=800E6, # Pa
    young_modulus=193.103448275E9, # Pa
    poisson_ratio=0.33
)

const pure_copper = MaterialProperties(;
    yield_strength=70E6, # Pa
    young_modulus=110E9, # Pa
    poisson_ratio=0.34
)

function mechanical_technology(dd::IMAS.dd, what::Symbol)
    if what != :pl && getproperty(dd.build, what).technology.material == "copper"
        material = pure_copper
    else
        material = stainless_steel
    end
    setproperty!(dd.solid_mechanics.center_stack.properties.yield_strength, what, material.yield_strength)
    setproperty!(dd.solid_mechanics.center_stack.properties.poisson_ratio, what, material.poisson_ratio)
    return setproperty!(dd.solid_mechanics.center_stack.properties.young_modulus, what, material.young_modulus)
end

"""
    coil_technology(technology::Symbol, coil_type::Symbol)

Return coil parameters from technology and coil type [:oh, :tf, :pf_active]"
"""
function coil_technology(coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology}, technology::Symbol, coil_type::Symbol)
    if coil_type ∉ (:oh, :tf, :pf_active)
        error("Supported coil type are [:oh, :tf, :pf_active]")
    end

    if technology == :copper
        resize!(coil_tech.material, 2)
        coil_tech.temperature = 293.0
        coil_tech.material[1].name = "copper"
        coil_tech.material[1].composition = 0.8

        coil_tech.material[2].name = "vacuum"
        coil_tech.material[2].composition = 0.2

    elseif technology ∈ (:nb3sn, :nbti, :nb3sn_iter, :nb3sn_kdemo, :rebco)
        resize!(coil_tech.material, 4)
        if technology == :nb3sn
            coil_tech.temperature = 4.2
            coil_tech.material[1].name = "nb3sn"
            coil_tech.material[1].composition = 0.2

            coil_tech.material[2].name = "vacuum"
            coil_tech.material[2].composition = 0.1

            coil_tech.material[3].name = "steel"
            coil_tech.material[3].composition = 0.5

            coil_tech.material[4].name = "copper"
            coil_tech.material[4].composition = 0.2
        elseif technology == :nbti
            coil_tech.temperature = 4.2
            coil_tech.material[1].name = "nbti"
            coil_tech.material[1].composition = 0.15

            coil_tech.material[2].name = "vacuum"
            coil_tech.material[2].composition = 0.2 # from Supercond. Sci. Technol. 36 (2023) 075009

            coil_tech.material[3].name = "steel"
            coil_tech.material[3].composition = 0.5

            coil_tech.material[4].name = "copper"
            coil_tech.material[4].composition = 0.15
        elseif technology == :nb3sn_iter
            coil_tech.temperature = 4.2
            coil_tech.material[1].name = "nb3sn"
            coil_tech.material[1].composition = 0.2

            coil_tech.material[2].name = "vacuum"
            coil_tech.material[2].composition = 0.1

            coil_tech.material[3].name = "steel"
            coil_tech.material[3].composition = 0.5

            coil_tech.material[4].name = "copper"
            coil_tech.material[4].composition = 0.2
        elseif technology == :nb3sn_kdemo
            coil_tech.temperature = 4.2
            coil_tech.material[1].name = "nb3sn_kdemo"
            coil_tech.material[1].composition = 0.2

            coil_tech.material[2].name = "vacuum"
            coil_tech.material[2].composition = 0.1

            coil_tech.material[3].name = "steel"
            coil_tech.material[3].composition = 0.5

            coil_tech.material[4].name = "copper"
            coil_tech.material[4].composition = 0.2

            if coil_type == :tf
                coil_tech.material[1].name = "nb3sn_kdemo"
                coil_tech.material[1].composition = 0.12

                coil_tech.material[2].name = "vacuum"
                coil_tech.material[2].composition = 0.26 # from NF 55 (2015) 053027, Table 2

                coil_tech.material[3].name = "steel"
                coil_tech.material[3].composition = 0.5

                coil_tech.material[4].name = "copper"
                coil_tech.material[4].composition = 0.12 
            end
        else
            coil_tech.temperature = 4.2
            coil_tech.material[1].name = "rebco"
            coil_tech.material[1].composition = 0.2

            coil_tech.material[2].name = "vacuum"
            coil_tech.material[2].composition = 0.1

            coil_tech.material[3].name = "steel"
            coil_tech.material[3].composition = 0.5

            coil_tech.material[4].name = "copper"
            coil_tech.material[4].composition = 0.2

        end
    end

    if technology == :nb3sn_iter
        resize!(coil_tech.material, 4)
        coil_tech.material[1].name = "nb3sn_iter"

        coil_tech.material[2].name = "vacuum"
        coil_tech.material[2].composition = 0.1

        coil_tech.material[3].name = "steel"

        coil_tech.material[4].name = "copper"

        if coil_type == :oh
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05

            coil_tech.material[1].composition = 0.22
            coil_tech.material[3].composition = 0.46
            coil_tech.material[4].composition = 0.22
        elseif coil_type == :tf
            coil_tech.thermal_strain = -0.69
            coil_tech.JxB_strain = -0.13

            coil_tech.material[1].composition = 0.175
            coil_tech.material[3].composition = 0.55
            coil_tech.material[4].composition = 0.175
        elseif coil_type == :pf_active
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05

            coil_tech.material[1].composition = 0.22
            coil_tech.material[3].composition = 0.46
            coil_tech.material[4].composition = 0.22
        end
    end

    coil_tech.thermal_strain = 0.0
    coil_tech.JxB_strain = 0.0

    return coil_tech
end


function GAMBL_blanket(bm::IMAS.blanket__module)
    layers = resize!(bm.layer, 3)

    n = 1
    layers[n].name = "First wall"
    layers[n].material = "tungsten"
    layers[n].thickness = 0.02

    n = n + 1
    layers[n].name = "Breeder"
    layers[n].material = "lithium-lead"
    layers[n].thickness = 0.5

    n = n + 1
    layers[n].name = "Shield"
    layers[n].material = "tungsten"
    layers[n].thickness = 0.05

    return bm
end
