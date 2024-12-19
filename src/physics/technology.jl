document[Symbol("Physics technology")] = Symbol[]

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

"""
    mechanical_technology(dd::IMAS.dd, what::Symbol)

Set `yield_strength`, `poisson_ratio`, and `young_modulus` properties for `:pl`, `:oh` and `:tf`
"""
function mechanical_technology(dd::IMAS.dd, what::Symbol)
    @assert what in (:pl, :oh, :tf)
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
        coil_tech.material = "copper"
        coil_tech.temperature = 293.0
        coil_tech.fraction_steel = 0.0
        coil_tech.ratio_SC_to_copper = 0.0
        coil_tech.fraction_void = 0.2

    elseif technology ∈ (:nb3sn, :nbti, :nb3sn_iter, :nb3sn_kdemo, :rebco)
        if technology == :nb3sn
            coil_tech.temperature = 4.2
            coil_tech.material = "nb3sn"
            coil_tech.fraction_void = 0.1
        elseif technology == :nbti
            coil_tech.temperature = 4.2
            coil_tech.material = "nbti"
            coil_tech.fraction_void = 0.2 # from Supercond. Sci. Technol. 36 (2023) 075009
        elseif technology == :nb3sn_iter
            coil_tech.temperature = 4.2
            coil_tech.material = "nb3sn_iter"
            coil_tech.fraction_void = 0.1
        elseif technology == :nb3sn_kdemo
            coil_tech.temperature = 4.2
            coil_tech.material = "nb3sn_kdemo"
            if coil_type == :tf
                coil_tech.fraction_void = 0.26 # from NF 55 (2015) 053027, Table 2
            end
        else
            coil_tech.temperature = 4.2
            coil_tech.material = "rebco"
        end
        coil_tech.fraction_steel = 0.5
        coil_tech.ratio_SC_to_copper = 1.0
        coil_tech.fraction_void = 0.1
    end

    if technology == :nb3sn_iter
        if coil_type == :oh
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_steel = 0.46
        elseif coil_type == :tf
            coil_tech.thermal_strain = -0.69
            coil_tech.JxB_strain = -0.13
            coil_tech.fraction_steel = 0.55
        elseif coil_type == :pf_active
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_steel = 0.46
        end
    end

    coil_tech.thermal_strain = 0.0
    coil_tech.JxB_strain = 0.0

    return coil_tech
end

@compat public coil_technology
push!(document[Symbol("Physics technology")], :coil_technology)

"""
    fraction_conductor(coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})

returns the fraction of (super)conductor in a coil technology
"""
function fraction_conductor(coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})
    frac = 1.0 - coil_tech.fraction_steel - coil_tech.fraction_void # fraction of coil that is a conductor
    @assert frac > 0.0 "coil technology has no room for conductor (fraction_steel = $(coil_tech.fraction_steel), fraction_void = $(coil_tech.fraction_void))"
    if coil_tech.material == "copper"
        return frac
    else
        return frac * coil_tech.ratio_SC_to_copper / (1.0 + coil_tech.ratio_SC_to_copper) # fraction of coil that is superconductor
    end
end

@compat public fraction_conductor
push!(document[Symbol("Physics technology")], :fraction_conductor)

"""
    GAMBL_blanket(bm::IMAS.blanket__module)

Define layers for the `dd.blanket.module` for a GAMBL type blanket technology
"""
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

@compat public GAMBL_blanket
push!(document[Symbol("Physics technology")], :GAMBL_blanket)
