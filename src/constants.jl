import PhysicalConstants.CODATA2018 as PCs

"""
Named tuple with physics constants:

    μ_0 = 1.25663706212e-6
    c = 2.99792458e8
    ϵ_0 = 8.8541878128e-12
    k_B = 1.380649e-23
    e = 1.602176634e-19
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    m_n = 1.67492749804e-27
    atm = 101325.0
    m_u = 1.6605390666e-27
"""
const constants = (
    μ_0=float(Float64, PCs.μ_0).val,
    c=float(Float64, PCs.c_0).val,
    ϵ_0=float(Float64, PCs.ε_0).val,
    k_B=float(Float64, PCs.k_B).val,
    e=float(Float64, PCs.e).val,
    m_e=float(Float64, PCs.m_e).val,
    m_p=float(Float64, PCs.m_p).val,
    m_n=float(Float64, PCs.m_n).val,
    atm=float(Float64, PCs.atm).val,
    m_u=float(Float64, PCs.m_u).val
)

const index_2_name__core_transport__model = Dict(
    0 => :unspecified, #Unspecified transport type
    1 => :combined, #Combination of data from available transport models. Representation of the total transport in the system
    2 => :transport_solver, #Output from a transport solver
    3 => :background, #Background transport level, ad-hoc transport model not directly related to a physics model
    4 => :database, #Transport specified by a database entry external to the dynamic evolution of the plasma
    5 => :neoclassical, #Neoclassical
    6 => :anomalous, #Representation of turbulent transport
    19 => :mhd, #Transport arising from MHD frequency modes
    20 => :ntm, #Transport arising from the presence of NTMs
    21 => :sawteeth, #Transport arising from the presence of sawteeth
    22 => :elm_continuous, #Continuous ELM model --- gives the ELM averaged profile
    23 => :elm_resolved, #Time resolved ELM model
    24 => :pedestal, #Transport level to give edge pedestal
    25 => :unknown) #Unknown transport type

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.core_transport__model}
    return index_2_name__core_transport__model
end

const index_2_name__core_sources__source = Dict(
    0 => :unspecified, #Unspecified source type
    1 => :total, #Total source; combines all sources
    2 => :nbi, #Source from Neutral Beam Injection
    3 => :ec, #Sources from electron cyclotron heating and current drive
    4 => :lh, #Sources from lower hybrid heating and current drive
    5 => :ic, #Sources from heating at the ion cyclotron range of frequencies
    6 => :fusion, #Sources from fusion reactions, e.g. alpha particle heating
    7 => :ohmic, #Source from ohmic heating
    8 => :bremsstrahlung, #Source from bremsstrahlung; radiation losses are negative sources
    9 => :synchrotron_radiation, #Source from synchrotron radiation; radiation losses are negative sources
    10 => :line_radiation, #Source from line radiation; radiation losses are negative sources
    11 => :collisional_equipartition, #Collisional equipartition
    12 => :cold_neutrals, #Source of cold neutrals
    13 => :bootstrap_current, #Bootstrap current
    14 => :pellet, #Sources from injection
    100 => :auxiliary, #Source from auxiliary systems, e.g. heating and current drive systems
    101 => :ic_nbi, #A combination of the ic and nbi sources
    102 => :ic_fusion, #A combination of the ic and fusion sources
    103 => :ic_nbi_fusion, #A combination of the ic and fusion sources
    104 => :ec_lh, #A combination of the ec and lh sources
    105 => :ec_ic, #A combination of the ec and ic sources
    106 => :lh_ic, #A combination of the lh and ic sources
    107 => :ec_lh_ic, #A combination of the ec, lh and ic sources
    108 => :gas_puff, #Gas puff
    109 => :killer_gas_puff, #Killer gas puff
    200 => :radiation, #Total radiation source; radiation losses are negative sources
    201 => :cyclotron_radiation, #Source from cyclotron radiation; radiation losses are negative sources
    202 => :cyclotron_synchrotron_radiation, #Source from combined cyclotron and synchrotron radiation; radiation losses are negative sources
    203 => :impurity_radiation, #Line radiation and Bremsstrahlung source; radiation losses are negative sources.
    303 => :particles_to_wall, #Particle pumping by the wall; negative source for plasma and positive source for the wall
    304 => :particles_to_pump, #Particle pumping by external pump; negative source for plasma and positive source for the pump
    305 => :charge_exchange, #Source from charge exchange. Charge exchange losses are negative sources
    400 => :transport, #Source term related to transport processes
    401 => :neoclassical, #Source term related to neoclassical processes
    402 => :equipartition, #Equipartition due to collisions and turbulence
    403 => :turbulent_equipartition, #Turbulent equipartition
    501 => :runaways, #Source from run-away processes; includes both electron and ion run-away
    601 => :ionisation, #Source from ionisation processes (not accounting for charge exchange)
    602 => :recombination, #Source from recombination processes (not accounting for charge exchange)
    603 => :excitation, #Source from excitation processes
    801 => :database, #Source from database entry
    802 => :gaussian) #Artificial source with a gaussian profile

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:core_sources__source}
    return index_2_name__core_sources__source
end

const index_2_name__stability__collection = Dict(
    1 => :default_limits,
    11 => :beta_limits, #Run all beta limits
    12 => :current_limit, #Run all current limits
    13 => :density_limit, #Run all density limits
)

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.stability__collection}
    return index_2_name__stability__collection
end

const index_2_name__stability__model = Dict(
    0 => :force_fail, #Instantly causes the actor to fail
    # 100s: Beta Limit Models
    101 => :beta_troyon_1984, #Beta limit defined by `F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209`
    102 => :beta_troyon_1985, #Beta limit defined by
    103 => :beta_tuda_1985, #Beta limit defined by
    104 => :beta_bernard_1983, #Beta limit defined by
    105 => :beta_model_105, #Beta limit defined by
    # 200s: Current Limit Models
    201 => :model_201, #Current limit defined by
    # 300s: Density Limit Models
    301 => :model_301, #Density limit defined by
    # 400s: Shaping Limit Models
    401 => :model_401, #Elongation limit defined by
    # 900s: Stability Codes
    999 => :unknown) #Unknown model type

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.stability__model}
    return index_2_name__stability__model
end

const index_2_name__balance_of_plant__power_electric_plant_operation = Dict((k - 1) => item for (k, item) in enumerate([:total, :HCD, :plant, :cryostat, :tritium_handling, :pumping, :pf_active]))

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:balance_of_plant__power_electric_plant_operation__system}
    return index_2_name__balance_of_plant__power_electric_plant_operation
end

"""
    name_2_index(ids::Union{IDS, IDSvector})

Return dict of name to IMAS indentifier.index
"""
function name_2_index(ids::Union{IDS,IDSvector})
    return Dict(v => k for (k, v) in index_2_name(ids))
end

"""
    findfirst(identifier_name::Symbol, ids::IDSvector)

Return item from IDSvector based on `ids.identifier.index` of `index_2_name(ids)`
"""
function Base.findfirst(identifier_name::Symbol, ids::IDSvector)
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids)))")
    elseif :identifier in fieldnames(eltype(ids))
        index = findfirst(idx -> idx.identifier.index == i, ids)
    else
        index = findfirst(idx -> idx.index == i, ids)
    end
    if index === nothing
        return nothing
    else
        return ids[index]
    end
end

"""
    findall(identifier_name::Symbol, ids::IDSvector)

Return items from IDSvector based on `ids.identifier.index` of `index_2_name(ids)`
"""
function Base.findall(identifier_name::Symbol, ids::IDSvector)
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids)))")
    elseif :identifier in fieldnames(eltype(ids))
        indexes = findall(idx -> idx.identifier.index == i, ids)
    else
        indexes = findall(idx -> idx.index == i, ids)
    end
    return eltype(ids)[ids[index] for index in indexes]
end

"""
    Base.resize!(@nospecialize(ids::IDSvector{T}), identifier_name::Symbol, conditions::Pair{String}...; wipe::Bool=true, error_multiple_matches::Bool=true)::T where {T<:IDSvectorElement}

Resize ids if `identifier_name` is not found based on `ids.identifier.index` of `index_2_name(ids)` and a set of conditions are not met.

If wipe=true and an entry matching the condition is found, then the content of the matching IDS is emptied.

Either way, the IDS is populated with the conditions.

NOTE: `error_multiple_matches` will delete all extra entries matching the conditions.

Returns the selected IDS
"""
function Base.resize!(@nospecialize(ids::IDSvector{T}), identifier_name::Symbol, conditions::Pair{String}...; wipe::Bool=true, error_multiple_matches::Bool=true)::T where {T<:IDSvectorElement}
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids)))")
    elseif :identifier in fieldnames(eltype(ids))
        return resize!(ids, "identifier.index" => i, conditions...; wipe, error_multiple_matches)
    else
        return resize!(ids, "index" => i, conditions...; wipe, error_multiple_matches)
    end
end

"""
    Base.deleteat!(@nospecialize(ids::T), identifier_name::Symbol)::T where {T<:IDSvector}

Deletes all entries that match based on `ids.identifier.index` of `index_2_name(ids)`
"""
function Base.deleteat!(@nospecialize(ids::T), identifier_name::Symbol, conditions::Pair{String}...)::T where {T<:IDSvector}
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids)))")
    elseif :identifier in fieldnames(eltype(ids))
        return deleteat!(ids, "identifier.index" => i, conditions...)
    else
        return deleteat!(ids, "index" => i, conditions...)
    end
end
