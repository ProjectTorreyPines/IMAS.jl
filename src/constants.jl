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

"""
    name_2_index(ids::Union{IDS, IDSvector})

returns dict of name to IMAS indentifier.index
"""
function name_2_index(ids::Union{IDS,IDSvector})
    return Dict(v => k for (k, v) in index_2_name(ids))
end

"""
    findfirst(identifier_name::Symbol, haystack::IDSvector)

return item from IDSvector based on identifier.index of index_2_name
"""
function Base.findfirst(identifier_name::Symbol, haystack::IDSvector)
    i = name_2_index(haystack)[identifier_name]
    index = findfirst(idx -> idx.identifier.index == i, haystack)
    return haystack[index]
end

"""
    findall(identifier_name::Symbol, haystack::IDSvector)

return items from IDSvector based on identifier.index of index_2_name
"""
function Base.findall(identifier_name::Symbol, haystack::IDSvector)
    i = name_2_index(haystack)[identifier_name]
    indexes = findall(idx -> idx.identifier.index == i, haystack)
    return [haystack[index] for index in indexes]
end