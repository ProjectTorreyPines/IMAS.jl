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
    avog = 6.02214076e23
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
    m_u=float(Float64, PCs.m_u).val,
    avog=float(Float64, PCs.AvogadroConstant).val
)

const index_2_name__core_transport__model = Dict(
    0 => :unspecified, # Unspecified transport type
    1 => :combined, # Combination of data from available transport models. Representation of the total transport in the system
    2 => :transport_solver, # Output from a transport solver
    3 => :background, # Background transport level, ad-hoc transport model not directly related to a physics model
    4 => :database, # Transport specified by a database entry external to the dynamic evolution of the plasma
    5 => :neoclassical, # Neoclassical
    6 => :anomalous, # Representation of turbulent transport
    19 => :mhd, # Transport arising from MHD frequency modes
    20 => :ntm, # Transport arising from the presence of NTMs
    21 => :sawteeth, # Transport arising from the presence of sawteeth
    22 => :elm_continuous, # Continuous ELM model --- gives the ELM averaged profile
    23 => :elm_resolved, # Time resolved ELM model
    24 => :pedestal, # Transport level to give edge pedestal
    25 => :unknown) #Unknown transport type

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.core_transport__model}
    return index_2_name__core_transport__model
end

const index_2_name__core_sources__source = Dict(
    0 => :unspecified, # Unspecified source type
    1 => :total, # Total source; combines all sources
    2 => :nbi, # Source from Neutral Beam Injection
    3 => :ec, # Sources from electron cyclotron heating and current drive
    4 => :lh, # Sources from lower hybrid heating and current drive
    5 => :ic, # Sources from heating at the ion cyclotron range of frequencies
    6 => :fusion, # Sources from fusion reactions, e.g. alpha particle heating
    7 => :ohmic, # Source from ohmic heating
    8 => :bremsstrahlung, # Source from bremsstrahlung; radiation losses are negative sources
    9 => :synchrotron_radiation, # Source from synchrotron radiation; radiation losses are negative sources
    10 => :line_radiation, # Source from line radiation; radiation losses are negative sources
    11 => :collisional_equipartition, # Collisional equipartition
    12 => :cold_neutrals, # Source of cold neutrals
    13 => :bootstrap_current, # Bootstrap current
    14 => :pellet, # Sources from injection
    100 => :auxiliary, # Source from auxiliary systems, e.g. heating and current drive systems
    101 => :ic_nbi, # A combination of the ic and nbi sources
    102 => :ic_fusion, # A combination of the ic and fusion sources
    103 => :ic_nbi_fusion, # A combination of the ic and fusion sources
    104 => :ec_lh, # A combination of the ec and lh sources
    105 => :ec_ic, # A combination of the ec and ic sources
    106 => :lh_ic, # A combination of the lh and ic sources
    107 => :ec_lh_ic, # A combination of the ec, lh and ic sources
    108 => :gas_puff, # Gas puff
    109 => :killer_gas_puff, # Killer gas puff
    200 => :radiation, # Total radiation source; radiation losses are negative sources
    201 => :cyclotron_radiation, # Source from cyclotron radiation; radiation losses are negative sources
    202 => :cyclotron_synchrotron_radiation, # Source from combined cyclotron and synchrotron radiation; radiation losses are negative sources
    203 => :impurity_radiation, # Line radiation and Bremsstrahlung source; radiation losses are negative sources.
    303 => :particles_to_wall, # Particle pumping by the wall; negative source for plasma and positive source for the wall
    304 => :particles_to_pump, # Particle pumping by external pump; negative source for plasma and positive source for the pump
    305 => :charge_exchange, # Source from charge exchange. Charge exchange losses are negative sources
    400 => :transport, # Source term related to transport processes
    401 => :neoclassical, # Source term related to neoclassical processes
    402 => :equipartition, # Equipartition due to collisions and turbulence
    403 => :turbulent_equipartition, # Turbulent equipartition
    409 => :time_derivative, # Source based on d/dt term in transport equations
    501 => :runaways, # Source from run-away processes; includes both electron and ion run-away
    601 => :ionisation, # Source from ionisation processes (not accounting for charge exchange)
    602 => :recombination, # Source from recombination processes (not accounting for charge exchange)
    603 => :excitation, # Source from excitation processes
    801 => :database, # Source from database entry
    802 => :gaussian) #Artificial source with a gaussian profile

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:core_sources__source}
    return index_2_name__core_sources__source
end

const index_2_name__pellets__launcher___shape__type = Dict(
    1 => :spherical,
    2 => :cylindrical,
    3 => :rectangular)

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:pellets__launcher___shape}
    return index_2_name__pellets__launcher___shape__type
end

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:pellets__launcher___shape__type}
    return index_2_name__pellets__launcher___shape__type
end

const index_2_name__stability__collection = Dict(
    1 => :default_limits,
    11 => :beta_limits, # Run all beta limits
    12 => :current_limits, # Run all current limits
    13 => :density_limits #Run all density limits
)

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.stability__collection}
    return index_2_name__stability__collection
end

const index_2_name__stability__model = Dict(
    0 => :force_fail, # Instantly causes the actor to fail
    # 100s: Beta Limit Models
    101 => :beta_troyon_1984, # Beta limit defined by `F Troyon et al 1984 Plasma Phys. Control. Fusion 26 209`
    102 => :beta_troyon_1985,
    103 => :beta_tuda_1985,
    104 => :beta_bernard_1983,
    105 => :beta_model_105,
    # 200s: Current Limit Models
    201 => :q95_gt_2, # q(@rho=0.95) > 2.
    202 => :q08_gt_2, # q(@rho=0.8) > 2.
    # 300s: Density Limit Models
    301 => :gw_density, # Density limit defined by Greenwald fraction
    # 400s: Shaping Limit Models
    401 => :κ_controllability,
    # 900s: Stability Codes
    999 => :unknown) #Unknown model type

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.stability__model}
    return index_2_name__stability__model
end

const index_2_name__equilibrium__time_slice___profiles_2d___grid_type = Dict(
    1 => :rectangular, # Cylindrical R,Z ala eqdsk (R=dim1, Z=dim2). In this case the position arrays should not be filled since they are redundant with grid/dim1 and dim2.
    2 => :inverse, # Rhopolar_polar 2D polar coordinates (rho=dim1, theta=dim2) with magnetic axis as centre of grid; theta and values following the COCOS=11 convention; the polar angle is theta=atan2(z-zaxis,r-raxis)
    11 => :inverse_psi_straight_field_line, # Flux surface type with psi as radial label (dim1) and the straight-field line poloidal angle (mod(index,10)=1) (dim2); could be non-equidistant; magnetic axis as centre of grid; following the COCOS=11 convention
    12 => :inverse_psi_equal_arc, # Flux surface type with psi as radial label (dim1) and the equal arc poloidal angle (mod(index,10)=2) (dim2)
    13 => :inverse_psi_polar, # Flux surface type with psi as radial label (dim1) and the polar poloidal angle (mod(index,10)=3) (dim2); could be non-equidistant
    14 => :inverse_psi_straight_field_line_fourier, # Flux surface type with psi as radial label (dim1) and Fourier modes in the straight-field line poloidal angle (mod(index,10)=4) (dim2), could be non-equidistant; magnetic axis as centre of grid; following the COCOS=11 convention
    15 => :inverse_psi_equal_arc_fourier, # Flux surface type with psi as radial label (dim1) and Fourier modes in the equal arc poloidal angle (mod(index,10)=5) (dim2)
    16 => :inverse_psi_polar_fourier, # Flux surface type with psi as radial label (dim1) and Fourier modes in the polar poloidal angle (mod(index,10)=6) (dim2); could be non-equidistant
    21 => :inverse_rhopolnorm_straight_field_line, # Flux surface type with radial label sqrt[(psi-psi_axis)/(psi_edge-psi_axis)] (dim1) and the straight-field line poloidal angle (dim2)
    22 => :inverse_rhopolnorm_equal_arc, # Flux surface type with radial label sqrt[(psi-psi_axis)/(psi_edge-psi_axis)] (dim1) and the equal arc poloidal angle (dim2)
    23 => :inverse_rhopolnorm_polar, # Flux surface type with radial label sqrt[(psi-psi_axis)/(psi_edge-psi_axis)] (dim1) and the polar poloidal angle (dim2)
    24 => :inverse_rhopolnorm_straight_field_line_fourier, # Flux surface type with radial label sqrt[(psi-psi_axis)/(psi_edge-psi_axis)] (dim1) and Fourier modes in the straight-field line poloidal angle (dim2)
    25 => :inverse_rhopolnorm_equal_arc_fourier, # Flux surface type with radial label sqrt[(psi-psi_axis)/(psi_edge-psi_axis)] (dim1) and Fourier modes in the equal arc poloidal angle (dim2)
    26 => :inverse_rhopolnorm_polar_fourier, # Flux surface type with radial label sqrt[(psi-psi_axis)/(psi_edge-psi_axis)] (dim1) and Fourier modes in the polar poloidal angle (dim2)
    31 => :inverse_rhotornorm_straight_field_line, # Flux surface type with radial label sqrt[Phi/Phi_edge] (dim1) and the straight-field line poloidal angle (dim2)
    32 => :inverse_rhotornorm_equal_arc, # Flux surface type with radial label sqrt[Phi/Phi_edge] (dim1) and the equal arc poloidal angle (dim2)
    33 => :inverse_rhotornorm_polar, # Flux surface type with radial label sqrt[Phi/Phi_edge] (dim1) and the polar poloidal angle (dim2)
    34 => :inverse_rhotornorm_straight_field_line_fourier, # Flux surface type with radial label sqrt[Phi/Phi_edge] (dim1) and Fourier modes in the straight-field line poloidal angle (dim2)
    35 => :inverse_rhotornorm_equal_arc_fourier, # Flux surface type with radial label sqrt[Phi/Phi_edge] (dim1) and Fourier modes in the equal arc poloidal angle (dim2)
    36 => :inverse_rhotornorm_polar_fourier, # Flux surface type with radial label sqrt[Phi/Phi_edge] (dim1) and Fourier modes in the polar poloidal angle (dim2)
    41 => :inverse_rhopol_straight_field_line, # Flux surface type with radial label sqrt[psi-psi_axis] (dim1) and the straight-field line poloidal angle (dim2)
    42 => :inverse_rhopol_equal_arc, # Flux surface type with radial label sqrt[psi-psi_axis] (dim1) and the equal arc poloidal angle (dim2)
    43 => :inverse_rhopol_polar, # Flux surface type with radial label sqrt[psi-psi_axis] (dim1) and the polar poloidal angle (dim2)
    44 => :inverse_rhopol_straight_field_line_fourier, # Flux surface type with radial label sqrt[psi-psi_axis] (dim1) and Fourier modes in the straight-field line poloidal angle (dim2)
    45 => :inverse_rhopol_equal_arc_fourier, # Flux surface type with radial label sqrt[psi-psi_axis] (dim1) and Fourier modes in the equal arc poloidal angle (dim2)
    46 => :inverse_rhopol_polar_fourier, # Flux surface type with radial label sqrt[psi-psi_axis] (dim1) and Fourier modes in the polar poloidal angle (dim2)
    51 => :inverse_rhotor_straight_field_line, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and the straight-field line poloidal angle (dim2)
    52 => :inverse_rhotor_equal_arc, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and the equal arc poloidal angle (dim2)
    53 => :inverse_rhotor_polar, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and the polar poloidal angle (dim2)
    54 => :inverse_rhotor_straight_field_line_fourier, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and Fourier modes in the straight-field line poloidal angle (dim2)
    55 => :inverse_rhotor_equal_arc_fourier, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and Fourier modes in the equal arc poloidal angle (dim2)
    56 => :inverse_rhotor_polar_fourier, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and Fourier modes in the polar poloidal angle (dim2)
    57 => :inverse_rhotor_mxh, # Flux surface type with radial label sqrt[Phi/pi/B0] (dim1), Phi being toroidal flux, and MXH coefficients R0, Z0, ϵ, κ, c0, c[...], s[...] describing each flux surface (dim2)
    91 => :irregular_rz_na) # Irregular grid, thus give list of vertices in dim1(1:ndim1), dim2(1:ndim1) and then all fields are on values(1:ndim1,1)

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.equilibrium__time_slice___profiles_2d}
    return index_2_name__equilibrium__time_slice___profiles_2d___grid_type
end

const index_2_name__pf_active__coil___function = Dict(
    0 => :flux, # Generate flux (drive toroidal current)
    1 => :shaping, # Generate magnetic field for shaping
    2 => :vertical) # Generate magnetic field for vertical force balance

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.pf_active__coil___function}
    return index_2_name__pf_active__coil___function
end

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.pf_active__coil}
    return index_2_name__pf_active__coil___function
end

const index_2_name__pf_active__coil___element___geometry__geometry_type = Dict(
    1 => :outline,
    2 => :rectangle,
    3 => :oblique,
    4 => :arcs_of_circle,
    5 => :annulus,
    6 => :thick_line)

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.pf_active__coil___element___geometry}
    return index_2_name__pf_active__coil___element___geometry__geometry_type
end

const index_2_name__pf_passive__loop___element___geometry__geometry_type = Dict(
    1 => :outline,
    2 => :rectangle,
    3 => :oblique,
    4 => :arcs_of_circle,
    5 => :annulus,
    6 => :thick_line)

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:IMAS.pf_passive__loop___element___geometry}
    return index_2_name__pf_passive__loop___element___geometry__geometry_type
end

const index_2_name__balance_of_plant__power_electric_plant_operation =
    Dict((k - 1) => item for (k, item) in enumerate((:total, :HCD, :plant, :cryostat, :tritium_handling, :pumping, :pf_active)))

function index_2_name(ids::Union{T,IDSvector{T}}) where {T<:balance_of_plant__power_electric_plant_operation__system}
    return index_2_name__balance_of_plant__power_electric_plant_operation
end

# ============= #

"""
    identifier_index(@nospecialize(ids::IDS); error_on_missing::Bool=true)

Return ids.identifier.index`
"""
function identifier_index(@nospecialize(ids::IDS); error_on_missing::Bool=true)
    if :identifier in fieldnames(typeof(ids))
        if hasdata(ids.identifier, :index) || error_on_missing
            return ids.identifier.index
        else
            return nothing
        end
    elseif :grid_type in fieldnames(typeof(ids))
        if hasdata(ids.grid_type, :index) || error_on_missing
            return ids.grid_type.index
        else
            return nothing
        end
    elseif :type in fieldnames(typeof(ids))
        if hasdata(ids.type, :index) || error_on_missing
            return ids.type.index
        else
            return nothing
        end
    elseif :geometry_type in fieldnames(typeof(ids))
        if hasdata(ids, :geometry_type) || error_on_missing
            return ids.geometry_type
        else
            return nothing
        end
    elseif :index in fieldnames(typeof(ids))
        if hasdata(ids, :index) || error_on_missing
            return ids.index
        else
            return nothing
        end
    else
        error("$(ulocation(ids)) does not have a `.index` field")
    end
end

function identifier_index(@nospecialize(ids::IDS), default::Int)
    index = identifier_index(ids; error_on_missing=false)
    if index === nothing
        return default
    else
        return index
    end
end

"""
    identifier_name(@nospecialize(ids::IDS); error_on_missing::Bool=true)

Return name (Symbol) based on `index` of `index_2_name(ids)`
"""
function identifier_name(@nospecialize(ids::IDS); error_on_missing::Bool=true)
    index = identifier_index(ids; error_on_missing)
    if index === nothing
        return nothing
    end
    name = get(index_2_name(ids), index, nothing)
    if name === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids))). Possible options are $(collect(values(index_2_name(ids))))")
    end
    return name
end

function identifier_name(@nospecialize(ids::IDS), default::Symbol)
    name = identifier_name(ids; error_on_missing=false)
    if name === nothing
        return default
    else
        return name
    end
end

"""
    name_2_index(@nospecialize(ids::Union{IDS,IDSvector}))

Return dict of name to IMAS indentifier.index
"""
function name_2_index(@nospecialize(ids::Union{IDS,IDSvector}))
    return Dict(v => k for (k, v) in index_2_name(ids))
end

"""
    findfirst(identifier_name::Symbol, @nospecialize(ids::IDSvector))

Return item from IDSvector based on `index` of `index_2_name(ids)`
"""
function Base.findfirst(identifier_name::Symbol, @nospecialize(ids::IDSvector))
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids))). Possible options are $(collect(values(index_2_name(ids))))")
    elseif :grid_type in fieldnames(eltype(ids))
        index = findfirst(idx -> idx.grid_type.index == i, ids)
    elseif :type in fieldnames(eltype(ids))
        index = findfirst(idx -> idx.type.index == i, ids)
    elseif :function in fieldnames(eltype(ids))
        index = findfirst(idx -> any(func.index == i for func in idx.function), ids)
    elseif :geometry_type in fieldnames(eltype(ids))
        index = findfirst(idx -> idx.geometry_type == i, ids)
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
    findall(identifier_name::Symbol, @nospecialize(ids::IDSvector))

Return items from IDSvector based on `index` of `index_2_name(ids)`
"""
function Base.findall(identifier_name::Symbol, @nospecialize(ids::IDSvector))
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids))). Possible options are $(collect(values(index_2_name(ids))))")
    elseif :grid_type in fieldnames(eltype(ids))
        indexes = findall(idx -> idx.grid_type.index == i, ids)
    elseif :type in fieldnames(eltype(ids))
        indexes = findall(idx -> idx.type.index == i, ids)
    elseif :function in fieldnames(eltype(ids))
        indexes = findall(idx -> any(func.index == i for func in idx.function), ids)
    elseif :geometry_type in fieldnames(eltype(ids))
        indexes = findall(idx -> idx.geometry_type == i, ids)
    elseif :identifier in fieldnames(eltype(ids))
        indexes = findall(idx -> idx.identifier.index == i, ids)
    else
        indexes = findall(idx -> idx.index == i, ids)
    end
    return (ids[index] for index in indexes)
end

"""
    resize!(
        @nospecialize(ids::IDSvector{T}),
        identifier_name::Symbol,
        conditions::Pair{String}...;
        wipe::Bool=true,
        error_multiple_matches::Bool=true
    )::T where {T<:IDSvectorElement}

Resize ids if `identifier_name` is not found based on `index` of `index_2_name(ids)` and a set of conditions are not met.

If wipe=true and an entry matching the condition is found, then the content of the matching IDS is emptied.

Either way, the IDS is populated with the conditions.

NOTE: `error_multiple_matches` will delete all extra entries matching the conditions.

Returns the selected IDS
"""
function Base.resize!(
    @nospecialize(ids::IDSvector{T}),
    identifier_name::Symbol,
    conditions::Pair{String}...;
    wipe::Bool=true,
    error_multiple_matches::Bool=true
)::T where {T<:IDSvectorElement}
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids))). Possible options are $(collect(values(index_2_name(ids))))")
    elseif :grid_type in fieldnames(eltype(ids))
        return resize!(ids, "grid_type.index" => i, conditions...; wipe, error_multiple_matches)
    elseif :type in fieldnames(eltype(ids))
        return resize!(ids, "type.index" => i, conditions...; wipe, error_multiple_matches)
    elseif :geometry_type in fieldnames(eltype(ids))
        return resize!(ids, "geometry_type" => i, conditions...; wipe, error_multiple_matches)
    elseif :identifier in fieldnames(eltype(ids))
        return resize!(ids, "identifier.index" => i, conditions...; wipe, error_multiple_matches)
    else
        return resize!(ids, "index" => i, conditions...; wipe, error_multiple_matches)
    end
end

"""
    deleteat!(@nospecialize(ids::T), identifier_name::Symbol, conditions::Pair{String}...)::T where {T<:IDSvector}

Deletes all entries that match based on `index` of `index_2_name(ids)`
"""
function Base.deleteat!(@nospecialize(ids::T), identifier_name::Symbol, conditions::Pair{String}...)::T where {T<:IDSvector}
    i = get(name_2_index(ids), identifier_name, nothing)
    if i === nothing
        error("`$(repr(identifier_name))` is not a known identifier for dd.$(fs2u(eltype(ids))). Possible options are $(collect(values(index_2_name(ids))))")
    elseif :grid_type in fieldnames(eltype(ids))
        return deleteat!(ids, "grid_type.index" => i, conditions...)
    elseif :type in fieldnames(eltype(ids))
        return deleteat!(ids, "type.index" => i, conditions...)
    elseif :geometry_type in fieldnames(eltype(ids))
        return deleteat!(ids, "geometry_type" => i, conditions...)
    elseif :identifier in fieldnames(eltype(ids))
        return deleteat!(ids, "identifier.index" => i, conditions...)
    else
        return deleteat!(ids, "index" => i, conditions...)
    end
end
