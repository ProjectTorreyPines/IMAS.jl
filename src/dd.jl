using StructArrays:StructArray

Base.@kwdef mutable struct waves__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__magnetic_axis
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__e_field_n_tor__plus
    phase :: Array{Float64, 2} = zeros(Float64,(0,0))
    amplitude :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__e_field_n_tor__plus
    phase :: Array{Float64, 1} = zeros(Float64,(0))
    amplitude :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__e_field__normal
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__e_field__parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__b_field__parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__e_field_n_tor__minus
    phase :: Array{Float64, 2} = zeros(Float64,(0,0))
    amplitude :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{waves__coherent_wave__full_wave__grid__space__objects_per_dimension__object__boundary} = StructArray(waves__coherent_wave__full_wave__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__space__objects_per_dimension
    object :: StructArray{waves__coherent_wave__full_wave__grid__space__objects_per_dimension__object} = StructArray(waves__coherent_wave__full_wave__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__global_quantities__electrons
    power_fast_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    power_fast :: Float64 = 0.0
    power_thermal :: Float64 = 0.0
    power_thermal_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    distribution_assumption :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__e_field__bi_normal
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__e_field__parallel
    imaginary :: Array{Float64, 1} = zeros(Float64,(0))
    real :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__e_field_n_tor__parallel
    phase :: Array{Float64, 1} = zeros(Float64,(0))
    amplitude :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__wave_solver_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__global_quantities__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__e_field__plus
    imaginary :: Array{Float64, 1} = zeros(Float64,(0))
    real :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__e_field__minus
    real :: Array{Float64, 1} = zeros(Float64,(0))
    imaginary :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__power_flow_norm
    perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__identifier__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__identifier
    antenna_name :: String = ""
    index_in_antenna :: Int32 = 0
    type :: waves__coherent_wave__identifier__type = waves__coherent_wave__identifier__type()
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__electrons
    power_density_fast_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    power_density_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_thermal_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    power_density_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    vibrational_level :: Float64 = 0.0
    power_density_fast_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    vibrational_mode :: String = ""
    power_density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_thermal_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__e_field_n_tor__parallel
    phase :: Array{Float64, 2} = zeros(Float64,(0,0))
    amplitude :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__e_field_n_tor
    minus :: waves__coherent_wave__profiles_2d__e_field_n_tor__minus = waves__coherent_wave__profiles_2d__e_field_n_tor__minus()
    plus :: waves__coherent_wave__profiles_2d__e_field_n_tor__plus = waves__coherent_wave__profiles_2d__e_field_n_tor__plus()
    parallel :: waves__coherent_wave__profiles_2d__e_field_n_tor__parallel = waves__coherent_wave__profiles_2d__e_field_n_tor__parallel()
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__electrons
    power :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__k_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__ids_properties
    provider :: String = ""
    version_put :: waves__ids_properties__version_put = waves__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__e_field_n_tor__minus
    phase :: Array{Float64, 1} = zeros(Float64,(0))
    amplitude :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__e_field_n_tor
    minus :: waves__coherent_wave__profiles_1d__e_field_n_tor__minus = waves__coherent_wave__profiles_1d__e_field_n_tor__minus()
    plus :: waves__coherent_wave__profiles_1d__e_field_n_tor__plus = waves__coherent_wave__profiles_1d__e_field_n_tor__plus()
    parallel :: waves__coherent_wave__profiles_1d__e_field_n_tor__parallel = waves__coherent_wave__profiles_1d__e_field_n_tor__parallel()
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__grid_subset__element
    object :: StructArray{waves__coherent_wave__full_wave__grid__grid_subset__element__object} = StructArray(waves__coherent_wave__full_wave__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__e_field
    minus :: waves__coherent_wave__beam_tracing__beam__e_field__minus = waves__coherent_wave__beam_tracing__beam__e_field__minus()
    plus :: waves__coherent_wave__beam_tracing__beam__e_field__plus = waves__coherent_wave__beam_tracing__beam__e_field__plus()
    parallel :: waves__coherent_wave__beam_tracing__beam__e_field__parallel = waves__coherent_wave__beam_tracing__beam__e_field__parallel()
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__b_field__normal
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__position
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    theta :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__e_field__plus
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__grid_subset
    base :: StructArray{waves__coherent_wave__full_wave__grid__grid_subset__base} = StructArray(waves__coherent_wave__full_wave__grid__grid_subset__base() for k in 1:1)
    metric :: waves__coherent_wave__full_wave__grid__grid_subset__metric = waves__coherent_wave__full_wave__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: waves__coherent_wave__full_wave__grid__grid_subset__identifier = waves__coherent_wave__full_wave__grid__grid_subset__identifier()
    element :: StructArray{waves__coherent_wave__full_wave__grid__grid_subset__element} = StructArray(waves__coherent_wave__full_wave__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__e_field__minus
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__e_field
    parallel :: StructArray{waves__coherent_wave__full_wave__e_field__parallel} = StructArray(waves__coherent_wave__full_wave__e_field__parallel() for k in 1:1)
    bi_normal :: StructArray{waves__coherent_wave__full_wave__e_field__bi_normal} = StructArray(waves__coherent_wave__full_wave__e_field__bi_normal() for k in 1:1)
    minus :: StructArray{waves__coherent_wave__full_wave__e_field__minus} = StructArray(waves__coherent_wave__full_wave__e_field__minus() for k in 1:1)
    plus :: StructArray{waves__coherent_wave__full_wave__e_field__plus} = StructArray(waves__coherent_wave__full_wave__e_field__plus() for k in 1:1)
    normal :: StructArray{waves__coherent_wave__full_wave__e_field__normal} = StructArray(waves__coherent_wave__full_wave__e_field__normal() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    geometry_type :: waves__coherent_wave__full_wave__grid__space__geometry_type = waves__coherent_wave__full_wave__grid__space__geometry_type()
    identifier :: waves__coherent_wave__full_wave__grid__space__identifier = waves__coherent_wave__full_wave__grid__space__identifier()
    objects_per_dimension :: StructArray{waves__coherent_wave__full_wave__grid__space__objects_per_dimension} = StructArray(waves__coherent_wave__full_wave__grid__space__objects_per_dimension() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__grid
    grid_subset :: StructArray{waves__coherent_wave__full_wave__grid__grid_subset} = StructArray(waves__coherent_wave__full_wave__grid__grid_subset() for k in 1:1)
    space :: StructArray{waves__coherent_wave__full_wave__grid__space} = StructArray(waves__coherent_wave__full_wave__grid__space() for k in 1:1)
    identifier :: waves__coherent_wave__full_wave__grid__identifier = waves__coherent_wave__full_wave__grid__identifier()
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__ion__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    power_inside_thermal_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    electron_configuration :: String = ""
    vibrational_mode :: String = ""
    power_inside_fast_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_max :: Float64 = 0.0
    power_inside_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_fast_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_inside_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_thermal_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__ion
    label :: String = ""
    power_inside_thermal_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    state :: StructArray{waves__coherent_wave__profiles_1d__ion__state} = StructArray(waves__coherent_wave__profiles_1d__ion__state() for k in 1:1)
    power_inside_fast_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_ion :: Float64 = 0.0
    power_inside_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_fast_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    power_inside_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_thermal_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    element :: StructArray{waves__coherent_wave__profiles_1d__ion__element} = StructArray(waves__coherent_wave__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__electrons
    power_inside_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_fast_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    power_inside_thermal_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_inside_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_thermal_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_inside_fast_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__wave_vector
    k_tor :: Array{Float64, 1} = zeros(Float64,(0))
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    k_z :: Array{Float64, 1} = zeros(Float64,(0))
    n_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    k_r :: Array{Float64, 1} = zeros(Float64,(0))
    n_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    varying_n_tor :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    power_density_fast_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    power_density_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    multiple_states_flag :: Int32 = 0
    power_density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    state :: StructArray{waves__coherent_wave__profiles_2d__ion__state} = StructArray(waves__coherent_wave__profiles_2d__ion__state() for k in 1:1)
    power_density_thermal_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    element :: StructArray{waves__coherent_wave__profiles_2d__ion__element} = StructArray(waves__coherent_wave__profiles_2d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__global_quantities__ion__state
    label :: String = ""
    power_fast_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    electron_configuration :: String = ""
    power_thermal :: Float64 = 0.0
    vibrational_level :: Float64 = 0.0
    power_fast :: Float64 = 0.0
    vibrational_mode :: String = ""
    power_thermal_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__global_quantities__ion
    power_fast_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    label :: String = ""
    power_fast :: Float64 = 0.0
    power_thermal :: Float64 = 0.0
    z_ion :: Float64 = 0.0
    power_thermal_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    state :: StructArray{waves__coherent_wave__global_quantities__ion__state} = StructArray(waves__coherent_wave__global_quantities__ion__state() for k in 1:1)
    distribution_assumption :: Int32 = 0
    element :: StructArray{waves__coherent_wave__global_quantities__ion__element} = StructArray(waves__coherent_wave__global_quantities__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__global_quantities
    ion :: StructArray{waves__coherent_wave__global_quantities__ion} = StructArray(waves__coherent_wave__global_quantities__ion() for k in 1:1)
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    time :: Float64 = 0.0
    current_tor_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    frequency :: Float64 = 0.0
    current_tor :: Float64 = 0.0
    power :: Float64 = 0.0
    power_n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: waves__coherent_wave__global_quantities__electrons = waves__coherent_wave__global_quantities__electrons()
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__b_field__bi_normal
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave__b_field
    parallel :: StructArray{waves__coherent_wave__full_wave__b_field__parallel} = StructArray(waves__coherent_wave__full_wave__b_field__parallel() for k in 1:1)
    bi_normal :: StructArray{waves__coherent_wave__full_wave__b_field__bi_normal} = StructArray(waves__coherent_wave__full_wave__b_field__bi_normal() for k in 1:1)
    normal :: StructArray{waves__coherent_wave__full_wave__b_field__normal} = StructArray(waves__coherent_wave__full_wave__b_field__normal() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__full_wave
    time :: Float64 = 0.0
    b_field :: waves__coherent_wave__full_wave__b_field = waves__coherent_wave__full_wave__b_field()
    grid :: waves__coherent_wave__full_wave__grid = waves__coherent_wave__full_wave__grid()
    e_field :: waves__coherent_wave__full_wave__e_field = waves__coherent_wave__full_wave__e_field()
    k_perpendicular :: StructArray{waves__coherent_wave__full_wave__k_perpendicular} = StructArray(waves__coherent_wave__full_wave__k_perpendicular() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__code
    library :: StructArray{waves__code__library} = StructArray(waves__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    power :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    multiple_states_flag :: Int32 = 0
    state :: StructArray{waves__coherent_wave__beam_tracing__beam__ion__state} = StructArray(waves__coherent_wave__beam_tracing__beam__ion__state() for k in 1:1)
    power :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{waves__coherent_wave__beam_tracing__beam__ion__element} = StructArray(waves__coherent_wave__beam_tracing__beam__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing__beam
    ion :: StructArray{waves__coherent_wave__beam_tracing__beam__ion} = StructArray(waves__coherent_wave__beam_tracing__beam__ion() for k in 1:1)
    length :: Array{Float64, 1} = zeros(Float64,(0))
    power_initial :: Float64 = 0.0
    wave_vector :: waves__coherent_wave__beam_tracing__beam__wave_vector = waves__coherent_wave__beam_tracing__beam__wave_vector()
    position :: waves__coherent_wave__beam_tracing__beam__position = waves__coherent_wave__beam_tracing__beam__position()
    e_field :: waves__coherent_wave__beam_tracing__beam__e_field = waves__coherent_wave__beam_tracing__beam__e_field()
    electrons :: waves__coherent_wave__beam_tracing__beam__electrons = waves__coherent_wave__beam_tracing__beam__electrons()
    power_flow_norm :: waves__coherent_wave__beam_tracing__beam__power_flow_norm = waves__coherent_wave__beam_tracing__beam__power_flow_norm()
end

Base.@kwdef mutable struct waves__coherent_wave__beam_tracing
    beam :: StructArray{waves__coherent_wave__beam_tracing__beam} = StructArray(waves__coherent_wave__beam_tracing__beam() for k in 1:1)
    time :: Float64 = 0.0
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__grid__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d__grid
    psi :: Array{Float64, 2} = zeros(Float64,(0,0))
    volume :: Array{Float64, 2} = zeros(Float64,(0,0))
    area :: Array{Float64, 2} = zeros(Float64,(0,0))
    theta_straight :: Array{Float64, 2} = zeros(Float64,(0,0))
    theta_geometric :: Array{Float64, 2} = zeros(Float64,(0,0))
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    rho_tor_norm :: Array{Float64, 2} = zeros(Float64,(0,0))
    type :: waves__coherent_wave__profiles_2d__grid__type = waves__coherent_wave__profiles_2d__grid__type()
    rho_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_2d
    ion :: StructArray{waves__coherent_wave__profiles_2d__ion} = StructArray(waves__coherent_wave__profiles_2d__ion() for k in 1:1)
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    time :: Float64 = 0.0
    e_field_n_tor :: StructArray{waves__coherent_wave__profiles_2d__e_field_n_tor} = StructArray(waves__coherent_wave__profiles_2d__e_field_n_tor() for k in 1:1)
    power_density :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid :: waves__coherent_wave__profiles_2d__grid = waves__coherent_wave__profiles_2d__grid()
    power_density_n_tor :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    electrons :: waves__coherent_wave__profiles_2d__electrons = waves__coherent_wave__profiles_2d__electrons()
end

Base.@kwdef mutable struct waves__coherent_wave__profiles_1d
    time :: Float64 = 0.0
    e_field_n_tor :: StructArray{waves__coherent_wave__profiles_1d__e_field_n_tor} = StructArray(waves__coherent_wave__profiles_1d__e_field_n_tor() for k in 1:1)
    ion :: StructArray{waves__coherent_wave__profiles_1d__ion} = StructArray(waves__coherent_wave__profiles_1d__ion() for k in 1:1)
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    current_tor_inside_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_tor_inside :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    electrons :: waves__coherent_wave__profiles_1d__electrons = waves__coherent_wave__profiles_1d__electrons()
    current_parallel_density :: Array{Float64, 1} = zeros(Float64,(0))
    power_inside_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid :: waves__coherent_wave__profiles_1d__grid = waves__coherent_wave__profiles_1d__grid()
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
    power_density :: Array{Float64, 1} = zeros(Float64,(0))
    k_perpendicular :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_parallel_density_n_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct waves__coherent_wave
    global_quantities :: StructArray{waves__coherent_wave__global_quantities} = StructArray(waves__coherent_wave__global_quantities() for k in 1:1)
    beam_tracing :: StructArray{waves__coherent_wave__beam_tracing} = StructArray(waves__coherent_wave__beam_tracing() for k in 1:1)
    full_wave :: StructArray{waves__coherent_wave__full_wave} = StructArray(waves__coherent_wave__full_wave() for k in 1:1)
    profiles_1d :: StructArray{waves__coherent_wave__profiles_1d} = StructArray(waves__coherent_wave__profiles_1d() for k in 1:1)
    profiles_2d :: StructArray{waves__coherent_wave__profiles_2d} = StructArray(waves__coherent_wave__profiles_2d() for k in 1:1)
    identifier :: waves__coherent_wave__identifier = waves__coherent_wave__identifier()
    wave_solver_type :: waves__coherent_wave__wave_solver_type = waves__coherent_wave__wave_solver_type()
end

Base.@kwdef mutable struct waves
    coherent_wave :: StructArray{waves__coherent_wave} = StructArray(waves__coherent_wave() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: waves__ids_properties = waves__ids_properties()
    vacuum_toroidal_field :: waves__vacuum_toroidal_field = waves__vacuum_toroidal_field()
    code :: waves__code = waves__code()
    magnetic_axis :: waves__magnetic_axis = waves__magnetic_axis()
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__annular__centreline
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_2d__mobile__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__annular__outline_inner
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__grid_subset__element
    object :: StructArray{wall__description_ggd__grid_ggd__grid_subset__element__object} = StructArray(wall__description_ggd__grid_ggd__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__annular__outline_outer
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_ggd__ggd__power_density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct wall__description_2d__vessel__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct wall__description_ggd__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__temperature_reference
    data :: Float64 = 0.0
    description :: String = ""
end

Base.@kwdef mutable struct wall__description_2d__limiter__unit__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_2d__limiter__unit
    name :: String = ""
    resistivity :: Float64 = 0.0
    closed :: Int32 = 0
    outline :: wall__description_2d__limiter__unit__outline = wall__description_2d__limiter__unit__outline()
    phi_extensions :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct wall__global_quantities__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__global_quantities__electrons
    pumping_speed :: Array{Float64, 1} = zeros(Float64,(0))
    gas_puff :: Array{Float64, 1} = zeros(Float64,(0))
    particle_flux_from_plasma :: Array{Float64, 1} = zeros(Float64,(0))
    power_outer_target :: Array{Float64, 1} = zeros(Float64,(0))
    particle_flux_from_wall :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_inner_target :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_ggd__ggd__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct wall__description_ggd__ggd
    temperature :: StructArray{wall__description_ggd__ggd__temperature} = StructArray(wall__description_ggd__ggd__temperature() for k in 1:1)
    time :: Float64 = 0.0
    power_density :: StructArray{wall__description_ggd__ggd__power_density} = StructArray(wall__description_ggd__ggd__power_density() for k in 1:1)
end

Base.@kwdef mutable struct wall__first_wall_power_flux_peak
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__ids_properties
    provider :: String = ""
    version_put :: wall__ids_properties__version_put = wall__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_2d__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__element__outline
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__element__j_tor
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__element
    name :: String = ""
    j_tor :: wall__description_2d__vessel__unit__element__j_tor = wall__description_2d__vessel__unit__element__j_tor()
    resistivity :: Float64 = 0.0
    outline :: wall__description_2d__vessel__unit__element__outline = wall__description_2d__vessel__unit__element__outline()
    resistance :: Float64 = 0.0
end

Base.@kwdef mutable struct wall__global_quantities__neutral
    label :: String = ""
    sputtering_chemical_coefficient :: Array{Float64, 2} = zeros(Float64,(0,0))
    gas_puff :: Array{Float64, 1} = zeros(Float64,(0))
    recycling_particles_coefficient :: Array{Float64, 2} = zeros(Float64,(0,0))
    pumping_speed :: Array{Float64, 1} = zeros(Float64,(0))
    particle_flux_from_wall :: Array{Float64, 2} = zeros(Float64,(0,0))
    recycling_energy_coefficient :: Array{Float64, 2} = zeros(Float64,(0,0))
    wall_inventory :: Array{Float64, 1} = zeros(Float64,(0))
    particle_flux_from_plasma :: Array{Float64, 1} = zeros(Float64,(0))
    sputtering_physical_coefficient :: Array{Float64, 2} = zeros(Float64,(0,0))
    element :: StructArray{wall__global_quantities__neutral__element} = StructArray(wall__global_quantities__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct wall__global_quantities
    neutral :: StructArray{wall__global_quantities__neutral} = StructArray(wall__global_quantities__neutral() for k in 1:1)
    power_incident :: Array{Float64, 1} = zeros(Float64,(0))
    power_radiated :: Array{Float64, 1} = zeros(Float64,(0))
    power_inner_target_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    power_conducted :: Array{Float64, 1} = zeros(Float64,(0))
    power_convected :: Array{Float64, 1} = zeros(Float64,(0))
    current_tor :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: wall__global_quantities__electrons = wall__global_quantities__electrons()
    power_density_inner_target_max :: Array{Float64, 1} = zeros(Float64,(0))
    power_black_body :: Array{Float64, 1} = zeros(Float64,(0))
    power_recombination_neutrals :: Array{Float64, 1} = zeros(Float64,(0))
    power_to_cooling :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_outer_target_max :: Array{Float64, 1} = zeros(Float64,(0))
    power_recombination_plasma :: Array{Float64, 1} = zeros(Float64,(0))
    power_currents :: Array{Float64, 1} = zeros(Float64,(0))
    power_neutrals :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit__annular
    outline_inner :: wall__description_2d__vessel__unit__annular__outline_inner = wall__description_2d__vessel__unit__annular__outline_inner()
    centreline :: wall__description_2d__vessel__unit__annular__centreline = wall__description_2d__vessel__unit__annular__centreline()
    thickness :: Array{Float64, 1} = zeros(Float64,(0))
    resistivity :: Float64 = 0.0
    outline_outer :: wall__description_2d__vessel__unit__annular__outline_outer = wall__description_2d__vessel__unit__annular__outline_outer()
end

Base.@kwdef mutable struct wall__description_2d__vessel__unit
    name :: String = ""
    annular :: wall__description_2d__vessel__unit__annular = wall__description_2d__vessel__unit__annular()
    identifier :: String = ""
    element :: StructArray{wall__description_2d__vessel__unit__element} = StructArray(wall__description_2d__vessel__unit__element() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{wall__description_ggd__grid_ggd__space__objects_per_dimension__object__boundary} = StructArray(wall__description_ggd__grid_ggd__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__space__objects_per_dimension
    object :: StructArray{wall__description_ggd__grid_ggd__space__objects_per_dimension__object} = StructArray(wall__description_ggd__grid_ggd__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct wall__code
    library :: StructArray{wall__code__library} = StructArray(wall__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct wall__description_2d__mobile__unit__outline
    time :: Float64 = 0.0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct wall__description_2d__mobile__unit
    name :: String = ""
    resistivity :: Float64 = 0.0
    closed :: Int32 = 0
    phi_extensions :: Array{Float64, 2} = zeros(Float64,(0,0))
    outline :: StructArray{wall__description_2d__mobile__unit__outline} = StructArray(wall__description_2d__mobile__unit__outline() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_2d__mobile
    unit :: StructArray{wall__description_2d__mobile__unit} = StructArray(wall__description_2d__mobile__unit() for k in 1:1)
    type :: wall__description_2d__mobile__type = wall__description_2d__mobile__type()
end

Base.@kwdef mutable struct wall__description_2d__vessel
    unit :: StructArray{wall__description_2d__vessel__unit} = StructArray(wall__description_2d__vessel__unit() for k in 1:1)
    type :: wall__description_2d__vessel__type = wall__description_2d__vessel__type()
end

Base.@kwdef mutable struct wall__description_2d__limiter__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct wall__description_2d__limiter
    type :: wall__description_2d__limiter__type = wall__description_2d__limiter__type()
    unit :: StructArray{wall__description_2d__limiter__unit} = StructArray(wall__description_2d__limiter__unit() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_2d
    mobile :: wall__description_2d__mobile = wall__description_2d__mobile()
    limiter :: wall__description_2d__limiter = wall__description_2d__limiter()
    type :: wall__description_2d__type = wall__description_2d__type()
    vessel :: wall__description_2d__vessel = wall__description_2d__vessel()
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__grid_subset
    base :: StructArray{wall__description_ggd__grid_ggd__grid_subset__base} = StructArray(wall__description_ggd__grid_ggd__grid_subset__base() for k in 1:1)
    metric :: wall__description_ggd__grid_ggd__grid_subset__metric = wall__description_ggd__grid_ggd__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: wall__description_ggd__grid_ggd__grid_subset__identifier = wall__description_ggd__grid_ggd__grid_subset__identifier()
    element :: StructArray{wall__description_ggd__grid_ggd__grid_subset__element} = StructArray(wall__description_ggd__grid_ggd__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{wall__description_ggd__grid_ggd__space__objects_per_dimension} = StructArray(wall__description_ggd__grid_ggd__space__objects_per_dimension() for k in 1:1)
    geometry_type :: wall__description_ggd__grid_ggd__space__geometry_type = wall__description_ggd__grid_ggd__space__geometry_type()
    identifier :: wall__description_ggd__grid_ggd__space__identifier = wall__description_ggd__grid_ggd__space__identifier()
end

Base.@kwdef mutable struct wall__description_ggd__grid_ggd
    time :: Float64 = 0.0
    grid_subset :: StructArray{wall__description_ggd__grid_ggd__grid_subset} = StructArray(wall__description_ggd__grid_ggd__grid_subset() for k in 1:1)
    space :: StructArray{wall__description_ggd__grid_ggd__space} = StructArray(wall__description_ggd__grid_ggd__space() for k in 1:1)
    identifier :: wall__description_ggd__grid_ggd__identifier = wall__description_ggd__grid_ggd__identifier()
end

Base.@kwdef mutable struct wall__description_ggd
    grid_ggd :: StructArray{wall__description_ggd__grid_ggd} = StructArray(wall__description_ggd__grid_ggd() for k in 1:1)
    type :: wall__description_ggd__type = wall__description_ggd__type()
    ggd :: StructArray{wall__description_ggd__ggd} = StructArray(wall__description_ggd__ggd() for k in 1:1)
end

Base.@kwdef mutable struct wall
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: wall__ids_properties = wall__ids_properties()
    description_ggd :: StructArray{wall__description_ggd} = StructArray(wall__description_ggd() for k in 1:1)
    description_2d :: StructArray{wall__description_2d} = StructArray(wall__description_2d() for k in 1:1)
    first_wall_surface_area :: Float64 = 0.0
    code :: wall__code = wall__code()
    global_quantities :: wall__global_quantities = wall__global_quantities()
    temperature_reference :: wall__temperature_reference = wall__temperature_reference()
    first_wall_power_flux_peak :: wall__first_wall_power_flux_peak = wall__first_wall_power_flux_peak()
end

Base.@kwdef mutable struct turbulence__grid_2d_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct turbulence__profiles_2d__electrons
    density_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    temperature :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct turbulence__grid_2d
    time :: Float64 = 0.0
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct turbulence__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct turbulence__profiles_2d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct turbulence__profiles_2d__neutral
    label :: String = ""
    ion_index :: Int32 = 0
    density_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    temperature :: Array{Float64, 2} = zeros(Float64,(0,0))
    element :: StructArray{turbulence__profiles_2d__neutral__element} = StructArray(turbulence__profiles_2d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct turbulence__profiles_2d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct turbulence__profiles_2d__ion
    z_ion :: Float64 = 0.0
    label :: String = ""
    density_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    neutral_index :: Int32 = 0
    temperature :: Array{Float64, 2} = zeros(Float64,(0,0))
    element :: StructArray{turbulence__profiles_2d__ion__element} = StructArray(turbulence__profiles_2d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct turbulence__profiles_2d
    ion :: StructArray{turbulence__profiles_2d__ion} = StructArray(turbulence__profiles_2d__ion() for k in 1:1)
    time :: Float64 = 0.0
    neutral :: StructArray{turbulence__profiles_2d__neutral} = StructArray(turbulence__profiles_2d__neutral() for k in 1:1)
    electrons :: turbulence__profiles_2d__electrons = turbulence__profiles_2d__electrons()
end

Base.@kwdef mutable struct turbulence__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct turbulence__ids_properties
    provider :: String = ""
    version_put :: turbulence__ids_properties__version_put = turbulence__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct turbulence__code
    library :: StructArray{turbulence__code__library} = StructArray(turbulence__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct turbulence
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: turbulence__ids_properties = turbulence__ids_properties()
    grid_2d_type :: turbulence__grid_2d_type = turbulence__grid_2d_type()
    grid_2d :: StructArray{turbulence__grid_2d} = StructArray(turbulence__grid_2d() for k in 1:1)
    code :: turbulence__code = turbulence__code()
    profiles_2d :: StructArray{turbulence__profiles_2d} = StructArray(turbulence__profiles_2d() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__control_parameters__real0d
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__d_dt
    pressure_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__energy__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__boundary_condition__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__energy__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__time_step_min
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__energy
    delta_relative :: transport_solver_numerics__convergence__equations__ion__energy__delta_relative = transport_solver_numerics__convergence__equations__ion__energy__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__element
    object :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__element__object} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__energy__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__control_parameters__integer0d
    value :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__particles__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__time_step
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__electrons__energy__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__d_drho_tor_norm
    pressure_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__state__energy__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__state__energy
    delta_relative :: transport_solver_numerics__convergence__equations__ion__state__energy__delta_relative = transport_solver_numerics__convergence__equations__ion__state__energy__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__particles__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__momentum_tor__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__ion__energy__identifier = transport_solver_numerics__boundary_conditions_ggd__ion__energy__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__time_step_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__particles
    value :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Float64 = 0.0
    identifier :: transport_solver_numerics__boundary_conditions_1d__ion__particles__identifier = transport_solver_numerics__boundary_conditions_1d__ion__particles__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__electrons__particles__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__electrons__particles
    delta_relative :: transport_solver_numerics__convergence__equations__electrons__particles__delta_relative = transport_solver_numerics__convergence__equations__electrons__particles__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__state__energy__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__computation_mode
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__primary_quantity__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__current__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__current
    delta_relative :: transport_solver_numerics__convergence__equations__current__delta_relative = transport_solver_numerics__convergence__equations__current__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__d2_drho_tor_norm2
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__electrons__energy__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__electrons__particles__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__electrons__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__electrons__particles__identifier = transport_solver_numerics__boundary_conditions_ggd__electrons__particles__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__d2_drho_tor_norm2
    pressure_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__convergence__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__convergence
    delta_relative :: transport_solver_numerics__solver_1d__equation__convergence__delta_relative = transport_solver_numerics__solver_1d__equation__convergence__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__electrons__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__electrons__energy__identifier = transport_solver_numerics__boundary_conditions_ggd__electrons__energy__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__primary_coordinate
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__restart_files
    names :: Array{String, 1} = String[]
    descriptions :: Array{String, 1} = String[]
    time :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__state__energy
    value :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Float64 = 0.0
    identifier :: transport_solver_numerics__boundary_conditions_1d__ion__state__energy__identifier = transport_solver_numerics__boundary_conditions_1d__ion__state__energy__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__coefficient
    profile :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__electrons
    particles :: StructArray{transport_solver_numerics__boundary_conditions_ggd__electrons__particles} = StructArray(transport_solver_numerics__boundary_conditions_ggd__electrons__particles() for k in 1:1)
    energy :: StructArray{transport_solver_numerics__boundary_conditions_ggd__electrons__energy} = StructArray(transport_solver_numerics__boundary_conditions_ggd__electrons__energy() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__state__particles__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__state__particles
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: transport_solver_numerics__boundary_conditions_1d__ion__state__particles__identifier = transport_solver_numerics__boundary_conditions_1d__ion__state__particles__identifier()
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__electrons__particles__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__electrons__particles
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: transport_solver_numerics__boundary_conditions_1d__electrons__particles__identifier = transport_solver_numerics__boundary_conditions_1d__electrons__particles__identifier()
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__state__energy__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__state__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__ion__state__energy__identifier = transport_solver_numerics__boundary_conditions_ggd__ion__state__energy__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__boundary_condition
    position :: Float64 = 0.0
    value :: Array{Float64, 1} = zeros(Float64,(0))
    type :: transport_solver_numerics__solver_1d__equation__boundary_condition__type = transport_solver_numerics__solver_1d__equation__boundary_condition__type()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__current__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__state__d_dt
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__electrons__energy__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__electrons__energy
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: transport_solver_numerics__boundary_conditions_1d__electrons__energy__identifier = transport_solver_numerics__boundary_conditions_1d__electrons__energy__identifier()
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__electrons
    particles :: transport_solver_numerics__boundary_conditions_1d__electrons__particles = transport_solver_numerics__boundary_conditions_1d__electrons__particles()
    energy :: transport_solver_numerics__boundary_conditions_1d__electrons__energy = transport_solver_numerics__boundary_conditions_1d__electrons__energy()
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__electrons__d2_drho_tor_norm2
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__state__d_drho_tor_norm
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__state__d2_drho_tor_norm2
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__control_parameters
    integer0d :: StructArray{transport_solver_numerics__solver_1d__control_parameters__integer0d} = StructArray(transport_solver_numerics__solver_1d__control_parameters__integer0d() for k in 1:1)
    real0d :: StructArray{transport_solver_numerics__solver_1d__control_parameters__real0d} = StructArray(transport_solver_numerics__solver_1d__control_parameters__real0d() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__electrons__energy
    delta_relative :: transport_solver_numerics__convergence__equations__electrons__energy__delta_relative = transport_solver_numerics__convergence__equations__electrons__energy__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__electrons
    particles :: transport_solver_numerics__convergence__equations__electrons__particles = transport_solver_numerics__convergence__equations__electrons__particles()
    energy :: transport_solver_numerics__convergence__equations__electrons__energy = transport_solver_numerics__convergence__equations__electrons__energy()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__state__particles__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__state__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__ion__state__particles__identifier = transport_solver_numerics__boundary_conditions_ggd__ion__state__particles__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__state
    label :: String = ""
    particles :: StructArray{transport_solver_numerics__boundary_conditions_ggd__ion__state__particles} = StructArray(transport_solver_numerics__boundary_conditions_ggd__ion__state__particles() for k in 1:1)
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    energy :: StructArray{transport_solver_numerics__boundary_conditions_ggd__ion__state__energy} = StructArray(transport_solver_numerics__boundary_conditions_ggd__ion__state__energy() for k in 1:1)
    is_neutral :: Int32 = 0
    neutral_type :: transport_solver_numerics__boundary_conditions_ggd__ion__state__neutral_type = transport_solver_numerics__boundary_conditions_ggd__ion__state__neutral_type()
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation__primary_quantity
    profile :: Array{Float64, 1} = zeros(Float64,(0))
    d_dt_cr :: Array{Float64, 1} = zeros(Float64,(0))
    d2_dr2 :: Array{Float64, 1} = zeros(Float64,(0))
    d_dt_cphi :: Array{Float64, 1} = zeros(Float64,(0))
    ion_index :: Int32 = 0
    state_index :: Int32 = 0
    identifier :: transport_solver_numerics__solver_1d__equation__primary_quantity__identifier = transport_solver_numerics__solver_1d__equation__primary_quantity__identifier()
    d_dr :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    d_dt :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__electrons__d_drho_tor_norm
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__particles__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__particles
    delta_relative :: transport_solver_numerics__convergence__equations__ion__particles__delta_relative = transport_solver_numerics__convergence__equations__ion__particles__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__energy
    value :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Float64 = 0.0
    identifier :: transport_solver_numerics__boundary_conditions_1d__ion__energy__identifier = transport_solver_numerics__boundary_conditions_1d__ion__energy__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__state__particles__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__state__particles
    delta_relative :: transport_solver_numerics__convergence__equations__ion__state__particles__delta_relative = transport_solver_numerics__convergence__equations__ion__state__particles__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__ion__particles__identifier = transport_solver_numerics__boundary_conditions_ggd__ion__particles__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__current__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__current
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: transport_solver_numerics__boundary_conditions_1d__current__identifier = transport_solver_numerics__boundary_conditions_1d__current__identifier()
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__time_step
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__energy_ion_total__delta_relative
    expression :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__energy_ion_total
    delta_relative :: transport_solver_numerics__convergence__equations__energy_ion_total__delta_relative = transport_solver_numerics__convergence__equations__energy_ion_total__delta_relative()
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__ion
    label :: String = ""
    particles :: StructArray{transport_solver_numerics__boundary_conditions_ggd__ion__particles} = StructArray(transport_solver_numerics__boundary_conditions_ggd__ion__particles() for k in 1:1)
    z_n :: Float64 = 0.0
    energy :: StructArray{transport_solver_numerics__boundary_conditions_ggd__ion__energy} = StructArray(transport_solver_numerics__boundary_conditions_ggd__ion__energy() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    state :: StructArray{transport_solver_numerics__boundary_conditions_ggd__ion__state} = StructArray(transport_solver_numerics__boundary_conditions_ggd__ion__state() for k in 1:1)
    a :: Float64 = 0.0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension__object__boundary} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension
    object :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension__object} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__space__objects_per_dimension() for k in 1:1)
    geometry_type :: transport_solver_numerics__boundary_conditions_ggd__grid__space__geometry_type = transport_solver_numerics__boundary_conditions_ggd__grid__space__geometry_type()
    identifier :: transport_solver_numerics__boundary_conditions_ggd__grid__space__identifier = transport_solver_numerics__boundary_conditions_ggd__grid__space__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__energy_ion_total__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__energy_ion_total
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: transport_solver_numerics__boundary_conditions_1d__energy_ion_total__identifier = transport_solver_numerics__boundary_conditions_1d__energy_ion_total__identifier()
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__code
    library :: StructArray{transport_solver_numerics__code__library} = StructArray(transport_solver_numerics__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__momentum_tor
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: transport_solver_numerics__boundary_conditions_1d__momentum_tor__identifier = transport_solver_numerics__boundary_conditions_1d__momentum_tor__identifier()
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    is_neutral :: Int32 = 0
    z_min :: Float64 = 0.0
    d_drho_tor_norm :: transport_solver_numerics__derivatives_1d__ion__state__d_drho_tor_norm = transport_solver_numerics__derivatives_1d__ion__state__d_drho_tor_norm()
    electron_configuration :: String = ""
    d2_drho_tor_norm2 :: transport_solver_numerics__derivatives_1d__ion__state__d2_drho_tor_norm2 = transport_solver_numerics__derivatives_1d__ion__state__d2_drho_tor_norm2()
    vibrational_mode :: String = ""
    z_max :: Float64 = 0.0
    neutral_type :: transport_solver_numerics__derivatives_1d__ion__state__neutral_type = transport_solver_numerics__derivatives_1d__ion__state__neutral_type()
    d_dt :: transport_solver_numerics__derivatives_1d__ion__state__d_dt = transport_solver_numerics__derivatives_1d__ion__state__d_dt()
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    particles :: transport_solver_numerics__convergence__equations__ion__state__particles = transport_solver_numerics__convergence__equations__ion__state__particles()
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    is_neutral :: Int32 = 0
    neutral_type :: transport_solver_numerics__convergence__equations__ion__state__neutral_type = transport_solver_numerics__convergence__equations__ion__state__neutral_type()
    energy :: transport_solver_numerics__convergence__equations__ion__state__energy = transport_solver_numerics__convergence__equations__ion__state__energy()
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations__ion
    label :: String = ""
    z_n :: Float64 = 0.0
    particles :: transport_solver_numerics__convergence__equations__ion__particles = transport_solver_numerics__convergence__equations__ion__particles()
    multiple_states_flag :: Int32 = 0
    state :: StructArray{transport_solver_numerics__convergence__equations__ion__state} = StructArray(transport_solver_numerics__convergence__equations__ion__state() for k in 1:1)
    energy :: transport_solver_numerics__convergence__equations__ion__energy = transport_solver_numerics__convergence__equations__ion__energy()
    a :: Float64 = 0.0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__convergence__equations
    ion :: StructArray{transport_solver_numerics__convergence__equations__ion} = StructArray(transport_solver_numerics__convergence__equations__ion() for k in 1:1)
    current :: transport_solver_numerics__convergence__equations__current = transport_solver_numerics__convergence__equations__current()
    time :: Float64 = 0.0
    electrons :: transport_solver_numerics__convergence__equations__electrons = transport_solver_numerics__convergence__equations__electrons()
    energy_ion_total :: transport_solver_numerics__convergence__equations__energy_ion_total = transport_solver_numerics__convergence__equations__energy_ion_total()
end

Base.@kwdef mutable struct transport_solver_numerics__convergence
    equations :: StructArray{transport_solver_numerics__convergence__equations} = StructArray(transport_solver_numerics__convergence__equations() for k in 1:1)
    time_step :: transport_solver_numerics__convergence__time_step = transport_solver_numerics__convergence__time_step()
end

Base.@kwdef mutable struct transport_solver_numerics__ids_properties
    provider :: String = ""
    version_put :: transport_solver_numerics__ids_properties__version_put = transport_solver_numerics__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d__equation
    boundary_condition :: StructArray{transport_solver_numerics__solver_1d__equation__boundary_condition} = StructArray(transport_solver_numerics__solver_1d__equation__boundary_condition() for k in 1:1)
    convergence :: transport_solver_numerics__solver_1d__equation__convergence = transport_solver_numerics__solver_1d__equation__convergence()
    computation_mode :: transport_solver_numerics__solver_1d__equation__computation_mode = transport_solver_numerics__solver_1d__equation__computation_mode()
    primary_quantity :: transport_solver_numerics__solver_1d__equation__primary_quantity = transport_solver_numerics__solver_1d__equation__primary_quantity()
    coefficient :: StructArray{transport_solver_numerics__solver_1d__equation__coefficient} = StructArray(transport_solver_numerics__solver_1d__equation__coefficient() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__solver_1d
    time :: Float64 = 0.0
    equation :: StructArray{transport_solver_numerics__solver_1d__equation} = StructArray(transport_solver_numerics__solver_1d__equation() for k in 1:1)
    drho_tor_dt :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: transport_solver_numerics__solver_1d__grid = transport_solver_numerics__solver_1d__grid()
    d_dvolume_drho_tor_dt :: Array{Float64, 1} = zeros(Float64,(0))
    control_parameters :: transport_solver_numerics__solver_1d__control_parameters = transport_solver_numerics__solver_1d__control_parameters()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset
    base :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__base} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__base() for k in 1:1)
    metric :: transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__metric = transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__identifier = transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__identifier()
    element :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__element} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__grid
    grid_subset :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__grid_subset() for k in 1:1)
    space :: StructArray{transport_solver_numerics__boundary_conditions_ggd__grid__space} = StructArray(transport_solver_numerics__boundary_conditions_ggd__grid__space() for k in 1:1)
    identifier :: transport_solver_numerics__boundary_conditions_ggd__grid__identifier = transport_solver_numerics__boundary_conditions_ggd__grid__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__d_drho_tor_norm
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    particles :: transport_solver_numerics__boundary_conditions_1d__ion__state__particles = transport_solver_numerics__boundary_conditions_1d__ion__state__particles()
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    is_neutral :: Int32 = 0
    neutral_type :: transport_solver_numerics__boundary_conditions_1d__ion__state__neutral_type = transport_solver_numerics__boundary_conditions_1d__ion__state__neutral_type()
    energy :: transport_solver_numerics__boundary_conditions_1d__ion__state__energy = transport_solver_numerics__boundary_conditions_1d__ion__state__energy()
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d__ion
    label :: String = ""
    z_n :: Float64 = 0.0
    particles :: transport_solver_numerics__boundary_conditions_1d__ion__particles = transport_solver_numerics__boundary_conditions_1d__ion__particles()
    multiple_states_flag :: Int32 = 0
    state :: StructArray{transport_solver_numerics__boundary_conditions_1d__ion__state} = StructArray(transport_solver_numerics__boundary_conditions_1d__ion__state() for k in 1:1)
    energy :: transport_solver_numerics__boundary_conditions_1d__ion__energy = transport_solver_numerics__boundary_conditions_1d__ion__energy()
    a :: Float64 = 0.0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_1d
    ion :: StructArray{transport_solver_numerics__boundary_conditions_1d__ion} = StructArray(transport_solver_numerics__boundary_conditions_1d__ion() for k in 1:1)
    current :: transport_solver_numerics__boundary_conditions_1d__current = transport_solver_numerics__boundary_conditions_1d__current()
    time :: Float64 = 0.0
    momentum_tor :: transport_solver_numerics__boundary_conditions_1d__momentum_tor = transport_solver_numerics__boundary_conditions_1d__momentum_tor()
    electrons :: transport_solver_numerics__boundary_conditions_1d__electrons = transport_solver_numerics__boundary_conditions_1d__electrons()
    energy_ion_total :: transport_solver_numerics__boundary_conditions_1d__energy_ion_total = transport_solver_numerics__boundary_conditions_1d__energy_ion_total()
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__electrons__d_dt
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__electrons
    d2_drho_tor_norm2 :: transport_solver_numerics__derivatives_1d__electrons__d2_drho_tor_norm2 = transport_solver_numerics__derivatives_1d__electrons__d2_drho_tor_norm2()
    d_dt :: transport_solver_numerics__derivatives_1d__electrons__d_dt = transport_solver_numerics__derivatives_1d__electrons__d_dt()
    d_drho_tor_norm :: transport_solver_numerics__derivatives_1d__electrons__d_drho_tor_norm = transport_solver_numerics__derivatives_1d__electrons__d_drho_tor_norm()
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion__d_dt
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    z_n :: Float64 = 0.0
    d2_drho_tor_norm2 :: transport_solver_numerics__derivatives_1d__ion__d2_drho_tor_norm2 = transport_solver_numerics__derivatives_1d__ion__d2_drho_tor_norm2()
    multiple_states_flag :: Int32 = 0
    state :: StructArray{transport_solver_numerics__derivatives_1d__ion__state} = StructArray(transport_solver_numerics__derivatives_1d__ion__state() for k in 1:1)
    a :: Float64 = 0.0
    d_dt :: transport_solver_numerics__derivatives_1d__ion__d_dt = transport_solver_numerics__derivatives_1d__ion__d_dt()
    d_drho_tor_norm :: transport_solver_numerics__derivatives_1d__ion__d_drho_tor_norm = transport_solver_numerics__derivatives_1d__ion__d_drho_tor_norm()
end

Base.@kwdef mutable struct transport_solver_numerics__derivatives_1d
    dpsi_dt_cphi :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    dpsi_drho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    dpsi_dt_crho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    ion :: StructArray{transport_solver_numerics__derivatives_1d__ion} = StructArray(transport_solver_numerics__derivatives_1d__ion() for k in 1:1)
    dpsi_dt :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: transport_solver_numerics__derivatives_1d__electrons = transport_solver_numerics__derivatives_1d__electrons()
    d_drho_tor_norm :: transport_solver_numerics__derivatives_1d__d_drho_tor_norm = transport_solver_numerics__derivatives_1d__d_drho_tor_norm()
    d2_drho_tor_norm2 :: transport_solver_numerics__derivatives_1d__d2_drho_tor_norm2 = transport_solver_numerics__derivatives_1d__d2_drho_tor_norm2()
    drho_tor_dt :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: transport_solver_numerics__derivatives_1d__grid = transport_solver_numerics__derivatives_1d__grid()
    d2psi_drho_tor2 :: Array{Float64, 1} = zeros(Float64,(0))
    d_dvolume_drho_tor_dt :: Array{Float64, 1} = zeros(Float64,(0))
    d_dt :: transport_solver_numerics__derivatives_1d__d_dt = transport_solver_numerics__derivatives_1d__d_dt()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd__current
    grid_index :: Int32 = 0
    values :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_subset_index :: Int32 = 0
    identifier :: transport_solver_numerics__boundary_conditions_ggd__current__identifier = transport_solver_numerics__boundary_conditions_ggd__current__identifier()
end

Base.@kwdef mutable struct transport_solver_numerics__boundary_conditions_ggd
    ion :: StructArray{transport_solver_numerics__boundary_conditions_ggd__ion} = StructArray(transport_solver_numerics__boundary_conditions_ggd__ion() for k in 1:1)
    time :: Float64 = 0.0
    grid :: transport_solver_numerics__boundary_conditions_ggd__grid = transport_solver_numerics__boundary_conditions_ggd__grid()
    current :: StructArray{transport_solver_numerics__boundary_conditions_ggd__current} = StructArray(transport_solver_numerics__boundary_conditions_ggd__current() for k in 1:1)
    electrons :: transport_solver_numerics__boundary_conditions_ggd__electrons = transport_solver_numerics__boundary_conditions_ggd__electrons()
end

Base.@kwdef mutable struct transport_solver_numerics
    time :: Array{Float64, 1} = zeros(Float64,(0))
    restart_files :: StructArray{transport_solver_numerics__restart_files} = StructArray(transport_solver_numerics__restart_files() for k in 1:1)
    time_step_min :: transport_solver_numerics__time_step_min = transport_solver_numerics__time_step_min()
    code :: transport_solver_numerics__code = transport_solver_numerics__code()
    time_step_average :: transport_solver_numerics__time_step_average = transport_solver_numerics__time_step_average()
    solver :: transport_solver_numerics__solver = transport_solver_numerics__solver()
    primary_coordinate :: transport_solver_numerics__primary_coordinate = transport_solver_numerics__primary_coordinate()
    convergence :: transport_solver_numerics__convergence = transport_solver_numerics__convergence()
    boundary_conditions_ggd :: StructArray{transport_solver_numerics__boundary_conditions_ggd} = StructArray(transport_solver_numerics__boundary_conditions_ggd() for k in 1:1)
    ids_properties :: transport_solver_numerics__ids_properties = transport_solver_numerics__ids_properties()
    vacuum_toroidal_field :: transport_solver_numerics__vacuum_toroidal_field = transport_solver_numerics__vacuum_toroidal_field()
    boundary_conditions_1d :: StructArray{transport_solver_numerics__boundary_conditions_1d} = StructArray(transport_solver_numerics__boundary_conditions_1d() for k in 1:1)
    time_step :: transport_solver_numerics__time_step = transport_solver_numerics__time_step()
    derivatives_1d :: StructArray{transport_solver_numerics__derivatives_1d} = StructArray(transport_solver_numerics__derivatives_1d() for k in 1:1)
    solver_1d :: StructArray{transport_solver_numerics__solver_1d} = StructArray(transport_solver_numerics__solver_1d() for k in 1:1)
end

Base.@kwdef mutable struct thomson_scattering__equilibrium_id__data_entry
    pulse_type :: String = ""
    run :: Int32 = 0
    machine :: String = ""
    pulse :: Int32 = 0
    user :: String = ""
end

Base.@kwdef mutable struct thomson_scattering__channel__n_e
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct thomson_scattering__channel__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct thomson_scattering__channel__t_e
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct thomson_scattering__equilibrium_id
    name :: String = ""
    data_entry :: thomson_scattering__equilibrium_id__data_entry = thomson_scattering__equilibrium_id__data_entry()
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct thomson_scattering__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct thomson_scattering__channel__delta_position
    time :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct thomson_scattering__channel__distance_separatrix_midplane
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct thomson_scattering__channel
    name :: String = ""
    distance_separatrix_midplane :: thomson_scattering__channel__distance_separatrix_midplane = thomson_scattering__channel__distance_separatrix_midplane()
    t_e :: thomson_scattering__channel__t_e = thomson_scattering__channel__t_e()
    n_e :: thomson_scattering__channel__n_e = thomson_scattering__channel__n_e()
    position :: thomson_scattering__channel__position = thomson_scattering__channel__position()
    delta_position :: thomson_scattering__channel__delta_position = thomson_scattering__channel__delta_position()
    identifier :: String = ""
end

Base.@kwdef mutable struct thomson_scattering__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct thomson_scattering__ids_properties
    provider :: String = ""
    version_put :: thomson_scattering__ids_properties__version_put = thomson_scattering__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct thomson_scattering__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct thomson_scattering__code
    library :: StructArray{thomson_scattering__code__library} = StructArray(thomson_scattering__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct thomson_scattering
    equilibrium_id :: thomson_scattering__equilibrium_id = thomson_scattering__equilibrium_id()
    latency :: Float64 = 0.0
    ids_properties :: thomson_scattering__ids_properties = thomson_scattering__ids_properties()
    channel :: StructArray{thomson_scattering__channel} = StructArray(thomson_scattering__channel() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    midplane :: thomson_scattering__midplane = thomson_scattering__midplane()
    code :: thomson_scattering__code = thomson_scattering__code()
end

Base.@kwdef mutable struct tf__coil__conductor__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct tf__coil__conductor__cross_section
    delta_z :: Array{Float64, 1} = zeros(Float64,(0))
    delta_r :: Array{Float64, 1} = zeros(Float64,(0))
    delta_phi :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__field_map__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct tf__field_map__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct tf__coil__conductor__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__field_map__a_field_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct tf__coil__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__field_map__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct tf__field_map__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct tf__field_map__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{tf__field_map__grid__space__objects_per_dimension__object__boundary} = StructArray(tf__field_map__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct tf__field_map__grid__space__objects_per_dimension
    object :: StructArray{tf__field_map__grid__space__objects_per_dimension__object} = StructArray(tf__field_map__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct tf__field_map__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct tf__b_field_tor_vacuum_r
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__field_map__b_field_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct tf__field_map__b_field_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct tf__coil__conductor__elements__intermediate_points
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__coil__conductor__elements__start_points
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__field_map__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{tf__field_map__grid__space__objects_per_dimension} = StructArray(tf__field_map__grid__space__objects_per_dimension() for k in 1:1)
    identifier :: tf__field_map__grid__space__identifier = tf__field_map__grid__space__identifier()
    geometry_type :: tf__field_map__grid__space__geometry_type = tf__field_map__grid__space__geometry_type()
end

Base.@kwdef mutable struct tf__field_map__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct tf__field_map__a_field_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct tf__coil__conductor__elements__end_points
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct tf__field_map__a_field_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct tf__coil__conductor__elements__centres
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__ids_properties
    provider :: String = ""
    version_put :: tf__ids_properties__version_put = tf__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct tf__field_map__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct tf__delta_b_field_tor_vacuum_r
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__field_map__b_field_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct tf__coil__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct tf__coil__conductor__elements
    start_points :: tf__coil__conductor__elements__start_points = tf__coil__conductor__elements__start_points()
    names :: Array{String, 1} = String[]
    types :: Array{Int32, 1} = zeros(Int32,(0))
    end_points :: tf__coil__conductor__elements__end_points = tf__coil__conductor__elements__end_points()
    intermediate_points :: tf__coil__conductor__elements__intermediate_points = tf__coil__conductor__elements__intermediate_points()
    centres :: tf__coil__conductor__elements__centres = tf__coil__conductor__elements__centres()
end

Base.@kwdef mutable struct tf__coil__conductor
    voltage :: tf__coil__conductor__voltage = tf__coil__conductor__voltage()
    elements :: tf__coil__conductor__elements = tf__coil__conductor__elements()
    cross_section :: tf__coil__conductor__cross_section = tf__coil__conductor__cross_section()
    resistance :: Float64 = 0.0
    current :: tf__coil__conductor__current = tf__coil__conductor__current()
end

Base.@kwdef mutable struct tf__coil
    voltage :: tf__coil__voltage = tf__coil__voltage()
    turns :: Float64 = 0.0
    conductor :: StructArray{tf__coil__conductor} = StructArray(tf__coil__conductor() for k in 1:1)
    name :: String = ""
    resistance :: Float64 = 0.0
    identifier :: String = ""
    current :: tf__coil__current = tf__coil__current()
end

Base.@kwdef mutable struct tf__code
    library :: StructArray{tf__code__library} = StructArray(tf__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct tf__field_map__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct tf__field_map__grid__grid_subset__element
    object :: StructArray{tf__field_map__grid__grid_subset__element__object} = StructArray(tf__field_map__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct tf__field_map__grid__grid_subset
    base :: StructArray{tf__field_map__grid__grid_subset__base} = StructArray(tf__field_map__grid__grid_subset__base() for k in 1:1)
    metric :: tf__field_map__grid__grid_subset__metric = tf__field_map__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: tf__field_map__grid__grid_subset__identifier = tf__field_map__grid__grid_subset__identifier()
    element :: StructArray{tf__field_map__grid__grid_subset__element} = StructArray(tf__field_map__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct tf__field_map__grid
    grid_subset :: StructArray{tf__field_map__grid__grid_subset} = StructArray(tf__field_map__grid__grid_subset() for k in 1:1)
    space :: StructArray{tf__field_map__grid__space} = StructArray(tf__field_map__grid__space() for k in 1:1)
    identifier :: tf__field_map__grid__identifier = tf__field_map__grid__identifier()
end

Base.@kwdef mutable struct tf__field_map
    b_field_z :: StructArray{tf__field_map__b_field_z} = StructArray(tf__field_map__b_field_z() for k in 1:1)
    time :: Float64 = 0.0
    a_field_z :: StructArray{tf__field_map__a_field_z} = StructArray(tf__field_map__a_field_z() for k in 1:1)
    a_field_tor :: StructArray{tf__field_map__a_field_tor} = StructArray(tf__field_map__a_field_tor() for k in 1:1)
    a_field_r :: StructArray{tf__field_map__a_field_r} = StructArray(tf__field_map__a_field_r() for k in 1:1)
    grid :: tf__field_map__grid = tf__field_map__grid()
    b_field_tor :: StructArray{tf__field_map__b_field_tor} = StructArray(tf__field_map__b_field_tor() for k in 1:1)
    b_field_r :: StructArray{tf__field_map__b_field_r} = StructArray(tf__field_map__b_field_r() for k in 1:1)
end

Base.@kwdef mutable struct tf
    coils_n :: Int32 = 0
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: tf__code = tf__code()
    is_periodic :: Int32 = 0
    latency :: Float64 = 0.0
    delta_b_field_tor_vacuum_r :: tf__delta_b_field_tor_vacuum_r = tf__delta_b_field_tor_vacuum_r()
    field_map :: StructArray{tf__field_map} = StructArray(tf__field_map() for k in 1:1)
    ids_properties :: tf__ids_properties = tf__ids_properties()
    coil :: StructArray{tf__coil} = StructArray(tf__coil() for k in 1:1)
    r0 :: Float64 = 0.0
    b_field_tor_vacuum_r :: tf__b_field_tor_vacuum_r = tf__b_field_tor_vacuum_r()
end

Base.@kwdef mutable struct temporary__dynamic_float4d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
end

Base.@kwdef mutable struct temporary__dynamic_float2d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_string0d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float3d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct temporary__dynamic_integer2d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 2} = zeros(Int32,(0, 0))
end

Base.@kwdef mutable struct temporary__constant_integer1d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float6d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 6} = zeros(Float64,(0,0,0,0,0,0))
end

Base.@kwdef mutable struct temporary__dynamic_integer2d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float5d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_integer1d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_string0d
    value :: String = ""
    identifier :: temporary__constant_string0d__identifier = temporary__constant_string0d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_integer3d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float2d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct temporary__dynamic_float1d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_integer3d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 3} = zeros(Int32,(0, 0, 0))
end

Base.@kwdef mutable struct temporary__constant_integer1d
    value :: Array{Int32, 1} = zeros(Int32,(0))
    identifier :: temporary__constant_integer1d__identifier = temporary__constant_integer1d__identifier()
end

Base.@kwdef mutable struct temporary__constant_integer0d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_integer0d
    value :: Int32 = 0
    identifier :: temporary__constant_integer0d__identifier = temporary__constant_integer0d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float2d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float2d
    value :: Array{Float64, 2} = zeros(Float64,(0,0))
    identifier :: temporary__constant_float2d__identifier = temporary__constant_float2d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float3d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float6d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float6d
    value :: Array{Float64, 6} = zeros(Float64,(0,0,0,0,0,0))
    identifier :: temporary__constant_float6d__identifier = temporary__constant_float6d__identifier()
end

Base.@kwdef mutable struct temporary__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct temporary__dynamic_integer2d
    value :: temporary__dynamic_integer2d__value = temporary__dynamic_integer2d__value()
    identifier :: temporary__dynamic_integer2d__identifier = temporary__dynamic_integer2d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_float6d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float6d
    value :: temporary__dynamic_float6d__value = temporary__dynamic_float6d__value()
    identifier :: temporary__dynamic_float6d__identifier = temporary__dynamic_float6d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_float1d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct temporary__dynamic_float1d
    value :: temporary__dynamic_float1d__value = temporary__dynamic_float1d__value()
    identifier :: temporary__dynamic_float1d__identifier = temporary__dynamic_float1d__identifier()
end

Base.@kwdef mutable struct temporary__constant_integer2d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float3d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float3d
    value :: temporary__dynamic_float3d__value = temporary__dynamic_float3d__value()
    identifier :: temporary__dynamic_float3d__identifier = temporary__dynamic_float3d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float3d
    value :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    identifier :: temporary__constant_float3d__identifier = temporary__constant_float3d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float0d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float0d
    value :: Float64 = 0.0
    identifier :: temporary__constant_float0d__identifier = temporary__constant_float0d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float5d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float5d
    value :: Array{Float64, 5} = zeros(Float64,(0,0,0,0,0))
    identifier :: temporary__constant_float5d__identifier = temporary__constant_float5d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float4d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float4d
    value :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
    identifier :: temporary__constant_float4d__identifier = temporary__constant_float4d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_float4d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_float4d
    value :: temporary__dynamic_float4d__value = temporary__dynamic_float4d__value()
    identifier :: temporary__dynamic_float4d__identifier = temporary__dynamic_float4d__identifier()
end

Base.@kwdef mutable struct temporary__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct temporary__code
    library :: StructArray{temporary__code__library} = StructArray(temporary__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct temporary__constant_string1d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_string1d
    value :: Array{String, 1} = String[]
    identifier :: temporary__constant_string1d__identifier = temporary__constant_string1d__identifier()
end

Base.@kwdef mutable struct temporary__constant_integer2d
    value :: Array{Int32, 2} = zeros(Int32,(0, 0))
    identifier :: temporary__constant_integer2d__identifier = temporary__constant_integer2d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_float5d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 5} = zeros(Float64,(0,0,0,0,0))
end

Base.@kwdef mutable struct temporary__dynamic_float5d
    value :: temporary__dynamic_float5d__value = temporary__dynamic_float5d__value()
    identifier :: temporary__dynamic_float5d__identifier = temporary__dynamic_float5d__identifier()
end

Base.@kwdef mutable struct temporary__constant_float1d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_float1d
    value :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: temporary__constant_float1d__identifier = temporary__constant_float1d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_integer3d
    value :: temporary__dynamic_integer3d__value = temporary__dynamic_integer3d__value()
    identifier :: temporary__dynamic_integer3d__identifier = temporary__dynamic_integer3d__identifier()
end

Base.@kwdef mutable struct temporary__ids_properties
    provider :: String = ""
    version_put :: temporary__ids_properties__version_put = temporary__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct temporary__dynamic_integer1d__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct temporary__dynamic_integer1d
    value :: temporary__dynamic_integer1d__value = temporary__dynamic_integer1d__value()
    identifier :: temporary__dynamic_integer1d__identifier = temporary__dynamic_integer1d__identifier()
end

Base.@kwdef mutable struct temporary__dynamic_float2d
    value :: temporary__dynamic_float2d__value = temporary__dynamic_float2d__value()
    identifier :: temporary__dynamic_float2d__identifier = temporary__dynamic_float2d__identifier()
end

Base.@kwdef mutable struct temporary__constant_integer3d__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct temporary__constant_integer3d
    value :: Array{Int32, 3} = zeros(Int32,(0, 0, 0))
    identifier :: temporary__constant_integer3d__identifier = temporary__constant_integer3d__identifier()
end

Base.@kwdef mutable struct temporary
    dynamic_integer3d :: StructArray{temporary__dynamic_integer3d} = StructArray(temporary__dynamic_integer3d() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    constant_float6d :: StructArray{temporary__constant_float6d} = StructArray(temporary__constant_float6d() for k in 1:1)
    constant_integer1d :: StructArray{temporary__constant_integer1d} = StructArray(temporary__constant_integer1d() for k in 1:1)
    code :: temporary__code = temporary__code()
    dynamic_float2d :: StructArray{temporary__dynamic_float2d} = StructArray(temporary__dynamic_float2d() for k in 1:1)
    dynamic_float6d :: StructArray{temporary__dynamic_float6d} = StructArray(temporary__dynamic_float6d() for k in 1:1)
    dynamic_float4d :: StructArray{temporary__dynamic_float4d} = StructArray(temporary__dynamic_float4d() for k in 1:1)
    dynamic_integer1d :: StructArray{temporary__dynamic_integer1d} = StructArray(temporary__dynamic_integer1d() for k in 1:1)
    dynamic_float1d :: StructArray{temporary__dynamic_float1d} = StructArray(temporary__dynamic_float1d() for k in 1:1)
    constant_float4d :: StructArray{temporary__constant_float4d} = StructArray(temporary__constant_float4d() for k in 1:1)
    constant_float1d :: StructArray{temporary__constant_float1d} = StructArray(temporary__constant_float1d() for k in 1:1)
    constant_string1d :: StructArray{temporary__constant_string1d} = StructArray(temporary__constant_string1d() for k in 1:1)
    dynamic_float3d :: StructArray{temporary__dynamic_float3d} = StructArray(temporary__dynamic_float3d() for k in 1:1)
    constant_integer0d :: StructArray{temporary__constant_integer0d} = StructArray(temporary__constant_integer0d() for k in 1:1)
    constant_string0d :: StructArray{temporary__constant_string0d} = StructArray(temporary__constant_string0d() for k in 1:1)
    constant_float5d :: StructArray{temporary__constant_float5d} = StructArray(temporary__constant_float5d() for k in 1:1)
    ids_properties :: temporary__ids_properties = temporary__ids_properties()
    constant_float3d :: StructArray{temporary__constant_float3d} = StructArray(temporary__constant_float3d() for k in 1:1)
    dynamic_integer2d :: StructArray{temporary__dynamic_integer2d} = StructArray(temporary__dynamic_integer2d() for k in 1:1)
    constant_float0d :: StructArray{temporary__constant_float0d} = StructArray(temporary__constant_float0d() for k in 1:1)
    constant_float2d :: StructArray{temporary__constant_float2d} = StructArray(temporary__constant_float2d() for k in 1:1)
    constant_integer2d :: StructArray{temporary__constant_integer2d} = StructArray(temporary__constant_integer2d() for k in 1:1)
    constant_integer3d :: StructArray{temporary__constant_integer3d} = StructArray(temporary__constant_integer3d() for k in 1:1)
    dynamic_float5d :: StructArray{temporary__dynamic_float5d} = StructArray(temporary__dynamic_float5d() for k in 1:1)
end

Base.@kwdef mutable struct summary__gas_injection_rates__methane_carbon_13
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pedestal_top_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__k_perpendicular
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_thermal_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__dn_e_dt
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__oxygen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_electrons_thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__magnetic_shear_flag
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__global_quantities__v_loop
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__methane_deuterated
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__runaways__current
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__denergy_thermal_dt
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_mhd
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__carbon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__current
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__current
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_ohm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__tau_helium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_line
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__energy
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__energy_fast
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__ethane
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_radiated_inside_lcfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__energy_thermal_pedestal_electron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__krypton
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__offset
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_lh
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__heat_flux_i_decay_length
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__offset
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_height
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__volume_inside_pedestal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__species__label
    source :: String = ""
    value :: String = ""
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__stationary_phase_flag
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__global_quantities__li_mhd
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__configuration
    source :: String = ""
    value :: String = ""
end

Base.@kwdef mutable struct summary__global_quantities__volume
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__fusion_fluence
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__zeff
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__lithium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__geometric_axis_r
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__frequency
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__li
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__beam_beam
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__tau_energy_98
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__deuterium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__position__phi
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__global_quantities__h_98
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__deuterium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__disruption__time
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__current_bootstrap
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_ion_total_thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__frequency
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__wall__evaporation
    source :: String = ""
    value :: String = ""
end

Base.@kwdef mutable struct summary__gas_injection_rates__bottom
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__beam_thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__beam_power_fraction
    source :: String = ""
    value :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__t_i_average_decay_length
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__pedestal_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_additional
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__offset
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__species__a
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__line_average__n_i__deuterium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__power
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__beam_beam
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__n_i_total_decay_length
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__xenon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__strike_point_outer_z
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__helium_4
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__nustar_pedestal_top_electron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__hydrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_norm_mhd
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__disruption__time_radiated_power_max
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__global_quantities__fusion_gain
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__helium_4
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__beryllium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__position__r
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__zeff
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__ip
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pedestal_top_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__power_launched
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_synchrotron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__q_95
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__neon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__gap_limiter_wall
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__pressure_neutral
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__pedestal_height
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__current_ohm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_height
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__midplane
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__nitrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__resistance
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__t_e
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_nbi_co_injected_ratio
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__triangularity_upper
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__volume_inside_pedestal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__carbon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__harmonic
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__n_tor
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__tritium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__hydrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__n_e_decay_length
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__argon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__isotope_fraction_hydrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__disruption__mitigation_valve
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__n_parallel
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__b0
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__nustar_pedestal_top_electron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_ic
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__angle_tor
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__current
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__beam_current_fraction
    source :: String = ""
    value :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter
    alpha_critical :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical()
    t_e_pedestal_top_critical :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical()
    alpha_ratio :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_width
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__xenon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__disruption__vertical_displacement
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__global_quantities__denergy_diamagnetic_dt
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__neon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_width
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__elms__frequency
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__angle_pol
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__isotope_fraction_hydrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__helium_3
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__kicks__occurrence
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__kicks
    occurrence :: summary__kicks__occurrence = summary__kicks__occurrence()
end

Base.@kwdef mutable struct summary__volume_average__meff_hydrogenic
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__alpha_experimental
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__methane
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__dn_e_dt
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_nbi
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_ic
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__xenon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__limiter__material
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct summary__limiter
    material :: summary__limiter__material = summary__limiter__material()
end

Base.@kwdef mutable struct summary__gas_injection_rates__ammonia_deuterated
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__strike_point_outer_r
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__separatrix
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__energy_fast
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__beta_pol_mhd
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__tungsten
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__meff_hydrogenic
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__power_launched
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__magnetic_axis_r
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_ec
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__alpha_electron_pedestal_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_loss
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__tau_energy
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_mhd
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__energy_fast
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__helium_3
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__oxygen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__current
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__krypton
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__d_dpsi_norm_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_ec
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__current
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt
    beam_thermal :: summary__fusion__neutron_fluxes__dt__beam_thermal = summary__fusion__neutron_fluxes__dt__beam_thermal()
    beam_beam :: summary__fusion__neutron_fluxes__dt__beam_beam = summary__fusion__neutron_fluxes__dt__beam_beam()
    total :: summary__fusion__neutron_fluxes__dt__total = summary__fusion__neutron_fluxes__dt__total()
    thermal :: summary__fusion__neutron_fluxes__dt__thermal = summary__fusion__neutron_fluxes__dt__thermal()
end

Base.@kwdef mutable struct summary__boundary__elongation
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_e
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__coulomb_factor_pedestal_top
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_steady
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__h_mode
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__global_quantities__greenwald_fraction
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_nbi
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__helium_3
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__beam_thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_power_total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__carbon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__d_dpsi_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__wall__material
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct summary__wall
    material :: summary__wall__material = summary__wall__material()
    evaporation :: summary__wall__evaporation = summary__wall__evaporation()
end

Base.@kwdef mutable struct summary__global_quantities__current_non_inductive
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__tritium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__energy_thermal_pedestal_ion
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__strike_point_inner_r
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__geometric_axis_z
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__beam_beam
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt
    beam_thermal :: summary__fusion__neutron_fluxes__tt__beam_thermal = summary__fusion__neutron_fluxes__tt__beam_thermal()
    beam_beam :: summary__fusion__neutron_fluxes__tt__beam_beam = summary__fusion__neutron_fluxes__tt__beam_beam()
    total :: summary__fusion__neutron_fluxes__tt__total = summary__fusion__neutron_fluxes__tt__total()
    thermal :: summary__fusion__neutron_fluxes__tt__thermal = summary__fusion__neutron_fluxes__tt__thermal()
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__neon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__ethylene
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__oxygen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__offset
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__polarisation
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__line_average__t_i_average
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__top
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__power
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__power_launched
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__offset
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__t_i_average
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__elms__type
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__elms
    frequency :: summary__elms__frequency = summary__elms__frequency()
    type :: summary__elms__type = summary__elms__type()
end

Base.@kwdef mutable struct summary__global_quantities__tau_resistive
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_height
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__strike_point_configuration
    source :: String = ""
    value :: String = ""
end

Base.@kwdef mutable struct summary__scrape_off_layer__power_radiated
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__iron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__species__z_n
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__line_average__n_i__argon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__phase
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__heat_flux_e_decay_length
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__strike_point_inner_z
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__triangularity_lower
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__pedestal_width
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__magnetic_axis_z
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_width
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__e_field_plus_minus_ratio
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__tangency_radius
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__position__z
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__position
    phi :: summary__heating_current_drive__nbi__position__phi = summary__heating_current_drive__nbi__position__phi()
    r :: summary__heating_current_drive__nbi__position__r = summary__heating_current_drive__nbi__position__r()
    z :: summary__heating_current_drive__nbi__position__z = summary__heating_current_drive__nbi__position__z()
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__pedestal_height
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ec__power_launched
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__argon
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__t_e
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__tag
    name :: String = ""
    comment :: String = ""
end

Base.@kwdef mutable struct summary__runaways__particles
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__runaways
    particles :: summary__runaways__particles = summary__runaways__particles()
    current :: summary__runaways__current = summary__runaways__current()
end

Base.@kwdef mutable struct summary__global_quantities__energy_b_field_pol
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_diamagnetic
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__separatrix
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron
    pedestal_width :: summary__pedestal_fits__linear__pressure_electron__pedestal_width = summary__pedestal_fits__linear__pressure_electron__pedestal_width()
    d_dpsi_norm_max_position :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position = summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position()
    d_dpsi_norm_max :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max = summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max()
    offset :: summary__pedestal_fits__linear__pressure_electron__offset = summary__pedestal_fits__linear__pressure_electron__offset()
    d_dpsi_norm :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm = summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm()
    separatrix :: summary__pedestal_fits__linear__pressure_electron__separatrix = summary__pedestal_fits__linear__pressure_electron__separatrix()
    pedestal_position :: summary__pedestal_fits__linear__pressure_electron__pedestal_position = summary__pedestal_fits__linear__pressure_electron__pedestal_position()
    pedestal_height :: summary__pedestal_fits__linear__pressure_electron__pedestal_height = summary__pedestal_fits__linear__pressure_electron__pedestal_height()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__separatrix
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__d_dpsi_norm_max
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__frequency
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__type
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager
    alpha_critical :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical()
    t_e_pedestal_top_critical :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical()
    alpha_ratio :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability
    alpha_experimental :: summary__pedestal_fits__mtanh__stability__alpha_experimental = summary__pedestal_fits__mtanh__stability__alpha_experimental()
    bootstrap_current_hager :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager()
    bootstrap_current_sauter :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter()
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__pedestal_position
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__nitrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__beam_thermal
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__code
    library :: StructArray{summary__code__library} = StructArray(summary__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_width
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__disruption__time_half_ip
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__volume_average__n_e
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__silane
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__rmps__occurrence
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__rmps
    occurrence :: summary__rmps__occurrence = summary__rmps__occurrence()
end

Base.@kwdef mutable struct summary__global_quantities__power_radiated
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__offset
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_bremsstrahlung
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__helium_4
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__lithium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__propane
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__lh__power
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_fast_perpendicular
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__hydrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__power_radiated_outside_lcfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__current_alignment
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__beta_pol
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary__minor_radius
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__boundary
    strike_point_inner_r :: summary__boundary__strike_point_inner_r = summary__boundary__strike_point_inner_r()
    strike_point_outer_z :: summary__boundary__strike_point_outer_z = summary__boundary__strike_point_outer_z()
    geometric_axis_z :: summary__boundary__geometric_axis_z = summary__boundary__geometric_axis_z()
    magnetic_axis_z :: summary__boundary__magnetic_axis_z = summary__boundary__magnetic_axis_z()
    geometric_axis_r :: summary__boundary__geometric_axis_r = summary__boundary__geometric_axis_r()
    strike_point_configuration :: summary__boundary__strike_point_configuration = summary__boundary__strike_point_configuration()
    triangularity_upper :: summary__boundary__triangularity_upper = summary__boundary__triangularity_upper()
    gap_limiter_wall :: summary__boundary__gap_limiter_wall = summary__boundary__gap_limiter_wall()
    strike_point_inner_z :: summary__boundary__strike_point_inner_z = summary__boundary__strike_point_inner_z()
    triangularity_lower :: summary__boundary__triangularity_lower = summary__boundary__triangularity_lower()
    minor_radius :: summary__boundary__minor_radius = summary__boundary__minor_radius()
    strike_point_outer_r :: summary__boundary__strike_point_outer_r = summary__boundary__strike_point_outer_r()
    elongation :: summary__boundary__elongation = summary__boundary__elongation()
    type :: summary__boundary__type = summary__boundary__type()
    magnetic_axis_r :: summary__boundary__magnetic_axis_r = summary__boundary__magnetic_axis_r()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__beryllium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__disruption
    time :: summary__disruption__time = summary__disruption__time()
    vertical_displacement :: summary__disruption__vertical_displacement = summary__disruption__vertical_displacement()
    mitigation_valve :: summary__disruption__mitigation_valve = summary__disruption__mitigation_valve()
    time_half_ip :: summary__disruption__time_half_ip = summary__disruption__time_half_ip()
    time_radiated_power_max :: summary__disruption__time_radiated_power_max = summary__disruption__time_radiated_power_max()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e
    pedestal_width :: summary__pedestal_fits__mtanh__t_e__pedestal_width = summary__pedestal_fits__mtanh__t_e__pedestal_width()
    d_dpsi_norm_max_position :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position = summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position()
    d_dpsi_norm_max :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max = summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max()
    offset :: summary__pedestal_fits__mtanh__t_e__offset = summary__pedestal_fits__mtanh__t_e__offset()
    d_dpsi_norm :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm = summary__pedestal_fits__mtanh__t_e__d_dpsi_norm()
    pedestal_position :: summary__pedestal_fits__mtanh__t_e__pedestal_position = summary__pedestal_fits__mtanh__t_e__pedestal_position()
    pedestal_height :: summary__pedestal_fits__mtanh__t_e__pedestal_height = summary__pedestal_fits__mtanh__t_e__pedestal_height()
end

Base.@kwdef mutable struct summary__line_average__n_i__beryllium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i_total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__power
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron
    pedestal_width :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_width = summary__pedestal_fits__mtanh__pressure_electron__pedestal_width()
    d_dpsi_norm_max_position :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position = summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position()
    d_dpsi_norm_max :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max = summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max()
    offset :: summary__pedestal_fits__mtanh__pressure_electron__offset = summary__pedestal_fits__mtanh__pressure_electron__offset()
    d_dpsi_norm :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm = summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm()
    separatrix :: summary__pedestal_fits__mtanh__pressure_electron__separatrix = summary__pedestal_fits__mtanh__pressure_electron__separatrix()
    pedestal_position :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_position = summary__pedestal_fits__mtanh__pressure_electron__pedestal_position()
    pedestal_height :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_height = summary__pedestal_fits__mtanh__pressure_electron__pedestal_height()
end

Base.@kwdef mutable struct summary__heating_current_drive__lh
    power_launched :: summary__heating_current_drive__lh__power_launched = summary__heating_current_drive__lh__power_launched()
    energy_fast :: summary__heating_current_drive__lh__energy_fast = summary__heating_current_drive__lh__energy_fast()
    power :: summary__heating_current_drive__lh__power = summary__heating_current_drive__lh__power()
    frequency :: summary__heating_current_drive__lh__frequency = summary__heating_current_drive__lh__frequency()
    position :: summary__heating_current_drive__lh__position = summary__heating_current_drive__lh__position()
    n_parallel :: summary__heating_current_drive__lh__n_parallel = summary__heating_current_drive__lh__n_parallel()
    current :: summary__heating_current_drive__lh__current = summary__heating_current_drive__lh__current()
end

Base.@kwdef mutable struct summary__volume_average__n_i_total
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__species
    label :: summary__heating_current_drive__nbi__species__label = summary__heating_current_drive__nbi__species__label()
    z_n :: summary__heating_current_drive__nbi__species__z_n = summary__heating_current_drive__nbi__species__z_n()
    a :: summary__heating_current_drive__nbi__species__a = summary__heating_current_drive__nbi__species__a()
end

Base.@kwdef mutable struct summary__global_quantities__r0
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__separatrix
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates__impurity_seeding
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__pellets__occurrence
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__pellets
    occurrence :: summary__pellets__occurrence = summary__pellets__occurrence()
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_lh
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__nitrogen
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__iron
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer__t_e_decay_length
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__power
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__scrape_off_layer
    n_e_decay_length :: summary__scrape_off_layer__n_e_decay_length = summary__scrape_off_layer__n_e_decay_length()
    heat_flux_i_decay_length :: summary__scrape_off_layer__heat_flux_i_decay_length = summary__scrape_off_layer__heat_flux_i_decay_length()
    pressure_neutral :: summary__scrape_off_layer__pressure_neutral = summary__scrape_off_layer__pressure_neutral()
    heat_flux_e_decay_length :: summary__scrape_off_layer__heat_flux_e_decay_length = summary__scrape_off_layer__heat_flux_e_decay_length()
    power_radiated :: summary__scrape_off_layer__power_radiated = summary__scrape_off_layer__power_radiated()
    t_e_decay_length :: summary__scrape_off_layer__t_e_decay_length = summary__scrape_off_layer__t_e_decay_length()
    n_i_total_decay_length :: summary__scrape_off_layer__n_i_total_decay_length = summary__scrape_off_layer__n_i_total_decay_length()
    t_i_average_decay_length :: summary__scrape_off_layer__t_i_average_decay_length = summary__scrape_off_layer__t_i_average_decay_length()
end

Base.@kwdef mutable struct summary__line_average__n_i__tungsten
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__direction
    source :: String = ""
    value :: Int32 = 0
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__d_dpsi_norm
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e
    pedestal_width :: summary__pedestal_fits__linear__n_e__pedestal_width = summary__pedestal_fits__linear__n_e__pedestal_width()
    d_dpsi_norm_max :: summary__pedestal_fits__linear__n_e__d_dpsi_norm_max = summary__pedestal_fits__linear__n_e__d_dpsi_norm_max()
    offset :: summary__pedestal_fits__linear__n_e__offset = summary__pedestal_fits__linear__n_e__offset()
    d_dpsi_norm :: summary__pedestal_fits__linear__n_e__d_dpsi_norm = summary__pedestal_fits__linear__n_e__d_dpsi_norm()
    separatrix :: summary__pedestal_fits__linear__n_e__separatrix = summary__pedestal_fits__linear__n_e__separatrix()
    pedestal_position :: summary__pedestal_fits__linear__n_e__pedestal_position = summary__pedestal_fits__linear__n_e__pedestal_position()
    pedestal_height :: summary__pedestal_fits__linear__n_e__pedestal_height = summary__pedestal_fits__linear__n_e__pedestal_height()
end

Base.@kwdef mutable struct summary__gas_injection_rates__ammonia
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__gas_injection_rates
    krypton :: summary__gas_injection_rates__krypton = summary__gas_injection_rates__krypton()
    lithium :: summary__gas_injection_rates__lithium = summary__gas_injection_rates__lithium()
    methane_carbon_13 :: summary__gas_injection_rates__methane_carbon_13 = summary__gas_injection_rates__methane_carbon_13()
    neon :: summary__gas_injection_rates__neon = summary__gas_injection_rates__neon()
    ethylene :: summary__gas_injection_rates__ethylene = summary__gas_injection_rates__ethylene()
    tritium :: summary__gas_injection_rates__tritium = summary__gas_injection_rates__tritium()
    helium_3 :: summary__gas_injection_rates__helium_3 = summary__gas_injection_rates__helium_3()
    deuterium :: summary__gas_injection_rates__deuterium = summary__gas_injection_rates__deuterium()
    bottom :: summary__gas_injection_rates__bottom = summary__gas_injection_rates__bottom()
    helium_4 :: summary__gas_injection_rates__helium_4 = summary__gas_injection_rates__helium_4()
    oxygen :: summary__gas_injection_rates__oxygen = summary__gas_injection_rates__oxygen()
    xenon :: summary__gas_injection_rates__xenon = summary__gas_injection_rates__xenon()
    methane :: summary__gas_injection_rates__methane = summary__gas_injection_rates__methane()
    hydrogen :: summary__gas_injection_rates__hydrogen = summary__gas_injection_rates__hydrogen()
    carbon :: summary__gas_injection_rates__carbon = summary__gas_injection_rates__carbon()
    ammonia_deuterated :: summary__gas_injection_rates__ammonia_deuterated = summary__gas_injection_rates__ammonia_deuterated()
    total :: summary__gas_injection_rates__total = summary__gas_injection_rates__total()
    top :: summary__gas_injection_rates__top = summary__gas_injection_rates__top()
    midplane :: summary__gas_injection_rates__midplane = summary__gas_injection_rates__midplane()
    silane :: summary__gas_injection_rates__silane = summary__gas_injection_rates__silane()
    ethane :: summary__gas_injection_rates__ethane = summary__gas_injection_rates__ethane()
    ammonia :: summary__gas_injection_rates__ammonia = summary__gas_injection_rates__ammonia()
    propane :: summary__gas_injection_rates__propane = summary__gas_injection_rates__propane()
    nitrogen :: summary__gas_injection_rates__nitrogen = summary__gas_injection_rates__nitrogen()
    methane_deuterated :: summary__gas_injection_rates__methane_deuterated = summary__gas_injection_rates__methane_deuterated()
    impurity_seeding :: summary__gas_injection_rates__impurity_seeding = summary__gas_injection_rates__impurity_seeding()
    argon :: summary__gas_injection_rates__argon = summary__gas_injection_rates__argon()
    beryllium :: summary__gas_injection_rates__beryllium = summary__gas_injection_rates__beryllium()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_height
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e
    pedestal_width :: summary__pedestal_fits__mtanh__n_e__pedestal_width = summary__pedestal_fits__mtanh__n_e__pedestal_width()
    d_dpsi_norm_max_position :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position = summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position()
    d_dpsi_norm_max :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max = summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max()
    offset :: summary__pedestal_fits__mtanh__n_e__offset = summary__pedestal_fits__mtanh__n_e__offset()
    d_dpsi_norm :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm = summary__pedestal_fits__mtanh__n_e__d_dpsi_norm()
    separatrix :: summary__pedestal_fits__mtanh__n_e__separatrix = summary__pedestal_fits__mtanh__n_e__separatrix()
    pedestal_position :: summary__pedestal_fits__mtanh__n_e__pedestal_position = summary__pedestal_fits__mtanh__n_e__pedestal_position()
    pedestal_height :: summary__pedestal_fits__mtanh__n_e__pedestal_height = summary__pedestal_fits__mtanh__n_e__pedestal_height()
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh
    pressure_electron :: summary__pedestal_fits__mtanh__pressure_electron = summary__pedestal_fits__mtanh__pressure_electron()
    b_field_pol_pedestal_top_hfs :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs = summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs()
    alpha_electron_pedestal_max_position :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position = summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position()
    rhostar_pedestal_top_electron_hfs :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs = summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs()
    energy_thermal_pedestal_electron :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron = summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron()
    t_e :: summary__pedestal_fits__mtanh__t_e = summary__pedestal_fits__mtanh__t_e()
    beta_pol_pedestal_top_electron_hfs :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs = summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs()
    parameters :: Array{Float64, 1} = zeros(Float64,(0))
    rhostar_pedestal_top_electron_lfs :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs = summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs()
    stability :: summary__pedestal_fits__mtanh__stability = summary__pedestal_fits__mtanh__stability()
    beta_pol_pedestal_top_electron_lfs :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs = summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs()
    b_field_tor_pedestal_top_lfs :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs = summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs()
    volume_inside_pedestal :: summary__pedestal_fits__mtanh__volume_inside_pedestal = summary__pedestal_fits__mtanh__volume_inside_pedestal()
    b_field_pol_pedestal_top_average :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average = summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average()
    coulomb_factor_pedestal_top :: summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top = summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top()
    beta_pol_pedestal_top_electron_average :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average = summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average()
    b_field_pedestal_top_lfs :: summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs = summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs()
    b_field_pol_pedestal_top_lfs :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs = summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs()
    nustar_pedestal_top_electron :: summary__pedestal_fits__mtanh__nustar_pedestal_top_electron = summary__pedestal_fits__mtanh__nustar_pedestal_top_electron()
    alpha_electron_pedestal_max :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max = summary__pedestal_fits__mtanh__alpha_electron_pedestal_max()
    b_field_pedestal_top_hfs :: summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs = summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs()
    n_e :: summary__pedestal_fits__mtanh__n_e = summary__pedestal_fits__mtanh__n_e()
    energy_thermal_pedestal_ion :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion = summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion()
    rhostar_pedestal_top_electron_magnetic_axis :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis = summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis()
    b_field_tor_pedestal_top_hfs :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs = summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs()
end

Base.@kwdef mutable struct summary__heating_current_drive__ic__harmonic
    source :: String = ""
    value :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__ic
    harmonic :: summary__heating_current_drive__ic__harmonic = summary__heating_current_drive__ic__harmonic()
    e_field_plus_minus_ratio :: summary__heating_current_drive__ic__e_field_plus_minus_ratio = summary__heating_current_drive__ic__e_field_plus_minus_ratio()
    current :: summary__heating_current_drive__ic__current = summary__heating_current_drive__ic__current()
    n_tor :: summary__heating_current_drive__ic__n_tor = summary__heating_current_drive__ic__n_tor()
    power_launched :: summary__heating_current_drive__ic__power_launched = summary__heating_current_drive__ic__power_launched()
    energy_fast :: summary__heating_current_drive__ic__energy_fast = summary__heating_current_drive__ic__energy_fast()
    position :: summary__heating_current_drive__ic__position = summary__heating_current_drive__ic__position()
    power :: summary__heating_current_drive__ic__power = summary__heating_current_drive__ic__power()
    phase :: summary__heating_current_drive__ic__phase = summary__heating_current_drive__ic__phase()
    frequency :: summary__heating_current_drive__ic__frequency = summary__heating_current_drive__ic__frequency()
    k_perpendicular :: summary__heating_current_drive__ic__k_perpendicular = summary__heating_current_drive__ic__k_perpendicular()
end

Base.@kwdef mutable struct summary__global_quantities__ratio_tau_helium_fuel
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd
    beam_thermal :: summary__fusion__neutron_fluxes__dd__beam_thermal = summary__fusion__neutron_fluxes__dd__beam_thermal()
    beam_beam :: summary__fusion__neutron_fluxes__dd__beam_beam = summary__fusion__neutron_fluxes__dd__beam_beam()
    total :: summary__fusion__neutron_fluxes__dd__total = summary__fusion__neutron_fluxes__dd__total()
    thermal :: summary__fusion__neutron_fluxes__dd__thermal = summary__fusion__neutron_fluxes__dd__thermal()
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes
    tt :: summary__fusion__neutron_fluxes__tt = summary__fusion__neutron_fluxes__tt()
    total :: summary__fusion__neutron_fluxes__total = summary__fusion__neutron_fluxes__total()
    thermal :: summary__fusion__neutron_fluxes__thermal = summary__fusion__neutron_fluxes__thermal()
    dt :: summary__fusion__neutron_fluxes__dt = summary__fusion__neutron_fluxes__dt()
    dd :: summary__fusion__neutron_fluxes__dd = summary__fusion__neutron_fluxes__dd()
end

Base.@kwdef mutable struct summary__fusion
    neutron_power_total :: summary__fusion__neutron_power_total = summary__fusion__neutron_power_total()
    neutron_fluxes :: summary__fusion__neutron_fluxes = summary__fusion__neutron_fluxes()
    power :: summary__fusion__power = summary__fusion__power()
    current :: summary__fusion__current = summary__fusion__current()
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_average
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__volume_average__n_i__lithium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi__angle
    source :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi
    energy :: summary__heating_current_drive__nbi__energy = summary__heating_current_drive__nbi__energy()
    angle :: summary__heating_current_drive__nbi__angle = summary__heating_current_drive__nbi__angle()
    current :: summary__heating_current_drive__nbi__current = summary__heating_current_drive__nbi__current()
    beam_current_fraction :: summary__heating_current_drive__nbi__beam_current_fraction = summary__heating_current_drive__nbi__beam_current_fraction()
    power_launched :: summary__heating_current_drive__nbi__power_launched = summary__heating_current_drive__nbi__power_launched()
    position :: summary__heating_current_drive__nbi__position = summary__heating_current_drive__nbi__position()
    tangency_radius :: summary__heating_current_drive__nbi__tangency_radius = summary__heating_current_drive__nbi__tangency_radius()
    direction :: summary__heating_current_drive__nbi__direction = summary__heating_current_drive__nbi__direction()
    power :: summary__heating_current_drive__nbi__power = summary__heating_current_drive__nbi__power()
    beam_power_fraction :: summary__heating_current_drive__nbi__beam_power_fraction = summary__heating_current_drive__nbi__beam_power_fraction()
    species :: summary__heating_current_drive__nbi__species = summary__heating_current_drive__nbi__species()
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i__tritium
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__line_average__n_i
    lithium :: summary__line_average__n_i__lithium = summary__line_average__n_i__lithium()
    krypton :: summary__line_average__n_i__krypton = summary__line_average__n_i__krypton()
    neon :: summary__line_average__n_i__neon = summary__line_average__n_i__neon()
    tritium :: summary__line_average__n_i__tritium = summary__line_average__n_i__tritium()
    helium_3 :: summary__line_average__n_i__helium_3 = summary__line_average__n_i__helium_3()
    iron :: summary__line_average__n_i__iron = summary__line_average__n_i__iron()
    deuterium :: summary__line_average__n_i__deuterium = summary__line_average__n_i__deuterium()
    helium_4 :: summary__line_average__n_i__helium_4 = summary__line_average__n_i__helium_4()
    oxygen :: summary__line_average__n_i__oxygen = summary__line_average__n_i__oxygen()
    tungsten :: summary__line_average__n_i__tungsten = summary__line_average__n_i__tungsten()
    xenon :: summary__line_average__n_i__xenon = summary__line_average__n_i__xenon()
    hydrogen :: summary__line_average__n_i__hydrogen = summary__line_average__n_i__hydrogen()
    carbon :: summary__line_average__n_i__carbon = summary__line_average__n_i__carbon()
    nitrogen :: summary__line_average__n_i__nitrogen = summary__line_average__n_i__nitrogen()
    beryllium :: summary__line_average__n_i__beryllium = summary__line_average__n_i__beryllium()
    argon :: summary__line_average__n_i__argon = summary__line_average__n_i__argon()
end

Base.@kwdef mutable struct summary__line_average
    n_i :: summary__line_average__n_i = summary__line_average__n_i()
    meff_hydrogenic :: summary__line_average__meff_hydrogenic = summary__line_average__meff_hydrogenic()
    t_e :: summary__line_average__t_e = summary__line_average__t_e()
    n_e :: summary__line_average__n_e = summary__line_average__n_e()
    isotope_fraction_hydrogen :: summary__line_average__isotope_fraction_hydrogen = summary__line_average__isotope_fraction_hydrogen()
    t_i_average :: summary__line_average__t_i_average = summary__line_average__t_i_average()
    n_i_total :: summary__line_average__n_i_total = summary__line_average__n_i_total()
    zeff :: summary__line_average__zeff = summary__line_average__zeff()
    dn_e_dt :: summary__line_average__dn_e_dt = summary__line_average__dn_e_dt()
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__pedestal_width
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e
    pedestal_width :: summary__pedestal_fits__linear__t_e__pedestal_width = summary__pedestal_fits__linear__t_e__pedestal_width()
    d_dpsi_norm_max :: summary__pedestal_fits__linear__t_e__d_dpsi_norm_max = summary__pedestal_fits__linear__t_e__d_dpsi_norm_max()
    offset :: summary__pedestal_fits__linear__t_e__offset = summary__pedestal_fits__linear__t_e__offset()
    d_dpsi_norm :: summary__pedestal_fits__linear__t_e__d_dpsi_norm = summary__pedestal_fits__linear__t_e__d_dpsi_norm()
    pedestal_position :: summary__pedestal_fits__linear__t_e__pedestal_position = summary__pedestal_fits__linear__t_e__pedestal_position()
    pedestal_height :: summary__pedestal_fits__linear__t_e__pedestal_height = summary__pedestal_fits__linear__t_e__pedestal_height()
end

Base.@kwdef mutable struct summary__volume_average__n_i__krypton
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities__energy_fast_parallel
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__global_quantities
    power_radiated_outside_lcfs :: summary__global_quantities__power_radiated_outside_lcfs = summary__global_quantities__power_radiated_outside_lcfs()
    beta_tor :: summary__global_quantities__beta_tor = summary__global_quantities__beta_tor()
    power_ohm :: summary__global_quantities__power_ohm = summary__global_quantities__power_ohm()
    energy_b_field_pol :: summary__global_quantities__energy_b_field_pol = summary__global_quantities__energy_b_field_pol()
    tau_energy :: summary__global_quantities__tau_energy = summary__global_quantities__tau_energy()
    energy_electrons_thermal :: summary__global_quantities__energy_electrons_thermal = summary__global_quantities__energy_electrons_thermal()
    greenwald_fraction :: summary__global_quantities__greenwald_fraction = summary__global_quantities__greenwald_fraction()
    power_steady :: summary__global_quantities__power_steady = summary__global_quantities__power_steady()
    resistance :: summary__global_quantities__resistance = summary__global_quantities__resistance()
    beta_pol :: summary__global_quantities__beta_pol = summary__global_quantities__beta_pol()
    energy_thermal :: summary__global_quantities__energy_thermal = summary__global_quantities__energy_thermal()
    fusion_gain :: summary__global_quantities__fusion_gain = summary__global_quantities__fusion_gain()
    power_radiated_inside_lcfs :: summary__global_quantities__power_radiated_inside_lcfs = summary__global_quantities__power_radiated_inside_lcfs()
    beta_tor_norm :: summary__global_quantities__beta_tor_norm = summary__global_quantities__beta_tor_norm()
    beta_tor_thermal_norm :: summary__global_quantities__beta_tor_thermal_norm = summary__global_quantities__beta_tor_thermal_norm()
    tau_helium :: summary__global_quantities__tau_helium = summary__global_quantities__tau_helium()
    power_loss :: summary__global_quantities__power_loss = summary__global_quantities__power_loss()
    current_alignment :: summary__global_quantities__current_alignment = summary__global_quantities__current_alignment()
    energy_diamagnetic :: summary__global_quantities__energy_diamagnetic = summary__global_quantities__energy_diamagnetic()
    tau_energy_98 :: summary__global_quantities__tau_energy_98 = summary__global_quantities__tau_energy_98()
    li_mhd :: summary__global_quantities__li_mhd = summary__global_quantities__li_mhd()
    energy_fast_perpendicular :: summary__global_quantities__energy_fast_perpendicular = summary__global_quantities__energy_fast_perpendicular()
    q_95 :: summary__global_quantities__q_95 = summary__global_quantities__q_95()
    beta_tor_mhd :: summary__global_quantities__beta_tor_mhd = summary__global_quantities__beta_tor_mhd()
    volume :: summary__global_quantities__volume = summary__global_quantities__volume()
    b0 :: summary__global_quantities__b0 = summary__global_quantities__b0()
    ip :: summary__global_quantities__ip = summary__global_quantities__ip()
    beta_pol_mhd :: summary__global_quantities__beta_pol_mhd = summary__global_quantities__beta_pol_mhd()
    denergy_diamagnetic_dt :: summary__global_quantities__denergy_diamagnetic_dt = summary__global_quantities__denergy_diamagnetic_dt()
    h_mode :: summary__global_quantities__h_mode = summary__global_quantities__h_mode()
    power_bremsstrahlung :: summary__global_quantities__power_bremsstrahlung = summary__global_quantities__power_bremsstrahlung()
    h_98 :: summary__global_quantities__h_98 = summary__global_quantities__h_98()
    tau_resistive :: summary__global_quantities__tau_resistive = summary__global_quantities__tau_resistive()
    energy_total :: summary__global_quantities__energy_total = summary__global_quantities__energy_total()
    current_non_inductive :: summary__global_quantities__current_non_inductive = summary__global_quantities__current_non_inductive()
    current_ohm :: summary__global_quantities__current_ohm = summary__global_quantities__current_ohm()
    power_synchrotron :: summary__global_quantities__power_synchrotron = summary__global_quantities__power_synchrotron()
    v_loop :: summary__global_quantities__v_loop = summary__global_quantities__v_loop()
    energy_ion_total_thermal :: summary__global_quantities__energy_ion_total_thermal = summary__global_quantities__energy_ion_total_thermal()
    power_line :: summary__global_quantities__power_line = summary__global_quantities__power_line()
    current_bootstrap :: summary__global_quantities__current_bootstrap = summary__global_quantities__current_bootstrap()
    fusion_fluence :: summary__global_quantities__fusion_fluence = summary__global_quantities__fusion_fluence()
    energy_fast_parallel :: summary__global_quantities__energy_fast_parallel = summary__global_quantities__energy_fast_parallel()
    power_radiated :: summary__global_quantities__power_radiated = summary__global_quantities__power_radiated()
    denergy_thermal_dt :: summary__global_quantities__denergy_thermal_dt = summary__global_quantities__denergy_thermal_dt()
    energy_mhd :: summary__global_quantities__energy_mhd = summary__global_quantities__energy_mhd()
    li :: summary__global_quantities__li = summary__global_quantities__li()
    ratio_tau_helium_fuel :: summary__global_quantities__ratio_tau_helium_fuel = summary__global_quantities__ratio_tau_helium_fuel()
    r0 :: summary__global_quantities__r0 = summary__global_quantities__r0()
    beta_tor_norm_mhd :: summary__global_quantities__beta_tor_norm_mhd = summary__global_quantities__beta_tor_norm_mhd()
end

Base.@kwdef mutable struct summary__ids_properties
    provider :: String = ""
    version_put :: summary__ids_properties__version_put = summary__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis
    source :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct summary__pedestal_fits__linear
    pressure_electron :: summary__pedestal_fits__linear__pressure_electron = summary__pedestal_fits__linear__pressure_electron()
    b_field_pol_pedestal_top_hfs :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs = summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs()
    rhostar_pedestal_top_electron_hfs :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs = summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs()
    beta_pol_pedestal_top_electron_hfs :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs = summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs()
    t_e :: summary__pedestal_fits__linear__t_e = summary__pedestal_fits__linear__t_e()
    energy_thermal_pedestal_electron :: summary__pedestal_fits__linear__energy_thermal_pedestal_electron = summary__pedestal_fits__linear__energy_thermal_pedestal_electron()
    parameters :: Array{Float64, 1} = zeros(Float64,(0))
    rhostar_pedestal_top_electron_lfs :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs = summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs()
    beta_pol_pedestal_top_electron_lfs :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs = summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs()
    b_field_tor_pedestal_top_lfs :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs = summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs()
    volume_inside_pedestal :: summary__pedestal_fits__linear__volume_inside_pedestal = summary__pedestal_fits__linear__volume_inside_pedestal()
    b_field_pol_pedestal_top_average :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_average = summary__pedestal_fits__linear__b_field_pol_pedestal_top_average()
    coulomb_factor_pedestal_top :: summary__pedestal_fits__linear__coulomb_factor_pedestal_top = summary__pedestal_fits__linear__coulomb_factor_pedestal_top()
    beta_pol_pedestal_top_electron_average :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average = summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average()
    b_field_pedestal_top_lfs :: summary__pedestal_fits__linear__b_field_pedestal_top_lfs = summary__pedestal_fits__linear__b_field_pedestal_top_lfs()
    b_field_pol_pedestal_top_lfs :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs = summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs()
    nustar_pedestal_top_electron :: summary__pedestal_fits__linear__nustar_pedestal_top_electron = summary__pedestal_fits__linear__nustar_pedestal_top_electron()
    b_field_pedestal_top_hfs :: summary__pedestal_fits__linear__b_field_pedestal_top_hfs = summary__pedestal_fits__linear__b_field_pedestal_top_hfs()
    n_e :: summary__pedestal_fits__linear__n_e = summary__pedestal_fits__linear__n_e()
    energy_thermal_pedestal_ion :: summary__pedestal_fits__linear__energy_thermal_pedestal_ion = summary__pedestal_fits__linear__energy_thermal_pedestal_ion()
    rhostar_pedestal_top_electron_magnetic_axis :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis = summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis()
    b_field_tor_pedestal_top_hfs :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs = summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs()
end

Base.@kwdef mutable struct summary__pedestal_fits
    linear :: summary__pedestal_fits__linear = summary__pedestal_fits__linear()
    mtanh :: summary__pedestal_fits__mtanh = summary__pedestal_fits__mtanh()
end

Base.@kwdef mutable struct summary__volume_average__n_i
    krypton :: summary__volume_average__n_i__krypton = summary__volume_average__n_i__krypton()
    lithium :: summary__volume_average__n_i__lithium = summary__volume_average__n_i__lithium()
    neon :: summary__volume_average__n_i__neon = summary__volume_average__n_i__neon()
    tritium :: summary__volume_average__n_i__tritium = summary__volume_average__n_i__tritium()
    helium_3 :: summary__volume_average__n_i__helium_3 = summary__volume_average__n_i__helium_3()
    iron :: summary__volume_average__n_i__iron = summary__volume_average__n_i__iron()
    deuterium :: summary__volume_average__n_i__deuterium = summary__volume_average__n_i__deuterium()
    helium_4 :: summary__volume_average__n_i__helium_4 = summary__volume_average__n_i__helium_4()
    oxygen :: summary__volume_average__n_i__oxygen = summary__volume_average__n_i__oxygen()
    tungsten :: summary__volume_average__n_i__tungsten = summary__volume_average__n_i__tungsten()
    xenon :: summary__volume_average__n_i__xenon = summary__volume_average__n_i__xenon()
    hydrogen :: summary__volume_average__n_i__hydrogen = summary__volume_average__n_i__hydrogen()
    carbon :: summary__volume_average__n_i__carbon = summary__volume_average__n_i__carbon()
    nitrogen :: summary__volume_average__n_i__nitrogen = summary__volume_average__n_i__nitrogen()
    beryllium :: summary__volume_average__n_i__beryllium = summary__volume_average__n_i__beryllium()
    argon :: summary__volume_average__n_i__argon = summary__volume_average__n_i__argon()
end

Base.@kwdef mutable struct summary__volume_average
    n_i :: summary__volume_average__n_i = summary__volume_average__n_i()
    meff_hydrogenic :: summary__volume_average__meff_hydrogenic = summary__volume_average__meff_hydrogenic()
    t_e :: summary__volume_average__t_e = summary__volume_average__t_e()
    n_e :: summary__volume_average__n_e = summary__volume_average__n_e()
    isotope_fraction_hydrogen :: summary__volume_average__isotope_fraction_hydrogen = summary__volume_average__isotope_fraction_hydrogen()
    t_i_average :: summary__volume_average__t_i_average = summary__volume_average__t_i_average()
    n_i_total :: summary__volume_average__n_i_total = summary__volume_average__n_i_total()
    dn_e_dt :: summary__volume_average__dn_e_dt = summary__volume_average__dn_e_dt()
    zeff :: summary__volume_average__zeff = summary__volume_average__zeff()
end

Base.@kwdef mutable struct summary__heating_current_drive__ec
    polarisation :: summary__heating_current_drive__ec__polarisation = summary__heating_current_drive__ec__polarisation()
    angle_tor :: summary__heating_current_drive__ec__angle_tor = summary__heating_current_drive__ec__angle_tor()
    power_launched :: summary__heating_current_drive__ec__power_launched = summary__heating_current_drive__ec__power_launched()
    harmonic :: summary__heating_current_drive__ec__harmonic = summary__heating_current_drive__ec__harmonic()
    energy_fast :: summary__heating_current_drive__ec__energy_fast = summary__heating_current_drive__ec__energy_fast()
    power :: summary__heating_current_drive__ec__power = summary__heating_current_drive__ec__power()
    position :: summary__heating_current_drive__ec__position = summary__heating_current_drive__ec__position()
    frequency :: summary__heating_current_drive__ec__frequency = summary__heating_current_drive__ec__frequency()
    angle_pol :: summary__heating_current_drive__ec__angle_pol = summary__heating_current_drive__ec__angle_pol()
    current :: summary__heating_current_drive__ec__current = summary__heating_current_drive__ec__current()
end

Base.@kwdef mutable struct summary__heating_current_drive
    power_launched_ec :: summary__heating_current_drive__power_launched_ec = summary__heating_current_drive__power_launched_ec()
    ec :: StructArray{summary__heating_current_drive__ec} = StructArray(summary__heating_current_drive__ec() for k in 1:1)
    power_nbi :: summary__heating_current_drive__power_nbi = summary__heating_current_drive__power_nbi()
    power_ec :: summary__heating_current_drive__power_ec = summary__heating_current_drive__power_ec()
    power_lh :: summary__heating_current_drive__power_lh = summary__heating_current_drive__power_lh()
    ic :: StructArray{summary__heating_current_drive__ic} = StructArray(summary__heating_current_drive__ic() for k in 1:1)
    lh :: StructArray{summary__heating_current_drive__lh} = StructArray(summary__heating_current_drive__lh() for k in 1:1)
    power_ic :: summary__heating_current_drive__power_ic = summary__heating_current_drive__power_ic()
    power_launched_lh :: summary__heating_current_drive__power_launched_lh = summary__heating_current_drive__power_launched_lh()
    nbi :: StructArray{summary__heating_current_drive__nbi} = StructArray(summary__heating_current_drive__nbi() for k in 1:1)
    power_launched_nbi_co_injected_ratio :: summary__heating_current_drive__power_launched_nbi_co_injected_ratio = summary__heating_current_drive__power_launched_nbi_co_injected_ratio()
    power_additional :: summary__heating_current_drive__power_additional = summary__heating_current_drive__power_additional()
    power_launched_ic :: summary__heating_current_drive__power_launched_ic = summary__heating_current_drive__power_launched_ic()
    power_launched_nbi :: summary__heating_current_drive__power_launched_nbi = summary__heating_current_drive__power_launched_nbi()
end

Base.@kwdef mutable struct summary
    magnetic_shear_flag :: summary__magnetic_shear_flag = summary__magnetic_shear_flag()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: summary__code = summary__code()
    stationary_phase_flag :: summary__stationary_phase_flag = summary__stationary_phase_flag()
    gas_injection_rates :: summary__gas_injection_rates = summary__gas_injection_rates()
    boundary :: summary__boundary = summary__boundary()
    configuration :: summary__configuration = summary__configuration()
    pedestal_fits :: summary__pedestal_fits = summary__pedestal_fits()
    heating_current_drive :: summary__heating_current_drive = summary__heating_current_drive()
    global_quantities :: summary__global_quantities = summary__global_quantities()
    disruption :: summary__disruption = summary__disruption()
    fusion :: summary__fusion = summary__fusion()
    rmps :: summary__rmps = summary__rmps()
    ids_properties :: summary__ids_properties = summary__ids_properties()
    midplane :: summary__midplane = summary__midplane()
    pellets :: summary__pellets = summary__pellets()
    tag :: summary__tag = summary__tag()
    elms :: summary__elms = summary__elms()
    kicks :: summary__kicks = summary__kicks()
    line_average :: summary__line_average = summary__line_average()
    scrape_off_layer :: summary__scrape_off_layer = summary__scrape_off_layer()
    volume_average :: summary__volume_average = summary__volume_average()
    runaways :: summary__runaways = summary__runaways()
    wall :: summary__wall = summary__wall()
    limiter :: summary__limiter = summary__limiter()
    time_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__filter_window
    material :: String = ""
    thickness :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__t_e_proxy
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__crystal__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__t_i_proxy
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__camera__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__crystal__summit
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__crystal__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__code
    library :: StructArray{spectrometer_x_ray_crystal__code__library} = StructArray(spectrometer_x_ray_crystal__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__velocity_tor_proxy
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__ids_properties
    provider :: String = ""
    version_put :: spectrometer_x_ray_crystal__ids_properties__version_put = spectrometer_x_ray_crystal__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__frame
    photon_count :: Array{Float64, 2} = zeros(Float64,(0,0))
    time :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__camera__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__crystal
    angular_span_horizontal :: Float64 = 0.0
    height :: Float64 = 0.0
    material :: String = ""
    angular_span_vertical :: Float64 = 0.0
    summit :: spectrometer_x_ray_crystal__crystal__summit = spectrometer_x_ray_crystal__crystal__summit()
    angle_bragg :: Float64 = 0.0
    curvature_vertical :: Float64 = 0.0
    curvature_horizontal :: Float64 = 0.0
    width :: Float64 = 0.0
    x3_unit_vector :: spectrometer_x_ray_crystal__crystal__x3_unit_vector = spectrometer_x_ray_crystal__crystal__x3_unit_vector()
    geometry_type :: spectrometer_x_ray_crystal__crystal__geometry_type = spectrometer_x_ray_crystal__crystal__geometry_type()
    wavelength_bragg :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal__camera
    centre :: spectrometer_x_ray_crystal__camera__centre = spectrometer_x_ray_crystal__camera__centre()
    pixel_n :: Array{Int32, 1} = zeros(Int32,(0))
    x3_unit_vector :: spectrometer_x_ray_crystal__camera__x3_unit_vector = spectrometer_x_ray_crystal__camera__x3_unit_vector()
    pixel_dimensions :: Array{Float64, 1} = zeros(Float64,(0))
    camera_dimensions :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_x_ray_crystal
    time :: Array{Float64, 1} = zeros(Float64,(0))
    wavelength_frames :: Array{Float64, 2} = zeros(Float64,(0,0))
    code :: spectrometer_x_ray_crystal__code = spectrometer_x_ray_crystal__code()
    frame :: StructArray{spectrometer_x_ray_crystal__frame} = StructArray(spectrometer_x_ray_crystal__frame() for k in 1:1)
    energy_bound_lower :: Array{Float64, 2} = zeros(Float64,(0,0))
    t_e_proxy :: spectrometer_x_ray_crystal__t_e_proxy = spectrometer_x_ray_crystal__t_e_proxy()
    latency :: Float64 = 0.0
    integration_time :: Float64 = 0.0
    filter_window :: spectrometer_x_ray_crystal__filter_window = spectrometer_x_ray_crystal__filter_window()
    camera :: spectrometer_x_ray_crystal__camera = spectrometer_x_ray_crystal__camera()
    crystal :: spectrometer_x_ray_crystal__crystal = spectrometer_x_ray_crystal__crystal()
    ids_properties :: spectrometer_x_ray_crystal__ids_properties = spectrometer_x_ray_crystal__ids_properties()
    velocity_tor_proxy :: spectrometer_x_ray_crystal__velocity_tor_proxy = spectrometer_x_ray_crystal__velocity_tor_proxy()
    t_i_proxy :: spectrometer_x_ray_crystal__t_i_proxy = spectrometer_x_ray_crystal__t_i_proxy()
    energy_bound_upper :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_frames :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__line_radiances_adjusted
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__validity_timed
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__output_voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__line_intensities
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__calibrated_line_integrals
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__active_spatial_resolution__width
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer__processed_line__radiance
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector_image__circular
    ellipticity :: Float64 = 0.0
    radius :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__light_collection_efficiencies__positions
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__light_collection_efficiencies
    values :: Array{Float64, 1} = zeros(Float64,(0))
    positions :: spectrometer_visible__channel__light_collection_efficiencies__positions = spectrometer_visible__channel__light_collection_efficiencies__positions()
end

Base.@kwdef mutable struct spectrometer_visible__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__line_power_radiances
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer__intensity_spectrum
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__etendue_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer__processed_line__intensity
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__fibre_image__circular
    ellipticity :: Float64 = 0.0
    radius :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer__processed_line
    intensity :: spectrometer_visible__channel__grating_spectrometer__processed_line__intensity = spectrometer_visible__channel__grating_spectrometer__processed_line__intensity()
    label :: String = ""
    radiance :: spectrometer_visible__channel__grating_spectrometer__processed_line__radiance = spectrometer_visible__channel__grating_spectrometer__processed_line__radiance()
    wavelength_central :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__isotope_ratios__isotope__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct spectrometer_visible__code
    library :: StructArray{spectrometer_visible__code__library} = StructArray(spectrometer_visible__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__photoelectric_voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarizer__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarizer__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer__radiance_spectral
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__line_radiances
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarizer__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer__photon_count
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_visible__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct spectrometer_visible__channel__isotope_ratios__method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__aperture
    x1_width :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: spectrometer_visible__channel__aperture__centre = spectrometer_visible__channel__aperture__centre()
    outline :: spectrometer_visible__channel__aperture__outline = spectrometer_visible__channel__aperture__outline()
    radius :: Float64 = 0.0
    x3_unit_vector :: spectrometer_visible__channel__aperture__x3_unit_vector = spectrometer_visible__channel__aperture__x3_unit_vector()
    x2_unit_vector :: spectrometer_visible__channel__aperture__x2_unit_vector = spectrometer_visible__channel__aperture__x2_unit_vector()
    surface :: Float64 = 0.0
    geometry_type :: Int32 = 0
    x1_unit_vector :: spectrometer_visible__channel__aperture__x1_unit_vector = spectrometer_visible__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct spectrometer_visible__channel__isotope_ratios__isotope
    label :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    cold_neutrals_fraction :: Array{Float64, 1} = zeros(Float64,(0))
    hot_neutrals_fraction :: Array{Float64, 1} = zeros(Float64,(0))
    density_ratio :: Array{Float64, 1} = zeros(Float64,(0))
    hot_neutrals_temperature :: Array{Float64, 1} = zeros(Float64,(0))
    cold_neutrals_temperature :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{spectrometer_visible__channel__isotope_ratios__isotope__element} = StructArray(spectrometer_visible__channel__isotope_ratios__isotope__element() for k in 1:1)
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector_image__outline
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer__wavelength_calibration
    offset :: Float64 = 0.0
    gain :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__active_spatial_resolution__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__active_spatial_resolution
    time :: Float64 = 0.0
    centre :: spectrometer_visible__channel__active_spatial_resolution__centre = spectrometer_visible__channel__active_spatial_resolution__centre()
    width :: spectrometer_visible__channel__active_spatial_resolution__width = spectrometer_visible__channel__active_spatial_resolution__width()
end

Base.@kwdef mutable struct spectrometer_visible__channel__fibre_image__outline
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__fibre_image
    circular :: spectrometer_visible__channel__fibre_image__circular = spectrometer_visible__channel__fibre_image__circular()
    outline :: spectrometer_visible__channel__fibre_image__outline = spectrometer_visible__channel__fibre_image__outline()
    geometry_type :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarizer__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__filter_spectrometer
    raw_lines :: Array{Float64, 1} = zeros(Float64,(0))
    calibrated_lines :: Array{Float64, 1} = zeros(Float64,(0))
    line_power_radiances :: spectrometer_visible__channel__filter_spectrometer__line_power_radiances = spectrometer_visible__channel__filter_spectrometer__line_power_radiances()
    calibrated_line_integrals :: spectrometer_visible__channel__filter_spectrometer__calibrated_line_integrals = spectrometer_visible__channel__filter_spectrometer__calibrated_line_integrals()
    line_radiances_adjusted :: spectrometer_visible__channel__filter_spectrometer__line_radiances_adjusted = spectrometer_visible__channel__filter_spectrometer__line_radiances_adjusted()
    radiance_calibration_date :: String = ""
    line_intensities :: spectrometer_visible__channel__filter_spectrometer__line_intensities = spectrometer_visible__channel__filter_spectrometer__line_intensities()
    processed_lines :: Array{Float64, 1} = zeros(Float64,(0))
    photon_count :: spectrometer_visible__channel__filter_spectrometer__photon_count = spectrometer_visible__channel__filter_spectrometer__photon_count()
    output_voltage :: spectrometer_visible__channel__filter_spectrometer__output_voltage = spectrometer_visible__channel__filter_spectrometer__output_voltage()
    line_radiances :: spectrometer_visible__channel__filter_spectrometer__line_radiances = spectrometer_visible__channel__filter_spectrometer__line_radiances()
    exposure_time :: Float64 = 0.0
    radiance_calibration :: Float64 = 0.0
    photoelectric_voltage :: spectrometer_visible__channel__filter_spectrometer__photoelectric_voltage = spectrometer_visible__channel__filter_spectrometer__photoelectric_voltage()
    line_labels :: Array{String, 1} = String[]
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector_image
    circular :: spectrometer_visible__channel__detector_image__circular = spectrometer_visible__channel__detector_image__circular()
    outline :: spectrometer_visible__channel__detector_image__outline = spectrometer_visible__channel__detector_image__outline()
    geometry_type :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarization_spectroscopy
    temperature_cold_neutrals :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_lh_tor :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_hot_neutrals :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_cold_neutrals :: Array{Float64, 1} = zeros(Float64,(0))
    n_e :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_lh_r :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_modulus :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_lh_z :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_hot_neutrals :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__line_of_sight
    first_point :: spectrometer_visible__channel__line_of_sight__first_point = spectrometer_visible__channel__line_of_sight__first_point()
    second_point :: spectrometer_visible__channel__line_of_sight__second_point = spectrometer_visible__channel__line_of_sight__second_point()
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarizer__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__channel__polarizer
    surface :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    radius :: Float64 = 0.0
    centre :: spectrometer_visible__channel__polarizer__centre = spectrometer_visible__channel__polarizer__centre()
    x1_width :: Float64 = 0.0
    x2_unit_vector :: spectrometer_visible__channel__polarizer__x2_unit_vector = spectrometer_visible__channel__polarizer__x2_unit_vector()
    x3_unit_vector :: spectrometer_visible__channel__polarizer__x3_unit_vector = spectrometer_visible__channel__polarizer__x3_unit_vector()
    geometry_type :: Int32 = 0
    outline :: spectrometer_visible__channel__polarizer__outline = spectrometer_visible__channel__polarizer__outline()
    x1_unit_vector :: spectrometer_visible__channel__polarizer__x1_unit_vector = spectrometer_visible__channel__polarizer__x1_unit_vector()
end

Base.@kwdef mutable struct spectrometer_visible__channel__isotope_ratios
    method :: spectrometer_visible__channel__isotope_ratios__method = spectrometer_visible__channel__isotope_ratios__method()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    isotope :: StructArray{spectrometer_visible__channel__isotope_ratios__isotope} = StructArray(spectrometer_visible__channel__isotope_ratios__isotope() for k in 1:1)
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
    signal_to_noise :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__geometry_matrix__emission_grid__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__geometry_matrix__emission_grid
    grid_type :: spectrometer_visible__channel__geometry_matrix__emission_grid__grid_type = spectrometer_visible__channel__geometry_matrix__emission_grid__grid_type()
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
    dim3 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__geometry_matrix
    emission_grid :: spectrometer_visible__channel__geometry_matrix__emission_grid = spectrometer_visible__channel__geometry_matrix__emission_grid()
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    voxel_map :: Array{Int32, 3} = zeros(Int32,(0, 0, 0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_visible__ids_properties
    provider :: String = ""
    version_put :: spectrometer_visible__ids_properties__version_put = spectrometer_visible__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_visible__channel__grating_spectrometer
    processed_line :: StructArray{spectrometer_visible__channel__grating_spectrometer__processed_line} = StructArray(spectrometer_visible__channel__grating_spectrometer__processed_line() for k in 1:1)
    wavelength_calibration_date :: String = ""
    intensity_spectrum :: spectrometer_visible__channel__grating_spectrometer__intensity_spectrum = spectrometer_visible__channel__grating_spectrometer__intensity_spectrum()
    wavelength_calibration :: spectrometer_visible__channel__grating_spectrometer__wavelength_calibration = spectrometer_visible__channel__grating_spectrometer__wavelength_calibration()
    radiance_spectral :: spectrometer_visible__channel__grating_spectrometer__radiance_spectral = spectrometer_visible__channel__grating_spectrometer__radiance_spectral()
    radiance_calibration_date :: String = ""
    wavelengths :: Array{Float64, 1} = zeros(Float64,(0))
    slit_width :: Float64 = 0.0
    exposure_time :: Float64 = 0.0
    grating :: Float64 = 0.0
    radiance_calibration :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_visible__channel__detector
    surface :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: spectrometer_visible__channel__detector__centre = spectrometer_visible__channel__detector__centre()
    outline :: spectrometer_visible__channel__detector__outline = spectrometer_visible__channel__detector__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: spectrometer_visible__channel__detector__x3_unit_vector = spectrometer_visible__channel__detector__x3_unit_vector()
    x2_unit_vector :: spectrometer_visible__channel__detector__x2_unit_vector = spectrometer_visible__channel__detector__x2_unit_vector()
    geometry_type :: Int32 = 0
    radius :: Float64 = 0.0
    x1_unit_vector :: spectrometer_visible__channel__detector__x1_unit_vector = spectrometer_visible__channel__detector__x1_unit_vector()
end

Base.@kwdef mutable struct spectrometer_visible__channel
    polarizer :: spectrometer_visible__channel__polarizer = spectrometer_visible__channel__polarizer()
    object_observed :: String = ""
    polarization_spectroscopy :: spectrometer_visible__channel__polarization_spectroscopy = spectrometer_visible__channel__polarization_spectroscopy()
    validity :: Int32 = 0
    filter_spectrometer :: spectrometer_visible__channel__filter_spectrometer = spectrometer_visible__channel__filter_spectrometer()
    etendue :: Float64 = 0.0
    fibre_image :: spectrometer_visible__channel__fibre_image = spectrometer_visible__channel__fibre_image()
    detector :: spectrometer_visible__channel__detector = spectrometer_visible__channel__detector()
    geometry_matrix :: spectrometer_visible__channel__geometry_matrix = spectrometer_visible__channel__geometry_matrix()
    light_collection_efficiencies :: spectrometer_visible__channel__light_collection_efficiencies = spectrometer_visible__channel__light_collection_efficiencies()
    grating_spectrometer :: spectrometer_visible__channel__grating_spectrometer = spectrometer_visible__channel__grating_spectrometer()
    name :: String = ""
    polarizer_active :: Int32 = 0
    validity_timed :: spectrometer_visible__channel__validity_timed = spectrometer_visible__channel__validity_timed()
    active_spatial_resolution :: StructArray{spectrometer_visible__channel__active_spatial_resolution} = StructArray(spectrometer_visible__channel__active_spatial_resolution() for k in 1:1)
    etendue_method :: spectrometer_visible__channel__etendue_method = spectrometer_visible__channel__etendue_method()
    isotope_ratios :: spectrometer_visible__channel__isotope_ratios = spectrometer_visible__channel__isotope_ratios()
    line_of_sight :: spectrometer_visible__channel__line_of_sight = spectrometer_visible__channel__line_of_sight()
    detector_image :: spectrometer_visible__channel__detector_image = spectrometer_visible__channel__detector_image()
    aperture :: StructArray{spectrometer_visible__channel__aperture} = StructArray(spectrometer_visible__channel__aperture() for k in 1:1)
    type :: spectrometer_visible__channel__type = spectrometer_visible__channel__type()
end

Base.@kwdef mutable struct spectrometer_visible
    latency :: Float64 = 0.0
    channel :: StructArray{spectrometer_visible__channel} = StructArray(spectrometer_visible__channel() for k in 1:1)
    ids_properties :: spectrometer_visible__ids_properties = spectrometer_visible__ids_properties()
    detector_layout :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: spectrometer_visible__code = spectrometer_visible__code()
end

Base.@kwdef mutable struct spectrometer_uv__channel__processed_line__intensity
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__image_field__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__image_field__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__radiance_spectral
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__wavelength_calibration
    offset :: Float64 = 0.0
    gain :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__summit
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__line_of_sight__position_parameter
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__line_of_sight__moving_mode
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__validity_timed
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 2} = zeros(Int32,(0, 0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__ids_properties
    provider :: String = ""
    version_put :: spectrometer_uv__ids_properties__version_put = spectrometer_uv__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_uv__channel__intensity_spectrum
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_uv__channel__line_of_sight__second_point
    time :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__line_of_sight
    first_point :: spectrometer_uv__channel__line_of_sight__first_point = spectrometer_uv__channel__line_of_sight__first_point()
    moving_mode :: spectrometer_uv__channel__line_of_sight__moving_mode = spectrometer_uv__channel__line_of_sight__moving_mode()
    period :: Float64 = 0.0
    second_point :: spectrometer_uv__channel__line_of_sight__second_point = spectrometer_uv__channel__line_of_sight__second_point()
    amplitude_parameter :: Float64 = 0.0
    position_parameter :: spectrometer_uv__channel__line_of_sight__position_parameter = spectrometer_uv__channel__line_of_sight__position_parameter()
end

Base.@kwdef mutable struct spectrometer_uv__etendue_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__processed_line__radiance
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__processed_line
    label :: String = ""
    intensity :: spectrometer_uv__channel__processed_line__intensity = spectrometer_uv__channel__processed_line__intensity()
    radiance :: spectrometer_uv__channel__processed_line__radiance = spectrometer_uv__channel__processed_line__radiance()
    wavelength_central :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__image_field__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating__image_field
    curvature_radius :: Float64 = 0.0
    centre :: spectrometer_uv__channel__grating__image_field__centre = spectrometer_uv__channel__grating__image_field__centre()
    geometry_type :: spectrometer_uv__channel__grating__image_field__geometry_type = spectrometer_uv__channel__grating__image_field__geometry_type()
    x3_unit_vector :: spectrometer_uv__channel__grating__image_field__x3_unit_vector = spectrometer_uv__channel__grating__image_field__x3_unit_vector()
end

Base.@kwdef mutable struct spectrometer_uv__channel__grating
    curvature_radius :: Float64 = 0.0
    x2_unit_vector :: spectrometer_uv__channel__grating__x2_unit_vector = spectrometer_uv__channel__grating__x2_unit_vector()
    groove_density :: Float64 = 0.0
    x1_unit_vector :: spectrometer_uv__channel__grating__x1_unit_vector = spectrometer_uv__channel__grating__x1_unit_vector()
    summit :: spectrometer_uv__channel__grating__summit = spectrometer_uv__channel__grating__summit()
    outline :: spectrometer_uv__channel__grating__outline = spectrometer_uv__channel__grating__outline()
    centre :: spectrometer_uv__channel__grating__centre = spectrometer_uv__channel__grating__centre()
    geometry_type :: spectrometer_uv__channel__grating__geometry_type = spectrometer_uv__channel__grating__geometry_type()
    image_field :: spectrometer_uv__channel__grating__image_field = spectrometer_uv__channel__grating__image_field()
    x3_unit_vector :: spectrometer_uv__channel__grating__x3_unit_vector = spectrometer_uv__channel__grating__x3_unit_vector()
    type :: spectrometer_uv__channel__grating__type = spectrometer_uv__channel__grating__type()
end

Base.@kwdef mutable struct spectrometer_uv__channel__supply_high_voltage__voltage_set
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__supply_high_voltage
    voltage_set :: spectrometer_uv__channel__supply_high_voltage__voltage_set = spectrometer_uv__channel__supply_high_voltage__voltage_set()
    object :: String = ""
end

Base.@kwdef mutable struct spectrometer_uv__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__aperture
    x2_unit_vector :: spectrometer_uv__channel__aperture__x2_unit_vector = spectrometer_uv__channel__aperture__x2_unit_vector()
    x2_width :: Float64 = 0.0
    centre :: spectrometer_uv__channel__aperture__centre = spectrometer_uv__channel__aperture__centre()
    outline :: spectrometer_uv__channel__aperture__outline = spectrometer_uv__channel__aperture__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: spectrometer_uv__channel__aperture__x3_unit_vector = spectrometer_uv__channel__aperture__x3_unit_vector()
    surface :: Float64 = 0.0
    geometry_type :: Int32 = 0
    radius :: Float64 = 0.0
    x1_unit_vector :: spectrometer_uv__channel__aperture__x1_unit_vector = spectrometer_uv__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector_position_parameter
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__code
    library :: StructArray{spectrometer_uv__code__library} = StructArray(spectrometer_uv__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector
    x1_width :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    radius :: Float64 = 0.0
    outline :: spectrometer_uv__channel__detector__outline = spectrometer_uv__channel__detector__outline()
    centre :: spectrometer_uv__channel__detector__centre = spectrometer_uv__channel__detector__centre()
    x2_unit_vector :: spectrometer_uv__channel__detector__x2_unit_vector = spectrometer_uv__channel__detector__x2_unit_vector()
    surface :: Float64 = 0.0
    x3_unit_vector :: spectrometer_uv__channel__detector__x3_unit_vector = spectrometer_uv__channel__detector__x3_unit_vector()
    geometry_type :: Int32 = 0
    x1_unit_vector :: spectrometer_uv__channel__detector__x1_unit_vector = spectrometer_uv__channel__detector__x1_unit_vector()
end

Base.@kwdef mutable struct spectrometer_uv__channel__detector_layout
    detector_dimensions :: Array{Float64, 1} = zeros(Float64,(0))
    pixel_n :: Array{Int32, 1} = zeros(Int32,(0))
    pixel_dimensions :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv__channel
    processed_line :: StructArray{spectrometer_uv__channel__processed_line} = StructArray(spectrometer_uv__channel__processed_line() for k in 1:1)
    validity :: Int32 = 0
    wavelength_calibration_date :: String = ""
    detector_position_parameter :: spectrometer_uv__channel__detector_position_parameter = spectrometer_uv__channel__detector_position_parameter()
    detector :: spectrometer_uv__channel__detector = spectrometer_uv__channel__detector()
    validity_timed :: spectrometer_uv__channel__validity_timed = spectrometer_uv__channel__validity_timed()
    name :: String = ""
    intensity_spectrum :: spectrometer_uv__channel__intensity_spectrum = spectrometer_uv__channel__intensity_spectrum()
    wavelength_calibration :: spectrometer_uv__channel__wavelength_calibration = spectrometer_uv__channel__wavelength_calibration()
    radiance_spectral :: spectrometer_uv__channel__radiance_spectral = spectrometer_uv__channel__radiance_spectral()
    radiance_calibration_date :: String = ""
    wavelengths :: Array{Float64, 1} = zeros(Float64,(0))
    supply_high_voltage :: StructArray{spectrometer_uv__channel__supply_high_voltage} = StructArray(spectrometer_uv__channel__supply_high_voltage() for k in 1:1)
    exposure_time :: Float64 = 0.0
    detector_layout :: spectrometer_uv__channel__detector_layout = spectrometer_uv__channel__detector_layout()
    aperture :: StructArray{spectrometer_uv__channel__aperture} = StructArray(spectrometer_uv__channel__aperture() for k in 1:1)
    line_of_sight :: spectrometer_uv__channel__line_of_sight = spectrometer_uv__channel__line_of_sight()
    grating :: spectrometer_uv__channel__grating = spectrometer_uv__channel__grating()
    radiance_calibration :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct spectrometer_uv
    latency :: Float64 = 0.0
    channel :: StructArray{spectrometer_uv__channel} = StructArray(spectrometer_uv__channel() for k in 1:1)
    ids_properties :: spectrometer_uv__ids_properties = spectrometer_uv__ids_properties()
    etendue :: Float64 = 0.0
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: spectrometer_uv__code = spectrometer_uv__code()
    etendue_method :: spectrometer_uv__etendue_method = spectrometer_uv__etendue_method()
end

Base.@kwdef mutable struct spectrometer_mass__pressures_partial
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_mass__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct spectrometer_mass__code
    library :: StructArray{spectrometer_mass__code__library} = StructArray(spectrometer_mass__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct spectrometer_mass__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct spectrometer_mass__ids_properties
    provider :: String = ""
    version_put :: spectrometer_mass__ids_properties__version_put = spectrometer_mass__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct spectrometer_mass
    name :: String = ""
    pressures_partial :: spectrometer_mass__pressures_partial = spectrometer_mass__pressures_partial()
    latency :: Float64 = 0.0
    ids_properties :: spectrometer_mass__ids_properties = spectrometer_mass__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: spectrometer_mass__code = spectrometer_mass__code()
    identifier :: String = ""
    a :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct soft_x_rays__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct soft_x_rays__channel__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct soft_x_rays__channel__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__etendue_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct soft_x_rays__ids_properties
    provider :: String = ""
    version_put :: soft_x_rays__ids_properties__version_put = soft_x_rays__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct soft_x_rays__channel__brightness
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct soft_x_rays__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__energy_band
    detection_efficiency :: Array{Float64, 1} = zeros(Float64,(0))
    lower_bound :: Float64 = 0.0
    energies :: Array{Float64, 1} = zeros(Float64,(0))
    upper_bound :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct soft_x_rays__code
    library :: StructArray{soft_x_rays__code__library} = StructArray(soft_x_rays__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct soft_x_rays__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__power
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct soft_x_rays__channel__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__validity_timed
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct soft_x_rays__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__line_of_sight
    first_point :: soft_x_rays__channel__line_of_sight__first_point = soft_x_rays__channel__line_of_sight__first_point()
    second_point :: soft_x_rays__channel__line_of_sight__second_point = soft_x_rays__channel__line_of_sight__second_point()
end

Base.@kwdef mutable struct soft_x_rays__channel__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct soft_x_rays__channel__detector
    geometry_type :: Int32 = 0
    x2_width :: Float64 = 0.0
    centre :: soft_x_rays__channel__detector__centre = soft_x_rays__channel__detector__centre()
    outline :: soft_x_rays__channel__detector__outline = soft_x_rays__channel__detector__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: soft_x_rays__channel__detector__x3_unit_vector = soft_x_rays__channel__detector__x3_unit_vector()
    surface :: Float64 = 0.0
    x2_unit_vector :: soft_x_rays__channel__detector__x2_unit_vector = soft_x_rays__channel__detector__x2_unit_vector()
    radius :: Float64 = 0.0
    x1_unit_vector :: soft_x_rays__channel__detector__x1_unit_vector = soft_x_rays__channel__detector__x1_unit_vector()
end

Base.@kwdef mutable struct soft_x_rays__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct soft_x_rays__channel__aperture
    geometry_type :: Int32 = 0
    x2_width :: Float64 = 0.0
    centre :: soft_x_rays__channel__aperture__centre = soft_x_rays__channel__aperture__centre()
    outline :: soft_x_rays__channel__aperture__outline = soft_x_rays__channel__aperture__outline()
    x1_width :: Float64 = 0.0
    x2_unit_vector :: soft_x_rays__channel__aperture__x2_unit_vector = soft_x_rays__channel__aperture__x2_unit_vector()
    x3_unit_vector :: soft_x_rays__channel__aperture__x3_unit_vector = soft_x_rays__channel__aperture__x3_unit_vector()
    surface :: Float64 = 0.0
    radius :: Float64 = 0.0
    x1_unit_vector :: soft_x_rays__channel__aperture__x1_unit_vector = soft_x_rays__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct soft_x_rays__channel
    brightness :: soft_x_rays__channel__brightness = soft_x_rays__channel__brightness()
    validity :: Int32 = 0
    etendue :: Float64 = 0.0
    energy_band :: StructArray{soft_x_rays__channel__energy_band} = StructArray(soft_x_rays__channel__energy_band() for k in 1:1)
    validity_timed :: soft_x_rays__channel__validity_timed = soft_x_rays__channel__validity_timed()
    detector :: soft_x_rays__channel__detector = soft_x_rays__channel__detector()
    name :: String = ""
    etendue_method :: soft_x_rays__channel__etendue_method = soft_x_rays__channel__etendue_method()
    power :: soft_x_rays__channel__power = soft_x_rays__channel__power()
    aperture :: StructArray{soft_x_rays__channel__aperture} = StructArray(soft_x_rays__channel__aperture() for k in 1:1)
    line_of_sight :: soft_x_rays__channel__line_of_sight = soft_x_rays__channel__line_of_sight()
    identifier :: String = ""
end

Base.@kwdef mutable struct soft_x_rays
    time :: Array{Float64, 1} = zeros(Float64,(0))
    channel :: StructArray{soft_x_rays__channel} = StructArray(soft_x_rays__channel() for k in 1:1)
    ids_properties :: soft_x_rays__ids_properties = soft_x_rays__ids_properties()
    latency :: Float64 = 0.0
    code :: soft_x_rays__code = soft_x_rays__code()
end

Base.@kwdef mutable struct sdn__topic__signal__value
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct sdn__topic__signal__quality
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct sdn__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct sdn__ids_properties
    provider :: String = ""
    version_put :: sdn__ids_properties__version_put = sdn__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct sdn__topic__signal
    allocated_position :: Int32 = 0
    name :: String = ""
    value :: sdn__topic__signal__value = sdn__topic__signal__value()
    quality :: sdn__topic__signal__quality = sdn__topic__signal__quality()
    definition :: String = ""
end

Base.@kwdef mutable struct sdn__topic
    name :: String = ""
    signal :: StructArray{sdn__topic__signal} = StructArray(sdn__topic__signal() for k in 1:1)
end

Base.@kwdef mutable struct sdn__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct sdn__code
    library :: StructArray{sdn__code__library} = StructArray(sdn__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct sdn
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: sdn__ids_properties = sdn__ids_properties()
    topic :: StructArray{sdn__topic} = StructArray(sdn__topic() for k in 1:1)
    code :: sdn__code = sdn__code()
end

Base.@kwdef mutable struct sawteeth__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct sawteeth__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct sawteeth__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct sawteeth__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct sawteeth__profiles_1d
    time :: Float64 = 0.0
    p_i_total_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    psi_star_post_crash :: Array{Float64, 1} = zeros(Float64,(0))
    p_e_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    magnetic_shear :: Array{Float64, 1} = zeros(Float64,(0))
    j_total :: Array{Float64, 1} = zeros(Float64,(0))
    n_e_fast :: Array{Float64, 1} = zeros(Float64,(0))
    t_e :: Array{Float64, 1} = zeros(Float64,(0))
    j_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    t_i_average :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    j_non_inductive :: Array{Float64, 1} = zeros(Float64,(0))
    conductivity_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    q :: Array{Float64, 1} = zeros(Float64,(0))
    p_i_total :: Array{Float64, 1} = zeros(Float64,(0))
    p_e :: Array{Float64, 1} = zeros(Float64,(0))
    j_ohmic :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: sawteeth__profiles_1d__grid = sawteeth__profiles_1d__grid()
    j_bootstrap :: Array{Float64, 1} = zeros(Float64,(0))
    n_e :: Array{Float64, 1} = zeros(Float64,(0))
    p_i_total_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    p_e_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    psi_star_pre_crash :: Array{Float64, 1} = zeros(Float64,(0))
    zeff :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct sawteeth__code
    library :: StructArray{sawteeth__code__library} = StructArray(sawteeth__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct sawteeth__diagnostics
    previous_period :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm_q1 :: Array{Float64, 1} = zeros(Float64,(0))
    previous_crash_trigger :: Array{Int32, 1} = zeros(Int32,(0))
    rho_tor_norm_mixing :: Array{Float64, 1} = zeros(Float64,(0))
    previous_crash_time :: Array{Float64, 1} = zeros(Float64,(0))
    magnetic_shear_q1 :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm_inversion :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct sawteeth__ids_properties
    provider :: String = ""
    version_put :: sawteeth__ids_properties__version_put = sawteeth__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct sawteeth
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: sawteeth__ids_properties = sawteeth__ids_properties()
    vacuum_toroidal_field :: sawteeth__vacuum_toroidal_field = sawteeth__vacuum_toroidal_field()
    diagnostics :: sawteeth__diagnostics = sawteeth__diagnostics()
    code :: sawteeth__code = sawteeth__code()
    profiles_1d :: StructArray{sawteeth__profiles_1d} = StructArray(sawteeth__profiles_1d() for k in 1:1)
    crash_trigger :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct refractometer__channel__n_e_line
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct refractometer__channel__bandwidth__phase_quadrature
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct refractometer__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct refractometer__channel__bandwidth__n_e_line
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct refractometer__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct refractometer__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct refractometer__ids_properties
    provider :: String = ""
    version_put :: refractometer__ids_properties__version_put = refractometer__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct refractometer__channel__bandwidth
    time_detector :: Array{Float64, 2} = zeros(Float64,(0,0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
    n_e_line :: refractometer__channel__bandwidth__n_e_line = refractometer__channel__bandwidth__n_e_line()
    phase_quadrature :: refractometer__channel__bandwidth__phase_quadrature = refractometer__channel__bandwidth__phase_quadrature()
    phase :: Array{Float64, 1} = zeros(Float64,(0))
    q_component :: Array{Float64, 2} = zeros(Float64,(0,0))
    frequency_main :: Float64 = 0.0
    i_component :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct refractometer__channel__n_e_profile_approximation__formula
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct refractometer__channel__n_e_profile_approximation
    parameters :: Array{Float64, 2} = zeros(Float64,(0,0))
    formula :: refractometer__channel__n_e_profile_approximation__formula = refractometer__channel__n_e_profile_approximation__formula()
end

Base.@kwdef mutable struct refractometer__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct refractometer__channel__line_of_sight
    first_point :: refractometer__channel__line_of_sight__first_point = refractometer__channel__line_of_sight__first_point()
    second_point :: refractometer__channel__line_of_sight__second_point = refractometer__channel__line_of_sight__second_point()
end

Base.@kwdef mutable struct refractometer__channel
    n_e_profile_approximation :: refractometer__channel__n_e_profile_approximation = refractometer__channel__n_e_profile_approximation()
    name :: String = ""
    bandwidth :: StructArray{refractometer__channel__bandwidth} = StructArray(refractometer__channel__bandwidth() for k in 1:1)
    n_e_line :: refractometer__channel__n_e_line = refractometer__channel__n_e_line()
    line_of_sight :: refractometer__channel__line_of_sight = refractometer__channel__line_of_sight()
    mode :: String = ""
    identifier :: String = ""
end

Base.@kwdef mutable struct refractometer__code
    library :: StructArray{refractometer__code__library} = StructArray(refractometer__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct refractometer
    latency :: Float64 = 0.0
    channel :: StructArray{refractometer__channel} = StructArray(refractometer__channel() for k in 1:1)
    ids_properties :: refractometer__ids_properties = refractometer__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: refractometer__code = refractometer__code()
    type :: String = ""
end

Base.@kwdef mutable struct reflectometer_profile__channel__phase
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct reflectometer_profile__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct reflectometer_profile__position
    psi :: Array{Float64, 2} = zeros(Float64,(0,0))
    theta :: Array{Float64, 2} = zeros(Float64,(0,0))
    phi :: Array{Float64, 2} = zeros(Float64,(0,0))
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
    rho_tor_norm :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct reflectometer_profile__channel__line_of_sight_detection__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct reflectometer_profile__channel__line_of_sight_detection__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct reflectometer_profile__channel__line_of_sight_detection
    first_point :: reflectometer_profile__channel__line_of_sight_detection__first_point = reflectometer_profile__channel__line_of_sight_detection__first_point()
    second_point :: reflectometer_profile__channel__line_of_sight_detection__second_point = reflectometer_profile__channel__line_of_sight_detection__second_point()
end

Base.@kwdef mutable struct reflectometer_profile__n_e
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct reflectometer_profile__channel__line_of_sight_emission__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct reflectometer_profile__channel__n_e
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct reflectometer_profile__channel__line_of_sight_emission__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct reflectometer_profile__channel__line_of_sight_emission
    first_point :: reflectometer_profile__channel__line_of_sight_emission__first_point = reflectometer_profile__channel__line_of_sight_emission__first_point()
    second_point :: reflectometer_profile__channel__line_of_sight_emission__second_point = reflectometer_profile__channel__line_of_sight_emission__second_point()
end

Base.@kwdef mutable struct reflectometer_profile__ids_properties
    provider :: String = ""
    version_put :: reflectometer_profile__ids_properties__version_put = reflectometer_profile__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct reflectometer_profile__channel__position
    psi :: Array{Float64, 2} = zeros(Float64,(0,0))
    theta :: Array{Float64, 2} = zeros(Float64,(0,0))
    phi :: Array{Float64, 2} = zeros(Float64,(0,0))
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
    rho_tor_norm :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct reflectometer_profile__channel
    cut_off_frequency :: Array{Float64, 2} = zeros(Float64,(0,0))
    mode :: String = ""
    line_of_sight_detection :: reflectometer_profile__channel__line_of_sight_detection = reflectometer_profile__channel__line_of_sight_detection()
    sweep_time :: Float64 = 0.0
    name :: String = ""
    position :: reflectometer_profile__channel__position = reflectometer_profile__channel__position()
    frequencies :: Array{Float64, 1} = zeros(Float64,(0))
    n_e :: reflectometer_profile__channel__n_e = reflectometer_profile__channel__n_e()
    phase :: reflectometer_profile__channel__phase = reflectometer_profile__channel__phase()
    identifier :: String = ""
    line_of_sight_emission :: reflectometer_profile__channel__line_of_sight_emission = reflectometer_profile__channel__line_of_sight_emission()
end

Base.@kwdef mutable struct reflectometer_profile__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct reflectometer_profile__psi_normalization
    psi_boundary :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct reflectometer_profile__code
    library :: StructArray{reflectometer_profile__code__library} = StructArray(reflectometer_profile__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct reflectometer_profile
    psi_normalization :: reflectometer_profile__psi_normalization = reflectometer_profile__psi_normalization()
    latency :: Float64 = 0.0
    channel :: StructArray{reflectometer_profile__channel} = StructArray(reflectometer_profile__channel() for k in 1:1)
    ids_properties :: reflectometer_profile__ids_properties = reflectometer_profile__ids_properties()
    n_e :: reflectometer_profile__n_e = reflectometer_profile__n_e()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    position :: reflectometer_profile__position = reflectometer_profile__position()
    code :: reflectometer_profile__code = reflectometer_profile__code()
    type :: String = ""
end

Base.@kwdef mutable struct radiation__process__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__process__profiles_1d__electrons
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
    emissivity :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct radiation__grid_ggd__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__grid_ggd__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct radiation__process__profiles_1d__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__grid_ggd__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__ggd__electrons__emissivity
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct radiation__process__ggd__electrons
    emissivity :: StructArray{radiation__process__ggd__electrons__emissivity} = StructArray(radiation__process__ggd__electrons__emissivity() for k in 1:1)
end

Base.@kwdef mutable struct radiation__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__process__profiles_1d__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    z_min :: Float64 = 0.0
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
    emissivity :: Array{Float64, 1} = zeros(Float64,(0))
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__process__ggd__neutral__state__emissivity
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct radiation__grid_ggd__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__ggd__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct radiation__process__global_quantities__inside_vessel
    power_electrons :: Float64 = 0.0
    power_neutral_total :: Float64 = 0.0
    power_ion_total :: Float64 = 0.0
    power :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__grid_ggd__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__global_quantities__inside_lcfs
    power_electrons :: Float64 = 0.0
    power_neutral_total :: Float64 = 0.0
    power_ion_total :: Float64 = 0.0
    power :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__process__global_quantities
    time :: Float64 = 0.0
    inside_vessel :: radiation__process__global_quantities__inside_vessel = radiation__process__global_quantities__inside_vessel()
    inside_lcfs :: radiation__process__global_quantities__inside_lcfs = radiation__process__global_quantities__inside_lcfs()
end

Base.@kwdef mutable struct radiation__grid_ggd__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__ggd__neutral__emissivity
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct radiation__process__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__process__ggd__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__ids_properties
    provider :: String = ""
    version_put :: radiation__ids_properties__version_put = radiation__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__ggd__ion__emissivity
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct radiation__process__ggd__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__process__ggd__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    emissivity :: StructArray{radiation__process__ggd__neutral__state__emissivity} = StructArray(radiation__process__ggd__neutral__state__emissivity() for k in 1:1)
    vibrational_mode :: String = ""
    neutral_type :: radiation__process__ggd__neutral__state__neutral_type = radiation__process__ggd__neutral__state__neutral_type()
end

Base.@kwdef mutable struct radiation__process__ggd__neutral
    label :: String = ""
    emissivity :: StructArray{radiation__process__ggd__neutral__emissivity} = StructArray(radiation__process__ggd__neutral__emissivity() for k in 1:1)
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    state :: StructArray{radiation__process__ggd__neutral__state} = StructArray(radiation__process__ggd__neutral__state() for k in 1:1)
    element :: StructArray{radiation__process__ggd__neutral__element} = StructArray(radiation__process__ggd__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct radiation__grid_ggd__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct radiation__grid_ggd__grid_subset__element
    object :: StructArray{radiation__grid_ggd__grid_subset__element__object} = StructArray(radiation__grid_ggd__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct radiation__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct radiation__code
    library :: StructArray{radiation__code__library} = StructArray(radiation__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct radiation__process__profiles_1d__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    multiple_states_flag :: Int32 = 0
    state :: StructArray{radiation__process__profiles_1d__ion__state} = StructArray(radiation__process__profiles_1d__ion__state() for k in 1:1)
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    emissivity :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{radiation__process__profiles_1d__ion__element} = StructArray(radiation__process__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct radiation__process__profiles_1d__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    neutral_type :: radiation__process__profiles_1d__neutral__state__neutral_type = radiation__process__profiles_1d__neutral__state__neutral_type()
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
    emissivity :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct radiation__process__profiles_1d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__grid_ggd__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct radiation__grid_ggd__grid_subset
    base :: StructArray{radiation__grid_ggd__grid_subset__base} = StructArray(radiation__grid_ggd__grid_subset__base() for k in 1:1)
    metric :: radiation__grid_ggd__grid_subset__metric = radiation__grid_ggd__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: radiation__grid_ggd__grid_subset__identifier = radiation__grid_ggd__grid_subset__identifier()
    element :: StructArray{radiation__grid_ggd__grid_subset__element} = StructArray(radiation__grid_ggd__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct radiation__process__ggd__ion__state__emissivity
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct radiation__process__ggd__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    emissivity :: StructArray{radiation__process__ggd__ion__state__emissivity} = StructArray(radiation__process__ggd__ion__state__emissivity() for k in 1:1)
    vibrational_mode :: String = ""
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct radiation__process__ggd__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    emissivity :: StructArray{radiation__process__ggd__ion__emissivity} = StructArray(radiation__process__ggd__ion__emissivity() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    state :: StructArray{radiation__process__ggd__ion__state} = StructArray(radiation__process__ggd__ion__state() for k in 1:1)
    neutral_index :: Int32 = 0
    element :: StructArray{radiation__process__ggd__ion__element} = StructArray(radiation__process__ggd__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct radiation__grid_ggd__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{radiation__grid_ggd__space__objects_per_dimension__object__boundary} = StructArray(radiation__grid_ggd__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct radiation__grid_ggd__space__objects_per_dimension
    object :: StructArray{radiation__grid_ggd__space__objects_per_dimension__object} = StructArray(radiation__grid_ggd__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct radiation__grid_ggd__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{radiation__grid_ggd__space__objects_per_dimension} = StructArray(radiation__grid_ggd__space__objects_per_dimension() for k in 1:1)
    geometry_type :: radiation__grid_ggd__space__geometry_type = radiation__grid_ggd__space__geometry_type()
    identifier :: radiation__grid_ggd__space__identifier = radiation__grid_ggd__space__identifier()
end

Base.@kwdef mutable struct radiation__process__ggd
    ion :: StructArray{radiation__process__ggd__ion} = StructArray(radiation__process__ggd__ion() for k in 1:1)
    time :: Float64 = 0.0
    neutral :: StructArray{radiation__process__ggd__neutral} = StructArray(radiation__process__ggd__neutral() for k in 1:1)
    electrons :: radiation__process__ggd__electrons = radiation__process__ggd__electrons()
end

Base.@kwdef mutable struct radiation__process__profiles_1d__neutral
    label :: String = ""
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    state :: StructArray{radiation__process__profiles_1d__neutral__state} = StructArray(radiation__process__profiles_1d__neutral__state() for k in 1:1)
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
    emissivity :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{radiation__process__profiles_1d__neutral__element} = StructArray(radiation__process__profiles_1d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct radiation__process__profiles_1d
    ion :: StructArray{radiation__process__profiles_1d__ion} = StructArray(radiation__process__profiles_1d__ion() for k in 1:1)
    emissivity_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    power_inside_neutral_total :: Array{Float64, 1} = zeros(Float64,(0))
    neutral :: StructArray{radiation__process__profiles_1d__neutral} = StructArray(radiation__process__profiles_1d__neutral() for k in 1:1)
    power_inside_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: radiation__process__profiles_1d__grid = radiation__process__profiles_1d__grid()
    emissivity_neutral_total :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: radiation__process__profiles_1d__electrons = radiation__process__profiles_1d__electrons()
end

Base.@kwdef mutable struct radiation__process
    global_quantities :: StructArray{radiation__process__global_quantities} = StructArray(radiation__process__global_quantities() for k in 1:1)
    profiles_1d :: StructArray{radiation__process__profiles_1d} = StructArray(radiation__process__profiles_1d() for k in 1:1)
    identifier :: radiation__process__identifier = radiation__process__identifier()
    ggd :: StructArray{radiation__process__ggd} = StructArray(radiation__process__ggd() for k in 1:1)
end

Base.@kwdef mutable struct radiation__grid_ggd
    time :: Float64 = 0.0
    grid_subset :: StructArray{radiation__grid_ggd__grid_subset} = StructArray(radiation__grid_ggd__grid_subset() for k in 1:1)
    space :: StructArray{radiation__grid_ggd__space} = StructArray(radiation__grid_ggd__space() for k in 1:1)
    identifier :: radiation__grid_ggd__identifier = radiation__grid_ggd__identifier()
end

Base.@kwdef mutable struct radiation
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: radiation__ids_properties = radiation__ids_properties()
    process :: StructArray{radiation__process} = StructArray(radiation__process() for k in 1:1)
    vacuum_toroidal_field :: radiation__vacuum_toroidal_field = radiation__vacuum_toroidal_field()
    grid_ggd :: StructArray{radiation__grid_ggd} = StructArray(radiation__grid_ggd() for k in 1:1)
    code :: radiation__code = radiation__code()
end

Base.@kwdef mutable struct pulse_schedule__position_control__boundary_outline__z__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__valve__flow_rate__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__geometric_axis__z__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__power__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__triangularity_lower__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__frequency__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__power__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__elongation_upper__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__nbi__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__magnetic_axis__z__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__gap__value__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__deposition_rho_tor_norm__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__power
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ic__antenna__power__reference = pulse_schedule__ic__antenna__power__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__flux_control__i_plasma__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__zeff__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__power__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__triangularity__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__active_limiter_point__r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__event__acquisition_state
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__power__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__active_limiter_point__z__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__x_point__z__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__power
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__nbi__unit__power__reference = pulse_schedule__nbi__unit__power__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__triangularity_lower
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__triangularity_lower__reference = pulse_schedule__position_control__triangularity_lower__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__tf__b_field_tor_vacuum_r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__tf__b_field_tor_vacuum_r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__tf__b_field_tor_vacuum_r__reference = pulse_schedule__tf__b_field_tor_vacuum_r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__phase__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__active_limiter_point__z
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__active_limiter_point__z__reference = pulse_schedule__position_control__active_limiter_point__z__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__frequency__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__frequency
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ic__antenna__frequency__reference = pulse_schedule__ic__antenna__frequency__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__phase__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__n_e_line__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__ic__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__strike_point__z__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__triangularity_upper__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__x_point__z
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__x_point__z__reference = pulse_schedule__position_control__x_point__z__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__magnetic_axis__r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__magnetic_axis__r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__magnetic_axis__r__reference = pulse_schedule__position_control__magnetic_axis__r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__frequency__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__x_point__r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__x_point__r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__x_point__r__reference = pulse_schedule__position_control__x_point__r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__x_point
    r :: pulse_schedule__position_control__x_point__r = pulse_schedule__position_control__x_point__r()
    z :: pulse_schedule__position_control__x_point__z = pulse_schedule__position_control__x_point__z()
end

Base.@kwdef mutable struct pulse_schedule__tf__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__tf
    mode :: pulse_schedule__tf__mode = pulse_schedule__tf__mode()
    b_field_tor_vacuum_r :: pulse_schedule__tf__b_field_tor_vacuum_r = pulse_schedule__tf__b_field_tor_vacuum_r()
end

Base.@kwdef mutable struct pulse_schedule__position_control__boundary_outline__z
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__boundary_outline__z__reference = pulse_schedule__position_control__boundary_outline__z__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__geometric_axis__r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__deposition_rho_tor_norm
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ec__launcher__deposition_rho_tor_norm__reference = pulse_schedule__ec__launcher__deposition_rho_tor_norm__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__phase
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__lh__antenna__phase__reference = pulse_schedule__lh__antenna__phase__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__n_parallel__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__flux_control__beta_normal__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__energy__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__energy
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__nbi__unit__energy__reference = pulse_schedule__nbi__unit__energy__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__power_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__power
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ec__launcher__power__reference = pulse_schedule__ec__launcher__power__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__triangularity_upper
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__triangularity_upper__reference = pulse_schedule__position_control__triangularity_upper__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__density_control__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__valve__species__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct pulse_schedule__density_control__valve__species
    label :: String = ""
    fraction :: Float64 = 0.0
    element :: StructArray{pulse_schedule__density_control__valve__species__element} = StructArray(pulse_schedule__density_control__valve__species__element() for k in 1:1)
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__steering_angle_pol__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__n_h_over_n_d__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__n_h_over_n_d
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__density_control__n_h_over_n_d__reference = pulse_schedule__density_control__n_h_over_n_d__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__density_control__n_t_over_n_d__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__density_control__n_t_over_n_d
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__density_control__n_t_over_n_d__reference = pulse_schedule__density_control__n_t_over_n_d__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ec__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__flux_control__li_3__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__flux_control__li_3
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__flux_control__li_3__reference = pulse_schedule__flux_control__li_3__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__density_control__zeff
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__density_control__zeff__reference = pulse_schedule__density_control__zeff__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__minor_radius__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__minor_radius
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__minor_radius__reference = pulse_schedule__position_control__minor_radius__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__gap__value
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__gap__value__reference = pulse_schedule__position_control__gap__value__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__ids_properties
    provider :: String = ""
    version_put :: pulse_schedule__ids_properties__version_put = pulse_schedule__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__steering_angle_pol
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ec__launcher__steering_angle_pol__reference = pulse_schedule__ec__launcher__steering_angle_pol__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__geometric_axis__z
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__geometric_axis__z__reference = pulse_schedule__position_control__geometric_axis__z__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__event__acquisition_strategy
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__elongation__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__power_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__magnetic_axis__z
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__magnetic_axis__z__reference = pulse_schedule__position_control__magnetic_axis__z__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__frequency
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__lh__antenna__frequency__reference = pulse_schedule__lh__antenna__frequency__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__geometric_axis__r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__geometric_axis__r__reference = pulse_schedule__position_control__geometric_axis__r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__geometric_axis
    r :: pulse_schedule__position_control__geometric_axis__r = pulse_schedule__position_control__geometric_axis__r()
    z :: pulse_schedule__position_control__geometric_axis__z = pulse_schedule__position_control__geometric_axis__z()
end

Base.@kwdef mutable struct pulse_schedule__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__frequency
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ec__launcher__frequency__reference = pulse_schedule__ec__launcher__frequency__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__event__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__event
    duration :: Float64 = 0.0
    provider :: String = ""
    time_stamp :: Float64 = 0.0
    acquisition_state :: pulse_schedule__event__acquisition_state = pulse_schedule__event__acquisition_state()
    listeners :: Array{String, 1} = String[]
    identifier :: String = ""
    type :: pulse_schedule__event__type = pulse_schedule__event__type()
    acquisition_strategy :: pulse_schedule__event__acquisition_strategy = pulse_schedule__event__acquisition_strategy()
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna__phase
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ic__antenna__phase__reference = pulse_schedule__ic__antenna__phase__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__density_control__valve__flow_rate
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__density_control__valve__flow_rate__reference = pulse_schedule__density_control__valve__flow_rate__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__density_control__valve
    name :: String = ""
    flow_rate :: pulse_schedule__density_control__valve__flow_rate = pulse_schedule__density_control__valve__flow_rate()
    species :: StructArray{pulse_schedule__density_control__valve__species} = StructArray(pulse_schedule__density_control__valve__species() for k in 1:1)
    identifier :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__position_control__boundary_outline__r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__boundary_outline__r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__boundary_outline__r__reference = pulse_schedule__position_control__boundary_outline__r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__boundary_outline
    r :: pulse_schedule__position_control__boundary_outline__r = pulse_schedule__position_control__boundary_outline__r()
    z :: pulse_schedule__position_control__boundary_outline__z = pulse_schedule__position_control__boundary_outline__z()
end

Base.@kwdef mutable struct pulse_schedule__position_control__triangularity
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__triangularity__reference = pulse_schedule__position_control__triangularity__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__steering_angle_tor__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__steering_angle_tor
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__ec__launcher__steering_angle_tor__reference = pulse_schedule__ec__launcher__steering_angle_tor__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__species__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit__species
    label :: String = ""
    fraction :: Float64 = 0.0
    element :: StructArray{pulse_schedule__nbi__unit__species__element} = StructArray(pulse_schedule__nbi__unit__species__element() for k in 1:1)
end

Base.@kwdef mutable struct pulse_schedule__nbi__unit
    name :: String = ""
    power_type :: pulse_schedule__nbi__unit__power_type = pulse_schedule__nbi__unit__power_type()
    energy :: pulse_schedule__nbi__unit__energy = pulse_schedule__nbi__unit__energy()
    species :: StructArray{pulse_schedule__nbi__unit__species} = StructArray(pulse_schedule__nbi__unit__species() for k in 1:1)
    power :: pulse_schedule__nbi__unit__power = pulse_schedule__nbi__unit__power()
    identifier :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__nbi
    mode :: pulse_schedule__nbi__mode = pulse_schedule__nbi__mode()
    unit :: StructArray{pulse_schedule__nbi__unit} = StructArray(pulse_schedule__nbi__unit() for k in 1:1)
end

Base.@kwdef mutable struct pulse_schedule__position_control__elongation
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__elongation__reference = pulse_schedule__position_control__elongation__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__power
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__lh__antenna__power__reference = pulse_schedule__lh__antenna__power__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__magnetic_axis
    r :: pulse_schedule__position_control__magnetic_axis__r = pulse_schedule__position_control__magnetic_axis__r()
    z :: pulse_schedule__position_control__magnetic_axis__z = pulse_schedule__position_control__magnetic_axis__z()
end

Base.@kwdef mutable struct pulse_schedule__flux_control__i_plasma
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__flux_control__i_plasma__reference = pulse_schedule__flux_control__i_plasma__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__strike_point__z
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__strike_point__z__reference = pulse_schedule__position_control__strike_point__z__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__elongation_upper
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__elongation_upper__reference = pulse_schedule__position_control__elongation_upper__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__elongation_lower__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__elongation_lower
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__elongation_lower__reference = pulse_schedule__position_control__elongation_lower__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__power_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna__n_parallel
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__lh__antenna__n_parallel__reference = pulse_schedule__lh__antenna__n_parallel__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__lh__antenna
    name :: String = ""
    power_type :: pulse_schedule__lh__antenna__power_type = pulse_schedule__lh__antenna__power_type()
    phase :: pulse_schedule__lh__antenna__phase = pulse_schedule__lh__antenna__phase()
    power :: pulse_schedule__lh__antenna__power = pulse_schedule__lh__antenna__power()
    frequency :: pulse_schedule__lh__antenna__frequency = pulse_schedule__lh__antenna__frequency()
    identifier :: String = ""
    n_parallel :: pulse_schedule__lh__antenna__n_parallel = pulse_schedule__lh__antenna__n_parallel()
end

Base.@kwdef mutable struct pulse_schedule__lh
    mode :: pulse_schedule__lh__mode = pulse_schedule__lh__mode()
    antenna :: StructArray{pulse_schedule__lh__antenna} = StructArray(pulse_schedule__lh__antenna() for k in 1:1)
end

Base.@kwdef mutable struct pulse_schedule__code
    library :: StructArray{pulse_schedule__code__library} = StructArray(pulse_schedule__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__position_control__strike_point__r__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__position_control__strike_point__r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__strike_point__r__reference = pulse_schedule__position_control__strike_point__r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__strike_point
    r :: pulse_schedule__position_control__strike_point__r = pulse_schedule__position_control__strike_point__r()
    z :: pulse_schedule__position_control__strike_point__z = pulse_schedule__position_control__strike_point__z()
end

Base.@kwdef mutable struct pulse_schedule__position_control__active_limiter_point__r
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__position_control__active_limiter_point__r__reference = pulse_schedule__position_control__active_limiter_point__r__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__position_control__active_limiter_point
    r :: pulse_schedule__position_control__active_limiter_point__r = pulse_schedule__position_control__active_limiter_point__r()
    z :: pulse_schedule__position_control__active_limiter_point__z = pulse_schedule__position_control__active_limiter_point__z()
end

Base.@kwdef mutable struct pulse_schedule__position_control__gap
    name :: String = ""
    r :: Float64 = 0.0
    value :: pulse_schedule__position_control__gap__value = pulse_schedule__position_control__gap__value()
    z :: Float64 = 0.0
    angle :: Float64 = 0.0
    identifier :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__position_control
    mode :: pulse_schedule__position_control__mode = pulse_schedule__position_control__mode()
    magnetic_axis :: pulse_schedule__position_control__magnetic_axis = pulse_schedule__position_control__magnetic_axis()
    elongation_lower :: pulse_schedule__position_control__elongation_lower = pulse_schedule__position_control__elongation_lower()
    strike_point :: StructArray{pulse_schedule__position_control__strike_point} = StructArray(pulse_schedule__position_control__strike_point() for k in 1:1)
    x_point :: StructArray{pulse_schedule__position_control__x_point} = StructArray(pulse_schedule__position_control__x_point() for k in 1:1)
    gap :: StructArray{pulse_schedule__position_control__gap} = StructArray(pulse_schedule__position_control__gap() for k in 1:1)
    triangularity :: pulse_schedule__position_control__triangularity = pulse_schedule__position_control__triangularity()
    elongation_upper :: pulse_schedule__position_control__elongation_upper = pulse_schedule__position_control__elongation_upper()
    triangularity_upper :: pulse_schedule__position_control__triangularity_upper = pulse_schedule__position_control__triangularity_upper()
    triangularity_lower :: pulse_schedule__position_control__triangularity_lower = pulse_schedule__position_control__triangularity_lower()
    boundary_outline :: StructArray{pulse_schedule__position_control__boundary_outline} = StructArray(pulse_schedule__position_control__boundary_outline() for k in 1:1)
    minor_radius :: pulse_schedule__position_control__minor_radius = pulse_schedule__position_control__minor_radius()
    geometric_axis :: pulse_schedule__position_control__geometric_axis = pulse_schedule__position_control__geometric_axis()
    elongation :: pulse_schedule__position_control__elongation = pulse_schedule__position_control__elongation()
    active_limiter_point :: pulse_schedule__position_control__active_limiter_point = pulse_schedule__position_control__active_limiter_point()
end

Base.@kwdef mutable struct pulse_schedule__flux_control__loop_voltage__reference
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pulse_schedule__flux_control__loop_voltage
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__flux_control__loop_voltage__reference = pulse_schedule__flux_control__loop_voltage__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__flux_control__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher__power_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__ec__launcher
    steering_angle_pol :: pulse_schedule__ec__launcher__steering_angle_pol = pulse_schedule__ec__launcher__steering_angle_pol()
    name :: String = ""
    deposition_rho_tor_norm :: pulse_schedule__ec__launcher__deposition_rho_tor_norm = pulse_schedule__ec__launcher__deposition_rho_tor_norm()
    power_type :: pulse_schedule__ec__launcher__power_type = pulse_schedule__ec__launcher__power_type()
    frequency :: pulse_schedule__ec__launcher__frequency = pulse_schedule__ec__launcher__frequency()
    power :: pulse_schedule__ec__launcher__power = pulse_schedule__ec__launcher__power()
    identifier :: String = ""
    steering_angle_tor :: pulse_schedule__ec__launcher__steering_angle_tor = pulse_schedule__ec__launcher__steering_angle_tor()
end

Base.@kwdef mutable struct pulse_schedule__ec
    launcher :: StructArray{pulse_schedule__ec__launcher} = StructArray(pulse_schedule__ec__launcher() for k in 1:1)
    mode :: pulse_schedule__ec__mode = pulse_schedule__ec__mode()
end

Base.@kwdef mutable struct pulse_schedule__flux_control__beta_normal
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__flux_control__beta_normal__reference = pulse_schedule__flux_control__beta_normal__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__flux_control
    beta_normal :: pulse_schedule__flux_control__beta_normal = pulse_schedule__flux_control__beta_normal()
    mode :: pulse_schedule__flux_control__mode = pulse_schedule__flux_control__mode()
    loop_voltage :: pulse_schedule__flux_control__loop_voltage = pulse_schedule__flux_control__loop_voltage()
    li_3 :: pulse_schedule__flux_control__li_3 = pulse_schedule__flux_control__li_3()
    i_plasma :: pulse_schedule__flux_control__i_plasma = pulse_schedule__flux_control__i_plasma()
end

Base.@kwdef mutable struct pulse_schedule__density_control__n_e_line
    reference_type :: Int32 = 0
    reference_name :: String = ""
    reference :: pulse_schedule__density_control__n_e_line__reference = pulse_schedule__density_control__n_e_line__reference()
    envelope_type :: Int32 = 0
end

Base.@kwdef mutable struct pulse_schedule__density_control
    n_e_line :: pulse_schedule__density_control__n_e_line = pulse_schedule__density_control__n_e_line()
    n_h_over_n_d :: pulse_schedule__density_control__n_h_over_n_d = pulse_schedule__density_control__n_h_over_n_d()
    mode :: pulse_schedule__density_control__mode = pulse_schedule__density_control__mode()
    n_t_over_n_d :: pulse_schedule__density_control__n_t_over_n_d = pulse_schedule__density_control__n_t_over_n_d()
    valve :: StructArray{pulse_schedule__density_control__valve} = StructArray(pulse_schedule__density_control__valve() for k in 1:1)
    zeff :: pulse_schedule__density_control__zeff = pulse_schedule__density_control__zeff()
end

Base.@kwdef mutable struct pulse_schedule__ic__antenna
    name :: String = ""
    power_type :: pulse_schedule__ic__antenna__power_type = pulse_schedule__ic__antenna__power_type()
    phase :: pulse_schedule__ic__antenna__phase = pulse_schedule__ic__antenna__phase()
    frequency :: pulse_schedule__ic__antenna__frequency = pulse_schedule__ic__antenna__frequency()
    power :: pulse_schedule__ic__antenna__power = pulse_schedule__ic__antenna__power()
    identifier :: String = ""
end

Base.@kwdef mutable struct pulse_schedule__ic
    mode :: pulse_schedule__ic__mode = pulse_schedule__ic__mode()
    antenna :: StructArray{pulse_schedule__ic__antenna} = StructArray(pulse_schedule__ic__antenna() for k in 1:1)
end

Base.@kwdef mutable struct pulse_schedule
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: pulse_schedule__code = pulse_schedule__code()
    event :: StructArray{pulse_schedule__event} = StructArray(pulse_schedule__event() for k in 1:1)
    position_control :: pulse_schedule__position_control = pulse_schedule__position_control()
    tf :: pulse_schedule__tf = pulse_schedule__tf()
    ids_properties :: pulse_schedule__ids_properties = pulse_schedule__ids_properties()
    density_control :: pulse_schedule__density_control = pulse_schedule__density_control()
    flux_control :: pulse_schedule__flux_control = pulse_schedule__flux_control()
    lh :: pulse_schedule__lh = pulse_schedule__lh()
    ic :: pulse_schedule__ic = pulse_schedule__ic()
    ec :: pulse_schedule__ec = pulse_schedule__ec()
    nbi :: pulse_schedule__nbi = pulse_schedule__nbi()
end

Base.@kwdef mutable struct polarimeter__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct polarimeter__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct polarimeter__code
    library :: StructArray{polarimeter__code__library} = StructArray(polarimeter__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct polarimeter__channel__ellipticity
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct polarimeter__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct polarimeter__channel__faraday_angle
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct polarimeter__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct polarimeter__ids_properties
    provider :: String = ""
    version_put :: polarimeter__ids_properties__version_put = polarimeter__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct polarimeter__channel__line_of_sight__third_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct polarimeter__channel__line_of_sight
    first_point :: polarimeter__channel__line_of_sight__first_point = polarimeter__channel__line_of_sight__first_point()
    second_point :: polarimeter__channel__line_of_sight__second_point = polarimeter__channel__line_of_sight__second_point()
    third_point :: polarimeter__channel__line_of_sight__third_point = polarimeter__channel__line_of_sight__third_point()
end

Base.@kwdef mutable struct polarimeter__channel
    name :: String = ""
    polarisation_initial :: Float64 = 0.0
    wavelength :: Float64 = 0.0
    ellipticity :: polarimeter__channel__ellipticity = polarimeter__channel__ellipticity()
    line_of_sight :: polarimeter__channel__line_of_sight = polarimeter__channel__line_of_sight()
    faraday_angle :: polarimeter__channel__faraday_angle = polarimeter__channel__faraday_angle()
    identifier :: String = ""
    ellipticity_initial :: Float64 = 0.0
end

Base.@kwdef mutable struct polarimeter
    latency :: Float64 = 0.0
    channel :: StructArray{polarimeter__channel} = StructArray(polarimeter__channel() for k in 1:1)
    ids_properties :: polarimeter__ids_properties = polarimeter__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: polarimeter__code = polarimeter__code()
end

Base.@kwdef mutable struct pf_passive__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct pf_passive__loop__element__geometry__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_passive__ids_properties
    provider :: String = ""
    version_put :: pf_passive__ids_properties__version_put = pf_passive__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct pf_passive__loop__element__geometry__rectangle
    height :: Float64 = 0.0
    r :: Float64 = 0.0
    width :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_passive__loop__element__geometry__arcs_of_circle
    curvature_radii :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_passive__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct pf_passive__code
    library :: StructArray{pf_passive__code__library} = StructArray(pf_passive__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct pf_passive__loop__element__geometry__oblique
    alpha :: Float64 = 0.0
    length_alpha :: Float64 = 0.0
    r :: Float64 = 0.0
    length_beta :: Float64 = 0.0
    z :: Float64 = 0.0
    beta :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_passive__loop__element__geometry
    arcs_of_circle :: pf_passive__loop__element__geometry__arcs_of_circle = pf_passive__loop__element__geometry__arcs_of_circle()
    rectangle :: pf_passive__loop__element__geometry__rectangle = pf_passive__loop__element__geometry__rectangle()
    outline :: pf_passive__loop__element__geometry__outline = pf_passive__loop__element__geometry__outline()
    geometry_type :: Int32 = 0
    oblique :: pf_passive__loop__element__geometry__oblique = pf_passive__loop__element__geometry__oblique()
end

Base.@kwdef mutable struct pf_passive__loop__element
    name :: String = ""
    turns_with_sign :: Float64 = 0.0
    area :: Float64 = 0.0
    geometry :: pf_passive__loop__element__geometry = pf_passive__loop__element__geometry()
    identifier :: String = ""
end

Base.@kwdef mutable struct pf_passive__loop
    current :: Array{Float64, 1} = zeros(Float64,(0))
    name :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    resistivity :: Float64 = 0.0
    resistance :: Float64 = 0.0
    element :: StructArray{pf_passive__loop__element} = StructArray(pf_passive__loop__element() for k in 1:1)
end

Base.@kwdef mutable struct pf_passive
    loop :: StructArray{pf_passive__loop} = StructArray(pf_passive__loop() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: pf_passive__ids_properties = pf_passive__ids_properties()
    code :: pf_passive__code = pf_passive__code()
end

Base.@kwdef mutable struct pf_active__coil__element__geometry__rectangle
    height :: Float64 = 0.0
    r :: Float64 = 0.0
    width :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_active__coil__b_field_max_timed
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__global_quantities
    psi_coils_list :: Array{Int32, 1} = zeros(Int32,(0))
    psi_coils_average :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__circuit__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__supply__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__supply__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__vertical_force__force
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__coil__element__geometry__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct pf_active__circuit__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct pf_active__ids_properties
    provider :: String = ""
    version_put :: pf_active__ids_properties__version_put = pf_active__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct pf_active__coil__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__vertical_force
    name :: String = ""
    combination :: Array{Float64, 1} = zeros(Float64,(0))
    limit_min :: Float64 = 0.0
    force :: pf_active__vertical_force__force = pf_active__vertical_force__force()
    limit_max :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_active__supply
    voltage_limit_min :: Float64 = 0.0
    energy_limit_max :: Float64 = 0.0
    filter_denominator :: Array{Float64, 1} = zeros(Float64,(0))
    current :: pf_active__supply__current = pf_active__supply__current()
    voltage_limit_max :: Float64 = 0.0
    name :: String = ""
    filter_numerator :: Array{Float64, 1} = zeros(Float64,(0))
    nonlinear_model :: String = ""
    current_limit_min :: Float64 = 0.0
    delay :: Float64 = 0.0
    resistance :: Float64 = 0.0
    voltage :: pf_active__supply__voltage = pf_active__supply__voltage()
    current_limiter_gain :: Float64 = 0.0
    identifier :: String = ""
    type :: Int32 = 0
    current_limit_max :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_active__coil__element__geometry__oblique
    alpha :: Float64 = 0.0
    length_alpha :: Float64 = 0.0
    r :: Float64 = 0.0
    length_beta :: Float64 = 0.0
    z :: Float64 = 0.0
    beta :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_active__radial_force__force
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__radial_force
    limit_min :: Float64 = 0.0
    combination :: Array{Float64, 1} = zeros(Float64,(0))
    name :: String = ""
    force :: pf_active__radial_force__force = pf_active__radial_force__force()
    limit_max :: Float64 = 0.0
end

Base.@kwdef mutable struct pf_active__coil__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__coil__element__geometry__arcs_of_circle
    curvature_radii :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pf_active__coil__element__geometry
    arcs_of_circle :: pf_active__coil__element__geometry__arcs_of_circle = pf_active__coil__element__geometry__arcs_of_circle()
    rectangle :: pf_active__coil__element__geometry__rectangle = pf_active__coil__element__geometry__rectangle()
    outline :: pf_active__coil__element__geometry__outline = pf_active__coil__element__geometry__outline()
    geometry_type :: Int32 = 0
    oblique :: pf_active__coil__element__geometry__oblique = pf_active__coil__element__geometry__oblique()
end

Base.@kwdef mutable struct pf_active__coil__element
    name :: String = ""
    turns_with_sign :: Float64 = 0.0
    area :: Float64 = 0.0
    geometry :: pf_active__coil__element__geometry = pf_active__coil__element__geometry()
    identifier :: String = ""
end

Base.@kwdef mutable struct pf_active__coil
    b_field_max :: Array{Float64, 1} = zeros(Float64,(0))
    energy_limit_max :: Float64 = 0.0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    current :: pf_active__coil__current = pf_active__coil__current()
    name :: String = ""
    b_field_max_timed :: pf_active__coil__b_field_max_timed = pf_active__coil__b_field_max_timed()
    resistance :: Float64 = 0.0
    voltage :: pf_active__coil__voltage = pf_active__coil__voltage()
    identifier :: String = ""
    current_limit_max :: Array{Float64, 2} = zeros(Float64,(0,0))
    element :: StructArray{pf_active__coil__element} = StructArray(pf_active__coil__element() for k in 1:1)
end

Base.@kwdef mutable struct pf_active__circuit
    current :: pf_active__circuit__current = pf_active__circuit__current()
    name :: String = ""
    connections :: Array{Int32, 2} = zeros(Int32,(0, 0))
    identifier :: String = ""
    type :: String = ""
    voltage :: pf_active__circuit__voltage = pf_active__circuit__voltage()
end

Base.@kwdef mutable struct pf_active__code
    library :: StructArray{pf_active__code__library} = StructArray(pf_active__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct pf_active
    vertical_force :: StructArray{pf_active__vertical_force} = StructArray(pf_active__vertical_force() for k in 1:1)
    latency :: Float64 = 0.0
    ids_properties :: pf_active__ids_properties = pf_active__ids_properties()
    coil :: StructArray{pf_active__coil} = StructArray(pf_active__coil() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: pf_active__code = pf_active__code()
    global_quantities :: pf_active__global_quantities = pf_active__global_quantities()
    radial_force :: StructArray{pf_active__radial_force} = StructArray(pf_active__radial_force() for k in 1:1)
    circuit :: StructArray{pf_active__circuit} = StructArray(pf_active__circuit() for k in 1:1)
    supply :: StructArray{pf_active__supply} = StructArray(pf_active__supply() for k in 1:1)
end

Base.@kwdef mutable struct pellets__time_slice__pellet__shape__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct pellets__time_slice__pellet__shape
    size :: Array{Float64, 1} = zeros(Float64,(0))
    type :: pellets__time_slice__pellet__shape__type = pellets__time_slice__pellet__shape__type()
end

Base.@kwdef mutable struct pellets__time_slice__pellet__path_geometry__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct pellets__time_slice__pellet__path_profiles__position
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pellets__time_slice__pellet__propellant_gas__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct pellets__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct pellets__ids_properties
    provider :: String = ""
    version_put :: pellets__ids_properties__version_put = pellets__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct pellets__time_slice__pellet__species
    label :: String = ""
    z_n :: Float64 = 0.0
    density :: Float64 = 0.0
    sublimation_energy :: Float64 = 0.0
    a :: Float64 = 0.0
    fraction :: Float64 = 0.0
end

Base.@kwdef mutable struct pellets__time_slice__pellet__path_profiles
    ablation_rate :: Array{Float64, 1} = zeros(Float64,(0))
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    distance :: Array{Float64, 1} = zeros(Float64,(0))
    ablated_particles :: Array{Float64, 1} = zeros(Float64,(0))
    t_e :: Array{Float64, 1} = zeros(Float64,(0))
    n_e :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm_drift :: Array{Float64, 1} = zeros(Float64,(0))
    position :: pellets__time_slice__pellet__path_profiles__position = pellets__time_slice__pellet__path_profiles__position()
    velocity :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct pellets__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct pellets__time_slice__pellet__path_geometry__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct pellets__time_slice__pellet__path_geometry
    first_point :: pellets__time_slice__pellet__path_geometry__first_point = pellets__time_slice__pellet__path_geometry__first_point()
    second_point :: pellets__time_slice__pellet__path_geometry__second_point = pellets__time_slice__pellet__path_geometry__second_point()
end

Base.@kwdef mutable struct pellets__time_slice__pellet__propellant_gas
    label :: String = ""
    molecules_n :: Float64 = 0.0
    element :: StructArray{pellets__time_slice__pellet__propellant_gas__element} = StructArray(pellets__time_slice__pellet__propellant_gas__element() for k in 1:1)
end

Base.@kwdef mutable struct pellets__time_slice__pellet
    path_geometry :: pellets__time_slice__pellet__path_geometry = pellets__time_slice__pellet__path_geometry()
    path_profiles :: pellets__time_slice__pellet__path_profiles = pellets__time_slice__pellet__path_profiles()
    shape :: pellets__time_slice__pellet__shape = pellets__time_slice__pellet__shape()
    propellant_gas :: pellets__time_slice__pellet__propellant_gas = pellets__time_slice__pellet__propellant_gas()
    velocity_initial :: Float64 = 0.0
    species :: StructArray{pellets__time_slice__pellet__species} = StructArray(pellets__time_slice__pellet__species() for k in 1:1)
end

Base.@kwdef mutable struct pellets__time_slice
    pellet :: StructArray{pellets__time_slice__pellet} = StructArray(pellets__time_slice__pellet() for k in 1:1)
    time :: Float64 = 0.0
end

Base.@kwdef mutable struct pellets__code
    library :: StructArray{pellets__code__library} = StructArray(pellets__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct pellets
    time_slice :: StructArray{pellets__time_slice} = StructArray(pellets__time_slice() for k in 1:1)
    latency :: Float64 = 0.0
    ids_properties :: pellets__ids_properties = pellets__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: pellets__code = pellets__code()
end

Base.@kwdef mutable struct numerics__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct numerics__ids_properties
    provider :: String = ""
    version_put :: numerics__ids_properties__version_put = numerics__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct numerics
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: numerics__ids_properties = numerics__ids_properties()
    time_step :: Array{Float64, 1} = zeros(Float64,(0))
    time_end :: Array{Float64, 1} = zeros(Float64,(0))
    time_start :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ntms__time_slice__mode__deltaw
    name :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct ntms__time_slice__mode__detailed_evolution__deltaw
    name :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ntms__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct ntms__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct ntms__ids_properties
    provider :: String = ""
    version_put :: ntms__ids_properties__version_put = ntms__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct ntms__time_slice__mode__torque
    name :: String = ""
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct ntms__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct ntms__time_slice__mode__detailed_evolution__torque
    name :: String = ""
    value :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ntms__time_slice__mode__detailed_evolution
    time_detailed :: Array{Float64, 1} = zeros(Float64,(0))
    calculation_method :: String = ""
    m_pol :: Int32 = 0
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    n_tor :: Int32 = 0
    dfrequency_dt :: Array{Float64, 1} = zeros(Float64,(0))
    torque :: StructArray{ntms__time_slice__mode__detailed_evolution__torque} = StructArray(ntms__time_slice__mode__detailed_evolution__torque() for k in 1:1)
    deltaw :: StructArray{ntms__time_slice__mode__detailed_evolution__deltaw} = StructArray(ntms__time_slice__mode__detailed_evolution__deltaw() for k in 1:1)
    dphase_dt :: Array{Float64, 1} = zeros(Float64,(0))
    width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    delta_diff :: Array{Float64, 2} = zeros(Float64,(0,0))
    dwidth_dt :: Array{Float64, 1} = zeros(Float64,(0))
    frequency :: Array{Float64, 1} = zeros(Float64,(0))
    phase :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ntms__code
    library :: StructArray{ntms__code__library} = StructArray(ntms__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct ntms__time_slice__mode__onset
    n_tor :: Int32 = 0
    time_onset :: Float64 = 0.0
    time_offset :: Float64 = 0.0
    m_pol :: Int32 = 0
    cause :: String = ""
    phase :: Float64 = 0.0
    width :: Float64 = 0.0
end

Base.@kwdef mutable struct ntms__time_slice__mode
    calculation_method :: String = ""
    m_pol :: Int32 = 0
    rho_tor :: Float64 = 0.0
    n_tor :: Int32 = 0
    dfrequency_dt :: Float64 = 0.0
    torque :: StructArray{ntms__time_slice__mode__torque} = StructArray(ntms__time_slice__mode__torque() for k in 1:1)
    detailed_evolution :: ntms__time_slice__mode__detailed_evolution = ntms__time_slice__mode__detailed_evolution()
    deltaw :: StructArray{ntms__time_slice__mode__deltaw} = StructArray(ntms__time_slice__mode__deltaw() for k in 1:1)
    dphase_dt :: Float64 = 0.0
    width :: Float64 = 0.0
    rho_tor_norm :: Float64 = 0.0
    delta_diff :: Array{Float64, 1} = zeros(Float64,(0))
    onset :: ntms__time_slice__mode__onset = ntms__time_slice__mode__onset()
    dwidth_dt :: Float64 = 0.0
    frequency :: Float64 = 0.0
    phase :: Float64 = 0.0
end

Base.@kwdef mutable struct ntms__time_slice
    time :: Float64 = 0.0
    mode :: StructArray{ntms__time_slice__mode} = StructArray(ntms__time_slice__mode() for k in 1:1)
end

Base.@kwdef mutable struct ntms
    time_slice :: StructArray{ntms__time_slice} = StructArray(ntms__time_slice() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: ntms__ids_properties = ntms__ids_properties()
    vacuum_toroidal_field :: ntms__vacuum_toroidal_field = ntms__vacuum_toroidal_field()
    code :: ntms__code = ntms__code()
end

Base.@kwdef mutable struct neutron_diagnostic__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__test_generator__frequency
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__field_of_view__direction_to_detector
    x :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    z :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    y :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__temperature_sensor__frequency
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__amplitude_peak
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__supply_low_voltage__voltage_set
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__test_generator__amplitude
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__supply_low_voltage__voltage_out
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__test_generator__shape
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__temperature_sensor__shape
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__temperature_sensor__amplitude
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__b_field_sensor__frequency
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__characteristics__reaction__mode
    count_limit_min :: Float64 = 0.0
    name :: String = ""
    count_limit_max :: Float64 = 0.0
    index :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__field_of_view__emission_grid
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__unit_source__radiation__reaction
    d2flux_drdz :: Array{Float64, 2} = zeros(Float64,(0,0))
    sensitivity :: Array{Float64, 2} = zeros(Float64,(0,0))
    reaction_rate :: Array{Float64, 2} = zeros(Float64,(0,0))
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct neutron_diagnostic__unit_source__radiation
    reaction :: StructArray{neutron_diagnostic__unit_source__radiation__reaction} = StructArray(neutron_diagnostic__unit_source__radiation__reaction() for k in 1:1)
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__b_field_sensor__amplitude
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__unit_source__position
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__unit_source
    position :: neutron_diagnostic__unit_source__position = neutron_diagnostic__unit_source__position()
    radiation :: StructArray{neutron_diagnostic__unit_source__radiation} = StructArray(neutron_diagnostic__unit_source__radiation() for k in 1:1)
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__field_of_view
    direction_to_detector :: neutron_diagnostic__detectors__field_of_view__direction_to_detector = neutron_diagnostic__detectors__field_of_view__direction_to_detector()
    emission_grid :: neutron_diagnostic__detectors__field_of_view__emission_grid = neutron_diagnostic__detectors__field_of_view__emission_grid()
    solid_angle :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__radiation
    converter_temperature :: Array{Float64, 1} = zeros(Float64,(0))
    converter_volume :: Float64 = 0.0
    converter_name :: String = ""
    converter_nuclear_density :: Float64 = 0.0
    index :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__characteristics__reaction
    error :: Float64 = 0.0
    probability_overlap :: Float64 = 0.0
    mode :: StructArray{neutron_diagnostic__characteristics__reaction__mode} = StructArray(neutron_diagnostic__characteristics__reaction__mode() for k in 1:1)
    index :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__characteristics
    reaction :: StructArray{neutron_diagnostic__characteristics__reaction} = StructArray(neutron_diagnostic__characteristics__reaction() for k in 1:1)
    pulse_length :: Float64 = 0.0
    dead_time :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__test_generator
    rise_time :: Float64 = 0.0
    shape :: neutron_diagnostic__detectors__test_generator__shape = neutron_diagnostic__detectors__test_generator__shape()
    fall_time :: Float64 = 0.0
    frequency :: neutron_diagnostic__detectors__test_generator__frequency = neutron_diagnostic__detectors__test_generator__frequency()
    amplitude :: neutron_diagnostic__detectors__test_generator__amplitude = neutron_diagnostic__detectors__test_generator__amplitude()
    power_switch :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__spectrum
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 2} = zeros(Int32,(0, 0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__synthetic_signals
    total_neutron_flux :: Array{Float64, 1} = zeros(Float64,(0))
    fusion_power :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__amplitude_raw
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__supply_high_voltage__voltage_out
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__adc
    discriminator_level_upper :: Int32 = 0
    input_range :: Float64 = 0.0
    discriminator_level_lower :: Int32 = 0
    sampling_rate :: Int32 = 0
    bias :: Float64 = 0.0
    power_switch :: Int32 = 0
    impedance :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__aperture
    radius :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: neutron_diagnostic__detectors__aperture__centre = neutron_diagnostic__detectors__aperture__centre()
    outline :: neutron_diagnostic__detectors__aperture__outline = neutron_diagnostic__detectors__aperture__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: neutron_diagnostic__detectors__aperture__x3_unit_vector = neutron_diagnostic__detectors__aperture__x3_unit_vector()
    x2_unit_vector :: neutron_diagnostic__detectors__aperture__x2_unit_vector = neutron_diagnostic__detectors__aperture__x2_unit_vector()
    surface :: Float64 = 0.0
    geometry_type :: Int32 = 0
    x1_unit_vector :: neutron_diagnostic__detectors__aperture__x1_unit_vector = neutron_diagnostic__detectors__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__energy_band
    detection_efficiency :: Array{Float64, 1} = zeros(Float64,(0))
    lower_bound :: Float64 = 0.0
    energies :: Array{Float64, 1} = zeros(Float64,(0))
    upper_bound :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__detector
    surface :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: neutron_diagnostic__detectors__detector__centre = neutron_diagnostic__detectors__detector__centre()
    outline :: neutron_diagnostic__detectors__detector__outline = neutron_diagnostic__detectors__detector__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: neutron_diagnostic__detectors__detector__x3_unit_vector = neutron_diagnostic__detectors__detector__x3_unit_vector()
    geometry_type :: Int32 = 0
    radius :: Float64 = 0.0
    x2_unit_vector :: neutron_diagnostic__detectors__detector__x2_unit_vector = neutron_diagnostic__detectors__detector__x2_unit_vector()
    x1_unit_vector :: neutron_diagnostic__detectors__detector__x1_unit_vector = neutron_diagnostic__detectors__detector__x1_unit_vector()
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__b_field_sensor__shape
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__b_field_sensor
    rise_time :: Float64 = 0.0
    shape :: neutron_diagnostic__detectors__b_field_sensor__shape = neutron_diagnostic__detectors__b_field_sensor__shape()
    fall_time :: Float64 = 0.0
    frequency :: neutron_diagnostic__detectors__b_field_sensor__frequency = neutron_diagnostic__detectors__b_field_sensor__frequency()
    amplitude :: neutron_diagnostic__detectors__b_field_sensor__amplitude = neutron_diagnostic__detectors__b_field_sensor__amplitude()
    power_switch :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__supply_high_voltage__voltage_set
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__supply_high_voltage
    voltage_set :: neutron_diagnostic__detectors__supply_high_voltage__voltage_set = neutron_diagnostic__detectors__supply_high_voltage__voltage_set()
    power_switch :: Int32 = 0
    voltage_out :: neutron_diagnostic__detectors__supply_high_voltage__voltage_out = neutron_diagnostic__detectors__supply_high_voltage__voltage_out()
end

Base.@kwdef mutable struct neutron_diagnostic__ids_properties
    provider :: String = ""
    version_put :: neutron_diagnostic__ids_properties__version_put = neutron_diagnostic__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__mode__counting
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__mode
    name :: String = ""
    counting :: neutron_diagnostic__detectors__mode__counting = neutron_diagnostic__detectors__mode__counting()
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__supply_low_voltage
    voltage_set :: neutron_diagnostic__detectors__supply_low_voltage__voltage_set = neutron_diagnostic__detectors__supply_low_voltage__voltage_set()
    power_switch :: Int32 = 0
    voltage_out :: neutron_diagnostic__detectors__supply_low_voltage__voltage_out = neutron_diagnostic__detectors__supply_low_voltage__voltage_out()
end

Base.@kwdef mutable struct neutron_diagnostic__detectors__temperature_sensor
    rise_time :: Float64 = 0.0
    shape :: neutron_diagnostic__detectors__temperature_sensor__shape = neutron_diagnostic__detectors__temperature_sensor__shape()
    fall_time :: Float64 = 0.0
    frequency :: neutron_diagnostic__detectors__temperature_sensor__frequency = neutron_diagnostic__detectors__temperature_sensor__frequency()
    amplitude :: neutron_diagnostic__detectors__temperature_sensor__amplitude = neutron_diagnostic__detectors__temperature_sensor__amplitude()
    power_switch :: Int32 = 0
end

Base.@kwdef mutable struct neutron_diagnostic__detectors
    end_time :: Float64 = 0.0
    field_of_view :: neutron_diagnostic__detectors__field_of_view = neutron_diagnostic__detectors__field_of_view()
    mode :: StructArray{neutron_diagnostic__detectors__mode} = StructArray(neutron_diagnostic__detectors__mode() for k in 1:1)
    energy_band :: StructArray{neutron_diagnostic__detectors__energy_band} = StructArray(neutron_diagnostic__detectors__energy_band() for k in 1:1)
    adc :: neutron_diagnostic__detectors__adc = neutron_diagnostic__detectors__adc()
    detector :: neutron_diagnostic__detectors__detector = neutron_diagnostic__detectors__detector()
    b_field_sensor :: neutron_diagnostic__detectors__b_field_sensor = neutron_diagnostic__detectors__b_field_sensor()
    temperature_sensor :: neutron_diagnostic__detectors__temperature_sensor = neutron_diagnostic__detectors__temperature_sensor()
    supply_high_voltage :: neutron_diagnostic__detectors__supply_high_voltage = neutron_diagnostic__detectors__supply_high_voltage()
    spectrum :: neutron_diagnostic__detectors__spectrum = neutron_diagnostic__detectors__spectrum()
    start_time :: Float64 = 0.0
    name :: String = ""
    spectrum_sampling_time :: Float64 = 0.0
    test_generator :: neutron_diagnostic__detectors__test_generator = neutron_diagnostic__detectors__test_generator()
    position :: neutron_diagnostic__detectors__position = neutron_diagnostic__detectors__position()
    supply_low_voltage :: neutron_diagnostic__detectors__supply_low_voltage = neutron_diagnostic__detectors__supply_low_voltage()
    amplitude_raw :: neutron_diagnostic__detectors__amplitude_raw = neutron_diagnostic__detectors__amplitude_raw()
    spectrum_total :: Array{Int32, 1} = zeros(Int32,(0))
    radiation :: StructArray{neutron_diagnostic__detectors__radiation} = StructArray(neutron_diagnostic__detectors__radiation() for k in 1:1)
    aperture :: StructArray{neutron_diagnostic__detectors__aperture} = StructArray(neutron_diagnostic__detectors__aperture() for k in 1:1)
    amplitude_peak :: neutron_diagnostic__detectors__amplitude_peak = neutron_diagnostic__detectors__amplitude_peak()
end

Base.@kwdef mutable struct neutron_diagnostic__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct neutron_diagnostic__code
    library :: StructArray{neutron_diagnostic__code__library} = StructArray(neutron_diagnostic__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct neutron_diagnostic
    detectors :: StructArray{neutron_diagnostic__detectors} = StructArray(neutron_diagnostic__detectors() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: neutron_diagnostic__ids_properties = neutron_diagnostic__ids_properties()
    latency :: Float64 = 0.0
    synthetic_signals :: neutron_diagnostic__synthetic_signals = neutron_diagnostic__synthetic_signals()
    unit_source :: StructArray{neutron_diagnostic__unit_source} = StructArray(neutron_diagnostic__unit_source() for k in 1:1)
    code :: neutron_diagnostic__code = neutron_diagnostic__code()
    characteristics :: neutron_diagnostic__characteristics = neutron_diagnostic__characteristics()
end

Base.@kwdef mutable struct nbi__unit__source__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__power_launched
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beam_power_fraction
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct nbi__unit__energy
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct nbi__unit__beam_current_fraction
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct nbi__unit__source__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct nbi__unit__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct nbi__unit__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__beamlets__positions
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__divergence_component
    horizontal :: Float64 = 0.0
    particles_fraction :: Float64 = 0.0
    vertical :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__source__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__species
    label :: String = ""
    z_n :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__focus
    width_min_horizontal :: Float64 = 0.0
    width_min_vertical :: Float64 = 0.0
    focal_length_vertical :: Float64 = 0.0
    focal_length_horizontal :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__tilting__delta_position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__source__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__beamlets
    tangency_radii :: Array{Float64, 1} = zeros(Float64,(0))
    power_fractions :: Array{Float64, 1} = zeros(Float64,(0))
    angles :: Array{Float64, 1} = zeros(Float64,(0))
    positions :: nbi__unit__beamlets_group__beamlets__positions = nbi__unit__beamlets_group__beamlets__positions()
end

Base.@kwdef mutable struct nbi__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct nbi__unit__source__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beamlets_group__tilting
    delta_angle :: Float64 = 0.0
    time :: Float64 = 0.0
    delta_position :: nbi__unit__beamlets_group__tilting__delta_position = nbi__unit__beamlets_group__tilting__delta_position()
    delta_tangency_radius :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__beamlets_group
    width_vertical :: Float64 = 0.0
    focus :: nbi__unit__beamlets_group__focus = nbi__unit__beamlets_group__focus()
    divergence_component :: StructArray{nbi__unit__beamlets_group__divergence_component} = StructArray(nbi__unit__beamlets_group__divergence_component() for k in 1:1)
    beamlets :: nbi__unit__beamlets_group__beamlets = nbi__unit__beamlets_group__beamlets()
    direction :: Int32 = 0
    tangency_radius :: Float64 = 0.0
    width_horizontal :: Float64 = 0.0
    tilting :: StructArray{nbi__unit__beamlets_group__tilting} = StructArray(nbi__unit__beamlets_group__tilting() for k in 1:1)
    position :: nbi__unit__beamlets_group__position = nbi__unit__beamlets_group__position()
    angle :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct nbi__ids_properties
    provider :: String = ""
    version_put :: nbi__ids_properties__version_put = nbi__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct nbi__unit__source
    x2_unit_vector :: nbi__unit__source__x2_unit_vector = nbi__unit__source__x2_unit_vector()
    x2_width :: Float64 = 0.0
    centre :: nbi__unit__source__centre = nbi__unit__source__centre()
    radius :: Float64 = 0.0
    x1_width :: Float64 = 0.0
    x3_unit_vector :: nbi__unit__source__x3_unit_vector = nbi__unit__source__x3_unit_vector()
    outline :: nbi__unit__source__outline = nbi__unit__source__outline()
    geometry_type :: Int32 = 0
    surface :: Float64 = 0.0
    x1_unit_vector :: nbi__unit__source__x1_unit_vector = nbi__unit__source__x1_unit_vector()
end

Base.@kwdef mutable struct nbi__unit__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct nbi__unit__aperture
    x1_width :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    radius :: Float64 = 0.0
    centre :: nbi__unit__aperture__centre = nbi__unit__aperture__centre()
    outline :: nbi__unit__aperture__outline = nbi__unit__aperture__outline()
    x2_unit_vector :: nbi__unit__aperture__x2_unit_vector = nbi__unit__aperture__x2_unit_vector()
    x3_unit_vector :: nbi__unit__aperture__x3_unit_vector = nbi__unit__aperture__x3_unit_vector()
    surface :: Float64 = 0.0
    geometry_type :: Int32 = 0
    x1_unit_vector :: nbi__unit__aperture__x1_unit_vector = nbi__unit__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct nbi__unit
    beam_current_fraction :: nbi__unit__beam_current_fraction = nbi__unit__beam_current_fraction()
    name :: String = ""
    source :: nbi__unit__source = nbi__unit__source()
    power_launched :: nbi__unit__power_launched = nbi__unit__power_launched()
    aperture :: StructArray{nbi__unit__aperture} = StructArray(nbi__unit__aperture() for k in 1:1)
    energy :: nbi__unit__energy = nbi__unit__energy()
    identifier :: String = ""
    beam_power_fraction :: nbi__unit__beam_power_fraction = nbi__unit__beam_power_fraction()
    beamlets_group :: StructArray{nbi__unit__beamlets_group} = StructArray(nbi__unit__beamlets_group() for k in 1:1)
    species :: nbi__unit__species = nbi__unit__species()
end

Base.@kwdef mutable struct nbi__code
    library :: StructArray{nbi__code__library} = StructArray(nbi__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct nbi
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: nbi__ids_properties = nbi__ids_properties()
    latency :: Float64 = 0.0
    code :: nbi__code = nbi__code()
    unit :: StructArray{nbi__unit} = StructArray(nbi__unit() for k in 1:1)
end

Base.@kwdef mutable struct mse__channel__polarisation_angle
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct mse__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct mse__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mse__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mse__ids_properties
    provider :: String = ""
    version_put :: mse__ids_properties__version_put = mse__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct mse__channel__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__active_spatial_resolution__width
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct mse__code
    library :: StructArray{mse__code__library} = StructArray(mse__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct mse__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__aperture
    radius :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: mse__channel__aperture__centre = mse__channel__aperture__centre()
    outline :: mse__channel__aperture__outline = mse__channel__aperture__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: mse__channel__aperture__x3_unit_vector = mse__channel__aperture__x3_unit_vector()
    x2_unit_vector :: mse__channel__aperture__x2_unit_vector = mse__channel__aperture__x2_unit_vector()
    surface :: Float64 = 0.0
    geometry_type :: Int32 = 0
    x1_unit_vector :: mse__channel__aperture__x1_unit_vector = mse__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct mse__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__line_of_sight
    first_point :: mse__channel__line_of_sight__first_point = mse__channel__line_of_sight__first_point()
    second_point :: mse__channel__line_of_sight__second_point = mse__channel__line_of_sight__second_point()
end

Base.@kwdef mutable struct mse__channel__detector
    geometry_type :: Int32 = 0
    x2_width :: Float64 = 0.0
    radius :: Float64 = 0.0
    centre :: mse__channel__detector__centre = mse__channel__detector__centre()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: mse__channel__detector__x3_unit_vector = mse__channel__detector__x3_unit_vector()
    x2_unit_vector :: mse__channel__detector__x2_unit_vector = mse__channel__detector__x2_unit_vector()
    outline :: mse__channel__detector__outline = mse__channel__detector__outline()
    surface :: Float64 = 0.0
    x1_unit_vector :: mse__channel__detector__x1_unit_vector = mse__channel__detector__x1_unit_vector()
end

Base.@kwdef mutable struct mse__channel__active_spatial_resolution__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct mse__channel__active_spatial_resolution
    time :: Float64 = 0.0
    geometric_coefficients :: Array{Float64, 1} = zeros(Float64,(0))
    centre :: mse__channel__active_spatial_resolution__centre = mse__channel__active_spatial_resolution__centre()
    width :: mse__channel__active_spatial_resolution__width = mse__channel__active_spatial_resolution__width()
end

Base.@kwdef mutable struct mse__channel
    polarisation_angle :: mse__channel__polarisation_angle = mse__channel__polarisation_angle()
    name :: String = ""
    line_of_sight :: mse__channel__line_of_sight = mse__channel__line_of_sight()
    active_spatial_resolution :: StructArray{mse__channel__active_spatial_resolution} = StructArray(mse__channel__active_spatial_resolution() for k in 1:1)
    aperture :: StructArray{mse__channel__aperture} = StructArray(mse__channel__aperture() for k in 1:1)
    detector :: mse__channel__detector = mse__channel__detector()
end

Base.@kwdef mutable struct mse
    time :: Array{Float64, 1} = zeros(Float64,(0))
    channel :: StructArray{mse__channel} = StructArray(mse__channel() for k in 1:1)
    ids_properties :: mse__ids_properties = mse__ids_properties()
    latency :: Float64 = 0.0
    code :: mse__code = mse__code()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__stress_maxwell
    imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__stress_reynolds
    real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__ntv
    imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate2
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__temperature_perturbed
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__displacement_perpendicular
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate3
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__mass_density_perturbed
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate2
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__psi_potential_perturbed
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate2
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate3
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__phi_potential_perturbed
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system__grid
    volume_element :: Array{Float64, 2} = zeros(Float64,(0,0))
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system__grid
    volume_element :: Array{Float64, 2} = zeros(Float64,(0,0))
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__grid
    volume_element :: Array{Float64, 2} = zeros(Float64,(0,0))
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__perturbation_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system
    jacobian :: Array{Float64, 2} = zeros(Float64,(0,0))
    tensor_contravariant :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
    grid_type :: mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system__grid_type = mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system__grid_type()
    grid :: mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system__grid = mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system__grid()
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
    tensor_covariant :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate2
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate1
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__pressure_perturbed
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate1
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__model_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__displacement_parallel
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate1
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate1
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate3
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed
    coordinate2 :: mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate2 = mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate2()
    coordinate1 :: mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate1 = mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate1()
    coordinate3 :: mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate3 = mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed__coordinate3()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate1
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__ballooning_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__equations
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__alfven_frequency_spectrum
    imaginary :: Array{Float64, 1} = zeros(Float64,(0))
    real :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mhd_linear__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct mhd_linear__code
    library :: StructArray{mhd_linear__code__library} = StructArray(mhd_linear__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate3
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed
    coordinate2 :: mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate2 = mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate2()
    coordinate1 :: mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate1 = mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate1()
    coordinate3 :: mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate3 = mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed__coordinate3()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate3
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__ids_properties
    provider :: String = ""
    version_put :: mhd_linear__ids_properties__version_put = mhd_linear__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate2
    imaginary :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_imaginary :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    real :: Array{Float64, 2} = zeros(Float64,(0,0))
    coefficients_real :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed
    coordinate2 :: mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate2 = mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate2()
    coordinate1 :: mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate1 = mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate1()
    coordinate3 :: mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate3 = mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed__coordinate3()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system
    jacobian :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_type :: mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system__grid_type = mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system__grid_type()
    tensor_contravariant :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid :: mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system__grid = mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system__grid()
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
    tensor_covariant :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed
    coordinate2 :: mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate2 = mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate2()
    coordinate1 :: mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate1 = mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate1()
    coordinate3 :: mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate3 = mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed__coordinate3()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__vacuum
    grid_type :: mhd_linear__time_slice__toroidal_mode__vacuum__grid_type = mhd_linear__time_slice__toroidal_mode__vacuum__grid_type()
    b_field_perturbed :: mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed = mhd_linear__time_slice__toroidal_mode__vacuum__b_field_perturbed()
    coordinate_system :: mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system = mhd_linear__time_slice__toroidal_mode__vacuum__coordinate_system()
    a_field_perturbed :: mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed = mhd_linear__time_slice__toroidal_mode__vacuum__a_field_perturbed()
    grid :: mhd_linear__time_slice__toroidal_mode__vacuum__grid = mhd_linear__time_slice__toroidal_mode__vacuum__grid()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__grid
    volume_element :: Array{Float64, 2} = zeros(Float64,(0,0))
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed
    coordinate2 :: mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate2 = mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate2()
    coordinate1 :: mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate1 = mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate1()
    coordinate3 :: mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate3 = mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed__coordinate3()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode__plasma
    displacement_perpendicular :: mhd_linear__time_slice__toroidal_mode__plasma__displacement_perpendicular = mhd_linear__time_slice__toroidal_mode__plasma__displacement_perpendicular()
    stress_maxwell :: mhd_linear__time_slice__toroidal_mode__plasma__stress_maxwell = mhd_linear__time_slice__toroidal_mode__plasma__stress_maxwell()
    alfven_frequency_spectrum :: StructArray{mhd_linear__time_slice__toroidal_mode__plasma__alfven_frequency_spectrum} = StructArray(mhd_linear__time_slice__toroidal_mode__plasma__alfven_frequency_spectrum() for k in 1:1)
    pressure_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__pressure_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__pressure_perturbed()
    tau_alfven :: Array{Float64, 1} = zeros(Float64,(0))
    phi_potential_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__phi_potential_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__phi_potential_perturbed()
    mass_density_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__mass_density_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__mass_density_perturbed()
    temperature_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__temperature_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__temperature_perturbed()
    ntv :: mhd_linear__time_slice__toroidal_mode__plasma__ntv = mhd_linear__time_slice__toroidal_mode__plasma__ntv()
    psi_potential_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__psi_potential_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__psi_potential_perturbed()
    a_field_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__a_field_perturbed()
    grid :: mhd_linear__time_slice__toroidal_mode__plasma__grid = mhd_linear__time_slice__toroidal_mode__plasma__grid()
    b_field_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__b_field_perturbed()
    coordinate_system :: mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system = mhd_linear__time_slice__toroidal_mode__plasma__coordinate_system()
    velocity_perturbed :: mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed = mhd_linear__time_slice__toroidal_mode__plasma__velocity_perturbed()
    tau_resistive :: Array{Float64, 1} = zeros(Float64,(0))
    grid_type :: mhd_linear__time_slice__toroidal_mode__plasma__grid_type = mhd_linear__time_slice__toroidal_mode__plasma__grid_type()
    stress_reynolds :: mhd_linear__time_slice__toroidal_mode__plasma__stress_reynolds = mhd_linear__time_slice__toroidal_mode__plasma__stress_reynolds()
    displacement_parallel :: mhd_linear__time_slice__toroidal_mode__plasma__displacement_parallel = mhd_linear__time_slice__toroidal_mode__plasma__displacement_parallel()
end

Base.@kwdef mutable struct mhd_linear__time_slice__toroidal_mode
    perturbation_type :: mhd_linear__time_slice__toroidal_mode__perturbation_type = mhd_linear__time_slice__toroidal_mode__perturbation_type()
    plasma :: mhd_linear__time_slice__toroidal_mode__plasma = mhd_linear__time_slice__toroidal_mode__plasma()
    ballooning_type :: mhd_linear__time_slice__toroidal_mode__ballooning_type = mhd_linear__time_slice__toroidal_mode__ballooning_type()
    energy_perturbed :: Float64 = 0.0
    n_tor :: Int32 = 0
    vacuum :: mhd_linear__time_slice__toroidal_mode__vacuum = mhd_linear__time_slice__toroidal_mode__vacuum()
    m_pol_dominant :: Float64 = 0.0
    radial_mode_number :: Float64 = 0.0
    growthrate :: Float64 = 0.0
    phase :: Float64 = 0.0
    frequency :: Float64 = 0.0
    amplitude_multiplier :: Float64 = 0.0
end

Base.@kwdef mutable struct mhd_linear__time_slice
    toroidal_mode :: StructArray{mhd_linear__time_slice__toroidal_mode} = StructArray(mhd_linear__time_slice__toroidal_mode() for k in 1:1)
    time :: Float64 = 0.0
end

Base.@kwdef mutable struct mhd_linear
    ideal_flag :: Int32 = 0
    time_slice :: StructArray{mhd_linear__time_slice} = StructArray(mhd_linear__time_slice() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: mhd_linear__ids_properties = mhd_linear__ids_properties()
    vacuum_toroidal_field :: mhd_linear__vacuum_toroidal_field = mhd_linear__vacuum_toroidal_field()
    equations :: mhd_linear__equations = mhd_linear__equations()
    fluids_n :: Int32 = 0
    model_type :: mhd_linear__model_type = mhd_linear__model_type()
    code :: mhd_linear__code = mhd_linear__code()
end

Base.@kwdef mutable struct mhd__ggd__zeff
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__mass_density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__velocity_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__phi_potential
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__j_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__grid_ggd__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd__grid_ggd__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd__ggd__a_field_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__electrons__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__electrons
    temperature :: StructArray{mhd__ggd__electrons__temperature} = StructArray(mhd__ggd__electrons__temperature() for k in 1:1)
end

Base.@kwdef mutable struct mhd__grid_ggd__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd__grid_ggd__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd__ggd__b_field_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__n_i_total
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__vorticity
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__j_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__velocity_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__grid_ggd__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd__ggd__b_field_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__velocity_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__b_field_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__j_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__grid_ggd__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{mhd__grid_ggd__space__objects_per_dimension__object__boundary} = StructArray(mhd__grid_ggd__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct mhd__grid_ggd__space__objects_per_dimension
    object :: StructArray{mhd__grid_ggd__space__objects_per_dimension__object} = StructArray(mhd__grid_ggd__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct mhd__grid_ggd__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd__ggd__t_i_average
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct mhd__grid_ggd__grid_subset__element
    object :: StructArray{mhd__grid_ggd__grid_subset__element__object} = StructArray(mhd__grid_ggd__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct mhd__ggd__velocity_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__ggd__psi
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__code
    library :: StructArray{mhd__code__library} = StructArray(mhd__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct mhd__ggd__a_field_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__grid_ggd__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct mhd__grid_ggd__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{mhd__grid_ggd__space__objects_per_dimension} = StructArray(mhd__grid_ggd__space__objects_per_dimension() for k in 1:1)
    identifier :: mhd__grid_ggd__space__identifier = mhd__grid_ggd__space__identifier()
    geometry_type :: mhd__grid_ggd__space__geometry_type = mhd__grid_ggd__space__geometry_type()
end

Base.@kwdef mutable struct mhd__ggd__a_field_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct mhd__grid_ggd__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct mhd__grid_ggd__grid_subset
    base :: StructArray{mhd__grid_ggd__grid_subset__base} = StructArray(mhd__grid_ggd__grid_subset__base() for k in 1:1)
    metric :: mhd__grid_ggd__grid_subset__metric = mhd__grid_ggd__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: mhd__grid_ggd__grid_subset__identifier = mhd__grid_ggd__grid_subset__identifier()
    element :: StructArray{mhd__grid_ggd__grid_subset__element} = StructArray(mhd__grid_ggd__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct mhd__grid_ggd
    time :: Float64 = 0.0
    grid_subset :: StructArray{mhd__grid_ggd__grid_subset} = StructArray(mhd__grid_ggd__grid_subset() for k in 1:1)
    space :: StructArray{mhd__grid_ggd__space} = StructArray(mhd__grid_ggd__space() for k in 1:1)
    identifier :: mhd__grid_ggd__identifier = mhd__grid_ggd__identifier()
end

Base.@kwdef mutable struct mhd__ggd
    b_field_z :: StructArray{mhd__ggd__b_field_z} = StructArray(mhd__ggd__b_field_z() for k in 1:1)
    time :: Float64 = 0.0
    a_field_z :: StructArray{mhd__ggd__a_field_z} = StructArray(mhd__ggd__a_field_z() for k in 1:1)
    velocity_parallel :: StructArray{mhd__ggd__velocity_parallel} = StructArray(mhd__ggd__velocity_parallel() for k in 1:1)
    psi :: StructArray{mhd__ggd__psi} = StructArray(mhd__ggd__psi() for k in 1:1)
    velocity_z :: StructArray{mhd__ggd__velocity_z} = StructArray(mhd__ggd__velocity_z() for k in 1:1)
    phi_potential :: StructArray{mhd__ggd__phi_potential} = StructArray(mhd__ggd__phi_potential() for k in 1:1)
    zeff :: StructArray{mhd__ggd__zeff} = StructArray(mhd__ggd__zeff() for k in 1:1)
    electrons :: mhd__ggd__electrons = mhd__ggd__electrons()
    j_r :: StructArray{mhd__ggd__j_r} = StructArray(mhd__ggd__j_r() for k in 1:1)
    a_field_tor :: StructArray{mhd__ggd__a_field_tor} = StructArray(mhd__ggd__a_field_tor() for k in 1:1)
    j_tor :: StructArray{mhd__ggd__j_tor} = StructArray(mhd__ggd__j_tor() for k in 1:1)
    velocity_tor :: StructArray{mhd__ggd__velocity_tor} = StructArray(mhd__ggd__velocity_tor() for k in 1:1)
    b_field_tor :: StructArray{mhd__ggd__b_field_tor} = StructArray(mhd__ggd__b_field_tor() for k in 1:1)
    t_i_average :: StructArray{mhd__ggd__t_i_average} = StructArray(mhd__ggd__t_i_average() for k in 1:1)
    vorticity :: StructArray{mhd__ggd__vorticity} = StructArray(mhd__ggd__vorticity() for k in 1:1)
    a_field_r :: StructArray{mhd__ggd__a_field_r} = StructArray(mhd__ggd__a_field_r() for k in 1:1)
    b_field_r :: StructArray{mhd__ggd__b_field_r} = StructArray(mhd__ggd__b_field_r() for k in 1:1)
    mass_density :: StructArray{mhd__ggd__mass_density} = StructArray(mhd__ggd__mass_density() for k in 1:1)
    j_z :: StructArray{mhd__ggd__j_z} = StructArray(mhd__ggd__j_z() for k in 1:1)
    velocity_r :: StructArray{mhd__ggd__velocity_r} = StructArray(mhd__ggd__velocity_r() for k in 1:1)
    n_i_total :: StructArray{mhd__ggd__n_i_total} = StructArray(mhd__ggd__n_i_total() for k in 1:1)
end

Base.@kwdef mutable struct mhd__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct mhd__ids_properties
    provider :: String = ""
    version_put :: mhd__ids_properties__version_put = mhd__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct mhd
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: mhd__ids_properties = mhd__ids_properties()
    grid_ggd :: StructArray{mhd__grid_ggd} = StructArray(mhd__grid_ggd() for k in 1:1)
    code :: mhd__code = mhd__code()
    ggd :: StructArray{mhd__ggd} = StructArray(mhd__ggd() for k in 1:1)
end

Base.@kwdef mutable struct magnetics__method__ip
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct magnetics__b_field_tor_probe__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__b_field_pol_probe__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__shunt__position__second_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__b_field_tor_probe__non_linear_response
    b_field_linear :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_non_linear :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct magnetics__shunt__position__first_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__shunt__position
    first_point :: magnetics__shunt__position__first_point = magnetics__shunt__position__first_point()
    second_point :: magnetics__shunt__position__second_point = magnetics__shunt__position__second_point()
end

Base.@kwdef mutable struct magnetics__diamagnetic_flux
    method_name :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct magnetics__flux_loop__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__method
    name :: String = ""
    ip :: magnetics__method__ip = magnetics__method__ip()
end

Base.@kwdef mutable struct magnetics__bpol_probe__non_linear_response
    b_field_linear :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_non_linear :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct magnetics__ip
    method_name :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct magnetics__rogowski_coil__measured_quantity
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct magnetics__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct magnetics__b_field_pol_probe__field
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__flux_loop__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__shunt__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__flux_loop__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct magnetics__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct magnetics__bpol_probe__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct magnetics__b_field_tor_probe__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct magnetics__b_field_pol_probe__non_linear_response
    b_field_linear :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_non_linear :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct magnetics__flux_loop__flux
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__flux_loop
    name :: String = ""
    indices_differential :: Array{Int32, 1} = zeros(Int32,(0))
    area :: Float64 = 0.0
    gm9 :: Float64 = 0.0
    flux :: magnetics__flux_loop__flux = magnetics__flux_loop__flux()
    position :: StructArray{magnetics__flux_loop__position} = StructArray(magnetics__flux_loop__position() for k in 1:1)
    identifier :: String = ""
    type :: magnetics__flux_loop__type = magnetics__flux_loop__type()
    voltage :: magnetics__flux_loop__voltage = magnetics__flux_loop__voltage()
end

Base.@kwdef mutable struct magnetics__bpol_probe__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__b_field_tor_probe__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__shunt
    tile_index :: Int32 = 0
    name :: String = ""
    target_index :: Int32 = 0
    position :: magnetics__shunt__position = magnetics__shunt__position()
    resistance :: Float64 = 0.0
    identifier :: String = ""
    divertor_index :: Int32 = 0
    voltage :: magnetics__shunt__voltage = magnetics__shunt__voltage()
end

Base.@kwdef mutable struct magnetics__rogowski_coil__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__b_field_pol_probe__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__b_field_tor_probe__field
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__bpol_probe__field
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__b_field_pol_probe__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct magnetics__b_field_pol_probe
    length :: Float64 = 0.0
    bandwidth_3db :: Array{Float64, 1} = zeros(Float64,(0))
    name :: String = ""
    field :: magnetics__b_field_pol_probe__field = magnetics__b_field_pol_probe__field()
    position :: magnetics__b_field_pol_probe__position = magnetics__b_field_pol_probe__position()
    turns :: Int32 = 0
    indices_differential :: Array{Int32, 1} = zeros(Int32,(0))
    area :: Float64 = 0.0
    toroidal_angle :: Float64 = 0.0
    poloidal_angle :: Float64 = 0.0
    voltage :: magnetics__b_field_pol_probe__voltage = magnetics__b_field_pol_probe__voltage()
    non_linear_response :: magnetics__b_field_pol_probe__non_linear_response = magnetics__b_field_pol_probe__non_linear_response()
    identifier :: String = ""
    type :: magnetics__b_field_pol_probe__type = magnetics__b_field_pol_probe__type()
end

Base.@kwdef mutable struct magnetics__rogowski_coil__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct magnetics__rogowski_coil
    indices_compound :: Array{Int32, 1} = zeros(Int32,(0))
    measured_quantity :: magnetics__rogowski_coil__measured_quantity = magnetics__rogowski_coil__measured_quantity()
    turns_per_metre :: Float64 = 0.0
    name :: String = ""
    area :: Float64 = 0.0
    position :: StructArray{magnetics__rogowski_coil__position} = StructArray(magnetics__rogowski_coil__position() for k in 1:1)
    identifier :: String = ""
    current :: magnetics__rogowski_coil__current = magnetics__rogowski_coil__current()
end

Base.@kwdef mutable struct magnetics__bpol_probe__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct magnetics__b_field_tor_probe
    length :: Float64 = 0.0
    bandwidth_3db :: Array{Float64, 1} = zeros(Float64,(0))
    name :: String = ""
    field :: magnetics__b_field_tor_probe__field = magnetics__b_field_tor_probe__field()
    position :: magnetics__b_field_tor_probe__position = magnetics__b_field_tor_probe__position()
    turns :: Int32 = 0
    indices_differential :: Array{Int32, 1} = zeros(Int32,(0))
    area :: Float64 = 0.0
    toroidal_angle :: Float64 = 0.0
    poloidal_angle :: Float64 = 0.0
    voltage :: magnetics__b_field_tor_probe__voltage = magnetics__b_field_tor_probe__voltage()
    non_linear_response :: magnetics__b_field_tor_probe__non_linear_response = magnetics__b_field_tor_probe__non_linear_response()
    identifier :: String = ""
    type :: magnetics__b_field_tor_probe__type = magnetics__b_field_tor_probe__type()
end

Base.@kwdef mutable struct magnetics__code
    library :: StructArray{magnetics__code__library} = StructArray(magnetics__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct magnetics__ids_properties
    provider :: String = ""
    version_put :: magnetics__ids_properties__version_put = magnetics__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct magnetics__bpol_probe
    length :: Float64 = 0.0
    bandwidth_3db :: Array{Float64, 1} = zeros(Float64,(0))
    name :: String = ""
    field :: magnetics__bpol_probe__field = magnetics__bpol_probe__field()
    position :: magnetics__bpol_probe__position = magnetics__bpol_probe__position()
    turns :: Int32 = 0
    indices_differential :: Array{Int32, 1} = zeros(Int32,(0))
    area :: Float64 = 0.0
    toroidal_angle :: Float64 = 0.0
    poloidal_angle :: Float64 = 0.0
    voltage :: magnetics__bpol_probe__voltage = magnetics__bpol_probe__voltage()
    non_linear_response :: magnetics__bpol_probe__non_linear_response = magnetics__bpol_probe__non_linear_response()
    identifier :: String = ""
    type :: magnetics__bpol_probe__type = magnetics__bpol_probe__type()
end

Base.@kwdef mutable struct magnetics
    b_field_tor_probe :: StructArray{magnetics__b_field_tor_probe} = StructArray(magnetics__b_field_tor_probe() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: magnetics__code = magnetics__code()
    latency :: Float64 = 0.0
    ip :: StructArray{magnetics__ip} = StructArray(magnetics__ip() for k in 1:1)
    method :: StructArray{magnetics__method} = StructArray(magnetics__method() for k in 1:1)
    ids_properties :: magnetics__ids_properties = magnetics__ids_properties()
    bpol_probe :: StructArray{magnetics__bpol_probe} = StructArray(magnetics__bpol_probe() for k in 1:1)
    diamagnetic_flux :: StructArray{magnetics__diamagnetic_flux} = StructArray(magnetics__diamagnetic_flux() for k in 1:1)
    b_field_pol_probe :: StructArray{magnetics__b_field_pol_probe} = StructArray(magnetics__b_field_pol_probe() for k in 1:1)
    rogowski_coil :: StructArray{magnetics__rogowski_coil} = StructArray(magnetics__rogowski_coil() for k in 1:1)
    shunt :: StructArray{magnetics__shunt} = StructArray(magnetics__shunt() for k in 1:1)
    flux_loop :: StructArray{magnetics__flux_loop} = StructArray(magnetics__flux_loop() for k in 1:1)
end

Base.@kwdef mutable struct lh_antennas__antenna__n_parallel_peak
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__power_forward
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct lh_antennas__antenna__pressure_tank
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__power_launched
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__n_e
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct lh_antennas__antenna__power_reflected
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__ids_properties
    provider :: String = ""
    version_put :: lh_antennas__ids_properties__version_put = lh_antennas__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct lh_antennas__antenna__reflection_coefficient
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__phase_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__position__phi
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__power_launched
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__position__z
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__row__position
    time :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct lh_antennas__antenna__position__r
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__position
    phi :: lh_antennas__antenna__position__phi = lh_antennas__antenna__position__phi()
    r :: lh_antennas__antenna__position__r = lh_antennas__antenna__position__r()
    z :: lh_antennas__antenna__position__z = lh_antennas__antenna__position__z()
    definition :: String = ""
end

Base.@kwdef mutable struct lh_antennas__power
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct lh_antennas__antenna__row
    name :: String = ""
    n_tor :: Array{Float64, 1} = zeros(Float64,(0))
    power_density_spectrum_2d :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
    n_pol :: Array{Float64, 1} = zeros(Float64,(0))
    position :: lh_antennas__antenna__row__position = lh_antennas__antenna__row__position()
    power_density_spectrum_1d :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct lh_antennas__code
    library :: StructArray{lh_antennas__code__library} = StructArray(lh_antennas__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct lh_antennas__reference_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct lh_antennas__antenna
    n_parallel_peak :: lh_antennas__antenna__n_parallel_peak = lh_antennas__antenna__n_parallel_peak()
    power_reflected :: lh_antennas__antenna__power_reflected = lh_antennas__antenna__power_reflected()
    phase_average :: lh_antennas__antenna__phase_average = lh_antennas__antenna__phase_average()
    name :: String = ""
    power_launched :: lh_antennas__antenna__power_launched = lh_antennas__antenna__power_launched()
    position :: lh_antennas__antenna__position = lh_antennas__antenna__position()
    power_forward :: lh_antennas__antenna__power_forward = lh_antennas__antenna__power_forward()
    row :: StructArray{lh_antennas__antenna__row} = StructArray(lh_antennas__antenna__row() for k in 1:1)
    pressure_tank :: lh_antennas__antenna__pressure_tank = lh_antennas__antenna__pressure_tank()
    n_e :: lh_antennas__antenna__n_e = lh_antennas__antenna__n_e()
    distance_to_antenna :: Array{Float64, 1} = zeros(Float64,(0))
    model_name :: String = ""
    frequency :: Float64 = 0.0
    reflection_coefficient :: lh_antennas__antenna__reflection_coefficient = lh_antennas__antenna__reflection_coefficient()
    identifier :: String = ""
end

Base.@kwdef mutable struct lh_antennas
    latency :: Float64 = 0.0
    reference_point :: lh_antennas__reference_point = lh_antennas__reference_point()
    power_launched :: lh_antennas__power_launched = lh_antennas__power_launched()
    ids_properties :: lh_antennas__ids_properties = lh_antennas__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    power :: lh_antennas__power = lh_antennas__power()
    code :: lh_antennas__code = lh_antennas__code()
    antenna :: StructArray{lh_antennas__antenna} = StructArray(lh_antennas__antenna() for k in 1:1)
end

Base.@kwdef mutable struct langmuir_probes__embedded__multi_temperature_fits__t_e
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__t_e
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__v_floating_sigma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__t_i
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__t_i
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__j_i_parallel
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__b_field_angle
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__j_i_saturation_skew
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__v_plasma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__v_plasma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__t_i_average
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__v_floating
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__distance_separatrix_midplane
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__heat_flux_parallel
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__t_e
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__j_i_parallel_sigma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__j_i_saturation_sigma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct langmuir_probes__embedded__ion_saturation_current
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__v_floating
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__surface_area_effective
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct langmuir_probes__embedded__multi_temperature_fits__t_i
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__multi_temperature_fits
    time :: Array{Float64, 1} = zeros(Float64,(0))
    t_e :: langmuir_probes__embedded__multi_temperature_fits__t_e = langmuir_probes__embedded__multi_temperature_fits__t_e()
    t_i :: langmuir_probes__embedded__multi_temperature_fits__t_i = langmuir_probes__embedded__multi_temperature_fits__t_i()
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__distance_x_point_z
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__j_i_saturation
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__t_e_average
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__distance_separatrix_midplane
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__ion_saturation_current
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__heat_flux_parallel
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__equilibrium_id__data_entry
    pulse_type :: String = ""
    run :: Int32 = 0
    machine :: String = ""
    pulse :: Int32 = 0
    user :: String = ""
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__j_i_kurtosis
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__j_i_sigma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__b_field_angle
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct langmuir_probes__ids_properties
    provider :: String = ""
    version_put :: langmuir_probes__ids_properties__version_put = langmuir_probes__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct langmuir_probes__embedded__n_e
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__j_i_saturation
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__j_i_saturation_kurtosis
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__j_i_skew
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__n_e
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded__j_i_parallel
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__equilibrium_id
    name :: String = ""
    data_entry :: langmuir_probes__equilibrium_id__data_entry = langmuir_probes__equilibrium_id__data_entry()
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct langmuir_probes__embedded__v_floating_sigma
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector__position
    validity :: Int32 = 0
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__mach_number_parallel
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__embedded
    v_floating_sigma :: langmuir_probes__embedded__v_floating_sigma = langmuir_probes__embedded__v_floating_sigma()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    j_i_saturation_skew :: langmuir_probes__embedded__j_i_saturation_skew = langmuir_probes__embedded__j_i_saturation_skew()
    surface_area :: Float64 = 0.0
    v_plasma :: langmuir_probes__embedded__v_plasma = langmuir_probes__embedded__v_plasma()
    v_floating :: langmuir_probes__embedded__v_floating = langmuir_probes__embedded__v_floating()
    j_i_saturation_kurtosis :: langmuir_probes__embedded__j_i_saturation_kurtosis = langmuir_probes__embedded__j_i_saturation_kurtosis()
    j_i_parallel_sigma :: langmuir_probes__embedded__j_i_parallel_sigma = langmuir_probes__embedded__j_i_parallel_sigma()
    name :: String = ""
    t_e :: langmuir_probes__embedded__t_e = langmuir_probes__embedded__t_e()
    j_i_saturation :: langmuir_probes__embedded__j_i_saturation = langmuir_probes__embedded__j_i_saturation()
    ion_saturation_current :: langmuir_probes__embedded__ion_saturation_current = langmuir_probes__embedded__ion_saturation_current()
    position :: langmuir_probes__embedded__position = langmuir_probes__embedded__position()
    j_i_saturation_sigma :: langmuir_probes__embedded__j_i_saturation_sigma = langmuir_probes__embedded__j_i_saturation_sigma()
    j_i_parallel :: langmuir_probes__embedded__j_i_parallel = langmuir_probes__embedded__j_i_parallel()
    surface_area_effective :: langmuir_probes__embedded__surface_area_effective = langmuir_probes__embedded__surface_area_effective()
    distance_separatrix_midplane :: langmuir_probes__embedded__distance_separatrix_midplane = langmuir_probes__embedded__distance_separatrix_midplane()
    b_field_angle :: langmuir_probes__embedded__b_field_angle = langmuir_probes__embedded__b_field_angle()
    t_i :: langmuir_probes__embedded__t_i = langmuir_probes__embedded__t_i()
    multi_temperature_fits :: StructArray{langmuir_probes__embedded__multi_temperature_fits} = StructArray(langmuir_probes__embedded__multi_temperature_fits() for k in 1:1)
    n_e :: langmuir_probes__embedded__n_e = langmuir_probes__embedded__n_e()
    heat_flux_parallel :: langmuir_probes__embedded__heat_flux_parallel = langmuir_probes__embedded__heat_flux_parallel()
    identifier :: String = ""
end

Base.@kwdef mutable struct langmuir_probes__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct langmuir_probes__code
    library :: StructArray{langmuir_probes__code__library} = StructArray(langmuir_probes__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__collector
    v_floating_sigma :: langmuir_probes__reciprocating__plunge__collector__v_floating_sigma = langmuir_probes__reciprocating__plunge__collector__v_floating_sigma()
    v_floating :: langmuir_probes__reciprocating__plunge__collector__v_floating = langmuir_probes__reciprocating__plunge__collector__v_floating()
    t_e :: langmuir_probes__reciprocating__plunge__collector__t_e = langmuir_probes__reciprocating__plunge__collector__t_e()
    j_i_saturation :: langmuir_probes__reciprocating__plunge__collector__j_i_saturation = langmuir_probes__reciprocating__plunge__collector__j_i_saturation()
    ion_saturation_current :: langmuir_probes__reciprocating__plunge__collector__ion_saturation_current = langmuir_probes__reciprocating__plunge__collector__ion_saturation_current()
    position :: langmuir_probes__reciprocating__plunge__collector__position = langmuir_probes__reciprocating__plunge__collector__position()
    j_i_parallel :: langmuir_probes__reciprocating__plunge__collector__j_i_parallel = langmuir_probes__reciprocating__plunge__collector__j_i_parallel()
    t_i :: langmuir_probes__reciprocating__plunge__collector__t_i = langmuir_probes__reciprocating__plunge__collector__t_i()
    j_i_kurtosis :: langmuir_probes__reciprocating__plunge__collector__j_i_kurtosis = langmuir_probes__reciprocating__plunge__collector__j_i_kurtosis()
    j_i_sigma :: langmuir_probes__reciprocating__plunge__collector__j_i_sigma = langmuir_probes__reciprocating__plunge__collector__j_i_sigma()
    j_i_skew :: langmuir_probes__reciprocating__plunge__collector__j_i_skew = langmuir_probes__reciprocating__plunge__collector__j_i_skew()
    heat_flux_parallel :: langmuir_probes__reciprocating__plunge__collector__heat_flux_parallel = langmuir_probes__reciprocating__plunge__collector__heat_flux_parallel()
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge__position_average
    validity :: Int32 = 0
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct langmuir_probes__reciprocating__plunge
    time :: Float64 = 0.0
    time_within_plunge :: Array{Float64, 1} = zeros(Float64,(0))
    mach_number_parallel :: langmuir_probes__reciprocating__plunge__mach_number_parallel = langmuir_probes__reciprocating__plunge__mach_number_parallel()
    v_plasma :: langmuir_probes__reciprocating__plunge__v_plasma = langmuir_probes__reciprocating__plunge__v_plasma()
    position_average :: langmuir_probes__reciprocating__plunge__position_average = langmuir_probes__reciprocating__plunge__position_average()
    t_i_average :: langmuir_probes__reciprocating__plunge__t_i_average = langmuir_probes__reciprocating__plunge__t_i_average()
    collector :: StructArray{langmuir_probes__reciprocating__plunge__collector} = StructArray(langmuir_probes__reciprocating__plunge__collector() for k in 1:1)
    distance_separatrix_midplane :: langmuir_probes__reciprocating__plunge__distance_separatrix_midplane = langmuir_probes__reciprocating__plunge__distance_separatrix_midplane()
    b_field_angle :: langmuir_probes__reciprocating__plunge__b_field_angle = langmuir_probes__reciprocating__plunge__b_field_angle()
    t_e_average :: langmuir_probes__reciprocating__plunge__t_e_average = langmuir_probes__reciprocating__plunge__t_e_average()
    n_e :: langmuir_probes__reciprocating__plunge__n_e = langmuir_probes__reciprocating__plunge__n_e()
    distance_x_point_z :: langmuir_probes__reciprocating__plunge__distance_x_point_z = langmuir_probes__reciprocating__plunge__distance_x_point_z()
end

Base.@kwdef mutable struct langmuir_probes__reciprocating
    name :: String = ""
    plunge :: StructArray{langmuir_probes__reciprocating__plunge} = StructArray(langmuir_probes__reciprocating__plunge() for k in 1:1)
    surface_area :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: String = ""
end

Base.@kwdef mutable struct langmuir_probes
    equilibrium_id :: langmuir_probes__equilibrium_id = langmuir_probes__equilibrium_id()
    reciprocating :: StructArray{langmuir_probes__reciprocating} = StructArray(langmuir_probes__reciprocating() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: langmuir_probes__ids_properties = langmuir_probes__ids_properties()
    latency :: Float64 = 0.0
    embedded :: StructArray{langmuir_probes__embedded} = StructArray(langmuir_probes__embedded() for k in 1:1)
    midplane :: langmuir_probes__midplane = langmuir_probes__midplane()
    code :: langmuir_probes__code = langmuir_probes__code()
end

Base.@kwdef mutable struct iron_core__segment__geometry__arcs_of_circle
    curvature_radii :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct iron_core__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct iron_core__segment__magnetisation_z
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct iron_core__segment__magnetisation_r
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct iron_core__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct iron_core__segment__geometry__rectangle
    height :: Float64 = 0.0
    r :: Float64 = 0.0
    width :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct iron_core__segment__geometry__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct iron_core__code
    library :: StructArray{iron_core__code__library} = StructArray(iron_core__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct iron_core__ids_properties
    provider :: String = ""
    version_put :: iron_core__ids_properties__version_put = iron_core__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct iron_core__segment__geometry__oblique
    alpha :: Float64 = 0.0
    length_alpha :: Float64 = 0.0
    r :: Float64 = 0.0
    length_beta :: Float64 = 0.0
    z :: Float64 = 0.0
    beta :: Float64 = 0.0
end

Base.@kwdef mutable struct iron_core__segment__geometry
    rectangle :: iron_core__segment__geometry__rectangle = iron_core__segment__geometry__rectangle()
    arcs_of_circle :: iron_core__segment__geometry__arcs_of_circle = iron_core__segment__geometry__arcs_of_circle()
    outline :: iron_core__segment__geometry__outline = iron_core__segment__geometry__outline()
    geometry_type :: Int32 = 0
    oblique :: iron_core__segment__geometry__oblique = iron_core__segment__geometry__oblique()
end

Base.@kwdef mutable struct iron_core__segment
    name :: String = ""
    geometry :: iron_core__segment__geometry = iron_core__segment__geometry()
    magnetisation_r :: iron_core__segment__magnetisation_r = iron_core__segment__magnetisation_r()
    b_field :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: String = ""
    magnetisation_z :: iron_core__segment__magnetisation_z = iron_core__segment__magnetisation_z()
    permeability_relative :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct iron_core
    segment :: StructArray{iron_core__segment} = StructArray(iron_core__segment() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: iron_core__ids_properties = iron_core__ids_properties()
    code :: iron_core__code = iron_core__code()
end

Base.@kwdef mutable struct interferometer__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct interferometer__channel__wavelength__phase_corrected
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct interferometer__n_e_volume_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct interferometer__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct interferometer__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct interferometer__ids_properties
    provider :: String = ""
    version_put :: interferometer__ids_properties__version_put = interferometer__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct interferometer__channel__n_e_line
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct interferometer__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct interferometer__channel__line_of_sight__third_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct interferometer__channel__line_of_sight
    first_point :: interferometer__channel__line_of_sight__first_point = interferometer__channel__line_of_sight__first_point()
    second_point :: interferometer__channel__line_of_sight__second_point = interferometer__channel__line_of_sight__second_point()
    third_point :: interferometer__channel__line_of_sight__third_point = interferometer__channel__line_of_sight__third_point()
end

Base.@kwdef mutable struct interferometer__channel__n_e_line_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct interferometer__channel__wavelength
    phase_corrected :: interferometer__channel__wavelength__phase_corrected = interferometer__channel__wavelength__phase_corrected()
    fringe_jump_correction :: Array{Int32, 1} = zeros(Int32,(0))
    phase_to_n_e_line :: Float64 = 0.0
    fringe_jump_correction_times :: Array{Float64, 1} = zeros(Float64,(0))
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct interferometer__code
    library :: StructArray{interferometer__code__library} = StructArray(interferometer__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct interferometer__electrons_n
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct interferometer__channel__path_length_variation
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct interferometer__channel
    wavelength :: StructArray{interferometer__channel__wavelength} = StructArray(interferometer__channel__wavelength() for k in 1:1)
    n_e_line_average :: interferometer__channel__n_e_line_average = interferometer__channel__n_e_line_average()
    name :: String = ""
    n_e_line :: interferometer__channel__n_e_line = interferometer__channel__n_e_line()
    line_of_sight :: interferometer__channel__line_of_sight = interferometer__channel__line_of_sight()
    path_length_variation :: interferometer__channel__path_length_variation = interferometer__channel__path_length_variation()
    identifier :: String = ""
end

Base.@kwdef mutable struct interferometer
    electrons_n :: interferometer__electrons_n = interferometer__electrons_n()
    latency :: Float64 = 0.0
    channel :: StructArray{interferometer__channel} = StructArray(interferometer__channel() for k in 1:1)
    ids_properties :: interferometer__ids_properties = interferometer__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: interferometer__code = interferometer__code()
    n_e_volume_average :: interferometer__n_e_volume_average = interferometer__n_e_volume_average()
end

Base.@kwdef mutable struct ic_antennas__antenna__power_reflected
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ic_antennas__antenna__power_forward
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ic_antennas__antenna__frequency
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ic_antennas__power_launched
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ic_antennas__antenna__power_launched
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ic_antennas__reference_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct ic_antennas__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct ic_antennas__ids_properties
    provider :: String = ""
    version_put :: ic_antennas__ids_properties__version_put = ic_antennas__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct ic_antennas__antenna__surface_current
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    time :: Float64 = 0.0
    m_pol :: Array{Int32, 1} = zeros(Int32,(0))
    spectrum :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct ic_antennas__antenna
    name :: String = ""
    surface_current :: StructArray{ic_antennas__antenna__surface_current} = StructArray(ic_antennas__antenna__surface_current() for k in 1:1)
    power_launched :: ic_antennas__antenna__power_launched = ic_antennas__antenna__power_launched()
    frequency :: ic_antennas__antenna__frequency = ic_antennas__antenna__frequency()
    power_forward :: ic_antennas__antenna__power_forward = ic_antennas__antenna__power_forward()
    identifier :: String = ""
    power_reflected :: ic_antennas__antenna__power_reflected = ic_antennas__antenna__power_reflected()
end

Base.@kwdef mutable struct ic_antennas__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct ic_antennas__code
    library :: StructArray{ic_antennas__code__library} = StructArray(ic_antennas__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct ic_antennas
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: ic_antennas__ids_properties = ic_antennas__ids_properties()
    reference_point :: ic_antennas__reference_point = ic_antennas__reference_point()
    power_launched :: ic_antennas__power_launched = ic_antennas__power_launched()
    latency :: Float64 = 0.0
    code :: ic_antennas__code = ic_antennas__code()
    antenna :: StructArray{ic_antennas__antenna} = StructArray(ic_antennas__antenna() for k in 1:1)
end

Base.@kwdef mutable struct hard_x_rays__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__radiance
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct hard_x_rays__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__etendue_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct hard_x_rays__channel__energy_band
    detection_efficiency :: Array{Float64, 1} = zeros(Float64,(0))
    lower_bound :: Float64 = 0.0
    energies :: Array{Float64, 1} = zeros(Float64,(0))
    upper_bound :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct hard_x_rays__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct hard_x_rays__channel__line_of_sight
    first_point :: hard_x_rays__channel__line_of_sight__first_point = hard_x_rays__channel__line_of_sight__first_point()
    second_point :: hard_x_rays__channel__line_of_sight__second_point = hard_x_rays__channel__line_of_sight__second_point()
end

Base.@kwdef mutable struct hard_x_rays__channel__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct hard_x_rays__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__detector
    geometry_type :: Int32 = 0
    x2_width :: Float64 = 0.0
    radius :: Float64 = 0.0
    outline :: hard_x_rays__channel__detector__outline = hard_x_rays__channel__detector__outline()
    x1_width :: Float64 = 0.0
    x3_unit_vector :: hard_x_rays__channel__detector__x3_unit_vector = hard_x_rays__channel__detector__x3_unit_vector()
    x2_unit_vector :: hard_x_rays__channel__detector__x2_unit_vector = hard_x_rays__channel__detector__x2_unit_vector()
    centre :: hard_x_rays__channel__detector__centre = hard_x_rays__channel__detector__centre()
    surface :: Float64 = 0.0
    x1_unit_vector :: hard_x_rays__channel__detector__x1_unit_vector = hard_x_rays__channel__detector__x1_unit_vector()
end

Base.@kwdef mutable struct hard_x_rays__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__code
    library :: StructArray{hard_x_rays__code__library} = StructArray(hard_x_rays__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct hard_x_rays__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct hard_x_rays__ids_properties
    provider :: String = ""
    version_put :: hard_x_rays__ids_properties__version_put = hard_x_rays__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct hard_x_rays__emissivity_profile_1d
    peak_position :: Array{Float64, 1} = zeros(Float64,(0))
    half_width_internal :: Array{Float64, 1} = zeros(Float64,(0))
    half_width_external :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
    lower_bound :: Float64 = 0.0
    upper_bound :: Float64 = 0.0
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
    emissivity :: Array{Float64, 2} = zeros(Float64,(0,0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct hard_x_rays__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct hard_x_rays__channel__aperture
    radius :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: hard_x_rays__channel__aperture__centre = hard_x_rays__channel__aperture__centre()
    outline :: hard_x_rays__channel__aperture__outline = hard_x_rays__channel__aperture__outline()
    x1_width :: Float64 = 0.0
    x2_unit_vector :: hard_x_rays__channel__aperture__x2_unit_vector = hard_x_rays__channel__aperture__x2_unit_vector()
    surface :: Float64 = 0.0
    x3_unit_vector :: hard_x_rays__channel__aperture__x3_unit_vector = hard_x_rays__channel__aperture__x3_unit_vector()
    geometry_type :: Int32 = 0
    x1_unit_vector :: hard_x_rays__channel__aperture__x1_unit_vector = hard_x_rays__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct hard_x_rays__channel
    name :: String = ""
    radiance :: hard_x_rays__channel__radiance = hard_x_rays__channel__radiance()
    aperture :: StructArray{hard_x_rays__channel__aperture} = StructArray(hard_x_rays__channel__aperture() for k in 1:1)
    line_of_sight :: hard_x_rays__channel__line_of_sight = hard_x_rays__channel__line_of_sight()
    etendue :: Float64 = 0.0
    energy_band :: StructArray{hard_x_rays__channel__energy_band} = StructArray(hard_x_rays__channel__energy_band() for k in 1:1)
    etendue_method :: hard_x_rays__channel__etendue_method = hard_x_rays__channel__etendue_method()
    identifier :: String = ""
    detector :: hard_x_rays__channel__detector = hard_x_rays__channel__detector()
end

Base.@kwdef mutable struct hard_x_rays
    latency :: Float64 = 0.0
    channel :: StructArray{hard_x_rays__channel} = StructArray(hard_x_rays__channel() for k in 1:1)
    emissivity_profile_1d :: StructArray{hard_x_rays__emissivity_profile_1d} = StructArray(hard_x_rays__emissivity_profile_1d() for k in 1:1)
    ids_properties :: hard_x_rays__ids_properties = hard_x_rays__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: hard_x_rays__code = hard_x_rays__code()
end

Base.@kwdef mutable struct gyrokinetics__collisions
    collisionality_norm :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct gyrokinetics__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct gyrokinetics__species_all
    debye_length_reference :: Float64 = 0.0
    beta_reference :: Float64 = 0.0
    shearing_rate_norm :: Float64 = 0.0
    velocity_tor_norm :: Float64 = 0.0
    zeff :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_gyrocenter
    energy_b_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_a_field_parallel :: Float64 = 0.0
    particles_a_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_a_field_parallel :: Float64 = 0.0
    energy_phi_potential :: Float64 = 0.0
    momentum_tor_parallel_phi_potential :: Float64 = 0.0
    particles_phi_potential :: Float64 = 0.0
    energy_a_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_phi_potential :: Float64 = 0.0
    particles_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_b_field_parallel :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct gyrokinetics__normalizing_quantities
    t_e :: Float64 = 0.0
    n_e :: Float64 = 0.0
    b_field_tor :: Float64 = 0.0
    r :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__code
    library :: StructArray{gyrokinetics__code__library} = StructArray(gyrokinetics__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct gyrokinetics__model
    include_b_field_parallel :: Int32 = 0
    include_full_curvature_drift :: Int32 = 0
    include_centrifugal_effects :: Int32 = 0
    collisions_momentum_conservation :: Int32 = 0
    collisions_finite_larmor_radius :: Int32 = 0
    non_linear_run :: Int32 = 0
    include_a_field_parallel :: Int32 = 0
    collisions_pitch_only :: Int32 = 0
    initial_value_run :: Int32 = 0
    collisions_energy_conservation :: Int32 = 0
    time_interval_norm :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_gyrocenter_rotating_frame
    energy_b_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_a_field_parallel :: Float64 = 0.0
    particles_a_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_a_field_parallel :: Float64 = 0.0
    energy_phi_potential :: Float64 = 0.0
    momentum_tor_parallel_phi_potential :: Float64 = 0.0
    particles_phi_potential :: Float64 = 0.0
    energy_a_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_phi_potential :: Float64 = 0.0
    particles_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_b_field_parallel :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__tag
    name :: String = ""
    comment :: String = ""
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_particle_rotating_frame
    energy_b_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_a_field_parallel :: Float64 = 0.0
    particles_a_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_a_field_parallel :: Float64 = 0.0
    energy_phi_potential :: Float64 = 0.0
    momentum_tor_parallel_phi_potential :: Float64 = 0.0
    particles_phi_potential :: Float64 = 0.0
    energy_a_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_phi_potential :: Float64 = 0.0
    particles_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_b_field_parallel :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__ids_properties
    provider :: String = ""
    version_put :: gyrokinetics__ids_properties__version_put = gyrokinetics__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct gyrokinetics__species
    density_log_gradient_norm :: Float64 = 0.0
    velocity_tor_gradient_norm :: Float64 = 0.0
    density_norm :: Float64 = 0.0
    charge_norm :: Float64 = 0.0
    mass_norm :: Float64 = 0.0
    temperature_log_gradient_norm :: Float64 = 0.0
    temperature_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__fluxes_integrated_norm
    energy_b_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_a_field_parallel :: Float64 = 0.0
    particles_a_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_a_field_parallel :: Float64 = 0.0
    energy_phi_potential :: Float64 = 0.0
    momentum_tor_parallel_phi_potential :: Float64 = 0.0
    particles_phi_potential :: Float64 = 0.0
    energy_a_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_phi_potential :: Float64 = 0.0
    particles_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_b_field_parallel :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments__moments_norm_particle
    density :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    v_parallel_energy_perpendicular :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    heat_flux_parallel :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    v_perpendicular_square_energy :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    pressure_parallel :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    j_parallel :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    pressure_perpendicular :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments__moments_norm_gyrocenter
    v_parallel_energy_perpendicular_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    pressure_perpendicular_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    j_parallel_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    v_perpendicular_square_energy_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    heat_flux_parallel_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    pressure_perpendicular :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    density_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    v_parallel_energy_perpendicular :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    pressure_parallel_gyroav :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    density :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    heat_flux_parallel :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    v_perpendicular_square_energy :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    pressure_parallel :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    j_parallel :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
end

Base.@kwdef mutable struct gyrokinetics__flux_surface
    shape_coefficients_c :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_gradient_norm :: Float64 = 0.0
    triangularity_upper :: Float64 = 0.0
    q :: Float64 = 0.0
    triangularity_lower :: Float64 = 0.0
    r_minor_norm :: Float64 = 0.0
    shape_coefficients_s :: Array{Float64, 1} = zeros(Float64,(0))
    ip_sign :: Float64 = 0.0
    magnetic_shear_r_minor :: Float64 = 0.0
    elongation :: Float64 = 0.0
    b_field_tor_sign :: Float64 = 0.0
    ds_dr_minor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    dc_dr_minor_norm :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_particle
    energy_b_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_a_field_parallel :: Float64 = 0.0
    particles_a_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_a_field_parallel :: Float64 = 0.0
    energy_phi_potential :: Float64 = 0.0
    momentum_tor_parallel_phi_potential :: Float64 = 0.0
    particles_phi_potential :: Float64 = 0.0
    energy_a_field_parallel :: Float64 = 0.0
    momentum_tor_parallel_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_phi_potential :: Float64 = 0.0
    particles_b_field_parallel :: Float64 = 0.0
    momentum_tor_perpendicular_b_field_parallel :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode__fluxes_moments
    fluxes_norm_particle_rotating_frame :: gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_particle_rotating_frame = gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_particle_rotating_frame()
    moments_norm_gyrocenter :: gyrokinetics__wavevector__eigenmode__fluxes_moments__moments_norm_gyrocenter = gyrokinetics__wavevector__eigenmode__fluxes_moments__moments_norm_gyrocenter()
    fluxes_norm_particle :: gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_particle = gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_particle()
    fluxes_norm_gyrocenter :: gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_gyrocenter = gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_gyrocenter()
    moments_norm_particle :: gyrokinetics__wavevector__eigenmode__fluxes_moments__moments_norm_particle = gyrokinetics__wavevector__eigenmode__fluxes_moments__moments_norm_particle()
    fluxes_norm_gyrocenter_rotating_frame :: gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_gyrocenter_rotating_frame = gyrokinetics__wavevector__eigenmode__fluxes_moments__fluxes_norm_gyrocenter_rotating_frame()
end

Base.@kwdef mutable struct gyrokinetics__wavevector__eigenmode
    growth_rate_norm :: Float64 = 0.0
    frequency_norm :: Float64 = 0.0
    phi_potential_perturbed_weight :: Array{Float64, 1} = zeros(Float64,(0))
    time_norm :: Array{Float64, 1} = zeros(Float64,(0))
    phi_potential_perturbed_norm :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    b_field_parallel_perturbed_parity :: Array{Float64, 1} = zeros(Float64,(0))
    fluxes_moments :: StructArray{gyrokinetics__wavevector__eigenmode__fluxes_moments} = StructArray(gyrokinetics__wavevector__eigenmode__fluxes_moments() for k in 1:1)
    a_field_parallel_perturbed_norm :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
    poloidal_angle :: Array{Float64, 1} = zeros(Float64,(0))
    growth_rate_tolerance :: Float64 = 0.0
    a_field_parallel_perturbed_weight :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_parallel_perturbed_weight :: Array{Float64, 1} = zeros(Float64,(0))
    a_field_parallel_perturbed_parity :: Array{Float64, 1} = zeros(Float64,(0))
    phi_potential_perturbed_parity :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_parallel_perturbed_norm :: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))
end

Base.@kwdef mutable struct gyrokinetics__wavevector
    eigenmode :: StructArray{gyrokinetics__wavevector__eigenmode} = StructArray(gyrokinetics__wavevector__eigenmode() for k in 1:1)
    poloidal_turns :: Int32 = 0
    radial_component_norm :: Float64 = 0.0
    binormal_component_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct gyrokinetics
    time :: Array{Float64, 1} = zeros(Float64,(0))
    model :: gyrokinetics__model = gyrokinetics__model()
    code :: gyrokinetics__code = gyrokinetics__code()
    normalizing_quantities :: gyrokinetics__normalizing_quantities = gyrokinetics__normalizing_quantities()
    tag :: StructArray{gyrokinetics__tag} = StructArray(gyrokinetics__tag() for k in 1:1)
    flux_surface :: gyrokinetics__flux_surface = gyrokinetics__flux_surface()
    ids_properties :: gyrokinetics__ids_properties = gyrokinetics__ids_properties()
    fluxes_integrated_norm :: StructArray{gyrokinetics__fluxes_integrated_norm} = StructArray(gyrokinetics__fluxes_integrated_norm() for k in 1:1)
    wavevector :: StructArray{gyrokinetics__wavevector} = StructArray(gyrokinetics__wavevector() for k in 1:1)
    collisions :: gyrokinetics__collisions = gyrokinetics__collisions()
    species :: StructArray{gyrokinetics__species} = StructArray(gyrokinetics__species() for k in 1:1)
    species_all :: gyrokinetics__species_all = gyrokinetics__species_all()
end

Base.@kwdef mutable struct gas_pumping__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct gas_pumping__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct gas_pumping__duct__species__flow_rate
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gas_pumping__duct__flow_rate
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gas_pumping__code
    library :: StructArray{gas_pumping__code__library} = StructArray(gas_pumping__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct gas_pumping__ids_properties
    provider :: String = ""
    version_put :: gas_pumping__ids_properties__version_put = gas_pumping__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct gas_pumping__duct__species__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct gas_pumping__duct__species
    label :: String = ""
    flow_rate :: gas_pumping__duct__species__flow_rate = gas_pumping__duct__species__flow_rate()
    element :: StructArray{gas_pumping__duct__species__element} = StructArray(gas_pumping__duct__species__element() for k in 1:1)
end

Base.@kwdef mutable struct gas_pumping__duct
    name :: String = ""
    flow_rate :: gas_pumping__duct__flow_rate = gas_pumping__duct__flow_rate()
    species :: StructArray{gas_pumping__duct__species} = StructArray(gas_pumping__duct__species() for k in 1:1)
    identifier :: String = ""
end

Base.@kwdef mutable struct gas_pumping
    duct :: StructArray{gas_pumping__duct} = StructArray(gas_pumping__duct() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: gas_pumping__ids_properties = gas_pumping__ids_properties()
    code :: gas_pumping__code = gas_pumping__code()
end

Base.@kwdef mutable struct gas_injection__pipe__flow_rate
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gas_injection__pipe__valve__flow_rate
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gas_injection__pipe__species__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct gas_injection__pipe__valve__electron_rate
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct gas_injection__pipe__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct gas_injection__pipe__species
    label :: String = ""
    fraction :: Float64 = 0.0
    element :: StructArray{gas_injection__pipe__species__element} = StructArray(gas_injection__pipe__species__element() for k in 1:1)
end

Base.@kwdef mutable struct gas_injection__pipe__valve__species__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct gas_injection__pipe__valve__species
    label :: String = ""
    fraction :: Float64 = 0.0
    element :: StructArray{gas_injection__pipe__valve__species__element} = StructArray(gas_injection__pipe__valve__species__element() for k in 1:1)
end

Base.@kwdef mutable struct gas_injection__pipe__exit_position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct gas_injection__pipe__valve
    flow_rate_max :: Float64 = 0.0
    name :: String = ""
    flow_rate :: gas_injection__pipe__valve__flow_rate = gas_injection__pipe__valve__flow_rate()
    species :: StructArray{gas_injection__pipe__valve__species} = StructArray(gas_injection__pipe__valve__species() for k in 1:1)
    identifier :: String = ""
    flow_rate_min :: Float64 = 0.0
    electron_rate :: gas_injection__pipe__valve__electron_rate = gas_injection__pipe__valve__electron_rate()
end

Base.@kwdef mutable struct gas_injection__pipe
    length :: Float64 = 0.0
    exit_position :: gas_injection__pipe__exit_position = gas_injection__pipe__exit_position()
    name :: String = ""
    flow_rate :: gas_injection__pipe__flow_rate = gas_injection__pipe__flow_rate()
    second_point :: gas_injection__pipe__second_point = gas_injection__pipe__second_point()
    valve :: StructArray{gas_injection__pipe__valve} = StructArray(gas_injection__pipe__valve() for k in 1:1)
    species :: StructArray{gas_injection__pipe__species} = StructArray(gas_injection__pipe__species() for k in 1:1)
    identifier :: String = ""
end

Base.@kwdef mutable struct gas_injection__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct gas_injection__code
    library :: StructArray{gas_injection__code__library} = StructArray(gas_injection__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct gas_injection__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct gas_injection__ids_properties
    provider :: String = ""
    version_put :: gas_injection__ids_properties__version_put = gas_injection__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct gas_injection
    latency :: Float64 = 0.0
    ids_properties :: gas_injection__ids_properties = gas_injection__ids_properties()
    pipe :: StructArray{gas_injection__pipe} = StructArray(gas_injection__pipe() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: gas_injection__code = gas_injection__code()
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__x_point__position_measured
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__b_field_r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_2d__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities__current_centre
    velocity_z :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__geometric_axis
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__strike_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_1d__geometric_axis
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__strike_point__position_measured
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__j_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__geometric_axis
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix__x_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__pressure
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__gap
    name :: String = ""
    r :: Float64 = 0.0
    value :: Float64 = 0.0
    identifier :: String = ""
    angle :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__active_limiter_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary} = StructArray(equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__objects_per_dimension
    object :: StructArray{equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object} = StructArray(equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__coordinate_system__grid_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__x_point__position_reconstructed
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__b_field_z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__j_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__phi
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__z
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__x_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__strike_point__position_reconstructed
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__strike_point
    chi_squared_z :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    position_measured :: equilibrium__time_slice__constraints__strike_point__position_measured = equilibrium__time_slice__constraints__strike_point__position_measured()
    time_measurement :: Float64 = 0.0
    chi_squared_r :: Float64 = 0.0
    position_reconstructed :: equilibrium__time_slice__constraints__strike_point__position_reconstructed = equilibrium__time_slice__constraints__strike_point__position_reconstructed()
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__bpol_probe
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__x_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__element
    object :: StructArray{equilibrium__time_slice__ggd__grid__grid_subset__element__object} = StructArray(equilibrium__time_slice__ggd__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct equilibrium__code
    library :: StructArray{equilibrium__code__library} = StructArray(equilibrium__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary} = StructArray(equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix__strike_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__lcfs
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities__q_min
    value :: Float64 = 0.0
    rho_tor_norm :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__flux_loop
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__r
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__closest_wall_point
    distance :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__faraday_angle
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__b_field_tor
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__diamagnetic_flux
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    measured :: Float64 = 0.0
    reconstructed :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__outline
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__ip
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__active_limiter_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__element
    object :: StructArray{equilibrium__grids_ggd__grid__grid_subset__element__object} = StructArray(equilibrium__grids_ggd__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__objects_per_dimension
    object :: StructArray{equilibrium__grids_ggd__grid__space__objects_per_dimension__object} = StructArray(equilibrium__grids_ggd__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    measured :: Float64 = 0.0
    reconstructed :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__iron_core_segment
    magnetisation_r :: equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r = equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r()
    magnetisation_z :: equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z = equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z()
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix
    psi :: Float64 = 0.0
    x_point :: StructArray{equilibrium__time_slice__boundary_secondary_separatrix__x_point} = StructArray(equilibrium__time_slice__boundary_secondary_separatrix__x_point() for k in 1:1)
    distance_inner_outer :: Float64 = 0.0
    outline :: equilibrium__time_slice__boundary_secondary_separatrix__outline = equilibrium__time_slice__boundary_secondary_separatrix__outline()
    strike_point :: StructArray{equilibrium__time_slice__boundary_secondary_separatrix__strike_point} = StructArray(equilibrium__time_slice__boundary_secondary_separatrix__strike_point() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__b_field_tor_vacuum_r
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix
    psi :: Float64 = 0.0
    elongation_lower :: Float64 = 0.0
    strike_point :: StructArray{equilibrium__time_slice__boundary_separatrix__strike_point} = StructArray(equilibrium__time_slice__boundary_separatrix__strike_point() for k in 1:1)
    x_point :: StructArray{equilibrium__time_slice__boundary_separatrix__x_point} = StructArray(equilibrium__time_slice__boundary_separatrix__x_point() for k in 1:1)
    triangularity :: Float64 = 0.0
    gap :: StructArray{equilibrium__time_slice__boundary_separatrix__gap} = StructArray(equilibrium__time_slice__boundary_separatrix__gap() for k in 1:1)
    elongation_upper :: Float64 = 0.0
    triangularity_upper :: Float64 = 0.0
    outline :: equilibrium__time_slice__boundary_separatrix__outline = equilibrium__time_slice__boundary_separatrix__outline()
    dr_dz_zero_point :: equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point = equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point()
    squareness_lower_outer :: Float64 = 0.0
    triangularity_lower :: Float64 = 0.0
    minor_radius :: Float64 = 0.0
    squareness_upper_inner :: Float64 = 0.0
    squareness_upper_outer :: Float64 = 0.0
    squareness_lower_inner :: Float64 = 0.0
    geometric_axis :: equilibrium__time_slice__boundary_separatrix__geometric_axis = equilibrium__time_slice__boundary_separatrix__geometric_axis()
    elongation :: Float64 = 0.0
    active_limiter_point :: equilibrium__time_slice__boundary_separatrix__active_limiter_point = equilibrium__time_slice__boundary_separatrix__active_limiter_point()
    closest_wall_point :: equilibrium__time_slice__boundary_separatrix__closest_wall_point = equilibrium__time_slice__boundary_separatrix__closest_wall_point()
    type :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__n_e_line
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    measured :: Float64 = 0.0
    reconstructed :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__n_e
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    measured :: Float64 = 0.0
    reconstructed :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{equilibrium__grids_ggd__grid__space__objects_per_dimension} = StructArray(equilibrium__grids_ggd__grid__space__objects_per_dimension() for k in 1:1)
    identifier :: equilibrium__grids_ggd__grid__space__identifier = equilibrium__grids_ggd__grid__space__identifier()
    geometry_type :: equilibrium__grids_ggd__grid__space__geometry_type = equilibrium__grids_ggd__grid__space__geometry_type()
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_2d__grid
    volume_element :: Array{Float64, 2} = zeros(Float64,(0,0))
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__mse_polarisation_angle
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__theta
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__convergence
    iterations_n :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__psi
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    geometry_type :: equilibrium__time_slice__ggd__grid__space__geometry_type = equilibrium__time_slice__ggd__grid__space__geometry_type()
    identifier :: equilibrium__time_slice__ggd__grid__space__identifier = equilibrium__time_slice__ggd__grid__space__identifier()
    objects_per_dimension :: StructArray{equilibrium__time_slice__ggd__grid__space__objects_per_dimension} = StructArray(equilibrium__time_slice__ggd__grid__space__objects_per_dimension() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__pf_current
    chi_squared :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    measured :: Float64 = 0.0
    reconstructed :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset
    base :: StructArray{equilibrium__time_slice__ggd__grid__grid_subset__base} = StructArray(equilibrium__time_slice__ggd__grid__grid_subset__base() for k in 1:1)
    metric :: equilibrium__time_slice__ggd__grid__grid_subset__metric = equilibrium__time_slice__ggd__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: equilibrium__time_slice__ggd__grid__grid_subset__identifier = equilibrium__time_slice__ggd__grid__grid_subset__identifier()
    element :: StructArray{equilibrium__time_slice__ggd__grid__grid_subset__element} = StructArray(equilibrium__time_slice__ggd__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid
    grid_subset :: StructArray{equilibrium__time_slice__ggd__grid__grid_subset} = StructArray(equilibrium__time_slice__ggd__grid__grid_subset() for k in 1:1)
    space :: StructArray{equilibrium__time_slice__ggd__grid__space} = StructArray(equilibrium__time_slice__ggd__grid__space() for k in 1:1)
    identifier :: equilibrium__time_slice__ggd__grid__identifier = equilibrium__time_slice__ggd__grid__identifier()
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_2d
    psi :: Array{Float64, 2} = zeros(Float64,(0,0))
    b_field_r :: Array{Float64, 2} = zeros(Float64,(0,0))
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    b_r :: Array{Float64, 2} = zeros(Float64,(0,0))
    theta :: Array{Float64, 2} = zeros(Float64,(0,0))
    b_field_z :: Array{Float64, 2} = zeros(Float64,(0,0))
    j_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    phi :: Array{Float64, 2} = zeros(Float64,(0,0))
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid :: equilibrium__time_slice__profiles_2d__grid = equilibrium__time_slice__profiles_2d__grid()
    b_z :: Array{Float64, 2} = zeros(Float64,(0,0))
    b_field_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_type :: equilibrium__time_slice__profiles_2d__grid_type = equilibrium__time_slice__profiles_2d__grid_type()
    j_parallel :: Array{Float64, 2} = zeros(Float64,(0,0))
    b_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__pf_passive_current
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    measured :: Float64 = 0.0
    reconstructed :: Float64 = 0.0
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__strike_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary
    psi :: Float64 = 0.0
    lcfs :: equilibrium__time_slice__boundary__lcfs = equilibrium__time_slice__boundary__lcfs()
    elongation_lower :: Float64 = 0.0
    strike_point :: StructArray{equilibrium__time_slice__boundary__strike_point} = StructArray(equilibrium__time_slice__boundary__strike_point() for k in 1:1)
    x_point :: StructArray{equilibrium__time_slice__boundary__x_point} = StructArray(equilibrium__time_slice__boundary__x_point() for k in 1:1)
    triangularity :: Float64 = 0.0
    elongation_upper :: Float64 = 0.0
    triangularity_upper :: Float64 = 0.0
    outline :: equilibrium__time_slice__boundary__outline = equilibrium__time_slice__boundary__outline()
    squareness_lower_outer :: Float64 = 0.0
    triangularity_lower :: Float64 = 0.0
    psi_norm :: Float64 = 0.0
    minor_radius :: Float64 = 0.0
    squareness_upper_inner :: Float64 = 0.0
    squareness_upper_outer :: Float64 = 0.0
    squareness_lower_inner :: Float64 = 0.0
    geometric_axis :: equilibrium__time_slice__boundary__geometric_axis = equilibrium__time_slice__boundary__geometric_axis()
    elongation :: Float64 = 0.0
    active_limiter_point :: equilibrium__time_slice__boundary__active_limiter_point = equilibrium__time_slice__boundary__active_limiter_point()
    b_flux_pol_norm :: Float64 = 0.0
    type :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset
    base :: StructArray{equilibrium__grids_ggd__grid__grid_subset__base} = StructArray(equilibrium__grids_ggd__grid__grid_subset__base() for k in 1:1)
    metric :: equilibrium__grids_ggd__grid__grid_subset__metric = equilibrium__grids_ggd__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: equilibrium__grids_ggd__grid__grid_subset__identifier = equilibrium__grids_ggd__grid__grid_subset__identifier()
    element :: StructArray{equilibrium__grids_ggd__grid__grid_subset__element} = StructArray(equilibrium__grids_ggd__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid
    grid_subset :: StructArray{equilibrium__grids_ggd__grid__grid_subset} = StructArray(equilibrium__grids_ggd__grid__grid_subset() for k in 1:1)
    space :: StructArray{equilibrium__grids_ggd__grid__space} = StructArray(equilibrium__grids_ggd__grid__space() for k in 1:1)
    identifier :: equilibrium__grids_ggd__grid__identifier = equilibrium__grids_ggd__grid__identifier()
end

Base.@kwdef mutable struct equilibrium__grids_ggd
    time :: Float64 = 0.0
    grid :: StructArray{equilibrium__grids_ggd__grid} = StructArray(equilibrium__grids_ggd__grid() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities__magnetic_axis
    r :: Float64 = 0.0
    b_field_tor :: Float64 = 0.0
    z :: Float64 = 0.0
    b_tor :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities
    ip :: Float64 = 0.0
    li_3 :: Float64 = 0.0
    beta_tor :: Float64 = 0.0
    surface :: Float64 = 0.0
    magnetic_axis :: equilibrium__time_slice__global_quantities__magnetic_axis = equilibrium__time_slice__global_quantities__magnetic_axis()
    energy_mhd :: Float64 = 0.0
    psi_boundary :: Float64 = 0.0
    length_pol :: Float64 = 0.0
    area :: Float64 = 0.0
    psi_external_average :: Float64 = 0.0
    q_95 :: Float64 = 0.0
    q_axis :: Float64 = 0.0
    psi_axis :: Float64 = 0.0
    w_mhd :: Float64 = 0.0
    volume :: Float64 = 0.0
    plasma_inductance :: Float64 = 0.0
    beta_pol :: Float64 = 0.0
    beta_normal :: Float64 = 0.0
    current_centre :: equilibrium__time_slice__global_quantities__current_centre = equilibrium__time_slice__global_quantities__current_centre()
    q_min :: equilibrium__time_slice__global_quantities__q_min = equilibrium__time_slice__global_quantities__q_min()
end

Base.@kwdef mutable struct equilibrium__ids_properties
    provider :: String = ""
    version_put :: equilibrium__ids_properties__version_put = equilibrium__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__x_point
    chi_squared_z :: Float64 = 0.0
    exact :: Int32 = 0
    weight :: Float64 = 0.0
    source :: String = ""
    position_measured :: equilibrium__time_slice__constraints__x_point__position_measured = equilibrium__time_slice__constraints__x_point__position_measured()
    time_measurement :: Float64 = 0.0
    chi_squared_r :: Float64 = 0.0
    position_reconstructed :: equilibrium__time_slice__constraints__x_point__position_reconstructed = equilibrium__time_slice__constraints__x_point__position_reconstructed()
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd
    b_field_z :: StructArray{equilibrium__time_slice__ggd__b_field_z} = StructArray(equilibrium__time_slice__ggd__b_field_z() for k in 1:1)
    psi :: StructArray{equilibrium__time_slice__ggd__psi} = StructArray(equilibrium__time_slice__ggd__psi() for k in 1:1)
    theta :: StructArray{equilibrium__time_slice__ggd__theta} = StructArray(equilibrium__time_slice__ggd__theta() for k in 1:1)
    z :: StructArray{equilibrium__time_slice__ggd__z} = StructArray(equilibrium__time_slice__ggd__z() for k in 1:1)
    phi :: StructArray{equilibrium__time_slice__ggd__phi} = StructArray(equilibrium__time_slice__ggd__phi() for k in 1:1)
    j_tor :: StructArray{equilibrium__time_slice__ggd__j_tor} = StructArray(equilibrium__time_slice__ggd__j_tor() for k in 1:1)
    grid :: equilibrium__time_slice__ggd__grid = equilibrium__time_slice__ggd__grid()
    b_field_tor :: StructArray{equilibrium__time_slice__ggd__b_field_tor} = StructArray(equilibrium__time_slice__ggd__b_field_tor() for k in 1:1)
    b_field_r :: StructArray{equilibrium__time_slice__ggd__b_field_r} = StructArray(equilibrium__time_slice__ggd__b_field_r() for k in 1:1)
    r :: StructArray{equilibrium__time_slice__ggd__r} = StructArray(equilibrium__time_slice__ggd__r() for k in 1:1)
    j_parallel :: StructArray{equilibrium__time_slice__ggd__j_parallel} = StructArray(equilibrium__time_slice__ggd__j_parallel() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__coordinate_system__grid
    volume_element :: Array{Float64, 2} = zeros(Float64,(0,0))
    dim2 :: Array{Float64, 1} = zeros(Float64,(0))
    dim1 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct equilibrium__time_slice__coordinate_system
    jacobian :: Array{Float64, 2} = zeros(Float64,(0,0))
    g13_covariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    g11_contravariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    g13_contravariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    r :: Array{Float64, 2} = zeros(Float64,(0,0))
    g12_contravariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    g22_contravariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    z :: Array{Float64, 2} = zeros(Float64,(0,0))
    g33_contravariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    g22_covariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    tensor_contravariant :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
    g33_covariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    g12_covariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid :: equilibrium__time_slice__coordinate_system__grid = equilibrium__time_slice__coordinate_system__grid()
    g23_covariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    g11_covariant :: Array{Float64, 2} = zeros(Float64,(0,0))
    tensor_covariant :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
    grid_type :: equilibrium__time_slice__coordinate_system__grid_type = equilibrium__time_slice__coordinate_system__grid_type()
    g23_contravariant :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__q__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__q
    chi_squared :: Float64 = 0.0
    weight :: Float64 = 0.0
    exact :: Int32 = 0
    source :: String = ""
    reconstructed :: Float64 = 0.0
    measured :: Float64 = 0.0
    position :: equilibrium__time_slice__constraints__q__position = equilibrium__time_slice__constraints__q__position()
    time_measurement :: Float64 = 0.0
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints
    flux_loop :: StructArray{equilibrium__time_slice__constraints__flux_loop} = StructArray(equilibrium__time_slice__constraints__flux_loop() for k in 1:1)
    n_e :: StructArray{equilibrium__time_slice__constraints__n_e} = StructArray(equilibrium__time_slice__constraints__n_e() for k in 1:1)
    ip :: equilibrium__time_slice__constraints__ip = equilibrium__time_slice__constraints__ip()
    n_e_line :: StructArray{equilibrium__time_slice__constraints__n_e_line} = StructArray(equilibrium__time_slice__constraints__n_e_line() for k in 1:1)
    pf_current :: StructArray{equilibrium__time_slice__constraints__pf_current} = StructArray(equilibrium__time_slice__constraints__pf_current() for k in 1:1)
    strike_point :: StructArray{equilibrium__time_slice__constraints__strike_point} = StructArray(equilibrium__time_slice__constraints__strike_point() for k in 1:1)
    x_point :: StructArray{equilibrium__time_slice__constraints__x_point} = StructArray(equilibrium__time_slice__constraints__x_point() for k in 1:1)
    iron_core_segment :: StructArray{equilibrium__time_slice__constraints__iron_core_segment} = StructArray(equilibrium__time_slice__constraints__iron_core_segment() for k in 1:1)
    pressure :: StructArray{equilibrium__time_slice__constraints__pressure} = StructArray(equilibrium__time_slice__constraints__pressure() for k in 1:1)
    diamagnetic_flux :: equilibrium__time_slice__constraints__diamagnetic_flux = equilibrium__time_slice__constraints__diamagnetic_flux()
    pf_passive_current :: StructArray{equilibrium__time_slice__constraints__pf_passive_current} = StructArray(equilibrium__time_slice__constraints__pf_passive_current() for k in 1:1)
    bpol_probe :: StructArray{equilibrium__time_slice__constraints__bpol_probe} = StructArray(equilibrium__time_slice__constraints__bpol_probe() for k in 1:1)
    q :: StructArray{equilibrium__time_slice__constraints__q} = StructArray(equilibrium__time_slice__constraints__q() for k in 1:1)
    mse_polarisation_angle :: StructArray{equilibrium__time_slice__constraints__mse_polarisation_angle} = StructArray(equilibrium__time_slice__constraints__mse_polarisation_angle() for k in 1:1)
    b_field_tor_vacuum_r :: equilibrium__time_slice__constraints__b_field_tor_vacuum_r = equilibrium__time_slice__constraints__b_field_tor_vacuum_r()
    faraday_angle :: StructArray{equilibrium__time_slice__constraints__faraday_angle} = StructArray(equilibrium__time_slice__constraints__faraday_angle() for k in 1:1)
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_1d
    dvolume_drho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_max :: Array{Float64, 1} = zeros(Float64,(0))
    gm9 :: Array{Float64, 1} = zeros(Float64,(0))
    dpsi_drho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    magnetic_shear :: Array{Float64, 1} = zeros(Float64,(0))
    b_average :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_min :: Array{Float64, 1} = zeros(Float64,(0))
    darea_dpsi :: Array{Float64, 1} = zeros(Float64,(0))
    gm3 :: Array{Float64, 1} = zeros(Float64,(0))
    squareness_upper_inner :: Array{Float64, 1} = zeros(Float64,(0))
    squareness_lower_inner :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    elongation :: Array{Float64, 1} = zeros(Float64,(0))
    beta_pol :: Array{Float64, 1} = zeros(Float64,(0))
    b_field_average :: Array{Float64, 1} = zeros(Float64,(0))
    j_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    gm6 :: Array{Float64, 1} = zeros(Float64,(0))
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    gm8 :: Array{Float64, 1} = zeros(Float64,(0))
    dpressure_dpsi :: Array{Float64, 1} = zeros(Float64,(0))
    triangularity_upper :: Array{Float64, 1} = zeros(Float64,(0))
    darea_drho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    trapped_fraction :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    dvolume_dpsi :: Array{Float64, 1} = zeros(Float64,(0))
    b_min :: Array{Float64, 1} = zeros(Float64,(0))
    f :: Array{Float64, 1} = zeros(Float64,(0))
    mass_density :: Array{Float64, 1} = zeros(Float64,(0))
    r_outboard :: Array{Float64, 1} = zeros(Float64,(0))
    gm4 :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    squareness_lower_outer :: Array{Float64, 1} = zeros(Float64,(0))
    triangularity_lower :: Array{Float64, 1} = zeros(Float64,(0))
    gm2 :: Array{Float64, 1} = zeros(Float64,(0))
    rho_volume_norm :: Array{Float64, 1} = zeros(Float64,(0))
    gm1 :: Array{Float64, 1} = zeros(Float64,(0))
    gm5 :: Array{Float64, 1} = zeros(Float64,(0))
    b_max :: Array{Float64, 1} = zeros(Float64,(0))
    f_df_dpsi :: Array{Float64, 1} = zeros(Float64,(0))
    j_tor :: Array{Float64, 1} = zeros(Float64,(0))
    r_inboard :: Array{Float64, 1} = zeros(Float64,(0))
    q :: Array{Float64, 1} = zeros(Float64,(0))
    gm7 :: Array{Float64, 1} = zeros(Float64,(0))
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    squareness_upper_outer :: Array{Float64, 1} = zeros(Float64,(0))
    geometric_axis :: equilibrium__time_slice__profiles_1d__geometric_axis = equilibrium__time_slice__profiles_1d__geometric_axis()
end

Base.@kwdef mutable struct equilibrium__time_slice
    time :: Float64 = 0.0
    ggd :: StructArray{equilibrium__time_slice__ggd} = StructArray(equilibrium__time_slice__ggd() for k in 1:1)
    boundary :: equilibrium__time_slice__boundary = equilibrium__time_slice__boundary()
    profiles_1d :: equilibrium__time_slice__profiles_1d = equilibrium__time_slice__profiles_1d()
    constraints :: equilibrium__time_slice__constraints = equilibrium__time_slice__constraints()
    global_quantities :: equilibrium__time_slice__global_quantities = equilibrium__time_slice__global_quantities()
    convergence :: equilibrium__time_slice__convergence = equilibrium__time_slice__convergence()
    coordinate_system :: equilibrium__time_slice__coordinate_system = equilibrium__time_slice__coordinate_system()
    boundary_secondary_separatrix :: equilibrium__time_slice__boundary_secondary_separatrix = equilibrium__time_slice__boundary_secondary_separatrix()
    profiles_2d :: StructArray{equilibrium__time_slice__profiles_2d} = StructArray(equilibrium__time_slice__profiles_2d() for k in 1:1)
    boundary_separatrix :: equilibrium__time_slice__boundary_separatrix = equilibrium__time_slice__boundary_separatrix()
end

Base.@kwdef mutable struct equilibrium
    time_slice :: StructArray{equilibrium__time_slice} = StructArray(equilibrium__time_slice() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: equilibrium__ids_properties = equilibrium__ids_properties()
    grids_ggd :: StructArray{equilibrium__grids_ggd} = StructArray(equilibrium__grids_ggd() for k in 1:1)
    vacuum_toroidal_field :: equilibrium__vacuum_toroidal_field = equilibrium__vacuum_toroidal_field()
    code :: equilibrium__code = equilibrium__code()
end

Base.@kwdef mutable struct em_coupling__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct em_coupling__ids_properties
    provider :: String = ""
    version_put :: em_coupling__ids_properties__version_put = em_coupling__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct em_coupling__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct em_coupling__code
    library :: StructArray{em_coupling__code__library} = StructArray(em_coupling__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct em_coupling
    poloidal_probes :: Array{String, 1} = String[]
    time :: Array{Float64, 1} = zeros(Float64,(0))
    field_probes_grid :: Array{Float64, 2} = zeros(Float64,(0,0))
    code :: em_coupling__code = em_coupling__code()
    field_probes_active :: Array{Float64, 2} = zeros(Float64,(0,0))
    mutual_active_active :: Array{Float64, 2} = zeros(Float64,(0,0))
    mutual_grid_passive :: Array{Float64, 2} = zeros(Float64,(0,0))
    mutual_grid_grid :: Array{Float64, 2} = zeros(Float64,(0,0))
    mutual_passive_active :: Array{Float64, 2} = zeros(Float64,(0,0))
    flux_loops :: Array{String, 1} = String[]
    mutual_loops_grid :: Array{Float64, 2} = zeros(Float64,(0,0))
    field_probes_passive :: Array{Float64, 2} = zeros(Float64,(0,0))
    ids_properties :: em_coupling__ids_properties = em_coupling__ids_properties()
    mutual_grid_active :: Array{Float64, 2} = zeros(Float64,(0,0))
    active_coils :: Array{String, 1} = String[]
    mutual_loops_active :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_points :: Array{String, 1} = String[]
    mutual_passive_passive :: Array{Float64, 2} = zeros(Float64,(0,0))
    passive_loops :: Array{String, 1} = String[]
    mutual_loops_passive :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__particles__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__grid_ggd__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__energy__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__particles__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__energy__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__energy__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__grid_ggd__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__neutral__particle_flux_integrated
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__energy__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__ion__particle_flux_integrated
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__power
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__energy__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__energy__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__power_ion_total
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__particles__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__energy__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__particles__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__grid_ggd__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{edge_transport__grid_ggd__space__objects_per_dimension__object__boundary} = StructArray(edge_transport__grid_ggd__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__energy__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__particles__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__energy_flux_max
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__particles__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__grid_ggd__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__electrons__power
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__particles__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__particles__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__grid_ggd__space__objects_per_dimension
    object :: StructArray{edge_transport__grid_ggd__space__objects_per_dimension__object} = StructArray(edge_transport__grid_ggd__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__grid_ggd__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__model__ggd__total_ion_energy__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__energy__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__momentum__flux_limiter
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__grid_ggd__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__momentum__flux_limiter
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__momentum__d
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__grid_ggd__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct edge_transport__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct edge_transport__model__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__model__ggd__total_ion_energy__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__momentum__d
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__total_ion_energy__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__electrons__particle_flux_integrated
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__momentum__v
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__energy__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__particles__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__code__output_flag
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__particles__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__energy__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__energy
    flux :: StructArray{edge_transport__model__ggd__neutral__state__energy__flux} = StructArray(edge_transport__model__ggd__neutral__state__energy__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__neutral__state__energy__d} = StructArray(edge_transport__model__ggd__neutral__state__energy__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__neutral__state__energy__v} = StructArray(edge_transport__model__ggd__neutral__state__energy__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__neutral__state__energy__flux_limiter} = StructArray(edge_transport__model__ggd__neutral__state__energy__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__momentum__v
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__momentum__flux_limiter
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__momentum__flux
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__energy__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__momentum__flux_limiter
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__particles__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__particles
    flux :: StructArray{edge_transport__model__ggd__ion__state__particles__flux} = StructArray(edge_transport__model__ggd__ion__state__particles__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__ion__state__particles__d} = StructArray(edge_transport__model__ggd__ion__state__particles__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__ion__state__particles__v} = StructArray(edge_transport__model__ggd__ion__state__particles__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__ion__state__particles__flux_limiter} = StructArray(edge_transport__model__ggd__ion__state__particles__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__energy__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__energy__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__momentum__flux
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__particles__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__particles__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__particles__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__particles
    flux :: StructArray{edge_transport__model__ggd__ion__particles__flux} = StructArray(edge_transport__model__ggd__ion__particles__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__ion__particles__d} = StructArray(edge_transport__model__ggd__ion__particles__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__ion__particles__v} = StructArray(edge_transport__model__ggd__ion__particles__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__ion__particles__flux_limiter} = StructArray(edge_transport__model__ggd__ion__particles__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__particles__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__particles__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__particles
    flux :: StructArray{edge_transport__model__ggd__neutral__particles__flux} = StructArray(edge_transport__model__ggd__neutral__particles__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__neutral__particles__d} = StructArray(edge_transport__model__ggd__neutral__particles__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__neutral__particles__v} = StructArray(edge_transport__model__ggd__neutral__particles__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__neutral__particles__flux_limiter} = StructArray(edge_transport__model__ggd__neutral__particles__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__particles__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__code
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: edge_transport__model__code__output_flag = edge_transport__model__code__output_flag()
    version :: String = ""
end

Base.@kwdef mutable struct edge_transport__grid_ggd__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    geometry_type :: edge_transport__grid_ggd__space__geometry_type = edge_transport__grid_ggd__space__geometry_type()
    objects_per_dimension :: StructArray{edge_transport__grid_ggd__space__objects_per_dimension} = StructArray(edge_transport__grid_ggd__space__objects_per_dimension() for k in 1:1)
    identifier :: edge_transport__grid_ggd__space__identifier = edge_transport__grid_ggd__space__identifier()
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__particles__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__energy__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__energy__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__momentum__flux
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__energy__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__energy
    flux :: StructArray{edge_transport__model__ggd__electrons__energy__flux} = StructArray(edge_transport__model__ggd__electrons__energy__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__electrons__energy__d} = StructArray(edge_transport__model__ggd__electrons__energy__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__electrons__energy__v} = StructArray(edge_transport__model__ggd__electrons__energy__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__electrons__energy__flux_limiter} = StructArray(edge_transport__model__ggd__electrons__energy__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__code
    library :: StructArray{edge_transport__code__library} = StructArray(edge_transport__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct edge_transport__grid_ggd__grid_subset__element
    object :: StructArray{edge_transport__grid_ggd__grid_subset__element__object} = StructArray(edge_transport__grid_ggd__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__particles__d
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__particles
    flux :: StructArray{edge_transport__model__ggd__neutral__state__particles__flux} = StructArray(edge_transport__model__ggd__neutral__state__particles__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__neutral__state__particles__d} = StructArray(edge_transport__model__ggd__neutral__state__particles__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__neutral__state__particles__v} = StructArray(edge_transport__model__ggd__neutral__state__particles__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__neutral__state__particles__flux_limiter} = StructArray(edge_transport__model__ggd__neutral__state__particles__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__total_ion_energy__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__total_ion_energy
    flux :: StructArray{edge_transport__model__ggd__total_ion_energy__flux} = StructArray(edge_transport__model__ggd__total_ion_energy__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__total_ion_energy__d} = StructArray(edge_transport__model__ggd__total_ion_energy__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__total_ion_energy__v} = StructArray(edge_transport__model__ggd__total_ion_energy__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__total_ion_energy__flux_limiter} = StructArray(edge_transport__model__ggd__total_ion_energy__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__momentum__flux_limiter
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__momentum__flux
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__momentum__d
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__momentum
    flux :: StructArray{edge_transport__model__ggd__ion__momentum__flux} = StructArray(edge_transport__model__ggd__ion__momentum__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__ion__momentum__d} = StructArray(edge_transport__model__ggd__ion__momentum__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__ion__momentum__v} = StructArray(edge_transport__model__ggd__ion__momentum__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__ion__momentum__flux_limiter} = StructArray(edge_transport__model__ggd__ion__momentum__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__particles__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons__particles
    flux :: StructArray{edge_transport__model__ggd__electrons__particles__flux} = StructArray(edge_transport__model__ggd__electrons__particles__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__electrons__particles__d} = StructArray(edge_transport__model__ggd__electrons__particles__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__electrons__particles__v} = StructArray(edge_transport__model__ggd__electrons__particles__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__electrons__particles__flux_limiter} = StructArray(edge_transport__model__ggd__electrons__particles__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__electrons
    particles :: edge_transport__model__ggd__electrons__particles = edge_transport__model__ggd__electrons__particles()
    energy :: edge_transport__model__ggd__electrons__energy = edge_transport__model__ggd__electrons__energy()
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    particle_flux_integrated :: StructArray{edge_transport__model__ggd_fast__ion__particle_flux_integrated} = StructArray(edge_transport__model__ggd_fast__ion__particle_flux_integrated() for k in 1:1)
    neutral_index :: Int32 = 0
    element :: StructArray{edge_transport__model__ggd_fast__ion__element} = StructArray(edge_transport__model__ggd_fast__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct edge_transport__ids_properties
    provider :: String = ""
    version_put :: edge_transport__ids_properties__version_put = edge_transport__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__grid_ggd__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__grid_ggd__grid_subset
    base :: StructArray{edge_transport__grid_ggd__grid_subset__base} = StructArray(edge_transport__grid_ggd__grid_subset__base() for k in 1:1)
    metric :: edge_transport__grid_ggd__grid_subset__metric = edge_transport__grid_ggd__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: edge_transport__grid_ggd__grid_subset__identifier = edge_transport__grid_ggd__grid_subset__identifier()
    element :: StructArray{edge_transport__grid_ggd__grid_subset__element} = StructArray(edge_transport__grid_ggd__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__momentum__v
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state__momentum
    flux :: StructArray{edge_transport__model__ggd__neutral__state__momentum__flux} = StructArray(edge_transport__model__ggd__neutral__state__momentum__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__neutral__state__momentum__d} = StructArray(edge_transport__model__ggd__neutral__state__momentum__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__neutral__state__momentum__v} = StructArray(edge_transport__model__ggd__neutral__state__momentum__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__neutral__state__momentum__flux_limiter} = StructArray(edge_transport__model__ggd__neutral__state__momentum__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    particles :: edge_transport__model__ggd__neutral__state__particles = edge_transport__model__ggd__neutral__state__particles()
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    momentum :: edge_transport__model__ggd__neutral__state__momentum = edge_transport__model__ggd__neutral__state__momentum()
    energy :: edge_transport__model__ggd__neutral__state__energy = edge_transport__model__ggd__neutral__state__energy()
    neutral_type :: edge_transport__model__ggd__neutral__state__neutral_type = edge_transport__model__ggd__neutral__state__neutral_type()
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__energy__v
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__energy
    flux :: StructArray{edge_transport__model__ggd__ion__state__energy__flux} = StructArray(edge_transport__model__ggd__ion__state__energy__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__ion__state__energy__d} = StructArray(edge_transport__model__ggd__ion__state__energy__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__ion__state__energy__v} = StructArray(edge_transport__model__ggd__ion__state__energy__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__ion__state__energy__flux_limiter} = StructArray(edge_transport__model__ggd__ion__state__energy__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__conductivity
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__energy__flux
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__energy
    flux :: StructArray{edge_transport__model__ggd__neutral__energy__flux} = StructArray(edge_transport__model__ggd__neutral__energy__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__neutral__energy__d} = StructArray(edge_transport__model__ggd__neutral__energy__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__neutral__energy__v} = StructArray(edge_transport__model__ggd__neutral__energy__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__neutral__energy__flux_limiter} = StructArray(edge_transport__model__ggd__neutral__energy__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__momentum__v
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__momentum__d
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__energy__flux_limiter
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__energy
    flux :: StructArray{edge_transport__model__ggd__ion__energy__flux} = StructArray(edge_transport__model__ggd__ion__energy__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__ion__energy__d} = StructArray(edge_transport__model__ggd__ion__energy__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__ion__energy__v} = StructArray(edge_transport__model__ggd__ion__energy__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__ion__energy__flux_limiter} = StructArray(edge_transport__model__ggd__ion__energy__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__neutral
    label :: String = ""
    particle_flux_integrated :: StructArray{edge_transport__model__ggd_fast__neutral__particle_flux_integrated} = StructArray(edge_transport__model__ggd_fast__neutral__particle_flux_integrated() for k in 1:1)
    ion_index :: Int32 = 0
    element :: StructArray{edge_transport__model__ggd_fast__neutral__element} = StructArray(edge_transport__model__ggd_fast__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__momentum__d
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__momentum
    flux :: StructArray{edge_transport__model__ggd__momentum__flux} = StructArray(edge_transport__model__ggd__momentum__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__momentum__d} = StructArray(edge_transport__model__ggd__momentum__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__momentum__v} = StructArray(edge_transport__model__ggd__momentum__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__momentum__flux_limiter} = StructArray(edge_transport__model__ggd__momentum__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__grid_ggd__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_transport__grid_ggd
    time :: Float64 = 0.0
    grid_subset :: StructArray{edge_transport__grid_ggd__grid_subset} = StructArray(edge_transport__grid_ggd__grid_subset() for k in 1:1)
    space :: StructArray{edge_transport__grid_ggd__space} = StructArray(edge_transport__grid_ggd__space() for k in 1:1)
    identifier :: edge_transport__grid_ggd__identifier = edge_transport__grid_ggd__identifier()
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__momentum__flux
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state__momentum
    flux :: StructArray{edge_transport__model__ggd__ion__state__momentum__flux} = StructArray(edge_transport__model__ggd__ion__state__momentum__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__ion__state__momentum__d} = StructArray(edge_transport__model__ggd__ion__state__momentum__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__ion__state__momentum__v} = StructArray(edge_transport__model__ggd__ion__state__momentum__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__ion__state__momentum__flux_limiter} = StructArray(edge_transport__model__ggd__ion__state__momentum__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    particles :: edge_transport__model__ggd__ion__state__particles = edge_transport__model__ggd__ion__state__particles()
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    momentum :: edge_transport__model__ggd__ion__state__momentum = edge_transport__model__ggd__ion__state__momentum()
    energy :: edge_transport__model__ggd__ion__state__energy = edge_transport__model__ggd__ion__state__energy()
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_transport__model__ggd__ion
    z_ion :: Float64 = 0.0
    label :: String = ""
    particles :: edge_transport__model__ggd__ion__particles = edge_transport__model__ggd__ion__particles()
    multiple_states_flag :: Int32 = 0
    momentum :: edge_transport__model__ggd__ion__momentum = edge_transport__model__ggd__ion__momentum()
    state :: StructArray{edge_transport__model__ggd__ion__state} = StructArray(edge_transport__model__ggd__ion__state() for k in 1:1)
    energy :: edge_transport__model__ggd__ion__energy = edge_transport__model__ggd__ion__energy()
    neutral_index :: Int32 = 0
    element :: StructArray{edge_transport__model__ggd__ion__element} = StructArray(edge_transport__model__ggd__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__momentum__v
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral__momentum
    flux :: StructArray{edge_transport__model__ggd__neutral__momentum__flux} = StructArray(edge_transport__model__ggd__neutral__momentum__flux() for k in 1:1)
    d :: StructArray{edge_transport__model__ggd__neutral__momentum__d} = StructArray(edge_transport__model__ggd__neutral__momentum__d() for k in 1:1)
    v :: StructArray{edge_transport__model__ggd__neutral__momentum__v} = StructArray(edge_transport__model__ggd__neutral__momentum__v() for k in 1:1)
    flux_limiter :: StructArray{edge_transport__model__ggd__neutral__momentum__flux_limiter} = StructArray(edge_transport__model__ggd__neutral__momentum__flux_limiter() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd__neutral
    label :: String = ""
    particles :: edge_transport__model__ggd__neutral__particles = edge_transport__model__ggd__neutral__particles()
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    momentum :: edge_transport__model__ggd__neutral__momentum = edge_transport__model__ggd__neutral__momentum()
    energy :: edge_transport__model__ggd__neutral__energy = edge_transport__model__ggd__neutral__energy()
    state :: StructArray{edge_transport__model__ggd__neutral__state} = StructArray(edge_transport__model__ggd__neutral__state() for k in 1:1)
    element :: StructArray{edge_transport__model__ggd__neutral__element} = StructArray(edge_transport__model__ggd__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd
    ion :: StructArray{edge_transport__model__ggd__ion} = StructArray(edge_transport__model__ggd__ion() for k in 1:1)
    time :: Float64 = 0.0
    total_ion_energy :: edge_transport__model__ggd__total_ion_energy = edge_transport__model__ggd__total_ion_energy()
    conductivity :: StructArray{edge_transport__model__ggd__conductivity} = StructArray(edge_transport__model__ggd__conductivity() for k in 1:1)
    neutral :: StructArray{edge_transport__model__ggd__neutral} = StructArray(edge_transport__model__ggd__neutral() for k in 1:1)
    momentum :: edge_transport__model__ggd__momentum = edge_transport__model__ggd__momentum()
    electrons :: edge_transport__model__ggd__electrons = edge_transport__model__ggd__electrons()
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast__electrons
    power :: StructArray{edge_transport__model__ggd_fast__electrons__power} = StructArray(edge_transport__model__ggd_fast__electrons__power() for k in 1:1)
    particle_flux_integrated :: StructArray{edge_transport__model__ggd_fast__electrons__particle_flux_integrated} = StructArray(edge_transport__model__ggd_fast__electrons__particle_flux_integrated() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport__model__ggd_fast
    ion :: StructArray{edge_transport__model__ggd_fast__ion} = StructArray(edge_transport__model__ggd_fast__ion() for k in 1:1)
    power :: StructArray{edge_transport__model__ggd_fast__power} = StructArray(edge_transport__model__ggd_fast__power() for k in 1:1)
    time :: Float64 = 0.0
    neutral :: StructArray{edge_transport__model__ggd_fast__neutral} = StructArray(edge_transport__model__ggd_fast__neutral() for k in 1:1)
    energy_flux_max :: StructArray{edge_transport__model__ggd_fast__energy_flux_max} = StructArray(edge_transport__model__ggd_fast__energy_flux_max() for k in 1:1)
    power_ion_total :: StructArray{edge_transport__model__ggd_fast__power_ion_total} = StructArray(edge_transport__model__ggd_fast__power_ion_total() for k in 1:1)
    electrons :: edge_transport__model__ggd_fast__electrons = edge_transport__model__ggd_fast__electrons()
end

Base.@kwdef mutable struct edge_transport__model
    ggd_fast :: StructArray{edge_transport__model__ggd_fast} = StructArray(edge_transport__model__ggd_fast() for k in 1:1)
    flux_multiplier :: Float64 = 0.0
    code :: edge_transport__model__code = edge_transport__model__code()
    identifier :: edge_transport__model__identifier = edge_transport__model__identifier()
    ggd :: StructArray{edge_transport__model__ggd} = StructArray(edge_transport__model__ggd() for k in 1:1)
end

Base.@kwdef mutable struct edge_transport
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: edge_transport__ids_properties = edge_transport__ids_properties()
    grid_ggd :: StructArray{edge_transport__grid_ggd} = StructArray(edge_transport__grid_ggd() for k in 1:1)
    midplane :: edge_transport__midplane = edge_transport__midplane()
    code :: edge_transport__code = edge_transport__code()
    model :: StructArray{edge_transport__model} = StructArray(edge_transport__model() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__source__ggd__electrons__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__total_ion_energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__grid_ggd__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__ggd_fast__ion__power
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__grid_ggd__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__grid_ggd__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__grid_ggd__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__species__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__grid_ggd__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__grid_ggd__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__grid_ggd__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{edge_sources__grid_ggd__space__objects_per_dimension__object__boundary} = StructArray(edge_sources__grid_ggd__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__grid_ggd__space__objects_per_dimension
    object :: StructArray{edge_sources__grid_ggd__space__objects_per_dimension__object} = StructArray(edge_sources__grid_ggd__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__grid_ggd__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{edge_sources__grid_ggd__space__objects_per_dimension} = StructArray(edge_sources__grid_ggd__space__objects_per_dimension() for k in 1:1)
    identifier :: edge_sources__grid_ggd__space__identifier = edge_sources__grid_ggd__space__identifier()
    geometry_type :: edge_sources__grid_ggd__space__geometry_type = edge_sources__grid_ggd__space__geometry_type()
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__state__momentum
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct edge_sources__grid_ggd__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct edge_sources__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__species__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__species__ion
    label :: String = ""
    element :: StructArray{edge_sources__source__species__ion__element} = StructArray(edge_sources__source__species__ion__element() for k in 1:1)
    z_ion :: Float64 = 0.0
    state :: edge_sources__source__species__ion__state = edge_sources__source__species__ion__state()
end

Base.@kwdef mutable struct edge_sources__source__ggd_fast__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__ggd_fast__ion
    power :: StructArray{edge_sources__source__ggd_fast__ion__power} = StructArray(edge_sources__source__ggd_fast__ion__power() for k in 1:1)
    z_ion :: Float64 = 0.0
    label :: String = ""
    neutral_index :: Int32 = 0
    element :: StructArray{edge_sources__source__ggd_fast__ion__element} = StructArray(edge_sources__source__ggd_fast__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__source__ggd_fast
    ion :: StructArray{edge_sources__source__ggd_fast__ion} = StructArray(edge_sources__source__ggd_fast__ion() for k in 1:1)
    time :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__momentum
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__current
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__grid_ggd__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__species__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__state__momentum
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__grid_ggd__grid_subset__element
    object :: StructArray{edge_sources__grid_ggd__grid_subset__element__object} = StructArray(edge_sources__grid_ggd__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__momentum
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__state__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__species__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__state__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__electrons__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__state__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__state
    label :: String = ""
    particles :: StructArray{edge_sources__source__ggd__ion__state__particles} = StructArray(edge_sources__source__ggd__ion__state__particles() for k in 1:1)
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    momentum :: StructArray{edge_sources__source__ggd__ion__state__momentum} = StructArray(edge_sources__source__ggd__ion__state__momentum() for k in 1:1)
    vibrational_mode :: String = ""
    energy :: StructArray{edge_sources__source__ggd__ion__state__energy} = StructArray(edge_sources__source__ggd__ion__state__energy() for k in 1:1)
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_sources__source__ggd__momentum
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_sources__grid_ggd__grid_subset
    base :: StructArray{edge_sources__grid_ggd__grid_subset__base} = StructArray(edge_sources__grid_ggd__grid_subset__base() for k in 1:1)
    metric :: edge_sources__grid_ggd__grid_subset__metric = edge_sources__grid_ggd__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: edge_sources__grid_ggd__grid_subset__identifier = edge_sources__grid_ggd__grid_subset__identifier()
    element :: StructArray{edge_sources__grid_ggd__grid_subset__element} = StructArray(edge_sources__grid_ggd__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__grid_ggd
    time :: Float64 = 0.0
    grid_subset :: StructArray{edge_sources__grid_ggd__grid_subset} = StructArray(edge_sources__grid_ggd__grid_subset() for k in 1:1)
    space :: StructArray{edge_sources__grid_ggd__space} = StructArray(edge_sources__grid_ggd__space() for k in 1:1)
    identifier :: edge_sources__grid_ggd__identifier = edge_sources__grid_ggd__identifier()
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__ggd__electrons
    particles :: StructArray{edge_sources__source__ggd__electrons__particles} = StructArray(edge_sources__source__ggd__electrons__particles() for k in 1:1)
    energy :: StructArray{edge_sources__source__ggd__electrons__energy} = StructArray(edge_sources__source__ggd__electrons__energy() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct edge_sources__code
    library :: StructArray{edge_sources__code__library} = StructArray(edge_sources__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__state__energy
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_sources__source__species__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__species__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    neutral_type :: edge_sources__source__species__neutral__state__neutral_type = edge_sources__source__species__neutral__state__neutral_type()
end

Base.@kwdef mutable struct edge_sources__source__species__neutral
    element :: StructArray{edge_sources__source__species__neutral__element} = StructArray(edge_sources__source__species__neutral__element() for k in 1:1)
    label :: String = ""
    state :: edge_sources__source__species__neutral__state = edge_sources__source__species__neutral__state()
end

Base.@kwdef mutable struct edge_sources__source__species
    ion :: edge_sources__source__species__ion = edge_sources__source__species__ion()
    type :: edge_sources__source__species__type = edge_sources__source__species__type()
    neutral :: edge_sources__source__species__neutral = edge_sources__source__species__neutral()
end

Base.@kwdef mutable struct edge_sources__ids_properties
    provider :: String = ""
    version_put :: edge_sources__ids_properties__version_put = edge_sources__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral__state
    label :: String = ""
    particles :: StructArray{edge_sources__source__ggd__neutral__state__particles} = StructArray(edge_sources__source__ggd__neutral__state__particles() for k in 1:1)
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    momentum :: StructArray{edge_sources__source__ggd__neutral__state__momentum} = StructArray(edge_sources__source__ggd__neutral__state__momentum() for k in 1:1)
    vibrational_mode :: String = ""
    energy :: StructArray{edge_sources__source__ggd__neutral__state__energy} = StructArray(edge_sources__source__ggd__neutral__state__energy() for k in 1:1)
    neutral_type :: edge_sources__source__ggd__neutral__state__neutral_type = edge_sources__source__ggd__neutral__state__neutral_type()
end

Base.@kwdef mutable struct edge_sources__source__ggd__neutral
    label :: String = ""
    particles :: StructArray{edge_sources__source__ggd__neutral__particles} = StructArray(edge_sources__source__ggd__neutral__particles() for k in 1:1)
    momentum :: StructArray{edge_sources__source__ggd__neutral__momentum} = StructArray(edge_sources__source__ggd__neutral__momentum() for k in 1:1)
    ion_index :: Int32 = 0
    energy :: StructArray{edge_sources__source__ggd__neutral__energy} = StructArray(edge_sources__source__ggd__neutral__energy() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    state :: StructArray{edge_sources__source__ggd__neutral__state} = StructArray(edge_sources__source__ggd__neutral__state() for k in 1:1)
    element :: StructArray{edge_sources__source__ggd__neutral__element} = StructArray(edge_sources__source__ggd__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__source__ggd__ion
    z_ion :: Float64 = 0.0
    particles :: StructArray{edge_sources__source__ggd__ion__particles} = StructArray(edge_sources__source__ggd__ion__particles() for k in 1:1)
    label :: String = ""
    momentum :: StructArray{edge_sources__source__ggd__ion__momentum} = StructArray(edge_sources__source__ggd__ion__momentum() for k in 1:1)
    energy :: StructArray{edge_sources__source__ggd__ion__energy} = StructArray(edge_sources__source__ggd__ion__energy() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    state :: StructArray{edge_sources__source__ggd__ion__state} = StructArray(edge_sources__source__ggd__ion__state() for k in 1:1)
    neutral_index :: Int32 = 0
    element :: StructArray{edge_sources__source__ggd__ion__element} = StructArray(edge_sources__source__ggd__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources__source__ggd
    ion :: StructArray{edge_sources__source__ggd__ion} = StructArray(edge_sources__source__ggd__ion() for k in 1:1)
    time :: Float64 = 0.0
    momentum :: StructArray{edge_sources__source__ggd__momentum} = StructArray(edge_sources__source__ggd__momentum() for k in 1:1)
    total_ion_energy :: StructArray{edge_sources__source__ggd__total_ion_energy} = StructArray(edge_sources__source__ggd__total_ion_energy() for k in 1:1)
    neutral :: StructArray{edge_sources__source__ggd__neutral} = StructArray(edge_sources__source__ggd__neutral() for k in 1:1)
    current :: StructArray{edge_sources__source__ggd__current} = StructArray(edge_sources__source__ggd__current() for k in 1:1)
    electrons :: edge_sources__source__ggd__electrons = edge_sources__source__ggd__electrons()
end

Base.@kwdef mutable struct edge_sources__source
    species :: edge_sources__source__species = edge_sources__source__species()
    ggd_fast :: StructArray{edge_sources__source__ggd_fast} = StructArray(edge_sources__source__ggd_fast() for k in 1:1)
    identifier :: edge_sources__source__identifier = edge_sources__source__identifier()
    ggd :: StructArray{edge_sources__source__ggd} = StructArray(edge_sources__source__ggd() for k in 1:1)
end

Base.@kwdef mutable struct edge_sources
    source :: StructArray{edge_sources__source} = StructArray(edge_sources__source() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: edge_sources__ids_properties = edge_sources__ids_properties()
    grid_ggd :: StructArray{edge_sources__grid_ggd} = StructArray(edge_sources__grid_ggd() for k in 1:1)
    midplane :: edge_sources__midplane = edge_sources__midplane()
    code :: edge_sources__code = edge_sources__code()
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__ion__content
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__pressure_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__pressure_fast_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__ion__density
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__electrons__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__state__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{edge_profiles__grid_ggd__space__objects_per_dimension__object__boundary} = StructArray(edge_profiles__grid_ggd__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__space__objects_per_dimension
    object :: StructArray{edge_profiles__grid_ggd__space__objects_per_dimension__object} = StructArray(edge_profiles__grid_ggd__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__neutral__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__pressure_fast_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__pressure_fast_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__ion__temperature
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__pressure_fast_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__grid_subset__element
    object :: StructArray{edge_profiles__grid_ggd__grid_subset__element__object} = StructArray(edge_profiles__grid_ggd__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__velocity
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ids_properties
    provider :: String = ""
    version_put :: edge_profiles__ids_properties__version_put = edge_profiles__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__e_field
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__phi_potential
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__zeff_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__energy_density_kinetic
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__energy_thermal
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__pressure
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__z_square_average
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__density_fast
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__density_fast
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__energy_density_kinetic
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__pressure
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_parallel_viscosity
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_pfirsch_schlueter
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__pressure
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__pressure_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd__j_total
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct edge_profiles__code
    library :: StructArray{edge_profiles__code__library} = StructArray(edge_profiles__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct edge_profiles__ggd__j_anomalous
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__state__density_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: edge_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method = edge_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__temperature
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__energy_density_kinetic
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__temperature_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: edge_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method = edge_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method()
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__density_fast
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__z_average
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__grid_subset
    base :: StructArray{edge_profiles__grid_ggd__grid_subset__base} = StructArray(edge_profiles__grid_ggd__grid_subset__base() for k in 1:1)
    metric :: edge_profiles__grid_ggd__grid_subset__metric = edge_profiles__grid_ggd__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: edge_profiles__grid_ggd__grid_subset__identifier = edge_profiles__grid_ggd__grid_subset__identifier()
    element :: StructArray{edge_profiles__grid_ggd__grid_subset__element} = StructArray(edge_profiles__grid_ggd__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__ion
    temperature :: StructArray{edge_profiles__ggd_fast__ion__temperature} = StructArray(edge_profiles__ggd_fast__ion__temperature() for k in 1:1)
    label :: String = ""
    content :: StructArray{edge_profiles__ggd_fast__ion__content} = StructArray(edge_profiles__ggd_fast__ion__content() for k in 1:1)
    z_ion :: Float64 = 0.0
    density :: StructArray{edge_profiles__ggd_fast__ion__density} = StructArray(edge_profiles__ggd_fast__ion__density() for k in 1:1)
    neutral_index :: Int32 = 0
    element :: StructArray{edge_profiles__ggd_fast__ion__element} = StructArray(edge_profiles__ggd_fast__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__t_i_average_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: edge_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method = edge_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__pressure_fast_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__neutral__state__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__neutral__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    electron_configuration :: String = ""
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_mode :: String = ""
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: edge_profiles__profiles_1d__neutral__state__velocity = edge_profiles__profiles_1d__neutral__state__velocity()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_type :: edge_profiles__profiles_1d__neutral__state__neutral_type = edge_profiles__profiles_1d__neutral__state__neutral_type()
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__pressure_fast_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_diamagnetic
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_inertial
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__electrons__density
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__velocity_diamagnetic
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__density_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: edge_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method = edge_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__velocity
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__j_heat_viscosity
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__neutral
    label :: String = ""
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    ion_index :: Int32 = 0
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{edge_profiles__profiles_1d__neutral__state} = StructArray(edge_profiles__profiles_1d__neutral__state() for k in 1:1)
    velocity :: edge_profiles__profiles_1d__neutral__velocity = edge_profiles__profiles_1d__neutral__velocity()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{edge_profiles__profiles_1d__neutral__element} = StructArray(edge_profiles__profiles_1d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__pressure
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__pressure_fast_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__velocity
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__pressure_fast_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__velocity_exb
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__energy_density_kinetic
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__n_i_total_over_n_e
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__pressure_fast_parallel
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__electrons__temperature
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    value :: Float64 = 0.0
end

Base.@kwdef mutable struct edge_profiles__ggd_fast__electrons
    temperature :: StructArray{edge_profiles__ggd_fast__electrons__temperature} = StructArray(edge_profiles__ggd_fast__electrons__temperature() for k in 1:1)
    density :: StructArray{edge_profiles__ggd_fast__electrons__density} = StructArray(edge_profiles__ggd_fast__electrons__density() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__ggd_fast
    ion :: StructArray{edge_profiles__ggd_fast__ion} = StructArray(edge_profiles__ggd_fast__ion() for k in 1:1)
    time :: Float64 = 0.0
    energy_thermal :: StructArray{edge_profiles__ggd_fast__energy_thermal} = StructArray(edge_profiles__ggd_fast__energy_thermal() for k in 1:1)
    electrons :: edge_profiles__ggd_fast__electrons = edge_profiles__ggd_fast__electrons()
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__density
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__velocity_exb
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state__ionisation_potential
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__zeff_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: edge_profiles__profiles_1d__zeff_fit__time_measurement_slice_method = edge_profiles__profiles_1d__zeff_fit__time_measurement_slice_method()
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__t_i_average
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_ion_neutral_friction
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__density_fast
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__pressure_fast_perpendicular
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__pressure_thermal
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__velocity
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__state
    label :: String = ""
    velocity_exb :: StructArray{edge_profiles__ggd__ion__state__velocity_exb} = StructArray(edge_profiles__ggd__ion__state__velocity_exb() for k in 1:1)
    vibrational_level :: Float64 = 0.0
    ionisation_potential :: StructArray{edge_profiles__ggd__ion__state__ionisation_potential} = StructArray(edge_profiles__ggd__ion__state__ionisation_potential() for k in 1:1)
    velocity :: StructArray{edge_profiles__ggd__ion__state__velocity} = StructArray(edge_profiles__ggd__ion__state__velocity() for k in 1:1)
    pressure_fast_parallel :: StructArray{edge_profiles__ggd__ion__state__pressure_fast_parallel} = StructArray(edge_profiles__ggd__ion__state__pressure_fast_parallel() for k in 1:1)
    pressure :: StructArray{edge_profiles__ggd__ion__state__pressure} = StructArray(edge_profiles__ggd__ion__state__pressure() for k in 1:1)
    z_min :: Float64 = 0.0
    energy_density_kinetic :: StructArray{edge_profiles__ggd__ion__state__energy_density_kinetic} = StructArray(edge_profiles__ggd__ion__state__energy_density_kinetic() for k in 1:1)
    z_average :: StructArray{edge_profiles__ggd__ion__state__z_average} = StructArray(edge_profiles__ggd__ion__state__z_average() for k in 1:1)
    density_fast :: StructArray{edge_profiles__ggd__ion__state__density_fast} = StructArray(edge_profiles__ggd__ion__state__density_fast() for k in 1:1)
    pressure_fast_perpendicular :: StructArray{edge_profiles__ggd__ion__state__pressure_fast_perpendicular} = StructArray(edge_profiles__ggd__ion__state__pressure_fast_perpendicular() for k in 1:1)
    temperature :: StructArray{edge_profiles__ggd__ion__state__temperature} = StructArray(edge_profiles__ggd__ion__state__temperature() for k in 1:1)
    electron_configuration :: String = ""
    vibrational_mode :: String = ""
    z_square_average :: StructArray{edge_profiles__ggd__ion__state__z_square_average} = StructArray(edge_profiles__ggd__ion__state__z_square_average() for k in 1:1)
    z_max :: Float64 = 0.0
    velocity_diamagnetic :: StructArray{edge_profiles__ggd__ion__state__velocity_diamagnetic} = StructArray(edge_profiles__ggd__ion__state__velocity_diamagnetic() for k in 1:1)
    density :: StructArray{edge_profiles__ggd__ion__state__density} = StructArray(edge_profiles__ggd__ion__state__density() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct edge_profiles__grid_ggd__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    geometry_type :: edge_profiles__grid_ggd__space__geometry_type = edge_profiles__grid_ggd__space__geometry_type()
    identifier :: edge_profiles__grid_ggd__space__identifier = edge_profiles__grid_ggd__space__identifier()
    objects_per_dimension :: StructArray{edge_profiles__grid_ggd__space__objects_per_dimension} = StructArray(edge_profiles__grid_ggd__space__objects_per_dimension() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__grid_ggd
    time :: Float64 = 0.0
    grid_subset :: StructArray{edge_profiles__grid_ggd__grid_subset} = StructArray(edge_profiles__grid_ggd__grid_subset() for k in 1:1)
    space :: StructArray{edge_profiles__grid_ggd__space} = StructArray(edge_profiles__grid_ggd__space() for k in 1:1)
    identifier :: edge_profiles__grid_ggd__identifier = edge_profiles__grid_ggd__identifier()
end

Base.@kwdef mutable struct edge_profiles__ggd__ion__pressure
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__ion
    label :: String = ""
    velocity :: StructArray{edge_profiles__ggd__ion__velocity} = StructArray(edge_profiles__ggd__ion__velocity() for k in 1:1)
    pressure_fast_parallel :: StructArray{edge_profiles__ggd__ion__pressure_fast_parallel} = StructArray(edge_profiles__ggd__ion__pressure_fast_parallel() for k in 1:1)
    pressure :: StructArray{edge_profiles__ggd__ion__pressure} = StructArray(edge_profiles__ggd__ion__pressure() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    energy_density_kinetic :: StructArray{edge_profiles__ggd__ion__energy_density_kinetic} = StructArray(edge_profiles__ggd__ion__energy_density_kinetic() for k in 1:1)
    density_fast :: StructArray{edge_profiles__ggd__ion__density_fast} = StructArray(edge_profiles__ggd__ion__density_fast() for k in 1:1)
    pressure_fast_perpendicular :: StructArray{edge_profiles__ggd__ion__pressure_fast_perpendicular} = StructArray(edge_profiles__ggd__ion__pressure_fast_perpendicular() for k in 1:1)
    temperature :: StructArray{edge_profiles__ggd__ion__temperature} = StructArray(edge_profiles__ggd__ion__temperature() for k in 1:1)
    neutral_index :: Int32 = 0
    state :: StructArray{edge_profiles__ggd__ion__state} = StructArray(edge_profiles__ggd__ion__state() for k in 1:1)
    z_ion :: Float64 = 0.0
    density :: StructArray{edge_profiles__ggd__ion__density} = StructArray(edge_profiles__ggd__ion__density() for k in 1:1)
    element :: StructArray{edge_profiles__ggd__ion__element} = StructArray(edge_profiles__ggd__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__density_fast
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons__velocity
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__electrons
    temperature :: StructArray{edge_profiles__ggd__electrons__temperature} = StructArray(edge_profiles__ggd__electrons__temperature() for k in 1:1)
    velocity :: StructArray{edge_profiles__ggd__electrons__velocity} = StructArray(edge_profiles__ggd__electrons__velocity() for k in 1:1)
    pressure_fast_parallel :: StructArray{edge_profiles__ggd__electrons__pressure_fast_parallel} = StructArray(edge_profiles__ggd__electrons__pressure_fast_parallel() for k in 1:1)
    pressure :: StructArray{edge_profiles__ggd__electrons__pressure} = StructArray(edge_profiles__ggd__electrons__pressure() for k in 1:1)
    density :: StructArray{edge_profiles__ggd__electrons__density} = StructArray(edge_profiles__ggd__electrons__density() for k in 1:1)
    density_fast :: StructArray{edge_profiles__ggd__electrons__density_fast} = StructArray(edge_profiles__ggd__electrons__density_fast() for k in 1:1)
    pressure_fast_perpendicular :: StructArray{edge_profiles__ggd__electrons__pressure_fast_perpendicular} = StructArray(edge_profiles__ggd__electrons__pressure_fast_perpendicular() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__ggd__zeff
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct edge_profiles__ggd__j_perpendicular_viscosity
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    rotation_frequency_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    electron_configuration :: String = ""
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_mode :: String = ""
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    z_average_square_1d :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: edge_profiles__profiles_1d__ion__state__velocity = edge_profiles__profiles_1d__ion__state__velocity()
    z_average :: Float64 = 0.0
    z_max :: Float64 = 0.0
    z_square_average :: Float64 = 0.0
    ionisation_potential :: Float64 = 0.0
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    z_average_1d :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: edge_profiles__profiles_1d__ion__state__density_fit = edge_profiles__profiles_1d__ion__state__density_fit()
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__ion
    label :: String = ""
    rotation_frequency_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_validity :: Int32 = 0
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    z_ion_1d :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    density_validity :: Int32 = 0
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{edge_profiles__profiles_1d__ion__state} = StructArray(edge_profiles__profiles_1d__ion__state() for k in 1:1)
    velocity :: edge_profiles__profiles_1d__ion__velocity = edge_profiles__profiles_1d__ion__velocity()
    z_ion :: Float64 = 0.0
    temperature_fit :: edge_profiles__profiles_1d__ion__temperature_fit = edge_profiles__profiles_1d__ion__temperature_fit()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: edge_profiles__profiles_1d__ion__density_fit = edge_profiles__profiles_1d__ion__density_fit()
    z_ion_square_1d :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{edge_profiles__profiles_1d__ion__element} = StructArray(edge_profiles__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__electrons__density_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: edge_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method = edge_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__electrons__temperature_fit
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    source :: Array{String, 1} = String[]
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: edge_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method = edge_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method()
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__profiles_1d__electrons
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_validity :: Int32 = 0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    density_validity :: Int32 = 0
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: edge_profiles__profiles_1d__electrons__velocity = edge_profiles__profiles_1d__electrons__velocity()
    temperature_fit :: edge_profiles__profiles_1d__electrons__temperature_fit = edge_profiles__profiles_1d__electrons__temperature_fit()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    collisionality_norm :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: edge_profiles__profiles_1d__electrons__density_fit = edge_profiles__profiles_1d__electrons__density_fit()
end

Base.@kwdef mutable struct edge_profiles__profiles_1d
    pressure_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    t_i_average_fit :: edge_profiles__profiles_1d__t_i_average_fit = edge_profiles__profiles_1d__t_i_average_fit()
    neutral :: StructArray{edge_profiles__profiles_1d__neutral} = StructArray(edge_profiles__profiles_1d__neutral() for k in 1:1)
    n_i_thermal_total :: Array{Float64, 1} = zeros(Float64,(0))
    magnetic_shear :: Array{Float64, 1} = zeros(Float64,(0))
    ion :: StructArray{edge_profiles__profiles_1d__ion} = StructArray(edge_profiles__profiles_1d__ion() for k in 1:1)
    j_total :: Array{Float64, 1} = zeros(Float64,(0))
    rotation_frequency_tor_sonic :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    j_tor :: Array{Float64, 1} = zeros(Float64,(0))
    current_parallel_inside :: Array{Float64, 1} = zeros(Float64,(0))
    j_non_inductive :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: edge_profiles__profiles_1d__electrons = edge_profiles__profiles_1d__electrons()
    q :: Array{Float64, 1} = zeros(Float64,(0))
    conductivity_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    t_i_average :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    j_ohmic :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: edge_profiles__profiles_1d__grid = edge_profiles__profiles_1d__grid()
    phi_potential :: Array{Float64, 1} = zeros(Float64,(0))
    j_bootstrap :: Array{Float64, 1} = zeros(Float64,(0))
    zeff_fit :: edge_profiles__profiles_1d__zeff_fit = edge_profiles__profiles_1d__zeff_fit()
    pressure_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    e_field :: edge_profiles__profiles_1d__e_field = edge_profiles__profiles_1d__e_field()
    zeff :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state__velocity_diamagnetic
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral__state
    label :: String = ""
    velocity_exb :: StructArray{edge_profiles__ggd__neutral__state__velocity_exb} = StructArray(edge_profiles__ggd__neutral__state__velocity_exb() for k in 1:1)
    vibrational_level :: Float64 = 0.0
    pressure_fast_parallel :: StructArray{edge_profiles__ggd__neutral__state__pressure_fast_parallel} = StructArray(edge_profiles__ggd__neutral__state__pressure_fast_parallel() for k in 1:1)
    velocity :: StructArray{edge_profiles__ggd__neutral__state__velocity} = StructArray(edge_profiles__ggd__neutral__state__velocity() for k in 1:1)
    pressure :: StructArray{edge_profiles__ggd__neutral__state__pressure} = StructArray(edge_profiles__ggd__neutral__state__pressure() for k in 1:1)
    energy_density_kinetic :: StructArray{edge_profiles__ggd__neutral__state__energy_density_kinetic} = StructArray(edge_profiles__ggd__neutral__state__energy_density_kinetic() for k in 1:1)
    density_fast :: StructArray{edge_profiles__ggd__neutral__state__density_fast} = StructArray(edge_profiles__ggd__neutral__state__density_fast() for k in 1:1)
    pressure_fast_perpendicular :: StructArray{edge_profiles__ggd__neutral__state__pressure_fast_perpendicular} = StructArray(edge_profiles__ggd__neutral__state__pressure_fast_perpendicular() for k in 1:1)
    temperature :: StructArray{edge_profiles__ggd__neutral__state__temperature} = StructArray(edge_profiles__ggd__neutral__state__temperature() for k in 1:1)
    electron_configuration :: String = ""
    vibrational_mode :: String = ""
    velocity_diamagnetic :: StructArray{edge_profiles__ggd__neutral__state__velocity_diamagnetic} = StructArray(edge_profiles__ggd__neutral__state__velocity_diamagnetic() for k in 1:1)
    density :: StructArray{edge_profiles__ggd__neutral__state__density} = StructArray(edge_profiles__ggd__neutral__state__density() for k in 1:1)
    neutral_type :: edge_profiles__ggd__neutral__state__neutral_type = edge_profiles__ggd__neutral__state__neutral_type()
end

Base.@kwdef mutable struct edge_profiles__ggd__neutral
    label :: String = ""
    velocity :: StructArray{edge_profiles__ggd__neutral__velocity} = StructArray(edge_profiles__ggd__neutral__velocity() for k in 1:1)
    pressure_fast_parallel :: StructArray{edge_profiles__ggd__neutral__pressure_fast_parallel} = StructArray(edge_profiles__ggd__neutral__pressure_fast_parallel() for k in 1:1)
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    pressure :: StructArray{edge_profiles__ggd__neutral__pressure} = StructArray(edge_profiles__ggd__neutral__pressure() for k in 1:1)
    energy_density_kinetic :: StructArray{edge_profiles__ggd__neutral__energy_density_kinetic} = StructArray(edge_profiles__ggd__neutral__energy_density_kinetic() for k in 1:1)
    density_fast :: StructArray{edge_profiles__ggd__neutral__density_fast} = StructArray(edge_profiles__ggd__neutral__density_fast() for k in 1:1)
    pressure_fast_perpendicular :: StructArray{edge_profiles__ggd__neutral__pressure_fast_perpendicular} = StructArray(edge_profiles__ggd__neutral__pressure_fast_perpendicular() for k in 1:1)
    temperature :: StructArray{edge_profiles__ggd__neutral__temperature} = StructArray(edge_profiles__ggd__neutral__temperature() for k in 1:1)
    state :: StructArray{edge_profiles__ggd__neutral__state} = StructArray(edge_profiles__ggd__neutral__state() for k in 1:1)
    density :: StructArray{edge_profiles__ggd__neutral__density} = StructArray(edge_profiles__ggd__neutral__density() for k in 1:1)
    element :: StructArray{edge_profiles__ggd__neutral__element} = StructArray(edge_profiles__ggd__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles__ggd__e_field
    diamagnetic_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    parallel_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    radial_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid_index :: Int32 = 0
    grid_subset_index :: Int32 = 0
    poloidal_coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct edge_profiles__ggd
    time :: Float64 = 0.0
    neutral :: StructArray{edge_profiles__ggd__neutral} = StructArray(edge_profiles__ggd__neutral() for k in 1:1)
    ion :: StructArray{edge_profiles__ggd__ion} = StructArray(edge_profiles__ggd__ion() for k in 1:1)
    pressure_perpendicular :: StructArray{edge_profiles__ggd__pressure_perpendicular} = StructArray(edge_profiles__ggd__pressure_perpendicular() for k in 1:1)
    j_pfirsch_schlueter :: StructArray{edge_profiles__ggd__j_pfirsch_schlueter} = StructArray(edge_profiles__ggd__j_pfirsch_schlueter() for k in 1:1)
    j_total :: StructArray{edge_profiles__ggd__j_total} = StructArray(edge_profiles__ggd__j_total() for k in 1:1)
    pressure_thermal :: StructArray{edge_profiles__ggd__pressure_thermal} = StructArray(edge_profiles__ggd__pressure_thermal() for k in 1:1)
    j_ion_neutral_friction :: StructArray{edge_profiles__ggd__j_ion_neutral_friction} = StructArray(edge_profiles__ggd__j_ion_neutral_friction() for k in 1:1)
    zeff :: StructArray{edge_profiles__ggd__zeff} = StructArray(edge_profiles__ggd__zeff() for k in 1:1)
    n_i_total_over_n_e :: StructArray{edge_profiles__ggd__n_i_total_over_n_e} = StructArray(edge_profiles__ggd__n_i_total_over_n_e() for k in 1:1)
    j_diamagnetic :: StructArray{edge_profiles__ggd__j_diamagnetic} = StructArray(edge_profiles__ggd__j_diamagnetic() for k in 1:1)
    electrons :: edge_profiles__ggd__electrons = edge_profiles__ggd__electrons()
    phi_potential :: StructArray{edge_profiles__ggd__phi_potential} = StructArray(edge_profiles__ggd__phi_potential() for k in 1:1)
    j_anomalous :: StructArray{edge_profiles__ggd__j_anomalous} = StructArray(edge_profiles__ggd__j_anomalous() for k in 1:1)
    t_i_average :: StructArray{edge_profiles__ggd__t_i_average} = StructArray(edge_profiles__ggd__t_i_average() for k in 1:1)
    j_inertial :: StructArray{edge_profiles__ggd__j_inertial} = StructArray(edge_profiles__ggd__j_inertial() for k in 1:1)
    j_heat_viscosity :: StructArray{edge_profiles__ggd__j_heat_viscosity} = StructArray(edge_profiles__ggd__j_heat_viscosity() for k in 1:1)
    pressure_parallel :: StructArray{edge_profiles__ggd__pressure_parallel} = StructArray(edge_profiles__ggd__pressure_parallel() for k in 1:1)
    j_parallel_viscosity :: StructArray{edge_profiles__ggd__j_parallel_viscosity} = StructArray(edge_profiles__ggd__j_parallel_viscosity() for k in 1:1)
    j_perpendicular_viscosity :: StructArray{edge_profiles__ggd__j_perpendicular_viscosity} = StructArray(edge_profiles__ggd__j_perpendicular_viscosity() for k in 1:1)
    j_parallel :: StructArray{edge_profiles__ggd__j_parallel} = StructArray(edge_profiles__ggd__j_parallel() for k in 1:1)
    e_field :: StructArray{edge_profiles__ggd__e_field} = StructArray(edge_profiles__ggd__e_field() for k in 1:1)
end

Base.@kwdef mutable struct edge_profiles
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: edge_profiles__ids_properties = edge_profiles__ids_properties()
    ggd_fast :: StructArray{edge_profiles__ggd_fast} = StructArray(edge_profiles__ggd_fast() for k in 1:1)
    vacuum_toroidal_field :: edge_profiles__vacuum_toroidal_field = edge_profiles__vacuum_toroidal_field()
    grid_ggd :: StructArray{edge_profiles__grid_ggd} = StructArray(edge_profiles__grid_ggd() for k in 1:1)
    midplane :: edge_profiles__midplane = edge_profiles__midplane()
    code :: edge_profiles__code = edge_profiles__code()
    profiles_1d :: StructArray{edge_profiles__profiles_1d} = StructArray(edge_profiles__profiles_1d() for k in 1:1)
    ggd :: StructArray{edge_profiles__ggd} = StructArray(edge_profiles__ggd() for k in 1:1)
end

Base.@kwdef mutable struct ece__channel__optical_depth
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__channel__beam__phase__curvature
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct ece__polarizer__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct ece__polarizer__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct ece__polarizer__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct ece__channel__position
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    theta :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ece__channel__frequency
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct ece__channel__beam__spot__size
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct ece__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct ece__psi_normalization
    psi_boundary :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ece__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct ece__polarizer__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct ece__polarizer
    polarization_angle :: Float64 = 0.0
    centre :: ece__polarizer__centre = ece__polarizer__centre()
    radius :: Float64 = 0.0
    x2_unit_vector :: ece__polarizer__x2_unit_vector = ece__polarizer__x2_unit_vector()
    x3_unit_vector :: ece__polarizer__x3_unit_vector = ece__polarizer__x3_unit_vector()
    x1_unit_vector :: ece__polarizer__x1_unit_vector = ece__polarizer__x1_unit_vector()
end

Base.@kwdef mutable struct ece__channel__beam__spot__angle
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ece__channel__beam__spot
    size :: ece__channel__beam__spot__size = ece__channel__beam__spot__size()
    angle :: ece__channel__beam__spot__angle = ece__channel__beam__spot__angle()
end

Base.@kwdef mutable struct ece__channel__harmonic
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct ece__line_of_sight
    first_point :: ece__line_of_sight__first_point = ece__line_of_sight__first_point()
    second_point :: ece__line_of_sight__second_point = ece__line_of_sight__second_point()
end

Base.@kwdef mutable struct ece__channel__beam__phase__angle
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ece__channel__t_e
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__t_e_central
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__ids_properties
    provider :: String = ""
    version_put :: ece__ids_properties__version_put = ece__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct ece__code
    library :: StructArray{ece__code__library} = StructArray(ece__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__channel__delta_position_suprathermal
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    theta :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ece__channel__beam__phase
    curvature :: ece__channel__beam__phase__curvature = ece__channel__beam__phase__curvature()
    angle :: ece__channel__beam__phase__angle = ece__channel__beam__phase__angle()
end

Base.@kwdef mutable struct ece__channel__beam
    spot :: ece__channel__beam__spot = ece__channel__beam__spot()
    phase :: ece__channel__beam__phase = ece__channel__beam__phase()
end

Base.@kwdef mutable struct ece__channel__t_e_voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ece__channel
    if_bandwidth :: Float64 = 0.0
    time :: Array{Float64, 1} = zeros(Float64,(0))
    harmonic :: ece__channel__harmonic = ece__channel__harmonic()
    name :: String = ""
    t_e :: ece__channel__t_e = ece__channel__t_e()
    position :: ece__channel__position = ece__channel__position()
    t_e_voltage :: ece__channel__t_e_voltage = ece__channel__t_e_voltage()
    beam :: ece__channel__beam = ece__channel__beam()
    optical_depth :: ece__channel__optical_depth = ece__channel__optical_depth()
    delta_position_suprathermal :: ece__channel__delta_position_suprathermal = ece__channel__delta_position_suprathermal()
    frequency :: ece__channel__frequency = ece__channel__frequency()
    identifier :: String = ""
end

Base.@kwdef mutable struct ece
    psi_normalization :: ece__psi_normalization = ece__psi_normalization()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    channel :: StructArray{ece__channel} = StructArray(ece__channel() for k in 1:1)
    ids_properties :: ece__ids_properties = ece__ids_properties()
    t_e_central :: ece__t_e_central = ece__t_e_central()
    line_of_sight :: ece__line_of_sight = ece__line_of_sight()
    latency :: Float64 = 0.0
    code :: ece__code = ece__code()
    polarizer :: StructArray{ece__polarizer} = StructArray(ece__polarizer() for k in 1:1)
end

Base.@kwdef mutable struct ec_launchers__launcher__beam__phase__curvature
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct ec_launchers__launcher__steering_angle_pol
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__launcher__steering_angle_tor
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__launcher__beam__spot__angle
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct ec_launchers__launcher__frequency
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct ec_launchers__launcher__beam__phase__angle
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__launcher__beam__phase
    curvature :: ec_launchers__launcher__beam__phase__curvature = ec_launchers__launcher__beam__phase__curvature()
    angle :: ec_launchers__launcher__beam__phase__angle = ec_launchers__launcher__beam__phase__angle()
end

Base.@kwdef mutable struct ec_launchers__launcher__power_launched
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__launcher__launching_position
    time :: Array{Float64, 1} = zeros(Float64,(0))
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r_limit_min :: Float64 = 0.0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    r_limit_max :: Float64 = 0.0
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct ec_launchers__code
    library :: StructArray{ec_launchers__code__library} = StructArray(ec_launchers__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct ec_launchers__launcher__mode
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct ec_launchers__launcher__beam__spot__size
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct ec_launchers__launcher__beam__spot
    size :: ec_launchers__launcher__beam__spot__size = ec_launchers__launcher__beam__spot__size()
    angle :: ec_launchers__launcher__beam__spot__angle = ec_launchers__launcher__beam__spot__angle()
end

Base.@kwdef mutable struct ec_launchers__launcher__beam
    spot :: ec_launchers__launcher__beam__spot = ec_launchers__launcher__beam__spot()
    phase :: ec_launchers__launcher__beam__phase = ec_launchers__launcher__beam__phase()
end

Base.@kwdef mutable struct ec_launchers__launcher
    steering_angle_pol :: ec_launchers__launcher__steering_angle_pol = ec_launchers__launcher__steering_angle_pol()
    name :: String = ""
    launching_position :: ec_launchers__launcher__launching_position = ec_launchers__launcher__launching_position()
    power_launched :: ec_launchers__launcher__power_launched = ec_launchers__launcher__power_launched()
    mode :: ec_launchers__launcher__mode = ec_launchers__launcher__mode()
    beam :: ec_launchers__launcher__beam = ec_launchers__launcher__beam()
    frequency :: ec_launchers__launcher__frequency = ec_launchers__launcher__frequency()
    identifier :: String = ""
    steering_angle_tor :: ec_launchers__launcher__steering_angle_tor = ec_launchers__launcher__steering_angle_tor()
end

Base.@kwdef mutable struct ec_launchers__ids_properties
    provider :: String = ""
    version_put :: ec_launchers__ids_properties__version_put = ec_launchers__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct ec_launchers
    launcher :: StructArray{ec_launchers__launcher} = StructArray(ec_launchers__launcher() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: ec_launchers__ids_properties = ec_launchers__ids_properties()
    latency :: Float64 = 0.0
    code :: ec_launchers__code = ec_launchers__code()
end

Base.@kwdef mutable struct divertors__divertor__target__power_incident_fraction
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_conducted
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_radiated
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_recombination_neutrals
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__wetted_area
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_incident
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__current_incident
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__wetted_area
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_convected
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_neutrals
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct divertors__divertor__target__power_convected
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_black_body
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__tile__surface_outline
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_radiated
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_currents
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_flux_peak
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_recombination_plasma
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__code
    library :: StructArray{divertors__code__library} = StructArray(divertors__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct divertors__divertor__target__tile__current_incident
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_recombination_neutrals
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct divertors__divertor__power_neutrals
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__tilt_angle_pol
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__flux_expansion
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_currents
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_conducted
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__two_point_model
    time :: Float64 = 0.0
    n_e_target :: Float64 = 0.0
    sol_heat_decay_length :: Float64 = 0.0
    t_e_target :: Float64 = 0.0
    sol_heat_spreading_length :: Float64 = 0.0
end

Base.@kwdef mutable struct divertors__ids_properties
    provider :: String = ""
    version_put :: divertors__ids_properties__version_put = divertors__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct divertors__divertor__power_recombination_plasma
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__current_incident
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__power_black_body
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__target__power_incident
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__divertor__particle_flux_recycled_total
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct divertors__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct divertors__divertor__target__tile
    name :: String = ""
    surface_outline :: divertors__divertor__target__tile__surface_outline = divertors__divertor__target__tile__surface_outline()
    surface_area :: Float64 = 0.0
    shunt_index :: Int32 = 0
    current_incident :: divertors__divertor__target__tile__current_incident = divertors__divertor__target__tile__current_incident()
    identifier :: String = ""
end

Base.@kwdef mutable struct divertors__divertor__target
    temperature_limit_max :: Float64 = 0.0
    flux_expansion :: divertors__divertor__target__flux_expansion = divertors__divertor__target__flux_expansion()
    extension_z :: Float64 = 0.0
    t_e_target_sputtering_limit_max :: Float64 = 0.0
    power_incident :: divertors__divertor__target__power_incident = divertors__divertor__target__power_incident()
    power_radiated :: divertors__divertor__target__power_radiated = divertors__divertor__target__power_radiated()
    two_point_model :: StructArray{divertors__divertor__target__two_point_model} = StructArray(divertors__divertor__target__two_point_model() for k in 1:1)
    power_conducted :: divertors__divertor__target__power_conducted = divertors__divertor__target__power_conducted()
    name :: String = ""
    wetted_area :: divertors__divertor__target__wetted_area = divertors__divertor__target__wetted_area()
    tile :: StructArray{divertors__divertor__target__tile} = StructArray(divertors__divertor__target__tile() for k in 1:1)
    power_convected :: divertors__divertor__target__power_convected = divertors__divertor__target__power_convected()
    power_black_body :: divertors__divertor__target__power_black_body = divertors__divertor__target__power_black_body()
    power_flux_peak :: divertors__divertor__target__power_flux_peak = divertors__divertor__target__power_flux_peak()
    tilt_angle_pol :: divertors__divertor__target__tilt_angle_pol = divertors__divertor__target__tilt_angle_pol()
    power_recombination_neutrals :: divertors__divertor__target__power_recombination_neutrals = divertors__divertor__target__power_recombination_neutrals()
    current_incident :: divertors__divertor__target__current_incident = divertors__divertor__target__current_incident()
    heat_flux_steady_limit_max :: Float64 = 0.0
    power_recombination_plasma :: divertors__divertor__target__power_recombination_plasma = divertors__divertor__target__power_recombination_plasma()
    power_incident_fraction :: divertors__divertor__target__power_incident_fraction = divertors__divertor__target__power_incident_fraction()
    power_currents :: divertors__divertor__target__power_currents = divertors__divertor__target__power_currents()
    power_neutrals :: divertors__divertor__target__power_neutrals = divertors__divertor__target__power_neutrals()
    identifier :: String = ""
    extension_r :: Float64 = 0.0
end

Base.@kwdef mutable struct divertors__divertor
    target :: StructArray{divertors__divertor__target} = StructArray(divertors__divertor__target() for k in 1:1)
    power_incident :: divertors__divertor__power_incident = divertors__divertor__power_incident()
    power_radiated :: divertors__divertor__power_radiated = divertors__divertor__power_radiated()
    power_conducted :: divertors__divertor__power_conducted = divertors__divertor__power_conducted()
    name :: String = ""
    particle_flux_recycled_total :: divertors__divertor__particle_flux_recycled_total = divertors__divertor__particle_flux_recycled_total()
    wetted_area :: divertors__divertor__wetted_area = divertors__divertor__wetted_area()
    power_convected :: divertors__divertor__power_convected = divertors__divertor__power_convected()
    power_black_body :: divertors__divertor__power_black_body = divertors__divertor__power_black_body()
    power_recombination_neutrals :: divertors__divertor__power_recombination_neutrals = divertors__divertor__power_recombination_neutrals()
    current_incident :: divertors__divertor__current_incident = divertors__divertor__current_incident()
    power_recombination_plasma :: divertors__divertor__power_recombination_plasma = divertors__divertor__power_recombination_plasma()
    power_currents :: divertors__divertor__power_currents = divertors__divertor__power_currents()
    power_neutrals :: divertors__divertor__power_neutrals = divertors__divertor__power_neutrals()
    identifier :: String = ""
end

Base.@kwdef mutable struct divertors
    divertor :: StructArray{divertors__divertor} = StructArray(divertors__divertor() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: divertors__ids_properties = divertors__ids_properties()
    midplane :: divertors__midplane = divertors__midplane()
    code :: divertors__code = divertors__code()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__source__identifier__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__source__identifier__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__grid_subset__element
    object :: StructArray{distributions__distribution__ggd__grid__grid_subset__element__object} = StructArray(distributions__distribution__ggd__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__source__identifier__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__markers__orbit_integrals
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    values :: Array{Complex{Float64}, 5} = zeros(Complex{Float64},(0,0,0,0,0))
    m_pol :: Array{Int32, 1} = zeros(Int32,(0))
    bounce_harmonics :: Array{Int32, 1} = zeros(Int32,(0))
    expressions :: Array{String, 1} = String[]
end

Base.@kwdef mutable struct distributions__distribution__process__nbi_energy
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__co_passing__collisions__electrons
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__grid__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    theta_straight :: Array{Float64, 1} = zeros(Float64,(0))
    theta_geometric :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    type :: distributions__distribution__profiles_2d__grid__type = distributions__distribution__profiles_2d__grid__type()
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct distributions__distribution__process__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__wave__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__thermalisation
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct distributions__distribution__species__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__collisions__electrons
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__species__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__counter_passing__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__magnetic_axis
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__co_passing__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__ggd__expansion__grid_subset
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__ggd__expansion
    grid_subset :: StructArray{distributions__distribution__ggd__expansion__grid_subset} = StructArray(distributions__distribution__ggd__expansion__grid_subset() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct distributions__distribution__process__reactant_energy
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__process
    reactant_energy :: distributions__distribution__process__reactant_energy = distributions__distribution__process__reactant_energy()
    nbi_beamlets_group :: Int32 = 0
    nbi_unit :: Int32 = 0
    type :: distributions__distribution__process__type = distributions__distribution__process__type()
    nbi_energy :: distributions__distribution__process__nbi_energy = distributions__distribution__process__nbi_energy()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__thermalisation
    particles :: Float64 = 0.0
    power :: Float64 = 0.0
    torque :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__source__identifier
    process_index :: Int32 = 0
    type :: distributions__distribution__profiles_1d__co_passing__source__identifier__type = distributions__distribution__profiles_1d__co_passing__source__identifier__type()
    wave_index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__source
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: distributions__distribution__profiles_1d__co_passing__source__identifier = distributions__distribution__profiles_1d__co_passing__source__identifier()
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__collisions__ion__state
    power_fast :: Float64 = 0.0
    label :: String = ""
    electron_configuration :: String = ""
    power_thermal :: Float64 = 0.0
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Float64 = 0.0
    torque_thermal_tor :: Float64 = 0.0
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{distributions__distribution__ggd__grid__space__objects_per_dimension__object__boundary} = StructArray(distributions__distribution__ggd__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__space__objects_per_dimension
    object :: StructArray{distributions__distribution__ggd__grid__space__objects_per_dimension__object} = StructArray(distributions__distribution__ggd__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    objects_per_dimension :: StructArray{distributions__distribution__ggd__grid__space__objects_per_dimension} = StructArray(distributions__distribution__ggd__grid__space__objects_per_dimension() for k in 1:1)
    geometry_type :: distributions__distribution__ggd__grid__space__geometry_type = distributions__distribution__ggd__grid__space__geometry_type()
    identifier :: distributions__distribution__ggd__grid__space__identifier = distributions__distribution__ggd__grid__space__identifier()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__source__identifier
    process_index :: Int32 = 0
    type :: distributions__distribution__profiles_1d__counter_passing__source__identifier__type = distributions__distribution__profiles_1d__counter_passing__source__identifier__type()
    wave_index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__source
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: distributions__distribution__profiles_1d__counter_passing__source__identifier = distributions__distribution__profiles_1d__counter_passing__source__identifier()
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__collisions__electrons
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__collisions__ion__state
    label :: String = ""
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__collisions__electrons
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__markers__orbit_integrals_instant
    values :: Array{Complex{Float64}, 3} = zeros(Complex{Float64},(0,0,0))
    time_orbit :: Array{Float64, 1} = zeros(Float64,(0))
    expressions :: Array{String, 1} = String[]
end

Base.@kwdef mutable struct distributions__distribution__species__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct distributions__code
    library :: StructArray{distributions__code__library} = StructArray(distributions__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__trapped__collisions__ion__state
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    label :: String = ""
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__markers__coordinate_identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__markers
    toroidal_mode :: Int32 = 0
    time :: Float64 = 0.0
    weights :: Array{Float64, 1} = zeros(Float64,(0))
    orbit_integrals_instant :: distributions__distribution__markers__orbit_integrals_instant = distributions__distribution__markers__orbit_integrals_instant()
    orbit_integrals :: distributions__distribution__markers__orbit_integrals = distributions__distribution__markers__orbit_integrals()
    positions :: Array{Float64, 2} = zeros(Float64,(0,0))
    coordinate_identifier :: StructArray{distributions__distribution__markers__coordinate_identifier} = StructArray(distributions__distribution__markers__coordinate_identifier() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__wave
    antenna_name :: String = ""
    index_in_antenna :: Int32 = 0
    type :: distributions__distribution__wave__type = distributions__distribution__wave__type()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__collisions__electrons
    power_fast :: Float64 = 0.0
    power_thermal :: Float64 = 0.0
    torque_fast_tor :: Float64 = 0.0
    torque_thermal_tor :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__counter_passing__collisions__electrons
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__collisions__ion__state
    label :: String = ""
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__collisions__ion
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    label :: String = ""
    z_ion :: Float64 = 0.0
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{distributions__distribution__profiles_1d__co_passing__collisions__ion__state} = StructArray(distributions__distribution__profiles_1d__co_passing__collisions__ion__state() for k in 1:1)
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    element :: StructArray{distributions__distribution__profiles_1d__co_passing__collisions__ion__element} = StructArray(distributions__distribution__profiles_1d__co_passing__collisions__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing__collisions
    ion :: StructArray{distributions__distribution__profiles_1d__co_passing__collisions__ion} = StructArray(distributions__distribution__profiles_1d__co_passing__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_1d__co_passing__collisions__electrons = distributions__distribution__profiles_1d__co_passing__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__co_passing
    source :: StructArray{distributions__distribution__profiles_1d__co_passing__source} = StructArray(distributions__distribution__profiles_1d__co_passing__source() for k in 1:1)
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    collisions :: distributions__distribution__profiles_1d__co_passing__collisions = distributions__distribution__profiles_1d__co_passing__collisions()
    current_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    current_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_tor_j_radial :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__collisions__electrons
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid__grid_subset
    base :: StructArray{distributions__distribution__ggd__grid__grid_subset__base} = StructArray(distributions__distribution__ggd__grid__grid_subset__base() for k in 1:1)
    metric :: distributions__distribution__ggd__grid__grid_subset__metric = distributions__distribution__ggd__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: distributions__distribution__ggd__grid__grid_subset__identifier = distributions__distribution__ggd__grid__grid_subset__identifier()
    element :: StructArray{distributions__distribution__ggd__grid__grid_subset__element} = StructArray(distributions__distribution__ggd__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__ggd__grid
    grid_subset :: StructArray{distributions__distribution__ggd__grid__grid_subset} = StructArray(distributions__distribution__ggd__grid__grid_subset() for k in 1:1)
    space :: StructArray{distributions__distribution__ggd__grid__space} = StructArray(distributions__distribution__ggd__grid__space() for k in 1:1)
    identifier :: distributions__distribution__ggd__grid__identifier = distributions__distribution__ggd__grid__identifier()
end

Base.@kwdef mutable struct distributions__distribution__ggd
    time :: Float64 = 0.0
    grid :: distributions__distribution__ggd__grid = distributions__distribution__ggd__grid()
    expansion :: StructArray{distributions__distribution__ggd__expansion} = StructArray(distributions__distribution__ggd__expansion() for k in 1:1)
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__species__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    neutral_type :: distributions__distribution__species__neutral__state__neutral_type = distributions__distribution__species__neutral__state__neutral_type()
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__source__identifier__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__collisions__ion__state
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    label :: String = ""
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__source__identifier
    process_index :: Int32 = 0
    type :: distributions__distribution__global_quantities__source__identifier__type = distributions__distribution__global_quantities__source__identifier__type()
    wave_index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__source
    particles :: Float64 = 0.0
    torque_tor :: Float64 = 0.0
    power :: Float64 = 0.0
    identifier :: distributions__distribution__global_quantities__source__identifier = distributions__distribution__global_quantities__source__identifier()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__counter_passing__collisions__ion__state
    label :: String = ""
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__ids_properties
    provider :: String = ""
    version_put :: distributions__ids_properties__version_put = distributions__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__trapped__collisions__electrons
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__collisions__ion__state
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    label :: String = ""
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__collisions__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    state :: StructArray{distributions__distribution__profiles_1d__trapped__collisions__ion__state} = StructArray(distributions__distribution__profiles_1d__trapped__collisions__ion__state() for k in 1:1)
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    element :: StructArray{distributions__distribution__profiles_1d__trapped__collisions__ion__element} = StructArray(distributions__distribution__profiles_1d__trapped__collisions__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__collisions
    ion :: StructArray{distributions__distribution__profiles_1d__trapped__collisions__ion} = StructArray(distributions__distribution__profiles_1d__trapped__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_1d__trapped__collisions__electrons = distributions__distribution__profiles_1d__trapped__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__collisions__ion
    label :: String = ""
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{distributions__distribution__profiles_1d__collisions__ion__element} = StructArray(distributions__distribution__profiles_1d__collisions__ion__element() for k in 1:1)
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    state :: StructArray{distributions__distribution__profiles_1d__collisions__ion__state} = StructArray(distributions__distribution__profiles_1d__collisions__ion__state() for k in 1:1)
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__collisions
    ion :: StructArray{distributions__distribution__profiles_1d__collisions__ion} = StructArray(distributions__distribution__profiles_1d__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_1d__collisions__electrons = distributions__distribution__profiles_1d__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__co_passing__collisions__ion__state
    label :: String = ""
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__co_passing__collisions__ion
    label :: String = ""
    z_ion :: Float64 = 0.0
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    multiple_states_flag :: Int32 = 0
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    state :: StructArray{distributions__distribution__profiles_2d__co_passing__collisions__ion__state} = StructArray(distributions__distribution__profiles_2d__co_passing__collisions__ion__state() for k in 1:1)
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    neutral_index :: Int32 = 0
    element :: StructArray{distributions__distribution__profiles_2d__co_passing__collisions__ion__element} = StructArray(distributions__distribution__profiles_2d__co_passing__collisions__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__co_passing__collisions
    ion :: StructArray{distributions__distribution__profiles_2d__co_passing__collisions__ion} = StructArray(distributions__distribution__profiles_2d__co_passing__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_2d__co_passing__collisions__electrons = distributions__distribution__profiles_2d__co_passing__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__collisions__electrons
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__species__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__species__ion
    label :: String = ""
    element :: StructArray{distributions__distribution__species__ion__element} = StructArray(distributions__distribution__species__ion__element() for k in 1:1)
    z_ion :: Float64 = 0.0
    state :: distributions__distribution__species__ion__state = distributions__distribution__species__ion__state()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__source__identifier
    process_index :: Int32 = 0
    type :: distributions__distribution__profiles_1d__trapped__source__identifier__type = distributions__distribution__profiles_1d__trapped__source__identifier__type()
    wave_index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped__source
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: distributions__distribution__profiles_1d__trapped__source__identifier = distributions__distribution__profiles_1d__trapped__source__identifier()
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__trapped
    source :: StructArray{distributions__distribution__profiles_1d__trapped__source} = StructArray(distributions__distribution__profiles_1d__trapped__source() for k in 1:1)
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    collisions :: distributions__distribution__profiles_1d__trapped__collisions = distributions__distribution__profiles_1d__trapped__collisions()
    current_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    torque_tor_j_radial :: Array{Float64, 1} = zeros(Float64,(0))
    current_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__collisions__ion
    z_ion :: Float64 = 0.0
    label :: String = ""
    power_fast :: Array{Float64, 1} = zeros(Float64,(0))
    power_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    torque_thermal_tor :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{distributions__distribution__profiles_1d__counter_passing__collisions__ion__state} = StructArray(distributions__distribution__profiles_1d__counter_passing__collisions__ion__state() for k in 1:1)
    torque_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    element :: StructArray{distributions__distribution__profiles_1d__counter_passing__collisions__ion__element} = StructArray(distributions__distribution__profiles_1d__counter_passing__collisions__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing__collisions
    ion :: StructArray{distributions__distribution__profiles_1d__counter_passing__collisions__ion} = StructArray(distributions__distribution__profiles_1d__counter_passing__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_1d__counter_passing__collisions__electrons = distributions__distribution__profiles_1d__counter_passing__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__counter_passing
    source :: StructArray{distributions__distribution__profiles_1d__counter_passing__source} = StructArray(distributions__distribution__profiles_1d__counter_passing__source() for k in 1:1)
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density :: Array{Float64, 1} = zeros(Float64,(0))
    collisions :: distributions__distribution__profiles_1d__counter_passing__collisions = distributions__distribution__profiles_1d__counter_passing__collisions()
    current_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    torque_tor_j_radial :: Array{Float64, 1} = zeros(Float64,(0))
    current_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__collisions__ion
    z_ion :: Float64 = 0.0
    label :: String = ""
    power_fast :: Float64 = 0.0
    power_thermal :: Float64 = 0.0
    torque_fast_tor :: Float64 = 0.0
    torque_thermal_tor :: Float64 = 0.0
    state :: StructArray{distributions__distribution__global_quantities__collisions__ion__state} = StructArray(distributions__distribution__global_quantities__collisions__ion__state() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    neutral_index :: Int32 = 0
    element :: StructArray{distributions__distribution__global_quantities__collisions__ion__element} = StructArray(distributions__distribution__global_quantities__collisions__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__global_quantities__collisions
    ion :: StructArray{distributions__distribution__global_quantities__collisions__ion} = StructArray(distributions__distribution__global_quantities__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__global_quantities__collisions__electrons = distributions__distribution__global_quantities__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__global_quantities
    time :: Float64 = 0.0
    particles_fast_n :: Float64 = 0.0
    energy_fast_parallel :: Float64 = 0.0
    energy :: Float64 = 0.0
    particles_n :: Float64 = 0.0
    torque_tor_j_radial :: Float64 = 0.0
    energy_fast :: Float64 = 0.0
    current_tor :: Float64 = 0.0
    source :: StructArray{distributions__distribution__global_quantities__source} = StructArray(distributions__distribution__global_quantities__source() for k in 1:1)
    thermalisation :: distributions__distribution__global_quantities__thermalisation = distributions__distribution__global_quantities__thermalisation()
    collisions :: distributions__distribution__global_quantities__collisions = distributions__distribution__global_quantities__collisions()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__counter_passing__collisions__ion
    element :: StructArray{distributions__distribution__profiles_2d__counter_passing__collisions__ion__element} = StructArray(distributions__distribution__profiles_2d__counter_passing__collisions__ion__element() for k in 1:1)
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    label :: String = ""
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    state :: StructArray{distributions__distribution__profiles_2d__counter_passing__collisions__ion__state} = StructArray(distributions__distribution__profiles_2d__counter_passing__collisions__ion__state() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    neutral_index :: Int32 = 0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__counter_passing__collisions
    ion :: StructArray{distributions__distribution__profiles_2d__counter_passing__collisions__ion} = StructArray(distributions__distribution__profiles_2d__counter_passing__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_2d__counter_passing__collisions__electrons = distributions__distribution__profiles_2d__counter_passing__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__counter_passing
    pressure :: Array{Float64, 2} = zeros(Float64,(0,0))
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    collisions :: distributions__distribution__profiles_2d__counter_passing__collisions = distributions__distribution__profiles_2d__counter_passing__collisions()
    current_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_tor_j_radial :: Array{Float64, 2} = zeros(Float64,(0,0))
    density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast_parallel :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__trapped__collisions__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__trapped__collisions__ion
    label :: String = ""
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_ion :: Float64 = 0.0
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    state :: StructArray{distributions__distribution__profiles_2d__trapped__collisions__ion__state} = StructArray(distributions__distribution__profiles_2d__trapped__collisions__ion__state() for k in 1:1)
    multiple_states_flag :: Int32 = 0
    neutral_index :: Int32 = 0
    element :: StructArray{distributions__distribution__profiles_2d__trapped__collisions__ion__element} = StructArray(distributions__distribution__profiles_2d__trapped__collisions__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__trapped__collisions
    ion :: StructArray{distributions__distribution__profiles_2d__trapped__collisions__ion} = StructArray(distributions__distribution__profiles_2d__trapped__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_2d__trapped__collisions__electrons = distributions__distribution__profiles_2d__trapped__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__trapped
    pressure :: Array{Float64, 2} = zeros(Float64,(0,0))
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    collisions :: distributions__distribution__profiles_2d__trapped__collisions = distributions__distribution__profiles_2d__trapped__collisions()
    current_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_tor_j_radial :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast_parallel :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__collisions__ion__state
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    label :: String = ""
    electron_configuration :: String = ""
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__collisions__ion
    element :: StructArray{distributions__distribution__profiles_2d__collisions__ion__element} = StructArray(distributions__distribution__profiles_2d__collisions__ion__element() for k in 1:1)
    power_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    label :: String = ""
    power_thermal :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    multiple_states_flag :: Int32 = 0
    state :: StructArray{distributions__distribution__profiles_2d__collisions__ion__state} = StructArray(distributions__distribution__profiles_2d__collisions__ion__state() for k in 1:1)
    torque_thermal_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    neutral_index :: Int32 = 0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__collisions
    ion :: StructArray{distributions__distribution__profiles_2d__collisions__ion} = StructArray(distributions__distribution__profiles_2d__collisions__ion() for k in 1:1)
    electrons :: distributions__distribution__profiles_2d__collisions__electrons = distributions__distribution__profiles_2d__collisions__electrons()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__fast_filter__method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__fast_filter
    method :: distributions__distribution__profiles_1d__fast_filter__method = distributions__distribution__profiles_1d__fast_filter__method()
    energy :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__species__neutral
    label :: String = ""
    element :: StructArray{distributions__distribution__species__neutral__element} = StructArray(distributions__distribution__species__neutral__element() for k in 1:1)
    state :: distributions__distribution__species__neutral__state = distributions__distribution__species__neutral__state()
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d__co_passing
    pressure :: Array{Float64, 2} = zeros(Float64,(0,0))
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    collisions :: distributions__distribution__profiles_2d__co_passing__collisions = distributions__distribution__profiles_2d__co_passing__collisions()
    current_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_tor_j_radial :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast_parallel :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_2d
    co_passing :: distributions__distribution__profiles_2d__co_passing = distributions__distribution__profiles_2d__co_passing()
    time :: Float64 = 0.0
    trapped :: distributions__distribution__profiles_2d__trapped = distributions__distribution__profiles_2d__trapped()
    pressure_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_fast_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    torque_tor_j_radial :: Array{Float64, 2} = zeros(Float64,(0,0))
    current_tor :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure :: Array{Float64, 2} = zeros(Float64,(0,0))
    pressure_fast_parallel :: Array{Float64, 2} = zeros(Float64,(0,0))
    grid :: distributions__distribution__profiles_2d__grid = distributions__distribution__profiles_2d__grid()
    counter_passing :: distributions__distribution__profiles_2d__counter_passing = distributions__distribution__profiles_2d__counter_passing()
    collisions :: distributions__distribution__profiles_2d__collisions = distributions__distribution__profiles_2d__collisions()
    density :: Array{Float64, 2} = zeros(Float64,(0,0))
    density_fast :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__source__identifier__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__source__identifier
    process_index :: Int32 = 0
    type :: distributions__distribution__profiles_1d__source__identifier__type = distributions__distribution__profiles_1d__source__identifier__type()
    wave_index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__source
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    identifier :: distributions__distribution__profiles_1d__source__identifier = distributions__distribution__profiles_1d__source__identifier()
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution__species__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distributions__distribution__species
    ion :: distributions__distribution__species__ion = distributions__distribution__species__ion()
    type :: distributions__distribution__species__type = distributions__distribution__species__type()
    neutral :: distributions__distribution__species__neutral = distributions__distribution__species__neutral()
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct distributions__distribution__profiles_1d
    co_passing :: distributions__distribution__profiles_1d__co_passing = distributions__distribution__profiles_1d__co_passing()
    time :: Float64 = 0.0
    trapped :: distributions__distribution__profiles_1d__trapped = distributions__distribution__profiles_1d__trapped()
    pressure_fast :: Array{Float64, 1} = zeros(Float64,(0))
    fast_filter :: distributions__distribution__profiles_1d__fast_filter = distributions__distribution__profiles_1d__fast_filter()
    current_fast_tor :: Array{Float64, 1} = zeros(Float64,(0))
    torque_tor_j_radial :: Array{Float64, 1} = zeros(Float64,(0))
    current_tor :: Array{Float64, 1} = zeros(Float64,(0))
    source :: StructArray{distributions__distribution__profiles_1d__source} = StructArray(distributions__distribution__profiles_1d__source() for k in 1:1)
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: distributions__distribution__profiles_1d__grid = distributions__distribution__profiles_1d__grid()
    thermalisation :: distributions__distribution__profiles_1d__thermalisation = distributions__distribution__profiles_1d__thermalisation()
    counter_passing :: distributions__distribution__profiles_1d__counter_passing = distributions__distribution__profiles_1d__counter_passing()
    collisions :: distributions__distribution__profiles_1d__collisions = distributions__distribution__profiles_1d__collisions()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distributions__distribution
    ggd :: StructArray{distributions__distribution__ggd} = StructArray(distributions__distribution__ggd() for k in 1:1)
    gyro_type :: Int32 = 0
    wave :: StructArray{distributions__distribution__wave} = StructArray(distributions__distribution__wave() for k in 1:1)
    is_delta_f :: Int32 = 0
    process :: StructArray{distributions__distribution__process} = StructArray(distributions__distribution__process() for k in 1:1)
    global_quantities :: StructArray{distributions__distribution__global_quantities} = StructArray(distributions__distribution__global_quantities() for k in 1:1)
    profiles_1d :: StructArray{distributions__distribution__profiles_1d} = StructArray(distributions__distribution__profiles_1d() for k in 1:1)
    profiles_2d :: StructArray{distributions__distribution__profiles_2d} = StructArray(distributions__distribution__profiles_2d() for k in 1:1)
    markers :: StructArray{distributions__distribution__markers} = StructArray(distributions__distribution__markers() for k in 1:1)
    species :: distributions__distribution__species = distributions__distribution__species()
end

Base.@kwdef mutable struct distributions
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: distributions__ids_properties = distributions__ids_properties()
    vacuum_toroidal_field :: distributions__vacuum_toroidal_field = distributions__vacuum_toroidal_field()
    distribution :: StructArray{distributions__distribution} = StructArray(distributions__distribution() for k in 1:1)
    code :: distributions__code = distributions__code()
    magnetic_axis :: distributions__magnetic_axis = distributions__magnetic_axis()
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__grid_subset__metric
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct distribution_sources__source__species__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__source__process__nbi_energy
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__space__objects_per_dimension__object__boundary
    neighbours :: Array{Int32, 1} = zeros(Int32,(0))
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__markers__coordinate_identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__grid_subset__base
    jacobian :: Array{Float64, 1} = zeros(Float64,(0))
    tensor_contravariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    tensor_covariant :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct distribution_sources__source__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__source__profiles_1d
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    grid :: distribution_sources__source__profiles_1d__grid = distribution_sources__source__profiles_1d__grid()
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distribution_sources__source__global_quantities__shinethrough
    particles :: Float64 = 0.0
    torque_tor :: Float64 = 0.0
    power :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__source__global_quantities
    particles :: Float64 = 0.0
    time :: Float64 = 0.0
    shinethrough :: distribution_sources__source__global_quantities__shinethrough = distribution_sources__source__global_quantities__shinethrough()
    torque_tor :: Float64 = 0.0
    power :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__source__markers__orbit_integrals
    n_tor :: Array{Int32, 1} = zeros(Int32,(0))
    values :: Array{Complex{Float64}, 5} = zeros(Complex{Float64},(0,0,0,0,0))
    m_pol :: Array{Int32, 1} = zeros(Int32,(0))
    bounce_harmonics :: Array{Int32, 1} = zeros(Int32,(0))
    expressions :: Array{String, 1} = String[]
end

Base.@kwdef mutable struct distribution_sources__source__species__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__magnetic_axis
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct distribution_sources__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct distribution_sources__code
    library :: StructArray{distribution_sources__code__library} = StructArray(distribution_sources__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct distribution_sources__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__source__species__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct distribution_sources__source__species__ion
    state :: distribution_sources__source__species__ion__state = distribution_sources__source__species__ion__state()
    z_ion :: Float64 = 0.0
    label :: String = ""
    element :: StructArray{distribution_sources__source__species__ion__element} = StructArray(distribution_sources__source__species__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__ids_properties
    provider :: String = ""
    version_put :: distribution_sources__ids_properties__version_put = distribution_sources__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__markers__orbit_integrals_instant
    values :: Array{Complex{Float64}, 3} = zeros(Complex{Float64},(0,0,0))
    time_orbit :: Array{Float64, 1} = zeros(Float64,(0))
    expressions :: Array{String, 1} = String[]
end

Base.@kwdef mutable struct distribution_sources__source__markers
    toroidal_mode :: Int32 = 0
    time :: Float64 = 0.0
    weights :: Array{Float64, 1} = zeros(Float64,(0))
    orbit_integrals_instant :: distribution_sources__source__markers__orbit_integrals_instant = distribution_sources__source__markers__orbit_integrals_instant()
    orbit_integrals :: distribution_sources__source__markers__orbit_integrals = distribution_sources__source__markers__orbit_integrals()
    positions :: Array{Float64, 2} = zeros(Float64,(0,0))
    coordinate_identifier :: StructArray{distribution_sources__source__markers__coordinate_identifier} = StructArray(distribution_sources__source__markers__coordinate_identifier() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__source__species__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__process__reactant_energy
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__space__objects_per_dimension__object
    nodes :: Array{Int32, 1} = zeros(Int32,(0))
    measure :: Float64 = 0.0
    geometry :: Array{Float64, 1} = zeros(Float64,(0))
    boundary :: StructArray{distribution_sources__source__ggd__grid__space__objects_per_dimension__object__boundary} = StructArray(distribution_sources__source__ggd__grid__space__objects_per_dimension__object__boundary() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__space__objects_per_dimension
    object :: StructArray{distribution_sources__source__ggd__grid__space__objects_per_dimension__object} = StructArray(distribution_sources__source__ggd__grid__space__objects_per_dimension__object() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__space__geometry_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__particles
    grid_index :: Int32 = 0
    values :: Array{Float64, 1} = zeros(Float64,(0))
    grid_subset_index :: Int32 = 0
    coefficients :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__space__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__space
    coordinates_type :: Array{Int32, 1} = zeros(Int32,(0))
    geometry_type :: distribution_sources__source__ggd__grid__space__geometry_type = distribution_sources__source__ggd__grid__space__geometry_type()
    identifier :: distribution_sources__source__ggd__grid__space__identifier = distribution_sources__source__ggd__grid__space__identifier()
    objects_per_dimension :: StructArray{distribution_sources__source__ggd__grid__space__objects_per_dimension} = StructArray(distribution_sources__source__ggd__grid__space__objects_per_dimension() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__source__process__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__process
    reactant_energy :: distribution_sources__source__process__reactant_energy = distribution_sources__source__process__reactant_energy()
    nbi_beamlets_group :: Int32 = 0
    nbi_unit :: Int32 = 0
    type :: distribution_sources__source__process__type = distribution_sources__source__process__type()
    nbi_energy :: distribution_sources__source__process__nbi_energy = distribution_sources__source__process__nbi_energy()
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__grid_subset__element__object
    dimension :: Int32 = 0
    space :: Int32 = 0
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__grid_subset__element
    object :: StructArray{distribution_sources__source__ggd__grid__grid_subset__element__object} = StructArray(distribution_sources__source__ggd__grid__grid_subset__element__object() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__grid_subset__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid__grid_subset
    base :: StructArray{distribution_sources__source__ggd__grid__grid_subset__base} = StructArray(distribution_sources__source__ggd__grid__grid_subset__base() for k in 1:1)
    metric :: distribution_sources__source__ggd__grid__grid_subset__metric = distribution_sources__source__ggd__grid__grid_subset__metric()
    dimension :: Int32 = 0
    identifier :: distribution_sources__source__ggd__grid__grid_subset__identifier = distribution_sources__source__ggd__grid__grid_subset__identifier()
    element :: StructArray{distribution_sources__source__ggd__grid__grid_subset__element} = StructArray(distribution_sources__source__ggd__grid__grid_subset__element() for k in 1:1)
end

Base.@kwdef mutable struct distribution_sources__source__ggd__grid
    grid_subset :: StructArray{distribution_sources__source__ggd__grid__grid_subset} = StructArray(distribution_sources__source__ggd__grid__grid_subset() for k in 1:1)
    space :: StructArray{distribution_sources__source__ggd__grid__space} = StructArray(distribution_sources__source__ggd__grid__space() for k in 1:1)
    identifier :: distribution_sources__source__ggd__grid__identifier = distribution_sources__source__ggd__grid__identifier()
end

Base.@kwdef mutable struct distribution_sources__source__ggd
    particles :: StructArray{distribution_sources__source__ggd__particles} = StructArray(distribution_sources__source__ggd__particles() for k in 1:1)
    time :: Float64 = 0.0
    grid :: distribution_sources__source__ggd__grid = distribution_sources__source__ggd__grid()
    discrete :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct distribution_sources__source__species__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct distribution_sources__source__species__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    neutral_type :: distribution_sources__source__species__neutral__state__neutral_type = distribution_sources__source__species__neutral__state__neutral_type()
end

Base.@kwdef mutable struct distribution_sources__source__species__neutral
    element :: StructArray{distribution_sources__source__species__neutral__element} = StructArray(distribution_sources__source__species__neutral__element() for k in 1:1)
    label :: String = ""
    state :: distribution_sources__source__species__neutral__state = distribution_sources__source__species__neutral__state()
end

Base.@kwdef mutable struct distribution_sources__source__species
    ion :: distribution_sources__source__species__ion = distribution_sources__source__species__ion()
    type :: distribution_sources__source__species__type = distribution_sources__source__species__type()
    neutral :: distribution_sources__source__species__neutral = distribution_sources__source__species__neutral()
end

Base.@kwdef mutable struct distribution_sources__source
    ggd :: StructArray{distribution_sources__source__ggd} = StructArray(distribution_sources__source__ggd() for k in 1:1)
    gyro_type :: Int32 = 0
    process :: StructArray{distribution_sources__source__process} = StructArray(distribution_sources__source__process() for k in 1:1)
    global_quantities :: StructArray{distribution_sources__source__global_quantities} = StructArray(distribution_sources__source__global_quantities() for k in 1:1)
    profiles_1d :: StructArray{distribution_sources__source__profiles_1d} = StructArray(distribution_sources__source__profiles_1d() for k in 1:1)
    markers :: StructArray{distribution_sources__source__markers} = StructArray(distribution_sources__source__markers() for k in 1:1)
    species :: distribution_sources__source__species = distribution_sources__source__species()
end

Base.@kwdef mutable struct distribution_sources
    source :: StructArray{distribution_sources__source} = StructArray(distribution_sources__source() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: distribution_sources__ids_properties = distribution_sources__ids_properties()
    vacuum_toroidal_field :: distribution_sources__vacuum_toroidal_field = distribution_sources__vacuum_toroidal_field()
    code :: distribution_sources__code = distribution_sources__code()
    magnetic_axis :: distribution_sources__magnetic_axis = distribution_sources__magnetic_axis()
end

Base.@kwdef mutable struct disruption__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct disruption__ids_properties
    provider :: String = ""
    version_put :: disruption__ids_properties__version_put = disruption__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct disruption__halo_currents__area__end_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct disruption__halo_currents__area__start_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct disruption__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct disruption__code
    library :: StructArray{disruption__code__library} = StructArray(disruption__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct disruption__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct disruption__halo_currents__active_wall_point
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct disruption__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct disruption__halo_currents__area
    start_point :: disruption__halo_currents__area__start_point = disruption__halo_currents__area__start_point()
    current_halo_pol :: Float64 = 0.0
    end_point :: disruption__halo_currents__area__end_point = disruption__halo_currents__area__end_point()
end

Base.@kwdef mutable struct disruption__halo_currents
    time :: Float64 = 0.0
    active_wall_point :: disruption__halo_currents__active_wall_point = disruption__halo_currents__active_wall_point()
    area :: StructArray{disruption__halo_currents__area} = StructArray(disruption__halo_currents__area() for k in 1:1)
end

Base.@kwdef mutable struct disruption__profiles_1d
    power_density_radiative_losses :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    power_density_conductive_losses :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: disruption__profiles_1d__grid = disruption__profiles_1d__grid()
    j_runaways :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct disruption__global_quantities
    energy_radiated_electrons_impurities_halo :: Array{Float64, 1} = zeros(Float64,(0))
    psi_halo_boundary :: Array{Float64, 1} = zeros(Float64,(0))
    power_ohm :: Array{Float64, 1} = zeros(Float64,(0))
    current_halo_tor :: Array{Float64, 1} = zeros(Float64,(0))
    energy_ohm_halo :: Array{Float64, 1} = zeros(Float64,(0))
    power_ohm_halo :: Array{Float64, 1} = zeros(Float64,(0))
    energy_parallel_halo :: Array{Float64, 1} = zeros(Float64,(0))
    power_radiated_electrons_impurities_halo :: Array{Float64, 1} = zeros(Float64,(0))
    current_halo_pol :: Array{Float64, 1} = zeros(Float64,(0))
    energy_radiated_electrons_impurities :: Array{Float64, 1} = zeros(Float64,(0))
    energy_ohm :: Array{Float64, 1} = zeros(Float64,(0))
    power_parallel_halo :: Array{Float64, 1} = zeros(Float64,(0))
    power_radiated_electrons_impurities :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct disruption
    halo_currents :: StructArray{disruption__halo_currents} = StructArray(disruption__halo_currents() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    vacuum_toroidal_field :: disruption__vacuum_toroidal_field = disruption__vacuum_toroidal_field()
    ids_properties :: disruption__ids_properties = disruption__ids_properties()
    code :: disruption__code = disruption__code()
    global_quantities :: disruption__global_quantities = disruption__global_quantities()
    profiles_1d :: StructArray{disruption__profiles_1d} = StructArray(disruption__profiles_1d() for k in 1:1)
end

Base.@kwdef mutable struct dataset_fair__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct dataset_fair__ids_properties
    provider :: String = ""
    version_put :: dataset_fair__ids_properties__version_put = dataset_fair__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct dataset_fair
    rights_holder :: String = ""
    is_referenced_by :: Array{String, 1} = String[]
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: dataset_fair__ids_properties = dataset_fair__ids_properties()
    valid :: String = ""
    identifier :: String = ""
    replaces :: String = ""
    license :: String = ""
    is_replaced_by :: String = ""
end

Base.@kwdef mutable struct dataset_description__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct dataset_description__pulse_time_end_epoch
    nanoseconds :: Int32 = 0
    seconds :: Int32 = 0
end

Base.@kwdef mutable struct dataset_description__data_entry
    pulse_type :: String = ""
    run :: Int32 = 0
    machine :: String = ""
    pulse :: Int32 = 0
    user :: String = ""
end

Base.@kwdef mutable struct dataset_description__parent_entry
    pulse_type :: String = ""
    run :: Int32 = 0
    machine :: String = ""
    pulse :: Int32 = 0
    user :: String = ""
end

Base.@kwdef mutable struct dataset_description__ids_properties
    provider :: String = ""
    version_put :: dataset_description__ids_properties__version_put = dataset_description__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct dataset_description__pulse_time_begin_epoch
    nanoseconds :: Int32 = 0
    seconds :: Int32 = 0
end

Base.@kwdef mutable struct dataset_description__simulation
    time_ended :: String = ""
    time_begun :: String = ""
    time_current :: Float64 = 0.0
    time_restart :: Float64 = 0.0
    workflow :: String = ""
    comment_after :: String = ""
    time_step :: Float64 = 0.0
    time_begin :: Float64 = 0.0
    comment_before :: String = ""
    time_end :: Float64 = 0.0
end

Base.@kwdef mutable struct dataset_description
    pulse_time_begin_epoch :: dataset_description__pulse_time_begin_epoch = dataset_description__pulse_time_begin_epoch()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: dataset_description__ids_properties = dataset_description__ids_properties()
    imas_version :: String = ""
    parent_entry :: dataset_description__parent_entry = dataset_description__parent_entry()
    dd_version :: String = ""
    simulation :: dataset_description__simulation = dataset_description__simulation()
    pulse_time_end_epoch :: dataset_description__pulse_time_end_epoch = dataset_description__pulse_time_end_epoch()
    data_entry :: dataset_description__data_entry = dataset_description__data_entry()
    pulse_time_begin :: String = ""
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__annular__centreline
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__element__j_tor
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__element__outline
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__annular__outline_outer
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__annular__outline_outer
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__element__j_tor
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__element__outline
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__ids_properties
    provider :: String = ""
    version_put :: cryostat__ids_properties__version_put = cryostat__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct cryostat__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__annular__outline_inner
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__element
    name :: String = ""
    j_tor :: cryostat__description_2d__thermal_shield__unit__element__j_tor = cryostat__description_2d__thermal_shield__unit__element__j_tor()
    resistivity :: Float64 = 0.0
    outline :: cryostat__description_2d__thermal_shield__unit__element__outline = cryostat__description_2d__thermal_shield__unit__element__outline()
    resistance :: Float64 = 0.0
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__annular__centreline
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit__annular
    outline_inner :: cryostat__description_2d__thermal_shield__unit__annular__outline_inner = cryostat__description_2d__thermal_shield__unit__annular__outline_inner()
    centreline :: cryostat__description_2d__thermal_shield__unit__annular__centreline = cryostat__description_2d__thermal_shield__unit__annular__centreline()
    thickness :: Array{Float64, 1} = zeros(Float64,(0))
    resistivity :: Float64 = 0.0
    outline_outer :: cryostat__description_2d__thermal_shield__unit__annular__outline_outer = cryostat__description_2d__thermal_shield__unit__annular__outline_outer()
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__element
    name :: String = ""
    j_tor :: cryostat__description_2d__cryostat__unit__element__j_tor = cryostat__description_2d__cryostat__unit__element__j_tor()
    resistivity :: Float64 = 0.0
    outline :: cryostat__description_2d__cryostat__unit__element__outline = cryostat__description_2d__cryostat__unit__element__outline()
    resistance :: Float64 = 0.0
end

Base.@kwdef mutable struct cryostat__code
    library :: StructArray{cryostat__code__library} = StructArray(cryostat__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield__unit
    name :: String = ""
    annular :: cryostat__description_2d__thermal_shield__unit__annular = cryostat__description_2d__thermal_shield__unit__annular()
    identifier :: String = ""
    element :: StructArray{cryostat__description_2d__thermal_shield__unit__element} = StructArray(cryostat__description_2d__thermal_shield__unit__element() for k in 1:1)
end

Base.@kwdef mutable struct cryostat__description_2d__thermal_shield
    unit :: StructArray{cryostat__description_2d__thermal_shield__unit} = StructArray(cryostat__description_2d__thermal_shield__unit() for k in 1:1)
    type :: cryostat__description_2d__thermal_shield__type = cryostat__description_2d__thermal_shield__type()
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__annular__outline_inner
    closed :: Int32 = 0
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit__annular
    outline_inner :: cryostat__description_2d__cryostat__unit__annular__outline_inner = cryostat__description_2d__cryostat__unit__annular__outline_inner()
    centreline :: cryostat__description_2d__cryostat__unit__annular__centreline = cryostat__description_2d__cryostat__unit__annular__centreline()
    thickness :: Array{Float64, 1} = zeros(Float64,(0))
    resistivity :: Float64 = 0.0
    outline_outer :: cryostat__description_2d__cryostat__unit__annular__outline_outer = cryostat__description_2d__cryostat__unit__annular__outline_outer()
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__unit
    name :: String = ""
    annular :: cryostat__description_2d__cryostat__unit__annular = cryostat__description_2d__cryostat__unit__annular()
    identifier :: String = ""
    element :: StructArray{cryostat__description_2d__cryostat__unit__element} = StructArray(cryostat__description_2d__cryostat__unit__element() for k in 1:1)
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct cryostat__description_2d__cryostat
    type :: cryostat__description_2d__cryostat__type = cryostat__description_2d__cryostat__type()
    unit :: StructArray{cryostat__description_2d__cryostat__unit} = StructArray(cryostat__description_2d__cryostat__unit() for k in 1:1)
end

Base.@kwdef mutable struct cryostat__description_2d
    cryostat :: cryostat__description_2d__cryostat = cryostat__description_2d__cryostat()
    thermal_shield :: cryostat__description_2d__thermal_shield = cryostat__description_2d__thermal_shield()
end

Base.@kwdef mutable struct cryostat
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: cryostat__ids_properties = cryostat__ids_properties()
    description_2d :: StructArray{cryostat__description_2d} = StructArray(cryostat__description_2d() for k in 1:1)
    code :: cryostat__code = cryostat__code()
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__momentum__poloidal
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral__state__particles
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral__particles
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__energy
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__particles
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__momentum__poloidal
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__electrons__energy
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral__energy
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_transport__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct core_transport__model__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_transport__model__code__output_flag
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral__state__energy
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__code
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: core_transport__model__code__output_flag = core_transport__model__code__output_flag()
    version :: String = ""
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__momentum__toroidal
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__total_ion_energy
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__momentum__parallel
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__momentum__toroidal
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__energy
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__electrons__particles
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__electrons
    particles :: core_transport__model__profiles_1d__electrons__particles = core_transport__model__profiles_1d__electrons__particles()
    energy :: core_transport__model__profiles_1d__electrons__energy = core_transport__model__profiles_1d__electrons__energy()
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__momentum_tor
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__grid_d
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__momentum__diamagnetic
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__grid_flux
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__momentum__parallel
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__grid_v
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_transport__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct core_transport__code
    library :: StructArray{core_transport__code__library} = StructArray(core_transport__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__particles
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__momentum__diamagnetic
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__momentum__radial
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__momentum
    parallel :: core_transport__model__profiles_1d__ion__momentum__parallel = core_transport__model__profiles_1d__ion__momentum__parallel()
    toroidal :: core_transport__model__profiles_1d__ion__momentum__toroidal = core_transport__model__profiles_1d__ion__momentum__toroidal()
    diamagnetic :: core_transport__model__profiles_1d__ion__momentum__diamagnetic = core_transport__model__profiles_1d__ion__momentum__diamagnetic()
    radial :: core_transport__model__profiles_1d__ion__momentum__radial = core_transport__model__profiles_1d__ion__momentum__radial()
    poloidal :: core_transport__model__profiles_1d__ion__momentum__poloidal = core_transport__model__profiles_1d__ion__momentum__poloidal()
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    particles :: core_transport__model__profiles_1d__neutral__state__particles = core_transport__model__profiles_1d__neutral__state__particles()
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    energy :: core_transport__model__profiles_1d__neutral__state__energy = core_transport__model__profiles_1d__neutral__state__energy()
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__neutral
    label :: String = ""
    particles :: core_transport__model__profiles_1d__neutral__particles = core_transport__model__profiles_1d__neutral__particles()
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    state :: StructArray{core_transport__model__profiles_1d__neutral__state} = StructArray(core_transport__model__profiles_1d__neutral__state() for k in 1:1)
    energy :: core_transport__model__profiles_1d__neutral__energy = core_transport__model__profiles_1d__neutral__energy()
    element :: StructArray{core_transport__model__profiles_1d__neutral__element} = StructArray(core_transport__model__profiles_1d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__momentum__radial
    v :: Array{Float64, 1} = zeros(Float64,(0))
    flow_damping_rate :: Array{Float64, 1} = zeros(Float64,(0))
    flux :: Array{Float64, 1} = zeros(Float64,(0))
    d :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state__momentum
    parallel :: core_transport__model__profiles_1d__ion__state__momentum__parallel = core_transport__model__profiles_1d__ion__state__momentum__parallel()
    toroidal :: core_transport__model__profiles_1d__ion__state__momentum__toroidal = core_transport__model__profiles_1d__ion__state__momentum__toroidal()
    diamagnetic :: core_transport__model__profiles_1d__ion__state__momentum__diamagnetic = core_transport__model__profiles_1d__ion__state__momentum__diamagnetic()
    radial :: core_transport__model__profiles_1d__ion__state__momentum__radial = core_transport__model__profiles_1d__ion__state__momentum__radial()
    poloidal :: core_transport__model__profiles_1d__ion__state__momentum__poloidal = core_transport__model__profiles_1d__ion__state__momentum__poloidal()
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion__state
    label :: String = ""
    particles :: core_transport__model__profiles_1d__ion__state__particles = core_transport__model__profiles_1d__ion__state__particles()
    vibrational_level :: Float64 = 0.0
    is_neutral :: Int32 = 0
    momentum :: core_transport__model__profiles_1d__ion__state__momentum = core_transport__model__profiles_1d__ion__state__momentum()
    energy :: core_transport__model__profiles_1d__ion__state__energy = core_transport__model__profiles_1d__ion__state__energy()
    z_min :: Float64 = 0.0
    electron_configuration :: String = ""
    vibrational_mode :: String = ""
    z_max :: Float64 = 0.0
    neutral_type :: core_transport__model__profiles_1d__ion__state__neutral_type = core_transport__model__profiles_1d__ion__state__neutral_type()
end

Base.@kwdef mutable struct core_transport__model__profiles_1d__ion
    z_ion :: Float64 = 0.0
    label :: String = ""
    particles :: core_transport__model__profiles_1d__ion__particles = core_transport__model__profiles_1d__ion__particles()
    multiple_states_flag :: Int32 = 0
    momentum :: core_transport__model__profiles_1d__ion__momentum = core_transport__model__profiles_1d__ion__momentum()
    state :: StructArray{core_transport__model__profiles_1d__ion__state} = StructArray(core_transport__model__profiles_1d__ion__state() for k in 1:1)
    energy :: core_transport__model__profiles_1d__ion__energy = core_transport__model__profiles_1d__ion__energy()
    neutral_index :: Int32 = 0
    element :: StructArray{core_transport__model__profiles_1d__ion__element} = StructArray(core_transport__model__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct core_transport__model__profiles_1d
    time :: Float64 = 0.0
    grid_v :: core_transport__model__profiles_1d__grid_v = core_transport__model__profiles_1d__grid_v()
    neutral :: StructArray{core_transport__model__profiles_1d__neutral} = StructArray(core_transport__model__profiles_1d__neutral() for k in 1:1)
    ion :: StructArray{core_transport__model__profiles_1d__ion} = StructArray(core_transport__model__profiles_1d__ion() for k in 1:1)
    momentum_tor :: core_transport__model__profiles_1d__momentum_tor = core_transport__model__profiles_1d__momentum_tor()
    electrons :: core_transport__model__profiles_1d__electrons = core_transport__model__profiles_1d__electrons()
    conductivity_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_radial :: Array{Float64, 1} = zeros(Float64,(0))
    grid_d :: core_transport__model__profiles_1d__grid_d = core_transport__model__profiles_1d__grid_d()
    total_ion_energy :: core_transport__model__profiles_1d__total_ion_energy = core_transport__model__profiles_1d__total_ion_energy()
    grid_flux :: core_transport__model__profiles_1d__grid_flux = core_transport__model__profiles_1d__grid_flux()
end

Base.@kwdef mutable struct core_transport__model
    flux_multiplier :: Float64 = 0.0
    code :: core_transport__model__code = core_transport__model__code()
    profiles_1d :: StructArray{core_transport__model__profiles_1d} = StructArray(core_transport__model__profiles_1d() for k in 1:1)
    identifier :: core_transport__model__identifier = core_transport__model__identifier()
    comment :: String = ""
end

Base.@kwdef mutable struct core_transport__ids_properties
    provider :: String = ""
    version_put :: core_transport__ids_properties__version_put = core_transport__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct core_transport
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: core_transport__ids_properties = core_transport__ids_properties()
    vacuum_toroidal_field :: core_transport__vacuum_toroidal_field = core_transport__vacuum_toroidal_field()
    code :: core_transport__code = core_transport__code()
    model :: StructArray{core_transport__model} = StructArray(core_transport__model() for k in 1:1)
end

Base.@kwdef mutable struct core_sources__source__species__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__electrons__particles_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_sources__source__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_sources__source__species__ion__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    z_min :: Float64 = 0.0
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__electrons__energy_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__momentum__toroidal_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__energy_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__state__particles_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__state__energy_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__total_ion_energy_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__global_quantities__electrons
    particles :: Float64 = 0.0
    power :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__particles_decomposed
    explicit_part :: Array{Float64, 1} = zeros(Float64,(0))
    implicit_part :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__electrons
    particles_inside :: Array{Float64, 1} = zeros(Float64,(0))
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    particles_decomposed :: core_sources__source__profiles_1d__electrons__particles_decomposed = core_sources__source__profiles_1d__electrons__particles_decomposed()
    energy_decomposed :: core_sources__source__profiles_1d__electrons__energy_decomposed = core_sources__source__profiles_1d__electrons__energy_decomposed()
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    power_inside :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    neutral_type :: core_sources__source__profiles_1d__neutral__state__neutral_type = core_sources__source__profiles_1d__neutral__state__neutral_type()
    energy :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__momentum
    toroidal_decomposed :: core_sources__source__profiles_1d__ion__momentum__toroidal_decomposed = core_sources__source__profiles_1d__ion__momentum__toroidal_decomposed()
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__species__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_sources__source__global_quantities
    time :: Float64 = 0.0
    current_parallel :: Float64 = 0.0
    total_ion_power :: Float64 = 0.0
    torque_tor :: Float64 = 0.0
    power :: Float64 = 0.0
    total_ion_particles :: Float64 = 0.0
    electrons :: core_sources__source__global_quantities__electrons = core_sources__source__global_quantities__electrons()
end

Base.@kwdef mutable struct core_sources__source__species__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__species__ion
    state :: core_sources__source__species__ion__state = core_sources__source__species__ion__state()
    z_ion :: Float64 = 0.0
    label :: String = ""
    element :: StructArray{core_sources__source__species__ion__element} = StructArray(core_sources__source__species__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct core_sources__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion__state
    label :: String = ""
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_level :: Float64 = 0.0
    is_neutral :: Int32 = 0
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    electron_configuration :: String = ""
    particles_decomposed :: core_sources__source__profiles_1d__ion__state__particles_decomposed = core_sources__source__profiles_1d__ion__state__particles_decomposed()
    vibrational_mode :: String = ""
    z_max :: Float64 = 0.0
    energy_decomposed :: core_sources__source__profiles_1d__ion__state__energy_decomposed = core_sources__source__profiles_1d__ion__state__energy_decomposed()
    neutral_type :: core_sources__source__profiles_1d__ion__state__neutral_type = core_sources__source__profiles_1d__ion__state__neutral_type()
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__ion
    label :: String = ""
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    momentum :: core_sources__source__profiles_1d__ion__momentum = core_sources__source__profiles_1d__ion__momentum()
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    neutral_index :: Int32 = 0
    particles_decomposed :: core_sources__source__profiles_1d__ion__particles_decomposed = core_sources__source__profiles_1d__ion__particles_decomposed()
    state :: StructArray{core_sources__source__profiles_1d__ion__state} = StructArray(core_sources__source__profiles_1d__ion__state() for k in 1:1)
    z_ion :: Float64 = 0.0
    energy_decomposed :: core_sources__source__profiles_1d__ion__energy_decomposed = core_sources__source__profiles_1d__ion__energy_decomposed()
    element :: StructArray{core_sources__source__profiles_1d__ion__element} = StructArray(core_sources__source__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct core_sources__code
    library :: StructArray{core_sources__code__library} = StructArray(core_sources__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct core_sources__source__species__neutral__state
    label :: String = ""
    electron_configuration :: String = ""
    vibrational_level :: Float64 = 0.0
    vibrational_mode :: String = ""
    neutral_type :: core_sources__source__species__neutral__state__neutral_type = core_sources__source__species__neutral__state__neutral_type()
end

Base.@kwdef mutable struct core_sources__source__species__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_sources__source__profiles_1d__neutral
    label :: String = ""
    particles :: Array{Float64, 1} = zeros(Float64,(0))
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    state :: StructArray{core_sources__source__profiles_1d__neutral__state} = StructArray(core_sources__source__profiles_1d__neutral__state() for k in 1:1)
    energy :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{core_sources__source__profiles_1d__neutral__element} = StructArray(core_sources__source__profiles_1d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct core_sources__source__profiles_1d
    time :: Float64 = 0.0
    neutral :: StructArray{core_sources__source__profiles_1d__neutral} = StructArray(core_sources__source__profiles_1d__neutral() for k in 1:1)
    torque_tor_inside :: Array{Float64, 1} = zeros(Float64,(0))
    ion :: StructArray{core_sources__source__profiles_1d__ion} = StructArray(core_sources__source__profiles_1d__ion() for k in 1:1)
    total_ion_energy_decomposed :: core_sources__source__profiles_1d__total_ion_energy_decomposed = core_sources__source__profiles_1d__total_ion_energy_decomposed()
    current_parallel_inside :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: core_sources__source__profiles_1d__electrons = core_sources__source__profiles_1d__electrons()
    conductivity_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor_j_cross_b_field :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: core_sources__source__profiles_1d__grid = core_sources__source__profiles_1d__grid()
    total_ion_energy :: Array{Float64, 1} = zeros(Float64,(0))
    j_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    total_ion_power_inside :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_sources__source__species__neutral
    state :: core_sources__source__species__neutral__state = core_sources__source__species__neutral__state()
    label :: String = ""
    element :: StructArray{core_sources__source__species__neutral__element} = StructArray(core_sources__source__species__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct core_sources__source__species
    ion :: core_sources__source__species__ion = core_sources__source__species__ion()
    type :: core_sources__source__species__type = core_sources__source__species__type()
    neutral :: core_sources__source__species__neutral = core_sources__source__species__neutral()
end

Base.@kwdef mutable struct core_sources__source
    global_quantities :: StructArray{core_sources__source__global_quantities} = StructArray(core_sources__source__global_quantities() for k in 1:1)
    profiles_1d :: StructArray{core_sources__source__profiles_1d} = StructArray(core_sources__source__profiles_1d() for k in 1:1)
    identifier :: core_sources__source__identifier = core_sources__source__identifier()
    species :: core_sources__source__species = core_sources__source__species()
end

Base.@kwdef mutable struct core_sources__ids_properties
    provider :: String = ""
    version_put :: core_sources__ids_properties__version_put = core_sources__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct core_sources
    source :: StructArray{core_sources__source} = StructArray(core_sources__source() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: core_sources__ids_properties = core_sources__ids_properties()
    vacuum_toroidal_field :: core_sources__vacuum_toroidal_field = core_sources__vacuum_toroidal_field()
    code :: core_sources__code = core_sources__code()
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_profiles__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__state__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state__density_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method = core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__e_field
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__ids_properties
    provider :: String = ""
    version_put :: core_profiles__ids_properties__version_put = core_profiles__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__temperature_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    rotation_frequency_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    electron_configuration :: String = ""
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_mode :: String = ""
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    z_average_square_1d :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: core_profiles__profiles_1d__ion__state__velocity = core_profiles__profiles_1d__ion__state__velocity()
    z_average :: Float64 = 0.0
    z_max :: Float64 = 0.0
    z_square_average :: Float64 = 0.0
    ionisation_potential :: Float64 = 0.0
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    z_average_1d :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: core_profiles__profiles_1d__ion__state__density_fit = core_profiles__profiles_1d__ion__state__density_fit()
end

Base.@kwdef mutable struct core_profiles__global_quantities
    beta_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    resistive_psi_losses :: Array{Float64, 1} = zeros(Float64,(0))
    t_i_average_peaking :: Array{Float64, 1} = zeros(Float64,(0))
    li_3 :: Array{Float64, 1} = zeros(Float64,(0))
    ip :: Array{Float64, 1} = zeros(Float64,(0))
    t_e_peaking :: Array{Float64, 1} = zeros(Float64,(0))
    beta_tor :: Array{Float64, 1} = zeros(Float64,(0))
    z_eff_resistive :: Array{Float64, 1} = zeros(Float64,(0))
    energy_diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    ejima :: Array{Float64, 1} = zeros(Float64,(0))
    li :: Array{Float64, 1} = zeros(Float64,(0))
    v_loop :: Array{Float64, 1} = zeros(Float64,(0))
    current_non_inductive :: Array{Float64, 1} = zeros(Float64,(0))
    beta_pol :: Array{Float64, 1} = zeros(Float64,(0))
    current_bootstrap :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    electron_configuration :: String = ""
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_mode :: String = ""
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: core_profiles__profiles_1d__neutral__state__velocity = core_profiles__profiles_1d__neutral__state__velocity()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_type :: core_profiles__profiles_1d__neutral__state__neutral_type = core_profiles__profiles_1d__neutral__state__neutral_type()
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral
    label :: String = ""
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    ion_index :: Int32 = 0
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{core_profiles__profiles_1d__neutral__state} = StructArray(core_profiles__profiles_1d__neutral__state() for k in 1:1)
    velocity :: core_profiles__profiles_1d__neutral__velocity = core_profiles__profiles_1d__neutral__velocity()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{core_profiles__profiles_1d__neutral__element} = StructArray(core_profiles__profiles_1d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct core_profiles__profiles_1d__zeff_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method = core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__density_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method = core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__temperature_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__t_i_average_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method = core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion
    label :: String = ""
    rotation_frequency_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_validity :: Int32 = 0
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    z_ion_1d :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    neutral_index :: Int32 = 0
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    density_validity :: Int32 = 0
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{core_profiles__profiles_1d__ion__state} = StructArray(core_profiles__profiles_1d__ion__state() for k in 1:1)
    velocity :: core_profiles__profiles_1d__ion__velocity = core_profiles__profiles_1d__ion__velocity()
    z_ion :: Float64 = 0.0
    temperature_fit :: core_profiles__profiles_1d__ion__temperature_fit = core_profiles__profiles_1d__ion__temperature_fit()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: core_profiles__profiles_1d__ion__density_fit = core_profiles__profiles_1d__ion__density_fit()
    z_ion_square_1d :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{core_profiles__profiles_1d__ion__element} = StructArray(core_profiles__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct core_profiles__code
    library :: StructArray{core_profiles__code__library} = StructArray(core_profiles__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__density_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method = core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_validity :: Int32 = 0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    density_validity :: Int32 = 0
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: core_profiles__profiles_1d__electrons__velocity = core_profiles__profiles_1d__electrons__velocity()
    temperature_fit :: core_profiles__profiles_1d__electrons__temperature_fit = core_profiles__profiles_1d__electrons__temperature_fit()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    collisionality_norm :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: core_profiles__profiles_1d__electrons__density_fit = core_profiles__profiles_1d__electrons__density_fit()
end

Base.@kwdef mutable struct core_profiles__profiles_1d
    pressure_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    t_i_average_fit :: core_profiles__profiles_1d__t_i_average_fit = core_profiles__profiles_1d__t_i_average_fit()
    neutral :: StructArray{core_profiles__profiles_1d__neutral} = StructArray(core_profiles__profiles_1d__neutral() for k in 1:1)
    n_i_thermal_total :: Array{Float64, 1} = zeros(Float64,(0))
    magnetic_shear :: Array{Float64, 1} = zeros(Float64,(0))
    ion :: StructArray{core_profiles__profiles_1d__ion} = StructArray(core_profiles__profiles_1d__ion() for k in 1:1)
    j_total :: Array{Float64, 1} = zeros(Float64,(0))
    rotation_frequency_tor_sonic :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    j_tor :: Array{Float64, 1} = zeros(Float64,(0))
    t_i_average :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: core_profiles__profiles_1d__electrons = core_profiles__profiles_1d__electrons()
    q :: Array{Float64, 1} = zeros(Float64,(0))
    j_non_inductive :: Array{Float64, 1} = zeros(Float64,(0))
    conductivity_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    current_parallel_inside :: Array{Float64, 1} = zeros(Float64,(0))
    j_ohmic :: Array{Float64, 1} = zeros(Float64,(0))
    phi_potential :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: core_profiles__profiles_1d__grid = core_profiles__profiles_1d__grid()
    j_bootstrap :: Array{Float64, 1} = zeros(Float64,(0))
    zeff_fit :: core_profiles__profiles_1d__zeff_fit = core_profiles__profiles_1d__zeff_fit()
    pressure_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    e_field :: core_profiles__profiles_1d__e_field = core_profiles__profiles_1d__e_field()
    zeff :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_profiles
    time :: Array{Float64, 1} = zeros(Float64,(0))
    vacuum_toroidal_field :: core_profiles__vacuum_toroidal_field = core_profiles__vacuum_toroidal_field()
    ids_properties :: core_profiles__ids_properties = core_profiles__ids_properties()
    code :: core_profiles__code = core_profiles__code()
    profiles_1d :: StructArray{core_profiles__profiles_1d} = StructArray(core_profiles__profiles_1d() for k in 1:1)
    global_quantities :: core_profiles__global_quantities = core_profiles__global_quantities()
end

Base.@kwdef mutable struct core_instant_changes__change__identifier
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__electrons__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__electrons__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__temperature_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__e_field
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__neutral__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__density_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__ion__density_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__ion__density_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__ids_properties
    provider :: String = ""
    version_put :: core_instant_changes__ids_properties__version_put = core_instant_changes__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__neutral__state__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__electrons__temperature_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__electrons__temperature_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__electrons__temperature_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__electrons__temperature_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__electrons__density_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__electrons__density_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__electrons__density_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__t_i_average_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__vacuum_toroidal_field
    b0 :: Array{Float64, 1} = zeros(Float64,(0))
    r0 :: Float64 = 0.0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__neutral__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__state__density_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__state__density_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__ion__state__density_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__ion__state__density_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__electrons
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_validity :: Int32 = 0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    density_validity :: Int32 = 0
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: core_instant_changes__change__profiles_1d__electrons__velocity = core_instant_changes__change__profiles_1d__electrons__velocity()
    temperature_fit :: core_instant_changes__change__profiles_1d__electrons__temperature_fit = core_instant_changes__change__profiles_1d__electrons__temperature_fit()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    collisionality_norm :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: core_instant_changes__change__profiles_1d__electrons__density_fit = core_instant_changes__change__profiles_1d__electrons__density_fit()
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__neutral__state__neutral_type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__neutral__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    electron_configuration :: String = ""
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_mode :: String = ""
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: core_instant_changes__change__profiles_1d__neutral__state__velocity = core_instant_changes__change__profiles_1d__neutral__state__velocity()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_type :: core_instant_changes__change__profiles_1d__neutral__state__neutral_type = core_instant_changes__change__profiles_1d__neutral__state__neutral_type()
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__neutral
    label :: String = ""
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    ion_index :: Int32 = 0
    multiple_states_flag :: Int32 = 0
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{core_instant_changes__change__profiles_1d__neutral__state} = StructArray(core_instant_changes__change__profiles_1d__neutral__state() for k in 1:1)
    velocity :: core_instant_changes__change__profiles_1d__neutral__velocity = core_instant_changes__change__profiles_1d__neutral__velocity()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{core_instant_changes__change__profiles_1d__neutral__element} = StructArray(core_instant_changes__change__profiles_1d__neutral__element() for k in 1:1)
end

Base.@kwdef mutable struct core_instant_changes__code
    library :: StructArray{core_instant_changes__code__library} = StructArray(core_instant_changes__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__zeff_fit__time_measurement_slice_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__zeff_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__zeff_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__zeff_fit__time_measurement_slice_method()
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__temperature_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__ion__temperature_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__ion__temperature_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__state__velocity
    parallel :: Array{Float64, 1} = zeros(Float64,(0))
    toroidal :: Array{Float64, 1} = zeros(Float64,(0))
    diamagnetic :: Array{Float64, 1} = zeros(Float64,(0))
    radial :: Array{Float64, 1} = zeros(Float64,(0))
    poloidal :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__t_i_average_fit
    chi_squared :: Array{Float64, 1} = zeros(Float64,(0))
    weight :: Array{Float64, 1} = zeros(Float64,(0))
    parameters :: String = ""
    source :: Array{String, 1} = String[]
    measured :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_slice_method :: core_instant_changes__change__profiles_1d__t_i_average_fit__time_measurement_slice_method = core_instant_changes__change__profiles_1d__t_i_average_fit__time_measurement_slice_method()
    reconstructed :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    time_measurement_width :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__grid
    psi :: Array{Float64, 1} = zeros(Float64,(0))
    psi_boundary :: Float64 = 0.0
    volume :: Array{Float64, 1} = zeros(Float64,(0))
    area :: Array{Float64, 1} = zeros(Float64,(0))
    rho_pol_norm :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor_norm :: Array{Float64, 1} = zeros(Float64,(0))
    surface :: Array{Float64, 1} = zeros(Float64,(0))
    rho_tor :: Array{Float64, 1} = zeros(Float64,(0))
    psi_magnetic_axis :: Float64 = 0.0
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion__state
    label :: String = ""
    vibrational_level :: Float64 = 0.0
    rotation_frequency_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    z_min :: Float64 = 0.0
    electron_configuration :: String = ""
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    vibrational_mode :: String = ""
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    z_average_square_1d :: Array{Float64, 1} = zeros(Float64,(0))
    velocity :: core_instant_changes__change__profiles_1d__ion__state__velocity = core_instant_changes__change__profiles_1d__ion__state__velocity()
    z_max :: Float64 = 0.0
    z_average :: Float64 = 0.0
    z_square_average :: Float64 = 0.0
    ionisation_potential :: Float64 = 0.0
    density :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    z_average_1d :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: core_instant_changes__change__profiles_1d__ion__state__density_fit = core_instant_changes__change__profiles_1d__ion__state__density_fit()
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d__ion
    label :: String = ""
    rotation_frequency_tor :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_tor :: Array{Float64, 1} = zeros(Float64,(0))
    temperature_validity :: Int32 = 0
    temperature :: Array{Float64, 1} = zeros(Float64,(0))
    z_ion_1d :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    multiple_states_flag :: Int32 = 0
    pressure_fast_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    neutral_index :: Int32 = 0
    pressure :: Array{Float64, 1} = zeros(Float64,(0))
    density_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    density_validity :: Int32 = 0
    pressure_fast_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    state :: StructArray{core_instant_changes__change__profiles_1d__ion__state} = StructArray(core_instant_changes__change__profiles_1d__ion__state() for k in 1:1)
    velocity :: core_instant_changes__change__profiles_1d__ion__velocity = core_instant_changes__change__profiles_1d__ion__velocity()
    z_ion :: Float64 = 0.0
    temperature_fit :: core_instant_changes__change__profiles_1d__ion__temperature_fit = core_instant_changes__change__profiles_1d__ion__temperature_fit()
    density :: Array{Float64, 1} = zeros(Float64,(0))
    velocity_pol :: Array{Float64, 1} = zeros(Float64,(0))
    density_fast :: Array{Float64, 1} = zeros(Float64,(0))
    density_fit :: core_instant_changes__change__profiles_1d__ion__density_fit = core_instant_changes__change__profiles_1d__ion__density_fit()
    z_ion_square_1d :: Array{Float64, 1} = zeros(Float64,(0))
    element :: StructArray{core_instant_changes__change__profiles_1d__ion__element} = StructArray(core_instant_changes__change__profiles_1d__ion__element() for k in 1:1)
end

Base.@kwdef mutable struct core_instant_changes__change__profiles_1d
    pressure_ion_total :: Array{Float64, 1} = zeros(Float64,(0))
    time :: Float64 = 0.0
    t_i_average_fit :: core_instant_changes__change__profiles_1d__t_i_average_fit = core_instant_changes__change__profiles_1d__t_i_average_fit()
    neutral :: StructArray{core_instant_changes__change__profiles_1d__neutral} = StructArray(core_instant_changes__change__profiles_1d__neutral() for k in 1:1)
    n_i_thermal_total :: Array{Float64, 1} = zeros(Float64,(0))
    magnetic_shear :: Array{Float64, 1} = zeros(Float64,(0))
    ion :: StructArray{core_instant_changes__change__profiles_1d__ion} = StructArray(core_instant_changes__change__profiles_1d__ion() for k in 1:1)
    j_total :: Array{Float64, 1} = zeros(Float64,(0))
    rotation_frequency_tor_sonic :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_thermal :: Array{Float64, 1} = zeros(Float64,(0))
    j_tor :: Array{Float64, 1} = zeros(Float64,(0))
    current_parallel_inside :: Array{Float64, 1} = zeros(Float64,(0))
    e_field_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    j_non_inductive :: Array{Float64, 1} = zeros(Float64,(0))
    pressure_perpendicular :: Array{Float64, 1} = zeros(Float64,(0))
    electrons :: core_instant_changes__change__profiles_1d__electrons = core_instant_changes__change__profiles_1d__electrons()
    conductivity_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    momentum_tor :: Array{Float64, 1} = zeros(Float64,(0))
    q :: Array{Float64, 1} = zeros(Float64,(0))
    t_i_average :: Array{Float64, 1} = zeros(Float64,(0))
    j_ohmic :: Array{Float64, 1} = zeros(Float64,(0))
    grid :: core_instant_changes__change__profiles_1d__grid = core_instant_changes__change__profiles_1d__grid()
    phi_potential :: Array{Float64, 1} = zeros(Float64,(0))
    j_bootstrap :: Array{Float64, 1} = zeros(Float64,(0))
    zeff_fit :: core_instant_changes__change__profiles_1d__zeff_fit = core_instant_changes__change__profiles_1d__zeff_fit()
    pressure_parallel :: Array{Float64, 1} = zeros(Float64,(0))
    e_field :: core_instant_changes__change__profiles_1d__e_field = core_instant_changes__change__profiles_1d__e_field()
    zeff :: Array{Float64, 1} = zeros(Float64,(0))
    n_i_total_over_n_e :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct core_instant_changes__change
    profiles_1d :: StructArray{core_instant_changes__change__profiles_1d} = StructArray(core_instant_changes__change__profiles_1d() for k in 1:1)
    identifier :: core_instant_changes__change__identifier = core_instant_changes__change__identifier()
end

Base.@kwdef mutable struct core_instant_changes
    change :: StructArray{core_instant_changes__change} = StructArray(core_instant_changes__change() for k in 1:1)
    time :: Array{Float64, 1} = zeros(Float64,(0))
    vacuum_toroidal_field :: core_instant_changes__vacuum_toroidal_field = core_instant_changes__vacuum_toroidal_field()
    ids_properties :: core_instant_changes__ids_properties = core_instant_changes__ids_properties()
    code :: core_instant_changes__code = core_instant_changes__code()
end

Base.@kwdef mutable struct controllers__linear_controller__outputs
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct controllers__nonlinear_controller__inputs
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__pid__i
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__nonlinear_controller__outputs
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__statespace__b
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__statespace__d
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct controllers__linear_controller__statespace__c
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__pid__d
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__pid__p
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__nonlinear_controller
    input_names :: Array{String, 1} = String[]
    name :: String = ""
    outputs :: controllers__nonlinear_controller__outputs = controllers__nonlinear_controller__outputs()
    controller_class :: String = ""
    output_names :: Array{String, 1} = String[]
    description :: String = ""
    inputs :: controllers__nonlinear_controller__inputs = controllers__nonlinear_controller__inputs()
end

Base.@kwdef mutable struct controllers__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct controllers__ids_properties
    provider :: String = ""
    version_put :: controllers__ids_properties__version_put = controllers__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct controllers__linear_controller__inputs
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__statespace__deltat
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct controllers__code
    library :: StructArray{controllers__code__library} = StructArray(controllers__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct controllers__linear_controller__pid__tau
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct controllers__linear_controller__pid
    i :: controllers__linear_controller__pid__i = controllers__linear_controller__pid__i()
    tau :: controllers__linear_controller__pid__tau = controllers__linear_controller__pid__tau()
    p :: controllers__linear_controller__pid__p = controllers__linear_controller__pid__p()
    d :: controllers__linear_controller__pid__d = controllers__linear_controller__pid__d()
end

Base.@kwdef mutable struct controllers__linear_controller__statespace__a
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 3} = zeros(Float64,(0,0,0))
end

Base.@kwdef mutable struct controllers__linear_controller__statespace
    d :: controllers__linear_controller__statespace__d = controllers__linear_controller__statespace__d()
    state_names :: Array{String, 1} = String[]
    c :: controllers__linear_controller__statespace__c = controllers__linear_controller__statespace__c()
    b :: controllers__linear_controller__statespace__b = controllers__linear_controller__statespace__b()
    a :: controllers__linear_controller__statespace__a = controllers__linear_controller__statespace__a()
    deltat :: controllers__linear_controller__statespace__deltat = controllers__linear_controller__statespace__deltat()
end

Base.@kwdef mutable struct controllers__linear_controller
    input_names :: Array{String, 1} = String[]
    controller_class :: String = ""
    outputs :: controllers__linear_controller__outputs = controllers__linear_controller__outputs()
    pid :: controllers__linear_controller__pid = controllers__linear_controller__pid()
    output_names :: Array{String, 1} = String[]
    name :: String = ""
    statespace :: controllers__linear_controller__statespace = controllers__linear_controller__statespace()
    description :: String = ""
    inputs :: controllers__linear_controller__inputs = controllers__linear_controller__inputs()
end

Base.@kwdef mutable struct controllers
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: controllers__ids_properties = controllers__ids_properties()
    nonlinear_controller :: StructArray{controllers__nonlinear_controller} = StructArray(controllers__nonlinear_controller() for k in 1:1)
    linear_controller :: StructArray{controllers__linear_controller} = StructArray(controllers__linear_controller() for k in 1:1)
    code :: controllers__code = controllers__code()
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__current
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__elements__centres
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__elements__end_points
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__elements__start_points
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__voltage
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__elements__intermediate_points
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__elements
    start_points :: coils_non_axisymmetric__coil__conductor__elements__start_points = coils_non_axisymmetric__coil__conductor__elements__start_points()
    names :: Array{String, 1} = String[]
    types :: Array{Int32, 1} = zeros(Int32,(0))
    end_points :: coils_non_axisymmetric__coil__conductor__elements__end_points = coils_non_axisymmetric__coil__conductor__elements__end_points()
    intermediate_points :: coils_non_axisymmetric__coil__conductor__elements__intermediate_points = coils_non_axisymmetric__coil__conductor__elements__intermediate_points()
    centres :: coils_non_axisymmetric__coil__conductor__elements__centres = coils_non_axisymmetric__coil__conductor__elements__centres()
end

Base.@kwdef mutable struct coils_non_axisymmetric__code
    library :: StructArray{coils_non_axisymmetric__code__library} = StructArray(coils_non_axisymmetric__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor__cross_section
    delta_r :: Array{Float64, 1} = zeros(Float64,(0))
    delta_z :: Array{Float64, 1} = zeros(Float64,(0))
    delta_phi :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct coils_non_axisymmetric__ids_properties
    provider :: String = ""
    version_put :: coils_non_axisymmetric__ids_properties__version_put = coils_non_axisymmetric__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil__conductor
    current :: coils_non_axisymmetric__coil__conductor__current = coils_non_axisymmetric__coil__conductor__current()
    elements :: coils_non_axisymmetric__coil__conductor__elements = coils_non_axisymmetric__coil__conductor__elements()
    cross_section :: coils_non_axisymmetric__coil__conductor__cross_section = coils_non_axisymmetric__coil__conductor__cross_section()
    resistance :: Float64 = 0.0
    voltage :: coils_non_axisymmetric__coil__conductor__voltage = coils_non_axisymmetric__coil__conductor__voltage()
end

Base.@kwdef mutable struct coils_non_axisymmetric__coil
    current :: coils_non_axisymmetric__coil__current = coils_non_axisymmetric__coil__current()
    name :: String = ""
    conductor :: StructArray{coils_non_axisymmetric__coil__conductor} = StructArray(coils_non_axisymmetric__coil__conductor() for k in 1:1)
    turns :: Float64 = 0.0
    resistance :: Float64 = 0.0
    identifier :: String = ""
    voltage :: coils_non_axisymmetric__coil__voltage = coils_non_axisymmetric__coil__voltage()
end

Base.@kwdef mutable struct coils_non_axisymmetric
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: coils_non_axisymmetric__ids_properties = coils_non_axisymmetric__ids_properties()
    coil :: StructArray{coils_non_axisymmetric__coil} = StructArray(coils_non_axisymmetric__coil() for k in 1:1)
    latency :: Float64 = 0.0
    code :: coils_non_axisymmetric__code = coils_non_axisymmetric__code()
end

Base.@kwdef mutable struct charge_exchange__channel__bes__radiances
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct charge_exchange__channel__t_i_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__etendue_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__processed_line__radiance
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__ion__n_i_over_n_e_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__momentum_tor
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__ion__t_i_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__processed_line__width
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__ion_fast__radiance
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__ion__velocity_pol
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__zeff_line_average_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct charge_exchange__channel__position__z
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__momentum_tor_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct charge_exchange__channel__ion__velocity_tor
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct charge_exchange__channel__ion__n_i_over_n_e
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__bes__doppler_shift
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__t_i_average_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__zeff_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct charge_exchange__channel__bes__lorentz_shift
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__bes
    label :: String = ""
    lorentz_shift :: charge_exchange__channel__bes__lorentz_shift = charge_exchange__channel__bes__lorentz_shift()
    z_n :: Float64 = 0.0
    transition_wavelength :: Float64 = 0.0
    radiances :: charge_exchange__channel__bes__radiances = charge_exchange__channel__bes__radiances()
    doppler_shift :: charge_exchange__channel__bes__doppler_shift = charge_exchange__channel__bes__doppler_shift()
    a :: Float64 = 0.0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct charge_exchange__channel__ion__velocity_pol_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__radiance_spectral
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__intensity_spectrum
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct charge_exchange__channel__ion__t_i
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__zeff
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__zeff_line_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__processed_line__shift
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__position__r
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__aperture
    surface :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: charge_exchange__aperture__centre = charge_exchange__aperture__centre()
    outline :: charge_exchange__aperture__outline = charge_exchange__aperture__outline()
    radius :: Float64 = 0.0
    x3_unit_vector :: charge_exchange__aperture__x3_unit_vector = charge_exchange__aperture__x3_unit_vector()
    x2_unit_vector :: charge_exchange__aperture__x2_unit_vector = charge_exchange__aperture__x2_unit_vector()
    x1_width :: Float64 = 0.0
    geometry_type :: Int32 = 0
    x1_unit_vector :: charge_exchange__aperture__x1_unit_vector = charge_exchange__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__processed_line__intensity
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__radiance_continuum
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 2} = zeros(Float64,(0,0))
end

Base.@kwdef mutable struct charge_exchange__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct charge_exchange__ids_properties
    provider :: String = ""
    version_put :: charge_exchange__ids_properties__version_put = charge_exchange__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__position__phi
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct charge_exchange__channel__position
    phi :: charge_exchange__channel__position__phi = charge_exchange__channel__position__phi()
    r :: charge_exchange__channel__position__r = charge_exchange__channel__position__r()
    z :: charge_exchange__channel__position__z = charge_exchange__channel__position__z()
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum__processed_line
    intensity :: charge_exchange__channel__spectrum__processed_line__intensity = charge_exchange__channel__spectrum__processed_line__intensity()
    label :: String = ""
    radiance :: charge_exchange__channel__spectrum__processed_line__radiance = charge_exchange__channel__spectrum__processed_line__radiance()
    wavelength_central :: Float64 = 0.0
    shift :: charge_exchange__channel__spectrum__processed_line__shift = charge_exchange__channel__spectrum__processed_line__shift()
    width :: charge_exchange__channel__spectrum__processed_line__width = charge_exchange__channel__spectrum__processed_line__width()
end

Base.@kwdef mutable struct charge_exchange__channel__spectrum
    processed_line :: StructArray{charge_exchange__channel__spectrum__processed_line} = StructArray(charge_exchange__channel__spectrum__processed_line() for k in 1:1)
    wavelength_calibration_date :: String = ""
    intensity_spectrum :: charge_exchange__channel__spectrum__intensity_spectrum = charge_exchange__channel__spectrum__intensity_spectrum()
    radiance_spectral :: charge_exchange__channel__spectrum__radiance_spectral = charge_exchange__channel__spectrum__radiance_spectral()
    radiance_calibration_date :: String = ""
    wavelengths :: Array{Float64, 1} = zeros(Float64,(0))
    slit_width :: Float64 = 0.0
    exposure_time :: Float64 = 0.0
    grating :: Float64 = 0.0
    radiance_calibration :: Array{Float64, 1} = zeros(Float64,(0))
    radiance_continuum :: charge_exchange__channel__spectrum__radiance_continuum = charge_exchange__channel__spectrum__radiance_continuum()
end

Base.@kwdef mutable struct charge_exchange__channel__ion__velocity_tor_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__ion
    n_i_over_n_e :: charge_exchange__channel__ion__n_i_over_n_e = charge_exchange__channel__ion__n_i_over_n_e()
    label :: String = ""
    z_n :: Float64 = 0.0
    velocity_tor :: charge_exchange__channel__ion__velocity_tor = charge_exchange__channel__ion__velocity_tor()
    t_i_method :: charge_exchange__channel__ion__t_i_method = charge_exchange__channel__ion__t_i_method()
    a :: Float64 = 0.0
    velocity_tor_method :: charge_exchange__channel__ion__velocity_tor_method = charge_exchange__channel__ion__velocity_tor_method()
    velocity_pol_method :: charge_exchange__channel__ion__velocity_pol_method = charge_exchange__channel__ion__velocity_pol_method()
    t_i :: charge_exchange__channel__ion__t_i = charge_exchange__channel__ion__t_i()
    n_i_over_n_e_method :: charge_exchange__channel__ion__n_i_over_n_e_method = charge_exchange__channel__ion__n_i_over_n_e_method()
    z_ion :: Float64 = 0.0
    velocity_pol :: charge_exchange__channel__ion__velocity_pol = charge_exchange__channel__ion__velocity_pol()
end

Base.@kwdef mutable struct charge_exchange__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct charge_exchange__code
    library :: StructArray{charge_exchange__code__library} = StructArray(charge_exchange__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct charge_exchange__channel__ion_fast__radiance_spectral_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct charge_exchange__channel__ion_fast
    label :: String = ""
    radiance_spectral_method :: charge_exchange__channel__ion_fast__radiance_spectral_method = charge_exchange__channel__ion_fast__radiance_spectral_method()
    z_n :: Float64 = 0.0
    transition_wavelength :: Float64 = 0.0
    radiance :: charge_exchange__channel__ion_fast__radiance = charge_exchange__channel__ion_fast__radiance()
    a :: Float64 = 0.0
    z_ion :: Float64 = 0.0
end

Base.@kwdef mutable struct charge_exchange__channel
    ion :: StructArray{charge_exchange__channel__ion} = StructArray(charge_exchange__channel__ion() for k in 1:1)
    name :: String = ""
    ion_fast :: StructArray{charge_exchange__channel__ion_fast} = StructArray(charge_exchange__channel__ion_fast() for k in 1:1)
    position :: charge_exchange__channel__position = charge_exchange__channel__position()
    zeff_line_average_method :: charge_exchange__channel__zeff_line_average_method = charge_exchange__channel__zeff_line_average_method()
    momentum_tor :: charge_exchange__channel__momentum_tor = charge_exchange__channel__momentum_tor()
    t_i_average :: charge_exchange__channel__t_i_average = charge_exchange__channel__t_i_average()
    zeff_method :: charge_exchange__channel__zeff_method = charge_exchange__channel__zeff_method()
    t_i_average_method :: charge_exchange__channel__t_i_average_method = charge_exchange__channel__t_i_average_method()
    zeff_line_average :: charge_exchange__channel__zeff_line_average = charge_exchange__channel__zeff_line_average()
    bes :: charge_exchange__channel__bes = charge_exchange__channel__bes()
    spectrum :: StructArray{charge_exchange__channel__spectrum} = StructArray(charge_exchange__channel__spectrum() for k in 1:1)
    momentum_tor_method :: charge_exchange__channel__momentum_tor_method = charge_exchange__channel__momentum_tor_method()
    identifier :: String = ""
    zeff :: charge_exchange__channel__zeff = charge_exchange__channel__zeff()
end

Base.@kwdef mutable struct charge_exchange
    aperture :: charge_exchange__aperture = charge_exchange__aperture()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    channel :: StructArray{charge_exchange__channel} = StructArray(charge_exchange__channel() for k in 1:1)
    ids_properties :: charge_exchange__ids_properties = charge_exchange__ids_properties()
    etendue :: Float64 = 0.0
    latency :: Float64 = 0.0
    code :: charge_exchange__code = charge_exchange__code()
    etendue_method :: charge_exchange__etendue_method = charge_exchange__etendue_method()
end

Base.@kwdef mutable struct camera_visible__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct camera_visible__channel__detector__frame
    time :: Float64 = 0.0
    radiance :: Array{Float64, 2} = zeros(Float64,(0,0))
    image_raw :: Array{Int32, 2} = zeros(Int32,(0, 0))
end

Base.@kwdef mutable struct camera_visible__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct camera_visible__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct camera_visible__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct camera_visible__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct camera_visible__code
    library :: StructArray{camera_visible__code__library} = StructArray(camera_visible__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct camera_visible__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct camera_visible__ids_properties
    provider :: String = ""
    version_put :: camera_visible__ids_properties__version_put = camera_visible__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct camera_visible__channel__detector__geometry_matrix__emission_grid
    phi :: Array{Float64, 1} = zeros(Float64,(0))
    r :: Array{Float64, 1} = zeros(Float64,(0))
    z :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct camera_visible__channel__detector__geometry_matrix
    emission_grid :: camera_visible__channel__detector__geometry_matrix__emission_grid = camera_visible__channel__detector__geometry_matrix__emission_grid()
    data :: Array{Float64, 5} = zeros(Float64,(0,0,0,0,0))
    voxel_map :: Array{Int32, 3} = zeros(Int32,(0, 0, 0))
end

Base.@kwdef mutable struct camera_visible__channel__detector
    pixel_to_alpha :: Array{Float64, 1} = zeros(Float64,(0))
    wavelength_upper :: Float64 = 0.0
    frame :: StructArray{camera_visible__channel__detector__frame} = StructArray(camera_visible__channel__detector__frame() for k in 1:1)
    exposure_time :: Float64 = 0.0
    pixel_to_beta :: Array{Float64, 1} = zeros(Float64,(0))
    counts_to_radiance :: Array{Float64, 2} = zeros(Float64,(0,0))
    wavelength_lower :: Float64 = 0.0
    noise :: Float64 = 0.0
    geometry_matrix :: camera_visible__channel__detector__geometry_matrix = camera_visible__channel__detector__geometry_matrix()
end

Base.@kwdef mutable struct camera_visible__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct camera_visible__channel__aperture
    x1_width :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: camera_visible__channel__aperture__centre = camera_visible__channel__aperture__centre()
    outline :: camera_visible__channel__aperture__outline = camera_visible__channel__aperture__outline()
    radius :: Float64 = 0.0
    x3_unit_vector :: camera_visible__channel__aperture__x3_unit_vector = camera_visible__channel__aperture__x3_unit_vector()
    surface :: Float64 = 0.0
    x2_unit_vector :: camera_visible__channel__aperture__x2_unit_vector = camera_visible__channel__aperture__x2_unit_vector()
    geometry_type :: Int32 = 0
    x1_unit_vector :: camera_visible__channel__aperture__x1_unit_vector = camera_visible__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct camera_visible__channel
    name :: String = ""
    viewing_angle_alpha_bounds :: Array{Float64, 1} = zeros(Float64,(0))
    viewing_angle_beta_bounds :: Array{Float64, 1} = zeros(Float64,(0))
    aperture :: StructArray{camera_visible__channel__aperture} = StructArray(camera_visible__channel__aperture() for k in 1:1)
    detector :: StructArray{camera_visible__channel__detector} = StructArray(camera_visible__channel__detector() for k in 1:1)
end

Base.@kwdef mutable struct camera_visible
    name :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    channel :: StructArray{camera_visible__channel} = StructArray(camera_visible__channel() for k in 1:1)
    ids_properties :: camera_visible__ids_properties = camera_visible__ids_properties()
    latency :: Float64 = 0.0
    code :: camera_visible__code = camera_visible__code()
end

Base.@kwdef mutable struct camera_ir__frame_analysis
    time :: Float64 = 0.0
    distance_separatrix_midplane :: Array{Float64, 1} = zeros(Float64,(0))
    sol_heat_decay_length :: Float64 = 0.0
    power_flux_parallel :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct camera_ir__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct camera_ir__ids_properties
    provider :: String = ""
    version_put :: camera_ir__ids_properties__version_put = camera_ir__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct camera_ir__calibration
    transmission_window :: Array{Int32, 2} = zeros(Int32,(0, 0))
    luminance_to_temperature :: Array{Int32, 2} = zeros(Int32,(0, 0))
    optical_temperature :: Array{Int32, 2} = zeros(Int32,(0, 0))
    transmission_mirror :: Array{Int32, 2} = zeros(Int32,(0, 0))
    transmission_barrel :: Array{Int32, 2} = zeros(Int32,(0, 0))
end

Base.@kwdef mutable struct camera_ir__midplane
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct camera_ir__frame
    time :: Float64 = 0.0
    image_raw :: Array{Int32, 2} = zeros(Int32,(0, 0))
end

Base.@kwdef mutable struct camera_ir__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct camera_ir__code
    library :: StructArray{camera_ir__code__library} = StructArray(camera_ir__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct camera_ir
    calibration :: camera_ir__calibration = camera_ir__calibration()
    name :: String = ""
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: camera_ir__ids_properties = camera_ir__ids_properties()
    frame_analysis :: StructArray{camera_ir__frame_analysis} = StructArray(camera_ir__frame_analysis() for k in 1:1)
    latency :: Float64 = 0.0
    midplane :: camera_ir__midplane = camera_ir__midplane()
    code :: camera_ir__code = camera_ir__code()
    frame :: StructArray{camera_ir__frame} = StructArray(camera_ir__frame() for k in 1:1)
end

Base.@kwdef mutable struct calorimetry__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct calorimetry__group__component__temperature_out
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct calorimetry__cooling_loop__mass_flow
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__ids_properties
    provider :: String = ""
    version_put :: calorimetry__ids_properties__version_put = calorimetry__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct calorimetry__group__component__transit_time
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__code
    library :: StructArray{calorimetry__code__library} = StructArray(calorimetry__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__group__component__power
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__cooling_loop__temperature_out
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__group__component__temperature_in
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__group__component__energy_total
    validity :: Int32 = 0
    data :: Float64 = 0.0
end

Base.@kwdef mutable struct calorimetry__cooling_loop__temperature_in
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__group__component__mass_flow
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__group__component__energy_cumulated
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct calorimetry__group__component
    name :: String = ""
    mass_flow :: calorimetry__group__component__mass_flow = calorimetry__group__component__mass_flow()
    energy_total :: calorimetry__group__component__energy_total = calorimetry__group__component__energy_total()
    temperature_in :: calorimetry__group__component__temperature_in = calorimetry__group__component__temperature_in()
    transit_time :: calorimetry__group__component__transit_time = calorimetry__group__component__transit_time()
    temperature_out :: calorimetry__group__component__temperature_out = calorimetry__group__component__temperature_out()
    energy_cumulated :: calorimetry__group__component__energy_cumulated = calorimetry__group__component__energy_cumulated()
    power :: calorimetry__group__component__power = calorimetry__group__component__power()
    identifier :: String = ""
end

Base.@kwdef mutable struct calorimetry__group
    name :: String = ""
    component :: StructArray{calorimetry__group__component} = StructArray(calorimetry__group__component() for k in 1:1)
    identifier :: String = ""
end

Base.@kwdef mutable struct calorimetry__cooling_loop
    name :: String = ""
    temperature_in :: calorimetry__cooling_loop__temperature_in = calorimetry__cooling_loop__temperature_in()
    mass_flow :: calorimetry__cooling_loop__mass_flow = calorimetry__cooling_loop__mass_flow()
    temperature_out :: calorimetry__cooling_loop__temperature_out = calorimetry__cooling_loop__temperature_out()
    identifier :: String = ""
end

Base.@kwdef mutable struct calorimetry
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: calorimetry__ids_properties = calorimetry__ids_properties()
    latency :: Float64 = 0.0
    cooling_loop :: StructArray{calorimetry__cooling_loop} = StructArray(calorimetry__cooling_loop() for k in 1:1)
    code :: calorimetry__code = calorimetry__code()
    group :: StructArray{calorimetry__group} = StructArray(calorimetry__group() for k in 1:1)
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__radiance_spectral
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__intensity
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct bremsstrahlung_visible__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct bremsstrahlung_visible__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bremsstrahlung_visible__ids_properties
    provider :: String = ""
    version_put :: bremsstrahlung_visible__ids_properties__version_put = bremsstrahlung_visible__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__zeff_line_average
    time :: Array{Float64, 1} = zeros(Float64,(0))
    validity :: Int32 = 0
    data :: Array{Float64, 1} = zeros(Float64,(0))
    validity_timed :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__line_of_sight
    first_point :: bremsstrahlung_visible__channel__line_of_sight__first_point = bremsstrahlung_visible__channel__line_of_sight__first_point()
    second_point :: bremsstrahlung_visible__channel__line_of_sight__second_point = bremsstrahlung_visible__channel__line_of_sight__second_point()
end

Base.@kwdef mutable struct bremsstrahlung_visible__code
    library :: StructArray{bremsstrahlung_visible__code__library} = StructArray(bremsstrahlung_visible__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel__filter
    wavelengths :: Array{Float64, 1} = zeros(Float64,(0))
    detection_efficiency :: Array{Float64, 1} = zeros(Float64,(0))
    wavelength_lower :: Float64 = 0.0
    wavelength_upper :: Float64 = 0.0
end

Base.@kwdef mutable struct bremsstrahlung_visible__channel
    intensity :: bremsstrahlung_visible__channel__intensity = bremsstrahlung_visible__channel__intensity()
    name :: String = ""
    line_of_sight :: bremsstrahlung_visible__channel__line_of_sight = bremsstrahlung_visible__channel__line_of_sight()
    filter :: bremsstrahlung_visible__channel__filter = bremsstrahlung_visible__channel__filter()
    radiance_spectral :: bremsstrahlung_visible__channel__radiance_spectral = bremsstrahlung_visible__channel__radiance_spectral()
    zeff_line_average :: bremsstrahlung_visible__channel__zeff_line_average = bremsstrahlung_visible__channel__zeff_line_average()
end

Base.@kwdef mutable struct bremsstrahlung_visible
    latency :: Float64 = 0.0
    channel :: StructArray{bremsstrahlung_visible__channel} = StructArray(bremsstrahlung_visible__channel() for k in 1:1)
    ids_properties :: bremsstrahlung_visible__ids_properties = bremsstrahlung_visible__ids_properties()
    time :: Array{Float64, 1} = zeros(Float64,(0))
    code :: bremsstrahlung_visible__code = bremsstrahlung_visible__code()
end

Base.@kwdef mutable struct bolometer__channel__validity_timed
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct bolometer__channel__line_of_sight__first_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__detector__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__detector__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct bolometer__channel__detector__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__detector__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__line_of_sight__second_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__aperture__centre
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__aperture__x1_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__detector__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__detector
    geometry_type :: Int32 = 0
    x2_width :: Float64 = 0.0
    centre :: bolometer__channel__detector__centre = bolometer__channel__detector__centre()
    outline :: bolometer__channel__detector__outline = bolometer__channel__detector__outline()
    x1_width :: Float64 = 0.0
    x2_unit_vector :: bolometer__channel__detector__x2_unit_vector = bolometer__channel__detector__x2_unit_vector()
    surface :: Float64 = 0.0
    x3_unit_vector :: bolometer__channel__detector__x3_unit_vector = bolometer__channel__detector__x3_unit_vector()
    radius :: Float64 = 0.0
    x1_unit_vector :: bolometer__channel__detector__x1_unit_vector = bolometer__channel__detector__x1_unit_vector()
end

Base.@kwdef mutable struct bolometer__channel__aperture__x2_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__line_of_sight__third_point
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__aperture__x3_unit_vector
    x :: Float64 = 0.0
    z :: Float64 = 0.0
    y :: Float64 = 0.0
end

Base.@kwdef mutable struct bolometer__channel__etendue_method
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct bolometer__channel__line_of_sight
    first_point :: bolometer__channel__line_of_sight__first_point = bolometer__channel__line_of_sight__first_point()
    second_point :: bolometer__channel__line_of_sight__second_point = bolometer__channel__line_of_sight__second_point()
    third_point :: bolometer__channel__line_of_sight__third_point = bolometer__channel__line_of_sight__third_point()
end

Base.@kwdef mutable struct bolometer__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct bolometer__ids_properties
    provider :: String = ""
    version_put :: bolometer__ids_properties__version_put = bolometer__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct bolometer__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct bolometer__code
    library :: StructArray{bolometer__code__library} = StructArray(bolometer__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
    version :: String = ""
end

Base.@kwdef mutable struct bolometer__channel__power
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct bolometer__channel__aperture__outline
    x1 :: Array{Float64, 1} = zeros(Float64,(0))
    x2 :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct bolometer__channel__aperture
    x1_width :: Float64 = 0.0
    x2_width :: Float64 = 0.0
    centre :: bolometer__channel__aperture__centre = bolometer__channel__aperture__centre()
    outline :: bolometer__channel__aperture__outline = bolometer__channel__aperture__outline()
    radius :: Float64 = 0.0
    x3_unit_vector :: bolometer__channel__aperture__x3_unit_vector = bolometer__channel__aperture__x3_unit_vector()
    x2_unit_vector :: bolometer__channel__aperture__x2_unit_vector = bolometer__channel__aperture__x2_unit_vector()
    surface :: Float64 = 0.0
    geometry_type :: Int32 = 0
    x1_unit_vector :: bolometer__channel__aperture__x1_unit_vector = bolometer__channel__aperture__x1_unit_vector()
end

Base.@kwdef mutable struct bolometer__channel
    name :: String = ""
    validity :: Int32 = 0
    aperture :: StructArray{bolometer__channel__aperture} = StructArray(bolometer__channel__aperture() for k in 1:1)
    line_of_sight :: bolometer__channel__line_of_sight = bolometer__channel__line_of_sight()
    etendue :: Float64 = 0.0
    power :: bolometer__channel__power = bolometer__channel__power()
    etendue_method :: bolometer__channel__etendue_method = bolometer__channel__etendue_method()
    identifier :: String = ""
    validity_timed :: bolometer__channel__validity_timed = bolometer__channel__validity_timed()
    detector :: bolometer__channel__detector = bolometer__channel__detector()
end

Base.@kwdef mutable struct bolometer
    time :: Array{Float64, 1} = zeros(Float64,(0))
    channel :: StructArray{bolometer__channel} = StructArray(bolometer__channel() for k in 1:1)
    ids_properties :: bolometer__ids_properties = bolometer__ids_properties()
    latency :: Float64 = 0.0
    power_radiated_validity :: Array{Int32, 1} = zeros(Int32,(0))
    power_radiated_inside_lcfs :: Array{Float64, 1} = zeros(Float64,(0))
    power_radiated_total :: Array{Float64, 1} = zeros(Float64,(0))
    code :: bolometer__code = bolometer__code()
end

Base.@kwdef mutable struct barometry__gauge__pressure
    time :: Array{Float64, 1} = zeros(Float64,(0))
    data :: Array{Float64, 1} = zeros(Float64,(0))
end

Base.@kwdef mutable struct barometry__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct barometry__ids_properties
    provider :: String = ""
    version_put :: barometry__ids_properties__version_put = barometry__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct barometry__gauge__position
    phi :: Float64 = 0.0
    r :: Float64 = 0.0
    z :: Float64 = 0.0
end

Base.@kwdef mutable struct barometry__gauge__type
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct barometry__gauge
    name :: String = ""
    pressure :: barometry__gauge__pressure = barometry__gauge__pressure()
    calibration_coefficient :: Float64 = 0.0
    position :: barometry__gauge__position = barometry__gauge__position()
    type :: barometry__gauge__type = barometry__gauge__type()
end

Base.@kwdef mutable struct barometry__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct barometry__code
    library :: StructArray{barometry__code__library} = StructArray(barometry__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct barometry
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: barometry__ids_properties = barometry__ids_properties()
    latency :: Float64 = 0.0
    gauge :: StructArray{barometry__gauge} = StructArray(barometry__gauge() for k in 1:1)
    code :: barometry__code = barometry__code()
end

Base.@kwdef mutable struct amns_data__release__data_entry
    run :: Int32 = 0
    shot :: Int32 = 0
    description :: String = ""
end

Base.@kwdef mutable struct amns_data__process__products__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct amns_data__process__charge_state
    table_6d :: Array{Float64, 6} = zeros(Float64,(0,0,0,0,0,0))
    label :: String = ""
    table_2d :: Array{Float64, 2} = zeros(Float64,(0,0))
    table_1d :: Array{Float64, 1} = zeros(Float64,(0))
    table_4d :: Array{Float64, 4} = zeros(Float64,(0,0,0,0))
    table_0d :: Float64 = 0.0
    table_3d :: Array{Float64, 3} = zeros(Float64,(0,0,0))
    z_min :: Float64 = 0.0
    table_5d :: Array{Float64, 5} = zeros(Float64,(0,0,0,0,0))
    z_max :: Float64 = 0.0
end

Base.@kwdef mutable struct amns_data__ids_properties__version_put
    access_layer_language :: String = ""
    data_dictionary :: String = ""
    access_layer :: String = ""
end

Base.@kwdef mutable struct amns_data__ids_properties
    provider :: String = ""
    version_put :: amns_data__ids_properties__version_put = amns_data__ids_properties__version_put()
    homogeneous_time :: Int32 = 0
    source :: String = ""
    creation_date :: String = ""
    comment :: String = ""
    occurrence :: Int32 = 0
end

Base.@kwdef mutable struct amns_data__code__library
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
end

Base.@kwdef mutable struct amns_data__process__products__role
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct amns_data__coordinate_system__coordinate
    units :: String = ""
    label :: String = ""
    values :: Array{Float64, 1} = zeros(Float64,(0))
    extrapolation_type :: Array{Int32, 1} = zeros(Int32,(0))
    interpolation_type :: Int32 = 0
    value_labels :: Array{String, 1} = String[]
    transformation :: Int32 = 0
    spacing :: Int32 = 0
end

Base.@kwdef mutable struct amns_data__coordinate_system
    coordinate :: StructArray{amns_data__coordinate_system__coordinate} = StructArray(amns_data__coordinate_system__coordinate() for k in 1:1)
end

Base.@kwdef mutable struct amns_data__process__products
    label :: String = ""
    role :: amns_data__process__products__role = amns_data__process__products__role()
    multiplicity :: Float64 = 0.0
    mass :: Float64 = 0.0
    metastable_label :: String = ""
    relative_charge :: Int32 = 0
    charge :: Float64 = 0.0
    metastable :: Array{Int32, 1} = zeros(Int32,(0))
    element :: StructArray{amns_data__process__products__element} = StructArray(amns_data__process__products__element() for k in 1:1)
end

Base.@kwdef mutable struct amns_data__process__reactants__element
    atoms_n :: Int32 = 0
    z_n :: Float64 = 0.0
    multiplicity :: Float64 = 0.0
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct amns_data__release
    date :: String = ""
    description :: String = ""
    data_entry :: StructArray{amns_data__release__data_entry} = StructArray(amns_data__release__data_entry() for k in 1:1)
end

Base.@kwdef mutable struct amns_data__code
    library :: StructArray{amns_data__code__library} = StructArray(amns_data__code__library() for k in 1:1)
    name :: String = ""
    parameters :: String = ""
    commit :: String = ""
    repository :: String = ""
    version :: String = ""
    output_flag :: Array{Int32, 1} = zeros(Int32,(0))
end

Base.@kwdef mutable struct amns_data__process__reactants__role
    name :: String = ""
    description :: String = ""
    index :: Int32 = 0
end

Base.@kwdef mutable struct amns_data__process__reactants
    label :: String = ""
    role :: amns_data__process__reactants__role = amns_data__process__reactants__role()
    multiplicity :: Float64 = 0.0
    mass :: Float64 = 0.0
    metastable_label :: String = ""
    charge :: Float64 = 0.0
    relative_charge :: Int32 = 0
    metastable :: Array{Int32, 1} = zeros(Int32,(0))
    element :: StructArray{amns_data__process__reactants__element} = StructArray(amns_data__process__reactants__element() for k in 1:1)
end

Base.@kwdef mutable struct amns_data__process
    label :: String = ""
    result_label :: String = ""
    charge_state :: StructArray{amns_data__process__charge_state} = StructArray(amns_data__process__charge_state() for k in 1:1)
    products :: StructArray{amns_data__process__products} = StructArray(amns_data__process__products() for k in 1:1)
    result_units :: String = ""
    coordinate_index :: Int32 = 0
    citation :: String = ""
    result_transformation :: Int32 = 0
    provider :: String = ""
    source :: String = ""
    table_dimension :: Int32 = 0
    reactants :: StructArray{amns_data__process__reactants} = StructArray(amns_data__process__reactants() for k in 1:1)
end

Base.@kwdef mutable struct amns_data
    z_n :: Float64 = 0.0
    time :: Array{Float64, 1} = zeros(Float64,(0))
    ids_properties :: amns_data__ids_properties = amns_data__ids_properties()
    process :: StructArray{amns_data__process} = StructArray(amns_data__process() for k in 1:1)
    code :: amns_data__code = amns_data__code()
    release :: StructArray{amns_data__release} = StructArray(amns_data__release() for k in 1:1)
    coordinate_system :: StructArray{amns_data__coordinate_system} = StructArray(amns_data__coordinate_system() for k in 1:1)
    a :: Float64 = 0.0
end

Base.@kwdef mutable struct dd
    gyrokinetics :: gyrokinetics = gyrokinetics()
    numerics :: numerics = numerics()
    edge_transport :: edge_transport = edge_transport()
    hard_x_rays :: hard_x_rays = hard_x_rays()
    pulse_schedule :: pulse_schedule = pulse_schedule()
    waves :: waves = waves()
    core_sources :: core_sources = core_sources()
    coils_non_axisymmetric :: coils_non_axisymmetric = coils_non_axisymmetric()
    interferometer :: interferometer = interferometer()
    spectrometer_uv :: spectrometer_uv = spectrometer_uv()
    dataset_description :: dataset_description = dataset_description()
    ece :: ece = ece()
    refractometer :: refractometer = refractometer()
    edge_sources :: edge_sources = edge_sources()
    reflectometer_profile :: reflectometer_profile = reflectometer_profile()
    edge_profiles :: edge_profiles = edge_profiles()
    core_instant_changes :: core_instant_changes = core_instant_changes()
    pellets :: pellets = pellets()
    neutron_diagnostic :: neutron_diagnostic = neutron_diagnostic()
    mhd :: mhd = mhd()
    em_coupling :: em_coupling = em_coupling()
    pf_passive :: pf_passive = pf_passive()
    mhd_linear :: mhd_linear = mhd_linear()
    camera_visible :: camera_visible = camera_visible()
    sawteeth :: sawteeth = sawteeth()
    polarimeter :: polarimeter = polarimeter()
    gas_pumping :: gas_pumping = gas_pumping()
    iron_core :: iron_core = iron_core()
    disruption :: disruption = disruption()
    ic_antennas :: ic_antennas = ic_antennas()
    core_transport :: core_transport = core_transport()
    amns_data :: amns_data = amns_data()
    spectrometer_x_ray_crystal :: spectrometer_x_ray_crystal = spectrometer_x_ray_crystal()
    temporary :: temporary = temporary()
    camera_ir :: camera_ir = camera_ir()
    charge_exchange :: charge_exchange = charge_exchange()
    thomson_scattering :: thomson_scattering = thomson_scattering()
    nbi :: nbi = nbi()
    divertors :: divertors = divertors()
    ntms :: ntms = ntms()
    pf_active :: pf_active = pf_active()
    controllers :: controllers = controllers()
    turbulence :: turbulence = turbulence()
    soft_x_rays :: soft_x_rays = soft_x_rays()
    barometry :: barometry = barometry()
    distributions :: distributions = distributions()
    tf :: tf = tf()
    magnetics :: magnetics = magnetics()
    equilibrium :: equilibrium = equilibrium()
    spectrometer_mass :: spectrometer_mass = spectrometer_mass()
    spectrometer_visible :: spectrometer_visible = spectrometer_visible()
    langmuir_probes :: langmuir_probes = langmuir_probes()
    core_profiles :: core_profiles = core_profiles()
    ec_launchers :: ec_launchers = ec_launchers()
    radiation :: radiation = radiation()
    transport_solver_numerics :: transport_solver_numerics = transport_solver_numerics()
    cryostat :: cryostat = cryostat()
    bremsstrahlung_visible :: bremsstrahlung_visible = bremsstrahlung_visible()
    distribution_sources :: distribution_sources = distribution_sources()
    bolometer :: bolometer = bolometer()
    dataset_fair :: dataset_fair = dataset_fair()
    calorimetry :: calorimetry = calorimetry()
    lh_antennas :: lh_antennas = lh_antennas()
    summary :: summary = summary()
    gas_injection :: gas_injection = gas_injection()
    mse :: mse = mse()
    sdn :: sdn = sdn()
    wall :: wall = wall()
end

