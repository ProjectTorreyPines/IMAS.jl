include("functionarrays.jl")

conversion_types = Union{Real, String, Array{Integer, N} where N, Array{Real, N} where N, Array{String, N} where N}

Base.@kwdef mutable struct wall__temperature_reference <: FDS
    var"data" :: Union{Missing, Real, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__temperature_reference(var"data"=missing, var"description"=missing, _parent=WeakRef(missing))
        fds = new(var"data", var"description", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String, Function} = missing
    var"data_dictionary" :: Union{Missing, String, Function} = missing
    var"access_layer" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__ids_properties <: FDS
    var"provider" :: Union{Missing, String, Function} = missing
    var"version_put" :: wall__ids_properties__version_put = wall__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Integer, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"creation_date" :: Union{Missing, String, Function} = missing
    var"comment" :: Union{Missing, String, Function} = missing
    var"occurrence" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__ids_properties(var"provider"=missing, var"version_put"=wall__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral___element <: FDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function} = missing
    var"z_n" :: Union{Missing, Real, Function} = missing
    var"multiplicity" :: Union{Missing, Real, Function} = missing
    var"a" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities__neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        fds = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral <: FDSvectorElement
    var"label" :: Union{Missing, String, Function} = missing
    var"sputtering_chemical_coefficient" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"gas_puff" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"recycling_particles_coefficient" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"pumping_speed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"particle_flux_from_wall" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"recycling_energy_coefficient" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"wall_inventory" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"particle_flux_from_plasma" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"sputtering_physical_coefficient" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"element" :: FDSvector{T} where {T<:wall__global_quantities__neutral___element} = FDSvector(wall__global_quantities__neutral___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities__neutral(var"label"=missing, var"sputtering_chemical_coefficient"=missing, var"gas_puff"=missing, var"recycling_particles_coefficient"=missing, var"pumping_speed"=missing, var"particle_flux_from_wall"=missing, var"recycling_energy_coefficient"=missing, var"wall_inventory"=missing, var"particle_flux_from_plasma"=missing, var"sputtering_physical_coefficient"=missing, var"element"=FDSvector(wall__global_quantities__neutral___element[]), _parent=WeakRef(missing))
        fds = new(var"label", var"sputtering_chemical_coefficient", var"gas_puff", var"recycling_particles_coefficient", var"pumping_speed", var"particle_flux_from_wall", var"recycling_energy_coefficient", var"wall_inventory", var"particle_flux_from_plasma", var"sputtering_physical_coefficient", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__global_quantities__electrons <: FDS
    var"particle_flux_from_plasma" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gas_puff" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_outer_target" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pumping_speed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"particle_flux_from_wall" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"power_inner_target" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities__electrons(var"particle_flux_from_plasma"=missing, var"gas_puff"=missing, var"power_outer_target"=missing, var"pumping_speed"=missing, var"particle_flux_from_wall"=missing, var"power_inner_target"=missing, _parent=WeakRef(missing))
        fds = new(var"particle_flux_from_plasma", var"gas_puff", var"power_outer_target", var"pumping_speed", var"particle_flux_from_wall", var"power_inner_target", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__global_quantities <: FDS
    var"neutral" :: FDSvector{T} where {T<:wall__global_quantities__neutral} = FDSvector(wall__global_quantities__neutral[])
    var"power_incident" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_radiated" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_inner_target_ion_total" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"temperature" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_conducted" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_convected" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"current_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"electrons" :: wall__global_quantities__electrons = wall__global_quantities__electrons()
    var"power_density_inner_target_max" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_black_body" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_recombination_neutrals" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_to_cooling" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_density_outer_target_max" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_recombination_plasma" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_currents" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"power_neutrals" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities(var"neutral"=FDSvector(wall__global_quantities__neutral[]), var"power_incident"=missing, var"power_radiated"=missing, var"power_inner_target_ion_total"=missing, var"temperature"=missing, var"power_conducted"=missing, var"power_convected"=missing, var"current_tor"=missing, var"electrons"=wall__global_quantities__electrons(), var"power_density_inner_target_max"=missing, var"power_black_body"=missing, var"power_recombination_neutrals"=missing, var"power_to_cooling"=missing, var"power_density_outer_target_max"=missing, var"power_recombination_plasma"=missing, var"power_currents"=missing, var"power_neutrals"=missing, _parent=WeakRef(missing))
        fds = new(var"neutral", var"power_incident", var"power_radiated", var"power_inner_target_ion_total", var"temperature", var"power_conducted", var"power_convected", var"current_tor", var"electrons", var"power_density_inner_target_max", var"power_black_body", var"power_recombination_neutrals", var"power_to_cooling", var"power_density_outer_target_max", var"power_recombination_plasma", var"power_currents", var"power_neutrals", _parent)
        assign_expressions(fds)
        setfield!(fds.neutral, :_parent, WeakRef(fds))
        setfield!(fds.electrons, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__first_wall_power_flux_peak <: FDS
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"data" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__first_wall_power_flux_peak(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        fds = new(var"time", var"data", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary <: FDSvectorElement
    var"neighbours" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"neighbours", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object <: FDSvectorElement
    var"nodes" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"measure" :: Union{Missing, Real, Function} = missing
    var"geometry" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"boundary" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        fds = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_expressions(fds)
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension(var"object"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_expressions(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___geometry_type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space <: FDSvectorElement
    var"coordinates_type" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"geometry_type" :: wall__description_ggd___grid_ggd___space___geometry_type = wall__description_ggd___grid_ggd___space___geometry_type()
    var"identifier" :: wall__description_ggd___grid_ggd___space___identifier = wall__description_ggd___grid_ggd___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space(var"coordinates_type"=missing, var"geometry_type"=wall__description_ggd___grid_ggd___space___geometry_type(), var"identifier"=wall__description_ggd___grid_ggd___space___identifier(), var"objects_per_dimension"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[]), _parent=WeakRef(missing))
        fds = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(fds)
        setfield!(fds.geometry_type, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.objects_per_dimension, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element___object <: FDSvectorElement
    var"dimension" :: Union{Missing, Integer, Function} = missing
    var"space" :: Union{Missing, Integer, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"dimension", var"space", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element___object} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___element(var"object"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_expressions(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___base <: FDSvectorElement
    var"jacobian" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset <: FDSvectorElement
    var"base" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___base} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[])
    var"metric" :: wall__description_ggd___grid_ggd___grid_subset___metric = wall__description_ggd___grid_ggd___grid_subset___metric()
    var"dimension" :: Union{Missing, Integer, Function} = missing
    var"identifier" :: wall__description_ggd___grid_ggd___grid_subset___identifier = wall__description_ggd___grid_ggd___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset(var"base"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[]), var"metric"=wall__description_ggd___grid_ggd___grid_subset___metric(), var"dimension"=missing, var"identifier"=wall__description_ggd___grid_ggd___grid_subset___identifier(), var"element"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[]), _parent=WeakRef(missing))
        fds = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.base, :_parent, WeakRef(fds))
        setfield!(fds.metric, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd <: FDSvectorElement
    var"time" :: Union{Missing, Real, Function} = missing
    var"grid_subset" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset} = FDSvector(wall__description_ggd___grid_ggd___grid_subset[])
    var"space" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space} = FDSvector(wall__description_ggd___grid_ggd___space[])
    var"identifier" :: wall__description_ggd___grid_ggd___identifier = wall__description_ggd___grid_ggd___identifier()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd(var"time"=missing, var"grid_subset"=FDSvector(wall__description_ggd___grid_ggd___grid_subset[]), var"space"=FDSvector(wall__description_ggd___grid_ggd___space[]), var"identifier"=wall__description_ggd___grid_ggd___identifier(), _parent=WeakRef(missing))
        fds = new(var"time", var"grid_subset", var"space", var"identifier", _parent)
        assign_expressions(fds)
        setfield!(fds.grid_subset, :_parent, WeakRef(fds))
        setfield!(fds.space, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___temperature <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd___temperature(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___power_density <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd___power_density(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd <: FDSvectorElement
    var"temperature" :: FDSvector{T} where {T<:wall__description_ggd___ggd___temperature} = FDSvector(wall__description_ggd___ggd___temperature[])
    var"time" :: Union{Missing, Real, Function} = missing
    var"power_density" :: FDSvector{T} where {T<:wall__description_ggd___ggd___power_density} = FDSvector(wall__description_ggd___ggd___power_density[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd(var"temperature"=FDSvector(wall__description_ggd___ggd___temperature[]), var"time"=missing, var"power_density"=FDSvector(wall__description_ggd___ggd___power_density[]), _parent=WeakRef(missing))
        fds = new(var"temperature", var"time", var"power_density", _parent)
        assign_expressions(fds)
        setfield!(fds.temperature, :_parent, WeakRef(fds))
        setfield!(fds.power_density, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd <: FDSvectorElement
    var"grid_ggd" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd} = FDSvector(wall__description_ggd___grid_ggd[])
    var"type" :: wall__description_ggd___type = wall__description_ggd___type()
    var"ggd" :: FDSvector{T} where {T<:wall__description_ggd___ggd} = FDSvector(wall__description_ggd___ggd[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd(var"grid_ggd"=FDSvector(wall__description_ggd___grid_ggd[]), var"type"=wall__description_ggd___type(), var"ggd"=FDSvector(wall__description_ggd___ggd[]), _parent=WeakRef(missing))
        fds = new(var"grid_ggd", var"type", var"ggd", _parent)
        assign_expressions(fds)
        setfield!(fds.grid_ggd, :_parent, WeakRef(fds))
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.ggd, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___outline <: FDS
    var"closed" :: Union{Missing, Integer, Function} = missing
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element___outline(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___j_tor <: FDS
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"data" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element___j_tor(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        fds = new(var"time", var"data", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"j_tor" :: wall__description_2d___vessel__unit___element___j_tor = wall__description_2d___vessel__unit___element___j_tor()
    var"resistivity" :: Union{Missing, Real, Function} = missing
    var"outline" :: wall__description_2d___vessel__unit___element___outline = wall__description_2d___vessel__unit___element___outline()
    var"resistance" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element(var"name"=missing, var"j_tor"=wall__description_2d___vessel__unit___element___j_tor(), var"resistivity"=missing, var"outline"=wall__description_2d___vessel__unit___element___outline(), var"resistance"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"j_tor", var"resistivity", var"outline", var"resistance", _parent)
        assign_expressions(fds)
        setfield!(fds.j_tor, :_parent, WeakRef(fds))
        setfield!(fds.outline, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_outer <: FDS
    var"closed" :: Union{Missing, Integer, Function} = missing
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular__outline_outer(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_inner <: FDS
    var"closed" :: Union{Missing, Integer, Function} = missing
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular__outline_inner(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__centreline <: FDS
    var"closed" :: Union{Missing, Integer, Function} = missing
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular__centreline(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular <: FDS
    var"outline_inner" :: wall__description_2d___vessel__unit___annular__outline_inner = wall__description_2d___vessel__unit___annular__outline_inner()
    var"centreline" :: wall__description_2d___vessel__unit___annular__centreline = wall__description_2d___vessel__unit___annular__centreline()
    var"thickness" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"resistivity" :: Union{Missing, Real, Function} = missing
    var"outline_outer" :: wall__description_2d___vessel__unit___annular__outline_outer = wall__description_2d___vessel__unit___annular__outline_outer()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular(var"outline_inner"=wall__description_2d___vessel__unit___annular__outline_inner(), var"centreline"=wall__description_2d___vessel__unit___annular__centreline(), var"thickness"=missing, var"resistivity"=missing, var"outline_outer"=wall__description_2d___vessel__unit___annular__outline_outer(), _parent=WeakRef(missing))
        fds = new(var"outline_inner", var"centreline", var"thickness", var"resistivity", var"outline_outer", _parent)
        assign_expressions(fds)
        setfield!(fds.outline_inner, :_parent, WeakRef(fds))
        setfield!(fds.centreline, :_parent, WeakRef(fds))
        setfield!(fds.outline_outer, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"annular" :: wall__description_2d___vessel__unit___annular = wall__description_2d___vessel__unit___annular()
    var"identifier" :: Union{Missing, String, Function} = missing
    var"element" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit___element} = FDSvector(wall__description_2d___vessel__unit___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit(var"name"=missing, var"annular"=wall__description_2d___vessel__unit___annular(), var"identifier"=missing, var"element"=FDSvector(wall__description_2d___vessel__unit___element[]), _parent=WeakRef(missing))
        fds = new(var"name", var"annular", var"identifier", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.annular, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel <: FDS
    var"type" :: wall__description_2d___vessel__type = wall__description_2d___vessel__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit} = FDSvector(wall__description_2d___vessel__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel(var"type"=wall__description_2d___vessel__type(), var"unit"=FDSvector(wall__description_2d___vessel__unit[]), _parent=WeakRef(missing))
        fds = new(var"type", var"unit", _parent)
        assign_expressions(fds)
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.unit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit___outline <: FDSvectorElement
    var"time" :: Union{Missing, Real, Function} = missing
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__unit___outline(var"time"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"time", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"resistivity" :: Union{Missing, Real, Function} = missing
    var"closed" :: Union{Missing, Integer, Function} = missing
    var"phi_extensions" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"outline" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit___outline} = FDSvector(wall__description_2d___mobile__unit___outline[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"phi_extensions"=missing, var"outline"=FDSvector(wall__description_2d___mobile__unit___outline[]), _parent=WeakRef(missing))
        fds = new(var"name", var"resistivity", var"closed", var"phi_extensions", var"outline", _parent)
        assign_expressions(fds)
        setfield!(fds.outline, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile <: FDS
    var"type" :: wall__description_2d___mobile__type = wall__description_2d___mobile__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit} = FDSvector(wall__description_2d___mobile__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile(var"type"=wall__description_2d___mobile__type(), var"unit"=FDSvector(wall__description_2d___mobile__unit[]), _parent=WeakRef(missing))
        fds = new(var"type", var"unit", _parent)
        assign_expressions(fds)
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.unit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit___outline <: FDS
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__unit___outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"resistivity" :: Union{Missing, Real, Function} = missing
    var"closed" :: Union{Missing, Integer, Function} = missing
    var"outline" :: wall__description_2d___limiter__unit___outline = wall__description_2d___limiter__unit___outline()
    var"phi_extensions" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"outline"=wall__description_2d___limiter__unit___outline(), var"phi_extensions"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"resistivity", var"closed", var"outline", var"phi_extensions", _parent)
        assign_expressions(fds)
        setfield!(fds.outline, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter <: FDS
    var"type" :: wall__description_2d___limiter__type = wall__description_2d___limiter__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___limiter__unit} = FDSvector(wall__description_2d___limiter__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter(var"type"=wall__description_2d___limiter__type(), var"unit"=FDSvector(wall__description_2d___limiter__unit[]), _parent=WeakRef(missing))
        fds = new(var"type", var"unit", _parent)
        assign_expressions(fds)
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.unit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d <: FDSvectorElement
    var"mobile" :: wall__description_2d___mobile = wall__description_2d___mobile()
    var"limiter" :: wall__description_2d___limiter = wall__description_2d___limiter()
    var"type" :: wall__description_2d___type = wall__description_2d___type()
    var"vessel" :: wall__description_2d___vessel = wall__description_2d___vessel()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d(var"mobile"=wall__description_2d___mobile(), var"limiter"=wall__description_2d___limiter(), var"type"=wall__description_2d___type(), var"vessel"=wall__description_2d___vessel(), _parent=WeakRef(missing))
        fds = new(var"mobile", var"limiter", var"type", var"vessel", _parent)
        assign_expressions(fds)
        setfield!(fds.mobile, :_parent, WeakRef(fds))
        setfield!(fds.limiter, :_parent, WeakRef(fds))
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.vessel, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__code <: FDS
    var"library" :: FDSvector{T} where {T<:wall__code__library} = FDSvector(wall__code__library[])
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"output_flag" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__code(var"library"=FDSvector(wall__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall <: FDS
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ids_properties" :: wall__ids_properties = wall__ids_properties()
    var"description_ggd" :: FDSvector{T} where {T<:wall__description_ggd} = FDSvector(wall__description_ggd[])
    var"description_2d" :: FDSvector{T} where {T<:wall__description_2d} = FDSvector(wall__description_2d[])
    var"first_wall_surface_area" :: Union{Missing, Real, Function} = missing
    var"code" :: wall__code = wall__code()
    var"global_quantities" :: wall__global_quantities = wall__global_quantities()
    var"temperature_reference" :: wall__temperature_reference = wall__temperature_reference()
    var"first_wall_power_flux_peak" :: wall__first_wall_power_flux_peak = wall__first_wall_power_flux_peak()
    _parent :: WeakRef = WeakRef(missing)
    function wall(var"time"=missing, var"ids_properties"=wall__ids_properties(), var"description_ggd"=FDSvector(wall__description_ggd[]), var"description_2d"=FDSvector(wall__description_2d[]), var"first_wall_surface_area"=missing, var"code"=wall__code(), var"global_quantities"=wall__global_quantities(), var"temperature_reference"=wall__temperature_reference(), var"first_wall_power_flux_peak"=wall__first_wall_power_flux_peak(), _parent=WeakRef(missing))
        fds = new(var"time", var"ids_properties", var"description_ggd", var"description_2d", var"first_wall_surface_area", var"code", var"global_quantities", var"temperature_reference", var"first_wall_power_flux_peak", _parent)
        assign_expressions(fds)
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.description_ggd, :_parent, WeakRef(fds))
        setfield!(fds.description_2d, :_parent, WeakRef(fds))
        setfield!(fds.code, :_parent, WeakRef(fds))
        setfield!(fds.global_quantities, :_parent, WeakRef(fds))
        setfield!(fds.temperature_reference, :_parent, WeakRef(fds))
        setfield!(fds.first_wall_power_flux_peak, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__wall__material <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__wall__material(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__wall__evaporation <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__wall__evaporation(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__wall <: FDS
    var"material" :: summary__wall__material = summary__wall__material()
    var"evaporation" :: summary__wall__evaporation = summary__wall__evaporation()
    _parent :: WeakRef = WeakRef(missing)
    function summary__wall(var"material"=summary__wall__material(), var"evaporation"=summary__wall__evaporation(), _parent=WeakRef(missing))
        fds = new(var"material", var"evaporation", _parent)
        assign_expressions(fds)
        setfield!(fds.material, :_parent, WeakRef(fds))
        setfield!(fds.evaporation, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_i <: FDS
    var"krypton" :: summary__volume_average__n_i__krypton = summary__volume_average__n_i__krypton()
    var"lithium" :: summary__volume_average__n_i__lithium = summary__volume_average__n_i__lithium()
    var"neon" :: summary__volume_average__n_i__neon = summary__volume_average__n_i__neon()
    var"tritium" :: summary__volume_average__n_i__tritium = summary__volume_average__n_i__tritium()
    var"helium_3" :: summary__volume_average__n_i__helium_3 = summary__volume_average__n_i__helium_3()
    var"deuterium" :: summary__volume_average__n_i__deuterium = summary__volume_average__n_i__deuterium()
    var"iron" :: summary__volume_average__n_i__iron = summary__volume_average__n_i__iron()
    var"helium_4" :: summary__volume_average__n_i__helium_4 = summary__volume_average__n_i__helium_4()
    var"oxygen" :: summary__volume_average__n_i__oxygen = summary__volume_average__n_i__oxygen()
    var"tungsten" :: summary__volume_average__n_i__tungsten = summary__volume_average__n_i__tungsten()
    var"xenon" :: summary__volume_average__n_i__xenon = summary__volume_average__n_i__xenon()
    var"hydrogen" :: summary__volume_average__n_i__hydrogen = summary__volume_average__n_i__hydrogen()
    var"carbon" :: summary__volume_average__n_i__carbon = summary__volume_average__n_i__carbon()
    var"nitrogen" :: summary__volume_average__n_i__nitrogen = summary__volume_average__n_i__nitrogen()
    var"beryllium" :: summary__volume_average__n_i__beryllium = summary__volume_average__n_i__beryllium()
    var"argon" :: summary__volume_average__n_i__argon = summary__volume_average__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_i(var"krypton"=summary__volume_average__n_i__krypton(), var"lithium"=summary__volume_average__n_i__lithium(), var"neon"=summary__volume_average__n_i__neon(), var"tritium"=summary__volume_average__n_i__tritium(), var"helium_3"=summary__volume_average__n_i__helium_3(), var"deuterium"=summary__volume_average__n_i__deuterium(), var"iron"=summary__volume_average__n_i__iron(), var"helium_4"=summary__volume_average__n_i__helium_4(), var"oxygen"=summary__volume_average__n_i__oxygen(), var"tungsten"=summary__volume_average__n_i__tungsten(), var"xenon"=summary__volume_average__n_i__xenon(), var"hydrogen"=summary__volume_average__n_i__hydrogen(), var"carbon"=summary__volume_average__n_i__carbon(), var"nitrogen"=summary__volume_average__n_i__nitrogen(), var"beryllium"=summary__volume_average__n_i__beryllium(), var"argon"=summary__volume_average__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__meff_hydrogenic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__meff_hydrogenic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__isotope_fraction_hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__isotope_fraction_hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average__dn_e_dt <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average__dn_e_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__volume_average <: FDS
    var"meff_hydrogenic" :: summary__volume_average__meff_hydrogenic = summary__volume_average__meff_hydrogenic()
    var"n_i" :: summary__volume_average__n_i = summary__volume_average__n_i()
    var"t_e" :: summary__volume_average__t_e = summary__volume_average__t_e()
    var"n_e" :: summary__volume_average__n_e = summary__volume_average__n_e()
    var"isotope_fraction_hydrogen" :: summary__volume_average__isotope_fraction_hydrogen = summary__volume_average__isotope_fraction_hydrogen()
    var"t_i_average" :: summary__volume_average__t_i_average = summary__volume_average__t_i_average()
    var"n_i_total" :: summary__volume_average__n_i_total = summary__volume_average__n_i_total()
    var"dn_e_dt" :: summary__volume_average__dn_e_dt = summary__volume_average__dn_e_dt()
    var"zeff" :: summary__volume_average__zeff = summary__volume_average__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__volume_average(var"meff_hydrogenic"=summary__volume_average__meff_hydrogenic(), var"n_i"=summary__volume_average__n_i(), var"t_e"=summary__volume_average__t_e(), var"n_e"=summary__volume_average__n_e(), var"isotope_fraction_hydrogen"=summary__volume_average__isotope_fraction_hydrogen(), var"t_i_average"=summary__volume_average__t_i_average(), var"n_i_total"=summary__volume_average__n_i_total(), var"dn_e_dt"=summary__volume_average__dn_e_dt(), var"zeff"=summary__volume_average__zeff(), _parent=WeakRef(missing))
        fds = new(var"meff_hydrogenic", var"n_i", var"t_e", var"n_e", var"isotope_fraction_hydrogen", var"t_i_average", var"n_i_total", var"dn_e_dt", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.meff_hydrogenic, :_parent, WeakRef(fds))
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.isotope_fraction_hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.dn_e_dt, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__tag <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"comment" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__tag(var"name"=missing, var"comment"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"comment", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__stationary_phase_flag <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__stationary_phase_flag(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__t_i_average_decay_length <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__t_i_average_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__t_e_decay_length <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__t_e_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__pressure_neutral <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__pressure_neutral(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__power_radiated <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__power_radiated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__n_i_total_decay_length <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__n_i_total_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__n_e_decay_length <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__n_e_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__heat_flux_i_decay_length <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__heat_flux_i_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer__heat_flux_e_decay_length <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer__heat_flux_e_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__scrape_off_layer <: FDS
    var"heat_flux_i_decay_length" :: summary__scrape_off_layer__heat_flux_i_decay_length = summary__scrape_off_layer__heat_flux_i_decay_length()
    var"n_e_decay_length" :: summary__scrape_off_layer__n_e_decay_length = summary__scrape_off_layer__n_e_decay_length()
    var"pressure_neutral" :: summary__scrape_off_layer__pressure_neutral = summary__scrape_off_layer__pressure_neutral()
    var"heat_flux_e_decay_length" :: summary__scrape_off_layer__heat_flux_e_decay_length = summary__scrape_off_layer__heat_flux_e_decay_length()
    var"power_radiated" :: summary__scrape_off_layer__power_radiated = summary__scrape_off_layer__power_radiated()
    var"n_i_total_decay_length" :: summary__scrape_off_layer__n_i_total_decay_length = summary__scrape_off_layer__n_i_total_decay_length()
    var"t_e_decay_length" :: summary__scrape_off_layer__t_e_decay_length = summary__scrape_off_layer__t_e_decay_length()
    var"t_i_average_decay_length" :: summary__scrape_off_layer__t_i_average_decay_length = summary__scrape_off_layer__t_i_average_decay_length()
    _parent :: WeakRef = WeakRef(missing)
    function summary__scrape_off_layer(var"heat_flux_i_decay_length"=summary__scrape_off_layer__heat_flux_i_decay_length(), var"n_e_decay_length"=summary__scrape_off_layer__n_e_decay_length(), var"pressure_neutral"=summary__scrape_off_layer__pressure_neutral(), var"heat_flux_e_decay_length"=summary__scrape_off_layer__heat_flux_e_decay_length(), var"power_radiated"=summary__scrape_off_layer__power_radiated(), var"n_i_total_decay_length"=summary__scrape_off_layer__n_i_total_decay_length(), var"t_e_decay_length"=summary__scrape_off_layer__t_e_decay_length(), var"t_i_average_decay_length"=summary__scrape_off_layer__t_i_average_decay_length(), _parent=WeakRef(missing))
        fds = new(var"heat_flux_i_decay_length", var"n_e_decay_length", var"pressure_neutral", var"heat_flux_e_decay_length", var"power_radiated", var"n_i_total_decay_length", var"t_e_decay_length", var"t_i_average_decay_length", _parent)
        assign_expressions(fds)
        setfield!(fds.heat_flux_i_decay_length, :_parent, WeakRef(fds))
        setfield!(fds.n_e_decay_length, :_parent, WeakRef(fds))
        setfield!(fds.pressure_neutral, :_parent, WeakRef(fds))
        setfield!(fds.heat_flux_e_decay_length, :_parent, WeakRef(fds))
        setfield!(fds.power_radiated, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total_decay_length, :_parent, WeakRef(fds))
        setfield!(fds.t_e_decay_length, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average_decay_length, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__runaways__particles <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__runaways__particles(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__runaways__current <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__runaways__current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__runaways <: FDS
    var"particles" :: summary__runaways__particles = summary__runaways__particles()
    var"current" :: summary__runaways__current = summary__runaways__current()
    _parent :: WeakRef = WeakRef(missing)
    function summary__runaways(var"particles"=summary__runaways__particles(), var"current"=summary__runaways__current(), _parent=WeakRef(missing))
        fds = new(var"particles", var"current", _parent)
        assign_expressions(fds)
        setfield!(fds.particles, :_parent, WeakRef(fds))
        setfield!(fds.current, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__rmps__occurrence <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__rmps__occurrence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__rmps <: FDS
    var"occurrence" :: summary__rmps__occurrence = summary__rmps__occurrence()
    _parent :: WeakRef = WeakRef(missing)
    function summary__rmps(var"occurrence"=summary__rmps__occurrence(), _parent=WeakRef(missing))
        fds = new(var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.occurrence, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pellets__occurrence <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pellets__occurrence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pellets <: FDS
    var"occurrence" :: summary__pellets__occurrence = summary__pellets__occurrence()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pellets(var"occurrence"=summary__pellets__occurrence(), _parent=WeakRef(missing))
        fds = new(var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.occurrence, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__volume_inside_pedestal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__volume_inside_pedestal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_width <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_height <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__offset <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__t_e <: FDS
    var"pedestal_width" :: summary__pedestal_fits__mtanh__t_e__pedestal_width = summary__pedestal_fits__mtanh__t_e__pedestal_width()
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position = summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position()
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max = summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max()
    var"offset" :: summary__pedestal_fits__mtanh__t_e__offset = summary__pedestal_fits__mtanh__t_e__offset()
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm = summary__pedestal_fits__mtanh__t_e__d_dpsi_norm()
    var"pedestal_position" :: summary__pedestal_fits__mtanh__t_e__pedestal_position = summary__pedestal_fits__mtanh__t_e__pedestal_position()
    var"pedestal_height" :: summary__pedestal_fits__mtanh__t_e__pedestal_height = summary__pedestal_fits__mtanh__t_e__pedestal_height()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__t_e(var"pedestal_width"=summary__pedestal_fits__mtanh__t_e__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__mtanh__t_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm(), var"pedestal_position"=summary__pedestal_fits__mtanh__t_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__mtanh__t_e__pedestal_height(), _parent=WeakRef(missing))
        fds = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(fds)
        setfield!(fds.pedestal_width, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max_position, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max, :_parent, WeakRef(fds))
        setfield!(fds.offset, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_position, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_height, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter <: FDS
    var"alpha_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical()
    var"t_e_pedestal_top_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical()
    var"alpha_ratio" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter(var"alpha_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical(), var"t_e_pedestal_top_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical(), var"alpha_ratio"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio(), _parent=WeakRef(missing))
        fds = new(var"alpha_critical", var"t_e_pedestal_top_critical", var"alpha_ratio", _parent)
        assign_expressions(fds)
        setfield!(fds.alpha_critical, :_parent, WeakRef(fds))
        setfield!(fds.t_e_pedestal_top_critical, :_parent, WeakRef(fds))
        setfield!(fds.alpha_ratio, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager <: FDS
    var"alpha_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical()
    var"t_e_pedestal_top_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical()
    var"alpha_ratio" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager(var"alpha_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical(), var"t_e_pedestal_top_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical(), var"alpha_ratio"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio(), _parent=WeakRef(missing))
        fds = new(var"alpha_critical", var"t_e_pedestal_top_critical", var"alpha_ratio", _parent)
        assign_expressions(fds)
        setfield!(fds.alpha_critical, :_parent, WeakRef(fds))
        setfield!(fds.t_e_pedestal_top_critical, :_parent, WeakRef(fds))
        setfield!(fds.alpha_ratio, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability__alpha_experimental <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability__alpha_experimental(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__stability <: FDS
    var"alpha_experimental" :: summary__pedestal_fits__mtanh__stability__alpha_experimental = summary__pedestal_fits__mtanh__stability__alpha_experimental()
    var"bootstrap_current_hager" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager = summary__pedestal_fits__mtanh__stability__bootstrap_current_hager()
    var"bootstrap_current_sauter" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter = summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__stability(var"alpha_experimental"=summary__pedestal_fits__mtanh__stability__alpha_experimental(), var"bootstrap_current_hager"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager(), var"bootstrap_current_sauter"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter(), _parent=WeakRef(missing))
        fds = new(var"alpha_experimental", var"bootstrap_current_hager", var"bootstrap_current_sauter", _parent)
        assign_expressions(fds)
        setfield!(fds.alpha_experimental, :_parent, WeakRef(fds))
        setfield!(fds.bootstrap_current_hager, :_parent, WeakRef(fds))
        setfield!(fds.bootstrap_current_sauter, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__separatrix <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_width <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_height <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__offset <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__pressure_electron <: FDS
    var"pedestal_width" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_width = summary__pedestal_fits__mtanh__pressure_electron__pedestal_width()
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position = summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position()
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max = summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max()
    var"offset" :: summary__pedestal_fits__mtanh__pressure_electron__offset = summary__pedestal_fits__mtanh__pressure_electron__offset()
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm = summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm()
    var"separatrix" :: summary__pedestal_fits__mtanh__pressure_electron__separatrix = summary__pedestal_fits__mtanh__pressure_electron__separatrix()
    var"pedestal_position" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_position = summary__pedestal_fits__mtanh__pressure_electron__pedestal_position()
    var"pedestal_height" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_height = summary__pedestal_fits__mtanh__pressure_electron__pedestal_height()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__pressure_electron(var"pedestal_width"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__mtanh__pressure_electron__offset(), var"d_dpsi_norm"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__mtanh__pressure_electron__separatrix(), var"pedestal_position"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_height(), _parent=WeakRef(missing))
        fds = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(fds)
        setfield!(fds.pedestal_width, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max_position, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max, :_parent, WeakRef(fds))
        setfield!(fds.offset, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm, :_parent, WeakRef(fds))
        setfield!(fds.separatrix, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_position, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_height, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__nustar_pedestal_top_electron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__nustar_pedestal_top_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__separatrix <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_width <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_height <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__offset <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__n_e <: FDS
    var"pedestal_width" :: summary__pedestal_fits__mtanh__n_e__pedestal_width = summary__pedestal_fits__mtanh__n_e__pedestal_width()
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position = summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position()
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max = summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max()
    var"offset" :: summary__pedestal_fits__mtanh__n_e__offset = summary__pedestal_fits__mtanh__n_e__offset()
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm = summary__pedestal_fits__mtanh__n_e__d_dpsi_norm()
    var"separatrix" :: summary__pedestal_fits__mtanh__n_e__separatrix = summary__pedestal_fits__mtanh__n_e__separatrix()
    var"pedestal_position" :: summary__pedestal_fits__mtanh__n_e__pedestal_position = summary__pedestal_fits__mtanh__n_e__pedestal_position()
    var"pedestal_height" :: summary__pedestal_fits__mtanh__n_e__pedestal_height = summary__pedestal_fits__mtanh__n_e__pedestal_height()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__n_e(var"pedestal_width"=summary__pedestal_fits__mtanh__n_e__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__mtanh__n_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__mtanh__n_e__separatrix(), var"pedestal_position"=summary__pedestal_fits__mtanh__n_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__mtanh__n_e__pedestal_height(), _parent=WeakRef(missing))
        fds = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(fds)
        setfield!(fds.pedestal_width, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max_position, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max, :_parent, WeakRef(fds))
        setfield!(fds.offset, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm, :_parent, WeakRef(fds))
        setfield!(fds.separatrix, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_position, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_height, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh__alpha_electron_pedestal_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh__alpha_electron_pedestal_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__mtanh <: FDS
    var"pressure_electron" :: summary__pedestal_fits__mtanh__pressure_electron = summary__pedestal_fits__mtanh__pressure_electron()
    var"b_field_pol_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs = summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs()
    var"alpha_electron_pedestal_max_position" :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position = summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position()
    var"rhostar_pedestal_top_electron_hfs" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs = summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs()
    var"beta_pol_pedestal_top_electron_hfs" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs = summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs()
    var"energy_thermal_pedestal_electron" :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron = summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron()
    var"parameters" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"t_e" :: summary__pedestal_fits__mtanh__t_e = summary__pedestal_fits__mtanh__t_e()
    var"rhostar_pedestal_top_electron_lfs" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs = summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs()
    var"beta_pol_pedestal_top_electron_lfs" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs = summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs()
    var"stability" :: summary__pedestal_fits__mtanh__stability = summary__pedestal_fits__mtanh__stability()
    var"b_field_tor_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs = summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs()
    var"volume_inside_pedestal" :: summary__pedestal_fits__mtanh__volume_inside_pedestal = summary__pedestal_fits__mtanh__volume_inside_pedestal()
    var"b_field_pol_pedestal_top_average" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average = summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average()
    var"coulomb_factor_pedestal_top" :: summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top = summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top()
    var"beta_pol_pedestal_top_electron_average" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average = summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average()
    var"b_field_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs = summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs()
    var"b_field_pol_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs = summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs()
    var"nustar_pedestal_top_electron" :: summary__pedestal_fits__mtanh__nustar_pedestal_top_electron = summary__pedestal_fits__mtanh__nustar_pedestal_top_electron()
    var"alpha_electron_pedestal_max" :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max = summary__pedestal_fits__mtanh__alpha_electron_pedestal_max()
    var"b_field_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs = summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs()
    var"n_e" :: summary__pedestal_fits__mtanh__n_e = summary__pedestal_fits__mtanh__n_e()
    var"energy_thermal_pedestal_ion" :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion = summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion()
    var"b_field_tor_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs = summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs()
    var"rhostar_pedestal_top_electron_magnetic_axis" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis = summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__mtanh(var"pressure_electron"=summary__pedestal_fits__mtanh__pressure_electron(), var"b_field_pol_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs(), var"alpha_electron_pedestal_max_position"=summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position(), var"rhostar_pedestal_top_electron_hfs"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs(), var"beta_pol_pedestal_top_electron_hfs"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs(), var"energy_thermal_pedestal_electron"=summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron(), var"parameters"=missing, var"t_e"=summary__pedestal_fits__mtanh__t_e(), var"rhostar_pedestal_top_electron_lfs"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs(), var"beta_pol_pedestal_top_electron_lfs"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs(), var"stability"=summary__pedestal_fits__mtanh__stability(), var"b_field_tor_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs(), var"volume_inside_pedestal"=summary__pedestal_fits__mtanh__volume_inside_pedestal(), var"b_field_pol_pedestal_top_average"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average(), var"coulomb_factor_pedestal_top"=summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top(), var"beta_pol_pedestal_top_electron_average"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average(), var"b_field_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs(), var"b_field_pol_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs(), var"nustar_pedestal_top_electron"=summary__pedestal_fits__mtanh__nustar_pedestal_top_electron(), var"alpha_electron_pedestal_max"=summary__pedestal_fits__mtanh__alpha_electron_pedestal_max(), var"b_field_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs(), var"n_e"=summary__pedestal_fits__mtanh__n_e(), var"energy_thermal_pedestal_ion"=summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion(), var"b_field_tor_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs(), var"rhostar_pedestal_top_electron_magnetic_axis"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis(), _parent=WeakRef(missing))
        fds = new(var"pressure_electron", var"b_field_pol_pedestal_top_hfs", var"alpha_electron_pedestal_max_position", var"rhostar_pedestal_top_electron_hfs", var"beta_pol_pedestal_top_electron_hfs", var"energy_thermal_pedestal_electron", var"parameters", var"t_e", var"rhostar_pedestal_top_electron_lfs", var"beta_pol_pedestal_top_electron_lfs", var"stability", var"b_field_tor_pedestal_top_lfs", var"volume_inside_pedestal", var"b_field_pol_pedestal_top_average", var"coulomb_factor_pedestal_top", var"beta_pol_pedestal_top_electron_average", var"b_field_pedestal_top_lfs", var"b_field_pol_pedestal_top_lfs", var"nustar_pedestal_top_electron", var"alpha_electron_pedestal_max", var"b_field_pedestal_top_hfs", var"n_e", var"energy_thermal_pedestal_ion", var"b_field_tor_pedestal_top_hfs", var"rhostar_pedestal_top_electron_magnetic_axis", _parent)
        assign_expressions(fds)
        setfield!(fds.pressure_electron, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pol_pedestal_top_hfs, :_parent, WeakRef(fds))
        setfield!(fds.alpha_electron_pedestal_max_position, :_parent, WeakRef(fds))
        setfield!(fds.rhostar_pedestal_top_electron_hfs, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_pedestal_top_electron_hfs, :_parent, WeakRef(fds))
        setfield!(fds.energy_thermal_pedestal_electron, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.rhostar_pedestal_top_electron_lfs, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_pedestal_top_electron_lfs, :_parent, WeakRef(fds))
        setfield!(fds.stability, :_parent, WeakRef(fds))
        setfield!(fds.b_field_tor_pedestal_top_lfs, :_parent, WeakRef(fds))
        setfield!(fds.volume_inside_pedestal, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pol_pedestal_top_average, :_parent, WeakRef(fds))
        setfield!(fds.coulomb_factor_pedestal_top, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_pedestal_top_electron_average, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pedestal_top_lfs, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pol_pedestal_top_lfs, :_parent, WeakRef(fds))
        setfield!(fds.nustar_pedestal_top_electron, :_parent, WeakRef(fds))
        setfield!(fds.alpha_electron_pedestal_max, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pedestal_top_hfs, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.energy_thermal_pedestal_ion, :_parent, WeakRef(fds))
        setfield!(fds.b_field_tor_pedestal_top_hfs, :_parent, WeakRef(fds))
        setfield!(fds.rhostar_pedestal_top_electron_magnetic_axis, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__volume_inside_pedestal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__volume_inside_pedestal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__pedestal_width <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__pedestal_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__pedestal_height <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__offset <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__d_dpsi_norm_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e__d_dpsi_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__t_e <: FDS
    var"pedestal_width" :: summary__pedestal_fits__linear__t_e__pedestal_width = summary__pedestal_fits__linear__t_e__pedestal_width()
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__t_e__d_dpsi_norm_max = summary__pedestal_fits__linear__t_e__d_dpsi_norm_max()
    var"offset" :: summary__pedestal_fits__linear__t_e__offset = summary__pedestal_fits__linear__t_e__offset()
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__t_e__d_dpsi_norm = summary__pedestal_fits__linear__t_e__d_dpsi_norm()
    var"pedestal_position" :: summary__pedestal_fits__linear__t_e__pedestal_position = summary__pedestal_fits__linear__t_e__pedestal_position()
    var"pedestal_height" :: summary__pedestal_fits__linear__t_e__pedestal_height = summary__pedestal_fits__linear__t_e__pedestal_height()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__t_e(var"pedestal_width"=summary__pedestal_fits__linear__t_e__pedestal_width(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__t_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__t_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__linear__t_e__d_dpsi_norm(), var"pedestal_position"=summary__pedestal_fits__linear__t_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__linear__t_e__pedestal_height(), _parent=WeakRef(missing))
        fds = new(var"pedestal_width", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(fds)
        setfield!(fds.pedestal_width, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max, :_parent, WeakRef(fds))
        setfield!(fds.offset, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_position, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_height, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__separatrix <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_width <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_height <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__offset <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__pressure_electron <: FDS
    var"pedestal_width" :: summary__pedestal_fits__linear__pressure_electron__pedestal_width = summary__pedestal_fits__linear__pressure_electron__pedestal_width()
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position = summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position()
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max = summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max()
    var"offset" :: summary__pedestal_fits__linear__pressure_electron__offset = summary__pedestal_fits__linear__pressure_electron__offset()
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm = summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm()
    var"separatrix" :: summary__pedestal_fits__linear__pressure_electron__separatrix = summary__pedestal_fits__linear__pressure_electron__separatrix()
    var"pedestal_position" :: summary__pedestal_fits__linear__pressure_electron__pedestal_position = summary__pedestal_fits__linear__pressure_electron__pedestal_position()
    var"pedestal_height" :: summary__pedestal_fits__linear__pressure_electron__pedestal_height = summary__pedestal_fits__linear__pressure_electron__pedestal_height()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__pressure_electron(var"pedestal_width"=summary__pedestal_fits__linear__pressure_electron__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__pressure_electron__offset(), var"d_dpsi_norm"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__linear__pressure_electron__separatrix(), var"pedestal_position"=summary__pedestal_fits__linear__pressure_electron__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__linear__pressure_electron__pedestal_height(), _parent=WeakRef(missing))
        fds = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(fds)
        setfield!(fds.pedestal_width, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max_position, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max, :_parent, WeakRef(fds))
        setfield!(fds.offset, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm, :_parent, WeakRef(fds))
        setfield!(fds.separatrix, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_position, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_height, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__nustar_pedestal_top_electron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__nustar_pedestal_top_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__separatrix <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__pedestal_width <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__pedestal_position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__pedestal_height <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__offset <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__d_dpsi_norm_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e__d_dpsi_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__n_e <: FDS
    var"pedestal_width" :: summary__pedestal_fits__linear__n_e__pedestal_width = summary__pedestal_fits__linear__n_e__pedestal_width()
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__n_e__d_dpsi_norm_max = summary__pedestal_fits__linear__n_e__d_dpsi_norm_max()
    var"offset" :: summary__pedestal_fits__linear__n_e__offset = summary__pedestal_fits__linear__n_e__offset()
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__n_e__d_dpsi_norm = summary__pedestal_fits__linear__n_e__d_dpsi_norm()
    var"separatrix" :: summary__pedestal_fits__linear__n_e__separatrix = summary__pedestal_fits__linear__n_e__separatrix()
    var"pedestal_position" :: summary__pedestal_fits__linear__n_e__pedestal_position = summary__pedestal_fits__linear__n_e__pedestal_position()
    var"pedestal_height" :: summary__pedestal_fits__linear__n_e__pedestal_height = summary__pedestal_fits__linear__n_e__pedestal_height()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__n_e(var"pedestal_width"=summary__pedestal_fits__linear__n_e__pedestal_width(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__n_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__n_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__linear__n_e__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__linear__n_e__separatrix(), var"pedestal_position"=summary__pedestal_fits__linear__n_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__linear__n_e__pedestal_height(), _parent=WeakRef(missing))
        fds = new(var"pedestal_width", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(fds)
        setfield!(fds.pedestal_width, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm_max, :_parent, WeakRef(fds))
        setfield!(fds.offset, :_parent, WeakRef(fds))
        setfield!(fds.d_dpsi_norm, :_parent, WeakRef(fds))
        setfield!(fds.separatrix, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_position, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_height, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__energy_thermal_pedestal_ion <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__energy_thermal_pedestal_ion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__energy_thermal_pedestal_electron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__energy_thermal_pedestal_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__coulomb_factor_pedestal_top <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__coulomb_factor_pedestal_top(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_pol_pedestal_top_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pedestal_top_lfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear__b_field_pedestal_top_hfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear__b_field_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits__linear <: FDS
    var"pressure_electron" :: summary__pedestal_fits__linear__pressure_electron = summary__pedestal_fits__linear__pressure_electron()
    var"b_field_pol_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs = summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs()
    var"rhostar_pedestal_top_electron_hfs" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs = summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs()
    var"beta_pol_pedestal_top_electron_hfs" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs = summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs()
    var"energy_thermal_pedestal_electron" :: summary__pedestal_fits__linear__energy_thermal_pedestal_electron = summary__pedestal_fits__linear__energy_thermal_pedestal_electron()
    var"parameters" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"t_e" :: summary__pedestal_fits__linear__t_e = summary__pedestal_fits__linear__t_e()
    var"rhostar_pedestal_top_electron_lfs" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs = summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs()
    var"beta_pol_pedestal_top_electron_lfs" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs = summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs()
    var"b_field_tor_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs = summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs()
    var"volume_inside_pedestal" :: summary__pedestal_fits__linear__volume_inside_pedestal = summary__pedestal_fits__linear__volume_inside_pedestal()
    var"b_field_pol_pedestal_top_average" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_average = summary__pedestal_fits__linear__b_field_pol_pedestal_top_average()
    var"coulomb_factor_pedestal_top" :: summary__pedestal_fits__linear__coulomb_factor_pedestal_top = summary__pedestal_fits__linear__coulomb_factor_pedestal_top()
    var"beta_pol_pedestal_top_electron_average" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average = summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average()
    var"b_field_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_pedestal_top_lfs = summary__pedestal_fits__linear__b_field_pedestal_top_lfs()
    var"b_field_pol_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs = summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs()
    var"nustar_pedestal_top_electron" :: summary__pedestal_fits__linear__nustar_pedestal_top_electron = summary__pedestal_fits__linear__nustar_pedestal_top_electron()
    var"b_field_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_pedestal_top_hfs = summary__pedestal_fits__linear__b_field_pedestal_top_hfs()
    var"n_e" :: summary__pedestal_fits__linear__n_e = summary__pedestal_fits__linear__n_e()
    var"energy_thermal_pedestal_ion" :: summary__pedestal_fits__linear__energy_thermal_pedestal_ion = summary__pedestal_fits__linear__energy_thermal_pedestal_ion()
    var"b_field_tor_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs = summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs()
    var"rhostar_pedestal_top_electron_magnetic_axis" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis = summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits__linear(var"pressure_electron"=summary__pedestal_fits__linear__pressure_electron(), var"b_field_pol_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs(), var"rhostar_pedestal_top_electron_hfs"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs(), var"beta_pol_pedestal_top_electron_hfs"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs(), var"energy_thermal_pedestal_electron"=summary__pedestal_fits__linear__energy_thermal_pedestal_electron(), var"parameters"=missing, var"t_e"=summary__pedestal_fits__linear__t_e(), var"rhostar_pedestal_top_electron_lfs"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs(), var"beta_pol_pedestal_top_electron_lfs"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs(), var"b_field_tor_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs(), var"volume_inside_pedestal"=summary__pedestal_fits__linear__volume_inside_pedestal(), var"b_field_pol_pedestal_top_average"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_average(), var"coulomb_factor_pedestal_top"=summary__pedestal_fits__linear__coulomb_factor_pedestal_top(), var"beta_pol_pedestal_top_electron_average"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average(), var"b_field_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_pedestal_top_lfs(), var"b_field_pol_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs(), var"nustar_pedestal_top_electron"=summary__pedestal_fits__linear__nustar_pedestal_top_electron(), var"b_field_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_pedestal_top_hfs(), var"n_e"=summary__pedestal_fits__linear__n_e(), var"energy_thermal_pedestal_ion"=summary__pedestal_fits__linear__energy_thermal_pedestal_ion(), var"b_field_tor_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs(), var"rhostar_pedestal_top_electron_magnetic_axis"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis(), _parent=WeakRef(missing))
        fds = new(var"pressure_electron", var"b_field_pol_pedestal_top_hfs", var"rhostar_pedestal_top_electron_hfs", var"beta_pol_pedestal_top_electron_hfs", var"energy_thermal_pedestal_electron", var"parameters", var"t_e", var"rhostar_pedestal_top_electron_lfs", var"beta_pol_pedestal_top_electron_lfs", var"b_field_tor_pedestal_top_lfs", var"volume_inside_pedestal", var"b_field_pol_pedestal_top_average", var"coulomb_factor_pedestal_top", var"beta_pol_pedestal_top_electron_average", var"b_field_pedestal_top_lfs", var"b_field_pol_pedestal_top_lfs", var"nustar_pedestal_top_electron", var"b_field_pedestal_top_hfs", var"n_e", var"energy_thermal_pedestal_ion", var"b_field_tor_pedestal_top_hfs", var"rhostar_pedestal_top_electron_magnetic_axis", _parent)
        assign_expressions(fds)
        setfield!(fds.pressure_electron, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pol_pedestal_top_hfs, :_parent, WeakRef(fds))
        setfield!(fds.rhostar_pedestal_top_electron_hfs, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_pedestal_top_electron_hfs, :_parent, WeakRef(fds))
        setfield!(fds.energy_thermal_pedestal_electron, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.rhostar_pedestal_top_electron_lfs, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_pedestal_top_electron_lfs, :_parent, WeakRef(fds))
        setfield!(fds.b_field_tor_pedestal_top_lfs, :_parent, WeakRef(fds))
        setfield!(fds.volume_inside_pedestal, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pol_pedestal_top_average, :_parent, WeakRef(fds))
        setfield!(fds.coulomb_factor_pedestal_top, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_pedestal_top_electron_average, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pedestal_top_lfs, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pol_pedestal_top_lfs, :_parent, WeakRef(fds))
        setfield!(fds.nustar_pedestal_top_electron, :_parent, WeakRef(fds))
        setfield!(fds.b_field_pedestal_top_hfs, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.energy_thermal_pedestal_ion, :_parent, WeakRef(fds))
        setfield!(fds.b_field_tor_pedestal_top_hfs, :_parent, WeakRef(fds))
        setfield!(fds.rhostar_pedestal_top_electron_magnetic_axis, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__pedestal_fits <: FDS
    var"linear" :: summary__pedestal_fits__linear = summary__pedestal_fits__linear()
    var"mtanh" :: summary__pedestal_fits__mtanh = summary__pedestal_fits__mtanh()
    _parent :: WeakRef = WeakRef(missing)
    function summary__pedestal_fits(var"linear"=summary__pedestal_fits__linear(), var"mtanh"=summary__pedestal_fits__mtanh(), _parent=WeakRef(missing))
        fds = new(var"linear", var"mtanh", _parent)
        assign_expressions(fds)
        setfield!(fds.linear, :_parent, WeakRef(fds))
        setfield!(fds.mtanh, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__midplane <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__midplane(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__magnetic_shear_flag <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__magnetic_shear_flag(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__velocity_tor <: FDS
    var"krypton" :: summary__local__separatrix__velocity_tor__krypton = summary__local__separatrix__velocity_tor__krypton()
    var"lithium" :: summary__local__separatrix__velocity_tor__lithium = summary__local__separatrix__velocity_tor__lithium()
    var"neon" :: summary__local__separatrix__velocity_tor__neon = summary__local__separatrix__velocity_tor__neon()
    var"tritium" :: summary__local__separatrix__velocity_tor__tritium = summary__local__separatrix__velocity_tor__tritium()
    var"helium_3" :: summary__local__separatrix__velocity_tor__helium_3 = summary__local__separatrix__velocity_tor__helium_3()
    var"deuterium" :: summary__local__separatrix__velocity_tor__deuterium = summary__local__separatrix__velocity_tor__deuterium()
    var"iron" :: summary__local__separatrix__velocity_tor__iron = summary__local__separatrix__velocity_tor__iron()
    var"helium_4" :: summary__local__separatrix__velocity_tor__helium_4 = summary__local__separatrix__velocity_tor__helium_4()
    var"oxygen" :: summary__local__separatrix__velocity_tor__oxygen = summary__local__separatrix__velocity_tor__oxygen()
    var"tungsten" :: summary__local__separatrix__velocity_tor__tungsten = summary__local__separatrix__velocity_tor__tungsten()
    var"xenon" :: summary__local__separatrix__velocity_tor__xenon = summary__local__separatrix__velocity_tor__xenon()
    var"hydrogen" :: summary__local__separatrix__velocity_tor__hydrogen = summary__local__separatrix__velocity_tor__hydrogen()
    var"carbon" :: summary__local__separatrix__velocity_tor__carbon = summary__local__separatrix__velocity_tor__carbon()
    var"nitrogen" :: summary__local__separatrix__velocity_tor__nitrogen = summary__local__separatrix__velocity_tor__nitrogen()
    var"beryllium" :: summary__local__separatrix__velocity_tor__beryllium = summary__local__separatrix__velocity_tor__beryllium()
    var"argon" :: summary__local__separatrix__velocity_tor__argon = summary__local__separatrix__velocity_tor__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__velocity_tor(var"krypton"=summary__local__separatrix__velocity_tor__krypton(), var"lithium"=summary__local__separatrix__velocity_tor__lithium(), var"neon"=summary__local__separatrix__velocity_tor__neon(), var"tritium"=summary__local__separatrix__velocity_tor__tritium(), var"helium_3"=summary__local__separatrix__velocity_tor__helium_3(), var"deuterium"=summary__local__separatrix__velocity_tor__deuterium(), var"iron"=summary__local__separatrix__velocity_tor__iron(), var"helium_4"=summary__local__separatrix__velocity_tor__helium_4(), var"oxygen"=summary__local__separatrix__velocity_tor__oxygen(), var"tungsten"=summary__local__separatrix__velocity_tor__tungsten(), var"xenon"=summary__local__separatrix__velocity_tor__xenon(), var"hydrogen"=summary__local__separatrix__velocity_tor__hydrogen(), var"carbon"=summary__local__separatrix__velocity_tor__carbon(), var"nitrogen"=summary__local__separatrix__velocity_tor__nitrogen(), var"beryllium"=summary__local__separatrix__velocity_tor__beryllium(), var"argon"=summary__local__separatrix__velocity_tor__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__q <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__position <: FDS
    var"psi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__position(var"psi"=missing, var"rho_tor_norm"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"rho_tor_norm", var"rho_tor", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_i <: FDS
    var"krypton" :: summary__local__separatrix__n_i__krypton = summary__local__separatrix__n_i__krypton()
    var"lithium" :: summary__local__separatrix__n_i__lithium = summary__local__separatrix__n_i__lithium()
    var"neon" :: summary__local__separatrix__n_i__neon = summary__local__separatrix__n_i__neon()
    var"tritium" :: summary__local__separatrix__n_i__tritium = summary__local__separatrix__n_i__tritium()
    var"helium_3" :: summary__local__separatrix__n_i__helium_3 = summary__local__separatrix__n_i__helium_3()
    var"deuterium" :: summary__local__separatrix__n_i__deuterium = summary__local__separatrix__n_i__deuterium()
    var"iron" :: summary__local__separatrix__n_i__iron = summary__local__separatrix__n_i__iron()
    var"helium_4" :: summary__local__separatrix__n_i__helium_4 = summary__local__separatrix__n_i__helium_4()
    var"oxygen" :: summary__local__separatrix__n_i__oxygen = summary__local__separatrix__n_i__oxygen()
    var"tungsten" :: summary__local__separatrix__n_i__tungsten = summary__local__separatrix__n_i__tungsten()
    var"xenon" :: summary__local__separatrix__n_i__xenon = summary__local__separatrix__n_i__xenon()
    var"hydrogen" :: summary__local__separatrix__n_i__hydrogen = summary__local__separatrix__n_i__hydrogen()
    var"carbon" :: summary__local__separatrix__n_i__carbon = summary__local__separatrix__n_i__carbon()
    var"nitrogen" :: summary__local__separatrix__n_i__nitrogen = summary__local__separatrix__n_i__nitrogen()
    var"beryllium" :: summary__local__separatrix__n_i__beryllium = summary__local__separatrix__n_i__beryllium()
    var"argon" :: summary__local__separatrix__n_i__argon = summary__local__separatrix__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_i(var"krypton"=summary__local__separatrix__n_i__krypton(), var"lithium"=summary__local__separatrix__n_i__lithium(), var"neon"=summary__local__separatrix__n_i__neon(), var"tritium"=summary__local__separatrix__n_i__tritium(), var"helium_3"=summary__local__separatrix__n_i__helium_3(), var"deuterium"=summary__local__separatrix__n_i__deuterium(), var"iron"=summary__local__separatrix__n_i__iron(), var"helium_4"=summary__local__separatrix__n_i__helium_4(), var"oxygen"=summary__local__separatrix__n_i__oxygen(), var"tungsten"=summary__local__separatrix__n_i__tungsten(), var"xenon"=summary__local__separatrix__n_i__xenon(), var"hydrogen"=summary__local__separatrix__n_i__hydrogen(), var"carbon"=summary__local__separatrix__n_i__carbon(), var"nitrogen"=summary__local__separatrix__n_i__nitrogen(), var"beryllium"=summary__local__separatrix__n_i__beryllium(), var"argon"=summary__local__separatrix__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__momentum_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__magnetic_shear <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix__e_field_parallel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__separatrix <: FDS
    var"n_i" :: summary__local__separatrix__n_i = summary__local__separatrix__n_i()
    var"velocity_tor" :: summary__local__separatrix__velocity_tor = summary__local__separatrix__velocity_tor()
    var"magnetic_shear" :: summary__local__separatrix__magnetic_shear = summary__local__separatrix__magnetic_shear()
    var"t_e" :: summary__local__separatrix__t_e = summary__local__separatrix__t_e()
    var"t_i_average" :: summary__local__separatrix__t_i_average = summary__local__separatrix__t_i_average()
    var"position" :: summary__local__separatrix__position = summary__local__separatrix__position()
    var"e_field_parallel" :: summary__local__separatrix__e_field_parallel = summary__local__separatrix__e_field_parallel()
    var"momentum_tor" :: summary__local__separatrix__momentum_tor = summary__local__separatrix__momentum_tor()
    var"n_i_total" :: summary__local__separatrix__n_i_total = summary__local__separatrix__n_i_total()
    var"q" :: summary__local__separatrix__q = summary__local__separatrix__q()
    var"n_e" :: summary__local__separatrix__n_e = summary__local__separatrix__n_e()
    var"zeff" :: summary__local__separatrix__zeff = summary__local__separatrix__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__separatrix(var"n_i"=summary__local__separatrix__n_i(), var"velocity_tor"=summary__local__separatrix__velocity_tor(), var"magnetic_shear"=summary__local__separatrix__magnetic_shear(), var"t_e"=summary__local__separatrix__t_e(), var"t_i_average"=summary__local__separatrix__t_i_average(), var"position"=summary__local__separatrix__position(), var"e_field_parallel"=summary__local__separatrix__e_field_parallel(), var"momentum_tor"=summary__local__separatrix__momentum_tor(), var"n_i_total"=summary__local__separatrix__n_i_total(), var"q"=summary__local__separatrix__q(), var"n_e"=summary__local__separatrix__n_e(), var"zeff"=summary__local__separatrix__zeff(), _parent=WeakRef(missing))
        fds = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"position", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"q", var"n_e", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.velocity_tor, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_shear, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.e_field_parallel, :_parent, WeakRef(fds))
        setfield!(fds.momentum_tor, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.q, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__r_eff_norm_2_3__plateau_factor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__r_eff_norm_2_3__plateau_factor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__r_eff_norm_2_3__iota <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__r_eff_norm_2_3__iota(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__r_eff_norm_2_3__effective_helical_ripple <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__r_eff_norm_2_3__effective_helical_ripple(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__r_eff_norm_2_3 <: FDS
    var"plateau_factor" :: summary__local__r_eff_norm_2_3__plateau_factor = summary__local__r_eff_norm_2_3__plateau_factor()
    var"iota" :: summary__local__r_eff_norm_2_3__iota = summary__local__r_eff_norm_2_3__iota()
    var"effective_helical_ripple" :: summary__local__r_eff_norm_2_3__effective_helical_ripple = summary__local__r_eff_norm_2_3__effective_helical_ripple()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__r_eff_norm_2_3(var"plateau_factor"=summary__local__r_eff_norm_2_3__plateau_factor(), var"iota"=summary__local__r_eff_norm_2_3__iota(), var"effective_helical_ripple"=summary__local__r_eff_norm_2_3__effective_helical_ripple(), _parent=WeakRef(missing))
        fds = new(var"plateau_factor", var"iota", var"effective_helical_ripple", _parent)
        assign_expressions(fds)
        setfield!(fds.plateau_factor, :_parent, WeakRef(fds))
        setfield!(fds.iota, :_parent, WeakRef(fds))
        setfield!(fds.effective_helical_ripple, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__velocity_tor <: FDS
    var"krypton" :: summary__local__pedestal__velocity_tor__krypton = summary__local__pedestal__velocity_tor__krypton()
    var"lithium" :: summary__local__pedestal__velocity_tor__lithium = summary__local__pedestal__velocity_tor__lithium()
    var"neon" :: summary__local__pedestal__velocity_tor__neon = summary__local__pedestal__velocity_tor__neon()
    var"tritium" :: summary__local__pedestal__velocity_tor__tritium = summary__local__pedestal__velocity_tor__tritium()
    var"helium_3" :: summary__local__pedestal__velocity_tor__helium_3 = summary__local__pedestal__velocity_tor__helium_3()
    var"deuterium" :: summary__local__pedestal__velocity_tor__deuterium = summary__local__pedestal__velocity_tor__deuterium()
    var"iron" :: summary__local__pedestal__velocity_tor__iron = summary__local__pedestal__velocity_tor__iron()
    var"helium_4" :: summary__local__pedestal__velocity_tor__helium_4 = summary__local__pedestal__velocity_tor__helium_4()
    var"oxygen" :: summary__local__pedestal__velocity_tor__oxygen = summary__local__pedestal__velocity_tor__oxygen()
    var"tungsten" :: summary__local__pedestal__velocity_tor__tungsten = summary__local__pedestal__velocity_tor__tungsten()
    var"xenon" :: summary__local__pedestal__velocity_tor__xenon = summary__local__pedestal__velocity_tor__xenon()
    var"hydrogen" :: summary__local__pedestal__velocity_tor__hydrogen = summary__local__pedestal__velocity_tor__hydrogen()
    var"carbon" :: summary__local__pedestal__velocity_tor__carbon = summary__local__pedestal__velocity_tor__carbon()
    var"nitrogen" :: summary__local__pedestal__velocity_tor__nitrogen = summary__local__pedestal__velocity_tor__nitrogen()
    var"beryllium" :: summary__local__pedestal__velocity_tor__beryllium = summary__local__pedestal__velocity_tor__beryllium()
    var"argon" :: summary__local__pedestal__velocity_tor__argon = summary__local__pedestal__velocity_tor__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__velocity_tor(var"krypton"=summary__local__pedestal__velocity_tor__krypton(), var"lithium"=summary__local__pedestal__velocity_tor__lithium(), var"neon"=summary__local__pedestal__velocity_tor__neon(), var"tritium"=summary__local__pedestal__velocity_tor__tritium(), var"helium_3"=summary__local__pedestal__velocity_tor__helium_3(), var"deuterium"=summary__local__pedestal__velocity_tor__deuterium(), var"iron"=summary__local__pedestal__velocity_tor__iron(), var"helium_4"=summary__local__pedestal__velocity_tor__helium_4(), var"oxygen"=summary__local__pedestal__velocity_tor__oxygen(), var"tungsten"=summary__local__pedestal__velocity_tor__tungsten(), var"xenon"=summary__local__pedestal__velocity_tor__xenon(), var"hydrogen"=summary__local__pedestal__velocity_tor__hydrogen(), var"carbon"=summary__local__pedestal__velocity_tor__carbon(), var"nitrogen"=summary__local__pedestal__velocity_tor__nitrogen(), var"beryllium"=summary__local__pedestal__velocity_tor__beryllium(), var"argon"=summary__local__pedestal__velocity_tor__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__q <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__position <: FDS
    var"psi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__position(var"psi"=missing, var"rho_tor_norm"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"rho_tor_norm", var"rho_tor", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_i <: FDS
    var"krypton" :: summary__local__pedestal__n_i__krypton = summary__local__pedestal__n_i__krypton()
    var"lithium" :: summary__local__pedestal__n_i__lithium = summary__local__pedestal__n_i__lithium()
    var"neon" :: summary__local__pedestal__n_i__neon = summary__local__pedestal__n_i__neon()
    var"tritium" :: summary__local__pedestal__n_i__tritium = summary__local__pedestal__n_i__tritium()
    var"helium_3" :: summary__local__pedestal__n_i__helium_3 = summary__local__pedestal__n_i__helium_3()
    var"deuterium" :: summary__local__pedestal__n_i__deuterium = summary__local__pedestal__n_i__deuterium()
    var"iron" :: summary__local__pedestal__n_i__iron = summary__local__pedestal__n_i__iron()
    var"helium_4" :: summary__local__pedestal__n_i__helium_4 = summary__local__pedestal__n_i__helium_4()
    var"oxygen" :: summary__local__pedestal__n_i__oxygen = summary__local__pedestal__n_i__oxygen()
    var"tungsten" :: summary__local__pedestal__n_i__tungsten = summary__local__pedestal__n_i__tungsten()
    var"xenon" :: summary__local__pedestal__n_i__xenon = summary__local__pedestal__n_i__xenon()
    var"hydrogen" :: summary__local__pedestal__n_i__hydrogen = summary__local__pedestal__n_i__hydrogen()
    var"carbon" :: summary__local__pedestal__n_i__carbon = summary__local__pedestal__n_i__carbon()
    var"nitrogen" :: summary__local__pedestal__n_i__nitrogen = summary__local__pedestal__n_i__nitrogen()
    var"beryllium" :: summary__local__pedestal__n_i__beryllium = summary__local__pedestal__n_i__beryllium()
    var"argon" :: summary__local__pedestal__n_i__argon = summary__local__pedestal__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_i(var"krypton"=summary__local__pedestal__n_i__krypton(), var"lithium"=summary__local__pedestal__n_i__lithium(), var"neon"=summary__local__pedestal__n_i__neon(), var"tritium"=summary__local__pedestal__n_i__tritium(), var"helium_3"=summary__local__pedestal__n_i__helium_3(), var"deuterium"=summary__local__pedestal__n_i__deuterium(), var"iron"=summary__local__pedestal__n_i__iron(), var"helium_4"=summary__local__pedestal__n_i__helium_4(), var"oxygen"=summary__local__pedestal__n_i__oxygen(), var"tungsten"=summary__local__pedestal__n_i__tungsten(), var"xenon"=summary__local__pedestal__n_i__xenon(), var"hydrogen"=summary__local__pedestal__n_i__hydrogen(), var"carbon"=summary__local__pedestal__n_i__carbon(), var"nitrogen"=summary__local__pedestal__n_i__nitrogen(), var"beryllium"=summary__local__pedestal__n_i__beryllium(), var"argon"=summary__local__pedestal__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__momentum_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__magnetic_shear <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal__e_field_parallel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__pedestal <: FDS
    var"n_i" :: summary__local__pedestal__n_i = summary__local__pedestal__n_i()
    var"velocity_tor" :: summary__local__pedestal__velocity_tor = summary__local__pedestal__velocity_tor()
    var"magnetic_shear" :: summary__local__pedestal__magnetic_shear = summary__local__pedestal__magnetic_shear()
    var"t_e" :: summary__local__pedestal__t_e = summary__local__pedestal__t_e()
    var"t_i_average" :: summary__local__pedestal__t_i_average = summary__local__pedestal__t_i_average()
    var"position" :: summary__local__pedestal__position = summary__local__pedestal__position()
    var"e_field_parallel" :: summary__local__pedestal__e_field_parallel = summary__local__pedestal__e_field_parallel()
    var"momentum_tor" :: summary__local__pedestal__momentum_tor = summary__local__pedestal__momentum_tor()
    var"n_i_total" :: summary__local__pedestal__n_i_total = summary__local__pedestal__n_i_total()
    var"q" :: summary__local__pedestal__q = summary__local__pedestal__q()
    var"n_e" :: summary__local__pedestal__n_e = summary__local__pedestal__n_e()
    var"zeff" :: summary__local__pedestal__zeff = summary__local__pedestal__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__pedestal(var"n_i"=summary__local__pedestal__n_i(), var"velocity_tor"=summary__local__pedestal__velocity_tor(), var"magnetic_shear"=summary__local__pedestal__magnetic_shear(), var"t_e"=summary__local__pedestal__t_e(), var"t_i_average"=summary__local__pedestal__t_i_average(), var"position"=summary__local__pedestal__position(), var"e_field_parallel"=summary__local__pedestal__e_field_parallel(), var"momentum_tor"=summary__local__pedestal__momentum_tor(), var"n_i_total"=summary__local__pedestal__n_i_total(), var"q"=summary__local__pedestal__q(), var"n_e"=summary__local__pedestal__n_e(), var"zeff"=summary__local__pedestal__zeff(), _parent=WeakRef(missing))
        fds = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"position", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"q", var"n_e", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.velocity_tor, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_shear, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.e_field_parallel, :_parent, WeakRef(fds))
        setfield!(fds.momentum_tor, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.q, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__velocity_tor <: FDS
    var"krypton" :: summary__local__magnetic_axis__velocity_tor__krypton = summary__local__magnetic_axis__velocity_tor__krypton()
    var"lithium" :: summary__local__magnetic_axis__velocity_tor__lithium = summary__local__magnetic_axis__velocity_tor__lithium()
    var"neon" :: summary__local__magnetic_axis__velocity_tor__neon = summary__local__magnetic_axis__velocity_tor__neon()
    var"tritium" :: summary__local__magnetic_axis__velocity_tor__tritium = summary__local__magnetic_axis__velocity_tor__tritium()
    var"helium_3" :: summary__local__magnetic_axis__velocity_tor__helium_3 = summary__local__magnetic_axis__velocity_tor__helium_3()
    var"deuterium" :: summary__local__magnetic_axis__velocity_tor__deuterium = summary__local__magnetic_axis__velocity_tor__deuterium()
    var"iron" :: summary__local__magnetic_axis__velocity_tor__iron = summary__local__magnetic_axis__velocity_tor__iron()
    var"helium_4" :: summary__local__magnetic_axis__velocity_tor__helium_4 = summary__local__magnetic_axis__velocity_tor__helium_4()
    var"oxygen" :: summary__local__magnetic_axis__velocity_tor__oxygen = summary__local__magnetic_axis__velocity_tor__oxygen()
    var"tungsten" :: summary__local__magnetic_axis__velocity_tor__tungsten = summary__local__magnetic_axis__velocity_tor__tungsten()
    var"xenon" :: summary__local__magnetic_axis__velocity_tor__xenon = summary__local__magnetic_axis__velocity_tor__xenon()
    var"hydrogen" :: summary__local__magnetic_axis__velocity_tor__hydrogen = summary__local__magnetic_axis__velocity_tor__hydrogen()
    var"carbon" :: summary__local__magnetic_axis__velocity_tor__carbon = summary__local__magnetic_axis__velocity_tor__carbon()
    var"nitrogen" :: summary__local__magnetic_axis__velocity_tor__nitrogen = summary__local__magnetic_axis__velocity_tor__nitrogen()
    var"beryllium" :: summary__local__magnetic_axis__velocity_tor__beryllium = summary__local__magnetic_axis__velocity_tor__beryllium()
    var"argon" :: summary__local__magnetic_axis__velocity_tor__argon = summary__local__magnetic_axis__velocity_tor__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__velocity_tor(var"krypton"=summary__local__magnetic_axis__velocity_tor__krypton(), var"lithium"=summary__local__magnetic_axis__velocity_tor__lithium(), var"neon"=summary__local__magnetic_axis__velocity_tor__neon(), var"tritium"=summary__local__magnetic_axis__velocity_tor__tritium(), var"helium_3"=summary__local__magnetic_axis__velocity_tor__helium_3(), var"deuterium"=summary__local__magnetic_axis__velocity_tor__deuterium(), var"iron"=summary__local__magnetic_axis__velocity_tor__iron(), var"helium_4"=summary__local__magnetic_axis__velocity_tor__helium_4(), var"oxygen"=summary__local__magnetic_axis__velocity_tor__oxygen(), var"tungsten"=summary__local__magnetic_axis__velocity_tor__tungsten(), var"xenon"=summary__local__magnetic_axis__velocity_tor__xenon(), var"hydrogen"=summary__local__magnetic_axis__velocity_tor__hydrogen(), var"carbon"=summary__local__magnetic_axis__velocity_tor__carbon(), var"nitrogen"=summary__local__magnetic_axis__velocity_tor__nitrogen(), var"beryllium"=summary__local__magnetic_axis__velocity_tor__beryllium(), var"argon"=summary__local__magnetic_axis__velocity_tor__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__q <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__position <: FDS
    var"psi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__position(var"psi"=missing, var"r"=missing, var"rho_tor_norm"=missing, var"z"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"r", var"rho_tor_norm", var"z", var"rho_tor", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_i <: FDS
    var"krypton" :: summary__local__magnetic_axis__n_i__krypton = summary__local__magnetic_axis__n_i__krypton()
    var"lithium" :: summary__local__magnetic_axis__n_i__lithium = summary__local__magnetic_axis__n_i__lithium()
    var"neon" :: summary__local__magnetic_axis__n_i__neon = summary__local__magnetic_axis__n_i__neon()
    var"tritium" :: summary__local__magnetic_axis__n_i__tritium = summary__local__magnetic_axis__n_i__tritium()
    var"helium_3" :: summary__local__magnetic_axis__n_i__helium_3 = summary__local__magnetic_axis__n_i__helium_3()
    var"deuterium" :: summary__local__magnetic_axis__n_i__deuterium = summary__local__magnetic_axis__n_i__deuterium()
    var"iron" :: summary__local__magnetic_axis__n_i__iron = summary__local__magnetic_axis__n_i__iron()
    var"helium_4" :: summary__local__magnetic_axis__n_i__helium_4 = summary__local__magnetic_axis__n_i__helium_4()
    var"oxygen" :: summary__local__magnetic_axis__n_i__oxygen = summary__local__magnetic_axis__n_i__oxygen()
    var"tungsten" :: summary__local__magnetic_axis__n_i__tungsten = summary__local__magnetic_axis__n_i__tungsten()
    var"xenon" :: summary__local__magnetic_axis__n_i__xenon = summary__local__magnetic_axis__n_i__xenon()
    var"hydrogen" :: summary__local__magnetic_axis__n_i__hydrogen = summary__local__magnetic_axis__n_i__hydrogen()
    var"carbon" :: summary__local__magnetic_axis__n_i__carbon = summary__local__magnetic_axis__n_i__carbon()
    var"nitrogen" :: summary__local__magnetic_axis__n_i__nitrogen = summary__local__magnetic_axis__n_i__nitrogen()
    var"beryllium" :: summary__local__magnetic_axis__n_i__beryllium = summary__local__magnetic_axis__n_i__beryllium()
    var"argon" :: summary__local__magnetic_axis__n_i__argon = summary__local__magnetic_axis__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_i(var"krypton"=summary__local__magnetic_axis__n_i__krypton(), var"lithium"=summary__local__magnetic_axis__n_i__lithium(), var"neon"=summary__local__magnetic_axis__n_i__neon(), var"tritium"=summary__local__magnetic_axis__n_i__tritium(), var"helium_3"=summary__local__magnetic_axis__n_i__helium_3(), var"deuterium"=summary__local__magnetic_axis__n_i__deuterium(), var"iron"=summary__local__magnetic_axis__n_i__iron(), var"helium_4"=summary__local__magnetic_axis__n_i__helium_4(), var"oxygen"=summary__local__magnetic_axis__n_i__oxygen(), var"tungsten"=summary__local__magnetic_axis__n_i__tungsten(), var"xenon"=summary__local__magnetic_axis__n_i__xenon(), var"hydrogen"=summary__local__magnetic_axis__n_i__hydrogen(), var"carbon"=summary__local__magnetic_axis__n_i__carbon(), var"nitrogen"=summary__local__magnetic_axis__n_i__nitrogen(), var"beryllium"=summary__local__magnetic_axis__n_i__beryllium(), var"argon"=summary__local__magnetic_axis__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__momentum_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__magnetic_shear <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__e_field_parallel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis__b_field <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis__b_field(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__magnetic_axis <: FDS
    var"n_i" :: summary__local__magnetic_axis__n_i = summary__local__magnetic_axis__n_i()
    var"velocity_tor" :: summary__local__magnetic_axis__velocity_tor = summary__local__magnetic_axis__velocity_tor()
    var"magnetic_shear" :: summary__local__magnetic_axis__magnetic_shear = summary__local__magnetic_axis__magnetic_shear()
    var"t_e" :: summary__local__magnetic_axis__t_e = summary__local__magnetic_axis__t_e()
    var"t_i_average" :: summary__local__magnetic_axis__t_i_average = summary__local__magnetic_axis__t_i_average()
    var"b_field" :: summary__local__magnetic_axis__b_field = summary__local__magnetic_axis__b_field()
    var"q" :: summary__local__magnetic_axis__q = summary__local__magnetic_axis__q()
    var"e_field_parallel" :: summary__local__magnetic_axis__e_field_parallel = summary__local__magnetic_axis__e_field_parallel()
    var"momentum_tor" :: summary__local__magnetic_axis__momentum_tor = summary__local__magnetic_axis__momentum_tor()
    var"n_i_total" :: summary__local__magnetic_axis__n_i_total = summary__local__magnetic_axis__n_i_total()
    var"position" :: summary__local__magnetic_axis__position = summary__local__magnetic_axis__position()
    var"n_e" :: summary__local__magnetic_axis__n_e = summary__local__magnetic_axis__n_e()
    var"zeff" :: summary__local__magnetic_axis__zeff = summary__local__magnetic_axis__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__magnetic_axis(var"n_i"=summary__local__magnetic_axis__n_i(), var"velocity_tor"=summary__local__magnetic_axis__velocity_tor(), var"magnetic_shear"=summary__local__magnetic_axis__magnetic_shear(), var"t_e"=summary__local__magnetic_axis__t_e(), var"t_i_average"=summary__local__magnetic_axis__t_i_average(), var"b_field"=summary__local__magnetic_axis__b_field(), var"q"=summary__local__magnetic_axis__q(), var"e_field_parallel"=summary__local__magnetic_axis__e_field_parallel(), var"momentum_tor"=summary__local__magnetic_axis__momentum_tor(), var"n_i_total"=summary__local__magnetic_axis__n_i_total(), var"position"=summary__local__magnetic_axis__position(), var"n_e"=summary__local__magnetic_axis__n_e(), var"zeff"=summary__local__magnetic_axis__zeff(), _parent=WeakRef(missing))
        fds = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"b_field", var"q", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"position", var"n_e", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.velocity_tor, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_shear, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.b_field, :_parent, WeakRef(fds))
        setfield!(fds.q, :_parent, WeakRef(fds))
        setfield!(fds.e_field_parallel, :_parent, WeakRef(fds))
        setfield!(fds.momentum_tor, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__power_flux_peak <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__power_flux_peak(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__name <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__name(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_i <: FDS
    var"krypton" :: summary__local__limiter__n_i__krypton = summary__local__limiter__n_i__krypton()
    var"lithium" :: summary__local__limiter__n_i__lithium = summary__local__limiter__n_i__lithium()
    var"neon" :: summary__local__limiter__n_i__neon = summary__local__limiter__n_i__neon()
    var"tritium" :: summary__local__limiter__n_i__tritium = summary__local__limiter__n_i__tritium()
    var"helium_3" :: summary__local__limiter__n_i__helium_3 = summary__local__limiter__n_i__helium_3()
    var"deuterium" :: summary__local__limiter__n_i__deuterium = summary__local__limiter__n_i__deuterium()
    var"iron" :: summary__local__limiter__n_i__iron = summary__local__limiter__n_i__iron()
    var"helium_4" :: summary__local__limiter__n_i__helium_4 = summary__local__limiter__n_i__helium_4()
    var"oxygen" :: summary__local__limiter__n_i__oxygen = summary__local__limiter__n_i__oxygen()
    var"tungsten" :: summary__local__limiter__n_i__tungsten = summary__local__limiter__n_i__tungsten()
    var"xenon" :: summary__local__limiter__n_i__xenon = summary__local__limiter__n_i__xenon()
    var"hydrogen" :: summary__local__limiter__n_i__hydrogen = summary__local__limiter__n_i__hydrogen()
    var"carbon" :: summary__local__limiter__n_i__carbon = summary__local__limiter__n_i__carbon()
    var"nitrogen" :: summary__local__limiter__n_i__nitrogen = summary__local__limiter__n_i__nitrogen()
    var"beryllium" :: summary__local__limiter__n_i__beryllium = summary__local__limiter__n_i__beryllium()
    var"argon" :: summary__local__limiter__n_i__argon = summary__local__limiter__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_i(var"krypton"=summary__local__limiter__n_i__krypton(), var"lithium"=summary__local__limiter__n_i__lithium(), var"neon"=summary__local__limiter__n_i__neon(), var"tritium"=summary__local__limiter__n_i__tritium(), var"helium_3"=summary__local__limiter__n_i__helium_3(), var"deuterium"=summary__local__limiter__n_i__deuterium(), var"iron"=summary__local__limiter__n_i__iron(), var"helium_4"=summary__local__limiter__n_i__helium_4(), var"oxygen"=summary__local__limiter__n_i__oxygen(), var"tungsten"=summary__local__limiter__n_i__tungsten(), var"xenon"=summary__local__limiter__n_i__xenon(), var"hydrogen"=summary__local__limiter__n_i__hydrogen(), var"carbon"=summary__local__limiter__n_i__carbon(), var"nitrogen"=summary__local__limiter__n_i__nitrogen(), var"beryllium"=summary__local__limiter__n_i__beryllium(), var"argon"=summary__local__limiter__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter__flux_expansion <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter__flux_expansion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__limiter <: FDS
    var"flux_expansion" :: summary__local__limiter__flux_expansion = summary__local__limiter__flux_expansion()
    var"n_i" :: summary__local__limiter__n_i = summary__local__limiter__n_i()
    var"name" :: summary__local__limiter__name = summary__local__limiter__name()
    var"power_flux_peak" :: summary__local__limiter__power_flux_peak = summary__local__limiter__power_flux_peak()
    var"t_e" :: summary__local__limiter__t_e = summary__local__limiter__t_e()
    var"n_e" :: summary__local__limiter__n_e = summary__local__limiter__n_e()
    var"t_i_average" :: summary__local__limiter__t_i_average = summary__local__limiter__t_i_average()
    var"n_i_total" :: summary__local__limiter__n_i_total = summary__local__limiter__n_i_total()
    var"zeff" :: summary__local__limiter__zeff = summary__local__limiter__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__limiter(var"flux_expansion"=summary__local__limiter__flux_expansion(), var"n_i"=summary__local__limiter__n_i(), var"name"=summary__local__limiter__name(), var"power_flux_peak"=summary__local__limiter__power_flux_peak(), var"t_e"=summary__local__limiter__t_e(), var"n_e"=summary__local__limiter__n_e(), var"t_i_average"=summary__local__limiter__t_i_average(), var"n_i_total"=summary__local__limiter__n_i_total(), var"zeff"=summary__local__limiter__zeff(), _parent=WeakRef(missing))
        fds = new(var"flux_expansion", var"n_i", var"name", var"power_flux_peak", var"t_e", var"n_e", var"t_i_average", var"n_i_total", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.flux_expansion, :_parent, WeakRef(fds))
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.name, :_parent, WeakRef(fds))
        setfield!(fds.power_flux_peak, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__velocity_tor <: FDS
    var"krypton" :: summary__local__itb__velocity_tor__krypton = summary__local__itb__velocity_tor__krypton()
    var"lithium" :: summary__local__itb__velocity_tor__lithium = summary__local__itb__velocity_tor__lithium()
    var"neon" :: summary__local__itb__velocity_tor__neon = summary__local__itb__velocity_tor__neon()
    var"tritium" :: summary__local__itb__velocity_tor__tritium = summary__local__itb__velocity_tor__tritium()
    var"helium_3" :: summary__local__itb__velocity_tor__helium_3 = summary__local__itb__velocity_tor__helium_3()
    var"deuterium" :: summary__local__itb__velocity_tor__deuterium = summary__local__itb__velocity_tor__deuterium()
    var"iron" :: summary__local__itb__velocity_tor__iron = summary__local__itb__velocity_tor__iron()
    var"helium_4" :: summary__local__itb__velocity_tor__helium_4 = summary__local__itb__velocity_tor__helium_4()
    var"oxygen" :: summary__local__itb__velocity_tor__oxygen = summary__local__itb__velocity_tor__oxygen()
    var"tungsten" :: summary__local__itb__velocity_tor__tungsten = summary__local__itb__velocity_tor__tungsten()
    var"xenon" :: summary__local__itb__velocity_tor__xenon = summary__local__itb__velocity_tor__xenon()
    var"hydrogen" :: summary__local__itb__velocity_tor__hydrogen = summary__local__itb__velocity_tor__hydrogen()
    var"carbon" :: summary__local__itb__velocity_tor__carbon = summary__local__itb__velocity_tor__carbon()
    var"nitrogen" :: summary__local__itb__velocity_tor__nitrogen = summary__local__itb__velocity_tor__nitrogen()
    var"beryllium" :: summary__local__itb__velocity_tor__beryllium = summary__local__itb__velocity_tor__beryllium()
    var"argon" :: summary__local__itb__velocity_tor__argon = summary__local__itb__velocity_tor__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__velocity_tor(var"krypton"=summary__local__itb__velocity_tor__krypton(), var"lithium"=summary__local__itb__velocity_tor__lithium(), var"neon"=summary__local__itb__velocity_tor__neon(), var"tritium"=summary__local__itb__velocity_tor__tritium(), var"helium_3"=summary__local__itb__velocity_tor__helium_3(), var"deuterium"=summary__local__itb__velocity_tor__deuterium(), var"iron"=summary__local__itb__velocity_tor__iron(), var"helium_4"=summary__local__itb__velocity_tor__helium_4(), var"oxygen"=summary__local__itb__velocity_tor__oxygen(), var"tungsten"=summary__local__itb__velocity_tor__tungsten(), var"xenon"=summary__local__itb__velocity_tor__xenon(), var"hydrogen"=summary__local__itb__velocity_tor__hydrogen(), var"carbon"=summary__local__itb__velocity_tor__carbon(), var"nitrogen"=summary__local__itb__velocity_tor__nitrogen(), var"beryllium"=summary__local__itb__velocity_tor__beryllium(), var"argon"=summary__local__itb__velocity_tor__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__q <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__position <: FDS
    var"psi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__position(var"psi"=missing, var"rho_tor_norm"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"rho_tor_norm", var"rho_tor", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_i <: FDS
    var"krypton" :: summary__local__itb__n_i__krypton = summary__local__itb__n_i__krypton()
    var"lithium" :: summary__local__itb__n_i__lithium = summary__local__itb__n_i__lithium()
    var"neon" :: summary__local__itb__n_i__neon = summary__local__itb__n_i__neon()
    var"tritium" :: summary__local__itb__n_i__tritium = summary__local__itb__n_i__tritium()
    var"helium_3" :: summary__local__itb__n_i__helium_3 = summary__local__itb__n_i__helium_3()
    var"deuterium" :: summary__local__itb__n_i__deuterium = summary__local__itb__n_i__deuterium()
    var"iron" :: summary__local__itb__n_i__iron = summary__local__itb__n_i__iron()
    var"helium_4" :: summary__local__itb__n_i__helium_4 = summary__local__itb__n_i__helium_4()
    var"oxygen" :: summary__local__itb__n_i__oxygen = summary__local__itb__n_i__oxygen()
    var"tungsten" :: summary__local__itb__n_i__tungsten = summary__local__itb__n_i__tungsten()
    var"xenon" :: summary__local__itb__n_i__xenon = summary__local__itb__n_i__xenon()
    var"hydrogen" :: summary__local__itb__n_i__hydrogen = summary__local__itb__n_i__hydrogen()
    var"carbon" :: summary__local__itb__n_i__carbon = summary__local__itb__n_i__carbon()
    var"nitrogen" :: summary__local__itb__n_i__nitrogen = summary__local__itb__n_i__nitrogen()
    var"beryllium" :: summary__local__itb__n_i__beryllium = summary__local__itb__n_i__beryllium()
    var"argon" :: summary__local__itb__n_i__argon = summary__local__itb__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_i(var"krypton"=summary__local__itb__n_i__krypton(), var"lithium"=summary__local__itb__n_i__lithium(), var"neon"=summary__local__itb__n_i__neon(), var"tritium"=summary__local__itb__n_i__tritium(), var"helium_3"=summary__local__itb__n_i__helium_3(), var"deuterium"=summary__local__itb__n_i__deuterium(), var"iron"=summary__local__itb__n_i__iron(), var"helium_4"=summary__local__itb__n_i__helium_4(), var"oxygen"=summary__local__itb__n_i__oxygen(), var"tungsten"=summary__local__itb__n_i__tungsten(), var"xenon"=summary__local__itb__n_i__xenon(), var"hydrogen"=summary__local__itb__n_i__hydrogen(), var"carbon"=summary__local__itb__n_i__carbon(), var"nitrogen"=summary__local__itb__n_i__nitrogen(), var"beryllium"=summary__local__itb__n_i__beryllium(), var"argon"=summary__local__itb__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__momentum_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__magnetic_shear <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb__e_field_parallel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__itb <: FDS
    var"n_i" :: summary__local__itb__n_i = summary__local__itb__n_i()
    var"velocity_tor" :: summary__local__itb__velocity_tor = summary__local__itb__velocity_tor()
    var"magnetic_shear" :: summary__local__itb__magnetic_shear = summary__local__itb__magnetic_shear()
    var"t_e" :: summary__local__itb__t_e = summary__local__itb__t_e()
    var"t_i_average" :: summary__local__itb__t_i_average = summary__local__itb__t_i_average()
    var"position" :: summary__local__itb__position = summary__local__itb__position()
    var"e_field_parallel" :: summary__local__itb__e_field_parallel = summary__local__itb__e_field_parallel()
    var"momentum_tor" :: summary__local__itb__momentum_tor = summary__local__itb__momentum_tor()
    var"n_i_total" :: summary__local__itb__n_i_total = summary__local__itb__n_i_total()
    var"q" :: summary__local__itb__q = summary__local__itb__q()
    var"n_e" :: summary__local__itb__n_e = summary__local__itb__n_e()
    var"zeff" :: summary__local__itb__zeff = summary__local__itb__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__itb(var"n_i"=summary__local__itb__n_i(), var"velocity_tor"=summary__local__itb__velocity_tor(), var"magnetic_shear"=summary__local__itb__magnetic_shear(), var"t_e"=summary__local__itb__t_e(), var"t_i_average"=summary__local__itb__t_i_average(), var"position"=summary__local__itb__position(), var"e_field_parallel"=summary__local__itb__e_field_parallel(), var"momentum_tor"=summary__local__itb__momentum_tor(), var"n_i_total"=summary__local__itb__n_i_total(), var"q"=summary__local__itb__q(), var"n_e"=summary__local__itb__n_e(), var"zeff"=summary__local__itb__zeff(), _parent=WeakRef(missing))
        fds = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"position", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"q", var"n_e", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.velocity_tor, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_shear, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.e_field_parallel, :_parent, WeakRef(fds))
        setfield!(fds.momentum_tor, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.q, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___power_flux_peak <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___power_flux_peak(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___name <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___name(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_i <: FDS
    var"krypton" :: summary__local__divertor_plate___n_i__krypton = summary__local__divertor_plate___n_i__krypton()
    var"lithium" :: summary__local__divertor_plate___n_i__lithium = summary__local__divertor_plate___n_i__lithium()
    var"neon" :: summary__local__divertor_plate___n_i__neon = summary__local__divertor_plate___n_i__neon()
    var"tritium" :: summary__local__divertor_plate___n_i__tritium = summary__local__divertor_plate___n_i__tritium()
    var"helium_3" :: summary__local__divertor_plate___n_i__helium_3 = summary__local__divertor_plate___n_i__helium_3()
    var"deuterium" :: summary__local__divertor_plate___n_i__deuterium = summary__local__divertor_plate___n_i__deuterium()
    var"iron" :: summary__local__divertor_plate___n_i__iron = summary__local__divertor_plate___n_i__iron()
    var"helium_4" :: summary__local__divertor_plate___n_i__helium_4 = summary__local__divertor_plate___n_i__helium_4()
    var"oxygen" :: summary__local__divertor_plate___n_i__oxygen = summary__local__divertor_plate___n_i__oxygen()
    var"tungsten" :: summary__local__divertor_plate___n_i__tungsten = summary__local__divertor_plate___n_i__tungsten()
    var"xenon" :: summary__local__divertor_plate___n_i__xenon = summary__local__divertor_plate___n_i__xenon()
    var"hydrogen" :: summary__local__divertor_plate___n_i__hydrogen = summary__local__divertor_plate___n_i__hydrogen()
    var"carbon" :: summary__local__divertor_plate___n_i__carbon = summary__local__divertor_plate___n_i__carbon()
    var"nitrogen" :: summary__local__divertor_plate___n_i__nitrogen = summary__local__divertor_plate___n_i__nitrogen()
    var"beryllium" :: summary__local__divertor_plate___n_i__beryllium = summary__local__divertor_plate___n_i__beryllium()
    var"argon" :: summary__local__divertor_plate___n_i__argon = summary__local__divertor_plate___n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_i(var"krypton"=summary__local__divertor_plate___n_i__krypton(), var"lithium"=summary__local__divertor_plate___n_i__lithium(), var"neon"=summary__local__divertor_plate___n_i__neon(), var"tritium"=summary__local__divertor_plate___n_i__tritium(), var"helium_3"=summary__local__divertor_plate___n_i__helium_3(), var"deuterium"=summary__local__divertor_plate___n_i__deuterium(), var"iron"=summary__local__divertor_plate___n_i__iron(), var"helium_4"=summary__local__divertor_plate___n_i__helium_4(), var"oxygen"=summary__local__divertor_plate___n_i__oxygen(), var"tungsten"=summary__local__divertor_plate___n_i__tungsten(), var"xenon"=summary__local__divertor_plate___n_i__xenon(), var"hydrogen"=summary__local__divertor_plate___n_i__hydrogen(), var"carbon"=summary__local__divertor_plate___n_i__carbon(), var"nitrogen"=summary__local__divertor_plate___n_i__nitrogen(), var"beryllium"=summary__local__divertor_plate___n_i__beryllium(), var"argon"=summary__local__divertor_plate___n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate___flux_expansion <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate___flux_expansion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__local__divertor_plate <: FDSvectorElement
    var"flux_expansion" :: summary__local__divertor_plate___flux_expansion = summary__local__divertor_plate___flux_expansion()
    var"n_i" :: summary__local__divertor_plate___n_i = summary__local__divertor_plate___n_i()
    var"name" :: summary__local__divertor_plate___name = summary__local__divertor_plate___name()
    var"power_flux_peak" :: summary__local__divertor_plate___power_flux_peak = summary__local__divertor_plate___power_flux_peak()
    var"t_e" :: summary__local__divertor_plate___t_e = summary__local__divertor_plate___t_e()
    var"n_e" :: summary__local__divertor_plate___n_e = summary__local__divertor_plate___n_e()
    var"t_i_average" :: summary__local__divertor_plate___t_i_average = summary__local__divertor_plate___t_i_average()
    var"n_i_total" :: summary__local__divertor_plate___n_i_total = summary__local__divertor_plate___n_i_total()
    var"zeff" :: summary__local__divertor_plate___zeff = summary__local__divertor_plate___zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local__divertor_plate(var"flux_expansion"=summary__local__divertor_plate___flux_expansion(), var"n_i"=summary__local__divertor_plate___n_i(), var"name"=summary__local__divertor_plate___name(), var"power_flux_peak"=summary__local__divertor_plate___power_flux_peak(), var"t_e"=summary__local__divertor_plate___t_e(), var"n_e"=summary__local__divertor_plate___n_e(), var"t_i_average"=summary__local__divertor_plate___t_i_average(), var"n_i_total"=summary__local__divertor_plate___n_i_total(), var"zeff"=summary__local__divertor_plate___zeff(), _parent=WeakRef(missing))
        fds = new(var"flux_expansion", var"n_i", var"name", var"power_flux_peak", var"t_e", var"n_e", var"t_i_average", var"n_i_total", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.flux_expansion, :_parent, WeakRef(fds))
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.name, :_parent, WeakRef(fds))
        setfield!(fds.power_flux_peak, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__local <: FDS
    var"itb" :: summary__local__itb = summary__local__itb()
    var"r_eff_norm_2_3" :: summary__local__r_eff_norm_2_3 = summary__local__r_eff_norm_2_3()
    var"limiter" :: summary__local__limiter = summary__local__limiter()
    var"divertor_plate" :: FDSvector{T} where {T<:summary__local__divertor_plate} = FDSvector(summary__local__divertor_plate[])
    var"magnetic_axis" :: summary__local__magnetic_axis = summary__local__magnetic_axis()
    var"pedestal" :: summary__local__pedestal = summary__local__pedestal()
    var"separatrix" :: summary__local__separatrix = summary__local__separatrix()
    _parent :: WeakRef = WeakRef(missing)
    function summary__local(var"itb"=summary__local__itb(), var"r_eff_norm_2_3"=summary__local__r_eff_norm_2_3(), var"limiter"=summary__local__limiter(), var"divertor_plate"=FDSvector(summary__local__divertor_plate[]), var"magnetic_axis"=summary__local__magnetic_axis(), var"pedestal"=summary__local__pedestal(), var"separatrix"=summary__local__separatrix(), _parent=WeakRef(missing))
        fds = new(var"itb", var"r_eff_norm_2_3", var"limiter", var"divertor_plate", var"magnetic_axis", var"pedestal", var"separatrix", _parent)
        assign_expressions(fds)
        setfield!(fds.itb, :_parent, WeakRef(fds))
        setfield!(fds.r_eff_norm_2_3, :_parent, WeakRef(fds))
        setfield!(fds.limiter, :_parent, WeakRef(fds))
        setfield!(fds.divertor_plate, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_axis, :_parent, WeakRef(fds))
        setfield!(fds.pedestal, :_parent, WeakRef(fds))
        setfield!(fds.separatrix, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__zeff <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__t_i_average <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__t_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__tungsten <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__iron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_i <: FDS
    var"krypton" :: summary__line_average__n_i__krypton = summary__line_average__n_i__krypton()
    var"lithium" :: summary__line_average__n_i__lithium = summary__line_average__n_i__lithium()
    var"neon" :: summary__line_average__n_i__neon = summary__line_average__n_i__neon()
    var"tritium" :: summary__line_average__n_i__tritium = summary__line_average__n_i__tritium()
    var"helium_3" :: summary__line_average__n_i__helium_3 = summary__line_average__n_i__helium_3()
    var"deuterium" :: summary__line_average__n_i__deuterium = summary__line_average__n_i__deuterium()
    var"iron" :: summary__line_average__n_i__iron = summary__line_average__n_i__iron()
    var"helium_4" :: summary__line_average__n_i__helium_4 = summary__line_average__n_i__helium_4()
    var"oxygen" :: summary__line_average__n_i__oxygen = summary__line_average__n_i__oxygen()
    var"tungsten" :: summary__line_average__n_i__tungsten = summary__line_average__n_i__tungsten()
    var"xenon" :: summary__line_average__n_i__xenon = summary__line_average__n_i__xenon()
    var"hydrogen" :: summary__line_average__n_i__hydrogen = summary__line_average__n_i__hydrogen()
    var"carbon" :: summary__line_average__n_i__carbon = summary__line_average__n_i__carbon()
    var"nitrogen" :: summary__line_average__n_i__nitrogen = summary__line_average__n_i__nitrogen()
    var"beryllium" :: summary__line_average__n_i__beryllium = summary__line_average__n_i__beryllium()
    var"argon" :: summary__line_average__n_i__argon = summary__line_average__n_i__argon()
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_i(var"krypton"=summary__line_average__n_i__krypton(), var"lithium"=summary__line_average__n_i__lithium(), var"neon"=summary__line_average__n_i__neon(), var"tritium"=summary__line_average__n_i__tritium(), var"helium_3"=summary__line_average__n_i__helium_3(), var"deuterium"=summary__line_average__n_i__deuterium(), var"iron"=summary__line_average__n_i__iron(), var"helium_4"=summary__line_average__n_i__helium_4(), var"oxygen"=summary__line_average__n_i__oxygen(), var"tungsten"=summary__line_average__n_i__tungsten(), var"xenon"=summary__line_average__n_i__xenon(), var"hydrogen"=summary__line_average__n_i__hydrogen(), var"carbon"=summary__line_average__n_i__carbon(), var"nitrogen"=summary__line_average__n_i__nitrogen(), var"beryllium"=summary__line_average__n_i__beryllium(), var"argon"=summary__line_average__n_i__argon(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.iron, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.tungsten, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__n_e <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__meff_hydrogenic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__meff_hydrogenic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__isotope_fraction_hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__isotope_fraction_hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average__dn_e_dt <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average__dn_e_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__line_average <: FDS
    var"meff_hydrogenic" :: summary__line_average__meff_hydrogenic = summary__line_average__meff_hydrogenic()
    var"n_i" :: summary__line_average__n_i = summary__line_average__n_i()
    var"t_e" :: summary__line_average__t_e = summary__line_average__t_e()
    var"n_e" :: summary__line_average__n_e = summary__line_average__n_e()
    var"isotope_fraction_hydrogen" :: summary__line_average__isotope_fraction_hydrogen = summary__line_average__isotope_fraction_hydrogen()
    var"t_i_average" :: summary__line_average__t_i_average = summary__line_average__t_i_average()
    var"n_i_total" :: summary__line_average__n_i_total = summary__line_average__n_i_total()
    var"dn_e_dt" :: summary__line_average__dn_e_dt = summary__line_average__dn_e_dt()
    var"zeff" :: summary__line_average__zeff = summary__line_average__zeff()
    _parent :: WeakRef = WeakRef(missing)
    function summary__line_average(var"meff_hydrogenic"=summary__line_average__meff_hydrogenic(), var"n_i"=summary__line_average__n_i(), var"t_e"=summary__line_average__t_e(), var"n_e"=summary__line_average__n_e(), var"isotope_fraction_hydrogen"=summary__line_average__isotope_fraction_hydrogen(), var"t_i_average"=summary__line_average__t_i_average(), var"n_i_total"=summary__line_average__n_i_total(), var"dn_e_dt"=summary__line_average__dn_e_dt(), var"zeff"=summary__line_average__zeff(), _parent=WeakRef(missing))
        fds = new(var"meff_hydrogenic", var"n_i", var"t_e", var"n_e", var"isotope_fraction_hydrogen", var"t_i_average", var"n_i_total", var"dn_e_dt", var"zeff", _parent)
        assign_expressions(fds)
        setfield!(fds.meff_hydrogenic, :_parent, WeakRef(fds))
        setfield!(fds.n_i, :_parent, WeakRef(fds))
        setfield!(fds.t_e, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.isotope_fraction_hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.t_i_average, :_parent, WeakRef(fds))
        setfield!(fds.n_i_total, :_parent, WeakRef(fds))
        setfield!(fds.dn_e_dt, :_parent, WeakRef(fds))
        setfield!(fds.zeff, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__limiter__material <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__limiter__material(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__limiter <: FDS
    var"material" :: summary__limiter__material = summary__limiter__material()
    _parent :: WeakRef = WeakRef(missing)
    function summary__limiter(var"material"=summary__limiter__material(), _parent=WeakRef(missing))
        fds = new(var"material", _parent)
        assign_expressions(fds)
        setfield!(fds.material, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__kicks__occurrence <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__kicks__occurrence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__kicks <: FDS
    var"occurrence" :: summary__kicks__occurrence = summary__kicks__occurrence()
    _parent :: WeakRef = WeakRef(missing)
    function summary__kicks(var"occurrence"=summary__kicks__occurrence(), _parent=WeakRef(missing))
        fds = new(var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.occurrence, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String, Function} = missing
    var"data_dictionary" :: Union{Missing, String, Function} = missing
    var"access_layer" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__ids_properties <: FDS
    var"provider" :: Union{Missing, String, Function} = missing
    var"version_put" :: summary__ids_properties__version_put = summary__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Integer, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"creation_date" :: Union{Missing, String, Function} = missing
    var"comment" :: Union{Missing, String, Function} = missing
    var"occurrence" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__ids_properties(var"provider"=missing, var"version_put"=summary__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_nbi <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_nbi(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_lh <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_lh(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_nbi_co_injected_ratio <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_launched_nbi_co_injected_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_nbi <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_launched_nbi(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_lh <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_launched_lh(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_ic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_launched_ic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_launched_ec <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_launched_ec(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_ic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_ic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_ec <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_ec(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__power_additional <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__power_additional(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___tangency_radius <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___tangency_radius(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___species__z_n <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___species__z_n(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___species__label <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___species__label(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___species__a <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___species__a(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___species <: FDS
    var"label" :: summary__heating_current_drive__nbi___species__label = summary__heating_current_drive__nbi___species__label()
    var"z_n" :: summary__heating_current_drive__nbi___species__z_n = summary__heating_current_drive__nbi___species__z_n()
    var"a" :: summary__heating_current_drive__nbi___species__a = summary__heating_current_drive__nbi___species__a()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___species(var"label"=summary__heating_current_drive__nbi___species__label(), var"z_n"=summary__heating_current_drive__nbi___species__z_n(), var"a"=summary__heating_current_drive__nbi___species__a(), _parent=WeakRef(missing))
        fds = new(var"label", var"z_n", var"a", _parent)
        assign_expressions(fds)
        setfield!(fds.label, :_parent, WeakRef(fds))
        setfield!(fds.z_n, :_parent, WeakRef(fds))
        setfield!(fds.a, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___power_launched <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___power <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___position__z <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___position__z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___position__r <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___position__r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___position__phi <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___position__phi(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___position <: FDS
    var"phi" :: summary__heating_current_drive__nbi___position__phi = summary__heating_current_drive__nbi___position__phi()
    var"r" :: summary__heating_current_drive__nbi___position__r = summary__heating_current_drive__nbi___position__r()
    var"z" :: summary__heating_current_drive__nbi___position__z = summary__heating_current_drive__nbi___position__z()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___position(var"phi"=summary__heating_current_drive__nbi___position__phi(), var"r"=summary__heating_current_drive__nbi___position__r(), var"z"=summary__heating_current_drive__nbi___position__z(), _parent=WeakRef(missing))
        fds = new(var"phi", var"r", var"z", _parent)
        assign_expressions(fds)
        setfield!(fds.phi, :_parent, WeakRef(fds))
        setfield!(fds.r, :_parent, WeakRef(fds))
        setfield!(fds.z, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___energy <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___energy(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___direction <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___direction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___current <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___beam_power_fraction <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___beam_power_fraction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___beam_current_fraction <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___beam_current_fraction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi___angle <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi___angle(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__nbi <: FDSvectorElement
    var"energy" :: summary__heating_current_drive__nbi___energy = summary__heating_current_drive__nbi___energy()
    var"angle" :: summary__heating_current_drive__nbi___angle = summary__heating_current_drive__nbi___angle()
    var"current" :: summary__heating_current_drive__nbi___current = summary__heating_current_drive__nbi___current()
    var"beam_current_fraction" :: summary__heating_current_drive__nbi___beam_current_fraction = summary__heating_current_drive__nbi___beam_current_fraction()
    var"power_launched" :: summary__heating_current_drive__nbi___power_launched = summary__heating_current_drive__nbi___power_launched()
    var"position" :: summary__heating_current_drive__nbi___position = summary__heating_current_drive__nbi___position()
    var"tangency_radius" :: summary__heating_current_drive__nbi___tangency_radius = summary__heating_current_drive__nbi___tangency_radius()
    var"power" :: summary__heating_current_drive__nbi___power = summary__heating_current_drive__nbi___power()
    var"direction" :: summary__heating_current_drive__nbi___direction = summary__heating_current_drive__nbi___direction()
    var"beam_power_fraction" :: summary__heating_current_drive__nbi___beam_power_fraction = summary__heating_current_drive__nbi___beam_power_fraction()
    var"species" :: summary__heating_current_drive__nbi___species = summary__heating_current_drive__nbi___species()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__nbi(var"energy"=summary__heating_current_drive__nbi___energy(), var"angle"=summary__heating_current_drive__nbi___angle(), var"current"=summary__heating_current_drive__nbi___current(), var"beam_current_fraction"=summary__heating_current_drive__nbi___beam_current_fraction(), var"power_launched"=summary__heating_current_drive__nbi___power_launched(), var"position"=summary__heating_current_drive__nbi___position(), var"tangency_radius"=summary__heating_current_drive__nbi___tangency_radius(), var"power"=summary__heating_current_drive__nbi___power(), var"direction"=summary__heating_current_drive__nbi___direction(), var"beam_power_fraction"=summary__heating_current_drive__nbi___beam_power_fraction(), var"species"=summary__heating_current_drive__nbi___species(), _parent=WeakRef(missing))
        fds = new(var"energy", var"angle", var"current", var"beam_current_fraction", var"power_launched", var"position", var"tangency_radius", var"power", var"direction", var"beam_power_fraction", var"species", _parent)
        assign_expressions(fds)
        setfield!(fds.energy, :_parent, WeakRef(fds))
        setfield!(fds.angle, :_parent, WeakRef(fds))
        setfield!(fds.current, :_parent, WeakRef(fds))
        setfield!(fds.beam_current_fraction, :_parent, WeakRef(fds))
        setfield!(fds.power_launched, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.tangency_radius, :_parent, WeakRef(fds))
        setfield!(fds.power, :_parent, WeakRef(fds))
        setfield!(fds.direction, :_parent, WeakRef(fds))
        setfield!(fds.beam_power_fraction, :_parent, WeakRef(fds))
        setfield!(fds.species, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___power_launched <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___power <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___n_parallel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___n_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___frequency <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___energy_fast <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___energy_fast(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh___current <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__lh <: FDSvectorElement
    var"power_launched" :: summary__heating_current_drive__lh___power_launched = summary__heating_current_drive__lh___power_launched()
    var"energy_fast" :: summary__heating_current_drive__lh___energy_fast = summary__heating_current_drive__lh___energy_fast()
    var"frequency" :: summary__heating_current_drive__lh___frequency = summary__heating_current_drive__lh___frequency()
    var"position" :: summary__heating_current_drive__lh___position = summary__heating_current_drive__lh___position()
    var"power" :: summary__heating_current_drive__lh___power = summary__heating_current_drive__lh___power()
    var"n_parallel" :: summary__heating_current_drive__lh___n_parallel = summary__heating_current_drive__lh___n_parallel()
    var"current" :: summary__heating_current_drive__lh___current = summary__heating_current_drive__lh___current()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__lh(var"power_launched"=summary__heating_current_drive__lh___power_launched(), var"energy_fast"=summary__heating_current_drive__lh___energy_fast(), var"frequency"=summary__heating_current_drive__lh___frequency(), var"position"=summary__heating_current_drive__lh___position(), var"power"=summary__heating_current_drive__lh___power(), var"n_parallel"=summary__heating_current_drive__lh___n_parallel(), var"current"=summary__heating_current_drive__lh___current(), _parent=WeakRef(missing))
        fds = new(var"power_launched", var"energy_fast", var"frequency", var"position", var"power", var"n_parallel", var"current", _parent)
        assign_expressions(fds)
        setfield!(fds.power_launched, :_parent, WeakRef(fds))
        setfield!(fds.energy_fast, :_parent, WeakRef(fds))
        setfield!(fds.frequency, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.power, :_parent, WeakRef(fds))
        setfield!(fds.n_parallel, :_parent, WeakRef(fds))
        setfield!(fds.current, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___power_launched <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___power <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___phase <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___phase(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___n_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___n_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___k_perpendicular <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___k_perpendicular(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___harmonic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___harmonic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___frequency <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___energy_fast <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___energy_fast(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___e_field_plus_minus_ratio <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___e_field_plus_minus_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic___current <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ic <: FDSvectorElement
    var"harmonic" :: summary__heating_current_drive__ic___harmonic = summary__heating_current_drive__ic___harmonic()
    var"e_field_plus_minus_ratio" :: summary__heating_current_drive__ic___e_field_plus_minus_ratio = summary__heating_current_drive__ic___e_field_plus_minus_ratio()
    var"current" :: summary__heating_current_drive__ic___current = summary__heating_current_drive__ic___current()
    var"n_tor" :: summary__heating_current_drive__ic___n_tor = summary__heating_current_drive__ic___n_tor()
    var"power_launched" :: summary__heating_current_drive__ic___power_launched = summary__heating_current_drive__ic___power_launched()
    var"energy_fast" :: summary__heating_current_drive__ic___energy_fast = summary__heating_current_drive__ic___energy_fast()
    var"position" :: summary__heating_current_drive__ic___position = summary__heating_current_drive__ic___position()
    var"power" :: summary__heating_current_drive__ic___power = summary__heating_current_drive__ic___power()
    var"phase" :: summary__heating_current_drive__ic___phase = summary__heating_current_drive__ic___phase()
    var"frequency" :: summary__heating_current_drive__ic___frequency = summary__heating_current_drive__ic___frequency()
    var"k_perpendicular" :: summary__heating_current_drive__ic___k_perpendicular = summary__heating_current_drive__ic___k_perpendicular()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ic(var"harmonic"=summary__heating_current_drive__ic___harmonic(), var"e_field_plus_minus_ratio"=summary__heating_current_drive__ic___e_field_plus_minus_ratio(), var"current"=summary__heating_current_drive__ic___current(), var"n_tor"=summary__heating_current_drive__ic___n_tor(), var"power_launched"=summary__heating_current_drive__ic___power_launched(), var"energy_fast"=summary__heating_current_drive__ic___energy_fast(), var"position"=summary__heating_current_drive__ic___position(), var"power"=summary__heating_current_drive__ic___power(), var"phase"=summary__heating_current_drive__ic___phase(), var"frequency"=summary__heating_current_drive__ic___frequency(), var"k_perpendicular"=summary__heating_current_drive__ic___k_perpendicular(), _parent=WeakRef(missing))
        fds = new(var"harmonic", var"e_field_plus_minus_ratio", var"current", var"n_tor", var"power_launched", var"energy_fast", var"position", var"power", var"phase", var"frequency", var"k_perpendicular", _parent)
        assign_expressions(fds)
        setfield!(fds.harmonic, :_parent, WeakRef(fds))
        setfield!(fds.e_field_plus_minus_ratio, :_parent, WeakRef(fds))
        setfield!(fds.current, :_parent, WeakRef(fds))
        setfield!(fds.n_tor, :_parent, WeakRef(fds))
        setfield!(fds.power_launched, :_parent, WeakRef(fds))
        setfield!(fds.energy_fast, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.power, :_parent, WeakRef(fds))
        setfield!(fds.phase, :_parent, WeakRef(fds))
        setfield!(fds.frequency, :_parent, WeakRef(fds))
        setfield!(fds.k_perpendicular, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___power_launched <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___power <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___position <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___polarisation <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___polarisation(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___harmonic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___harmonic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___frequency <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___energy_fast <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___energy_fast(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___current <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___angle_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___angle_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec___angle_pol <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec___angle_pol(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive__ec <: FDSvectorElement
    var"angle_tor" :: summary__heating_current_drive__ec___angle_tor = summary__heating_current_drive__ec___angle_tor()
    var"polarisation" :: summary__heating_current_drive__ec___polarisation = summary__heating_current_drive__ec___polarisation()
    var"power_launched" :: summary__heating_current_drive__ec___power_launched = summary__heating_current_drive__ec___power_launched()
    var"harmonic" :: summary__heating_current_drive__ec___harmonic = summary__heating_current_drive__ec___harmonic()
    var"energy_fast" :: summary__heating_current_drive__ec___energy_fast = summary__heating_current_drive__ec___energy_fast()
    var"frequency" :: summary__heating_current_drive__ec___frequency = summary__heating_current_drive__ec___frequency()
    var"position" :: summary__heating_current_drive__ec___position = summary__heating_current_drive__ec___position()
    var"power" :: summary__heating_current_drive__ec___power = summary__heating_current_drive__ec___power()
    var"angle_pol" :: summary__heating_current_drive__ec___angle_pol = summary__heating_current_drive__ec___angle_pol()
    var"current" :: summary__heating_current_drive__ec___current = summary__heating_current_drive__ec___current()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive__ec(var"angle_tor"=summary__heating_current_drive__ec___angle_tor(), var"polarisation"=summary__heating_current_drive__ec___polarisation(), var"power_launched"=summary__heating_current_drive__ec___power_launched(), var"harmonic"=summary__heating_current_drive__ec___harmonic(), var"energy_fast"=summary__heating_current_drive__ec___energy_fast(), var"frequency"=summary__heating_current_drive__ec___frequency(), var"position"=summary__heating_current_drive__ec___position(), var"power"=summary__heating_current_drive__ec___power(), var"angle_pol"=summary__heating_current_drive__ec___angle_pol(), var"current"=summary__heating_current_drive__ec___current(), _parent=WeakRef(missing))
        fds = new(var"angle_tor", var"polarisation", var"power_launched", var"harmonic", var"energy_fast", var"frequency", var"position", var"power", var"angle_pol", var"current", _parent)
        assign_expressions(fds)
        setfield!(fds.angle_tor, :_parent, WeakRef(fds))
        setfield!(fds.polarisation, :_parent, WeakRef(fds))
        setfield!(fds.power_launched, :_parent, WeakRef(fds))
        setfield!(fds.harmonic, :_parent, WeakRef(fds))
        setfield!(fds.energy_fast, :_parent, WeakRef(fds))
        setfield!(fds.frequency, :_parent, WeakRef(fds))
        setfield!(fds.position, :_parent, WeakRef(fds))
        setfield!(fds.power, :_parent, WeakRef(fds))
        setfield!(fds.angle_pol, :_parent, WeakRef(fds))
        setfield!(fds.current, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__heating_current_drive <: FDS
    var"power_launched_ec" :: summary__heating_current_drive__power_launched_ec = summary__heating_current_drive__power_launched_ec()
    var"ec" :: FDSvector{T} where {T<:summary__heating_current_drive__ec} = FDSvector(summary__heating_current_drive__ec[])
    var"power_nbi" :: summary__heating_current_drive__power_nbi = summary__heating_current_drive__power_nbi()
    var"power_ec" :: summary__heating_current_drive__power_ec = summary__heating_current_drive__power_ec()
    var"power_lh" :: summary__heating_current_drive__power_lh = summary__heating_current_drive__power_lh()
    var"ic" :: FDSvector{T} where {T<:summary__heating_current_drive__ic} = FDSvector(summary__heating_current_drive__ic[])
    var"lh" :: FDSvector{T} where {T<:summary__heating_current_drive__lh} = FDSvector(summary__heating_current_drive__lh[])
    var"power_ic" :: summary__heating_current_drive__power_ic = summary__heating_current_drive__power_ic()
    var"power_launched_lh" :: summary__heating_current_drive__power_launched_lh = summary__heating_current_drive__power_launched_lh()
    var"nbi" :: FDSvector{T} where {T<:summary__heating_current_drive__nbi} = FDSvector(summary__heating_current_drive__nbi[])
    var"power_launched_nbi_co_injected_ratio" :: summary__heating_current_drive__power_launched_nbi_co_injected_ratio = summary__heating_current_drive__power_launched_nbi_co_injected_ratio()
    var"power_additional" :: summary__heating_current_drive__power_additional = summary__heating_current_drive__power_additional()
    var"power_launched_ic" :: summary__heating_current_drive__power_launched_ic = summary__heating_current_drive__power_launched_ic()
    var"power_launched_nbi" :: summary__heating_current_drive__power_launched_nbi = summary__heating_current_drive__power_launched_nbi()
    _parent :: WeakRef = WeakRef(missing)
    function summary__heating_current_drive(var"power_launched_ec"=summary__heating_current_drive__power_launched_ec(), var"ec"=FDSvector(summary__heating_current_drive__ec[]), var"power_nbi"=summary__heating_current_drive__power_nbi(), var"power_ec"=summary__heating_current_drive__power_ec(), var"power_lh"=summary__heating_current_drive__power_lh(), var"ic"=FDSvector(summary__heating_current_drive__ic[]), var"lh"=FDSvector(summary__heating_current_drive__lh[]), var"power_ic"=summary__heating_current_drive__power_ic(), var"power_launched_lh"=summary__heating_current_drive__power_launched_lh(), var"nbi"=FDSvector(summary__heating_current_drive__nbi[]), var"power_launched_nbi_co_injected_ratio"=summary__heating_current_drive__power_launched_nbi_co_injected_ratio(), var"power_additional"=summary__heating_current_drive__power_additional(), var"power_launched_ic"=summary__heating_current_drive__power_launched_ic(), var"power_launched_nbi"=summary__heating_current_drive__power_launched_nbi(), _parent=WeakRef(missing))
        fds = new(var"power_launched_ec", var"ec", var"power_nbi", var"power_ec", var"power_lh", var"ic", var"lh", var"power_ic", var"power_launched_lh", var"nbi", var"power_launched_nbi_co_injected_ratio", var"power_additional", var"power_launched_ic", var"power_launched_nbi", _parent)
        assign_expressions(fds)
        setfield!(fds.power_launched_ec, :_parent, WeakRef(fds))
        setfield!(fds.ec, :_parent, WeakRef(fds))
        setfield!(fds.power_nbi, :_parent, WeakRef(fds))
        setfield!(fds.power_ec, :_parent, WeakRef(fds))
        setfield!(fds.power_lh, :_parent, WeakRef(fds))
        setfield!(fds.ic, :_parent, WeakRef(fds))
        setfield!(fds.lh, :_parent, WeakRef(fds))
        setfield!(fds.power_ic, :_parent, WeakRef(fds))
        setfield!(fds.power_launched_lh, :_parent, WeakRef(fds))
        setfield!(fds.nbi, :_parent, WeakRef(fds))
        setfield!(fds.power_launched_nbi_co_injected_ratio, :_parent, WeakRef(fds))
        setfield!(fds.power_additional, :_parent, WeakRef(fds))
        setfield!(fds.power_launched_ic, :_parent, WeakRef(fds))
        setfield!(fds.power_launched_nbi, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__volume <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__volume(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__v_loop <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__v_loop(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__tau_resistive <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__tau_resistive(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__tau_helium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__tau_helium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__tau_energy_98 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__tau_energy_98(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__tau_energy <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__tau_energy(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__resistance <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__resistance(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__ratio_tau_helium_fuel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__ratio_tau_helium_fuel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__r0 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__r0(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__q_95 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__q_95(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_synchrotron <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_synchrotron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_steady <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_steady(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_radiated_outside_lcfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_radiated_outside_lcfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_radiated_inside_lcfs <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_radiated_inside_lcfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_radiated <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_radiated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_ohm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_ohm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_loss <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_loss(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_line <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_line(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__power_bremsstrahlung <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__power_bremsstrahlung(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__li_mhd <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__li_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__li <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__li(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__ip <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__ip(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__h_mode <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__h_mode(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__h_98 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__h_98(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__greenwald_fraction <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__greenwald_fraction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__fusion_gain <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__fusion_gain(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__fusion_fluence <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__fusion_fluence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_mhd <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_ion_total_thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_ion_total_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_fast_perpendicular <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_fast_perpendicular(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_fast_parallel <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_fast_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_electrons_thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_electrons_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_diamagnetic <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_diamagnetic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__energy_b_field_pol <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__energy_b_field_pol(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__denergy_thermal_dt <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__denergy_thermal_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__denergy_diamagnetic_dt <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__denergy_diamagnetic_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__current_ohm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__current_ohm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__current_non_inductive <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__current_non_inductive(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__current_bootstrap <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__current_bootstrap(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__current_alignment <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__current_alignment(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_thermal_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_tor_thermal_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_norm_mhd <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_tor_norm_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_norm <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_tor_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor_mhd <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_tor_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_tor <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_pol_mhd <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_pol_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__beta_pol <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__beta_pol(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities__b0 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities__b0(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__global_quantities <: FDS
    var"power_radiated_outside_lcfs" :: summary__global_quantities__power_radiated_outside_lcfs = summary__global_quantities__power_radiated_outside_lcfs()
    var"beta_tor" :: summary__global_quantities__beta_tor = summary__global_quantities__beta_tor()
    var"power_ohm" :: summary__global_quantities__power_ohm = summary__global_quantities__power_ohm()
    var"energy_b_field_pol" :: summary__global_quantities__energy_b_field_pol = summary__global_quantities__energy_b_field_pol()
    var"tau_energy" :: summary__global_quantities__tau_energy = summary__global_quantities__tau_energy()
    var"energy_electrons_thermal" :: summary__global_quantities__energy_electrons_thermal = summary__global_quantities__energy_electrons_thermal()
    var"greenwald_fraction" :: summary__global_quantities__greenwald_fraction = summary__global_quantities__greenwald_fraction()
    var"power_steady" :: summary__global_quantities__power_steady = summary__global_quantities__power_steady()
    var"resistance" :: summary__global_quantities__resistance = summary__global_quantities__resistance()
    var"beta_pol" :: summary__global_quantities__beta_pol = summary__global_quantities__beta_pol()
    var"energy_thermal" :: summary__global_quantities__energy_thermal = summary__global_quantities__energy_thermal()
    var"fusion_gain" :: summary__global_quantities__fusion_gain = summary__global_quantities__fusion_gain()
    var"power_radiated_inside_lcfs" :: summary__global_quantities__power_radiated_inside_lcfs = summary__global_quantities__power_radiated_inside_lcfs()
    var"beta_tor_norm" :: summary__global_quantities__beta_tor_norm = summary__global_quantities__beta_tor_norm()
    var"beta_tor_thermal_norm" :: summary__global_quantities__beta_tor_thermal_norm = summary__global_quantities__beta_tor_thermal_norm()
    var"tau_helium" :: summary__global_quantities__tau_helium = summary__global_quantities__tau_helium()
    var"power_loss" :: summary__global_quantities__power_loss = summary__global_quantities__power_loss()
    var"current_alignment" :: summary__global_quantities__current_alignment = summary__global_quantities__current_alignment()
    var"energy_diamagnetic" :: summary__global_quantities__energy_diamagnetic = summary__global_quantities__energy_diamagnetic()
    var"tau_energy_98" :: summary__global_quantities__tau_energy_98 = summary__global_quantities__tau_energy_98()
    var"li_mhd" :: summary__global_quantities__li_mhd = summary__global_quantities__li_mhd()
    var"energy_fast_perpendicular" :: summary__global_quantities__energy_fast_perpendicular = summary__global_quantities__energy_fast_perpendicular()
    var"q_95" :: summary__global_quantities__q_95 = summary__global_quantities__q_95()
    var"beta_tor_mhd" :: summary__global_quantities__beta_tor_mhd = summary__global_quantities__beta_tor_mhd()
    var"volume" :: summary__global_quantities__volume = summary__global_quantities__volume()
    var"b0" :: summary__global_quantities__b0 = summary__global_quantities__b0()
    var"ip" :: summary__global_quantities__ip = summary__global_quantities__ip()
    var"beta_pol_mhd" :: summary__global_quantities__beta_pol_mhd = summary__global_quantities__beta_pol_mhd()
    var"denergy_diamagnetic_dt" :: summary__global_quantities__denergy_diamagnetic_dt = summary__global_quantities__denergy_diamagnetic_dt()
    var"h_mode" :: summary__global_quantities__h_mode = summary__global_quantities__h_mode()
    var"power_bremsstrahlung" :: summary__global_quantities__power_bremsstrahlung = summary__global_quantities__power_bremsstrahlung()
    var"h_98" :: summary__global_quantities__h_98 = summary__global_quantities__h_98()
    var"tau_resistive" :: summary__global_quantities__tau_resistive = summary__global_quantities__tau_resistive()
    var"current_non_inductive" :: summary__global_quantities__current_non_inductive = summary__global_quantities__current_non_inductive()
    var"current_ohm" :: summary__global_quantities__current_ohm = summary__global_quantities__current_ohm()
    var"energy_total" :: summary__global_quantities__energy_total = summary__global_quantities__energy_total()
    var"power_synchrotron" :: summary__global_quantities__power_synchrotron = summary__global_quantities__power_synchrotron()
    var"v_loop" :: summary__global_quantities__v_loop = summary__global_quantities__v_loop()
    var"energy_ion_total_thermal" :: summary__global_quantities__energy_ion_total_thermal = summary__global_quantities__energy_ion_total_thermal()
    var"power_line" :: summary__global_quantities__power_line = summary__global_quantities__power_line()
    var"current_bootstrap" :: summary__global_quantities__current_bootstrap = summary__global_quantities__current_bootstrap()
    var"fusion_fluence" :: summary__global_quantities__fusion_fluence = summary__global_quantities__fusion_fluence()
    var"energy_fast_parallel" :: summary__global_quantities__energy_fast_parallel = summary__global_quantities__energy_fast_parallel()
    var"power_radiated" :: summary__global_quantities__power_radiated = summary__global_quantities__power_radiated()
    var"denergy_thermal_dt" :: summary__global_quantities__denergy_thermal_dt = summary__global_quantities__denergy_thermal_dt()
    var"energy_mhd" :: summary__global_quantities__energy_mhd = summary__global_quantities__energy_mhd()
    var"li" :: summary__global_quantities__li = summary__global_quantities__li()
    var"ratio_tau_helium_fuel" :: summary__global_quantities__ratio_tau_helium_fuel = summary__global_quantities__ratio_tau_helium_fuel()
    var"r0" :: summary__global_quantities__r0 = summary__global_quantities__r0()
    var"beta_tor_norm_mhd" :: summary__global_quantities__beta_tor_norm_mhd = summary__global_quantities__beta_tor_norm_mhd()
    _parent :: WeakRef = WeakRef(missing)
    function summary__global_quantities(var"power_radiated_outside_lcfs"=summary__global_quantities__power_radiated_outside_lcfs(), var"beta_tor"=summary__global_quantities__beta_tor(), var"power_ohm"=summary__global_quantities__power_ohm(), var"energy_b_field_pol"=summary__global_quantities__energy_b_field_pol(), var"tau_energy"=summary__global_quantities__tau_energy(), var"energy_electrons_thermal"=summary__global_quantities__energy_electrons_thermal(), var"greenwald_fraction"=summary__global_quantities__greenwald_fraction(), var"power_steady"=summary__global_quantities__power_steady(), var"resistance"=summary__global_quantities__resistance(), var"beta_pol"=summary__global_quantities__beta_pol(), var"energy_thermal"=summary__global_quantities__energy_thermal(), var"fusion_gain"=summary__global_quantities__fusion_gain(), var"power_radiated_inside_lcfs"=summary__global_quantities__power_radiated_inside_lcfs(), var"beta_tor_norm"=summary__global_quantities__beta_tor_norm(), var"beta_tor_thermal_norm"=summary__global_quantities__beta_tor_thermal_norm(), var"tau_helium"=summary__global_quantities__tau_helium(), var"power_loss"=summary__global_quantities__power_loss(), var"current_alignment"=summary__global_quantities__current_alignment(), var"energy_diamagnetic"=summary__global_quantities__energy_diamagnetic(), var"tau_energy_98"=summary__global_quantities__tau_energy_98(), var"li_mhd"=summary__global_quantities__li_mhd(), var"energy_fast_perpendicular"=summary__global_quantities__energy_fast_perpendicular(), var"q_95"=summary__global_quantities__q_95(), var"beta_tor_mhd"=summary__global_quantities__beta_tor_mhd(), var"volume"=summary__global_quantities__volume(), var"b0"=summary__global_quantities__b0(), var"ip"=summary__global_quantities__ip(), var"beta_pol_mhd"=summary__global_quantities__beta_pol_mhd(), var"denergy_diamagnetic_dt"=summary__global_quantities__denergy_diamagnetic_dt(), var"h_mode"=summary__global_quantities__h_mode(), var"power_bremsstrahlung"=summary__global_quantities__power_bremsstrahlung(), var"h_98"=summary__global_quantities__h_98(), var"tau_resistive"=summary__global_quantities__tau_resistive(), var"current_non_inductive"=summary__global_quantities__current_non_inductive(), var"current_ohm"=summary__global_quantities__current_ohm(), var"energy_total"=summary__global_quantities__energy_total(), var"power_synchrotron"=summary__global_quantities__power_synchrotron(), var"v_loop"=summary__global_quantities__v_loop(), var"energy_ion_total_thermal"=summary__global_quantities__energy_ion_total_thermal(), var"power_line"=summary__global_quantities__power_line(), var"current_bootstrap"=summary__global_quantities__current_bootstrap(), var"fusion_fluence"=summary__global_quantities__fusion_fluence(), var"energy_fast_parallel"=summary__global_quantities__energy_fast_parallel(), var"power_radiated"=summary__global_quantities__power_radiated(), var"denergy_thermal_dt"=summary__global_quantities__denergy_thermal_dt(), var"energy_mhd"=summary__global_quantities__energy_mhd(), var"li"=summary__global_quantities__li(), var"ratio_tau_helium_fuel"=summary__global_quantities__ratio_tau_helium_fuel(), var"r0"=summary__global_quantities__r0(), var"beta_tor_norm_mhd"=summary__global_quantities__beta_tor_norm_mhd(), _parent=WeakRef(missing))
        fds = new(var"power_radiated_outside_lcfs", var"beta_tor", var"power_ohm", var"energy_b_field_pol", var"tau_energy", var"energy_electrons_thermal", var"greenwald_fraction", var"power_steady", var"resistance", var"beta_pol", var"energy_thermal", var"fusion_gain", var"power_radiated_inside_lcfs", var"beta_tor_norm", var"beta_tor_thermal_norm", var"tau_helium", var"power_loss", var"current_alignment", var"energy_diamagnetic", var"tau_energy_98", var"li_mhd", var"energy_fast_perpendicular", var"q_95", var"beta_tor_mhd", var"volume", var"b0", var"ip", var"beta_pol_mhd", var"denergy_diamagnetic_dt", var"h_mode", var"power_bremsstrahlung", var"h_98", var"tau_resistive", var"current_non_inductive", var"current_ohm", var"energy_total", var"power_synchrotron", var"v_loop", var"energy_ion_total_thermal", var"power_line", var"current_bootstrap", var"fusion_fluence", var"energy_fast_parallel", var"power_radiated", var"denergy_thermal_dt", var"energy_mhd", var"li", var"ratio_tau_helium_fuel", var"r0", var"beta_tor_norm_mhd", _parent)
        assign_expressions(fds)
        setfield!(fds.power_radiated_outside_lcfs, :_parent, WeakRef(fds))
        setfield!(fds.beta_tor, :_parent, WeakRef(fds))
        setfield!(fds.power_ohm, :_parent, WeakRef(fds))
        setfield!(fds.energy_b_field_pol, :_parent, WeakRef(fds))
        setfield!(fds.tau_energy, :_parent, WeakRef(fds))
        setfield!(fds.energy_electrons_thermal, :_parent, WeakRef(fds))
        setfield!(fds.greenwald_fraction, :_parent, WeakRef(fds))
        setfield!(fds.power_steady, :_parent, WeakRef(fds))
        setfield!(fds.resistance, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol, :_parent, WeakRef(fds))
        setfield!(fds.energy_thermal, :_parent, WeakRef(fds))
        setfield!(fds.fusion_gain, :_parent, WeakRef(fds))
        setfield!(fds.power_radiated_inside_lcfs, :_parent, WeakRef(fds))
        setfield!(fds.beta_tor_norm, :_parent, WeakRef(fds))
        setfield!(fds.beta_tor_thermal_norm, :_parent, WeakRef(fds))
        setfield!(fds.tau_helium, :_parent, WeakRef(fds))
        setfield!(fds.power_loss, :_parent, WeakRef(fds))
        setfield!(fds.current_alignment, :_parent, WeakRef(fds))
        setfield!(fds.energy_diamagnetic, :_parent, WeakRef(fds))
        setfield!(fds.tau_energy_98, :_parent, WeakRef(fds))
        setfield!(fds.li_mhd, :_parent, WeakRef(fds))
        setfield!(fds.energy_fast_perpendicular, :_parent, WeakRef(fds))
        setfield!(fds.q_95, :_parent, WeakRef(fds))
        setfield!(fds.beta_tor_mhd, :_parent, WeakRef(fds))
        setfield!(fds.volume, :_parent, WeakRef(fds))
        setfield!(fds.b0, :_parent, WeakRef(fds))
        setfield!(fds.ip, :_parent, WeakRef(fds))
        setfield!(fds.beta_pol_mhd, :_parent, WeakRef(fds))
        setfield!(fds.denergy_diamagnetic_dt, :_parent, WeakRef(fds))
        setfield!(fds.h_mode, :_parent, WeakRef(fds))
        setfield!(fds.power_bremsstrahlung, :_parent, WeakRef(fds))
        setfield!(fds.h_98, :_parent, WeakRef(fds))
        setfield!(fds.tau_resistive, :_parent, WeakRef(fds))
        setfield!(fds.current_non_inductive, :_parent, WeakRef(fds))
        setfield!(fds.current_ohm, :_parent, WeakRef(fds))
        setfield!(fds.energy_total, :_parent, WeakRef(fds))
        setfield!(fds.power_synchrotron, :_parent, WeakRef(fds))
        setfield!(fds.v_loop, :_parent, WeakRef(fds))
        setfield!(fds.energy_ion_total_thermal, :_parent, WeakRef(fds))
        setfield!(fds.power_line, :_parent, WeakRef(fds))
        setfield!(fds.current_bootstrap, :_parent, WeakRef(fds))
        setfield!(fds.fusion_fluence, :_parent, WeakRef(fds))
        setfield!(fds.energy_fast_parallel, :_parent, WeakRef(fds))
        setfield!(fds.power_radiated, :_parent, WeakRef(fds))
        setfield!(fds.denergy_thermal_dt, :_parent, WeakRef(fds))
        setfield!(fds.energy_mhd, :_parent, WeakRef(fds))
        setfield!(fds.li, :_parent, WeakRef(fds))
        setfield!(fds.ratio_tau_helium_fuel, :_parent, WeakRef(fds))
        setfield!(fds.r0, :_parent, WeakRef(fds))
        setfield!(fds.beta_tor_norm_mhd, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__xenon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__tritium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__top <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__top(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__silane <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__silane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__propane <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__propane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__oxygen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__nitrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__neon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__midplane <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__midplane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__methane_deuterated <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__methane_deuterated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__methane_carbon_13 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__methane_carbon_13(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__methane <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__methane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__lithium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__krypton <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__impurity_seeding <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__impurity_seeding(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__hydrogen <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__helium_4 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__helium_3 <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__ethylene <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__ethylene(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__ethane <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__ethane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__deuterium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__carbon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__bottom <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__bottom(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__beryllium <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__argon <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__ammonia_deuterated <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__ammonia_deuterated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates__ammonia <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates__ammonia(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__gas_injection_rates <: FDS
    var"krypton" :: summary__gas_injection_rates__krypton = summary__gas_injection_rates__krypton()
    var"lithium" :: summary__gas_injection_rates__lithium = summary__gas_injection_rates__lithium()
    var"methane_carbon_13" :: summary__gas_injection_rates__methane_carbon_13 = summary__gas_injection_rates__methane_carbon_13()
    var"neon" :: summary__gas_injection_rates__neon = summary__gas_injection_rates__neon()
    var"ethylene" :: summary__gas_injection_rates__ethylene = summary__gas_injection_rates__ethylene()
    var"tritium" :: summary__gas_injection_rates__tritium = summary__gas_injection_rates__tritium()
    var"helium_3" :: summary__gas_injection_rates__helium_3 = summary__gas_injection_rates__helium_3()
    var"deuterium" :: summary__gas_injection_rates__deuterium = summary__gas_injection_rates__deuterium()
    var"bottom" :: summary__gas_injection_rates__bottom = summary__gas_injection_rates__bottom()
    var"helium_4" :: summary__gas_injection_rates__helium_4 = summary__gas_injection_rates__helium_4()
    var"oxygen" :: summary__gas_injection_rates__oxygen = summary__gas_injection_rates__oxygen()
    var"xenon" :: summary__gas_injection_rates__xenon = summary__gas_injection_rates__xenon()
    var"methane" :: summary__gas_injection_rates__methane = summary__gas_injection_rates__methane()
    var"hydrogen" :: summary__gas_injection_rates__hydrogen = summary__gas_injection_rates__hydrogen()
    var"carbon" :: summary__gas_injection_rates__carbon = summary__gas_injection_rates__carbon()
    var"ammonia_deuterated" :: summary__gas_injection_rates__ammonia_deuterated = summary__gas_injection_rates__ammonia_deuterated()
    var"total" :: summary__gas_injection_rates__total = summary__gas_injection_rates__total()
    var"top" :: summary__gas_injection_rates__top = summary__gas_injection_rates__top()
    var"midplane" :: summary__gas_injection_rates__midplane = summary__gas_injection_rates__midplane()
    var"ethane" :: summary__gas_injection_rates__ethane = summary__gas_injection_rates__ethane()
    var"silane" :: summary__gas_injection_rates__silane = summary__gas_injection_rates__silane()
    var"ammonia" :: summary__gas_injection_rates__ammonia = summary__gas_injection_rates__ammonia()
    var"methane_deuterated" :: summary__gas_injection_rates__methane_deuterated = summary__gas_injection_rates__methane_deuterated()
    var"nitrogen" :: summary__gas_injection_rates__nitrogen = summary__gas_injection_rates__nitrogen()
    var"propane" :: summary__gas_injection_rates__propane = summary__gas_injection_rates__propane()
    var"beryllium" :: summary__gas_injection_rates__beryllium = summary__gas_injection_rates__beryllium()
    var"argon" :: summary__gas_injection_rates__argon = summary__gas_injection_rates__argon()
    var"impurity_seeding" :: summary__gas_injection_rates__impurity_seeding = summary__gas_injection_rates__impurity_seeding()
    _parent :: WeakRef = WeakRef(missing)
    function summary__gas_injection_rates(var"krypton"=summary__gas_injection_rates__krypton(), var"lithium"=summary__gas_injection_rates__lithium(), var"methane_carbon_13"=summary__gas_injection_rates__methane_carbon_13(), var"neon"=summary__gas_injection_rates__neon(), var"ethylene"=summary__gas_injection_rates__ethylene(), var"tritium"=summary__gas_injection_rates__tritium(), var"helium_3"=summary__gas_injection_rates__helium_3(), var"deuterium"=summary__gas_injection_rates__deuterium(), var"bottom"=summary__gas_injection_rates__bottom(), var"helium_4"=summary__gas_injection_rates__helium_4(), var"oxygen"=summary__gas_injection_rates__oxygen(), var"xenon"=summary__gas_injection_rates__xenon(), var"methane"=summary__gas_injection_rates__methane(), var"hydrogen"=summary__gas_injection_rates__hydrogen(), var"carbon"=summary__gas_injection_rates__carbon(), var"ammonia_deuterated"=summary__gas_injection_rates__ammonia_deuterated(), var"total"=summary__gas_injection_rates__total(), var"top"=summary__gas_injection_rates__top(), var"midplane"=summary__gas_injection_rates__midplane(), var"ethane"=summary__gas_injection_rates__ethane(), var"silane"=summary__gas_injection_rates__silane(), var"ammonia"=summary__gas_injection_rates__ammonia(), var"methane_deuterated"=summary__gas_injection_rates__methane_deuterated(), var"nitrogen"=summary__gas_injection_rates__nitrogen(), var"propane"=summary__gas_injection_rates__propane(), var"beryllium"=summary__gas_injection_rates__beryllium(), var"argon"=summary__gas_injection_rates__argon(), var"impurity_seeding"=summary__gas_injection_rates__impurity_seeding(), _parent=WeakRef(missing))
        fds = new(var"krypton", var"lithium", var"methane_carbon_13", var"neon", var"ethylene", var"tritium", var"helium_3", var"deuterium", var"bottom", var"helium_4", var"oxygen", var"xenon", var"methane", var"hydrogen", var"carbon", var"ammonia_deuterated", var"total", var"top", var"midplane", var"ethane", var"silane", var"ammonia", var"methane_deuterated", var"nitrogen", var"propane", var"beryllium", var"argon", var"impurity_seeding", _parent)
        assign_expressions(fds)
        setfield!(fds.krypton, :_parent, WeakRef(fds))
        setfield!(fds.lithium, :_parent, WeakRef(fds))
        setfield!(fds.methane_carbon_13, :_parent, WeakRef(fds))
        setfield!(fds.neon, :_parent, WeakRef(fds))
        setfield!(fds.ethylene, :_parent, WeakRef(fds))
        setfield!(fds.tritium, :_parent, WeakRef(fds))
        setfield!(fds.helium_3, :_parent, WeakRef(fds))
        setfield!(fds.deuterium, :_parent, WeakRef(fds))
        setfield!(fds.bottom, :_parent, WeakRef(fds))
        setfield!(fds.helium_4, :_parent, WeakRef(fds))
        setfield!(fds.oxygen, :_parent, WeakRef(fds))
        setfield!(fds.xenon, :_parent, WeakRef(fds))
        setfield!(fds.methane, :_parent, WeakRef(fds))
        setfield!(fds.hydrogen, :_parent, WeakRef(fds))
        setfield!(fds.carbon, :_parent, WeakRef(fds))
        setfield!(fds.ammonia_deuterated, :_parent, WeakRef(fds))
        setfield!(fds.total, :_parent, WeakRef(fds))
        setfield!(fds.top, :_parent, WeakRef(fds))
        setfield!(fds.midplane, :_parent, WeakRef(fds))
        setfield!(fds.ethane, :_parent, WeakRef(fds))
        setfield!(fds.silane, :_parent, WeakRef(fds))
        setfield!(fds.ammonia, :_parent, WeakRef(fds))
        setfield!(fds.methane_deuterated, :_parent, WeakRef(fds))
        setfield!(fds.nitrogen, :_parent, WeakRef(fds))
        setfield!(fds.propane, :_parent, WeakRef(fds))
        setfield!(fds.beryllium, :_parent, WeakRef(fds))
        setfield!(fds.argon, :_parent, WeakRef(fds))
        setfield!(fds.impurity_seeding, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__power <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_power_total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_power_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__tt__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__tt__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__beam_thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__tt__beam_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt__beam_beam <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__tt__beam_beam(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__tt <: FDS
    var"beam_thermal" :: summary__fusion__neutron_fluxes__tt__beam_thermal = summary__fusion__neutron_fluxes__tt__beam_thermal()
    var"beam_beam" :: summary__fusion__neutron_fluxes__tt__beam_beam = summary__fusion__neutron_fluxes__tt__beam_beam()
    var"total" :: summary__fusion__neutron_fluxes__tt__total = summary__fusion__neutron_fluxes__tt__total()
    var"thermal" :: summary__fusion__neutron_fluxes__tt__thermal = summary__fusion__neutron_fluxes__tt__thermal()
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__tt(var"beam_thermal"=summary__fusion__neutron_fluxes__tt__beam_thermal(), var"beam_beam"=summary__fusion__neutron_fluxes__tt__beam_beam(), var"total"=summary__fusion__neutron_fluxes__tt__total(), var"thermal"=summary__fusion__neutron_fluxes__tt__thermal(), _parent=WeakRef(missing))
        fds = new(var"beam_thermal", var"beam_beam", var"total", var"thermal", _parent)
        assign_expressions(fds)
        setfield!(fds.beam_thermal, :_parent, WeakRef(fds))
        setfield!(fds.beam_beam, :_parent, WeakRef(fds))
        setfield!(fds.total, :_parent, WeakRef(fds))
        setfield!(fds.thermal, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dt__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dt__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__beam_thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dt__beam_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt__beam_beam <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dt__beam_beam(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dt <: FDS
    var"beam_thermal" :: summary__fusion__neutron_fluxes__dt__beam_thermal = summary__fusion__neutron_fluxes__dt__beam_thermal()
    var"beam_beam" :: summary__fusion__neutron_fluxes__dt__beam_beam = summary__fusion__neutron_fluxes__dt__beam_beam()
    var"total" :: summary__fusion__neutron_fluxes__dt__total = summary__fusion__neutron_fluxes__dt__total()
    var"thermal" :: summary__fusion__neutron_fluxes__dt__thermal = summary__fusion__neutron_fluxes__dt__thermal()
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dt(var"beam_thermal"=summary__fusion__neutron_fluxes__dt__beam_thermal(), var"beam_beam"=summary__fusion__neutron_fluxes__dt__beam_beam(), var"total"=summary__fusion__neutron_fluxes__dt__total(), var"thermal"=summary__fusion__neutron_fluxes__dt__thermal(), _parent=WeakRef(missing))
        fds = new(var"beam_thermal", var"beam_beam", var"total", var"thermal", _parent)
        assign_expressions(fds)
        setfield!(fds.beam_thermal, :_parent, WeakRef(fds))
        setfield!(fds.beam_beam, :_parent, WeakRef(fds))
        setfield!(fds.total, :_parent, WeakRef(fds))
        setfield!(fds.thermal, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__total <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dd__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dd__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__beam_thermal <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dd__beam_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd__beam_beam <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dd__beam_beam(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes__dd <: FDS
    var"beam_thermal" :: summary__fusion__neutron_fluxes__dd__beam_thermal = summary__fusion__neutron_fluxes__dd__beam_thermal()
    var"beam_beam" :: summary__fusion__neutron_fluxes__dd__beam_beam = summary__fusion__neutron_fluxes__dd__beam_beam()
    var"total" :: summary__fusion__neutron_fluxes__dd__total = summary__fusion__neutron_fluxes__dd__total()
    var"thermal" :: summary__fusion__neutron_fluxes__dd__thermal = summary__fusion__neutron_fluxes__dd__thermal()
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes__dd(var"beam_thermal"=summary__fusion__neutron_fluxes__dd__beam_thermal(), var"beam_beam"=summary__fusion__neutron_fluxes__dd__beam_beam(), var"total"=summary__fusion__neutron_fluxes__dd__total(), var"thermal"=summary__fusion__neutron_fluxes__dd__thermal(), _parent=WeakRef(missing))
        fds = new(var"beam_thermal", var"beam_beam", var"total", var"thermal", _parent)
        assign_expressions(fds)
        setfield!(fds.beam_thermal, :_parent, WeakRef(fds))
        setfield!(fds.beam_beam, :_parent, WeakRef(fds))
        setfield!(fds.total, :_parent, WeakRef(fds))
        setfield!(fds.thermal, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__neutron_fluxes <: FDS
    var"tt" :: summary__fusion__neutron_fluxes__tt = summary__fusion__neutron_fluxes__tt()
    var"total" :: summary__fusion__neutron_fluxes__total = summary__fusion__neutron_fluxes__total()
    var"thermal" :: summary__fusion__neutron_fluxes__thermal = summary__fusion__neutron_fluxes__thermal()
    var"dt" :: summary__fusion__neutron_fluxes__dt = summary__fusion__neutron_fluxes__dt()
    var"dd" :: summary__fusion__neutron_fluxes__dd = summary__fusion__neutron_fluxes__dd()
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__neutron_fluxes(var"tt"=summary__fusion__neutron_fluxes__tt(), var"total"=summary__fusion__neutron_fluxes__total(), var"thermal"=summary__fusion__neutron_fluxes__thermal(), var"dt"=summary__fusion__neutron_fluxes__dt(), var"dd"=summary__fusion__neutron_fluxes__dd(), _parent=WeakRef(missing))
        fds = new(var"tt", var"total", var"thermal", var"dt", var"dd", _parent)
        assign_expressions(fds)
        setfield!(fds.tt, :_parent, WeakRef(fds))
        setfield!(fds.total, :_parent, WeakRef(fds))
        setfield!(fds.thermal, :_parent, WeakRef(fds))
        setfield!(fds.dt, :_parent, WeakRef(fds))
        setfield!(fds.dd, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion__current <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion__current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__fusion <: FDS
    var"neutron_power_total" :: summary__fusion__neutron_power_total = summary__fusion__neutron_power_total()
    var"neutron_fluxes" :: summary__fusion__neutron_fluxes = summary__fusion__neutron_fluxes()
    var"power" :: summary__fusion__power = summary__fusion__power()
    var"current" :: summary__fusion__current = summary__fusion__current()
    _parent :: WeakRef = WeakRef(missing)
    function summary__fusion(var"neutron_power_total"=summary__fusion__neutron_power_total(), var"neutron_fluxes"=summary__fusion__neutron_fluxes(), var"power"=summary__fusion__power(), var"current"=summary__fusion__current(), _parent=WeakRef(missing))
        fds = new(var"neutron_power_total", var"neutron_fluxes", var"power", var"current", _parent)
        assign_expressions(fds)
        setfield!(fds.neutron_power_total, :_parent, WeakRef(fds))
        setfield!(fds.neutron_fluxes, :_parent, WeakRef(fds))
        setfield!(fds.power, :_parent, WeakRef(fds))
        setfield!(fds.current, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__elms__type <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__elms__type(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__elms__frequency <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__elms__frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__elms <: FDS
    var"frequency" :: summary__elms__frequency = summary__elms__frequency()
    var"type" :: summary__elms__type = summary__elms__type()
    _parent :: WeakRef = WeakRef(missing)
    function summary__elms(var"frequency"=summary__elms__frequency(), var"type"=summary__elms__type(), _parent=WeakRef(missing))
        fds = new(var"frequency", var"type", _parent)
        assign_expressions(fds)
        setfield!(fds.frequency, :_parent, WeakRef(fds))
        setfield!(fds.type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__disruption__vertical_displacement <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__disruption__vertical_displacement(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__disruption__time_radiated_power_max <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__disruption__time_radiated_power_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__disruption__time_half_ip <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__disruption__time_half_ip(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__disruption__time <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__disruption__time(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__disruption__mitigation_valve <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__disruption__mitigation_valve(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__disruption <: FDS
    var"time" :: summary__disruption__time = summary__disruption__time()
    var"mitigation_valve" :: summary__disruption__mitigation_valve = summary__disruption__mitigation_valve()
    var"vertical_displacement" :: summary__disruption__vertical_displacement = summary__disruption__vertical_displacement()
    var"time_half_ip" :: summary__disruption__time_half_ip = summary__disruption__time_half_ip()
    var"time_radiated_power_max" :: summary__disruption__time_radiated_power_max = summary__disruption__time_radiated_power_max()
    _parent :: WeakRef = WeakRef(missing)
    function summary__disruption(var"time"=summary__disruption__time(), var"mitigation_valve"=summary__disruption__mitigation_valve(), var"vertical_displacement"=summary__disruption__vertical_displacement(), var"time_half_ip"=summary__disruption__time_half_ip(), var"time_radiated_power_max"=summary__disruption__time_radiated_power_max(), _parent=WeakRef(missing))
        fds = new(var"time", var"mitigation_valve", var"vertical_displacement", var"time_half_ip", var"time_radiated_power_max", _parent)
        assign_expressions(fds)
        setfield!(fds.time, :_parent, WeakRef(fds))
        setfield!(fds.mitigation_valve, :_parent, WeakRef(fds))
        setfield!(fds.vertical_displacement, :_parent, WeakRef(fds))
        setfield!(fds.time_half_ip, :_parent, WeakRef(fds))
        setfield!(fds.time_radiated_power_max, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__configuration <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__configuration(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__code <: FDS
    var"library" :: FDSvector{T} where {T<:summary__code__library} = FDSvector(summary__code__library[])
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"output_flag" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__code(var"library"=FDSvector(summary__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__type <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__type(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__triangularity_upper <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__triangularity_upper(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__triangularity_lower <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__triangularity_lower(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__strike_point_outer_z <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__strike_point_outer_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__strike_point_outer_r <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__strike_point_outer_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__strike_point_inner_z <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__strike_point_inner_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__strike_point_inner_r <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__strike_point_inner_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__strike_point_configuration <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__strike_point_configuration(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__minor_radius <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__minor_radius(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__magnetic_axis_z <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__magnetic_axis_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__magnetic_axis_r <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__magnetic_axis_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__geometric_axis_z <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__geometric_axis_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__geometric_axis_r <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__geometric_axis_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__gap_limiter_wall <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__gap_limiter_wall(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary__elongation <: FDS
    var"source" :: Union{Missing, String, Function} = missing
    var"value" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary__elongation(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        fds = new(var"source", var"value", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct summary__boundary <: FDS
    var"strike_point_inner_r" :: summary__boundary__strike_point_inner_r = summary__boundary__strike_point_inner_r()
    var"strike_point_outer_z" :: summary__boundary__strike_point_outer_z = summary__boundary__strike_point_outer_z()
    var"geometric_axis_z" :: summary__boundary__geometric_axis_z = summary__boundary__geometric_axis_z()
    var"magnetic_axis_z" :: summary__boundary__magnetic_axis_z = summary__boundary__magnetic_axis_z()
    var"geometric_axis_r" :: summary__boundary__geometric_axis_r = summary__boundary__geometric_axis_r()
    var"strike_point_configuration" :: summary__boundary__strike_point_configuration = summary__boundary__strike_point_configuration()
    var"triangularity_upper" :: summary__boundary__triangularity_upper = summary__boundary__triangularity_upper()
    var"gap_limiter_wall" :: summary__boundary__gap_limiter_wall = summary__boundary__gap_limiter_wall()
    var"strike_point_inner_z" :: summary__boundary__strike_point_inner_z = summary__boundary__strike_point_inner_z()
    var"triangularity_lower" :: summary__boundary__triangularity_lower = summary__boundary__triangularity_lower()
    var"minor_radius" :: summary__boundary__minor_radius = summary__boundary__minor_radius()
    var"strike_point_outer_r" :: summary__boundary__strike_point_outer_r = summary__boundary__strike_point_outer_r()
    var"elongation" :: summary__boundary__elongation = summary__boundary__elongation()
    var"type" :: summary__boundary__type = summary__boundary__type()
    var"magnetic_axis_r" :: summary__boundary__magnetic_axis_r = summary__boundary__magnetic_axis_r()
    _parent :: WeakRef = WeakRef(missing)
    function summary__boundary(var"strike_point_inner_r"=summary__boundary__strike_point_inner_r(), var"strike_point_outer_z"=summary__boundary__strike_point_outer_z(), var"geometric_axis_z"=summary__boundary__geometric_axis_z(), var"magnetic_axis_z"=summary__boundary__magnetic_axis_z(), var"geometric_axis_r"=summary__boundary__geometric_axis_r(), var"strike_point_configuration"=summary__boundary__strike_point_configuration(), var"triangularity_upper"=summary__boundary__triangularity_upper(), var"gap_limiter_wall"=summary__boundary__gap_limiter_wall(), var"strike_point_inner_z"=summary__boundary__strike_point_inner_z(), var"triangularity_lower"=summary__boundary__triangularity_lower(), var"minor_radius"=summary__boundary__minor_radius(), var"strike_point_outer_r"=summary__boundary__strike_point_outer_r(), var"elongation"=summary__boundary__elongation(), var"type"=summary__boundary__type(), var"magnetic_axis_r"=summary__boundary__magnetic_axis_r(), _parent=WeakRef(missing))
        fds = new(var"strike_point_inner_r", var"strike_point_outer_z", var"geometric_axis_z", var"magnetic_axis_z", var"geometric_axis_r", var"strike_point_configuration", var"triangularity_upper", var"gap_limiter_wall", var"strike_point_inner_z", var"triangularity_lower", var"minor_radius", var"strike_point_outer_r", var"elongation", var"type", var"magnetic_axis_r", _parent)
        assign_expressions(fds)
        setfield!(fds.strike_point_inner_r, :_parent, WeakRef(fds))
        setfield!(fds.strike_point_outer_z, :_parent, WeakRef(fds))
        setfield!(fds.geometric_axis_z, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_axis_z, :_parent, WeakRef(fds))
        setfield!(fds.geometric_axis_r, :_parent, WeakRef(fds))
        setfield!(fds.strike_point_configuration, :_parent, WeakRef(fds))
        setfield!(fds.triangularity_upper, :_parent, WeakRef(fds))
        setfield!(fds.gap_limiter_wall, :_parent, WeakRef(fds))
        setfield!(fds.strike_point_inner_z, :_parent, WeakRef(fds))
        setfield!(fds.triangularity_lower, :_parent, WeakRef(fds))
        setfield!(fds.minor_radius, :_parent, WeakRef(fds))
        setfield!(fds.strike_point_outer_r, :_parent, WeakRef(fds))
        setfield!(fds.elongation, :_parent, WeakRef(fds))
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_axis_r, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct summary <: FDS
    var"local" :: summary__local = summary__local()
    var"wall" :: summary__wall = summary__wall()
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"code" :: summary__code = summary__code()
    var"stationary_phase_flag" :: summary__stationary_phase_flag = summary__stationary_phase_flag()
    var"gas_injection_rates" :: summary__gas_injection_rates = summary__gas_injection_rates()
    var"configuration" :: summary__configuration = summary__configuration()
    var"boundary" :: summary__boundary = summary__boundary()
    var"pedestal_fits" :: summary__pedestal_fits = summary__pedestal_fits()
    var"heating_current_drive" :: summary__heating_current_drive = summary__heating_current_drive()
    var"global_quantities" :: summary__global_quantities = summary__global_quantities()
    var"disruption" :: summary__disruption = summary__disruption()
    var"fusion" :: summary__fusion = summary__fusion()
    var"rmps" :: summary__rmps = summary__rmps()
    var"ids_properties" :: summary__ids_properties = summary__ids_properties()
    var"midplane" :: summary__midplane = summary__midplane()
    var"pellets" :: summary__pellets = summary__pellets()
    var"elms" :: summary__elms = summary__elms()
    var"kicks" :: summary__kicks = summary__kicks()
    var"tag" :: summary__tag = summary__tag()
    var"line_average" :: summary__line_average = summary__line_average()
    var"scrape_off_layer" :: summary__scrape_off_layer = summary__scrape_off_layer()
    var"limiter" :: summary__limiter = summary__limiter()
    var"runaways" :: summary__runaways = summary__runaways()
    var"time_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"magnetic_shear_flag" :: summary__magnetic_shear_flag = summary__magnetic_shear_flag()
    var"volume_average" :: summary__volume_average = summary__volume_average()
    _parent :: WeakRef = WeakRef(missing)
    function summary(var"local"=summary__local(), var"wall"=summary__wall(), var"time"=missing, var"code"=summary__code(), var"stationary_phase_flag"=summary__stationary_phase_flag(), var"gas_injection_rates"=summary__gas_injection_rates(), var"configuration"=summary__configuration(), var"boundary"=summary__boundary(), var"pedestal_fits"=summary__pedestal_fits(), var"heating_current_drive"=summary__heating_current_drive(), var"global_quantities"=summary__global_quantities(), var"disruption"=summary__disruption(), var"fusion"=summary__fusion(), var"rmps"=summary__rmps(), var"ids_properties"=summary__ids_properties(), var"midplane"=summary__midplane(), var"pellets"=summary__pellets(), var"elms"=summary__elms(), var"kicks"=summary__kicks(), var"tag"=summary__tag(), var"line_average"=summary__line_average(), var"scrape_off_layer"=summary__scrape_off_layer(), var"limiter"=summary__limiter(), var"runaways"=summary__runaways(), var"time_width"=missing, var"magnetic_shear_flag"=summary__magnetic_shear_flag(), var"volume_average"=summary__volume_average(), _parent=WeakRef(missing))
        fds = new(var"local", var"wall", var"time", var"code", var"stationary_phase_flag", var"gas_injection_rates", var"configuration", var"boundary", var"pedestal_fits", var"heating_current_drive", var"global_quantities", var"disruption", var"fusion", var"rmps", var"ids_properties", var"midplane", var"pellets", var"elms", var"kicks", var"tag", var"line_average", var"scrape_off_layer", var"limiter", var"runaways", var"time_width", var"magnetic_shear_flag", var"volume_average", _parent)
        assign_expressions(fds)
        setfield!(fds.local, :_parent, WeakRef(fds))
        setfield!(fds.wall, :_parent, WeakRef(fds))
        setfield!(fds.code, :_parent, WeakRef(fds))
        setfield!(fds.stationary_phase_flag, :_parent, WeakRef(fds))
        setfield!(fds.gas_injection_rates, :_parent, WeakRef(fds))
        setfield!(fds.configuration, :_parent, WeakRef(fds))
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        setfield!(fds.pedestal_fits, :_parent, WeakRef(fds))
        setfield!(fds.heating_current_drive, :_parent, WeakRef(fds))
        setfield!(fds.global_quantities, :_parent, WeakRef(fds))
        setfield!(fds.disruption, :_parent, WeakRef(fds))
        setfield!(fds.fusion, :_parent, WeakRef(fds))
        setfield!(fds.rmps, :_parent, WeakRef(fds))
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.midplane, :_parent, WeakRef(fds))
        setfield!(fds.pellets, :_parent, WeakRef(fds))
        setfield!(fds.elms, :_parent, WeakRef(fds))
        setfield!(fds.kicks, :_parent, WeakRef(fds))
        setfield!(fds.tag, :_parent, WeakRef(fds))
        setfield!(fds.line_average, :_parent, WeakRef(fds))
        setfield!(fds.scrape_off_layer, :_parent, WeakRef(fds))
        setfield!(fds.limiter, :_parent, WeakRef(fds))
        setfield!(fds.runaways, :_parent, WeakRef(fds))
        setfield!(fds.magnetic_shear_flag, :_parent, WeakRef(fds))
        setfield!(fds.volume_average, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"r0" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        fds = new(var"b0", var"r0", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid_type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_2d___grid_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid <: FDS
    var"volume_element" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"dim2" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dim1" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_2d___grid(var"volume_element"=missing, var"dim2"=missing, var"dim1"=missing, _parent=WeakRef(missing))
        fds = new(var"volume_element", var"dim2", var"dim1", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d <: FDSvectorElement
    var"psi" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"b_field_r" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"r" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"b_r" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"theta" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"b_field_z" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"j_tor" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"phi" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"b_field_tor" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"b_z" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"grid" :: equilibrium__time_slice___profiles_2d___grid = equilibrium__time_slice___profiles_2d___grid()
    var"grid_type" :: equilibrium__time_slice___profiles_2d___grid_type = equilibrium__time_slice___profiles_2d___grid_type()
    var"j_parallel" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"b_tor" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_2d(var"psi"=missing, var"b_field_r"=missing, var"r"=missing, var"b_r"=missing, var"theta"=missing, var"b_field_z"=missing, var"j_tor"=missing, var"phi"=missing, var"z"=missing, var"b_field_tor"=missing, var"b_z"=missing, var"grid"=equilibrium__time_slice___profiles_2d___grid(), var"grid_type"=equilibrium__time_slice___profiles_2d___grid_type(), var"j_parallel"=missing, var"b_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"b_field_r", var"r", var"b_r", var"theta", var"b_field_z", var"j_tor", var"phi", var"z", var"b_field_tor", var"b_z", var"grid", var"grid_type", var"j_parallel", var"b_tor", _parent)
        assign_expressions(fds)
        setfield!(fds.grid, :_parent, WeakRef(fds))
        setfield!(fds.grid_type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d__geometric_axis <: FDS
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_1d__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d <: FDS
    var"b_field_max" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dvolume_drho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm9" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dpsi_drho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"surface" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"magnetic_shear" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"b_average" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"b_field_min" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"darea_dpsi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm3" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"squareness_upper_inner" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"squareness_lower_inner" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"elongation" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"beta_pol" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"b_field_average" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"j_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm6" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"psi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm8" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dpressure_dpsi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"triangularity_upper" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"darea_drho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"area" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"trapped_fraction" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"volume" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dvolume_dpsi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"b_min" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"f" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"mass_density" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"r_outboard" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm4" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"phi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"squareness_lower_outer" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"triangularity_lower" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm2" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_volume_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm1" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm5" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"b_max" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"f_df_dpsi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"j_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"r_inboard" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"q" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"gm7" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"squareness_upper_outer" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"geometric_axis" :: equilibrium__time_slice___profiles_1d__geometric_axis = equilibrium__time_slice___profiles_1d__geometric_axis()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_1d(var"b_field_max"=missing, var"dvolume_drho_tor"=missing, var"gm9"=missing, var"dpsi_drho_tor"=missing, var"surface"=missing, var"rho_tor"=missing, var"magnetic_shear"=missing, var"b_average"=missing, var"b_field_min"=missing, var"darea_dpsi"=missing, var"gm3"=missing, var"squareness_upper_inner"=missing, var"squareness_lower_inner"=missing, var"rho_tor_norm"=missing, var"elongation"=missing, var"beta_pol"=missing, var"b_field_average"=missing, var"j_parallel"=missing, var"gm6"=missing, var"psi"=missing, var"gm8"=missing, var"dpressure_dpsi"=missing, var"triangularity_upper"=missing, var"darea_drho_tor"=missing, var"area"=missing, var"trapped_fraction"=missing, var"volume"=missing, var"dvolume_dpsi"=missing, var"b_min"=missing, var"f"=missing, var"mass_density"=missing, var"r_outboard"=missing, var"gm4"=missing, var"phi"=missing, var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"gm2"=missing, var"rho_volume_norm"=missing, var"gm1"=missing, var"gm5"=missing, var"b_max"=missing, var"f_df_dpsi"=missing, var"j_tor"=missing, var"r_inboard"=missing, var"q"=missing, var"gm7"=missing, var"pressure"=missing, var"squareness_upper_outer"=missing, var"geometric_axis"=equilibrium__time_slice___profiles_1d__geometric_axis(), _parent=WeakRef(missing))
        fds = new(var"b_field_max", var"dvolume_drho_tor", var"gm9", var"dpsi_drho_tor", var"surface", var"rho_tor", var"magnetic_shear", var"b_average", var"b_field_min", var"darea_dpsi", var"gm3", var"squareness_upper_inner", var"squareness_lower_inner", var"rho_tor_norm", var"elongation", var"beta_pol", var"b_field_average", var"j_parallel", var"gm6", var"psi", var"gm8", var"dpressure_dpsi", var"triangularity_upper", var"darea_drho_tor", var"area", var"trapped_fraction", var"volume", var"dvolume_dpsi", var"b_min", var"f", var"mass_density", var"r_outboard", var"gm4", var"phi", var"squareness_lower_outer", var"triangularity_lower", var"gm2", var"rho_volume_norm", var"gm1", var"gm5", var"b_max", var"f_df_dpsi", var"j_tor", var"r_inboard", var"q", var"gm7", var"pressure", var"squareness_upper_outer", var"geometric_axis", _parent)
        assign_expressions(fds)
        setfield!(fds.geometric_axis, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__q_min <: FDS
    var"value" :: Union{Missing, Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities__q_min(var"value"=missing, var"rho_tor_norm"=missing, _parent=WeakRef(missing))
        fds = new(var"value", var"rho_tor_norm", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__magnetic_axis <: FDS
    var"b_field_tor" :: Union{Missing, Real, Function} = missing
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    var"b_tor" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities__magnetic_axis(var"b_field_tor"=missing, var"r"=missing, var"z"=missing, var"b_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"b_field_tor", var"r", var"z", var"b_tor", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__current_centre <: FDS
    var"velocity_z" :: Union{Missing, Real, Function} = missing
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities__current_centre(var"velocity_z"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"velocity_z", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities <: FDS
    var"ip" :: Union{Missing, Real, Function} = missing
    var"li_3" :: Union{Missing, Real, Function} = missing
    var"beta_tor" :: Union{Missing, Real, Function} = missing
    var"surface" :: Union{Missing, Real, Function} = missing
    var"magnetic_axis" :: equilibrium__time_slice___global_quantities__magnetic_axis = equilibrium__time_slice___global_quantities__magnetic_axis()
    var"energy_mhd" :: Union{Missing, Real, Function} = missing
    var"psi_boundary" :: Union{Missing, Real, Function} = missing
    var"length_pol" :: Union{Missing, Real, Function} = missing
    var"area" :: Union{Missing, Real, Function} = missing
    var"psi_external_average" :: Union{Missing, Real, Function} = missing
    var"q_95" :: Union{Missing, Real, Function} = missing
    var"q_axis" :: Union{Missing, Real, Function} = missing
    var"psi_axis" :: Union{Missing, Real, Function} = missing
    var"w_mhd" :: Union{Missing, Real, Function} = missing
    var"volume" :: Union{Missing, Real, Function} = missing
    var"plasma_inductance" :: Union{Missing, Real, Function} = missing
    var"beta_pol" :: Union{Missing, Real, Function} = missing
    var"beta_normal" :: Union{Missing, Real, Function} = missing
    var"current_centre" :: equilibrium__time_slice___global_quantities__current_centre = equilibrium__time_slice___global_quantities__current_centre()
    var"q_min" :: equilibrium__time_slice___global_quantities__q_min = equilibrium__time_slice___global_quantities__q_min()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities(var"ip"=missing, var"li_3"=missing, var"beta_tor"=missing, var"surface"=missing, var"magnetic_axis"=equilibrium__time_slice___global_quantities__magnetic_axis(), var"energy_mhd"=missing, var"psi_boundary"=missing, var"length_pol"=missing, var"area"=missing, var"psi_external_average"=missing, var"q_95"=missing, var"q_axis"=missing, var"psi_axis"=missing, var"w_mhd"=missing, var"volume"=missing, var"plasma_inductance"=missing, var"beta_pol"=missing, var"beta_normal"=missing, var"current_centre"=equilibrium__time_slice___global_quantities__current_centre(), var"q_min"=equilibrium__time_slice___global_quantities__q_min(), _parent=WeakRef(missing))
        fds = new(var"ip", var"li_3", var"beta_tor", var"surface", var"magnetic_axis", var"energy_mhd", var"psi_boundary", var"length_pol", var"area", var"psi_external_average", var"q_95", var"q_axis", var"psi_axis", var"w_mhd", var"volume", var"plasma_inductance", var"beta_pol", var"beta_normal", var"current_centre", var"q_min", _parent)
        assign_expressions(fds)
        setfield!(fds.magnetic_axis, :_parent, WeakRef(fds))
        setfield!(fds.current_centre, :_parent, WeakRef(fds))
        setfield!(fds.q_min, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___z <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___z(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___theta <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___theta(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___r <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___r(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___psi <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___psi(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___phi <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___phi(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_tor <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___j_tor(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_parallel <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___j_parallel(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary <: FDSvectorElement
    var"neighbours" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"neighbours", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object <: FDSvectorElement
    var"nodes" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"measure" :: Union{Missing, Real, Function} = missing
    var"geometry" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        fds = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_expressions(fds)
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_expressions(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___geometry_type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space <: FDSvectorElement
    var"coordinates_type" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"geometry_type" :: equilibrium__time_slice___ggd___grid__space___geometry_type = equilibrium__time_slice___ggd___grid__space___geometry_type()
    var"identifier" :: equilibrium__time_slice___ggd___grid__space___identifier = equilibrium__time_slice___ggd___grid__space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__time_slice___ggd___grid__space___geometry_type(), var"identifier"=equilibrium__time_slice___ggd___grid__space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[]), _parent=WeakRef(missing))
        fds = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(fds)
        setfield!(fds.geometry_type, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.objects_per_dimension, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element___object <: FDSvectorElement
    var"dimension" :: Union{Missing, Integer, Function} = missing
    var"space" :: Union{Missing, Integer, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"dimension", var"space", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element___object} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_expressions(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___base <: FDSvectorElement
    var"jacobian" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset <: FDSvectorElement
    var"base" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___base} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[])
    var"metric" :: equilibrium__time_slice___ggd___grid__grid_subset___metric = equilibrium__time_slice___ggd___grid__grid_subset___metric()
    var"dimension" :: Union{Missing, Integer, Function} = missing
    var"identifier" :: equilibrium__time_slice___ggd___grid__grid_subset___identifier = equilibrium__time_slice___ggd___grid__grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset(var"base"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[]), var"metric"=equilibrium__time_slice___ggd___grid__grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__time_slice___ggd___grid__grid_subset___identifier(), var"element"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[]), _parent=WeakRef(missing))
        fds = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.base, :_parent, WeakRef(fds))
        setfield!(fds.metric, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space} = FDSvector(equilibrium__time_slice___ggd___grid__space[])
    var"identifier" :: equilibrium__time_slice___ggd___grid__identifier = equilibrium__time_slice___ggd___grid__identifier()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid(var"grid_subset"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset[]), var"space"=FDSvector(equilibrium__time_slice___ggd___grid__space[]), var"identifier"=equilibrium__time_slice___ggd___grid__identifier(), _parent=WeakRef(missing))
        fds = new(var"grid_subset", var"space", var"identifier", _parent)
        assign_expressions(fds)
        setfield!(fds.grid_subset, :_parent, WeakRef(fds))
        setfield!(fds.space, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_z <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___b_field_z(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_tor <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___b_field_tor(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_r <: FDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function} = missing
    var"values" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid_subset_index" :: Union{Missing, Integer, Function} = missing
    var"coefficients" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___b_field_r(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd <: FDSvectorElement
    var"b_field_z" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___b_field_z} = FDSvector(equilibrium__time_slice___ggd___b_field_z[])
    var"psi" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___psi} = FDSvector(equilibrium__time_slice___ggd___psi[])
    var"theta" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___theta} = FDSvector(equilibrium__time_slice___ggd___theta[])
    var"z" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___z} = FDSvector(equilibrium__time_slice___ggd___z[])
    var"phi" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___phi} = FDSvector(equilibrium__time_slice___ggd___phi[])
    var"j_tor" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___j_tor} = FDSvector(equilibrium__time_slice___ggd___j_tor[])
    var"grid" :: equilibrium__time_slice___ggd___grid = equilibrium__time_slice___ggd___grid()
    var"b_field_tor" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___b_field_tor} = FDSvector(equilibrium__time_slice___ggd___b_field_tor[])
    var"b_field_r" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___b_field_r} = FDSvector(equilibrium__time_slice___ggd___b_field_r[])
    var"r" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___r} = FDSvector(equilibrium__time_slice___ggd___r[])
    var"j_parallel" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___j_parallel} = FDSvector(equilibrium__time_slice___ggd___j_parallel[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd(var"b_field_z"=FDSvector(equilibrium__time_slice___ggd___b_field_z[]), var"psi"=FDSvector(equilibrium__time_slice___ggd___psi[]), var"theta"=FDSvector(equilibrium__time_slice___ggd___theta[]), var"z"=FDSvector(equilibrium__time_slice___ggd___z[]), var"phi"=FDSvector(equilibrium__time_slice___ggd___phi[]), var"j_tor"=FDSvector(equilibrium__time_slice___ggd___j_tor[]), var"grid"=equilibrium__time_slice___ggd___grid(), var"b_field_tor"=FDSvector(equilibrium__time_slice___ggd___b_field_tor[]), var"b_field_r"=FDSvector(equilibrium__time_slice___ggd___b_field_r[]), var"r"=FDSvector(equilibrium__time_slice___ggd___r[]), var"j_parallel"=FDSvector(equilibrium__time_slice___ggd___j_parallel[]), _parent=WeakRef(missing))
        fds = new(var"b_field_z", var"psi", var"theta", var"z", var"phi", var"j_tor", var"grid", var"b_field_tor", var"b_field_r", var"r", var"j_parallel", _parent)
        assign_expressions(fds)
        setfield!(fds.b_field_z, :_parent, WeakRef(fds))
        setfield!(fds.psi, :_parent, WeakRef(fds))
        setfield!(fds.theta, :_parent, WeakRef(fds))
        setfield!(fds.z, :_parent, WeakRef(fds))
        setfield!(fds.phi, :_parent, WeakRef(fds))
        setfield!(fds.j_tor, :_parent, WeakRef(fds))
        setfield!(fds.grid, :_parent, WeakRef(fds))
        setfield!(fds.b_field_tor, :_parent, WeakRef(fds))
        setfield!(fds.b_field_r, :_parent, WeakRef(fds))
        setfield!(fds.r, :_parent, WeakRef(fds))
        setfield!(fds.j_parallel, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid_type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___coordinate_system__grid_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid <: FDS
    var"volume_element" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"dim2" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dim1" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___coordinate_system__grid(var"volume_element"=missing, var"dim2"=missing, var"dim1"=missing, _parent=WeakRef(missing))
        fds = new(var"volume_element", var"dim2", var"dim1", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system <: FDS
    var"jacobian" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g13_covariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g11_contravariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g13_contravariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"r" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g12_contravariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g22_contravariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g33_contravariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g22_covariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 4} where T<:Real, Function} = missing
    var"g12_covariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g33_covariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"grid" :: equilibrium__time_slice___coordinate_system__grid = equilibrium__time_slice___coordinate_system__grid()
    var"g23_covariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"g11_covariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 4} where T<:Real, Function} = missing
    var"g23_contravariant" :: Union{Missing, Array{T, 2} where T<:Real, Function} = missing
    var"grid_type" :: equilibrium__time_slice___coordinate_system__grid_type = equilibrium__time_slice___coordinate_system__grid_type()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___coordinate_system(var"jacobian"=missing, var"g13_covariant"=missing, var"g11_contravariant"=missing, var"g13_contravariant"=missing, var"r"=missing, var"g12_contravariant"=missing, var"g22_contravariant"=missing, var"z"=missing, var"g33_contravariant"=missing, var"g22_covariant"=missing, var"tensor_contravariant"=missing, var"g12_covariant"=missing, var"g33_covariant"=missing, var"grid"=equilibrium__time_slice___coordinate_system__grid(), var"g23_covariant"=missing, var"g11_covariant"=missing, var"tensor_covariant"=missing, var"g23_contravariant"=missing, var"grid_type"=equilibrium__time_slice___coordinate_system__grid_type(), _parent=WeakRef(missing))
        fds = new(var"jacobian", var"g13_covariant", var"g11_contravariant", var"g13_contravariant", var"r", var"g12_contravariant", var"g22_contravariant", var"z", var"g33_contravariant", var"g22_covariant", var"tensor_contravariant", var"g12_covariant", var"g33_covariant", var"grid", var"g23_covariant", var"g11_covariant", var"tensor_covariant", var"g23_contravariant", var"grid_type", _parent)
        assign_expressions(fds)
        setfield!(fds.grid, :_parent, WeakRef(fds))
        setfield!(fds.grid_type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___convergence <: FDS
    var"iterations_n" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___convergence(var"iterations_n"=missing, _parent=WeakRef(missing))
        fds = new(var"iterations_n", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_reconstructed <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__x_point___position_reconstructed(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_measured <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__x_point___position_measured(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point <: FDSvectorElement
    var"chi_squared_z" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"position_measured" :: equilibrium__time_slice___constraints__x_point___position_measured = equilibrium__time_slice___constraints__x_point___position_measured()
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    var"chi_squared_r" :: Union{Missing, Real, Function} = missing
    var"position_reconstructed" :: equilibrium__time_slice___constraints__x_point___position_reconstructed = equilibrium__time_slice___constraints__x_point___position_reconstructed()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__x_point(var"chi_squared_z"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"position_measured"=equilibrium__time_slice___constraints__x_point___position_measured(), var"time_measurement"=missing, var"chi_squared_r"=missing, var"position_reconstructed"=equilibrium__time_slice___constraints__x_point___position_reconstructed(), _parent=WeakRef(missing))
        fds = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        assign_expressions(fds)
        setfield!(fds.position_measured, :_parent, WeakRef(fds))
        setfield!(fds.position_reconstructed, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_reconstructed <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__strike_point___position_reconstructed(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_measured <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__strike_point___position_measured(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point <: FDSvectorElement
    var"chi_squared_z" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"position_measured" :: equilibrium__time_slice___constraints__strike_point___position_measured = equilibrium__time_slice___constraints__strike_point___position_measured()
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    var"chi_squared_r" :: Union{Missing, Real, Function} = missing
    var"position_reconstructed" :: equilibrium__time_slice___constraints__strike_point___position_reconstructed = equilibrium__time_slice___constraints__strike_point___position_reconstructed()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__strike_point(var"chi_squared_z"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"position_measured"=equilibrium__time_slice___constraints__strike_point___position_measured(), var"time_measurement"=missing, var"chi_squared_r"=missing, var"position_reconstructed"=equilibrium__time_slice___constraints__strike_point___position_reconstructed(), _parent=WeakRef(missing))
        fds = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        assign_expressions(fds)
        setfield!(fds.position_measured, :_parent, WeakRef(fds))
        setfield!(fds.position_reconstructed, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q___position <: FDS
    var"phi" :: Union{Missing, Real, Function} = missing
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__q___position(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"phi", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"position" :: equilibrium__time_slice___constraints__q___position = equilibrium__time_slice___constraints__q___position()
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__q(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"position"=equilibrium__time_slice___constraints__q___position(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"position", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.position, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pressure <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__pressure(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_passive_current <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__pf_passive_current(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_current <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__pf_current(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e_line <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__n_e_line(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__n_e(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__mse_polarisation_angle <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__mse_polarisation_angle(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z <: FDS
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r <: FDS
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment <: FDSvectorElement
    var"magnetisation_r" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r()
    var"magnetisation_z" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__iron_core_segment(var"magnetisation_r"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(), var"magnetisation_z"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(), _parent=WeakRef(missing))
        fds = new(var"magnetisation_r", var"magnetisation_z", _parent)
        assign_expressions(fds)
        setfield!(fds.magnetisation_r, :_parent, WeakRef(fds))
        setfield!(fds.magnetisation_z, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__ip <: FDS
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__ip(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__flux_loop <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__flux_loop(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__faraday_angle <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__faraday_angle(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__diamagnetic_flux <: FDS
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__diamagnetic_flux(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__bpol_probe <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__bpol_probe(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__b_field_tor_vacuum_r <: FDS
    var"chi_squared" :: Union{Missing, Real, Function} = missing
    var"exact" :: Union{Missing, Integer, Function} = missing
    var"weight" :: Union{Missing, Real, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"measured" :: Union{Missing, Real, Function} = missing
    var"reconstructed" :: Union{Missing, Real, Function} = missing
    var"time_measurement" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__b_field_tor_vacuum_r(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints <: FDS
    var"faraday_angle" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__faraday_angle} = FDSvector(equilibrium__time_slice___constraints__faraday_angle[])
    var"n_e" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__n_e} = FDSvector(equilibrium__time_slice___constraints__n_e[])
    var"ip" :: equilibrium__time_slice___constraints__ip = equilibrium__time_slice___constraints__ip()
    var"n_e_line" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__n_e_line} = FDSvector(equilibrium__time_slice___constraints__n_e_line[])
    var"pf_current" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__pf_current} = FDSvector(equilibrium__time_slice___constraints__pf_current[])
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__strike_point} = FDSvector(equilibrium__time_slice___constraints__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__x_point} = FDSvector(equilibrium__time_slice___constraints__x_point[])
    var"iron_core_segment" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__iron_core_segment} = FDSvector(equilibrium__time_slice___constraints__iron_core_segment[])
    var"pressure" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__pressure} = FDSvector(equilibrium__time_slice___constraints__pressure[])
    var"diamagnetic_flux" :: equilibrium__time_slice___constraints__diamagnetic_flux = equilibrium__time_slice___constraints__diamagnetic_flux()
    var"pf_passive_current" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__pf_passive_current} = FDSvector(equilibrium__time_slice___constraints__pf_passive_current[])
    var"bpol_probe" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__bpol_probe} = FDSvector(equilibrium__time_slice___constraints__bpol_probe[])
    var"mse_polarisation_angle" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__mse_polarisation_angle} = FDSvector(equilibrium__time_slice___constraints__mse_polarisation_angle[])
    var"q" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__q} = FDSvector(equilibrium__time_slice___constraints__q[])
    var"b_field_tor_vacuum_r" :: equilibrium__time_slice___constraints__b_field_tor_vacuum_r = equilibrium__time_slice___constraints__b_field_tor_vacuum_r()
    var"flux_loop" :: FDSvector{T} where {T<:equilibrium__time_slice___constraints__flux_loop} = FDSvector(equilibrium__time_slice___constraints__flux_loop[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints(var"faraday_angle"=FDSvector(equilibrium__time_slice___constraints__faraday_angle[]), var"n_e"=FDSvector(equilibrium__time_slice___constraints__n_e[]), var"ip"=equilibrium__time_slice___constraints__ip(), var"n_e_line"=FDSvector(equilibrium__time_slice___constraints__n_e_line[]), var"pf_current"=FDSvector(equilibrium__time_slice___constraints__pf_current[]), var"strike_point"=FDSvector(equilibrium__time_slice___constraints__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___constraints__x_point[]), var"iron_core_segment"=FDSvector(equilibrium__time_slice___constraints__iron_core_segment[]), var"pressure"=FDSvector(equilibrium__time_slice___constraints__pressure[]), var"diamagnetic_flux"=equilibrium__time_slice___constraints__diamagnetic_flux(), var"pf_passive_current"=FDSvector(equilibrium__time_slice___constraints__pf_passive_current[]), var"bpol_probe"=FDSvector(equilibrium__time_slice___constraints__bpol_probe[]), var"mse_polarisation_angle"=FDSvector(equilibrium__time_slice___constraints__mse_polarisation_angle[]), var"q"=FDSvector(equilibrium__time_slice___constraints__q[]), var"b_field_tor_vacuum_r"=equilibrium__time_slice___constraints__b_field_tor_vacuum_r(), var"flux_loop"=FDSvector(equilibrium__time_slice___constraints__flux_loop[]), _parent=WeakRef(missing))
        fds = new(var"faraday_angle", var"n_e", var"ip", var"n_e_line", var"pf_current", var"strike_point", var"x_point", var"iron_core_segment", var"pressure", var"diamagnetic_flux", var"pf_passive_current", var"bpol_probe", var"mse_polarisation_angle", var"q", var"b_field_tor_vacuum_r", var"flux_loop", _parent)
        assign_expressions(fds)
        setfield!(fds.faraday_angle, :_parent, WeakRef(fds))
        setfield!(fds.n_e, :_parent, WeakRef(fds))
        setfield!(fds.ip, :_parent, WeakRef(fds))
        setfield!(fds.n_e_line, :_parent, WeakRef(fds))
        setfield!(fds.pf_current, :_parent, WeakRef(fds))
        setfield!(fds.strike_point, :_parent, WeakRef(fds))
        setfield!(fds.x_point, :_parent, WeakRef(fds))
        setfield!(fds.iron_core_segment, :_parent, WeakRef(fds))
        setfield!(fds.pressure, :_parent, WeakRef(fds))
        setfield!(fds.diamagnetic_flux, :_parent, WeakRef(fds))
        setfield!(fds.pf_passive_current, :_parent, WeakRef(fds))
        setfield!(fds.bpol_probe, :_parent, WeakRef(fds))
        setfield!(fds.mse_polarisation_angle, :_parent, WeakRef(fds))
        setfield!(fds.q, :_parent, WeakRef(fds))
        setfield!(fds.b_field_tor_vacuum_r, :_parent, WeakRef(fds))
        setfield!(fds.flux_loop, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__x_point <: FDSvectorElement
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__strike_point <: FDSvectorElement
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__outline <: FDS
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__geometric_axis <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__gap <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"r" :: Union{Missing, Real, Function} = missing
    var"value" :: Union{Missing, Real, Function} = missing
    var"identifier" :: Union{Missing, String, Function} = missing
    var"angle" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__gap(var"name"=missing, var"r"=missing, var"value"=missing, var"identifier"=missing, var"angle"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"r", var"value", var"identifier", var"angle", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__closest_wall_point <: FDS
    var"distance" :: Union{Missing, Real, Function} = missing
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__closest_wall_point(var"distance"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"distance", var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__active_limiter_point <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__active_limiter_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix <: FDS
    var"psi" :: Union{Missing, Real, Function} = missing
    var"elongation_lower" :: Union{Missing, Real, Function} = missing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_separatrix__x_point[])
    var"gap" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__gap} = FDSvector(equilibrium__time_slice___boundary_separatrix__gap[])
    var"triangularity" :: Union{Missing, Real, Function} = missing
    var"elongation_upper" :: Union{Missing, Real, Function} = missing
    var"triangularity_upper" :: Union{Missing, Real, Function} = missing
    var"outline" :: equilibrium__time_slice___boundary_separatrix__outline = equilibrium__time_slice___boundary_separatrix__outline()
    var"dr_dz_zero_point" :: equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point = equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point()
    var"squareness_lower_outer" :: Union{Missing, Real, Function} = missing
    var"triangularity_lower" :: Union{Missing, Real, Function} = missing
    var"minor_radius" :: Union{Missing, Real, Function} = missing
    var"squareness_upper_inner" :: Union{Missing, Real, Function} = missing
    var"squareness_upper_outer" :: Union{Missing, Real, Function} = missing
    var"squareness_lower_inner" :: Union{Missing, Real, Function} = missing
    var"geometric_axis" :: equilibrium__time_slice___boundary_separatrix__geometric_axis = equilibrium__time_slice___boundary_separatrix__geometric_axis()
    var"elongation" :: Union{Missing, Real, Function} = missing
    var"active_limiter_point" :: equilibrium__time_slice___boundary_separatrix__active_limiter_point = equilibrium__time_slice___boundary_separatrix__active_limiter_point()
    var"closest_wall_point" :: equilibrium__time_slice___boundary_separatrix__closest_wall_point = equilibrium__time_slice___boundary_separatrix__closest_wall_point()
    var"type" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix(var"psi"=missing, var"elongation_lower"=missing, var"strike_point"=FDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___boundary_separatrix__x_point[]), var"gap"=FDSvector(equilibrium__time_slice___boundary_separatrix__gap[]), var"triangularity"=missing, var"elongation_upper"=missing, var"triangularity_upper"=missing, var"outline"=equilibrium__time_slice___boundary_separatrix__outline(), var"dr_dz_zero_point"=equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(), var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"minor_radius"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"squareness_lower_inner"=missing, var"geometric_axis"=equilibrium__time_slice___boundary_separatrix__geometric_axis(), var"elongation"=missing, var"active_limiter_point"=equilibrium__time_slice___boundary_separatrix__active_limiter_point(), var"closest_wall_point"=equilibrium__time_slice___boundary_separatrix__closest_wall_point(), var"type"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"elongation_lower", var"strike_point", var"x_point", var"gap", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"dr_dz_zero_point", var"squareness_lower_outer", var"triangularity_lower", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"closest_wall_point", var"type", _parent)
        assign_expressions(fds)
        setfield!(fds.strike_point, :_parent, WeakRef(fds))
        setfield!(fds.x_point, :_parent, WeakRef(fds))
        setfield!(fds.gap, :_parent, WeakRef(fds))
        setfield!(fds.outline, :_parent, WeakRef(fds))
        setfield!(fds.dr_dz_zero_point, :_parent, WeakRef(fds))
        setfield!(fds.geometric_axis, :_parent, WeakRef(fds))
        setfield!(fds.active_limiter_point, :_parent, WeakRef(fds))
        setfield!(fds.closest_wall_point, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__x_point <: FDSvectorElement
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__strike_point <: FDSvectorElement
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__outline <: FDS
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix <: FDS
    var"psi" :: Union{Missing, Real, Function} = missing
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[])
    var"distance_inner_outer" :: Union{Missing, Real, Function} = missing
    var"outline" :: equilibrium__time_slice___boundary_secondary_separatrix__outline = equilibrium__time_slice___boundary_secondary_separatrix__outline()
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix(var"psi"=missing, var"x_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[]), var"distance_inner_outer"=missing, var"outline"=equilibrium__time_slice___boundary_secondary_separatrix__outline(), var"strike_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[]), _parent=WeakRef(missing))
        fds = new(var"psi", var"x_point", var"distance_inner_outer", var"outline", var"strike_point", _parent)
        assign_expressions(fds)
        setfield!(fds.x_point, :_parent, WeakRef(fds))
        setfield!(fds.outline, :_parent, WeakRef(fds))
        setfield!(fds.strike_point, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__x_point <: FDSvectorElement
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__strike_point <: FDSvectorElement
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__outline <: FDS
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__lcfs <: FDS
    var"r" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__lcfs(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__geometric_axis <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__active_limiter_point <: FDS
    var"r" :: Union{Missing, Real, Function} = missing
    var"z" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__active_limiter_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary <: FDS
    var"psi" :: Union{Missing, Real, Function} = missing
    var"lcfs" :: equilibrium__time_slice___boundary__lcfs = equilibrium__time_slice___boundary__lcfs()
    var"elongation_lower" :: Union{Missing, Real, Function} = missing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary__strike_point} = FDSvector(equilibrium__time_slice___boundary__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary__x_point} = FDSvector(equilibrium__time_slice___boundary__x_point[])
    var"triangularity" :: Union{Missing, Real, Function} = missing
    var"elongation_upper" :: Union{Missing, Real, Function} = missing
    var"triangularity_upper" :: Union{Missing, Real, Function} = missing
    var"outline" :: equilibrium__time_slice___boundary__outline = equilibrium__time_slice___boundary__outline()
    var"squareness_lower_outer" :: Union{Missing, Real, Function} = missing
    var"triangularity_lower" :: Union{Missing, Real, Function} = missing
    var"psi_norm" :: Union{Missing, Real, Function} = missing
    var"minor_radius" :: Union{Missing, Real, Function} = missing
    var"squareness_upper_inner" :: Union{Missing, Real, Function} = missing
    var"squareness_upper_outer" :: Union{Missing, Real, Function} = missing
    var"squareness_lower_inner" :: Union{Missing, Real, Function} = missing
    var"geometric_axis" :: equilibrium__time_slice___boundary__geometric_axis = equilibrium__time_slice___boundary__geometric_axis()
    var"elongation" :: Union{Missing, Real, Function} = missing
    var"active_limiter_point" :: equilibrium__time_slice___boundary__active_limiter_point = equilibrium__time_slice___boundary__active_limiter_point()
    var"b_flux_pol_norm" :: Union{Missing, Real, Function} = missing
    var"type" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary(var"psi"=missing, var"lcfs"=equilibrium__time_slice___boundary__lcfs(), var"elongation_lower"=missing, var"strike_point"=FDSvector(equilibrium__time_slice___boundary__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___boundary__x_point[]), var"triangularity"=missing, var"elongation_upper"=missing, var"triangularity_upper"=missing, var"outline"=equilibrium__time_slice___boundary__outline(), var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"psi_norm"=missing, var"minor_radius"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"squareness_lower_inner"=missing, var"geometric_axis"=equilibrium__time_slice___boundary__geometric_axis(), var"elongation"=missing, var"active_limiter_point"=equilibrium__time_slice___boundary__active_limiter_point(), var"b_flux_pol_norm"=missing, var"type"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"lcfs", var"elongation_lower", var"strike_point", var"x_point", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"squareness_lower_outer", var"triangularity_lower", var"psi_norm", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"b_flux_pol_norm", var"type", _parent)
        assign_expressions(fds)
        setfield!(fds.lcfs, :_parent, WeakRef(fds))
        setfield!(fds.strike_point, :_parent, WeakRef(fds))
        setfield!(fds.x_point, :_parent, WeakRef(fds))
        setfield!(fds.outline, :_parent, WeakRef(fds))
        setfield!(fds.geometric_axis, :_parent, WeakRef(fds))
        setfield!(fds.active_limiter_point, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice <: FDSvectorElement
    var"time" :: Union{Missing, Real, Function} = missing
    var"ggd" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd} = FDSvector(equilibrium__time_slice___ggd[])
    var"profiles_1d" :: equilibrium__time_slice___profiles_1d = equilibrium__time_slice___profiles_1d()
    var"boundary" :: equilibrium__time_slice___boundary = equilibrium__time_slice___boundary()
    var"constraints" :: equilibrium__time_slice___constraints = equilibrium__time_slice___constraints()
    var"global_quantities" :: equilibrium__time_slice___global_quantities = equilibrium__time_slice___global_quantities()
    var"convergence" :: equilibrium__time_slice___convergence = equilibrium__time_slice___convergence()
    var"coordinate_system" :: equilibrium__time_slice___coordinate_system = equilibrium__time_slice___coordinate_system()
    var"boundary_secondary_separatrix" :: equilibrium__time_slice___boundary_secondary_separatrix = equilibrium__time_slice___boundary_secondary_separatrix()
    var"boundary_separatrix" :: equilibrium__time_slice___boundary_separatrix = equilibrium__time_slice___boundary_separatrix()
    var"profiles_2d" :: FDSvector{T} where {T<:equilibrium__time_slice___profiles_2d} = FDSvector(equilibrium__time_slice___profiles_2d[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice(var"time"=missing, var"ggd"=FDSvector(equilibrium__time_slice___ggd[]), var"profiles_1d"=equilibrium__time_slice___profiles_1d(), var"boundary"=equilibrium__time_slice___boundary(), var"constraints"=equilibrium__time_slice___constraints(), var"global_quantities"=equilibrium__time_slice___global_quantities(), var"convergence"=equilibrium__time_slice___convergence(), var"coordinate_system"=equilibrium__time_slice___coordinate_system(), var"boundary_secondary_separatrix"=equilibrium__time_slice___boundary_secondary_separatrix(), var"boundary_separatrix"=equilibrium__time_slice___boundary_separatrix(), var"profiles_2d"=FDSvector(equilibrium__time_slice___profiles_2d[]), _parent=WeakRef(missing))
        fds = new(var"time", var"ggd", var"profiles_1d", var"boundary", var"constraints", var"global_quantities", var"convergence", var"coordinate_system", var"boundary_secondary_separatrix", var"boundary_separatrix", var"profiles_2d", _parent)
        assign_expressions(fds)
        setfield!(fds.ggd, :_parent, WeakRef(fds))
        setfield!(fds.profiles_1d, :_parent, WeakRef(fds))
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        setfield!(fds.constraints, :_parent, WeakRef(fds))
        setfield!(fds.global_quantities, :_parent, WeakRef(fds))
        setfield!(fds.convergence, :_parent, WeakRef(fds))
        setfield!(fds.coordinate_system, :_parent, WeakRef(fds))
        setfield!(fds.boundary_secondary_separatrix, :_parent, WeakRef(fds))
        setfield!(fds.boundary_separatrix, :_parent, WeakRef(fds))
        setfield!(fds.profiles_2d, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String, Function} = missing
    var"data_dictionary" :: Union{Missing, String, Function} = missing
    var"access_layer" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__ids_properties <: FDS
    var"provider" :: Union{Missing, String, Function} = missing
    var"version_put" :: equilibrium__ids_properties__version_put = equilibrium__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Integer, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"creation_date" :: Union{Missing, String, Function} = missing
    var"comment" :: Union{Missing, String, Function} = missing
    var"occurrence" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__ids_properties(var"provider"=missing, var"version_put"=equilibrium__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary <: FDSvectorElement
    var"neighbours" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"neighbours", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object <: FDSvectorElement
    var"nodes" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"measure" :: Union{Missing, Real, Function} = missing
    var"geometry" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        fds = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_expressions(fds)
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension(var"object"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_expressions(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___geometry_type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space <: FDSvectorElement
    var"coordinates_type" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"geometry_type" :: equilibrium__grids_ggd___grid___space___geometry_type = equilibrium__grids_ggd___grid___space___geometry_type()
    var"identifier" :: equilibrium__grids_ggd___grid___space___identifier = equilibrium__grids_ggd___grid___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__grids_ggd___grid___space___geometry_type(), var"identifier"=equilibrium__grids_ggd___grid___space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[]), _parent=WeakRef(missing))
        fds = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(fds)
        setfield!(fds.geometry_type, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.objects_per_dimension, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element___object <: FDSvectorElement
    var"dimension" :: Union{Missing, Integer, Function} = missing
    var"space" :: Union{Missing, Integer, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"dimension", var"space", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element___object} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___element(var"object"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_expressions(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___base <: FDSvectorElement
    var"jacobian" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"tensor_contravariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    var"tensor_covariant" :: Union{Missing, Array{T, 3} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset <: FDSvectorElement
    var"base" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___base} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[])
    var"metric" :: equilibrium__grids_ggd___grid___grid_subset___metric = equilibrium__grids_ggd___grid___grid_subset___metric()
    var"dimension" :: Union{Missing, Integer, Function} = missing
    var"identifier" :: equilibrium__grids_ggd___grid___grid_subset___identifier = equilibrium__grids_ggd___grid___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset(var"base"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[]), var"metric"=equilibrium__grids_ggd___grid___grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__grids_ggd___grid___grid_subset___identifier(), var"element"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[]), _parent=WeakRef(missing))
        fds = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.base, :_parent, WeakRef(fds))
        setfield!(fds.metric, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid <: FDSvectorElement
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset} = FDSvector(equilibrium__grids_ggd___grid___grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space} = FDSvector(equilibrium__grids_ggd___grid___space[])
    var"identifier" :: equilibrium__grids_ggd___grid___identifier = equilibrium__grids_ggd___grid___identifier()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid(var"grid_subset"=FDSvector(equilibrium__grids_ggd___grid___grid_subset[]), var"space"=FDSvector(equilibrium__grids_ggd___grid___space[]), var"identifier"=equilibrium__grids_ggd___grid___identifier(), _parent=WeakRef(missing))
        fds = new(var"grid_subset", var"space", var"identifier", _parent)
        assign_expressions(fds)
        setfield!(fds.grid_subset, :_parent, WeakRef(fds))
        setfield!(fds.space, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd <: FDSvectorElement
    var"time" :: Union{Missing, Real, Function} = missing
    var"grid" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid} = FDSvector(equilibrium__grids_ggd___grid[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd(var"time"=missing, var"grid"=FDSvector(equilibrium__grids_ggd___grid[]), _parent=WeakRef(missing))
        fds = new(var"time", var"grid", _parent)
        assign_expressions(fds)
        setfield!(fds.grid, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__code <: FDS
    var"library" :: FDSvector{T} where {T<:equilibrium__code__library} = FDSvector(equilibrium__code__library[])
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"output_flag" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__code(var"library"=FDSvector(equilibrium__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium <: FDS
    var"time_slice" :: FDSvector{T} where {T<:equilibrium__time_slice} = FDSvector(equilibrium__time_slice[])
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ids_properties" :: equilibrium__ids_properties = equilibrium__ids_properties()
    var"grids_ggd" :: FDSvector{T} where {T<:equilibrium__grids_ggd} = FDSvector(equilibrium__grids_ggd[])
    var"vacuum_toroidal_field" :: equilibrium__vacuum_toroidal_field = equilibrium__vacuum_toroidal_field()
    var"code" :: equilibrium__code = equilibrium__code()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium(var"time_slice"=FDSvector(equilibrium__time_slice[]), var"time"=missing, var"ids_properties"=equilibrium__ids_properties(), var"grids_ggd"=FDSvector(equilibrium__grids_ggd[]), var"vacuum_toroidal_field"=equilibrium__vacuum_toroidal_field(), var"code"=equilibrium__code(), _parent=WeakRef(missing))
        fds = new(var"time_slice", var"time", var"ids_properties", var"grids_ggd", var"vacuum_toroidal_field", var"code", _parent)
        assign_expressions(fds)
        setfield!(fds.time_slice, :_parent, WeakRef(fds))
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.grids_ggd, :_parent, WeakRef(fds))
        setfield!(fds.vacuum_toroidal_field, :_parent, WeakRef(fds))
        setfield!(fds.code, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__simulation <: FDS
    var"time_ended" :: Union{Missing, String, Function} = missing
    var"time_begun" :: Union{Missing, String, Function} = missing
    var"time_current" :: Union{Missing, Real, Function} = missing
    var"time_restart" :: Union{Missing, Real, Function} = missing
    var"workflow" :: Union{Missing, String, Function} = missing
    var"comment_after" :: Union{Missing, String, Function} = missing
    var"time_begin" :: Union{Missing, Real, Function} = missing
    var"time_end" :: Union{Missing, Real, Function} = missing
    var"comment_before" :: Union{Missing, String, Function} = missing
    var"time_step" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__simulation(var"time_ended"=missing, var"time_begun"=missing, var"time_current"=missing, var"time_restart"=missing, var"workflow"=missing, var"comment_after"=missing, var"time_begin"=missing, var"time_end"=missing, var"comment_before"=missing, var"time_step"=missing, _parent=WeakRef(missing))
        fds = new(var"time_ended", var"time_begun", var"time_current", var"time_restart", var"workflow", var"comment_after", var"time_begin", var"time_end", var"comment_before", var"time_step", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__pulse_time_end_epoch <: FDS
    var"nanoseconds" :: Union{Missing, Integer, Function} = missing
    var"seconds" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__pulse_time_end_epoch(var"nanoseconds"=missing, var"seconds"=missing, _parent=WeakRef(missing))
        fds = new(var"nanoseconds", var"seconds", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__pulse_time_begin_epoch <: FDS
    var"nanoseconds" :: Union{Missing, Integer, Function} = missing
    var"seconds" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__pulse_time_begin_epoch(var"nanoseconds"=missing, var"seconds"=missing, _parent=WeakRef(missing))
        fds = new(var"nanoseconds", var"seconds", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__parent_entry <: FDS
    var"pulse_type" :: Union{Missing, String, Function} = missing
    var"run" :: Union{Missing, Integer, Function} = missing
    var"machine" :: Union{Missing, String, Function} = missing
    var"pulse" :: Union{Missing, Integer, Function} = missing
    var"user" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__parent_entry(var"pulse_type"=missing, var"run"=missing, var"machine"=missing, var"pulse"=missing, var"user"=missing, _parent=WeakRef(missing))
        fds = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String, Function} = missing
    var"data_dictionary" :: Union{Missing, String, Function} = missing
    var"access_layer" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__ids_properties <: FDS
    var"provider" :: Union{Missing, String, Function} = missing
    var"version_put" :: dataset_description__ids_properties__version_put = dataset_description__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Integer, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"creation_date" :: Union{Missing, String, Function} = missing
    var"comment" :: Union{Missing, String, Function} = missing
    var"occurrence" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__ids_properties(var"provider"=missing, var"version_put"=dataset_description__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__data_entry <: FDS
    var"pulse_type" :: Union{Missing, String, Function} = missing
    var"run" :: Union{Missing, Integer, Function} = missing
    var"machine" :: Union{Missing, String, Function} = missing
    var"pulse" :: Union{Missing, Integer, Function} = missing
    var"user" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__data_entry(var"pulse_type"=missing, var"run"=missing, var"machine"=missing, var"pulse"=missing, var"user"=missing, _parent=WeakRef(missing))
        fds = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description <: FDS
    var"pulse_time_begin_epoch" :: dataset_description__pulse_time_begin_epoch = dataset_description__pulse_time_begin_epoch()
    var"imas_version" :: Union{Missing, String, Function} = missing
    var"ids_properties" :: dataset_description__ids_properties = dataset_description__ids_properties()
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"dd_version" :: Union{Missing, String, Function} = missing
    var"parent_entry" :: dataset_description__parent_entry = dataset_description__parent_entry()
    var"simulation" :: dataset_description__simulation = dataset_description__simulation()
    var"pulse_time_end_epoch" :: dataset_description__pulse_time_end_epoch = dataset_description__pulse_time_end_epoch()
    var"data_entry" :: dataset_description__data_entry = dataset_description__data_entry()
    var"pulse_time_begin" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description(var"pulse_time_begin_epoch"=dataset_description__pulse_time_begin_epoch(), var"imas_version"=missing, var"ids_properties"=dataset_description__ids_properties(), var"time"=missing, var"dd_version"=missing, var"parent_entry"=dataset_description__parent_entry(), var"simulation"=dataset_description__simulation(), var"pulse_time_end_epoch"=dataset_description__pulse_time_end_epoch(), var"data_entry"=dataset_description__data_entry(), var"pulse_time_begin"=missing, _parent=WeakRef(missing))
        fds = new(var"pulse_time_begin_epoch", var"imas_version", var"ids_properties", var"time", var"dd_version", var"parent_entry", var"simulation", var"pulse_time_end_epoch", var"data_entry", var"pulse_time_begin", _parent)
        assign_expressions(fds)
        setfield!(fds.pulse_time_begin_epoch, :_parent, WeakRef(fds))
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.parent_entry, :_parent, WeakRef(fds))
        setfield!(fds.simulation, :_parent, WeakRef(fds))
        setfield!(fds.pulse_time_end_epoch, :_parent, WeakRef(fds))
        setfield!(fds.data_entry, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"r0" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        fds = new(var"b0", var"r0", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method = core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___zeff_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method = core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___t_i_average_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___velocity <: FDS
    var"parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"toroidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"radial" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"poloidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___velocity <: FDS
    var"parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"toroidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"radial" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"poloidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___state___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___neutral_type <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___state___neutral_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state <: FDSvectorElement
    var"label" :: Union{Missing, String, Function} = missing
    var"vibrational_level" :: Union{Missing, Real, Function} = missing
    var"temperature" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"electron_configuration" :: Union{Missing, String, Function} = missing
    var"pressure" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"vibrational_mode" :: Union{Missing, String, Function} = missing
    var"pressure_fast_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"velocity" :: core_profiles__profiles_1d___neutral___state___velocity = core_profiles__profiles_1d___neutral___state___velocity()
    var"density" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fast" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"neutral_type" :: core_profiles__profiles_1d___neutral___state___neutral_type = core_profiles__profiles_1d___neutral___state___neutral_type()
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___state(var"label"=missing, var"vibrational_level"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"electron_configuration"=missing, var"pressure"=missing, var"density_thermal"=missing, var"vibrational_mode"=missing, var"pressure_fast_parallel"=missing, var"velocity"=core_profiles__profiles_1d___neutral___state___velocity(), var"density"=missing, var"density_fast"=missing, var"neutral_type"=core_profiles__profiles_1d___neutral___state___neutral_type(), _parent=WeakRef(missing))
        fds = new(var"label", var"vibrational_level", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"velocity", var"density", var"density_fast", var"neutral_type", _parent)
        assign_expressions(fds)
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.neutral_type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___element <: FDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function} = missing
    var"z_n" :: Union{Missing, Real, Function} = missing
    var"multiplicity" :: Union{Missing, Real, Function} = missing
    var"a" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        fds = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral <: FDSvectorElement
    var"label" :: Union{Missing, String, Function} = missing
    var"temperature" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ion_index" :: Union{Missing, Integer, Function} = missing
    var"multiple_states_flag" :: Union{Missing, Integer, Function} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_fast_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral___state} = FDSvector(core_profiles__profiles_1d___neutral___state[])
    var"velocity" :: core_profiles__profiles_1d___neutral___velocity = core_profiles__profiles_1d___neutral___velocity()
    var"density" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fast" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral___element} = FDSvector(core_profiles__profiles_1d___neutral___element[])
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral(var"label"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"ion_index"=missing, var"multiple_states_flag"=missing, var"pressure_fast_perpendicular"=missing, var"pressure"=missing, var"density_thermal"=missing, var"pressure_fast_parallel"=missing, var"state"=FDSvector(core_profiles__profiles_1d___neutral___state[]), var"velocity"=core_profiles__profiles_1d___neutral___velocity(), var"density"=missing, var"density_fast"=missing, var"element"=FDSvector(core_profiles__profiles_1d___neutral___element[]), _parent=WeakRef(missing))
        fds = new(var"label", var"temperature", var"pressure_thermal", var"ion_index", var"multiple_states_flag", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"pressure_fast_parallel", var"state", var"velocity", var"density", var"density_fast", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.state, :_parent, WeakRef(fds))
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___velocity <: FDS
    var"parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"toroidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"radial" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"poloidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___temperature_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___velocity <: FDS
    var"parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"toroidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"radial" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"poloidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state___density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state <: FDSvectorElement
    var"label" :: Union{Missing, String, Function} = missing
    var"vibrational_level" :: Union{Missing, Real, Function} = missing
    var"rotation_frequency_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"temperature" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z_min" :: Union{Missing, Real, Function} = missing
    var"electron_configuration" :: Union{Missing, String, Function} = missing
    var"pressure" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"vibrational_mode" :: Union{Missing, String, Function} = missing
    var"pressure_fast_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z_average_square_1d" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"velocity" :: core_profiles__profiles_1d___ion___state___velocity = core_profiles__profiles_1d___ion___state___velocity()
    var"z_average" :: Union{Missing, Real, Function} = missing
    var"z_max" :: Union{Missing, Real, Function} = missing
    var"z_square_average" :: Union{Missing, Real, Function} = missing
    var"ionisation_potential" :: Union{Missing, Real, Function} = missing
    var"density" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fast" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z_average_1d" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fit" :: core_profiles__profiles_1d___ion___state___density_fit = core_profiles__profiles_1d___ion___state___density_fit()
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state(var"label"=missing, var"vibrational_level"=missing, var"rotation_frequency_tor"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"z_min"=missing, var"electron_configuration"=missing, var"pressure"=missing, var"density_thermal"=missing, var"vibrational_mode"=missing, var"pressure_fast_parallel"=missing, var"z_average_square_1d"=missing, var"velocity"=core_profiles__profiles_1d___ion___state___velocity(), var"z_average"=missing, var"z_max"=missing, var"z_square_average"=missing, var"ionisation_potential"=missing, var"density"=missing, var"density_fast"=missing, var"z_average_1d"=missing, var"density_fit"=core_profiles__profiles_1d___ion___state___density_fit(), _parent=WeakRef(missing))
        fds = new(var"label", var"vibrational_level", var"rotation_frequency_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"z_min", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"z_average_square_1d", var"velocity", var"z_average", var"z_max", var"z_square_average", var"ionisation_potential", var"density", var"density_fast", var"z_average_1d", var"density_fit", _parent)
        assign_expressions(fds)
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.density_fit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___element <: FDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function} = missing
    var"z_n" :: Union{Missing, Real, Function} = missing
    var"multiplicity" :: Union{Missing, Real, Function} = missing
    var"a" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        fds = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion <: FDSvectorElement
    var"label" :: Union{Missing, String, Function} = missing
    var"rotation_frequency_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"temperature_validity" :: Union{Missing, Integer, Function} = missing
    var"velocity_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"temperature" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z_ion_1d" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"multiple_states_flag" :: Union{Missing, Integer, Function} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"neutral_index" :: Union{Missing, Integer, Function} = missing
    var"pressure" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_validity" :: Union{Missing, Integer, Function} = missing
    var"pressure_fast_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion___state} = FDSvector(core_profiles__profiles_1d___ion___state[])
    var"velocity" :: core_profiles__profiles_1d___ion___velocity = core_profiles__profiles_1d___ion___velocity()
    var"z_ion" :: Union{Missing, Real, Function} = missing
    var"temperature_fit" :: core_profiles__profiles_1d___ion___temperature_fit = core_profiles__profiles_1d___ion___temperature_fit()
    var"density" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"velocity_pol" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fast" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fit" :: core_profiles__profiles_1d___ion___density_fit = core_profiles__profiles_1d___ion___density_fit()
    var"z_ion_square_1d" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion___element} = FDSvector(core_profiles__profiles_1d___ion___element[])
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion(var"label"=missing, var"rotation_frequency_tor"=missing, var"temperature_validity"=missing, var"velocity_tor"=missing, var"temperature"=missing, var"z_ion_1d"=missing, var"pressure_thermal"=missing, var"multiple_states_flag"=missing, var"pressure_fast_perpendicular"=missing, var"neutral_index"=missing, var"pressure"=missing, var"density_thermal"=missing, var"density_validity"=missing, var"pressure_fast_parallel"=missing, var"state"=FDSvector(core_profiles__profiles_1d___ion___state[]), var"velocity"=core_profiles__profiles_1d___ion___velocity(), var"z_ion"=missing, var"temperature_fit"=core_profiles__profiles_1d___ion___temperature_fit(), var"density"=missing, var"velocity_pol"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___ion___density_fit(), var"z_ion_square_1d"=missing, var"element"=FDSvector(core_profiles__profiles_1d___ion___element[]), _parent=WeakRef(missing))
        fds = new(var"label", var"rotation_frequency_tor", var"temperature_validity", var"velocity_tor", var"temperature", var"z_ion_1d", var"pressure_thermal", var"multiple_states_flag", var"pressure_fast_perpendicular", var"neutral_index", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"state", var"velocity", var"z_ion", var"temperature_fit", var"density", var"velocity_pol", var"density_fast", var"density_fit", var"z_ion_square_1d", var"element", _parent)
        assign_expressions(fds)
        setfield!(fds.state, :_parent, WeakRef(fds))
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.temperature_fit, :_parent, WeakRef(fds))
        setfield!(fds.density_fit, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___grid <: FDS
    var"psi" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"psi_boundary" :: Union{Missing, Real, Function} = missing
    var"volume" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"area" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_pol_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"surface" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"psi_magnetic_axis" :: Union{Missing, Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___grid(var"psi"=missing, var"psi_boundary"=missing, var"volume"=missing, var"area"=missing, var"rho_pol_norm"=missing, var"rho_tor_norm"=missing, var"surface"=missing, var"rho_tor"=missing, var"psi_magnetic_axis"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"psi_boundary", var"volume", var"area", var"rho_pol_norm", var"rho_tor_norm", var"surface", var"rho_tor", var"psi_magnetic_axis", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__velocity <: FDS
    var"parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"toroidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"radial" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"poloidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__temperature_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String, Function} = missing
    var"description" :: Union{Missing, String, Function} = missing
    var"index" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit <: FDS
    var"local" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"chi_squared" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"reconstructed" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_width" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rho_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"weight" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"source" :: Union{Missing, Array{T, 1} where T<:String, Function} = missing
    var"measured" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method = core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons <: FDS
    var"temperature_validity" :: Union{Missing, Integer, Function} = missing
    var"velocity_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"temperature" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_validity" :: Union{Missing, Integer, Function} = missing
    var"pressure_fast_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"velocity" :: core_profiles__profiles_1d___electrons__velocity = core_profiles__profiles_1d___electrons__velocity()
    var"temperature_fit" :: core_profiles__profiles_1d___electrons__temperature_fit = core_profiles__profiles_1d___electrons__temperature_fit()
    var"density" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"velocity_pol" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"collisionality_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fast" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"density_fit" :: core_profiles__profiles_1d___electrons__density_fit = core_profiles__profiles_1d___electrons__density_fit()
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons(var"temperature_validity"=missing, var"velocity_tor"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"pressure"=missing, var"density_thermal"=missing, var"density_validity"=missing, var"pressure_fast_parallel"=missing, var"velocity"=core_profiles__profiles_1d___electrons__velocity(), var"temperature_fit"=core_profiles__profiles_1d___electrons__temperature_fit(), var"density"=missing, var"velocity_pol"=missing, var"collisionality_norm"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___electrons__density_fit(), _parent=WeakRef(missing))
        fds = new(var"temperature_validity", var"velocity_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"velocity", var"temperature_fit", var"density", var"velocity_pol", var"collisionality_norm", var"density_fast", var"density_fit", _parent)
        assign_expressions(fds)
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.temperature_fit, :_parent, WeakRef(fds))
        setfield!(fds.density_fit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___e_field <: FDS
    var"parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"toroidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"radial" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"poloidal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___e_field(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d <: FDSvectorElement
    var"pressure_ion_total" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"time" :: Union{Missing, Real, Function} = missing
    var"t_i_average_fit" :: core_profiles__profiles_1d___t_i_average_fit = core_profiles__profiles_1d___t_i_average_fit()
    var"neutral" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral} = FDSvector(core_profiles__profiles_1d___neutral[])
    var"n_i_thermal_total" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"magnetic_shear" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ion" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion} = FDSvector(core_profiles__profiles_1d___ion[])
    var"j_total" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"rotation_frequency_tor_sonic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"pressure_thermal" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"j_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"current_parallel_inside" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"j_non_inductive" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"e_field_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"momentum_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"conductivity_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"electrons" :: core_profiles__profiles_1d___electrons = core_profiles__profiles_1d___electrons()
    var"pressure_perpendicular" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"q" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"t_i_average" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"j_ohmic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"grid" :: core_profiles__profiles_1d___grid = core_profiles__profiles_1d___grid()
    var"phi_potential" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"j_bootstrap" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"zeff_fit" :: core_profiles__profiles_1d___zeff_fit = core_profiles__profiles_1d___zeff_fit()
    var"pressure_parallel" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"e_field" :: core_profiles__profiles_1d___e_field = core_profiles__profiles_1d___e_field()
    var"zeff" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"n_i_total_over_n_e" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d(var"pressure_ion_total"=missing, var"time"=missing, var"t_i_average_fit"=core_profiles__profiles_1d___t_i_average_fit(), var"neutral"=FDSvector(core_profiles__profiles_1d___neutral[]), var"n_i_thermal_total"=missing, var"magnetic_shear"=missing, var"ion"=FDSvector(core_profiles__profiles_1d___ion[]), var"j_total"=missing, var"rotation_frequency_tor_sonic"=missing, var"pressure_thermal"=missing, var"j_tor"=missing, var"current_parallel_inside"=missing, var"j_non_inductive"=missing, var"e_field_parallel"=missing, var"momentum_tor"=missing, var"conductivity_parallel"=missing, var"electrons"=core_profiles__profiles_1d___electrons(), var"pressure_perpendicular"=missing, var"q"=missing, var"t_i_average"=missing, var"j_ohmic"=missing, var"grid"=core_profiles__profiles_1d___grid(), var"phi_potential"=missing, var"j_bootstrap"=missing, var"zeff_fit"=core_profiles__profiles_1d___zeff_fit(), var"pressure_parallel"=missing, var"e_field"=core_profiles__profiles_1d___e_field(), var"zeff"=missing, var"n_i_total_over_n_e"=missing, _parent=WeakRef(missing))
        fds = new(var"pressure_ion_total", var"time", var"t_i_average_fit", var"neutral", var"n_i_thermal_total", var"magnetic_shear", var"ion", var"j_total", var"rotation_frequency_tor_sonic", var"pressure_thermal", var"j_tor", var"current_parallel_inside", var"j_non_inductive", var"e_field_parallel", var"momentum_tor", var"conductivity_parallel", var"electrons", var"pressure_perpendicular", var"q", var"t_i_average", var"j_ohmic", var"grid", var"phi_potential", var"j_bootstrap", var"zeff_fit", var"pressure_parallel", var"e_field", var"zeff", var"n_i_total_over_n_e", _parent)
        assign_expressions(fds)
        setfield!(fds.t_i_average_fit, :_parent, WeakRef(fds))
        setfield!(fds.neutral, :_parent, WeakRef(fds))
        setfield!(fds.ion, :_parent, WeakRef(fds))
        setfield!(fds.electrons, :_parent, WeakRef(fds))
        setfield!(fds.grid, :_parent, WeakRef(fds))
        setfield!(fds.zeff_fit, :_parent, WeakRef(fds))
        setfield!(fds.e_field, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String, Function} = missing
    var"data_dictionary" :: Union{Missing, String, Function} = missing
    var"access_layer" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__ids_properties <: FDS
    var"provider" :: Union{Missing, String, Function} = missing
    var"version_put" :: core_profiles__ids_properties__version_put = core_profiles__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Integer, Function} = missing
    var"source" :: Union{Missing, String, Function} = missing
    var"creation_date" :: Union{Missing, String, Function} = missing
    var"comment" :: Union{Missing, String, Function} = missing
    var"occurrence" :: Union{Missing, Integer, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__ids_properties(var"provider"=missing, var"version_put"=core_profiles__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__global_quantities <: FDS
    var"beta_tor_norm" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"resistive_psi_losses" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ip" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"li_3" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"t_i_average_peaking" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"t_e_peaking" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"beta_tor" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"z_eff_resistive" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ejima" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"energy_diamagnetic" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"li" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"current_non_inductive" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"v_loop" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"beta_pol" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"current_bootstrap" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__global_quantities(var"beta_tor_norm"=missing, var"resistive_psi_losses"=missing, var"ip"=missing, var"li_3"=missing, var"t_i_average_peaking"=missing, var"t_e_peaking"=missing, var"beta_tor"=missing, var"z_eff_resistive"=missing, var"ejima"=missing, var"energy_diamagnetic"=missing, var"li"=missing, var"current_non_inductive"=missing, var"v_loop"=missing, var"beta_pol"=missing, var"current_bootstrap"=missing, _parent=WeakRef(missing))
        fds = new(var"beta_tor_norm", var"resistive_psi_losses", var"ip", var"li_3", var"t_i_average_peaking", var"t_e_peaking", var"beta_tor", var"z_eff_resistive", var"ejima", var"energy_diamagnetic", var"li", var"current_non_inductive", var"v_loop", var"beta_pol", var"current_bootstrap", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__code <: FDS
    var"library" :: FDSvector{T} where {T<:core_profiles__code__library} = FDSvector(core_profiles__code__library[])
    var"name" :: Union{Missing, String, Function} = missing
    var"parameters" :: Union{Missing, String, Function} = missing
    var"commit" :: Union{Missing, String, Function} = missing
    var"repository" :: Union{Missing, String, Function} = missing
    var"output_flag" :: Union{Missing, Array{T, 1} where T<:Integer, Function} = missing
    var"version" :: Union{Missing, String, Function} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__code(var"library"=FDSvector(core_profiles__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles <: FDS
    var"time" :: Union{Missing, Array{T, 1} where T<:Real, Function} = missing
    var"ids_properties" :: core_profiles__ids_properties = core_profiles__ids_properties()
    var"vacuum_toroidal_field" :: core_profiles__vacuum_toroidal_field = core_profiles__vacuum_toroidal_field()
    var"code" :: core_profiles__code = core_profiles__code()
    var"global_quantities" :: core_profiles__global_quantities = core_profiles__global_quantities()
    var"profiles_1d" :: FDSvector{T} where {T<:core_profiles__profiles_1d} = FDSvector(core_profiles__profiles_1d[])
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles(var"time"=missing, var"ids_properties"=core_profiles__ids_properties(), var"vacuum_toroidal_field"=core_profiles__vacuum_toroidal_field(), var"code"=core_profiles__code(), var"global_quantities"=core_profiles__global_quantities(), var"profiles_1d"=FDSvector(core_profiles__profiles_1d[]), _parent=WeakRef(missing))
        fds = new(var"time", var"ids_properties", var"vacuum_toroidal_field", var"code", var"global_quantities", var"profiles_1d", _parent)
        assign_expressions(fds)
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.vacuum_toroidal_field, :_parent, WeakRef(fds))
        setfield!(fds.code, :_parent, WeakRef(fds))
        setfield!(fds.global_quantities, :_parent, WeakRef(fds))
        setfield!(fds.profiles_1d, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct dd <: FDS
    var"summary" :: Union{Missing, summary} = summary()
    var"equilibrium" :: Union{Missing, equilibrium} = equilibrium()
    var"core_profiles" :: Union{Missing, core_profiles} = core_profiles()
    var"wall" :: Union{Missing, wall} = wall()
    var"dataset_description" :: Union{Missing, dataset_description} = dataset_description()
    _parent :: WeakRef = WeakRef(missing)
    function dd(var"summary"=summary(), var"equilibrium"=equilibrium(), var"core_profiles"=core_profiles(), var"wall"=wall(), var"dataset_description"=dataset_description(), _parent=WeakRef(missing))
        fds = new(var"summary", var"equilibrium", var"core_profiles", var"wall", var"dataset_description", _parent)
        assign_expressions(fds)
        setfield!(fds.summary, :_parent, WeakRef(fds))
        setfield!(fds.equilibrium, :_parent, WeakRef(fds))
        setfield!(fds.core_profiles, :_parent, WeakRef(fds))
        setfield!(fds.wall, :_parent, WeakRef(fds))
        setfield!(fds.dataset_description, :_parent, WeakRef(fds))
        return fds
    end
end

