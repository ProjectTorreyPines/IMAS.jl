include("functionarrays.jl")

conversion_types = Union{Real, String, Array{Integer, N} where N, Array{Real, N} where N, Array{String, N} where N}

mutable struct wall__temperature_reference <: IDS
    var"data" :: Union{Missing, Real, Function}
    var"description" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__temperature_reference(var"data"=missing, var"description"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"description", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: wall__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__ids_properties(var"provider"=missing, var"version_put"=wall__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__global_quantities__neutral___element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__global_quantities__neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__global_quantities__neutral <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"sputtering_chemical_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"gas_puff" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"recycling_particles_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"pumping_speed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particle_flux_from_wall" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"recycling_energy_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"wall_inventory" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particle_flux_from_plasma" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"sputtering_physical_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"element" :: IDSvector{T} where {T<:wall__global_quantities__neutral___element}
    _parent :: WeakRef
    function wall__global_quantities__neutral(var"label"=missing, var"sputtering_chemical_coefficient"=missing, var"gas_puff"=missing, var"recycling_particles_coefficient"=missing, var"pumping_speed"=missing, var"particle_flux_from_wall"=missing, var"recycling_energy_coefficient"=missing, var"wall_inventory"=missing, var"particle_flux_from_plasma"=missing, var"sputtering_physical_coefficient"=missing, var"element"=IDSvector(wall__global_quantities__neutral___element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"sputtering_chemical_coefficient", var"gas_puff", var"recycling_particles_coefficient", var"pumping_speed", var"particle_flux_from_wall", var"recycling_energy_coefficient", var"wall_inventory", var"particle_flux_from_plasma", var"sputtering_physical_coefficient", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__global_quantities__electrons <: IDS
    var"particle_flux_from_plasma" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gas_puff" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_outer_target" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pumping_speed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particle_flux_from_wall" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"power_inner_target" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__global_quantities__electrons(var"particle_flux_from_plasma"=missing, var"gas_puff"=missing, var"power_outer_target"=missing, var"pumping_speed"=missing, var"particle_flux_from_wall"=missing, var"power_inner_target"=missing, _parent=WeakRef(missing))
        ids = new(var"particle_flux_from_plasma", var"gas_puff", var"power_outer_target", var"pumping_speed", var"particle_flux_from_wall", var"power_inner_target", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__global_quantities <: IDS
    var"neutral" :: IDSvector{T} where {T<:wall__global_quantities__neutral}
    var"power_incident" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_radiated" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_inner_target_ion_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_conducted" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_convected" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electrons" :: wall__global_quantities__electrons
    var"power_density_inner_target_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_black_body" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_recombination_neutrals" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_to_cooling" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_density_outer_target_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_recombination_plasma" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_currents" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_neutrals" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__global_quantities(var"neutral"=IDSvector(wall__global_quantities__neutral[]), var"power_incident"=missing, var"power_radiated"=missing, var"power_inner_target_ion_total"=missing, var"temperature"=missing, var"power_conducted"=missing, var"power_convected"=missing, var"current_tor"=missing, var"electrons"=wall__global_quantities__electrons(), var"power_density_inner_target_max"=missing, var"power_black_body"=missing, var"power_recombination_neutrals"=missing, var"power_to_cooling"=missing, var"power_density_outer_target_max"=missing, var"power_recombination_plasma"=missing, var"power_currents"=missing, var"power_neutrals"=missing, _parent=WeakRef(missing))
        ids = new(var"neutral", var"power_incident", var"power_radiated", var"power_inner_target_ion_total", var"temperature", var"power_conducted", var"power_convected", var"current_tor", var"electrons", var"power_density_inner_target_max", var"power_black_body", var"power_recombination_neutrals", var"power_to_cooling", var"power_density_outer_target_max", var"power_recombination_plasma", var"power_currents", var"power_neutrals", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__first_wall_power_flux_peak <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__first_wall_power_flux_peak(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary <: IDSvectorElement
    var"neighbours" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"neighbours", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object <: IDSvectorElement
    var"nodes" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measure" :: Union{Missing, Real, Function}
    var"geometry" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"boundary" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=IDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        ids = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension <: IDSvectorElement
    var"object" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___objects_per_dimension(var"object"=IDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___geometry_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space <: IDSvectorElement
    var"coordinates_type" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"geometry_type" :: wall__description_ggd___grid_ggd___space___geometry_type
    var"identifier" :: wall__description_ggd___grid_ggd___space___identifier
    var"objects_per_dimension" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space(var"coordinates_type"=missing, var"geometry_type"=wall__description_ggd___grid_ggd___space___geometry_type(), var"identifier"=wall__description_ggd___grid_ggd___space___identifier(), var"objects_per_dimension"=IDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[]), _parent=WeakRef(missing))
        ids = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(ids)
        setfield!(ids.geometry_type, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.objects_per_dimension, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___metric <: IDS
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___element___object <: IDSvectorElement
    var"dimension" :: Union{Missing, Integer, Function}
    var"space" :: Union{Missing, Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"dimension", var"space", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___element <: IDSvectorElement
    var"object" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element___object}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___element(var"object"=IDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___base <: IDSvectorElement
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset <: IDSvectorElement
    var"base" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___base}
    var"metric" :: wall__description_ggd___grid_ggd___grid_subset___metric
    var"dimension" :: Union{Missing, Integer, Function}
    var"identifier" :: wall__description_ggd___grid_ggd___grid_subset___identifier
    var"element" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset(var"base"=IDSvector(wall__description_ggd___grid_ggd___grid_subset___base[]), var"metric"=wall__description_ggd___grid_ggd___grid_subset___metric(), var"dimension"=missing, var"identifier"=wall__description_ggd___grid_ggd___grid_subset___identifier(), var"element"=IDSvector(wall__description_ggd___grid_ggd___grid_subset___element[]), _parent=WeakRef(missing))
        ids = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.base, :_parent, WeakRef(ids))
        setfield!(ids.metric, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd <: IDSvectorElement
    var"time" :: Union{Missing, Real, Function}
    var"grid_subset" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset}
    var"space" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___space}
    var"identifier" :: wall__description_ggd___grid_ggd___identifier
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd(var"time"=missing, var"grid_subset"=IDSvector(wall__description_ggd___grid_ggd___grid_subset[]), var"space"=IDSvector(wall__description_ggd___grid_ggd___space[]), var"identifier"=wall__description_ggd___grid_ggd___identifier(), _parent=WeakRef(missing))
        ids = new(var"time", var"grid_subset", var"space", var"identifier", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_subset, :_parent, WeakRef(ids))
        setfield!(ids.space, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___ggd___temperature <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_ggd___ggd___temperature(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___ggd___power_density <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_ggd___ggd___power_density(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___ggd <: IDSvectorElement
    var"temperature" :: IDSvector{T} where {T<:wall__description_ggd___ggd___temperature}
    var"time" :: Union{Missing, Real, Function}
    var"power_density" :: IDSvector{T} where {T<:wall__description_ggd___ggd___power_density}
    _parent :: WeakRef
    function wall__description_ggd___ggd(var"temperature"=IDSvector(wall__description_ggd___ggd___temperature[]), var"time"=missing, var"power_density"=IDSvector(wall__description_ggd___ggd___power_density[]), _parent=WeakRef(missing))
        ids = new(var"temperature", var"time", var"power_density", _parent)
        assign_expressions(ids)
        setfield!(ids.temperature, :_parent, WeakRef(ids))
        setfield!(ids.power_density, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd <: IDSvectorElement
    var"grid_ggd" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd}
    var"type" :: wall__description_ggd___type
    var"ggd" :: IDSvector{T} where {T<:wall__description_ggd___ggd}
    _parent :: WeakRef
    function wall__description_ggd(var"grid_ggd"=IDSvector(wall__description_ggd___grid_ggd[]), var"type"=wall__description_ggd___type(), var"ggd"=IDSvector(wall__description_ggd___ggd[]), _parent=WeakRef(missing))
        ids = new(var"grid_ggd", var"type", var"ggd", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_ggd, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.ggd, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___element___outline <: IDS
    var"closed" :: Union{Missing, Integer, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___element___outline(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"closed", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___element___j_tor <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___element___j_tor(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___element <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"j_tor" :: wall__description_2d___vessel__unit___element___j_tor
    var"resistivity" :: Union{Missing, Real, Function}
    var"outline" :: wall__description_2d___vessel__unit___element___outline
    var"resistance" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___element(var"name"=missing, var"j_tor"=wall__description_2d___vessel__unit___element___j_tor(), var"resistivity"=missing, var"outline"=wall__description_2d___vessel__unit___element___outline(), var"resistance"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"j_tor", var"resistivity", var"outline", var"resistance", _parent)
        assign_expressions(ids)
        setfield!(ids.j_tor, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___annular__outline_outer <: IDS
    var"closed" :: Union{Missing, Integer, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___annular__outline_outer(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"closed", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___annular__outline_inner <: IDS
    var"closed" :: Union{Missing, Integer, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___annular__outline_inner(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"closed", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___annular__centreline <: IDS
    var"closed" :: Union{Missing, Integer, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___annular__centreline(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"closed", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___annular <: IDS
    var"outline_inner" :: wall__description_2d___vessel__unit___annular__outline_inner
    var"centreline" :: wall__description_2d___vessel__unit___annular__centreline
    var"thickness" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"resistivity" :: Union{Missing, Real, Function}
    var"outline_outer" :: wall__description_2d___vessel__unit___annular__outline_outer
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___annular(var"outline_inner"=wall__description_2d___vessel__unit___annular__outline_inner(), var"centreline"=wall__description_2d___vessel__unit___annular__centreline(), var"thickness"=missing, var"resistivity"=missing, var"outline_outer"=wall__description_2d___vessel__unit___annular__outline_outer(), _parent=WeakRef(missing))
        ids = new(var"outline_inner", var"centreline", var"thickness", var"resistivity", var"outline_outer", _parent)
        assign_expressions(ids)
        setfield!(ids.outline_inner, :_parent, WeakRef(ids))
        setfield!(ids.centreline, :_parent, WeakRef(ids))
        setfield!(ids.outline_outer, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"annular" :: wall__description_2d___vessel__unit___annular
    var"identifier" :: Union{Missing, String, Function}
    var"element" :: IDSvector{T} where {T<:wall__description_2d___vessel__unit___element}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit(var"name"=missing, var"annular"=wall__description_2d___vessel__unit___annular(), var"identifier"=missing, var"element"=IDSvector(wall__description_2d___vessel__unit___element[]), _parent=WeakRef(missing))
        ids = new(var"name", var"annular", var"identifier", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.annular, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___vessel__type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel <: IDS
    var"type" :: wall__description_2d___vessel__type
    var"unit" :: IDSvector{T} where {T<:wall__description_2d___vessel__unit}
    _parent :: WeakRef
    function wall__description_2d___vessel(var"type"=wall__description_2d___vessel__type(), var"unit"=IDSvector(wall__description_2d___vessel__unit[]), _parent=WeakRef(missing))
        ids = new(var"type", var"unit", _parent)
        assign_expressions(ids)
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.unit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_2d___type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___mobile__unit___outline <: IDSvectorElement
    var"time" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___mobile__unit___outline(var"time"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___mobile__unit <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"resistivity" :: Union{Missing, Real, Function}
    var"closed" :: Union{Missing, Integer, Function}
    var"phi_extensions" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"outline" :: IDSvector{T} where {T<:wall__description_2d___mobile__unit___outline}
    _parent :: WeakRef
    function wall__description_2d___mobile__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"phi_extensions"=missing, var"outline"=IDSvector(wall__description_2d___mobile__unit___outline[]), _parent=WeakRef(missing))
        ids = new(var"name", var"resistivity", var"closed", var"phi_extensions", var"outline", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___mobile__type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_2d___mobile__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___mobile <: IDS
    var"type" :: wall__description_2d___mobile__type
    var"unit" :: IDSvector{T} where {T<:wall__description_2d___mobile__unit}
    _parent :: WeakRef
    function wall__description_2d___mobile(var"type"=wall__description_2d___mobile__type(), var"unit"=IDSvector(wall__description_2d___mobile__unit[]), _parent=WeakRef(missing))
        ids = new(var"type", var"unit", _parent)
        assign_expressions(ids)
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.unit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___limiter__unit___outline <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___limiter__unit___outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___limiter__unit <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"resistivity" :: Union{Missing, Real, Function}
    var"closed" :: Union{Missing, Integer, Function}
    var"outline" :: wall__description_2d___limiter__unit___outline
    var"phi_extensions" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___limiter__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"outline"=wall__description_2d___limiter__unit___outline(), var"phi_extensions"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"resistivity", var"closed", var"outline", var"phi_extensions", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___limiter__type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_2d___limiter__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___limiter <: IDS
    var"type" :: wall__description_2d___limiter__type
    var"unit" :: IDSvector{T} where {T<:wall__description_2d___limiter__unit}
    _parent :: WeakRef
    function wall__description_2d___limiter(var"type"=wall__description_2d___limiter__type(), var"unit"=IDSvector(wall__description_2d___limiter__unit[]), _parent=WeakRef(missing))
        ids = new(var"type", var"unit", _parent)
        assign_expressions(ids)
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.unit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d <: IDSvectorElement
    var"mobile" :: wall__description_2d___mobile
    var"limiter" :: wall__description_2d___limiter
    var"type" :: wall__description_2d___type
    var"vessel" :: wall__description_2d___vessel
    _parent :: WeakRef
    function wall__description_2d(var"mobile"=wall__description_2d___mobile(), var"limiter"=wall__description_2d___limiter(), var"type"=wall__description_2d___type(), var"vessel"=wall__description_2d___vessel(), _parent=WeakRef(missing))
        ids = new(var"mobile", var"limiter", var"type", var"vessel", _parent)
        assign_expressions(ids)
        setfield!(ids.mobile, :_parent, WeakRef(ids))
        setfield!(ids.limiter, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.vessel, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__code <: IDS
    var"library" :: IDSvector{T} where {T<:wall__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__code(var"library"=IDSvector(wall__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ids_properties" :: wall__ids_properties
    var"description_ggd" :: IDSvector{T} where {T<:wall__description_ggd}
    var"description_2d" :: IDSvector{T} where {T<:wall__description_2d}
    var"first_wall_surface_area" :: Union{Missing, Real, Function}
    var"code" :: wall__code
    var"global_quantities" :: wall__global_quantities
    var"temperature_reference" :: wall__temperature_reference
    var"first_wall_power_flux_peak" :: wall__first_wall_power_flux_peak
    _parent :: WeakRef
    function wall(var"time"=missing, var"ids_properties"=wall__ids_properties(), var"description_ggd"=IDSvector(wall__description_ggd[]), var"description_2d"=IDSvector(wall__description_2d[]), var"first_wall_surface_area"=missing, var"code"=wall__code(), var"global_quantities"=wall__global_quantities(), var"temperature_reference"=wall__temperature_reference(), var"first_wall_power_flux_peak"=wall__first_wall_power_flux_peak(), _parent=WeakRef(missing))
        ids = new(var"time", var"ids_properties", var"description_ggd", var"description_2d", var"first_wall_surface_area", var"code", var"global_quantities", var"temperature_reference", var"first_wall_power_flux_peak", _parent)
        assign_expressions(ids)
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.description_ggd, :_parent, WeakRef(ids))
        setfield!(ids.description_2d, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.temperature_reference, :_parent, WeakRef(ids))
        setfield!(ids.first_wall_power_flux_peak, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__wall__material <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__wall__material(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__wall__evaporation <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__wall__evaporation(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__wall <: IDS
    var"material" :: summary__wall__material
    var"evaporation" :: summary__wall__evaporation
    _parent :: WeakRef
    function summary__wall(var"material"=summary__wall__material(), var"evaporation"=summary__wall__evaporation(), _parent=WeakRef(missing))
        ids = new(var"material", var"evaporation", _parent)
        assign_expressions(ids)
        setfield!(ids.material, :_parent, WeakRef(ids))
        setfield!(ids.evaporation, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__volume_average__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__n_i <: IDS
    var"krypton" :: summary__volume_average__n_i__krypton
    var"lithium" :: summary__volume_average__n_i__lithium
    var"neon" :: summary__volume_average__n_i__neon
    var"tritium" :: summary__volume_average__n_i__tritium
    var"helium_3" :: summary__volume_average__n_i__helium_3
    var"deuterium" :: summary__volume_average__n_i__deuterium
    var"iron" :: summary__volume_average__n_i__iron
    var"helium_4" :: summary__volume_average__n_i__helium_4
    var"oxygen" :: summary__volume_average__n_i__oxygen
    var"tungsten" :: summary__volume_average__n_i__tungsten
    var"xenon" :: summary__volume_average__n_i__xenon
    var"hydrogen" :: summary__volume_average__n_i__hydrogen
    var"carbon" :: summary__volume_average__n_i__carbon
    var"nitrogen" :: summary__volume_average__n_i__nitrogen
    var"beryllium" :: summary__volume_average__n_i__beryllium
    var"argon" :: summary__volume_average__n_i__argon
    _parent :: WeakRef
    function summary__volume_average__n_i(var"krypton"=summary__volume_average__n_i__krypton(), var"lithium"=summary__volume_average__n_i__lithium(), var"neon"=summary__volume_average__n_i__neon(), var"tritium"=summary__volume_average__n_i__tritium(), var"helium_3"=summary__volume_average__n_i__helium_3(), var"deuterium"=summary__volume_average__n_i__deuterium(), var"iron"=summary__volume_average__n_i__iron(), var"helium_4"=summary__volume_average__n_i__helium_4(), var"oxygen"=summary__volume_average__n_i__oxygen(), var"tungsten"=summary__volume_average__n_i__tungsten(), var"xenon"=summary__volume_average__n_i__xenon(), var"hydrogen"=summary__volume_average__n_i__hydrogen(), var"carbon"=summary__volume_average__n_i__carbon(), var"nitrogen"=summary__volume_average__n_i__nitrogen(), var"beryllium"=summary__volume_average__n_i__beryllium(), var"argon"=summary__volume_average__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__volume_average__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__meff_hydrogenic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__meff_hydrogenic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__isotope_fraction_hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__isotope_fraction_hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average__dn_e_dt <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__volume_average__dn_e_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__volume_average <: IDS
    var"meff_hydrogenic" :: summary__volume_average__meff_hydrogenic
    var"n_i" :: summary__volume_average__n_i
    var"t_e" :: summary__volume_average__t_e
    var"n_e" :: summary__volume_average__n_e
    var"isotope_fraction_hydrogen" :: summary__volume_average__isotope_fraction_hydrogen
    var"t_i_average" :: summary__volume_average__t_i_average
    var"n_i_total" :: summary__volume_average__n_i_total
    var"dn_e_dt" :: summary__volume_average__dn_e_dt
    var"zeff" :: summary__volume_average__zeff
    _parent :: WeakRef
    function summary__volume_average(var"meff_hydrogenic"=summary__volume_average__meff_hydrogenic(), var"n_i"=summary__volume_average__n_i(), var"t_e"=summary__volume_average__t_e(), var"n_e"=summary__volume_average__n_e(), var"isotope_fraction_hydrogen"=summary__volume_average__isotope_fraction_hydrogen(), var"t_i_average"=summary__volume_average__t_i_average(), var"n_i_total"=summary__volume_average__n_i_total(), var"dn_e_dt"=summary__volume_average__dn_e_dt(), var"zeff"=summary__volume_average__zeff(), _parent=WeakRef(missing))
        ids = new(var"meff_hydrogenic", var"n_i", var"t_e", var"n_e", var"isotope_fraction_hydrogen", var"t_i_average", var"n_i_total", var"dn_e_dt", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.meff_hydrogenic, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.isotope_fraction_hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.dn_e_dt, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__tag <: IDS
    var"name" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__tag(var"name"=missing, var"comment"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"comment", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__stationary_phase_flag <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__stationary_phase_flag(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__t_i_average_decay_length <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__t_i_average_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__t_e_decay_length <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__t_e_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__pressure_neutral <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__pressure_neutral(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__power_radiated <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__power_radiated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__n_i_total_decay_length <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__n_i_total_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__n_e_decay_length <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__n_e_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__heat_flux_i_decay_length <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__heat_flux_i_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer__heat_flux_e_decay_length <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__scrape_off_layer__heat_flux_e_decay_length(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__scrape_off_layer <: IDS
    var"heat_flux_i_decay_length" :: summary__scrape_off_layer__heat_flux_i_decay_length
    var"n_e_decay_length" :: summary__scrape_off_layer__n_e_decay_length
    var"pressure_neutral" :: summary__scrape_off_layer__pressure_neutral
    var"heat_flux_e_decay_length" :: summary__scrape_off_layer__heat_flux_e_decay_length
    var"power_radiated" :: summary__scrape_off_layer__power_radiated
    var"n_i_total_decay_length" :: summary__scrape_off_layer__n_i_total_decay_length
    var"t_e_decay_length" :: summary__scrape_off_layer__t_e_decay_length
    var"t_i_average_decay_length" :: summary__scrape_off_layer__t_i_average_decay_length
    _parent :: WeakRef
    function summary__scrape_off_layer(var"heat_flux_i_decay_length"=summary__scrape_off_layer__heat_flux_i_decay_length(), var"n_e_decay_length"=summary__scrape_off_layer__n_e_decay_length(), var"pressure_neutral"=summary__scrape_off_layer__pressure_neutral(), var"heat_flux_e_decay_length"=summary__scrape_off_layer__heat_flux_e_decay_length(), var"power_radiated"=summary__scrape_off_layer__power_radiated(), var"n_i_total_decay_length"=summary__scrape_off_layer__n_i_total_decay_length(), var"t_e_decay_length"=summary__scrape_off_layer__t_e_decay_length(), var"t_i_average_decay_length"=summary__scrape_off_layer__t_i_average_decay_length(), _parent=WeakRef(missing))
        ids = new(var"heat_flux_i_decay_length", var"n_e_decay_length", var"pressure_neutral", var"heat_flux_e_decay_length", var"power_radiated", var"n_i_total_decay_length", var"t_e_decay_length", var"t_i_average_decay_length", _parent)
        assign_expressions(ids)
        setfield!(ids.heat_flux_i_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.n_e_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.pressure_neutral, :_parent, WeakRef(ids))
        setfield!(ids.heat_flux_e_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.t_e_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average_decay_length, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__runaways__particles <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__runaways__particles(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__runaways__current <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__runaways__current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__runaways <: IDS
    var"particles" :: summary__runaways__particles
    var"current" :: summary__runaways__current
    _parent :: WeakRef
    function summary__runaways(var"particles"=summary__runaways__particles(), var"current"=summary__runaways__current(), _parent=WeakRef(missing))
        ids = new(var"particles", var"current", _parent)
        assign_expressions(ids)
        setfield!(ids.particles, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__rmps__occurrence <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__rmps__occurrence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__rmps <: IDS
    var"occurrence" :: summary__rmps__occurrence
    _parent :: WeakRef
    function summary__rmps(var"occurrence"=summary__rmps__occurrence(), _parent=WeakRef(missing))
        ids = new(var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.occurrence, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pellets__occurrence <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__pellets__occurrence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pellets <: IDS
    var"occurrence" :: summary__pellets__occurrence
    _parent :: WeakRef
    function summary__pellets(var"occurrence"=summary__pellets__occurrence(), _parent=WeakRef(missing))
        ids = new(var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.occurrence, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__volume_inside_pedestal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__volume_inside_pedestal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_width <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__pedestal_height <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__offset <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e__d_dpsi_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__t_e <: IDS
    var"pedestal_width" :: summary__pedestal_fits__mtanh__t_e__pedestal_width
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__mtanh__t_e__offset
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm
    var"pedestal_position" :: summary__pedestal_fits__mtanh__t_e__pedestal_position
    var"pedestal_height" :: summary__pedestal_fits__mtanh__t_e__pedestal_height
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e(var"pedestal_width"=summary__pedestal_fits__mtanh__t_e__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__mtanh__t_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm(), var"pedestal_position"=summary__pedestal_fits__mtanh__t_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__mtanh__t_e__pedestal_height(), _parent=WeakRef(missing))
        ids = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(ids)
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter <: IDS
    var"alpha_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical
    var"t_e_pedestal_top_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical
    var"alpha_ratio" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter(var"alpha_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical(), var"t_e_pedestal_top_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical(), var"alpha_ratio"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio(), _parent=WeakRef(missing))
        ids = new(var"alpha_critical", var"t_e_pedestal_top_critical", var"alpha_ratio", _parent)
        assign_expressions(ids)
        setfield!(ids.alpha_critical, :_parent, WeakRef(ids))
        setfield!(ids.t_e_pedestal_top_critical, :_parent, WeakRef(ids))
        setfield!(ids.alpha_ratio, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__bootstrap_current_hager <: IDS
    var"alpha_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical
    var"t_e_pedestal_top_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical
    var"alpha_ratio" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager(var"alpha_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical(), var"t_e_pedestal_top_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical(), var"alpha_ratio"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio(), _parent=WeakRef(missing))
        ids = new(var"alpha_critical", var"t_e_pedestal_top_critical", var"alpha_ratio", _parent)
        assign_expressions(ids)
        setfield!(ids.alpha_critical, :_parent, WeakRef(ids))
        setfield!(ids.t_e_pedestal_top_critical, :_parent, WeakRef(ids))
        setfield!(ids.alpha_ratio, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability__alpha_experimental <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__alpha_experimental(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__stability <: IDS
    var"alpha_experimental" :: summary__pedestal_fits__mtanh__stability__alpha_experimental
    var"bootstrap_current_hager" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager
    var"bootstrap_current_sauter" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability(var"alpha_experimental"=summary__pedestal_fits__mtanh__stability__alpha_experimental(), var"bootstrap_current_hager"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager(), var"bootstrap_current_sauter"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter(), _parent=WeakRef(missing))
        ids = new(var"alpha_experimental", var"bootstrap_current_hager", var"bootstrap_current_sauter", _parent)
        assign_expressions(ids)
        setfield!(ids.alpha_experimental, :_parent, WeakRef(ids))
        setfield!(ids.bootstrap_current_hager, :_parent, WeakRef(ids))
        setfield!(ids.bootstrap_current_sauter, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__separatrix <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_width <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__pedestal_height <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__offset <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__pressure_electron <: IDS
    var"pedestal_width" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_width
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__mtanh__pressure_electron__offset
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm
    var"separatrix" :: summary__pedestal_fits__mtanh__pressure_electron__separatrix
    var"pedestal_position" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_position
    var"pedestal_height" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_height
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron(var"pedestal_width"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__mtanh__pressure_electron__offset(), var"d_dpsi_norm"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__mtanh__pressure_electron__separatrix(), var"pedestal_position"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_height(), _parent=WeakRef(missing))
        ids = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(ids)
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__nustar_pedestal_top_electron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__nustar_pedestal_top_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__separatrix <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_width <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__pedestal_height <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__offset <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e__d_dpsi_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__n_e <: IDS
    var"pedestal_width" :: summary__pedestal_fits__mtanh__n_e__pedestal_width
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__mtanh__n_e__offset
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm
    var"separatrix" :: summary__pedestal_fits__mtanh__n_e__separatrix
    var"pedestal_position" :: summary__pedestal_fits__mtanh__n_e__pedestal_position
    var"pedestal_height" :: summary__pedestal_fits__mtanh__n_e__pedestal_height
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e(var"pedestal_width"=summary__pedestal_fits__mtanh__n_e__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__mtanh__n_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__mtanh__n_e__separatrix(), var"pedestal_position"=summary__pedestal_fits__mtanh__n_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__mtanh__n_e__pedestal_height(), _parent=WeakRef(missing))
        ids = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(ids)
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh__alpha_electron_pedestal_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__alpha_electron_pedestal_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__mtanh <: IDS
    var"pressure_electron" :: summary__pedestal_fits__mtanh__pressure_electron
    var"b_field_pol_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs
    var"alpha_electron_pedestal_max_position" :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position
    var"rhostar_pedestal_top_electron_hfs" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs
    var"beta_pol_pedestal_top_electron_hfs" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs
    var"energy_thermal_pedestal_electron" :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron
    var"parameters" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_e" :: summary__pedestal_fits__mtanh__t_e
    var"rhostar_pedestal_top_electron_lfs" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs
    var"beta_pol_pedestal_top_electron_lfs" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs
    var"stability" :: summary__pedestal_fits__mtanh__stability
    var"b_field_tor_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs
    var"volume_inside_pedestal" :: summary__pedestal_fits__mtanh__volume_inside_pedestal
    var"b_field_pol_pedestal_top_average" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average
    var"coulomb_factor_pedestal_top" :: summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top
    var"beta_pol_pedestal_top_electron_average" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average
    var"b_field_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs
    var"b_field_pol_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs
    var"nustar_pedestal_top_electron" :: summary__pedestal_fits__mtanh__nustar_pedestal_top_electron
    var"alpha_electron_pedestal_max" :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max
    var"b_field_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs
    var"n_e" :: summary__pedestal_fits__mtanh__n_e
    var"energy_thermal_pedestal_ion" :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion
    var"b_field_tor_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs
    var"rhostar_pedestal_top_electron_magnetic_axis" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh(var"pressure_electron"=summary__pedestal_fits__mtanh__pressure_electron(), var"b_field_pol_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs(), var"alpha_electron_pedestal_max_position"=summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position(), var"rhostar_pedestal_top_electron_hfs"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs(), var"beta_pol_pedestal_top_electron_hfs"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs(), var"energy_thermal_pedestal_electron"=summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron(), var"parameters"=missing, var"t_e"=summary__pedestal_fits__mtanh__t_e(), var"rhostar_pedestal_top_electron_lfs"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs(), var"beta_pol_pedestal_top_electron_lfs"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs(), var"stability"=summary__pedestal_fits__mtanh__stability(), var"b_field_tor_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs(), var"volume_inside_pedestal"=summary__pedestal_fits__mtanh__volume_inside_pedestal(), var"b_field_pol_pedestal_top_average"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average(), var"coulomb_factor_pedestal_top"=summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top(), var"beta_pol_pedestal_top_electron_average"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average(), var"b_field_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs(), var"b_field_pol_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs(), var"nustar_pedestal_top_electron"=summary__pedestal_fits__mtanh__nustar_pedestal_top_electron(), var"alpha_electron_pedestal_max"=summary__pedestal_fits__mtanh__alpha_electron_pedestal_max(), var"b_field_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs(), var"n_e"=summary__pedestal_fits__mtanh__n_e(), var"energy_thermal_pedestal_ion"=summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion(), var"b_field_tor_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs(), var"rhostar_pedestal_top_electron_magnetic_axis"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis(), _parent=WeakRef(missing))
        ids = new(var"pressure_electron", var"b_field_pol_pedestal_top_hfs", var"alpha_electron_pedestal_max_position", var"rhostar_pedestal_top_electron_hfs", var"beta_pol_pedestal_top_electron_hfs", var"energy_thermal_pedestal_electron", var"parameters", var"t_e", var"rhostar_pedestal_top_electron_lfs", var"beta_pol_pedestal_top_electron_lfs", var"stability", var"b_field_tor_pedestal_top_lfs", var"volume_inside_pedestal", var"b_field_pol_pedestal_top_average", var"coulomb_factor_pedestal_top", var"beta_pol_pedestal_top_electron_average", var"b_field_pedestal_top_lfs", var"b_field_pol_pedestal_top_lfs", var"nustar_pedestal_top_electron", var"alpha_electron_pedestal_max", var"b_field_pedestal_top_hfs", var"n_e", var"energy_thermal_pedestal_ion", var"b_field_tor_pedestal_top_hfs", var"rhostar_pedestal_top_electron_magnetic_axis", _parent)
        assign_expressions(ids)
        setfield!(ids.pressure_electron, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.alpha_electron_pedestal_max_position, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_electron, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.stability, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.volume_inside_pedestal, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_average, :_parent, WeakRef(ids))
        setfield!(ids.coulomb_factor_pedestal_top, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_average, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.nustar_pedestal_top_electron, :_parent, WeakRef(ids))
        setfield!(ids.alpha_electron_pedestal_max, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_ion, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_magnetic_axis, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__volume_inside_pedestal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__volume_inside_pedestal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e__pedestal_width <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e__pedestal_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e__pedestal_height <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e__offset <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e__d_dpsi_norm_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e__d_dpsi_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__t_e <: IDS
    var"pedestal_width" :: summary__pedestal_fits__linear__t_e__pedestal_width
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__t_e__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__linear__t_e__offset
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__t_e__d_dpsi_norm
    var"pedestal_position" :: summary__pedestal_fits__linear__t_e__pedestal_position
    var"pedestal_height" :: summary__pedestal_fits__linear__t_e__pedestal_height
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e(var"pedestal_width"=summary__pedestal_fits__linear__t_e__pedestal_width(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__t_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__t_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__linear__t_e__d_dpsi_norm(), var"pedestal_position"=summary__pedestal_fits__linear__t_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__linear__t_e__pedestal_height(), _parent=WeakRef(missing))
        ids = new(var"pedestal_width", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(ids)
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__separatrix <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_width <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__pedestal_height <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__offset <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__pressure_electron <: IDS
    var"pedestal_width" :: summary__pedestal_fits__linear__pressure_electron__pedestal_width
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__linear__pressure_electron__offset
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm
    var"separatrix" :: summary__pedestal_fits__linear__pressure_electron__separatrix
    var"pedestal_position" :: summary__pedestal_fits__linear__pressure_electron__pedestal_position
    var"pedestal_height" :: summary__pedestal_fits__linear__pressure_electron__pedestal_height
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron(var"pedestal_width"=summary__pedestal_fits__linear__pressure_electron__pedestal_width(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__pressure_electron__offset(), var"d_dpsi_norm"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__linear__pressure_electron__separatrix(), var"pedestal_position"=summary__pedestal_fits__linear__pressure_electron__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__linear__pressure_electron__pedestal_height(), _parent=WeakRef(missing))
        ids = new(var"pedestal_width", var"d_dpsi_norm_max_position", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(ids)
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__nustar_pedestal_top_electron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__nustar_pedestal_top_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__separatrix <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__separatrix(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__pedestal_width <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__pedestal_width(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__pedestal_position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__pedestal_position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__pedestal_height <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__pedestal_height(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__offset <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__offset(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__d_dpsi_norm_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__d_dpsi_norm_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e__d_dpsi_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e__d_dpsi_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__n_e <: IDS
    var"pedestal_width" :: summary__pedestal_fits__linear__n_e__pedestal_width
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__n_e__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__linear__n_e__offset
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__n_e__d_dpsi_norm
    var"separatrix" :: summary__pedestal_fits__linear__n_e__separatrix
    var"pedestal_position" :: summary__pedestal_fits__linear__n_e__pedestal_position
    var"pedestal_height" :: summary__pedestal_fits__linear__n_e__pedestal_height
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e(var"pedestal_width"=summary__pedestal_fits__linear__n_e__pedestal_width(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__n_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__n_e__offset(), var"d_dpsi_norm"=summary__pedestal_fits__linear__n_e__d_dpsi_norm(), var"separatrix"=summary__pedestal_fits__linear__n_e__separatrix(), var"pedestal_position"=summary__pedestal_fits__linear__n_e__pedestal_position(), var"pedestal_height"=summary__pedestal_fits__linear__n_e__pedestal_height(), _parent=WeakRef(missing))
        ids = new(var"pedestal_width", var"d_dpsi_norm_max", var"offset", var"d_dpsi_norm", var"separatrix", var"pedestal_position", var"pedestal_height", _parent)
        assign_expressions(ids)
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__energy_thermal_pedestal_ion <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__energy_thermal_pedestal_ion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__energy_thermal_pedestal_electron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__energy_thermal_pedestal_electron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__coulomb_factor_pedestal_top <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__coulomb_factor_pedestal_top(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_pol_pedestal_top_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_pol_pedestal_top_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_pedestal_top_lfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_pedestal_top_lfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear__b_field_pedestal_top_hfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__pedestal_fits__linear__b_field_pedestal_top_hfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__pedestal_fits__linear <: IDS
    var"pressure_electron" :: summary__pedestal_fits__linear__pressure_electron
    var"b_field_pol_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs
    var"rhostar_pedestal_top_electron_hfs" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs
    var"beta_pol_pedestal_top_electron_hfs" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs
    var"energy_thermal_pedestal_electron" :: summary__pedestal_fits__linear__energy_thermal_pedestal_electron
    var"parameters" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_e" :: summary__pedestal_fits__linear__t_e
    var"rhostar_pedestal_top_electron_lfs" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs
    var"beta_pol_pedestal_top_electron_lfs" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs
    var"b_field_tor_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs
    var"volume_inside_pedestal" :: summary__pedestal_fits__linear__volume_inside_pedestal
    var"b_field_pol_pedestal_top_average" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_average
    var"coulomb_factor_pedestal_top" :: summary__pedestal_fits__linear__coulomb_factor_pedestal_top
    var"beta_pol_pedestal_top_electron_average" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average
    var"b_field_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_pedestal_top_lfs
    var"b_field_pol_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs
    var"nustar_pedestal_top_electron" :: summary__pedestal_fits__linear__nustar_pedestal_top_electron
    var"b_field_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_pedestal_top_hfs
    var"n_e" :: summary__pedestal_fits__linear__n_e
    var"energy_thermal_pedestal_ion" :: summary__pedestal_fits__linear__energy_thermal_pedestal_ion
    var"b_field_tor_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs
    var"rhostar_pedestal_top_electron_magnetic_axis" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis
    _parent :: WeakRef
    function summary__pedestal_fits__linear(var"pressure_electron"=summary__pedestal_fits__linear__pressure_electron(), var"b_field_pol_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs(), var"rhostar_pedestal_top_electron_hfs"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs(), var"beta_pol_pedestal_top_electron_hfs"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs(), var"energy_thermal_pedestal_electron"=summary__pedestal_fits__linear__energy_thermal_pedestal_electron(), var"parameters"=missing, var"t_e"=summary__pedestal_fits__linear__t_e(), var"rhostar_pedestal_top_electron_lfs"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs(), var"beta_pol_pedestal_top_electron_lfs"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs(), var"b_field_tor_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs(), var"volume_inside_pedestal"=summary__pedestal_fits__linear__volume_inside_pedestal(), var"b_field_pol_pedestal_top_average"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_average(), var"coulomb_factor_pedestal_top"=summary__pedestal_fits__linear__coulomb_factor_pedestal_top(), var"beta_pol_pedestal_top_electron_average"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average(), var"b_field_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_pedestal_top_lfs(), var"b_field_pol_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs(), var"nustar_pedestal_top_electron"=summary__pedestal_fits__linear__nustar_pedestal_top_electron(), var"b_field_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_pedestal_top_hfs(), var"n_e"=summary__pedestal_fits__linear__n_e(), var"energy_thermal_pedestal_ion"=summary__pedestal_fits__linear__energy_thermal_pedestal_ion(), var"b_field_tor_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs(), var"rhostar_pedestal_top_electron_magnetic_axis"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis(), _parent=WeakRef(missing))
        ids = new(var"pressure_electron", var"b_field_pol_pedestal_top_hfs", var"rhostar_pedestal_top_electron_hfs", var"beta_pol_pedestal_top_electron_hfs", var"energy_thermal_pedestal_electron", var"parameters", var"t_e", var"rhostar_pedestal_top_electron_lfs", var"beta_pol_pedestal_top_electron_lfs", var"b_field_tor_pedestal_top_lfs", var"volume_inside_pedestal", var"b_field_pol_pedestal_top_average", var"coulomb_factor_pedestal_top", var"beta_pol_pedestal_top_electron_average", var"b_field_pedestal_top_lfs", var"b_field_pol_pedestal_top_lfs", var"nustar_pedestal_top_electron", var"b_field_pedestal_top_hfs", var"n_e", var"energy_thermal_pedestal_ion", var"b_field_tor_pedestal_top_hfs", var"rhostar_pedestal_top_electron_magnetic_axis", _parent)
        assign_expressions(ids)
        setfield!(ids.pressure_electron, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_electron, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.volume_inside_pedestal, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_average, :_parent, WeakRef(ids))
        setfield!(ids.coulomb_factor_pedestal_top, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_average, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.nustar_pedestal_top_electron, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_ion, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_magnetic_axis, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__pedestal_fits <: IDS
    var"linear" :: summary__pedestal_fits__linear
    var"mtanh" :: summary__pedestal_fits__mtanh
    _parent :: WeakRef
    function summary__pedestal_fits(var"linear"=summary__pedestal_fits__linear(), var"mtanh"=summary__pedestal_fits__mtanh(), _parent=WeakRef(missing))
        ids = new(var"linear", var"mtanh", _parent)
        assign_expressions(ids)
        setfield!(ids.linear, :_parent, WeakRef(ids))
        setfield!(ids.mtanh, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__midplane <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__midplane(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__magnetic_shear_flag <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__magnetic_shear_flag(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__velocity_tor <: IDS
    var"krypton" :: summary__local__separatrix__velocity_tor__krypton
    var"lithium" :: summary__local__separatrix__velocity_tor__lithium
    var"neon" :: summary__local__separatrix__velocity_tor__neon
    var"tritium" :: summary__local__separatrix__velocity_tor__tritium
    var"helium_3" :: summary__local__separatrix__velocity_tor__helium_3
    var"deuterium" :: summary__local__separatrix__velocity_tor__deuterium
    var"iron" :: summary__local__separatrix__velocity_tor__iron
    var"helium_4" :: summary__local__separatrix__velocity_tor__helium_4
    var"oxygen" :: summary__local__separatrix__velocity_tor__oxygen
    var"tungsten" :: summary__local__separatrix__velocity_tor__tungsten
    var"xenon" :: summary__local__separatrix__velocity_tor__xenon
    var"hydrogen" :: summary__local__separatrix__velocity_tor__hydrogen
    var"carbon" :: summary__local__separatrix__velocity_tor__carbon
    var"nitrogen" :: summary__local__separatrix__velocity_tor__nitrogen
    var"beryllium" :: summary__local__separatrix__velocity_tor__beryllium
    var"argon" :: summary__local__separatrix__velocity_tor__argon
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor(var"krypton"=summary__local__separatrix__velocity_tor__krypton(), var"lithium"=summary__local__separatrix__velocity_tor__lithium(), var"neon"=summary__local__separatrix__velocity_tor__neon(), var"tritium"=summary__local__separatrix__velocity_tor__tritium(), var"helium_3"=summary__local__separatrix__velocity_tor__helium_3(), var"deuterium"=summary__local__separatrix__velocity_tor__deuterium(), var"iron"=summary__local__separatrix__velocity_tor__iron(), var"helium_4"=summary__local__separatrix__velocity_tor__helium_4(), var"oxygen"=summary__local__separatrix__velocity_tor__oxygen(), var"tungsten"=summary__local__separatrix__velocity_tor__tungsten(), var"xenon"=summary__local__separatrix__velocity_tor__xenon(), var"hydrogen"=summary__local__separatrix__velocity_tor__hydrogen(), var"carbon"=summary__local__separatrix__velocity_tor__carbon(), var"nitrogen"=summary__local__separatrix__velocity_tor__nitrogen(), var"beryllium"=summary__local__separatrix__velocity_tor__beryllium(), var"argon"=summary__local__separatrix__velocity_tor__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__separatrix__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__q <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__position <: IDS
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__position(var"psi"=missing, var"rho_tor_norm"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"rho_tor_norm", var"rho_tor", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__n_i <: IDS
    var"krypton" :: summary__local__separatrix__n_i__krypton
    var"lithium" :: summary__local__separatrix__n_i__lithium
    var"neon" :: summary__local__separatrix__n_i__neon
    var"tritium" :: summary__local__separatrix__n_i__tritium
    var"helium_3" :: summary__local__separatrix__n_i__helium_3
    var"deuterium" :: summary__local__separatrix__n_i__deuterium
    var"iron" :: summary__local__separatrix__n_i__iron
    var"helium_4" :: summary__local__separatrix__n_i__helium_4
    var"oxygen" :: summary__local__separatrix__n_i__oxygen
    var"tungsten" :: summary__local__separatrix__n_i__tungsten
    var"xenon" :: summary__local__separatrix__n_i__xenon
    var"hydrogen" :: summary__local__separatrix__n_i__hydrogen
    var"carbon" :: summary__local__separatrix__n_i__carbon
    var"nitrogen" :: summary__local__separatrix__n_i__nitrogen
    var"beryllium" :: summary__local__separatrix__n_i__beryllium
    var"argon" :: summary__local__separatrix__n_i__argon
    _parent :: WeakRef
    function summary__local__separatrix__n_i(var"krypton"=summary__local__separatrix__n_i__krypton(), var"lithium"=summary__local__separatrix__n_i__lithium(), var"neon"=summary__local__separatrix__n_i__neon(), var"tritium"=summary__local__separatrix__n_i__tritium(), var"helium_3"=summary__local__separatrix__n_i__helium_3(), var"deuterium"=summary__local__separatrix__n_i__deuterium(), var"iron"=summary__local__separatrix__n_i__iron(), var"helium_4"=summary__local__separatrix__n_i__helium_4(), var"oxygen"=summary__local__separatrix__n_i__oxygen(), var"tungsten"=summary__local__separatrix__n_i__tungsten(), var"xenon"=summary__local__separatrix__n_i__xenon(), var"hydrogen"=summary__local__separatrix__n_i__hydrogen(), var"carbon"=summary__local__separatrix__n_i__carbon(), var"nitrogen"=summary__local__separatrix__n_i__nitrogen(), var"beryllium"=summary__local__separatrix__n_i__beryllium(), var"argon"=summary__local__separatrix__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__separatrix__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__momentum_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__magnetic_shear <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix__e_field_parallel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__separatrix <: IDS
    var"n_i" :: summary__local__separatrix__n_i
    var"velocity_tor" :: summary__local__separatrix__velocity_tor
    var"magnetic_shear" :: summary__local__separatrix__magnetic_shear
    var"t_e" :: summary__local__separatrix__t_e
    var"t_i_average" :: summary__local__separatrix__t_i_average
    var"position" :: summary__local__separatrix__position
    var"e_field_parallel" :: summary__local__separatrix__e_field_parallel
    var"momentum_tor" :: summary__local__separatrix__momentum_tor
    var"n_i_total" :: summary__local__separatrix__n_i_total
    var"q" :: summary__local__separatrix__q
    var"n_e" :: summary__local__separatrix__n_e
    var"zeff" :: summary__local__separatrix__zeff
    _parent :: WeakRef
    function summary__local__separatrix(var"n_i"=summary__local__separatrix__n_i(), var"velocity_tor"=summary__local__separatrix__velocity_tor(), var"magnetic_shear"=summary__local__separatrix__magnetic_shear(), var"t_e"=summary__local__separatrix__t_e(), var"t_i_average"=summary__local__separatrix__t_i_average(), var"position"=summary__local__separatrix__position(), var"e_field_parallel"=summary__local__separatrix__e_field_parallel(), var"momentum_tor"=summary__local__separatrix__momentum_tor(), var"n_i_total"=summary__local__separatrix__n_i_total(), var"q"=summary__local__separatrix__q(), var"n_e"=summary__local__separatrix__n_e(), var"zeff"=summary__local__separatrix__zeff(), _parent=WeakRef(missing))
        ids = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"position", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"q", var"n_e", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__r_eff_norm_2_3__plateau_factor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__r_eff_norm_2_3__plateau_factor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__r_eff_norm_2_3__iota <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__r_eff_norm_2_3__iota(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__r_eff_norm_2_3__effective_helical_ripple <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__r_eff_norm_2_3__effective_helical_ripple(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__r_eff_norm_2_3 <: IDS
    var"plateau_factor" :: summary__local__r_eff_norm_2_3__plateau_factor
    var"iota" :: summary__local__r_eff_norm_2_3__iota
    var"effective_helical_ripple" :: summary__local__r_eff_norm_2_3__effective_helical_ripple
    _parent :: WeakRef
    function summary__local__r_eff_norm_2_3(var"plateau_factor"=summary__local__r_eff_norm_2_3__plateau_factor(), var"iota"=summary__local__r_eff_norm_2_3__iota(), var"effective_helical_ripple"=summary__local__r_eff_norm_2_3__effective_helical_ripple(), _parent=WeakRef(missing))
        ids = new(var"plateau_factor", var"iota", var"effective_helical_ripple", _parent)
        assign_expressions(ids)
        setfield!(ids.plateau_factor, :_parent, WeakRef(ids))
        setfield!(ids.iota, :_parent, WeakRef(ids))
        setfield!(ids.effective_helical_ripple, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__pedestal__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__velocity_tor <: IDS
    var"krypton" :: summary__local__pedestal__velocity_tor__krypton
    var"lithium" :: summary__local__pedestal__velocity_tor__lithium
    var"neon" :: summary__local__pedestal__velocity_tor__neon
    var"tritium" :: summary__local__pedestal__velocity_tor__tritium
    var"helium_3" :: summary__local__pedestal__velocity_tor__helium_3
    var"deuterium" :: summary__local__pedestal__velocity_tor__deuterium
    var"iron" :: summary__local__pedestal__velocity_tor__iron
    var"helium_4" :: summary__local__pedestal__velocity_tor__helium_4
    var"oxygen" :: summary__local__pedestal__velocity_tor__oxygen
    var"tungsten" :: summary__local__pedestal__velocity_tor__tungsten
    var"xenon" :: summary__local__pedestal__velocity_tor__xenon
    var"hydrogen" :: summary__local__pedestal__velocity_tor__hydrogen
    var"carbon" :: summary__local__pedestal__velocity_tor__carbon
    var"nitrogen" :: summary__local__pedestal__velocity_tor__nitrogen
    var"beryllium" :: summary__local__pedestal__velocity_tor__beryllium
    var"argon" :: summary__local__pedestal__velocity_tor__argon
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor(var"krypton"=summary__local__pedestal__velocity_tor__krypton(), var"lithium"=summary__local__pedestal__velocity_tor__lithium(), var"neon"=summary__local__pedestal__velocity_tor__neon(), var"tritium"=summary__local__pedestal__velocity_tor__tritium(), var"helium_3"=summary__local__pedestal__velocity_tor__helium_3(), var"deuterium"=summary__local__pedestal__velocity_tor__deuterium(), var"iron"=summary__local__pedestal__velocity_tor__iron(), var"helium_4"=summary__local__pedestal__velocity_tor__helium_4(), var"oxygen"=summary__local__pedestal__velocity_tor__oxygen(), var"tungsten"=summary__local__pedestal__velocity_tor__tungsten(), var"xenon"=summary__local__pedestal__velocity_tor__xenon(), var"hydrogen"=summary__local__pedestal__velocity_tor__hydrogen(), var"carbon"=summary__local__pedestal__velocity_tor__carbon(), var"nitrogen"=summary__local__pedestal__velocity_tor__nitrogen(), var"beryllium"=summary__local__pedestal__velocity_tor__beryllium(), var"argon"=summary__local__pedestal__velocity_tor__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__pedestal__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__q <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__position <: IDS
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__position(var"psi"=missing, var"rho_tor_norm"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"rho_tor_norm", var"rho_tor", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__n_i <: IDS
    var"krypton" :: summary__local__pedestal__n_i__krypton
    var"lithium" :: summary__local__pedestal__n_i__lithium
    var"neon" :: summary__local__pedestal__n_i__neon
    var"tritium" :: summary__local__pedestal__n_i__tritium
    var"helium_3" :: summary__local__pedestal__n_i__helium_3
    var"deuterium" :: summary__local__pedestal__n_i__deuterium
    var"iron" :: summary__local__pedestal__n_i__iron
    var"helium_4" :: summary__local__pedestal__n_i__helium_4
    var"oxygen" :: summary__local__pedestal__n_i__oxygen
    var"tungsten" :: summary__local__pedestal__n_i__tungsten
    var"xenon" :: summary__local__pedestal__n_i__xenon
    var"hydrogen" :: summary__local__pedestal__n_i__hydrogen
    var"carbon" :: summary__local__pedestal__n_i__carbon
    var"nitrogen" :: summary__local__pedestal__n_i__nitrogen
    var"beryllium" :: summary__local__pedestal__n_i__beryllium
    var"argon" :: summary__local__pedestal__n_i__argon
    _parent :: WeakRef
    function summary__local__pedestal__n_i(var"krypton"=summary__local__pedestal__n_i__krypton(), var"lithium"=summary__local__pedestal__n_i__lithium(), var"neon"=summary__local__pedestal__n_i__neon(), var"tritium"=summary__local__pedestal__n_i__tritium(), var"helium_3"=summary__local__pedestal__n_i__helium_3(), var"deuterium"=summary__local__pedestal__n_i__deuterium(), var"iron"=summary__local__pedestal__n_i__iron(), var"helium_4"=summary__local__pedestal__n_i__helium_4(), var"oxygen"=summary__local__pedestal__n_i__oxygen(), var"tungsten"=summary__local__pedestal__n_i__tungsten(), var"xenon"=summary__local__pedestal__n_i__xenon(), var"hydrogen"=summary__local__pedestal__n_i__hydrogen(), var"carbon"=summary__local__pedestal__n_i__carbon(), var"nitrogen"=summary__local__pedestal__n_i__nitrogen(), var"beryllium"=summary__local__pedestal__n_i__beryllium(), var"argon"=summary__local__pedestal__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__pedestal__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__momentum_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__magnetic_shear <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal__e_field_parallel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__pedestal <: IDS
    var"n_i" :: summary__local__pedestal__n_i
    var"velocity_tor" :: summary__local__pedestal__velocity_tor
    var"magnetic_shear" :: summary__local__pedestal__magnetic_shear
    var"t_e" :: summary__local__pedestal__t_e
    var"t_i_average" :: summary__local__pedestal__t_i_average
    var"position" :: summary__local__pedestal__position
    var"e_field_parallel" :: summary__local__pedestal__e_field_parallel
    var"momentum_tor" :: summary__local__pedestal__momentum_tor
    var"n_i_total" :: summary__local__pedestal__n_i_total
    var"q" :: summary__local__pedestal__q
    var"n_e" :: summary__local__pedestal__n_e
    var"zeff" :: summary__local__pedestal__zeff
    _parent :: WeakRef
    function summary__local__pedestal(var"n_i"=summary__local__pedestal__n_i(), var"velocity_tor"=summary__local__pedestal__velocity_tor(), var"magnetic_shear"=summary__local__pedestal__magnetic_shear(), var"t_e"=summary__local__pedestal__t_e(), var"t_i_average"=summary__local__pedestal__t_i_average(), var"position"=summary__local__pedestal__position(), var"e_field_parallel"=summary__local__pedestal__e_field_parallel(), var"momentum_tor"=summary__local__pedestal__momentum_tor(), var"n_i_total"=summary__local__pedestal__n_i_total(), var"q"=summary__local__pedestal__q(), var"n_e"=summary__local__pedestal__n_e(), var"zeff"=summary__local__pedestal__zeff(), _parent=WeakRef(missing))
        ids = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"position", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"q", var"n_e", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__magnetic_axis__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__velocity_tor <: IDS
    var"krypton" :: summary__local__magnetic_axis__velocity_tor__krypton
    var"lithium" :: summary__local__magnetic_axis__velocity_tor__lithium
    var"neon" :: summary__local__magnetic_axis__velocity_tor__neon
    var"tritium" :: summary__local__magnetic_axis__velocity_tor__tritium
    var"helium_3" :: summary__local__magnetic_axis__velocity_tor__helium_3
    var"deuterium" :: summary__local__magnetic_axis__velocity_tor__deuterium
    var"iron" :: summary__local__magnetic_axis__velocity_tor__iron
    var"helium_4" :: summary__local__magnetic_axis__velocity_tor__helium_4
    var"oxygen" :: summary__local__magnetic_axis__velocity_tor__oxygen
    var"tungsten" :: summary__local__magnetic_axis__velocity_tor__tungsten
    var"xenon" :: summary__local__magnetic_axis__velocity_tor__xenon
    var"hydrogen" :: summary__local__magnetic_axis__velocity_tor__hydrogen
    var"carbon" :: summary__local__magnetic_axis__velocity_tor__carbon
    var"nitrogen" :: summary__local__magnetic_axis__velocity_tor__nitrogen
    var"beryllium" :: summary__local__magnetic_axis__velocity_tor__beryllium
    var"argon" :: summary__local__magnetic_axis__velocity_tor__argon
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor(var"krypton"=summary__local__magnetic_axis__velocity_tor__krypton(), var"lithium"=summary__local__magnetic_axis__velocity_tor__lithium(), var"neon"=summary__local__magnetic_axis__velocity_tor__neon(), var"tritium"=summary__local__magnetic_axis__velocity_tor__tritium(), var"helium_3"=summary__local__magnetic_axis__velocity_tor__helium_3(), var"deuterium"=summary__local__magnetic_axis__velocity_tor__deuterium(), var"iron"=summary__local__magnetic_axis__velocity_tor__iron(), var"helium_4"=summary__local__magnetic_axis__velocity_tor__helium_4(), var"oxygen"=summary__local__magnetic_axis__velocity_tor__oxygen(), var"tungsten"=summary__local__magnetic_axis__velocity_tor__tungsten(), var"xenon"=summary__local__magnetic_axis__velocity_tor__xenon(), var"hydrogen"=summary__local__magnetic_axis__velocity_tor__hydrogen(), var"carbon"=summary__local__magnetic_axis__velocity_tor__carbon(), var"nitrogen"=summary__local__magnetic_axis__velocity_tor__nitrogen(), var"beryllium"=summary__local__magnetic_axis__velocity_tor__beryllium(), var"argon"=summary__local__magnetic_axis__velocity_tor__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__magnetic_axis__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__q <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__position <: IDS
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__position(var"psi"=missing, var"r"=missing, var"rho_tor_norm"=missing, var"z"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"r", var"rho_tor_norm", var"z", var"rho_tor", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_i <: IDS
    var"krypton" :: summary__local__magnetic_axis__n_i__krypton
    var"lithium" :: summary__local__magnetic_axis__n_i__lithium
    var"neon" :: summary__local__magnetic_axis__n_i__neon
    var"tritium" :: summary__local__magnetic_axis__n_i__tritium
    var"helium_3" :: summary__local__magnetic_axis__n_i__helium_3
    var"deuterium" :: summary__local__magnetic_axis__n_i__deuterium
    var"iron" :: summary__local__magnetic_axis__n_i__iron
    var"helium_4" :: summary__local__magnetic_axis__n_i__helium_4
    var"oxygen" :: summary__local__magnetic_axis__n_i__oxygen
    var"tungsten" :: summary__local__magnetic_axis__n_i__tungsten
    var"xenon" :: summary__local__magnetic_axis__n_i__xenon
    var"hydrogen" :: summary__local__magnetic_axis__n_i__hydrogen
    var"carbon" :: summary__local__magnetic_axis__n_i__carbon
    var"nitrogen" :: summary__local__magnetic_axis__n_i__nitrogen
    var"beryllium" :: summary__local__magnetic_axis__n_i__beryllium
    var"argon" :: summary__local__magnetic_axis__n_i__argon
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i(var"krypton"=summary__local__magnetic_axis__n_i__krypton(), var"lithium"=summary__local__magnetic_axis__n_i__lithium(), var"neon"=summary__local__magnetic_axis__n_i__neon(), var"tritium"=summary__local__magnetic_axis__n_i__tritium(), var"helium_3"=summary__local__magnetic_axis__n_i__helium_3(), var"deuterium"=summary__local__magnetic_axis__n_i__deuterium(), var"iron"=summary__local__magnetic_axis__n_i__iron(), var"helium_4"=summary__local__magnetic_axis__n_i__helium_4(), var"oxygen"=summary__local__magnetic_axis__n_i__oxygen(), var"tungsten"=summary__local__magnetic_axis__n_i__tungsten(), var"xenon"=summary__local__magnetic_axis__n_i__xenon(), var"hydrogen"=summary__local__magnetic_axis__n_i__hydrogen(), var"carbon"=summary__local__magnetic_axis__n_i__carbon(), var"nitrogen"=summary__local__magnetic_axis__n_i__nitrogen(), var"beryllium"=summary__local__magnetic_axis__n_i__beryllium(), var"argon"=summary__local__magnetic_axis__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__magnetic_axis__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__momentum_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__magnetic_shear <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__e_field_parallel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis__b_field <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__b_field(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__magnetic_axis <: IDS
    var"n_i" :: summary__local__magnetic_axis__n_i
    var"velocity_tor" :: summary__local__magnetic_axis__velocity_tor
    var"magnetic_shear" :: summary__local__magnetic_axis__magnetic_shear
    var"t_e" :: summary__local__magnetic_axis__t_e
    var"t_i_average" :: summary__local__magnetic_axis__t_i_average
    var"b_field" :: summary__local__magnetic_axis__b_field
    var"q" :: summary__local__magnetic_axis__q
    var"e_field_parallel" :: summary__local__magnetic_axis__e_field_parallel
    var"momentum_tor" :: summary__local__magnetic_axis__momentum_tor
    var"n_i_total" :: summary__local__magnetic_axis__n_i_total
    var"position" :: summary__local__magnetic_axis__position
    var"n_e" :: summary__local__magnetic_axis__n_e
    var"zeff" :: summary__local__magnetic_axis__zeff
    _parent :: WeakRef
    function summary__local__magnetic_axis(var"n_i"=summary__local__magnetic_axis__n_i(), var"velocity_tor"=summary__local__magnetic_axis__velocity_tor(), var"magnetic_shear"=summary__local__magnetic_axis__magnetic_shear(), var"t_e"=summary__local__magnetic_axis__t_e(), var"t_i_average"=summary__local__magnetic_axis__t_i_average(), var"b_field"=summary__local__magnetic_axis__b_field(), var"q"=summary__local__magnetic_axis__q(), var"e_field_parallel"=summary__local__magnetic_axis__e_field_parallel(), var"momentum_tor"=summary__local__magnetic_axis__momentum_tor(), var"n_i_total"=summary__local__magnetic_axis__n_i_total(), var"position"=summary__local__magnetic_axis__position(), var"n_e"=summary__local__magnetic_axis__n_e(), var"zeff"=summary__local__magnetic_axis__zeff(), _parent=WeakRef(missing))
        ids = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"b_field", var"q", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"position", var"n_e", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.b_field, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__limiter__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__power_flux_peak <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__power_flux_peak(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__name <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__local__limiter__name(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__n_i <: IDS
    var"krypton" :: summary__local__limiter__n_i__krypton
    var"lithium" :: summary__local__limiter__n_i__lithium
    var"neon" :: summary__local__limiter__n_i__neon
    var"tritium" :: summary__local__limiter__n_i__tritium
    var"helium_3" :: summary__local__limiter__n_i__helium_3
    var"deuterium" :: summary__local__limiter__n_i__deuterium
    var"iron" :: summary__local__limiter__n_i__iron
    var"helium_4" :: summary__local__limiter__n_i__helium_4
    var"oxygen" :: summary__local__limiter__n_i__oxygen
    var"tungsten" :: summary__local__limiter__n_i__tungsten
    var"xenon" :: summary__local__limiter__n_i__xenon
    var"hydrogen" :: summary__local__limiter__n_i__hydrogen
    var"carbon" :: summary__local__limiter__n_i__carbon
    var"nitrogen" :: summary__local__limiter__n_i__nitrogen
    var"beryllium" :: summary__local__limiter__n_i__beryllium
    var"argon" :: summary__local__limiter__n_i__argon
    _parent :: WeakRef
    function summary__local__limiter__n_i(var"krypton"=summary__local__limiter__n_i__krypton(), var"lithium"=summary__local__limiter__n_i__lithium(), var"neon"=summary__local__limiter__n_i__neon(), var"tritium"=summary__local__limiter__n_i__tritium(), var"helium_3"=summary__local__limiter__n_i__helium_3(), var"deuterium"=summary__local__limiter__n_i__deuterium(), var"iron"=summary__local__limiter__n_i__iron(), var"helium_4"=summary__local__limiter__n_i__helium_4(), var"oxygen"=summary__local__limiter__n_i__oxygen(), var"tungsten"=summary__local__limiter__n_i__tungsten(), var"xenon"=summary__local__limiter__n_i__xenon(), var"hydrogen"=summary__local__limiter__n_i__hydrogen(), var"carbon"=summary__local__limiter__n_i__carbon(), var"nitrogen"=summary__local__limiter__n_i__nitrogen(), var"beryllium"=summary__local__limiter__n_i__beryllium(), var"argon"=summary__local__limiter__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__limiter__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter__flux_expansion <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__limiter__flux_expansion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__limiter <: IDS
    var"flux_expansion" :: summary__local__limiter__flux_expansion
    var"n_i" :: summary__local__limiter__n_i
    var"name" :: summary__local__limiter__name
    var"power_flux_peak" :: summary__local__limiter__power_flux_peak
    var"t_e" :: summary__local__limiter__t_e
    var"n_e" :: summary__local__limiter__n_e
    var"t_i_average" :: summary__local__limiter__t_i_average
    var"n_i_total" :: summary__local__limiter__n_i_total
    var"zeff" :: summary__local__limiter__zeff
    _parent :: WeakRef
    function summary__local__limiter(var"flux_expansion"=summary__local__limiter__flux_expansion(), var"n_i"=summary__local__limiter__n_i(), var"name"=summary__local__limiter__name(), var"power_flux_peak"=summary__local__limiter__power_flux_peak(), var"t_e"=summary__local__limiter__t_e(), var"n_e"=summary__local__limiter__n_e(), var"t_i_average"=summary__local__limiter__t_i_average(), var"n_i_total"=summary__local__limiter__n_i_total(), var"zeff"=summary__local__limiter__zeff(), _parent=WeakRef(missing))
        ids = new(var"flux_expansion", var"n_i", var"name", var"power_flux_peak", var"t_e", var"n_e", var"t_i_average", var"n_i_total", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.flux_expansion, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.name, :_parent, WeakRef(ids))
        setfield!(ids.power_flux_peak, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__itb__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__velocity_tor__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__velocity_tor <: IDS
    var"krypton" :: summary__local__itb__velocity_tor__krypton
    var"lithium" :: summary__local__itb__velocity_tor__lithium
    var"neon" :: summary__local__itb__velocity_tor__neon
    var"tritium" :: summary__local__itb__velocity_tor__tritium
    var"helium_3" :: summary__local__itb__velocity_tor__helium_3
    var"deuterium" :: summary__local__itb__velocity_tor__deuterium
    var"iron" :: summary__local__itb__velocity_tor__iron
    var"helium_4" :: summary__local__itb__velocity_tor__helium_4
    var"oxygen" :: summary__local__itb__velocity_tor__oxygen
    var"tungsten" :: summary__local__itb__velocity_tor__tungsten
    var"xenon" :: summary__local__itb__velocity_tor__xenon
    var"hydrogen" :: summary__local__itb__velocity_tor__hydrogen
    var"carbon" :: summary__local__itb__velocity_tor__carbon
    var"nitrogen" :: summary__local__itb__velocity_tor__nitrogen
    var"beryllium" :: summary__local__itb__velocity_tor__beryllium
    var"argon" :: summary__local__itb__velocity_tor__argon
    _parent :: WeakRef
    function summary__local__itb__velocity_tor(var"krypton"=summary__local__itb__velocity_tor__krypton(), var"lithium"=summary__local__itb__velocity_tor__lithium(), var"neon"=summary__local__itb__velocity_tor__neon(), var"tritium"=summary__local__itb__velocity_tor__tritium(), var"helium_3"=summary__local__itb__velocity_tor__helium_3(), var"deuterium"=summary__local__itb__velocity_tor__deuterium(), var"iron"=summary__local__itb__velocity_tor__iron(), var"helium_4"=summary__local__itb__velocity_tor__helium_4(), var"oxygen"=summary__local__itb__velocity_tor__oxygen(), var"tungsten"=summary__local__itb__velocity_tor__tungsten(), var"xenon"=summary__local__itb__velocity_tor__xenon(), var"hydrogen"=summary__local__itb__velocity_tor__hydrogen(), var"carbon"=summary__local__itb__velocity_tor__carbon(), var"nitrogen"=summary__local__itb__velocity_tor__nitrogen(), var"beryllium"=summary__local__itb__velocity_tor__beryllium(), var"argon"=summary__local__itb__velocity_tor__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__itb__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__q <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__q(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__position <: IDS
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__position(var"psi"=missing, var"rho_tor_norm"=missing, var"rho_tor"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"rho_tor_norm", var"rho_tor", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__n_i <: IDS
    var"krypton" :: summary__local__itb__n_i__krypton
    var"lithium" :: summary__local__itb__n_i__lithium
    var"neon" :: summary__local__itb__n_i__neon
    var"tritium" :: summary__local__itb__n_i__tritium
    var"helium_3" :: summary__local__itb__n_i__helium_3
    var"deuterium" :: summary__local__itb__n_i__deuterium
    var"iron" :: summary__local__itb__n_i__iron
    var"helium_4" :: summary__local__itb__n_i__helium_4
    var"oxygen" :: summary__local__itb__n_i__oxygen
    var"tungsten" :: summary__local__itb__n_i__tungsten
    var"xenon" :: summary__local__itb__n_i__xenon
    var"hydrogen" :: summary__local__itb__n_i__hydrogen
    var"carbon" :: summary__local__itb__n_i__carbon
    var"nitrogen" :: summary__local__itb__n_i__nitrogen
    var"beryllium" :: summary__local__itb__n_i__beryllium
    var"argon" :: summary__local__itb__n_i__argon
    _parent :: WeakRef
    function summary__local__itb__n_i(var"krypton"=summary__local__itb__n_i__krypton(), var"lithium"=summary__local__itb__n_i__lithium(), var"neon"=summary__local__itb__n_i__neon(), var"tritium"=summary__local__itb__n_i__tritium(), var"helium_3"=summary__local__itb__n_i__helium_3(), var"deuterium"=summary__local__itb__n_i__deuterium(), var"iron"=summary__local__itb__n_i__iron(), var"helium_4"=summary__local__itb__n_i__helium_4(), var"oxygen"=summary__local__itb__n_i__oxygen(), var"tungsten"=summary__local__itb__n_i__tungsten(), var"xenon"=summary__local__itb__n_i__xenon(), var"hydrogen"=summary__local__itb__n_i__hydrogen(), var"carbon"=summary__local__itb__n_i__carbon(), var"nitrogen"=summary__local__itb__n_i__nitrogen(), var"beryllium"=summary__local__itb__n_i__beryllium(), var"argon"=summary__local__itb__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__itb__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__momentum_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__momentum_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__magnetic_shear <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__magnetic_shear(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb__e_field_parallel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__e_field_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__itb <: IDS
    var"n_i" :: summary__local__itb__n_i
    var"velocity_tor" :: summary__local__itb__velocity_tor
    var"magnetic_shear" :: summary__local__itb__magnetic_shear
    var"t_e" :: summary__local__itb__t_e
    var"t_i_average" :: summary__local__itb__t_i_average
    var"position" :: summary__local__itb__position
    var"e_field_parallel" :: summary__local__itb__e_field_parallel
    var"momentum_tor" :: summary__local__itb__momentum_tor
    var"n_i_total" :: summary__local__itb__n_i_total
    var"q" :: summary__local__itb__q
    var"n_e" :: summary__local__itb__n_e
    var"zeff" :: summary__local__itb__zeff
    _parent :: WeakRef
    function summary__local__itb(var"n_i"=summary__local__itb__n_i(), var"velocity_tor"=summary__local__itb__velocity_tor(), var"magnetic_shear"=summary__local__itb__magnetic_shear(), var"t_e"=summary__local__itb__t_e(), var"t_i_average"=summary__local__itb__t_i_average(), var"position"=summary__local__itb__position(), var"e_field_parallel"=summary__local__itb__e_field_parallel(), var"momentum_tor"=summary__local__itb__momentum_tor(), var"n_i_total"=summary__local__itb__n_i_total(), var"q"=summary__local__itb__q(), var"n_e"=summary__local__itb__n_e(), var"zeff"=summary__local__itb__zeff(), _parent=WeakRef(missing))
        ids = new(var"n_i", var"velocity_tor", var"magnetic_shear", var"t_e", var"t_i_average", var"position", var"e_field_parallel", var"momentum_tor", var"n_i_total", var"q", var"n_e", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__divertor_plate___zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___power_flux_peak <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___power_flux_peak(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___name <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___name(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_i <: IDS
    var"krypton" :: summary__local__divertor_plate___n_i__krypton
    var"lithium" :: summary__local__divertor_plate___n_i__lithium
    var"neon" :: summary__local__divertor_plate___n_i__neon
    var"tritium" :: summary__local__divertor_plate___n_i__tritium
    var"helium_3" :: summary__local__divertor_plate___n_i__helium_3
    var"deuterium" :: summary__local__divertor_plate___n_i__deuterium
    var"iron" :: summary__local__divertor_plate___n_i__iron
    var"helium_4" :: summary__local__divertor_plate___n_i__helium_4
    var"oxygen" :: summary__local__divertor_plate___n_i__oxygen
    var"tungsten" :: summary__local__divertor_plate___n_i__tungsten
    var"xenon" :: summary__local__divertor_plate___n_i__xenon
    var"hydrogen" :: summary__local__divertor_plate___n_i__hydrogen
    var"carbon" :: summary__local__divertor_plate___n_i__carbon
    var"nitrogen" :: summary__local__divertor_plate___n_i__nitrogen
    var"beryllium" :: summary__local__divertor_plate___n_i__beryllium
    var"argon" :: summary__local__divertor_plate___n_i__argon
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i(var"krypton"=summary__local__divertor_plate___n_i__krypton(), var"lithium"=summary__local__divertor_plate___n_i__lithium(), var"neon"=summary__local__divertor_plate___n_i__neon(), var"tritium"=summary__local__divertor_plate___n_i__tritium(), var"helium_3"=summary__local__divertor_plate___n_i__helium_3(), var"deuterium"=summary__local__divertor_plate___n_i__deuterium(), var"iron"=summary__local__divertor_plate___n_i__iron(), var"helium_4"=summary__local__divertor_plate___n_i__helium_4(), var"oxygen"=summary__local__divertor_plate___n_i__oxygen(), var"tungsten"=summary__local__divertor_plate___n_i__tungsten(), var"xenon"=summary__local__divertor_plate___n_i__xenon(), var"hydrogen"=summary__local__divertor_plate___n_i__hydrogen(), var"carbon"=summary__local__divertor_plate___n_i__carbon(), var"nitrogen"=summary__local__divertor_plate___n_i__nitrogen(), var"beryllium"=summary__local__divertor_plate___n_i__beryllium(), var"argon"=summary__local__divertor_plate___n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local__divertor_plate___n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate___flux_expansion <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__divertor_plate___flux_expansion(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__local__divertor_plate <: IDSvectorElement
    var"flux_expansion" :: summary__local__divertor_plate___flux_expansion
    var"n_i" :: summary__local__divertor_plate___n_i
    var"name" :: summary__local__divertor_plate___name
    var"power_flux_peak" :: summary__local__divertor_plate___power_flux_peak
    var"t_e" :: summary__local__divertor_plate___t_e
    var"n_e" :: summary__local__divertor_plate___n_e
    var"t_i_average" :: summary__local__divertor_plate___t_i_average
    var"n_i_total" :: summary__local__divertor_plate___n_i_total
    var"zeff" :: summary__local__divertor_plate___zeff
    _parent :: WeakRef
    function summary__local__divertor_plate(var"flux_expansion"=summary__local__divertor_plate___flux_expansion(), var"n_i"=summary__local__divertor_plate___n_i(), var"name"=summary__local__divertor_plate___name(), var"power_flux_peak"=summary__local__divertor_plate___power_flux_peak(), var"t_e"=summary__local__divertor_plate___t_e(), var"n_e"=summary__local__divertor_plate___n_e(), var"t_i_average"=summary__local__divertor_plate___t_i_average(), var"n_i_total"=summary__local__divertor_plate___n_i_total(), var"zeff"=summary__local__divertor_plate___zeff(), _parent=WeakRef(missing))
        ids = new(var"flux_expansion", var"n_i", var"name", var"power_flux_peak", var"t_e", var"n_e", var"t_i_average", var"n_i_total", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.flux_expansion, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.name, :_parent, WeakRef(ids))
        setfield!(ids.power_flux_peak, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local <: IDS
    var"itb" :: summary__local__itb
    var"r_eff_norm_2_3" :: summary__local__r_eff_norm_2_3
    var"limiter" :: summary__local__limiter
    var"divertor_plate" :: IDSvector{T} where {T<:summary__local__divertor_plate}
    var"magnetic_axis" :: summary__local__magnetic_axis
    var"pedestal" :: summary__local__pedestal
    var"separatrix" :: summary__local__separatrix
    _parent :: WeakRef
    function summary__local(var"itb"=summary__local__itb(), var"r_eff_norm_2_3"=summary__local__r_eff_norm_2_3(), var"limiter"=summary__local__limiter(), var"divertor_plate"=IDSvector(summary__local__divertor_plate[]), var"magnetic_axis"=summary__local__magnetic_axis(), var"pedestal"=summary__local__pedestal(), var"separatrix"=summary__local__separatrix(), _parent=WeakRef(missing))
        ids = new(var"itb", var"r_eff_norm_2_3", var"limiter", var"divertor_plate", var"magnetic_axis", var"pedestal", var"separatrix", _parent)
        assign_expressions(ids)
        setfield!(ids.itb, :_parent, WeakRef(ids))
        setfield!(ids.r_eff_norm_2_3, :_parent, WeakRef(ids))
        setfield!(ids.limiter, :_parent, WeakRef(ids))
        setfield!(ids.divertor_plate, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis, :_parent, WeakRef(ids))
        setfield!(ids.pedestal, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__line_average__zeff <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__zeff(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__t_i_average <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__t_i_average(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__t_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__t_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__tungsten <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__tungsten(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__iron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__iron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_i__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__n_i <: IDS
    var"krypton" :: summary__line_average__n_i__krypton
    var"lithium" :: summary__line_average__n_i__lithium
    var"neon" :: summary__line_average__n_i__neon
    var"tritium" :: summary__line_average__n_i__tritium
    var"helium_3" :: summary__line_average__n_i__helium_3
    var"deuterium" :: summary__line_average__n_i__deuterium
    var"iron" :: summary__line_average__n_i__iron
    var"helium_4" :: summary__line_average__n_i__helium_4
    var"oxygen" :: summary__line_average__n_i__oxygen
    var"tungsten" :: summary__line_average__n_i__tungsten
    var"xenon" :: summary__line_average__n_i__xenon
    var"hydrogen" :: summary__line_average__n_i__hydrogen
    var"carbon" :: summary__line_average__n_i__carbon
    var"nitrogen" :: summary__line_average__n_i__nitrogen
    var"beryllium" :: summary__line_average__n_i__beryllium
    var"argon" :: summary__line_average__n_i__argon
    _parent :: WeakRef
    function summary__line_average__n_i(var"krypton"=summary__line_average__n_i__krypton(), var"lithium"=summary__line_average__n_i__lithium(), var"neon"=summary__line_average__n_i__neon(), var"tritium"=summary__line_average__n_i__tritium(), var"helium_3"=summary__line_average__n_i__helium_3(), var"deuterium"=summary__line_average__n_i__deuterium(), var"iron"=summary__line_average__n_i__iron(), var"helium_4"=summary__line_average__n_i__helium_4(), var"oxygen"=summary__line_average__n_i__oxygen(), var"tungsten"=summary__line_average__n_i__tungsten(), var"xenon"=summary__line_average__n_i__xenon(), var"hydrogen"=summary__line_average__n_i__hydrogen(), var"carbon"=summary__line_average__n_i__carbon(), var"nitrogen"=summary__line_average__n_i__nitrogen(), var"beryllium"=summary__line_average__n_i__beryllium(), var"argon"=summary__line_average__n_i__argon(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"neon", var"tritium", var"helium_3", var"deuterium", var"iron", var"helium_4", var"oxygen", var"tungsten", var"xenon", var"hydrogen", var"carbon", var"nitrogen", var"beryllium", var"argon", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__line_average__n_e <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__n_e(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__meff_hydrogenic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__meff_hydrogenic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__isotope_fraction_hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__isotope_fraction_hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average__dn_e_dt <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__line_average__dn_e_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__line_average <: IDS
    var"meff_hydrogenic" :: summary__line_average__meff_hydrogenic
    var"n_i" :: summary__line_average__n_i
    var"t_e" :: summary__line_average__t_e
    var"n_e" :: summary__line_average__n_e
    var"isotope_fraction_hydrogen" :: summary__line_average__isotope_fraction_hydrogen
    var"t_i_average" :: summary__line_average__t_i_average
    var"n_i_total" :: summary__line_average__n_i_total
    var"dn_e_dt" :: summary__line_average__dn_e_dt
    var"zeff" :: summary__line_average__zeff
    _parent :: WeakRef
    function summary__line_average(var"meff_hydrogenic"=summary__line_average__meff_hydrogenic(), var"n_i"=summary__line_average__n_i(), var"t_e"=summary__line_average__t_e(), var"n_e"=summary__line_average__n_e(), var"isotope_fraction_hydrogen"=summary__line_average__isotope_fraction_hydrogen(), var"t_i_average"=summary__line_average__t_i_average(), var"n_i_total"=summary__line_average__n_i_total(), var"dn_e_dt"=summary__line_average__dn_e_dt(), var"zeff"=summary__line_average__zeff(), _parent=WeakRef(missing))
        ids = new(var"meff_hydrogenic", var"n_i", var"t_e", var"n_e", var"isotope_fraction_hydrogen", var"t_i_average", var"n_i_total", var"dn_e_dt", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.meff_hydrogenic, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.isotope_fraction_hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.dn_e_dt, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__limiter__material <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__limiter__material(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__limiter <: IDS
    var"material" :: summary__limiter__material
    _parent :: WeakRef
    function summary__limiter(var"material"=summary__limiter__material(), _parent=WeakRef(missing))
        ids = new(var"material", _parent)
        assign_expressions(ids)
        setfield!(ids.material, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__kicks__occurrence <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__kicks__occurrence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__kicks <: IDS
    var"occurrence" :: summary__kicks__occurrence
    _parent :: WeakRef
    function summary__kicks(var"occurrence"=summary__kicks__occurrence(), _parent=WeakRef(missing))
        ids = new(var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.occurrence, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: summary__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__ids_properties(var"provider"=missing, var"version_put"=summary__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive__power_nbi <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_nbi(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_lh <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_lh(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_launched_nbi_co_injected_ratio <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_launched_nbi_co_injected_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_launched_nbi <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_launched_nbi(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_launched_lh <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_launched_lh(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_launched_ic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_launched_ic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_launched_ec <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_launched_ec(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_ic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_ic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_ec <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_ec(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__power_additional <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__power_additional(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___tangency_radius <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___tangency_radius(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___species__z_n <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___species__z_n(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___species__label <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___species__label(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___species__a <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___species__a(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___species <: IDS
    var"label" :: summary__heating_current_drive__nbi___species__label
    var"z_n" :: summary__heating_current_drive__nbi___species__z_n
    var"a" :: summary__heating_current_drive__nbi___species__a
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___species(var"label"=summary__heating_current_drive__nbi___species__label(), var"z_n"=summary__heating_current_drive__nbi___species__z_n(), var"a"=summary__heating_current_drive__nbi___species__a(), _parent=WeakRef(missing))
        ids = new(var"label", var"z_n", var"a", _parent)
        assign_expressions(ids)
        setfield!(ids.label, :_parent, WeakRef(ids))
        setfield!(ids.z_n, :_parent, WeakRef(ids))
        setfield!(ids.a, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___power_launched <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___power <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___position__z <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___position__z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___position__r <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___position__r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___position__phi <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___position__phi(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___position <: IDS
    var"phi" :: summary__heating_current_drive__nbi___position__phi
    var"r" :: summary__heating_current_drive__nbi___position__r
    var"z" :: summary__heating_current_drive__nbi___position__z
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___position(var"phi"=summary__heating_current_drive__nbi___position__phi(), var"r"=summary__heating_current_drive__nbi___position__r(), var"z"=summary__heating_current_drive__nbi___position__z(), _parent=WeakRef(missing))
        ids = new(var"phi", var"r", var"z", _parent)
        assign_expressions(ids)
        setfield!(ids.phi, :_parent, WeakRef(ids))
        setfield!(ids.r, :_parent, WeakRef(ids))
        setfield!(ids.z, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___energy <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___energy(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___direction <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___direction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___current <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___beam_power_fraction <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___beam_power_fraction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___beam_current_fraction <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___beam_current_fraction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi___angle <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___angle(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__nbi <: IDSvectorElement
    var"energy" :: summary__heating_current_drive__nbi___energy
    var"angle" :: summary__heating_current_drive__nbi___angle
    var"current" :: summary__heating_current_drive__nbi___current
    var"beam_current_fraction" :: summary__heating_current_drive__nbi___beam_current_fraction
    var"power_launched" :: summary__heating_current_drive__nbi___power_launched
    var"position" :: summary__heating_current_drive__nbi___position
    var"tangency_radius" :: summary__heating_current_drive__nbi___tangency_radius
    var"power" :: summary__heating_current_drive__nbi___power
    var"direction" :: summary__heating_current_drive__nbi___direction
    var"beam_power_fraction" :: summary__heating_current_drive__nbi___beam_power_fraction
    var"species" :: summary__heating_current_drive__nbi___species
    _parent :: WeakRef
    function summary__heating_current_drive__nbi(var"energy"=summary__heating_current_drive__nbi___energy(), var"angle"=summary__heating_current_drive__nbi___angle(), var"current"=summary__heating_current_drive__nbi___current(), var"beam_current_fraction"=summary__heating_current_drive__nbi___beam_current_fraction(), var"power_launched"=summary__heating_current_drive__nbi___power_launched(), var"position"=summary__heating_current_drive__nbi___position(), var"tangency_radius"=summary__heating_current_drive__nbi___tangency_radius(), var"power"=summary__heating_current_drive__nbi___power(), var"direction"=summary__heating_current_drive__nbi___direction(), var"beam_power_fraction"=summary__heating_current_drive__nbi___beam_power_fraction(), var"species"=summary__heating_current_drive__nbi___species(), _parent=WeakRef(missing))
        ids = new(var"energy", var"angle", var"current", var"beam_current_fraction", var"power_launched", var"position", var"tangency_radius", var"power", var"direction", var"beam_power_fraction", var"species", _parent)
        assign_expressions(ids)
        setfield!(ids.energy, :_parent, WeakRef(ids))
        setfield!(ids.angle, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.beam_current_fraction, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.tangency_radius, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.direction, :_parent, WeakRef(ids))
        setfield!(ids.beam_power_fraction, :_parent, WeakRef(ids))
        setfield!(ids.species, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___power_launched <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___power <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___n_parallel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___n_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___frequency <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___energy_fast <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___energy_fast(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh___current <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__lh___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__lh <: IDSvectorElement
    var"power_launched" :: summary__heating_current_drive__lh___power_launched
    var"energy_fast" :: summary__heating_current_drive__lh___energy_fast
    var"frequency" :: summary__heating_current_drive__lh___frequency
    var"position" :: summary__heating_current_drive__lh___position
    var"power" :: summary__heating_current_drive__lh___power
    var"n_parallel" :: summary__heating_current_drive__lh___n_parallel
    var"current" :: summary__heating_current_drive__lh___current
    _parent :: WeakRef
    function summary__heating_current_drive__lh(var"power_launched"=summary__heating_current_drive__lh___power_launched(), var"energy_fast"=summary__heating_current_drive__lh___energy_fast(), var"frequency"=summary__heating_current_drive__lh___frequency(), var"position"=summary__heating_current_drive__lh___position(), var"power"=summary__heating_current_drive__lh___power(), var"n_parallel"=summary__heating_current_drive__lh___n_parallel(), var"current"=summary__heating_current_drive__lh___current(), _parent=WeakRef(missing))
        ids = new(var"power_launched", var"energy_fast", var"frequency", var"position", var"power", var"n_parallel", var"current", _parent)
        assign_expressions(ids)
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast, :_parent, WeakRef(ids))
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.n_parallel, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___power_launched <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___power <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___phase <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___phase(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___n_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___n_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___k_perpendicular <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___k_perpendicular(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___harmonic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___harmonic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___frequency <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___energy_fast <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___energy_fast(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___e_field_plus_minus_ratio <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___e_field_plus_minus_ratio(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic___current <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ic___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ic <: IDSvectorElement
    var"harmonic" :: summary__heating_current_drive__ic___harmonic
    var"e_field_plus_minus_ratio" :: summary__heating_current_drive__ic___e_field_plus_minus_ratio
    var"current" :: summary__heating_current_drive__ic___current
    var"n_tor" :: summary__heating_current_drive__ic___n_tor
    var"power_launched" :: summary__heating_current_drive__ic___power_launched
    var"energy_fast" :: summary__heating_current_drive__ic___energy_fast
    var"position" :: summary__heating_current_drive__ic___position
    var"power" :: summary__heating_current_drive__ic___power
    var"phase" :: summary__heating_current_drive__ic___phase
    var"frequency" :: summary__heating_current_drive__ic___frequency
    var"k_perpendicular" :: summary__heating_current_drive__ic___k_perpendicular
    _parent :: WeakRef
    function summary__heating_current_drive__ic(var"harmonic"=summary__heating_current_drive__ic___harmonic(), var"e_field_plus_minus_ratio"=summary__heating_current_drive__ic___e_field_plus_minus_ratio(), var"current"=summary__heating_current_drive__ic___current(), var"n_tor"=summary__heating_current_drive__ic___n_tor(), var"power_launched"=summary__heating_current_drive__ic___power_launched(), var"energy_fast"=summary__heating_current_drive__ic___energy_fast(), var"position"=summary__heating_current_drive__ic___position(), var"power"=summary__heating_current_drive__ic___power(), var"phase"=summary__heating_current_drive__ic___phase(), var"frequency"=summary__heating_current_drive__ic___frequency(), var"k_perpendicular"=summary__heating_current_drive__ic___k_perpendicular(), _parent=WeakRef(missing))
        ids = new(var"harmonic", var"e_field_plus_minus_ratio", var"current", var"n_tor", var"power_launched", var"energy_fast", var"position", var"power", var"phase", var"frequency", var"k_perpendicular", _parent)
        assign_expressions(ids)
        setfield!(ids.harmonic, :_parent, WeakRef(ids))
        setfield!(ids.e_field_plus_minus_ratio, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.n_tor, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.phase, :_parent, WeakRef(ids))
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.k_perpendicular, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___power_launched <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___power_launched(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___power <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___position <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___position(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___polarisation <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___polarisation(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___harmonic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___harmonic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___frequency <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___energy_fast <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___energy_fast(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___current <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___angle_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___angle_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec___angle_pol <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__heating_current_drive__ec___angle_pol(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__heating_current_drive__ec <: IDSvectorElement
    var"angle_tor" :: summary__heating_current_drive__ec___angle_tor
    var"polarisation" :: summary__heating_current_drive__ec___polarisation
    var"power_launched" :: summary__heating_current_drive__ec___power_launched
    var"harmonic" :: summary__heating_current_drive__ec___harmonic
    var"energy_fast" :: summary__heating_current_drive__ec___energy_fast
    var"frequency" :: summary__heating_current_drive__ec___frequency
    var"position" :: summary__heating_current_drive__ec___position
    var"power" :: summary__heating_current_drive__ec___power
    var"angle_pol" :: summary__heating_current_drive__ec___angle_pol
    var"current" :: summary__heating_current_drive__ec___current
    _parent :: WeakRef
    function summary__heating_current_drive__ec(var"angle_tor"=summary__heating_current_drive__ec___angle_tor(), var"polarisation"=summary__heating_current_drive__ec___polarisation(), var"power_launched"=summary__heating_current_drive__ec___power_launched(), var"harmonic"=summary__heating_current_drive__ec___harmonic(), var"energy_fast"=summary__heating_current_drive__ec___energy_fast(), var"frequency"=summary__heating_current_drive__ec___frequency(), var"position"=summary__heating_current_drive__ec___position(), var"power"=summary__heating_current_drive__ec___power(), var"angle_pol"=summary__heating_current_drive__ec___angle_pol(), var"current"=summary__heating_current_drive__ec___current(), _parent=WeakRef(missing))
        ids = new(var"angle_tor", var"polarisation", var"power_launched", var"harmonic", var"energy_fast", var"frequency", var"position", var"power", var"angle_pol", var"current", _parent)
        assign_expressions(ids)
        setfield!(ids.angle_tor, :_parent, WeakRef(ids))
        setfield!(ids.polarisation, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
        setfield!(ids.harmonic, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast, :_parent, WeakRef(ids))
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.angle_pol, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive <: IDS
    var"power_launched_ec" :: summary__heating_current_drive__power_launched_ec
    var"ec" :: IDSvector{T} where {T<:summary__heating_current_drive__ec}
    var"power_nbi" :: summary__heating_current_drive__power_nbi
    var"power_ec" :: summary__heating_current_drive__power_ec
    var"power_lh" :: summary__heating_current_drive__power_lh
    var"ic" :: IDSvector{T} where {T<:summary__heating_current_drive__ic}
    var"lh" :: IDSvector{T} where {T<:summary__heating_current_drive__lh}
    var"power_ic" :: summary__heating_current_drive__power_ic
    var"power_launched_lh" :: summary__heating_current_drive__power_launched_lh
    var"nbi" :: IDSvector{T} where {T<:summary__heating_current_drive__nbi}
    var"power_launched_nbi_co_injected_ratio" :: summary__heating_current_drive__power_launched_nbi_co_injected_ratio
    var"power_additional" :: summary__heating_current_drive__power_additional
    var"power_launched_ic" :: summary__heating_current_drive__power_launched_ic
    var"power_launched_nbi" :: summary__heating_current_drive__power_launched_nbi
    _parent :: WeakRef
    function summary__heating_current_drive(var"power_launched_ec"=summary__heating_current_drive__power_launched_ec(), var"ec"=IDSvector(summary__heating_current_drive__ec[]), var"power_nbi"=summary__heating_current_drive__power_nbi(), var"power_ec"=summary__heating_current_drive__power_ec(), var"power_lh"=summary__heating_current_drive__power_lh(), var"ic"=IDSvector(summary__heating_current_drive__ic[]), var"lh"=IDSvector(summary__heating_current_drive__lh[]), var"power_ic"=summary__heating_current_drive__power_ic(), var"power_launched_lh"=summary__heating_current_drive__power_launched_lh(), var"nbi"=IDSvector(summary__heating_current_drive__nbi[]), var"power_launched_nbi_co_injected_ratio"=summary__heating_current_drive__power_launched_nbi_co_injected_ratio(), var"power_additional"=summary__heating_current_drive__power_additional(), var"power_launched_ic"=summary__heating_current_drive__power_launched_ic(), var"power_launched_nbi"=summary__heating_current_drive__power_launched_nbi(), _parent=WeakRef(missing))
        ids = new(var"power_launched_ec", var"ec", var"power_nbi", var"power_ec", var"power_lh", var"ic", var"lh", var"power_ic", var"power_launched_lh", var"nbi", var"power_launched_nbi_co_injected_ratio", var"power_additional", var"power_launched_ic", var"power_launched_nbi", _parent)
        assign_expressions(ids)
        setfield!(ids.power_launched_ec, :_parent, WeakRef(ids))
        setfield!(ids.ec, :_parent, WeakRef(ids))
        setfield!(ids.power_nbi, :_parent, WeakRef(ids))
        setfield!(ids.power_ec, :_parent, WeakRef(ids))
        setfield!(ids.power_lh, :_parent, WeakRef(ids))
        setfield!(ids.ic, :_parent, WeakRef(ids))
        setfield!(ids.lh, :_parent, WeakRef(ids))
        setfield!(ids.power_ic, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_lh, :_parent, WeakRef(ids))
        setfield!(ids.nbi, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_nbi_co_injected_ratio, :_parent, WeakRef(ids))
        setfield!(ids.power_additional, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_ic, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_nbi, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__global_quantities__volume <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__volume(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__v_loop <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__v_loop(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__tau_resistive <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__tau_resistive(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__tau_helium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__tau_helium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__tau_energy_98 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__tau_energy_98(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__tau_energy <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__tau_energy(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__resistance <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__resistance(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__ratio_tau_helium_fuel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__ratio_tau_helium_fuel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__r0 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__r0(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__q_95 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__q_95(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_synchrotron <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_synchrotron(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_steady <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_steady(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_radiated_outside_lcfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_radiated_outside_lcfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_radiated_inside_lcfs <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_radiated_inside_lcfs(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_radiated <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_radiated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_ohm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_ohm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_loss <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_loss(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_line <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_line(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__power_bremsstrahlung <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__power_bremsstrahlung(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__li_mhd <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__li_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__li <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__li(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__ip <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__ip(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__h_mode <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__global_quantities__h_mode(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__h_98 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__h_98(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__greenwald_fraction <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__greenwald_fraction(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__fusion_gain <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__fusion_gain(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__fusion_fluence <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__fusion_fluence(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_mhd <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_ion_total_thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_ion_total_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_fast_perpendicular <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_fast_perpendicular(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_fast_parallel <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_fast_parallel(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_electrons_thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_electrons_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_diamagnetic <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_diamagnetic(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__energy_b_field_pol <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__energy_b_field_pol(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__denergy_thermal_dt <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__denergy_thermal_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__denergy_diamagnetic_dt <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__denergy_diamagnetic_dt(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__current_ohm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__current_ohm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__current_non_inductive <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__current_non_inductive(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__current_bootstrap <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__current_bootstrap(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__current_alignment <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__current_alignment(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_tor_thermal_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_tor_thermal_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_tor_norm_mhd <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_tor_norm_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_tor_norm <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_tor_norm(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_tor_mhd <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_tor_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_tor <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_tor(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_pol_mhd <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_pol_mhd(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__beta_pol <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__beta_pol(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities__b0 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__global_quantities__b0(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__global_quantities <: IDS
    var"power_radiated_outside_lcfs" :: summary__global_quantities__power_radiated_outside_lcfs
    var"beta_tor" :: summary__global_quantities__beta_tor
    var"power_ohm" :: summary__global_quantities__power_ohm
    var"energy_b_field_pol" :: summary__global_quantities__energy_b_field_pol
    var"tau_energy" :: summary__global_quantities__tau_energy
    var"energy_electrons_thermal" :: summary__global_quantities__energy_electrons_thermal
    var"greenwald_fraction" :: summary__global_quantities__greenwald_fraction
    var"power_steady" :: summary__global_quantities__power_steady
    var"resistance" :: summary__global_quantities__resistance
    var"beta_pol" :: summary__global_quantities__beta_pol
    var"energy_thermal" :: summary__global_quantities__energy_thermal
    var"fusion_gain" :: summary__global_quantities__fusion_gain
    var"power_radiated_inside_lcfs" :: summary__global_quantities__power_radiated_inside_lcfs
    var"beta_tor_norm" :: summary__global_quantities__beta_tor_norm
    var"beta_tor_thermal_norm" :: summary__global_quantities__beta_tor_thermal_norm
    var"tau_helium" :: summary__global_quantities__tau_helium
    var"power_loss" :: summary__global_quantities__power_loss
    var"current_alignment" :: summary__global_quantities__current_alignment
    var"energy_diamagnetic" :: summary__global_quantities__energy_diamagnetic
    var"tau_energy_98" :: summary__global_quantities__tau_energy_98
    var"li_mhd" :: summary__global_quantities__li_mhd
    var"energy_fast_perpendicular" :: summary__global_quantities__energy_fast_perpendicular
    var"q_95" :: summary__global_quantities__q_95
    var"beta_tor_mhd" :: summary__global_quantities__beta_tor_mhd
    var"volume" :: summary__global_quantities__volume
    var"b0" :: summary__global_quantities__b0
    var"ip" :: summary__global_quantities__ip
    var"beta_pol_mhd" :: summary__global_quantities__beta_pol_mhd
    var"denergy_diamagnetic_dt" :: summary__global_quantities__denergy_diamagnetic_dt
    var"h_mode" :: summary__global_quantities__h_mode
    var"power_bremsstrahlung" :: summary__global_quantities__power_bremsstrahlung
    var"h_98" :: summary__global_quantities__h_98
    var"tau_resistive" :: summary__global_quantities__tau_resistive
    var"current_non_inductive" :: summary__global_quantities__current_non_inductive
    var"current_ohm" :: summary__global_quantities__current_ohm
    var"energy_total" :: summary__global_quantities__energy_total
    var"power_synchrotron" :: summary__global_quantities__power_synchrotron
    var"v_loop" :: summary__global_quantities__v_loop
    var"energy_ion_total_thermal" :: summary__global_quantities__energy_ion_total_thermal
    var"power_line" :: summary__global_quantities__power_line
    var"current_bootstrap" :: summary__global_quantities__current_bootstrap
    var"fusion_fluence" :: summary__global_quantities__fusion_fluence
    var"energy_fast_parallel" :: summary__global_quantities__energy_fast_parallel
    var"power_radiated" :: summary__global_quantities__power_radiated
    var"denergy_thermal_dt" :: summary__global_quantities__denergy_thermal_dt
    var"energy_mhd" :: summary__global_quantities__energy_mhd
    var"li" :: summary__global_quantities__li
    var"ratio_tau_helium_fuel" :: summary__global_quantities__ratio_tau_helium_fuel
    var"r0" :: summary__global_quantities__r0
    var"beta_tor_norm_mhd" :: summary__global_quantities__beta_tor_norm_mhd
    _parent :: WeakRef
    function summary__global_quantities(var"power_radiated_outside_lcfs"=summary__global_quantities__power_radiated_outside_lcfs(), var"beta_tor"=summary__global_quantities__beta_tor(), var"power_ohm"=summary__global_quantities__power_ohm(), var"energy_b_field_pol"=summary__global_quantities__energy_b_field_pol(), var"tau_energy"=summary__global_quantities__tau_energy(), var"energy_electrons_thermal"=summary__global_quantities__energy_electrons_thermal(), var"greenwald_fraction"=summary__global_quantities__greenwald_fraction(), var"power_steady"=summary__global_quantities__power_steady(), var"resistance"=summary__global_quantities__resistance(), var"beta_pol"=summary__global_quantities__beta_pol(), var"energy_thermal"=summary__global_quantities__energy_thermal(), var"fusion_gain"=summary__global_quantities__fusion_gain(), var"power_radiated_inside_lcfs"=summary__global_quantities__power_radiated_inside_lcfs(), var"beta_tor_norm"=summary__global_quantities__beta_tor_norm(), var"beta_tor_thermal_norm"=summary__global_quantities__beta_tor_thermal_norm(), var"tau_helium"=summary__global_quantities__tau_helium(), var"power_loss"=summary__global_quantities__power_loss(), var"current_alignment"=summary__global_quantities__current_alignment(), var"energy_diamagnetic"=summary__global_quantities__energy_diamagnetic(), var"tau_energy_98"=summary__global_quantities__tau_energy_98(), var"li_mhd"=summary__global_quantities__li_mhd(), var"energy_fast_perpendicular"=summary__global_quantities__energy_fast_perpendicular(), var"q_95"=summary__global_quantities__q_95(), var"beta_tor_mhd"=summary__global_quantities__beta_tor_mhd(), var"volume"=summary__global_quantities__volume(), var"b0"=summary__global_quantities__b0(), var"ip"=summary__global_quantities__ip(), var"beta_pol_mhd"=summary__global_quantities__beta_pol_mhd(), var"denergy_diamagnetic_dt"=summary__global_quantities__denergy_diamagnetic_dt(), var"h_mode"=summary__global_quantities__h_mode(), var"power_bremsstrahlung"=summary__global_quantities__power_bremsstrahlung(), var"h_98"=summary__global_quantities__h_98(), var"tau_resistive"=summary__global_quantities__tau_resistive(), var"current_non_inductive"=summary__global_quantities__current_non_inductive(), var"current_ohm"=summary__global_quantities__current_ohm(), var"energy_total"=summary__global_quantities__energy_total(), var"power_synchrotron"=summary__global_quantities__power_synchrotron(), var"v_loop"=summary__global_quantities__v_loop(), var"energy_ion_total_thermal"=summary__global_quantities__energy_ion_total_thermal(), var"power_line"=summary__global_quantities__power_line(), var"current_bootstrap"=summary__global_quantities__current_bootstrap(), var"fusion_fluence"=summary__global_quantities__fusion_fluence(), var"energy_fast_parallel"=summary__global_quantities__energy_fast_parallel(), var"power_radiated"=summary__global_quantities__power_radiated(), var"denergy_thermal_dt"=summary__global_quantities__denergy_thermal_dt(), var"energy_mhd"=summary__global_quantities__energy_mhd(), var"li"=summary__global_quantities__li(), var"ratio_tau_helium_fuel"=summary__global_quantities__ratio_tau_helium_fuel(), var"r0"=summary__global_quantities__r0(), var"beta_tor_norm_mhd"=summary__global_quantities__beta_tor_norm_mhd(), _parent=WeakRef(missing))
        ids = new(var"power_radiated_outside_lcfs", var"beta_tor", var"power_ohm", var"energy_b_field_pol", var"tau_energy", var"energy_electrons_thermal", var"greenwald_fraction", var"power_steady", var"resistance", var"beta_pol", var"energy_thermal", var"fusion_gain", var"power_radiated_inside_lcfs", var"beta_tor_norm", var"beta_tor_thermal_norm", var"tau_helium", var"power_loss", var"current_alignment", var"energy_diamagnetic", var"tau_energy_98", var"li_mhd", var"energy_fast_perpendicular", var"q_95", var"beta_tor_mhd", var"volume", var"b0", var"ip", var"beta_pol_mhd", var"denergy_diamagnetic_dt", var"h_mode", var"power_bremsstrahlung", var"h_98", var"tau_resistive", var"current_non_inductive", var"current_ohm", var"energy_total", var"power_synchrotron", var"v_loop", var"energy_ion_total_thermal", var"power_line", var"current_bootstrap", var"fusion_fluence", var"energy_fast_parallel", var"power_radiated", var"denergy_thermal_dt", var"energy_mhd", var"li", var"ratio_tau_helium_fuel", var"r0", var"beta_tor_norm_mhd", _parent)
        assign_expressions(ids)
        setfield!(ids.power_radiated_outside_lcfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor, :_parent, WeakRef(ids))
        setfield!(ids.power_ohm, :_parent, WeakRef(ids))
        setfield!(ids.energy_b_field_pol, :_parent, WeakRef(ids))
        setfield!(ids.tau_energy, :_parent, WeakRef(ids))
        setfield!(ids.energy_electrons_thermal, :_parent, WeakRef(ids))
        setfield!(ids.greenwald_fraction, :_parent, WeakRef(ids))
        setfield!(ids.power_steady, :_parent, WeakRef(ids))
        setfield!(ids.resistance, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal, :_parent, WeakRef(ids))
        setfield!(ids.fusion_gain, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated_inside_lcfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_norm, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_thermal_norm, :_parent, WeakRef(ids))
        setfield!(ids.tau_helium, :_parent, WeakRef(ids))
        setfield!(ids.power_loss, :_parent, WeakRef(ids))
        setfield!(ids.current_alignment, :_parent, WeakRef(ids))
        setfield!(ids.energy_diamagnetic, :_parent, WeakRef(ids))
        setfield!(ids.tau_energy_98, :_parent, WeakRef(ids))
        setfield!(ids.li_mhd, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast_perpendicular, :_parent, WeakRef(ids))
        setfield!(ids.q_95, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_mhd, :_parent, WeakRef(ids))
        setfield!(ids.volume, :_parent, WeakRef(ids))
        setfield!(ids.b0, :_parent, WeakRef(ids))
        setfield!(ids.ip, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_mhd, :_parent, WeakRef(ids))
        setfield!(ids.denergy_diamagnetic_dt, :_parent, WeakRef(ids))
        setfield!(ids.h_mode, :_parent, WeakRef(ids))
        setfield!(ids.power_bremsstrahlung, :_parent, WeakRef(ids))
        setfield!(ids.h_98, :_parent, WeakRef(ids))
        setfield!(ids.tau_resistive, :_parent, WeakRef(ids))
        setfield!(ids.current_non_inductive, :_parent, WeakRef(ids))
        setfield!(ids.current_ohm, :_parent, WeakRef(ids))
        setfield!(ids.energy_total, :_parent, WeakRef(ids))
        setfield!(ids.power_synchrotron, :_parent, WeakRef(ids))
        setfield!(ids.v_loop, :_parent, WeakRef(ids))
        setfield!(ids.energy_ion_total_thermal, :_parent, WeakRef(ids))
        setfield!(ids.power_line, :_parent, WeakRef(ids))
        setfield!(ids.current_bootstrap, :_parent, WeakRef(ids))
        setfield!(ids.fusion_fluence, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast_parallel, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated, :_parent, WeakRef(ids))
        setfield!(ids.denergy_thermal_dt, :_parent, WeakRef(ids))
        setfield!(ids.energy_mhd, :_parent, WeakRef(ids))
        setfield!(ids.li, :_parent, WeakRef(ids))
        setfield!(ids.ratio_tau_helium_fuel, :_parent, WeakRef(ids))
        setfield!(ids.r0, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_norm_mhd, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__gas_injection_rates__xenon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__xenon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__tritium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__tritium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__top <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__top(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__silane <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__silane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__propane <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__propane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__oxygen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__oxygen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__nitrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__nitrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__neon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__neon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__midplane <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__midplane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__methane_deuterated <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__methane_deuterated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__methane_carbon_13 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__methane_carbon_13(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__methane <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__methane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__lithium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__lithium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__krypton <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__krypton(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__impurity_seeding <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__impurity_seeding(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__hydrogen <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__hydrogen(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__helium_4 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__helium_4(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__helium_3 <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__helium_3(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__ethylene <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__ethylene(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__ethane <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__ethane(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__deuterium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__deuterium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__carbon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__carbon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__bottom <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__bottom(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__beryllium <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__beryllium(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__argon <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__argon(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__ammonia_deuterated <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__ammonia_deuterated(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates__ammonia <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__gas_injection_rates__ammonia(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__gas_injection_rates <: IDS
    var"krypton" :: summary__gas_injection_rates__krypton
    var"lithium" :: summary__gas_injection_rates__lithium
    var"methane_carbon_13" :: summary__gas_injection_rates__methane_carbon_13
    var"neon" :: summary__gas_injection_rates__neon
    var"ethylene" :: summary__gas_injection_rates__ethylene
    var"tritium" :: summary__gas_injection_rates__tritium
    var"helium_3" :: summary__gas_injection_rates__helium_3
    var"deuterium" :: summary__gas_injection_rates__deuterium
    var"bottom" :: summary__gas_injection_rates__bottom
    var"helium_4" :: summary__gas_injection_rates__helium_4
    var"oxygen" :: summary__gas_injection_rates__oxygen
    var"xenon" :: summary__gas_injection_rates__xenon
    var"methane" :: summary__gas_injection_rates__methane
    var"hydrogen" :: summary__gas_injection_rates__hydrogen
    var"carbon" :: summary__gas_injection_rates__carbon
    var"ammonia_deuterated" :: summary__gas_injection_rates__ammonia_deuterated
    var"total" :: summary__gas_injection_rates__total
    var"top" :: summary__gas_injection_rates__top
    var"midplane" :: summary__gas_injection_rates__midplane
    var"ethane" :: summary__gas_injection_rates__ethane
    var"silane" :: summary__gas_injection_rates__silane
    var"ammonia" :: summary__gas_injection_rates__ammonia
    var"methane_deuterated" :: summary__gas_injection_rates__methane_deuterated
    var"nitrogen" :: summary__gas_injection_rates__nitrogen
    var"propane" :: summary__gas_injection_rates__propane
    var"beryllium" :: summary__gas_injection_rates__beryllium
    var"argon" :: summary__gas_injection_rates__argon
    var"impurity_seeding" :: summary__gas_injection_rates__impurity_seeding
    _parent :: WeakRef
    function summary__gas_injection_rates(var"krypton"=summary__gas_injection_rates__krypton(), var"lithium"=summary__gas_injection_rates__lithium(), var"methane_carbon_13"=summary__gas_injection_rates__methane_carbon_13(), var"neon"=summary__gas_injection_rates__neon(), var"ethylene"=summary__gas_injection_rates__ethylene(), var"tritium"=summary__gas_injection_rates__tritium(), var"helium_3"=summary__gas_injection_rates__helium_3(), var"deuterium"=summary__gas_injection_rates__deuterium(), var"bottom"=summary__gas_injection_rates__bottom(), var"helium_4"=summary__gas_injection_rates__helium_4(), var"oxygen"=summary__gas_injection_rates__oxygen(), var"xenon"=summary__gas_injection_rates__xenon(), var"methane"=summary__gas_injection_rates__methane(), var"hydrogen"=summary__gas_injection_rates__hydrogen(), var"carbon"=summary__gas_injection_rates__carbon(), var"ammonia_deuterated"=summary__gas_injection_rates__ammonia_deuterated(), var"total"=summary__gas_injection_rates__total(), var"top"=summary__gas_injection_rates__top(), var"midplane"=summary__gas_injection_rates__midplane(), var"ethane"=summary__gas_injection_rates__ethane(), var"silane"=summary__gas_injection_rates__silane(), var"ammonia"=summary__gas_injection_rates__ammonia(), var"methane_deuterated"=summary__gas_injection_rates__methane_deuterated(), var"nitrogen"=summary__gas_injection_rates__nitrogen(), var"propane"=summary__gas_injection_rates__propane(), var"beryllium"=summary__gas_injection_rates__beryllium(), var"argon"=summary__gas_injection_rates__argon(), var"impurity_seeding"=summary__gas_injection_rates__impurity_seeding(), _parent=WeakRef(missing))
        ids = new(var"krypton", var"lithium", var"methane_carbon_13", var"neon", var"ethylene", var"tritium", var"helium_3", var"deuterium", var"bottom", var"helium_4", var"oxygen", var"xenon", var"methane", var"hydrogen", var"carbon", var"ammonia_deuterated", var"total", var"top", var"midplane", var"ethane", var"silane", var"ammonia", var"methane_deuterated", var"nitrogen", var"propane", var"beryllium", var"argon", var"impurity_seeding", _parent)
        assign_expressions(ids)
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.methane_carbon_13, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.ethylene, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.bottom, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
        setfield!(ids.methane, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.ammonia_deuterated, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.top, :_parent, WeakRef(ids))
        setfield!(ids.midplane, :_parent, WeakRef(ids))
        setfield!(ids.ethane, :_parent, WeakRef(ids))
        setfield!(ids.silane, :_parent, WeakRef(ids))
        setfield!(ids.ammonia, :_parent, WeakRef(ids))
        setfield!(ids.methane_deuterated, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.propane, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.impurity_seeding, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__fusion__power <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__power(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_power_total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_power_total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__tt__total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__tt__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__tt__thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__tt__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__tt__beam_thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__tt__beam_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__tt__beam_beam <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__tt__beam_beam(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__tt <: IDS
    var"beam_thermal" :: summary__fusion__neutron_fluxes__tt__beam_thermal
    var"beam_beam" :: summary__fusion__neutron_fluxes__tt__beam_beam
    var"total" :: summary__fusion__neutron_fluxes__tt__total
    var"thermal" :: summary__fusion__neutron_fluxes__tt__thermal
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__tt(var"beam_thermal"=summary__fusion__neutron_fluxes__tt__beam_thermal(), var"beam_beam"=summary__fusion__neutron_fluxes__tt__beam_beam(), var"total"=summary__fusion__neutron_fluxes__tt__total(), var"thermal"=summary__fusion__neutron_fluxes__tt__thermal(), _parent=WeakRef(missing))
        ids = new(var"beam_thermal", var"beam_beam", var"total", var"thermal", _parent)
        assign_expressions(ids)
        setfield!(ids.beam_thermal, :_parent, WeakRef(ids))
        setfield!(ids.beam_beam, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dt__total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dt__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dt__thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dt__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dt__beam_thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dt__beam_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dt__beam_beam <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dt__beam_beam(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dt <: IDS
    var"beam_thermal" :: summary__fusion__neutron_fluxes__dt__beam_thermal
    var"beam_beam" :: summary__fusion__neutron_fluxes__dt__beam_beam
    var"total" :: summary__fusion__neutron_fluxes__dt__total
    var"thermal" :: summary__fusion__neutron_fluxes__dt__thermal
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dt(var"beam_thermal"=summary__fusion__neutron_fluxes__dt__beam_thermal(), var"beam_beam"=summary__fusion__neutron_fluxes__dt__beam_beam(), var"total"=summary__fusion__neutron_fluxes__dt__total(), var"thermal"=summary__fusion__neutron_fluxes__dt__thermal(), _parent=WeakRef(missing))
        ids = new(var"beam_thermal", var"beam_beam", var"total", var"thermal", _parent)
        assign_expressions(ids)
        setfield!(ids.beam_thermal, :_parent, WeakRef(ids))
        setfield!(ids.beam_beam, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dd__total <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dd__total(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dd__thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dd__thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dd__beam_thermal <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dd__beam_thermal(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dd__beam_beam <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dd__beam_beam(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes__dd <: IDS
    var"beam_thermal" :: summary__fusion__neutron_fluxes__dd__beam_thermal
    var"beam_beam" :: summary__fusion__neutron_fluxes__dd__beam_beam
    var"total" :: summary__fusion__neutron_fluxes__dd__total
    var"thermal" :: summary__fusion__neutron_fluxes__dd__thermal
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dd(var"beam_thermal"=summary__fusion__neutron_fluxes__dd__beam_thermal(), var"beam_beam"=summary__fusion__neutron_fluxes__dd__beam_beam(), var"total"=summary__fusion__neutron_fluxes__dd__total(), var"thermal"=summary__fusion__neutron_fluxes__dd__thermal(), _parent=WeakRef(missing))
        ids = new(var"beam_thermal", var"beam_beam", var"total", var"thermal", _parent)
        assign_expressions(ids)
        setfield!(ids.beam_thermal, :_parent, WeakRef(ids))
        setfield!(ids.beam_beam, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes <: IDS
    var"tt" :: summary__fusion__neutron_fluxes__tt
    var"total" :: summary__fusion__neutron_fluxes__total
    var"thermal" :: summary__fusion__neutron_fluxes__thermal
    var"dt" :: summary__fusion__neutron_fluxes__dt
    var"dd" :: summary__fusion__neutron_fluxes__dd
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes(var"tt"=summary__fusion__neutron_fluxes__tt(), var"total"=summary__fusion__neutron_fluxes__total(), var"thermal"=summary__fusion__neutron_fluxes__thermal(), var"dt"=summary__fusion__neutron_fluxes__dt(), var"dd"=summary__fusion__neutron_fluxes__dd(), _parent=WeakRef(missing))
        ids = new(var"tt", var"total", var"thermal", var"dt", var"dd", _parent)
        assign_expressions(ids)
        setfield!(ids.tt, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        setfield!(ids.dt, :_parent, WeakRef(ids))
        setfield!(ids.dd, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__fusion__current <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__fusion__current(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__fusion <: IDS
    var"neutron_power_total" :: summary__fusion__neutron_power_total
    var"neutron_fluxes" :: summary__fusion__neutron_fluxes
    var"power" :: summary__fusion__power
    var"current" :: summary__fusion__current
    _parent :: WeakRef
    function summary__fusion(var"neutron_power_total"=summary__fusion__neutron_power_total(), var"neutron_fluxes"=summary__fusion__neutron_fluxes(), var"power"=summary__fusion__power(), var"current"=summary__fusion__current(), _parent=WeakRef(missing))
        ids = new(var"neutron_power_total", var"neutron_fluxes", var"power", var"current", _parent)
        assign_expressions(ids)
        setfield!(ids.neutron_power_total, :_parent, WeakRef(ids))
        setfield!(ids.neutron_fluxes, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__elms__type <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__elms__type(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__elms__frequency <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__elms__frequency(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__elms <: IDS
    var"frequency" :: summary__elms__frequency
    var"type" :: summary__elms__type
    _parent :: WeakRef
    function summary__elms(var"frequency"=summary__elms__frequency(), var"type"=summary__elms__type(), _parent=WeakRef(missing))
        ids = new(var"frequency", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__disruption__vertical_displacement <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__disruption__vertical_displacement(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__disruption__time_radiated_power_max <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__disruption__time_radiated_power_max(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__disruption__time_half_ip <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__disruption__time_half_ip(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__disruption__time <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function summary__disruption__time(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__disruption__mitigation_valve <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function summary__disruption__mitigation_valve(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__disruption <: IDS
    var"time" :: summary__disruption__time
    var"mitigation_valve" :: summary__disruption__mitigation_valve
    var"vertical_displacement" :: summary__disruption__vertical_displacement
    var"time_half_ip" :: summary__disruption__time_half_ip
    var"time_radiated_power_max" :: summary__disruption__time_radiated_power_max
    _parent :: WeakRef
    function summary__disruption(var"time"=summary__disruption__time(), var"mitigation_valve"=summary__disruption__mitigation_valve(), var"vertical_displacement"=summary__disruption__vertical_displacement(), var"time_half_ip"=summary__disruption__time_half_ip(), var"time_radiated_power_max"=summary__disruption__time_radiated_power_max(), _parent=WeakRef(missing))
        ids = new(var"time", var"mitigation_valve", var"vertical_displacement", var"time_half_ip", var"time_radiated_power_max", _parent)
        assign_expressions(ids)
        setfield!(ids.time, :_parent, WeakRef(ids))
        setfield!(ids.mitigation_valve, :_parent, WeakRef(ids))
        setfield!(ids.vertical_displacement, :_parent, WeakRef(ids))
        setfield!(ids.time_half_ip, :_parent, WeakRef(ids))
        setfield!(ids.time_radiated_power_max, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__configuration <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__configuration(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__code <: IDS
    var"library" :: IDSvector{T} where {T<:summary__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__code(var"library"=IDSvector(summary__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__boundary__type <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function summary__boundary__type(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__triangularity_upper <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__triangularity_upper(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__triangularity_lower <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__triangularity_lower(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__strike_point_outer_z <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__strike_point_outer_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__strike_point_outer_r <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__strike_point_outer_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__strike_point_inner_z <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__strike_point_inner_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__strike_point_inner_r <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__strike_point_inner_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__strike_point_configuration <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__boundary__strike_point_configuration(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__minor_radius <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__minor_radius(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__magnetic_axis_z <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__magnetic_axis_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__magnetic_axis_r <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__magnetic_axis_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__geometric_axis_z <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__geometric_axis_z(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__geometric_axis_r <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__geometric_axis_r(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__gap_limiter_wall <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__gap_limiter_wall(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary__elongation <: IDS
    var"source" :: Union{Missing, String, Function}
    var"value" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__boundary__elongation(var"source"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"source", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct summary__boundary <: IDS
    var"strike_point_inner_r" :: summary__boundary__strike_point_inner_r
    var"strike_point_outer_z" :: summary__boundary__strike_point_outer_z
    var"geometric_axis_z" :: summary__boundary__geometric_axis_z
    var"magnetic_axis_z" :: summary__boundary__magnetic_axis_z
    var"geometric_axis_r" :: summary__boundary__geometric_axis_r
    var"strike_point_configuration" :: summary__boundary__strike_point_configuration
    var"triangularity_upper" :: summary__boundary__triangularity_upper
    var"gap_limiter_wall" :: summary__boundary__gap_limiter_wall
    var"strike_point_inner_z" :: summary__boundary__strike_point_inner_z
    var"triangularity_lower" :: summary__boundary__triangularity_lower
    var"minor_radius" :: summary__boundary__minor_radius
    var"strike_point_outer_r" :: summary__boundary__strike_point_outer_r
    var"elongation" :: summary__boundary__elongation
    var"type" :: summary__boundary__type
    var"magnetic_axis_r" :: summary__boundary__magnetic_axis_r
    _parent :: WeakRef
    function summary__boundary(var"strike_point_inner_r"=summary__boundary__strike_point_inner_r(), var"strike_point_outer_z"=summary__boundary__strike_point_outer_z(), var"geometric_axis_z"=summary__boundary__geometric_axis_z(), var"magnetic_axis_z"=summary__boundary__magnetic_axis_z(), var"geometric_axis_r"=summary__boundary__geometric_axis_r(), var"strike_point_configuration"=summary__boundary__strike_point_configuration(), var"triangularity_upper"=summary__boundary__triangularity_upper(), var"gap_limiter_wall"=summary__boundary__gap_limiter_wall(), var"strike_point_inner_z"=summary__boundary__strike_point_inner_z(), var"triangularity_lower"=summary__boundary__triangularity_lower(), var"minor_radius"=summary__boundary__minor_radius(), var"strike_point_outer_r"=summary__boundary__strike_point_outer_r(), var"elongation"=summary__boundary__elongation(), var"type"=summary__boundary__type(), var"magnetic_axis_r"=summary__boundary__magnetic_axis_r(), _parent=WeakRef(missing))
        ids = new(var"strike_point_inner_r", var"strike_point_outer_z", var"geometric_axis_z", var"magnetic_axis_z", var"geometric_axis_r", var"strike_point_configuration", var"triangularity_upper", var"gap_limiter_wall", var"strike_point_inner_z", var"triangularity_lower", var"minor_radius", var"strike_point_outer_r", var"elongation", var"type", var"magnetic_axis_r", _parent)
        assign_expressions(ids)
        setfield!(ids.strike_point_inner_r, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_outer_z, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis_z, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis_z, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis_r, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_configuration, :_parent, WeakRef(ids))
        setfield!(ids.triangularity_upper, :_parent, WeakRef(ids))
        setfield!(ids.gap_limiter_wall, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_inner_z, :_parent, WeakRef(ids))
        setfield!(ids.triangularity_lower, :_parent, WeakRef(ids))
        setfield!(ids.minor_radius, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_outer_r, :_parent, WeakRef(ids))
        setfield!(ids.elongation, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis_r, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary <: IDS
    var"local" :: summary__local
    var"wall" :: summary__wall
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"code" :: summary__code
    var"stationary_phase_flag" :: summary__stationary_phase_flag
    var"gas_injection_rates" :: summary__gas_injection_rates
    var"configuration" :: summary__configuration
    var"boundary" :: summary__boundary
    var"pedestal_fits" :: summary__pedestal_fits
    var"heating_current_drive" :: summary__heating_current_drive
    var"global_quantities" :: summary__global_quantities
    var"disruption" :: summary__disruption
    var"fusion" :: summary__fusion
    var"rmps" :: summary__rmps
    var"ids_properties" :: summary__ids_properties
    var"midplane" :: summary__midplane
    var"pellets" :: summary__pellets
    var"elms" :: summary__elms
    var"kicks" :: summary__kicks
    var"tag" :: summary__tag
    var"line_average" :: summary__line_average
    var"scrape_off_layer" :: summary__scrape_off_layer
    var"limiter" :: summary__limiter
    var"runaways" :: summary__runaways
    var"time_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"magnetic_shear_flag" :: summary__magnetic_shear_flag
    var"volume_average" :: summary__volume_average
    _parent :: WeakRef
    function summary(var"local"=summary__local(), var"wall"=summary__wall(), var"time"=missing, var"code"=summary__code(), var"stationary_phase_flag"=summary__stationary_phase_flag(), var"gas_injection_rates"=summary__gas_injection_rates(), var"configuration"=summary__configuration(), var"boundary"=summary__boundary(), var"pedestal_fits"=summary__pedestal_fits(), var"heating_current_drive"=summary__heating_current_drive(), var"global_quantities"=summary__global_quantities(), var"disruption"=summary__disruption(), var"fusion"=summary__fusion(), var"rmps"=summary__rmps(), var"ids_properties"=summary__ids_properties(), var"midplane"=summary__midplane(), var"pellets"=summary__pellets(), var"elms"=summary__elms(), var"kicks"=summary__kicks(), var"tag"=summary__tag(), var"line_average"=summary__line_average(), var"scrape_off_layer"=summary__scrape_off_layer(), var"limiter"=summary__limiter(), var"runaways"=summary__runaways(), var"time_width"=missing, var"magnetic_shear_flag"=summary__magnetic_shear_flag(), var"volume_average"=summary__volume_average(), _parent=WeakRef(missing))
        ids = new(var"local", var"wall", var"time", var"code", var"stationary_phase_flag", var"gas_injection_rates", var"configuration", var"boundary", var"pedestal_fits", var"heating_current_drive", var"global_quantities", var"disruption", var"fusion", var"rmps", var"ids_properties", var"midplane", var"pellets", var"elms", var"kicks", var"tag", var"line_average", var"scrape_off_layer", var"limiter", var"runaways", var"time_width", var"magnetic_shear_flag", var"volume_average", _parent)
        assign_expressions(ids)
        setfield!(ids.local, :_parent, WeakRef(ids))
        setfield!(ids.wall, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        setfield!(ids.stationary_phase_flag, :_parent, WeakRef(ids))
        setfield!(ids.gas_injection_rates, :_parent, WeakRef(ids))
        setfield!(ids.configuration, :_parent, WeakRef(ids))
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_fits, :_parent, WeakRef(ids))
        setfield!(ids.heating_current_drive, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.disruption, :_parent, WeakRef(ids))
        setfield!(ids.fusion, :_parent, WeakRef(ids))
        setfield!(ids.rmps, :_parent, WeakRef(ids))
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.midplane, :_parent, WeakRef(ids))
        setfield!(ids.pellets, :_parent, WeakRef(ids))
        setfield!(ids.elms, :_parent, WeakRef(ids))
        setfield!(ids.kicks, :_parent, WeakRef(ids))
        setfield!(ids.tag, :_parent, WeakRef(ids))
        setfield!(ids.line_average, :_parent, WeakRef(ids))
        setfield!(ids.scrape_off_layer, :_parent, WeakRef(ids))
        setfield!(ids.limiter, :_parent, WeakRef(ids))
        setfield!(ids.runaways, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear_flag, :_parent, WeakRef(ids))
        setfield!(ids.volume_average, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct radial_build__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function radial_build__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct radial_build__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: radial_build__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function radial_build__ids_properties(var"provider"=missing, var"version_put"=radial_build__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct radial_build__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function radial_build__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct radial_build__code <: IDS
    var"library" :: IDSvector{T} where {T<:radial_build__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function radial_build__code(var"library"=IDSvector(radial_build__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct radial_build__center_stack <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"material" :: Union{Missing, String, Function}
    var"start" :: Union{Missing, Real, Function}
    var"end" :: Union{Missing, Real, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function radial_build__center_stack(var"name"=missing, var"material"=missing, var"start"=missing, var"end"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"material", var"start", var"end", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct radial_build <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"center_stack" :: IDSvector{T} where {T<:radial_build__center_stack}
    var"ids_properties" :: radial_build__ids_properties
    var"code" :: radial_build__code
    _parent :: WeakRef
    function radial_build(var"time"=missing, var"center_stack"=IDSvector(radial_build__center_stack[]), var"ids_properties"=radial_build__ids_properties(), var"code"=radial_build__code(), _parent=WeakRef(missing))
        ids = new(var"time", var"center_stack", var"ids_properties", var"code", _parent)
        assign_expressions(ids)
        setfield!(ids.center_stack, :_parent, WeakRef(ids))
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__vertical_force___force <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__vertical_force___force(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__vertical_force <: IDSvectorElement
    var"limit_min" :: Union{Missing, Real, Function}
    var"combination" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"name" :: Union{Missing, String, Function}
    var"force" :: pf_active__vertical_force___force
    var"limit_max" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__vertical_force(var"limit_min"=missing, var"combination"=missing, var"name"=missing, var"force"=pf_active__vertical_force___force(), var"limit_max"=missing, _parent=WeakRef(missing))
        ids = new(var"limit_min", var"combination", var"name", var"force", var"limit_max", _parent)
        assign_expressions(ids)
        setfield!(ids.force, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__supply___voltage <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__supply___voltage(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__supply___current <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__supply___current(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__supply <: IDSvectorElement
    var"voltage_limit_min" :: Union{Missing, Real, Function}
    var"energy_limit_max" :: Union{Missing, Real, Function}
    var"filter_denominator" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current" :: pf_active__supply___current
    var"voltage_limit_max" :: Union{Missing, Real, Function}
    var"name" :: Union{Missing, String, Function}
    var"filter_numerator" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"nonlinear_model" :: Union{Missing, String, Function}
    var"current_limit_min" :: Union{Missing, Real, Function}
    var"delay" :: Union{Missing, Real, Function}
    var"resistance" :: Union{Missing, Real, Function}
    var"voltage" :: pf_active__supply___voltage
    var"current_limiter_gain" :: Union{Missing, Real, Function}
    var"identifier" :: Union{Missing, String, Function}
    var"current_limit_max" :: Union{Missing, Real, Function}
    var"type" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function pf_active__supply(var"voltage_limit_min"=missing, var"energy_limit_max"=missing, var"filter_denominator"=missing, var"current"=pf_active__supply___current(), var"voltage_limit_max"=missing, var"name"=missing, var"filter_numerator"=missing, var"nonlinear_model"=missing, var"current_limit_min"=missing, var"delay"=missing, var"resistance"=missing, var"voltage"=pf_active__supply___voltage(), var"current_limiter_gain"=missing, var"identifier"=missing, var"current_limit_max"=missing, var"type"=missing, _parent=WeakRef(missing))
        ids = new(var"voltage_limit_min", var"energy_limit_max", var"filter_denominator", var"current", var"voltage_limit_max", var"name", var"filter_numerator", var"nonlinear_model", var"current_limit_min", var"delay", var"resistance", var"voltage", var"current_limiter_gain", var"identifier", var"current_limit_max", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__radial_force___force <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__radial_force___force(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__radial_force <: IDSvectorElement
    var"limit_min" :: Union{Missing, Real, Function}
    var"combination" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"name" :: Union{Missing, String, Function}
    var"force" :: pf_active__radial_force___force
    var"limit_max" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__radial_force(var"limit_min"=missing, var"combination"=missing, var"name"=missing, var"force"=pf_active__radial_force___force(), var"limit_max"=missing, _parent=WeakRef(missing))
        ids = new(var"limit_min", var"combination", var"name", var"force", var"limit_max", _parent)
        assign_expressions(ids)
        setfield!(ids.force, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function pf_active__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: pf_active__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function pf_active__ids_properties(var"provider"=missing, var"version_put"=pf_active__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__global_quantities <: IDS
    var"psi_coils_list" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"psi_coils_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__global_quantities(var"psi_coils_list"=missing, var"psi_coils_average"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"psi_coils_list", var"psi_coils_average", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___voltage <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___voltage(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___element___geometry__rectangle <: IDS
    var"height" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"width" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__coil___element___geometry__rectangle(var"height"=missing, var"r"=missing, var"width"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"height", var"r", var"width", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___element___geometry__outline <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___element___geometry__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___element___geometry__oblique <: IDS
    var"alpha" :: Union{Missing, Real, Function}
    var"length_alpha" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"length_beta" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    var"beta" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__coil___element___geometry__oblique(var"alpha"=missing, var"length_alpha"=missing, var"r"=missing, var"length_beta"=missing, var"z"=missing, var"beta"=missing, _parent=WeakRef(missing))
        ids = new(var"alpha", var"length_alpha", var"r", var"length_beta", var"z", var"beta", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___element___geometry__arcs_of_circle <: IDS
    var"curvature_radii" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___element___geometry__arcs_of_circle(var"curvature_radii"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"curvature_radii", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___element___geometry <: IDS
    var"arcs_of_circle" :: pf_active__coil___element___geometry__arcs_of_circle
    var"rectangle" :: pf_active__coil___element___geometry__rectangle
    var"outline" :: pf_active__coil___element___geometry__outline
    var"geometry_type" :: Union{Missing, Integer, Function}
    var"oblique" :: pf_active__coil___element___geometry__oblique
    _parent :: WeakRef
    function pf_active__coil___element___geometry(var"arcs_of_circle"=pf_active__coil___element___geometry__arcs_of_circle(), var"rectangle"=pf_active__coil___element___geometry__rectangle(), var"outline"=pf_active__coil___element___geometry__outline(), var"geometry_type"=missing, var"oblique"=pf_active__coil___element___geometry__oblique(), _parent=WeakRef(missing))
        ids = new(var"arcs_of_circle", var"rectangle", var"outline", var"geometry_type", var"oblique", _parent)
        assign_expressions(ids)
        setfield!(ids.arcs_of_circle, :_parent, WeakRef(ids))
        setfield!(ids.rectangle, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.oblique, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__coil___element <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"turns_with_sign" :: Union{Missing, Real, Function}
    var"area" :: Union{Missing, Real, Function}
    var"geometry" :: pf_active__coil___element___geometry
    var"identifier" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function pf_active__coil___element(var"name"=missing, var"turns_with_sign"=missing, var"area"=missing, var"geometry"=pf_active__coil___element___geometry(), var"identifier"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"turns_with_sign", var"area", var"geometry", var"identifier", _parent)
        assign_expressions(ids)
        setfield!(ids.geometry, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__coil___current <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___current(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___b_field_max_timed <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___b_field_max_timed(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil <: IDSvectorElement
    var"b_field_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"energy_limit_max" :: Union{Missing, Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current" :: pf_active__coil___current
    var"name" :: Union{Missing, String, Function}
    var"b_field_max_timed" :: pf_active__coil___b_field_max_timed
    var"resistance" :: Union{Missing, Real, Function}
    var"voltage" :: pf_active__coil___voltage
    var"identifier" :: Union{Missing, String, Function}
    var"current_limit_max" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"element" :: IDSvector{T} where {T<:pf_active__coil___element}
    _parent :: WeakRef
    function pf_active__coil(var"b_field_max"=missing, var"energy_limit_max"=missing, var"temperature"=missing, var"current"=pf_active__coil___current(), var"name"=missing, var"b_field_max_timed"=pf_active__coil___b_field_max_timed(), var"resistance"=missing, var"voltage"=pf_active__coil___voltage(), var"identifier"=missing, var"current_limit_max"=missing, var"element"=IDSvector(pf_active__coil___element[]), _parent=WeakRef(missing))
        ids = new(var"b_field_max", var"energy_limit_max", var"temperature", var"current", var"name", var"b_field_max_timed", var"resistance", var"voltage", var"identifier", var"current_limit_max", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.b_field_max_timed, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function pf_active__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__code <: IDS
    var"library" :: IDSvector{T} where {T<:pf_active__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function pf_active__code(var"library"=IDSvector(pf_active__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__circuit___voltage <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__circuit___voltage(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__circuit___current <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__circuit___current(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"data", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__circuit <: IDSvectorElement
    var"voltage" :: pf_active__circuit___voltage
    var"name" :: Union{Missing, String, Function}
    var"connections" :: Union{Missing, AbstractArray{T, 2} where T<:Integer, Function}
    var"identifier" :: Union{Missing, String, Function}
    var"type" :: Union{Missing, String, Function}
    var"current" :: pf_active__circuit___current
    _parent :: WeakRef
    function pf_active__circuit(var"voltage"=pf_active__circuit___voltage(), var"name"=missing, var"connections"=missing, var"identifier"=missing, var"type"=missing, var"current"=pf_active__circuit___current(), _parent=WeakRef(missing))
        ids = new(var"voltage", var"name", var"connections", var"identifier", var"type", var"current", _parent)
        assign_expressions(ids)
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active <: IDS
    var"vertical_force" :: IDSvector{T} where {T<:pf_active__vertical_force}
    var"latency" :: Union{Missing, Real, Function}
    var"ids_properties" :: pf_active__ids_properties
    var"coil" :: IDSvector{T} where {T<:pf_active__coil}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"code" :: pf_active__code
    var"global_quantities" :: pf_active__global_quantities
    var"radial_force" :: IDSvector{T} where {T<:pf_active__radial_force}
    var"circuit" :: IDSvector{T} where {T<:pf_active__circuit}
    var"supply" :: IDSvector{T} where {T<:pf_active__supply}
    _parent :: WeakRef
    function pf_active(var"vertical_force"=IDSvector(pf_active__vertical_force[]), var"latency"=missing, var"ids_properties"=pf_active__ids_properties(), var"coil"=IDSvector(pf_active__coil[]), var"time"=missing, var"code"=pf_active__code(), var"global_quantities"=pf_active__global_quantities(), var"radial_force"=IDSvector(pf_active__radial_force[]), var"circuit"=IDSvector(pf_active__circuit[]), var"supply"=IDSvector(pf_active__supply[]), _parent=WeakRef(missing))
        ids = new(var"vertical_force", var"latency", var"ids_properties", var"coil", var"time", var"code", var"global_quantities", var"radial_force", var"circuit", var"supply", _parent)
        assign_expressions(ids)
        setfield!(ids.vertical_force, :_parent, WeakRef(ids))
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.coil, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.radial_force, :_parent, WeakRef(ids))
        setfield!(ids.circuit, :_parent, WeakRef(ids))
        setfield!(ids.supply, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__vacuum_toroidal_field <: IDS
    var"b0" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r0" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        ids = new(var"b0", var"r0", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_2d___grid_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_2d___grid_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_2d___grid <: IDS
    var"volume_element" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"dim2" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dim1" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_2d___grid(var"volume_element"=missing, var"dim2"=missing, var"dim1"=missing, _parent=WeakRef(missing))
        ids = new(var"volume_element", var"dim2", var"dim1", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_2d <: IDSvectorElement
    var"psi" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_field_r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"theta" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_field_z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"j_tor" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"phi" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_field_tor" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid" :: equilibrium__time_slice___profiles_2d___grid
    var"grid_type" :: equilibrium__time_slice___profiles_2d___grid_type
    var"j_parallel" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_tor" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_2d(var"psi"=missing, var"b_field_r"=missing, var"r"=missing, var"b_r"=missing, var"theta"=missing, var"b_field_z"=missing, var"j_tor"=missing, var"phi"=missing, var"z"=missing, var"b_field_tor"=missing, var"b_z"=missing, var"grid"=equilibrium__time_slice___profiles_2d___grid(), var"grid_type"=equilibrium__time_slice___profiles_2d___grid_type(), var"j_parallel"=missing, var"b_tor"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"b_field_r", var"r", var"b_r", var"theta", var"b_field_z", var"j_tor", var"phi", var"z", var"b_field_tor", var"b_z", var"grid", var"grid_type", var"j_parallel", var"b_tor", _parent)
        assign_expressions(ids)
        setfield!(ids.grid, :_parent, WeakRef(ids))
        setfield!(ids.grid_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_1d__geometric_axis <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_1d__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_1d <: IDS
    var"b_field_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dvolume_drho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm9" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dpsi_drho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"surface" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"magnetic_shear" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_field_min" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"darea_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm3" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_upper_inner" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_lower_inner" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"elongation" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"beta_pol" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_field_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm6" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm8" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dpressure_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"triangularity_upper" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"darea_drho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"area" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"trapped_fraction" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dvolume_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_min" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"f" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"mass_density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r_outboard" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm4" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_lower_outer" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"triangularity_lower" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm2" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_volume_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm1" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm5" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"f_df_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r_inboard" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"q" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm7" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_upper_outer" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"geometric_axis" :: equilibrium__time_slice___profiles_1d__geometric_axis
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_1d(var"b_field_max"=missing, var"dvolume_drho_tor"=missing, var"gm9"=missing, var"dpsi_drho_tor"=missing, var"surface"=missing, var"rho_tor"=missing, var"magnetic_shear"=missing, var"b_average"=missing, var"b_field_min"=missing, var"darea_dpsi"=missing, var"gm3"=missing, var"squareness_upper_inner"=missing, var"squareness_lower_inner"=missing, var"rho_tor_norm"=missing, var"elongation"=missing, var"beta_pol"=missing, var"b_field_average"=missing, var"j_parallel"=missing, var"gm6"=missing, var"psi"=missing, var"gm8"=missing, var"dpressure_dpsi"=missing, var"triangularity_upper"=missing, var"darea_drho_tor"=missing, var"area"=missing, var"trapped_fraction"=missing, var"volume"=missing, var"dvolume_dpsi"=missing, var"b_min"=missing, var"f"=missing, var"mass_density"=missing, var"r_outboard"=missing, var"gm4"=missing, var"phi"=missing, var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"gm2"=missing, var"rho_volume_norm"=missing, var"gm1"=missing, var"gm5"=missing, var"b_max"=missing, var"f_df_dpsi"=missing, var"j_tor"=missing, var"r_inboard"=missing, var"q"=missing, var"gm7"=missing, var"pressure"=missing, var"squareness_upper_outer"=missing, var"geometric_axis"=equilibrium__time_slice___profiles_1d__geometric_axis(), _parent=WeakRef(missing))
        ids = new(var"b_field_max", var"dvolume_drho_tor", var"gm9", var"dpsi_drho_tor", var"surface", var"rho_tor", var"magnetic_shear", var"b_average", var"b_field_min", var"darea_dpsi", var"gm3", var"squareness_upper_inner", var"squareness_lower_inner", var"rho_tor_norm", var"elongation", var"beta_pol", var"b_field_average", var"j_parallel", var"gm6", var"psi", var"gm8", var"dpressure_dpsi", var"triangularity_upper", var"darea_drho_tor", var"area", var"trapped_fraction", var"volume", var"dvolume_dpsi", var"b_min", var"f", var"mass_density", var"r_outboard", var"gm4", var"phi", var"squareness_lower_outer", var"triangularity_lower", var"gm2", var"rho_volume_norm", var"gm1", var"gm5", var"b_max", var"f_df_dpsi", var"j_tor", var"r_inboard", var"q", var"gm7", var"pressure", var"squareness_upper_outer", var"geometric_axis", _parent)
        assign_expressions(ids)
        setfield!(ids.geometric_axis, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities__q_min <: IDS
    var"value" :: Union{Missing, Real, Function}
    var"rho_tor_norm" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities__q_min(var"value"=missing, var"rho_tor_norm"=missing, _parent=WeakRef(missing))
        ids = new(var"value", var"rho_tor_norm", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities__magnetic_axis <: IDS
    var"b_field_tor" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    var"b_tor" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities__magnetic_axis(var"b_field_tor"=missing, var"r"=missing, var"z"=missing, var"b_tor"=missing, _parent=WeakRef(missing))
        ids = new(var"b_field_tor", var"r", var"z", var"b_tor", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities__current_centre <: IDS
    var"velocity_z" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities__current_centre(var"velocity_z"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"velocity_z", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities <: IDS
    var"ip" :: Union{Missing, Real, Function}
    var"li_3" :: Union{Missing, Real, Function}
    var"beta_tor" :: Union{Missing, Real, Function}
    var"surface" :: Union{Missing, Real, Function}
    var"magnetic_axis" :: equilibrium__time_slice___global_quantities__magnetic_axis
    var"energy_mhd" :: Union{Missing, Real, Function}
    var"psi_boundary" :: Union{Missing, Real, Function}
    var"length_pol" :: Union{Missing, Real, Function}
    var"area" :: Union{Missing, Real, Function}
    var"psi_external_average" :: Union{Missing, Real, Function}
    var"q_95" :: Union{Missing, Real, Function}
    var"q_axis" :: Union{Missing, Real, Function}
    var"psi_axis" :: Union{Missing, Real, Function}
    var"w_mhd" :: Union{Missing, Real, Function}
    var"volume" :: Union{Missing, Real, Function}
    var"plasma_inductance" :: Union{Missing, Real, Function}
    var"beta_pol" :: Union{Missing, Real, Function}
    var"beta_normal" :: Union{Missing, Real, Function}
    var"current_centre" :: equilibrium__time_slice___global_quantities__current_centre
    var"q_min" :: equilibrium__time_slice___global_quantities__q_min
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities(var"ip"=missing, var"li_3"=missing, var"beta_tor"=missing, var"surface"=missing, var"magnetic_axis"=equilibrium__time_slice___global_quantities__magnetic_axis(), var"energy_mhd"=missing, var"psi_boundary"=missing, var"length_pol"=missing, var"area"=missing, var"psi_external_average"=missing, var"q_95"=missing, var"q_axis"=missing, var"psi_axis"=missing, var"w_mhd"=missing, var"volume"=missing, var"plasma_inductance"=missing, var"beta_pol"=missing, var"beta_normal"=missing, var"current_centre"=equilibrium__time_slice___global_quantities__current_centre(), var"q_min"=equilibrium__time_slice___global_quantities__q_min(), _parent=WeakRef(missing))
        ids = new(var"ip", var"li_3", var"beta_tor", var"surface", var"magnetic_axis", var"energy_mhd", var"psi_boundary", var"length_pol", var"area", var"psi_external_average", var"q_95", var"q_axis", var"psi_axis", var"w_mhd", var"volume", var"plasma_inductance", var"beta_pol", var"beta_normal", var"current_centre", var"q_min", _parent)
        assign_expressions(ids)
        setfield!(ids.magnetic_axis, :_parent, WeakRef(ids))
        setfield!(ids.current_centre, :_parent, WeakRef(ids))
        setfield!(ids.q_min, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___z <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___z(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___theta <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___theta(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___r <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___r(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___psi <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___psi(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___phi <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___phi(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___j_tor <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___j_tor(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___j_parallel <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___j_parallel(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary <: IDSvectorElement
    var"neighbours" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"neighbours", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object <: IDSvectorElement
    var"nodes" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measure" :: Union{Missing, Real, Function}
    var"geometry" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"boundary" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=IDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        ids = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension <: IDSvectorElement
    var"object" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension(var"object"=IDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__space___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__space___geometry_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__space <: IDSvectorElement
    var"coordinates_type" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"geometry_type" :: equilibrium__time_slice___ggd___grid__space___geometry_type
    var"identifier" :: equilibrium__time_slice___ggd___grid__space___identifier
    var"objects_per_dimension" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__time_slice___ggd___grid__space___geometry_type(), var"identifier"=equilibrium__time_slice___ggd___grid__space___identifier(), var"objects_per_dimension"=IDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[]), _parent=WeakRef(missing))
        ids = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(ids)
        setfield!(ids.geometry_type, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.objects_per_dimension, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__grid_subset___metric <: IDS
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__grid_subset___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element___object <: IDSvectorElement
    var"dimension" :: Union{Missing, Integer, Function}
    var"space" :: Union{Missing, Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"dimension", var"space", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element <: IDSvectorElement
    var"object" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element___object}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__grid_subset___element(var"object"=IDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__grid_subset___base <: IDSvectorElement
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid__grid_subset <: IDSvectorElement
    var"base" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___base}
    var"metric" :: equilibrium__time_slice___ggd___grid__grid_subset___metric
    var"dimension" :: Union{Missing, Integer, Function}
    var"identifier" :: equilibrium__time_slice___ggd___grid__grid_subset___identifier
    var"element" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid__grid_subset(var"base"=IDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[]), var"metric"=equilibrium__time_slice___ggd___grid__grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__time_slice___ggd___grid__grid_subset___identifier(), var"element"=IDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[]), _parent=WeakRef(missing))
        ids = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.base, :_parent, WeakRef(ids))
        setfield!(ids.metric, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___grid <: IDS
    var"grid_subset" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset}
    var"space" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space}
    var"identifier" :: equilibrium__time_slice___ggd___grid__identifier
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___grid(var"grid_subset"=IDSvector(equilibrium__time_slice___ggd___grid__grid_subset[]), var"space"=IDSvector(equilibrium__time_slice___ggd___grid__space[]), var"identifier"=equilibrium__time_slice___ggd___grid__identifier(), _parent=WeakRef(missing))
        ids = new(var"grid_subset", var"space", var"identifier", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_subset, :_parent, WeakRef(ids))
        setfield!(ids.space, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___b_field_z <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___b_field_z(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___b_field_tor <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___b_field_tor(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd___b_field_r <: IDSvectorElement
    var"grid_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd___b_field_r(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___ggd <: IDSvectorElement
    var"b_field_z" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___b_field_z}
    var"psi" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___psi}
    var"theta" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___theta}
    var"z" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___z}
    var"phi" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___phi}
    var"j_tor" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___j_tor}
    var"grid" :: equilibrium__time_slice___ggd___grid
    var"b_field_tor" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___b_field_tor}
    var"b_field_r" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___b_field_r}
    var"r" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___r}
    var"j_parallel" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd___j_parallel}
    _parent :: WeakRef
    function equilibrium__time_slice___ggd(var"b_field_z"=IDSvector(equilibrium__time_slice___ggd___b_field_z[]), var"psi"=IDSvector(equilibrium__time_slice___ggd___psi[]), var"theta"=IDSvector(equilibrium__time_slice___ggd___theta[]), var"z"=IDSvector(equilibrium__time_slice___ggd___z[]), var"phi"=IDSvector(equilibrium__time_slice___ggd___phi[]), var"j_tor"=IDSvector(equilibrium__time_slice___ggd___j_tor[]), var"grid"=equilibrium__time_slice___ggd___grid(), var"b_field_tor"=IDSvector(equilibrium__time_slice___ggd___b_field_tor[]), var"b_field_r"=IDSvector(equilibrium__time_slice___ggd___b_field_r[]), var"r"=IDSvector(equilibrium__time_slice___ggd___r[]), var"j_parallel"=IDSvector(equilibrium__time_slice___ggd___j_parallel[]), _parent=WeakRef(missing))
        ids = new(var"b_field_z", var"psi", var"theta", var"z", var"phi", var"j_tor", var"grid", var"b_field_tor", var"b_field_r", var"r", var"j_parallel", _parent)
        assign_expressions(ids)
        setfield!(ids.b_field_z, :_parent, WeakRef(ids))
        setfield!(ids.psi, :_parent, WeakRef(ids))
        setfield!(ids.theta, :_parent, WeakRef(ids))
        setfield!(ids.z, :_parent, WeakRef(ids))
        setfield!(ids.phi, :_parent, WeakRef(ids))
        setfield!(ids.j_tor, :_parent, WeakRef(ids))
        setfield!(ids.grid, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor, :_parent, WeakRef(ids))
        setfield!(ids.b_field_r, :_parent, WeakRef(ids))
        setfield!(ids.r, :_parent, WeakRef(ids))
        setfield!(ids.j_parallel, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___coordinate_system__grid_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___coordinate_system__grid_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___coordinate_system__grid <: IDS
    var"volume_element" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"dim2" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dim1" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___coordinate_system__grid(var"volume_element"=missing, var"dim2"=missing, var"dim1"=missing, _parent=WeakRef(missing))
        ids = new(var"volume_element", var"dim2", var"dim1", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___coordinate_system <: IDS
    var"jacobian" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g13_covariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g11_contravariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g13_contravariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g12_contravariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g22_contravariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g33_contravariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g22_covariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 4} where T<:Real, Function}
    var"g12_covariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g33_covariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid" :: equilibrium__time_slice___coordinate_system__grid
    var"g23_covariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"g11_covariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 4} where T<:Real, Function}
    var"g23_contravariant" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_type" :: equilibrium__time_slice___coordinate_system__grid_type
    _parent :: WeakRef
    function equilibrium__time_slice___coordinate_system(var"jacobian"=missing, var"g13_covariant"=missing, var"g11_contravariant"=missing, var"g13_contravariant"=missing, var"r"=missing, var"g12_contravariant"=missing, var"g22_contravariant"=missing, var"z"=missing, var"g33_contravariant"=missing, var"g22_covariant"=missing, var"tensor_contravariant"=missing, var"g12_covariant"=missing, var"g33_covariant"=missing, var"grid"=equilibrium__time_slice___coordinate_system__grid(), var"g23_covariant"=missing, var"g11_covariant"=missing, var"tensor_covariant"=missing, var"g23_contravariant"=missing, var"grid_type"=equilibrium__time_slice___coordinate_system__grid_type(), _parent=WeakRef(missing))
        ids = new(var"jacobian", var"g13_covariant", var"g11_contravariant", var"g13_contravariant", var"r", var"g12_contravariant", var"g22_contravariant", var"z", var"g33_contravariant", var"g22_covariant", var"tensor_contravariant", var"g12_covariant", var"g33_covariant", var"grid", var"g23_covariant", var"g11_covariant", var"tensor_covariant", var"g23_contravariant", var"grid_type", _parent)
        assign_expressions(ids)
        setfield!(ids.grid, :_parent, WeakRef(ids))
        setfield!(ids.grid_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___convergence <: IDS
    var"iterations_n" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___convergence(var"iterations_n"=missing, _parent=WeakRef(missing))
        ids = new(var"iterations_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__x_point___position_reconstructed <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__x_point___position_reconstructed(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__x_point___position_measured <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__x_point___position_measured(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__x_point <: IDSvectorElement
    var"chi_squared_z" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"position_measured" :: equilibrium__time_slice___constraints__x_point___position_measured
    var"time_measurement" :: Union{Missing, Real, Function}
    var"chi_squared_r" :: Union{Missing, Real, Function}
    var"position_reconstructed" :: equilibrium__time_slice___constraints__x_point___position_reconstructed
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__x_point(var"chi_squared_z"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"position_measured"=equilibrium__time_slice___constraints__x_point___position_measured(), var"time_measurement"=missing, var"chi_squared_r"=missing, var"position_reconstructed"=equilibrium__time_slice___constraints__x_point___position_reconstructed(), _parent=WeakRef(missing))
        ids = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        assign_expressions(ids)
        setfield!(ids.position_measured, :_parent, WeakRef(ids))
        setfield!(ids.position_reconstructed, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__strike_point___position_reconstructed <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__strike_point___position_reconstructed(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__strike_point___position_measured <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__strike_point___position_measured(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__strike_point <: IDSvectorElement
    var"chi_squared_z" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"position_measured" :: equilibrium__time_slice___constraints__strike_point___position_measured
    var"time_measurement" :: Union{Missing, Real, Function}
    var"chi_squared_r" :: Union{Missing, Real, Function}
    var"position_reconstructed" :: equilibrium__time_slice___constraints__strike_point___position_reconstructed
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__strike_point(var"chi_squared_z"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"position_measured"=equilibrium__time_slice___constraints__strike_point___position_measured(), var"time_measurement"=missing, var"chi_squared_r"=missing, var"position_reconstructed"=equilibrium__time_slice___constraints__strike_point___position_reconstructed(), _parent=WeakRef(missing))
        ids = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        assign_expressions(ids)
        setfield!(ids.position_measured, :_parent, WeakRef(ids))
        setfield!(ids.position_reconstructed, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__q___position <: IDS
    var"phi" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__q___position(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"phi", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__q <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"position" :: equilibrium__time_slice___constraints__q___position
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__q(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"position"=equilibrium__time_slice___constraints__q___position(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"position", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.position, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__pressure <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__pressure(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__pf_passive_current <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__pf_passive_current(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__pf_current <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__pf_current(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__n_e_line <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__n_e_line(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__n_e <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__n_e(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__mse_polarisation_angle <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__mse_polarisation_angle(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__iron_core_segment <: IDSvectorElement
    var"magnetisation_r" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r
    var"magnetisation_z" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__iron_core_segment(var"magnetisation_r"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(), var"magnetisation_z"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(), _parent=WeakRef(missing))
        ids = new(var"magnetisation_r", var"magnetisation_z", _parent)
        assign_expressions(ids)
        setfield!(ids.magnetisation_r, :_parent, WeakRef(ids))
        setfield!(ids.magnetisation_z, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__ip <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__ip(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__flux_loop <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__flux_loop(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__faraday_angle <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__faraday_angle(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__diamagnetic_flux <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__diamagnetic_flux(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__bpol_probe <: IDSvectorElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__bpol_probe(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__b_field_tor_vacuum_r <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"weight" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__b_field_tor_vacuum_r(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints <: IDS
    var"faraday_angle" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__faraday_angle}
    var"n_e" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__n_e}
    var"ip" :: equilibrium__time_slice___constraints__ip
    var"n_e_line" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__n_e_line}
    var"pf_current" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__pf_current}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__strike_point}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__x_point}
    var"iron_core_segment" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__iron_core_segment}
    var"pressure" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__pressure}
    var"diamagnetic_flux" :: equilibrium__time_slice___constraints__diamagnetic_flux
    var"pf_passive_current" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__pf_passive_current}
    var"bpol_probe" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__bpol_probe}
    var"mse_polarisation_angle" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__mse_polarisation_angle}
    var"q" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__q}
    var"b_field_tor_vacuum_r" :: equilibrium__time_slice___constraints__b_field_tor_vacuum_r
    var"flux_loop" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__flux_loop}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints(var"faraday_angle"=IDSvector(equilibrium__time_slice___constraints__faraday_angle[]), var"n_e"=IDSvector(equilibrium__time_slice___constraints__n_e[]), var"ip"=equilibrium__time_slice___constraints__ip(), var"n_e_line"=IDSvector(equilibrium__time_slice___constraints__n_e_line[]), var"pf_current"=IDSvector(equilibrium__time_slice___constraints__pf_current[]), var"strike_point"=IDSvector(equilibrium__time_slice___constraints__strike_point[]), var"x_point"=IDSvector(equilibrium__time_slice___constraints__x_point[]), var"iron_core_segment"=IDSvector(equilibrium__time_slice___constraints__iron_core_segment[]), var"pressure"=IDSvector(equilibrium__time_slice___constraints__pressure[]), var"diamagnetic_flux"=equilibrium__time_slice___constraints__diamagnetic_flux(), var"pf_passive_current"=IDSvector(equilibrium__time_slice___constraints__pf_passive_current[]), var"bpol_probe"=IDSvector(equilibrium__time_slice___constraints__bpol_probe[]), var"mse_polarisation_angle"=IDSvector(equilibrium__time_slice___constraints__mse_polarisation_angle[]), var"q"=IDSvector(equilibrium__time_slice___constraints__q[]), var"b_field_tor_vacuum_r"=equilibrium__time_slice___constraints__b_field_tor_vacuum_r(), var"flux_loop"=IDSvector(equilibrium__time_slice___constraints__flux_loop[]), _parent=WeakRef(missing))
        ids = new(var"faraday_angle", var"n_e", var"ip", var"n_e_line", var"pf_current", var"strike_point", var"x_point", var"iron_core_segment", var"pressure", var"diamagnetic_flux", var"pf_passive_current", var"bpol_probe", var"mse_polarisation_angle", var"q", var"b_field_tor_vacuum_r", var"flux_loop", _parent)
        assign_expressions(ids)
        setfield!(ids.faraday_angle, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.ip, :_parent, WeakRef(ids))
        setfield!(ids.n_e_line, :_parent, WeakRef(ids))
        setfield!(ids.pf_current, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        setfield!(ids.iron_core_segment, :_parent, WeakRef(ids))
        setfield!(ids.pressure, :_parent, WeakRef(ids))
        setfield!(ids.diamagnetic_flux, :_parent, WeakRef(ids))
        setfield!(ids.pf_passive_current, :_parent, WeakRef(ids))
        setfield!(ids.bpol_probe, :_parent, WeakRef(ids))
        setfield!(ids.mse_polarisation_angle, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_vacuum_r, :_parent, WeakRef(ids))
        setfield!(ids.flux_loop, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__x_point <: IDSvectorElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__strike_point <: IDSvectorElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__outline <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__geometric_axis <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__gap <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"r" :: Union{Missing, Real, Function}
    var"value" :: Union{Missing, Real, Function}
    var"identifier" :: Union{Missing, String, Function}
    var"angle" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__gap(var"name"=missing, var"r"=missing, var"value"=missing, var"identifier"=missing, var"angle"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"r", var"value", var"identifier", var"angle", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__closest_wall_point <: IDS
    var"distance" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__closest_wall_point(var"distance"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"distance", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__active_limiter_point <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__active_limiter_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix <: IDS
    var"psi" :: Union{Missing, Real, Function}
    var"elongation_lower" :: Union{Missing, Real, Function}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__strike_point}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__x_point}
    var"gap" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__gap}
    var"triangularity" :: Union{Missing, Real, Function}
    var"elongation_upper" :: Union{Missing, Real, Function}
    var"triangularity_upper" :: Union{Missing, Real, Function}
    var"outline" :: equilibrium__time_slice___boundary_separatrix__outline
    var"dr_dz_zero_point" :: equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point
    var"squareness_lower_outer" :: Union{Missing, Real, Function}
    var"triangularity_lower" :: Union{Missing, Real, Function}
    var"minor_radius" :: Union{Missing, Real, Function}
    var"squareness_upper_inner" :: Union{Missing, Real, Function}
    var"squareness_upper_outer" :: Union{Missing, Real, Function}
    var"squareness_lower_inner" :: Union{Missing, Real, Function}
    var"geometric_axis" :: equilibrium__time_slice___boundary_separatrix__geometric_axis
    var"elongation" :: Union{Missing, Real, Function}
    var"active_limiter_point" :: equilibrium__time_slice___boundary_separatrix__active_limiter_point
    var"closest_wall_point" :: equilibrium__time_slice___boundary_separatrix__closest_wall_point
    var"type" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix(var"psi"=missing, var"elongation_lower"=missing, var"strike_point"=IDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[]), var"x_point"=IDSvector(equilibrium__time_slice___boundary_separatrix__x_point[]), var"gap"=IDSvector(equilibrium__time_slice___boundary_separatrix__gap[]), var"triangularity"=missing, var"elongation_upper"=missing, var"triangularity_upper"=missing, var"outline"=equilibrium__time_slice___boundary_separatrix__outline(), var"dr_dz_zero_point"=equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(), var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"minor_radius"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"squareness_lower_inner"=missing, var"geometric_axis"=equilibrium__time_slice___boundary_separatrix__geometric_axis(), var"elongation"=missing, var"active_limiter_point"=equilibrium__time_slice___boundary_separatrix__active_limiter_point(), var"closest_wall_point"=equilibrium__time_slice___boundary_separatrix__closest_wall_point(), var"type"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"elongation_lower", var"strike_point", var"x_point", var"gap", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"dr_dz_zero_point", var"squareness_lower_outer", var"triangularity_lower", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"closest_wall_point", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        setfield!(ids.gap, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.dr_dz_zero_point, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis, :_parent, WeakRef(ids))
        setfield!(ids.active_limiter_point, :_parent, WeakRef(ids))
        setfield!(ids.closest_wall_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_secondary_separatrix__x_point <: IDSvectorElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_secondary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_secondary_separatrix__strike_point <: IDSvectorElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_secondary_separatrix__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_secondary_separatrix__outline <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_secondary_separatrix__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_secondary_separatrix <: IDS
    var"psi" :: Union{Missing, Real, Function}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__x_point}
    var"distance_inner_outer" :: Union{Missing, Real, Function}
    var"outline" :: equilibrium__time_slice___boundary_secondary_separatrix__outline
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__strike_point}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_secondary_separatrix(var"psi"=missing, var"x_point"=IDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[]), var"distance_inner_outer"=missing, var"outline"=equilibrium__time_slice___boundary_secondary_separatrix__outline(), var"strike_point"=IDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[]), _parent=WeakRef(missing))
        ids = new(var"psi", var"x_point", var"distance_inner_outer", var"outline", var"strike_point", _parent)
        assign_expressions(ids)
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__x_point <: IDSvectorElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__strike_point <: IDSvectorElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__outline <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__lcfs <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__lcfs(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__geometric_axis <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__active_limiter_point <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__active_limiter_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary <: IDS
    var"psi" :: Union{Missing, Real, Function}
    var"lcfs" :: equilibrium__time_slice___boundary__lcfs
    var"elongation_lower" :: Union{Missing, Real, Function}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary__strike_point}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary__x_point}
    var"triangularity" :: Union{Missing, Real, Function}
    var"elongation_upper" :: Union{Missing, Real, Function}
    var"triangularity_upper" :: Union{Missing, Real, Function}
    var"outline" :: equilibrium__time_slice___boundary__outline
    var"squareness_lower_outer" :: Union{Missing, Real, Function}
    var"triangularity_lower" :: Union{Missing, Real, Function}
    var"psi_norm" :: Union{Missing, Real, Function}
    var"minor_radius" :: Union{Missing, Real, Function}
    var"squareness_upper_inner" :: Union{Missing, Real, Function}
    var"squareness_upper_outer" :: Union{Missing, Real, Function}
    var"squareness_lower_inner" :: Union{Missing, Real, Function}
    var"geometric_axis" :: equilibrium__time_slice___boundary__geometric_axis
    var"elongation" :: Union{Missing, Real, Function}
    var"active_limiter_point" :: equilibrium__time_slice___boundary__active_limiter_point
    var"b_flux_pol_norm" :: Union{Missing, Real, Function}
    var"type" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary(var"psi"=missing, var"lcfs"=equilibrium__time_slice___boundary__lcfs(), var"elongation_lower"=missing, var"strike_point"=IDSvector(equilibrium__time_slice___boundary__strike_point[]), var"x_point"=IDSvector(equilibrium__time_slice___boundary__x_point[]), var"triangularity"=missing, var"elongation_upper"=missing, var"triangularity_upper"=missing, var"outline"=equilibrium__time_slice___boundary__outline(), var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"psi_norm"=missing, var"minor_radius"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"squareness_lower_inner"=missing, var"geometric_axis"=equilibrium__time_slice___boundary__geometric_axis(), var"elongation"=missing, var"active_limiter_point"=equilibrium__time_slice___boundary__active_limiter_point(), var"b_flux_pol_norm"=missing, var"type"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"lcfs", var"elongation_lower", var"strike_point", var"x_point", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"squareness_lower_outer", var"triangularity_lower", var"psi_norm", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"b_flux_pol_norm", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.lcfs, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis, :_parent, WeakRef(ids))
        setfield!(ids.active_limiter_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice <: IDSvectorElement
    var"time" :: Union{Missing, Real, Function}
    var"ggd" :: IDSvector{T} where {T<:equilibrium__time_slice___ggd}
    var"profiles_1d" :: equilibrium__time_slice___profiles_1d
    var"boundary" :: equilibrium__time_slice___boundary
    var"constraints" :: equilibrium__time_slice___constraints
    var"global_quantities" :: equilibrium__time_slice___global_quantities
    var"convergence" :: equilibrium__time_slice___convergence
    var"coordinate_system" :: equilibrium__time_slice___coordinate_system
    var"boundary_secondary_separatrix" :: equilibrium__time_slice___boundary_secondary_separatrix
    var"boundary_separatrix" :: equilibrium__time_slice___boundary_separatrix
    var"profiles_2d" :: IDSvector{T} where {T<:equilibrium__time_slice___profiles_2d}
    _parent :: WeakRef
    function equilibrium__time_slice(var"time"=missing, var"ggd"=IDSvector(equilibrium__time_slice___ggd[]), var"profiles_1d"=equilibrium__time_slice___profiles_1d(), var"boundary"=equilibrium__time_slice___boundary(), var"constraints"=equilibrium__time_slice___constraints(), var"global_quantities"=equilibrium__time_slice___global_quantities(), var"convergence"=equilibrium__time_slice___convergence(), var"coordinate_system"=equilibrium__time_slice___coordinate_system(), var"boundary_secondary_separatrix"=equilibrium__time_slice___boundary_secondary_separatrix(), var"boundary_separatrix"=equilibrium__time_slice___boundary_separatrix(), var"profiles_2d"=IDSvector(equilibrium__time_slice___profiles_2d[]), _parent=WeakRef(missing))
        ids = new(var"time", var"ggd", var"profiles_1d", var"boundary", var"constraints", var"global_quantities", var"convergence", var"coordinate_system", var"boundary_secondary_separatrix", var"boundary_separatrix", var"profiles_2d", _parent)
        assign_expressions(ids)
        setfield!(ids.ggd, :_parent, WeakRef(ids))
        setfield!(ids.profiles_1d, :_parent, WeakRef(ids))
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        setfield!(ids.constraints, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.convergence, :_parent, WeakRef(ids))
        setfield!(ids.coordinate_system, :_parent, WeakRef(ids))
        setfield!(ids.boundary_secondary_separatrix, :_parent, WeakRef(ids))
        setfield!(ids.boundary_separatrix, :_parent, WeakRef(ids))
        setfield!(ids.profiles_2d, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function equilibrium__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: equilibrium__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__ids_properties(var"provider"=missing, var"version_put"=equilibrium__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary <: IDSvectorElement
    var"neighbours" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"neighbours", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object <: IDSvectorElement
    var"nodes" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measure" :: Union{Missing, Real, Function}
    var"geometry" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"boundary" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=IDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        ids = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension <: IDSvectorElement
    var"object" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___space___objects_per_dimension(var"object"=IDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___space___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___space___geometry_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___space <: IDSvectorElement
    var"coordinates_type" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"geometry_type" :: equilibrium__grids_ggd___grid___space___geometry_type
    var"identifier" :: equilibrium__grids_ggd___grid___space___identifier
    var"objects_per_dimension" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__grids_ggd___grid___space___geometry_type(), var"identifier"=equilibrium__grids_ggd___grid___space___identifier(), var"objects_per_dimension"=IDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[]), _parent=WeakRef(missing))
        ids = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(ids)
        setfield!(ids.geometry_type, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.objects_per_dimension, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___grid_subset___metric <: IDS
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___grid_subset___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___grid_subset___element___object <: IDSvectorElement
    var"dimension" :: Union{Missing, Integer, Function}
    var"space" :: Union{Missing, Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"dimension", var"space", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___grid_subset___element <: IDSvectorElement
    var"object" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element___object}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___grid_subset___element(var"object"=IDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___grid_subset___base <: IDSvectorElement
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid___grid_subset <: IDSvectorElement
    var"base" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___base}
    var"metric" :: equilibrium__grids_ggd___grid___grid_subset___metric
    var"dimension" :: Union{Missing, Integer, Function}
    var"identifier" :: equilibrium__grids_ggd___grid___grid_subset___identifier
    var"element" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element}
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid___grid_subset(var"base"=IDSvector(equilibrium__grids_ggd___grid___grid_subset___base[]), var"metric"=equilibrium__grids_ggd___grid___grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__grids_ggd___grid___grid_subset___identifier(), var"element"=IDSvector(equilibrium__grids_ggd___grid___grid_subset___element[]), _parent=WeakRef(missing))
        ids = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.base, :_parent, WeakRef(ids))
        setfield!(ids.metric, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd___grid <: IDSvectorElement
    var"grid_subset" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset}
    var"space" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid___space}
    var"identifier" :: equilibrium__grids_ggd___grid___identifier
    _parent :: WeakRef
    function equilibrium__grids_ggd___grid(var"grid_subset"=IDSvector(equilibrium__grids_ggd___grid___grid_subset[]), var"space"=IDSvector(equilibrium__grids_ggd___grid___space[]), var"identifier"=equilibrium__grids_ggd___grid___identifier(), _parent=WeakRef(missing))
        ids = new(var"grid_subset", var"space", var"identifier", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_subset, :_parent, WeakRef(ids))
        setfield!(ids.space, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__grids_ggd <: IDSvectorElement
    var"time" :: Union{Missing, Real, Function}
    var"grid" :: IDSvector{T} where {T<:equilibrium__grids_ggd___grid}
    _parent :: WeakRef
    function equilibrium__grids_ggd(var"time"=missing, var"grid"=IDSvector(equilibrium__grids_ggd___grid[]), _parent=WeakRef(missing))
        ids = new(var"time", var"grid", _parent)
        assign_expressions(ids)
        setfield!(ids.grid, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function equilibrium__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__code <: IDS
    var"library" :: IDSvector{T} where {T<:equilibrium__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function equilibrium__code(var"library"=IDSvector(equilibrium__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium <: IDS
    var"time_slice" :: IDSvector{T} where {T<:equilibrium__time_slice}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ids_properties" :: equilibrium__ids_properties
    var"grids_ggd" :: IDSvector{T} where {T<:equilibrium__grids_ggd}
    var"vacuum_toroidal_field" :: equilibrium__vacuum_toroidal_field
    var"code" :: equilibrium__code
    _parent :: WeakRef
    function equilibrium(var"time_slice"=IDSvector(equilibrium__time_slice[]), var"time"=missing, var"ids_properties"=equilibrium__ids_properties(), var"grids_ggd"=IDSvector(equilibrium__grids_ggd[]), var"vacuum_toroidal_field"=equilibrium__vacuum_toroidal_field(), var"code"=equilibrium__code(), _parent=WeakRef(missing))
        ids = new(var"time_slice", var"time", var"ids_properties", var"grids_ggd", var"vacuum_toroidal_field", var"code", _parent)
        assign_expressions(ids)
        setfield!(ids.time_slice, :_parent, WeakRef(ids))
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.grids_ggd, :_parent, WeakRef(ids))
        setfield!(ids.vacuum_toroidal_field, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct dataset_description__simulation <: IDS
    var"time_ended" :: Union{Missing, String, Function}
    var"time_begun" :: Union{Missing, String, Function}
    var"time_current" :: Union{Missing, Real, Function}
    var"time_restart" :: Union{Missing, Real, Function}
    var"workflow" :: Union{Missing, String, Function}
    var"comment_after" :: Union{Missing, String, Function}
    var"time_begin" :: Union{Missing, Real, Function}
    var"time_end" :: Union{Missing, Real, Function}
    var"comment_before" :: Union{Missing, String, Function}
    var"time_step" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function dataset_description__simulation(var"time_ended"=missing, var"time_begun"=missing, var"time_current"=missing, var"time_restart"=missing, var"workflow"=missing, var"comment_after"=missing, var"time_begin"=missing, var"time_end"=missing, var"comment_before"=missing, var"time_step"=missing, _parent=WeakRef(missing))
        ids = new(var"time_ended", var"time_begun", var"time_current", var"time_restart", var"workflow", var"comment_after", var"time_begin", var"time_end", var"comment_before", var"time_step", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description__pulse_time_end_epoch <: IDS
    var"nanoseconds" :: Union{Missing, Integer, Function}
    var"seconds" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function dataset_description__pulse_time_end_epoch(var"nanoseconds"=missing, var"seconds"=missing, _parent=WeakRef(missing))
        ids = new(var"nanoseconds", var"seconds", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description__pulse_time_begin_epoch <: IDS
    var"nanoseconds" :: Union{Missing, Integer, Function}
    var"seconds" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function dataset_description__pulse_time_begin_epoch(var"nanoseconds"=missing, var"seconds"=missing, _parent=WeakRef(missing))
        ids = new(var"nanoseconds", var"seconds", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description__parent_entry <: IDS
    var"pulse_type" :: Union{Missing, String, Function}
    var"run" :: Union{Missing, Integer, Function}
    var"machine" :: Union{Missing, String, Function}
    var"pulse" :: Union{Missing, Integer, Function}
    var"user" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description__parent_entry(var"pulse_type"=missing, var"run"=missing, var"machine"=missing, var"pulse"=missing, var"user"=missing, _parent=WeakRef(missing))
        ids = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: dataset_description__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function dataset_description__ids_properties(var"provider"=missing, var"version_put"=dataset_description__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct dataset_description__data_entry <: IDS
    var"pulse_type" :: Union{Missing, String, Function}
    var"run" :: Union{Missing, Integer, Function}
    var"machine" :: Union{Missing, String, Function}
    var"pulse" :: Union{Missing, Integer, Function}
    var"user" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description__data_entry(var"pulse_type"=missing, var"run"=missing, var"machine"=missing, var"pulse"=missing, var"user"=missing, _parent=WeakRef(missing))
        ids = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description <: IDS
    var"pulse_time_begin_epoch" :: dataset_description__pulse_time_begin_epoch
    var"imas_version" :: Union{Missing, String, Function}
    var"ids_properties" :: dataset_description__ids_properties
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dd_version" :: Union{Missing, String, Function}
    var"parent_entry" :: dataset_description__parent_entry
    var"simulation" :: dataset_description__simulation
    var"pulse_time_end_epoch" :: dataset_description__pulse_time_end_epoch
    var"data_entry" :: dataset_description__data_entry
    var"pulse_time_begin" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description(var"pulse_time_begin_epoch"=dataset_description__pulse_time_begin_epoch(), var"imas_version"=missing, var"ids_properties"=dataset_description__ids_properties(), var"time"=missing, var"dd_version"=missing, var"parent_entry"=dataset_description__parent_entry(), var"simulation"=dataset_description__simulation(), var"pulse_time_end_epoch"=dataset_description__pulse_time_end_epoch(), var"data_entry"=dataset_description__data_entry(), var"pulse_time_begin"=missing, _parent=WeakRef(missing))
        ids = new(var"pulse_time_begin_epoch", var"imas_version", var"ids_properties", var"time", var"dd_version", var"parent_entry", var"simulation", var"pulse_time_end_epoch", var"data_entry", var"pulse_time_begin", _parent)
        assign_expressions(ids)
        setfield!(ids.pulse_time_begin_epoch, :_parent, WeakRef(ids))
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.parent_entry, :_parent, WeakRef(ids))
        setfield!(ids.simulation, :_parent, WeakRef(ids))
        setfield!(ids.pulse_time_end_epoch, :_parent, WeakRef(ids))
        setfield!(ids.data_entry, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__vacuum_toroidal_field <: IDS
    var"b0" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r0" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        ids = new(var"b0", var"r0", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_sources__source___species__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__neutral__state__neutral_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_sources__source___species__neutral__state__neutral_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__neutral__state <: IDS
    var"label" :: Union{Missing, String, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"neutral_type" :: core_sources__source___species__neutral__state__neutral_type
    _parent :: WeakRef
    function core_sources__source___species__neutral__state(var"label"=missing, var"electron_configuration"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, var"neutral_type"=core_sources__source___species__neutral__state__neutral_type(), _parent=WeakRef(missing))
        ids = new(var"label", var"electron_configuration", var"vibrational_level", var"vibrational_mode", var"neutral_type", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___species__neutral__element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__neutral__element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__neutral <: IDS
    var"label" :: Union{Missing, String, Function}
    var"state" :: core_sources__source___species__neutral__state
    var"element" :: IDSvector{T} where {T<:core_sources__source___species__neutral__element}
    _parent :: WeakRef
    function core_sources__source___species__neutral(var"label"=missing, var"state"=core_sources__source___species__neutral__state(), var"element"=IDSvector(core_sources__source___species__neutral__element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"state", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___species__ion__state <: IDS
    var"label" :: Union{Missing, String, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"z_min" :: Union{Missing, Real, Function}
    var"z_max" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__ion__state(var"label"=missing, var"electron_configuration"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, var"z_min"=missing, var"z_max"=missing, _parent=WeakRef(missing))
        ids = new(var"label", var"electron_configuration", var"vibrational_level", var"vibrational_mode", var"z_min", var"z_max", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__ion__element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__ion__element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__ion <: IDS
    var"label" :: Union{Missing, String, Function}
    var"state" :: core_sources__source___species__ion__state
    var"z_ion" :: Union{Missing, Real, Function}
    var"element" :: IDSvector{T} where {T<:core_sources__source___species__ion__element}
    _parent :: WeakRef
    function core_sources__source___species__ion(var"label"=missing, var"state"=core_sources__source___species__ion__state(), var"z_ion"=missing, var"element"=IDSvector(core_sources__source___species__ion__element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"state", var"z_ion", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___species <: IDS
    var"ion" :: core_sources__source___species__ion
    var"type" :: core_sources__source___species__type
    var"neutral" :: core_sources__source___species__neutral
    _parent :: WeakRef
    function core_sources__source___species(var"ion"=core_sources__source___species__ion(), var"type"=core_sources__source___species__type(), var"neutral"=core_sources__source___species__neutral(), _parent=WeakRef(missing))
        ids = new(var"ion", var"type", var"neutral", _parent)
        assign_expressions(ids)
        setfield!(ids.ion, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___total_ion_energy_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___total_ion_energy_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral___state___neutral_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral___state___neutral_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral___state <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"neutral_type" :: core_sources__source___profiles_1d___neutral___state___neutral_type
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral___state(var"label"=missing, var"electron_configuration"=missing, var"particles"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, var"neutral_type"=core_sources__source___profiles_1d___neutral___state___neutral_type(), var"energy"=missing, _parent=WeakRef(missing))
        ids = new(var"label", var"electron_configuration", var"particles", var"vibrational_level", var"vibrational_mode", var"neutral_type", var"energy", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral___element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ion_index" :: Union{Missing, Integer, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"state" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___neutral___state}
    var"element" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___neutral___element}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral(var"label"=missing, var"particles"=missing, var"ion_index"=missing, var"multiple_states_flag"=missing, var"energy"=missing, var"state"=IDSvector(core_sources__source___profiles_1d___neutral___state[]), var"element"=IDSvector(core_sources__source___profiles_1d___neutral___element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"particles", var"ion_index", var"multiple_states_flag", var"energy", var"state", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___state___particles_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___state___particles_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___state___neutral_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___state___neutral_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___state___energy_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___state___energy_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___state <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"is_neutral" :: Union{Missing, Integer, Function}
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_min" :: Union{Missing, Real, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"particles_decomposed" :: core_sources__source___profiles_1d___ion___state___particles_decomposed
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"z_max" :: Union{Missing, Real, Function}
    var"energy_decomposed" :: core_sources__source___profiles_1d___ion___state___energy_decomposed
    var"neutral_type" :: core_sources__source___profiles_1d___ion___state___neutral_type
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___state(var"label"=missing, var"particles"=missing, var"vibrational_level"=missing, var"is_neutral"=missing, var"energy"=missing, var"z_min"=missing, var"electron_configuration"=missing, var"particles_decomposed"=core_sources__source___profiles_1d___ion___state___particles_decomposed(), var"vibrational_mode"=missing, var"z_max"=missing, var"energy_decomposed"=core_sources__source___profiles_1d___ion___state___energy_decomposed(), var"neutral_type"=core_sources__source___profiles_1d___ion___state___neutral_type(), _parent=WeakRef(missing))
        ids = new(var"label", var"particles", var"vibrational_level", var"is_neutral", var"energy", var"z_min", var"electron_configuration", var"particles_decomposed", var"vibrational_mode", var"z_max", var"energy_decomposed", var"neutral_type", _parent)
        assign_expressions(ids)
        setfield!(ids.particles_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.energy_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___particles_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___particles_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___momentum__toroidal_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___momentum__toroidal_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___momentum <: IDS
    var"toroidal_decomposed" :: core_sources__source___profiles_1d___ion___momentum__toroidal_decomposed
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___momentum(var"toroidal_decomposed"=core_sources__source___profiles_1d___ion___momentum__toroidal_decomposed(), var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"toroidal_decomposed", var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        setfield!(ids.toroidal_decomposed, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___energy_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___energy_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion___element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"momentum" :: core_sources__source___profiles_1d___ion___momentum
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"neutral_index" :: Union{Missing, Integer, Function}
    var"particles_decomposed" :: core_sources__source___profiles_1d___ion___particles_decomposed
    var"state" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___ion___state}
    var"z_ion" :: Union{Missing, Real, Function}
    var"energy_decomposed" :: core_sources__source___profiles_1d___ion___energy_decomposed
    var"element" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___ion___element}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion(var"label"=missing, var"particles"=missing, var"momentum"=core_sources__source___profiles_1d___ion___momentum(), var"energy"=missing, var"multiple_states_flag"=missing, var"neutral_index"=missing, var"particles_decomposed"=core_sources__source___profiles_1d___ion___particles_decomposed(), var"state"=IDSvector(core_sources__source___profiles_1d___ion___state[]), var"z_ion"=missing, var"energy_decomposed"=core_sources__source___profiles_1d___ion___energy_decomposed(), var"element"=IDSvector(core_sources__source___profiles_1d___ion___element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"particles", var"momentum", var"energy", var"multiple_states_flag", var"neutral_index", var"particles_decomposed", var"state", var"z_ion", var"energy_decomposed", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.momentum, :_parent, WeakRef(ids))
        setfield!(ids.particles_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.energy_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___grid <: IDS
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_boundary" :: Union{Missing, Real, Function}
    var"volume" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"area" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_pol_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"surface" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_magnetic_axis" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___grid(var"psi"=missing, var"psi_boundary"=missing, var"volume"=missing, var"area"=missing, var"rho_pol_norm"=missing, var"rho_tor_norm"=missing, var"surface"=missing, var"rho_tor"=missing, var"psi_magnetic_axis"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"psi_boundary", var"volume", var"area", var"rho_pol_norm", var"rho_tor_norm", var"surface", var"rho_tor", var"psi_magnetic_axis", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___electrons__particles_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___electrons__particles_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___electrons__energy_decomposed <: IDS
    var"explicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"implicit_part" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___electrons__energy_decomposed(var"explicit_part"=missing, var"implicit_part"=missing, _parent=WeakRef(missing))
        ids = new(var"explicit_part", var"implicit_part", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___electrons <: IDS
    var"particles_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particles_decomposed" :: core_sources__source___profiles_1d___electrons__particles_decomposed
    var"energy_decomposed" :: core_sources__source___profiles_1d___electrons__energy_decomposed
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___electrons(var"particles_inside"=missing, var"particles"=missing, var"particles_decomposed"=core_sources__source___profiles_1d___electrons__particles_decomposed(), var"energy_decomposed"=core_sources__source___profiles_1d___electrons__energy_decomposed(), var"energy"=missing, var"power_inside"=missing, _parent=WeakRef(missing))
        ids = new(var"particles_inside", var"particles", var"particles_decomposed", var"energy_decomposed", var"energy", var"power_inside", _parent)
        assign_expressions(ids)
        setfield!(ids.particles_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.energy_decomposed, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d <: IDSvectorElement
    var"time" :: Union{Missing, Real, Function}
    var"neutral" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___neutral}
    var"torque_tor_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ion" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___ion}
    var"total_ion_energy_decomposed" :: core_sources__source___profiles_1d___total_ion_energy_decomposed
    var"current_parallel_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"momentum_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"conductivity_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electrons" :: core_sources__source___profiles_1d___electrons
    var"momentum_tor_j_cross_b_field" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid" :: core_sources__source___profiles_1d___grid
    var"total_ion_energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"total_ion_power_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d(var"time"=missing, var"neutral"=IDSvector(core_sources__source___profiles_1d___neutral[]), var"torque_tor_inside"=missing, var"ion"=IDSvector(core_sources__source___profiles_1d___ion[]), var"total_ion_energy_decomposed"=core_sources__source___profiles_1d___total_ion_energy_decomposed(), var"current_parallel_inside"=missing, var"momentum_tor"=missing, var"conductivity_parallel"=missing, var"electrons"=core_sources__source___profiles_1d___electrons(), var"momentum_tor_j_cross_b_field"=missing, var"grid"=core_sources__source___profiles_1d___grid(), var"total_ion_energy"=missing, var"j_parallel"=missing, var"total_ion_power_inside"=missing, _parent=WeakRef(missing))
        ids = new(var"time", var"neutral", var"torque_tor_inside", var"ion", var"total_ion_energy_decomposed", var"current_parallel_inside", var"momentum_tor", var"conductivity_parallel", var"electrons", var"momentum_tor_j_cross_b_field", var"grid", var"total_ion_energy", var"j_parallel", var"total_ion_power_inside", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        setfield!(ids.ion, :_parent, WeakRef(ids))
        setfield!(ids.total_ion_energy_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        setfield!(ids.grid, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___identifier <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_sources__source___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___global_quantities___electrons <: IDS
    var"particles" :: Union{Missing, Real, Function}
    var"power" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___global_quantities___electrons(var"particles"=missing, var"power"=missing, _parent=WeakRef(missing))
        ids = new(var"particles", var"power", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___global_quantities <: IDSvectorElement
    var"time" :: Union{Missing, Real, Function}
    var"current_parallel" :: Union{Missing, Real, Function}
    var"total_ion_power" :: Union{Missing, Real, Function}
    var"torque_tor" :: Union{Missing, Real, Function}
    var"power" :: Union{Missing, Real, Function}
    var"total_ion_particles" :: Union{Missing, Real, Function}
    var"electrons" :: core_sources__source___global_quantities___electrons
    _parent :: WeakRef
    function core_sources__source___global_quantities(var"time"=missing, var"current_parallel"=missing, var"total_ion_power"=missing, var"torque_tor"=missing, var"power"=missing, var"total_ion_particles"=missing, var"electrons"=core_sources__source___global_quantities___electrons(), _parent=WeakRef(missing))
        ids = new(var"time", var"current_parallel", var"total_ion_power", var"torque_tor", var"power", var"total_ion_particles", var"electrons", _parent)
        assign_expressions(ids)
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source <: IDSvectorElement
    var"global_quantities" :: IDSvector{T} where {T<:core_sources__source___global_quantities}
    var"identifier" :: core_sources__source___identifier
    var"profiles_1d" :: IDSvector{T} where {T<:core_sources__source___profiles_1d}
    var"species" :: core_sources__source___species
    _parent :: WeakRef
    function core_sources__source(var"global_quantities"=IDSvector(core_sources__source___global_quantities[]), var"identifier"=core_sources__source___identifier(), var"profiles_1d"=IDSvector(core_sources__source___profiles_1d[]), var"species"=core_sources__source___species(), _parent=WeakRef(missing))
        ids = new(var"global_quantities", var"identifier", var"profiles_1d", var"species", _parent)
        assign_expressions(ids)
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.profiles_1d, :_parent, WeakRef(ids))
        setfield!(ids.species, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: core_sources__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_sources__ids_properties(var"provider"=missing, var"version_put"=core_sources__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__code <: IDS
    var"library" :: IDSvector{T} where {T<:core_sources__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__code(var"library"=IDSvector(core_sources__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources <: IDS
    var"source" :: IDSvector{T} where {T<:core_sources__source}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ids_properties" :: core_sources__ids_properties
    var"vacuum_toroidal_field" :: core_sources__vacuum_toroidal_field
    var"code" :: core_sources__code
    _parent :: WeakRef
    function core_sources(var"source"=IDSvector(core_sources__source[]), var"time"=missing, var"ids_properties"=core_sources__ids_properties(), var"vacuum_toroidal_field"=core_sources__vacuum_toroidal_field(), var"code"=core_sources__code(), _parent=WeakRef(missing))
        ids = new(var"source", var"time", var"ids_properties", var"vacuum_toroidal_field", var"code", _parent)
        assign_expressions(ids)
        setfield!(ids.source, :_parent, WeakRef(ids))
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.vacuum_toroidal_field, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__vacuum_toroidal_field <: IDS
    var"b0" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r0" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        ids = new(var"b0", var"r0", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___zeff_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___zeff_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___t_i_average_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___t_i_average_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___velocity <: IDS
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___state___velocity <: IDS
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___state___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___state___neutral_type <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___state___neutral_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___state <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"velocity" :: core_profiles__profiles_1d___neutral___state___velocity
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"neutral_type" :: core_profiles__profiles_1d___neutral___state___neutral_type
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___state(var"label"=missing, var"vibrational_level"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"electron_configuration"=missing, var"pressure"=missing, var"density_thermal"=missing, var"vibrational_mode"=missing, var"pressure_fast_parallel"=missing, var"velocity"=core_profiles__profiles_1d___neutral___state___velocity(), var"density"=missing, var"density_fast"=missing, var"neutral_type"=core_profiles__profiles_1d___neutral___state___neutral_type(), _parent=WeakRef(missing))
        ids = new(var"label", var"vibrational_level", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"velocity", var"density", var"density_fast", var"neutral_type", _parent)
        assign_expressions(ids)
        setfield!(ids.velocity, :_parent, WeakRef(ids))
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ion_index" :: Union{Missing, Integer, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"state" :: IDSvector{T} where {T<:core_profiles__profiles_1d___neutral___state}
    var"velocity" :: core_profiles__profiles_1d___neutral___velocity
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"element" :: IDSvector{T} where {T<:core_profiles__profiles_1d___neutral___element}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral(var"label"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"ion_index"=missing, var"multiple_states_flag"=missing, var"pressure_fast_perpendicular"=missing, var"pressure"=missing, var"density_thermal"=missing, var"pressure_fast_parallel"=missing, var"state"=IDSvector(core_profiles__profiles_1d___neutral___state[]), var"velocity"=core_profiles__profiles_1d___neutral___velocity(), var"density"=missing, var"density_fast"=missing, var"element"=IDSvector(core_profiles__profiles_1d___neutral___element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"temperature", var"pressure_thermal", var"ion_index", var"multiple_states_flag", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"pressure_fast_parallel", var"state", var"velocity", var"density", var"density_fast", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.velocity, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___velocity <: IDS
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___temperature_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___temperature_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state___velocity <: IDS
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state___density_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state___density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"rotation_frequency_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_min" :: Union{Missing, Real, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_average_square_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"velocity" :: core_profiles__profiles_1d___ion___state___velocity
    var"z_average" :: Union{Missing, Real, Function}
    var"z_max" :: Union{Missing, Real, Function}
    var"z_square_average" :: Union{Missing, Real, Function}
    var"ionisation_potential" :: Union{Missing, Real, Function}
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_average_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fit" :: core_profiles__profiles_1d___ion___state___density_fit
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state(var"label"=missing, var"vibrational_level"=missing, var"rotation_frequency_tor"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"z_min"=missing, var"electron_configuration"=missing, var"pressure"=missing, var"density_thermal"=missing, var"vibrational_mode"=missing, var"pressure_fast_parallel"=missing, var"z_average_square_1d"=missing, var"velocity"=core_profiles__profiles_1d___ion___state___velocity(), var"z_average"=missing, var"z_max"=missing, var"z_square_average"=missing, var"ionisation_potential"=missing, var"density"=missing, var"density_fast"=missing, var"z_average_1d"=missing, var"density_fit"=core_profiles__profiles_1d___ion___state___density_fit(), _parent=WeakRef(missing))
        ids = new(var"label", var"vibrational_level", var"rotation_frequency_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"z_min", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"z_average_square_1d", var"velocity", var"z_average", var"z_max", var"z_square_average", var"ionisation_potential", var"density", var"density_fast", var"z_average_1d", var"density_fit", _parent)
        assign_expressions(ids)
        setfield!(ids.velocity, :_parent, WeakRef(ids))
        setfield!(ids.density_fit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___element <: IDSvectorElement
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    var"multiplicity" :: Union{Missing, Real, Function}
    var"a" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        ids = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___density_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion <: IDSvectorElement
    var"label" :: Union{Missing, String, Function}
    var"rotation_frequency_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature_validity" :: Union{Missing, Integer, Function}
    var"velocity_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_ion_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"neutral_index" :: Union{Missing, Integer, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_validity" :: Union{Missing, Integer, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"state" :: IDSvector{T} where {T<:core_profiles__profiles_1d___ion___state}
    var"velocity" :: core_profiles__profiles_1d___ion___velocity
    var"z_ion" :: Union{Missing, Real, Function}
    var"temperature_fit" :: core_profiles__profiles_1d___ion___temperature_fit
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"velocity_pol" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fit" :: core_profiles__profiles_1d___ion___density_fit
    var"z_ion_square_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"element" :: IDSvector{T} where {T<:core_profiles__profiles_1d___ion___element}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion(var"label"=missing, var"rotation_frequency_tor"=missing, var"temperature_validity"=missing, var"velocity_tor"=missing, var"temperature"=missing, var"z_ion_1d"=missing, var"pressure_thermal"=missing, var"multiple_states_flag"=missing, var"pressure_fast_perpendicular"=missing, var"neutral_index"=missing, var"pressure"=missing, var"density_thermal"=missing, var"density_validity"=missing, var"pressure_fast_parallel"=missing, var"state"=IDSvector(core_profiles__profiles_1d___ion___state[]), var"velocity"=core_profiles__profiles_1d___ion___velocity(), var"z_ion"=missing, var"temperature_fit"=core_profiles__profiles_1d___ion___temperature_fit(), var"density"=missing, var"velocity_pol"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___ion___density_fit(), var"z_ion_square_1d"=missing, var"element"=IDSvector(core_profiles__profiles_1d___ion___element[]), _parent=WeakRef(missing))
        ids = new(var"label", var"rotation_frequency_tor", var"temperature_validity", var"velocity_tor", var"temperature", var"z_ion_1d", var"pressure_thermal", var"multiple_states_flag", var"pressure_fast_perpendicular", var"neutral_index", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"state", var"velocity", var"z_ion", var"temperature_fit", var"density", var"velocity_pol", var"density_fast", var"density_fit", var"z_ion_square_1d", var"element", _parent)
        assign_expressions(ids)
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.velocity, :_parent, WeakRef(ids))
        setfield!(ids.temperature_fit, :_parent, WeakRef(ids))
        setfield!(ids.density_fit, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___grid <: IDS
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_boundary" :: Union{Missing, Real, Function}
    var"volume" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"area" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_pol_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"surface" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_magnetic_axis" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___grid(var"psi"=missing, var"psi_boundary"=missing, var"volume"=missing, var"area"=missing, var"rho_pol_norm"=missing, var"rho_tor_norm"=missing, var"surface"=missing, var"rho_tor"=missing, var"psi_magnetic_axis"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"psi_boundary", var"volume", var"area", var"rho_pol_norm", var"rho_tor_norm", var"surface", var"rho_tor", var"psi_magnetic_axis", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__velocity <: IDS
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__temperature_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__temperature_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method <: IDS
    var"name" :: Union{Missing, String, Function}
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"description", var"index", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__density_fit <: IDS
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        ids = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons <: IDS
    var"temperature_validity" :: Union{Missing, Integer, Function}
    var"velocity_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_validity" :: Union{Missing, Integer, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"velocity" :: core_profiles__profiles_1d___electrons__velocity
    var"temperature_fit" :: core_profiles__profiles_1d___electrons__temperature_fit
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"velocity_pol" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"collisionality_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fit" :: core_profiles__profiles_1d___electrons__density_fit
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons(var"temperature_validity"=missing, var"velocity_tor"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"pressure"=missing, var"density_thermal"=missing, var"density_validity"=missing, var"pressure_fast_parallel"=missing, var"velocity"=core_profiles__profiles_1d___electrons__velocity(), var"temperature_fit"=core_profiles__profiles_1d___electrons__temperature_fit(), var"density"=missing, var"velocity_pol"=missing, var"collisionality_norm"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___electrons__density_fit(), _parent=WeakRef(missing))
        ids = new(var"temperature_validity", var"velocity_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"velocity", var"temperature_fit", var"density", var"velocity_pol", var"collisionality_norm", var"density_fast", var"density_fit", _parent)
        assign_expressions(ids)
        setfield!(ids.velocity, :_parent, WeakRef(ids))
        setfield!(ids.temperature_fit, :_parent, WeakRef(ids))
        setfield!(ids.density_fit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___e_field <: IDS
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___e_field(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        ids = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d <: IDSvectorElement
    var"pressure_ion_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, Real, Function}
    var"t_i_average_fit" :: core_profiles__profiles_1d___t_i_average_fit
    var"neutral" :: IDSvector{T} where {T<:core_profiles__profiles_1d___neutral}
    var"n_i_thermal_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"magnetic_shear" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ion" :: IDSvector{T} where {T<:core_profiles__profiles_1d___ion}
    var"j_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rotation_frequency_tor_sonic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_parallel_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_non_inductive" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"e_field_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"momentum_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"conductivity_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electrons" :: core_profiles__profiles_1d___electrons
    var"pressure_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"q" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_i_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_ohmic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"grid" :: core_profiles__profiles_1d___grid
    var"phi_potential" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_bootstrap" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"zeff_fit" :: core_profiles__profiles_1d___zeff_fit
    var"pressure_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"e_field" :: core_profiles__profiles_1d___e_field
    var"zeff" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"n_i_total_over_n_e" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d(var"pressure_ion_total"=missing, var"time"=missing, var"t_i_average_fit"=core_profiles__profiles_1d___t_i_average_fit(), var"neutral"=IDSvector(core_profiles__profiles_1d___neutral[]), var"n_i_thermal_total"=missing, var"magnetic_shear"=missing, var"ion"=IDSvector(core_profiles__profiles_1d___ion[]), var"j_total"=missing, var"rotation_frequency_tor_sonic"=missing, var"pressure_thermal"=missing, var"j_tor"=missing, var"current_parallel_inside"=missing, var"j_non_inductive"=missing, var"e_field_parallel"=missing, var"momentum_tor"=missing, var"conductivity_parallel"=missing, var"electrons"=core_profiles__profiles_1d___electrons(), var"pressure_perpendicular"=missing, var"q"=missing, var"t_i_average"=missing, var"j_ohmic"=missing, var"grid"=core_profiles__profiles_1d___grid(), var"phi_potential"=missing, var"j_bootstrap"=missing, var"zeff_fit"=core_profiles__profiles_1d___zeff_fit(), var"pressure_parallel"=missing, var"e_field"=core_profiles__profiles_1d___e_field(), var"zeff"=missing, var"n_i_total_over_n_e"=missing, _parent=WeakRef(missing))
        ids = new(var"pressure_ion_total", var"time", var"t_i_average_fit", var"neutral", var"n_i_thermal_total", var"magnetic_shear", var"ion", var"j_total", var"rotation_frequency_tor_sonic", var"pressure_thermal", var"j_tor", var"current_parallel_inside", var"j_non_inductive", var"e_field_parallel", var"momentum_tor", var"conductivity_parallel", var"electrons", var"pressure_perpendicular", var"q", var"t_i_average", var"j_ohmic", var"grid", var"phi_potential", var"j_bootstrap", var"zeff_fit", var"pressure_parallel", var"e_field", var"zeff", var"n_i_total_over_n_e", _parent)
        assign_expressions(ids)
        setfield!(ids.t_i_average_fit, :_parent, WeakRef(ids))
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        setfield!(ids.ion, :_parent, WeakRef(ids))
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        setfield!(ids.grid, :_parent, WeakRef(ids))
        setfield!(ids.zeff_fit, :_parent, WeakRef(ids))
        setfield!(ids.e_field, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__ids_properties__version_put <: IDS
    var"access_layer_language" :: Union{Missing, String, Function}
    var"data_dictionary" :: Union{Missing, String, Function}
    var"access_layer" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        ids = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__ids_properties <: IDS
    var"provider" :: Union{Missing, String, Function}
    var"version_put" :: core_profiles__ids_properties__version_put
    var"homogeneous_time" :: Union{Missing, Integer, Function}
    var"source" :: Union{Missing, String, Function}
    var"creation_date" :: Union{Missing, String, Function}
    var"comment" :: Union{Missing, String, Function}
    var"occurrence" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__ids_properties(var"provider"=missing, var"version_put"=core_profiles__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        ids = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_expressions(ids)
        setfield!(ids.version_put, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__global_quantities <: IDS
    var"beta_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"resistive_psi_losses" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ip" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"li_3" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_i_average_peaking" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_e_peaking" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"beta_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_eff_resistive" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ejima" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"energy_diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"li" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_non_inductive" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"v_loop" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"beta_pol" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_bootstrap" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__global_quantities(var"beta_tor_norm"=missing, var"resistive_psi_losses"=missing, var"ip"=missing, var"li_3"=missing, var"t_i_average_peaking"=missing, var"t_e_peaking"=missing, var"beta_tor"=missing, var"z_eff_resistive"=missing, var"ejima"=missing, var"energy_diamagnetic"=missing, var"li"=missing, var"current_non_inductive"=missing, var"v_loop"=missing, var"beta_pol"=missing, var"current_bootstrap"=missing, _parent=WeakRef(missing))
        ids = new(var"beta_tor_norm", var"resistive_psi_losses", var"ip", var"li_3", var"t_i_average_peaking", var"t_e_peaking", var"beta_tor", var"z_eff_resistive", var"ejima", var"energy_diamagnetic", var"li", var"current_non_inductive", var"v_loop", var"beta_pol", var"current_bootstrap", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__code__library <: IDSvectorElement
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__code <: IDS
    var"library" :: IDSvector{T} where {T<:core_profiles__code__library}
    var"name" :: Union{Missing, String, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"commit" :: Union{Missing, String, Function}
    var"repository" :: Union{Missing, String, Function}
    var"output_flag" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"version" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__code(var"library"=IDSvector(core_profiles__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        ids = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_expressions(ids)
        setfield!(ids.library, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ids_properties" :: core_profiles__ids_properties
    var"vacuum_toroidal_field" :: core_profiles__vacuum_toroidal_field
    var"code" :: core_profiles__code
    var"global_quantities" :: core_profiles__global_quantities
    var"profiles_1d" :: IDSvector{T} where {T<:core_profiles__profiles_1d}
    _parent :: WeakRef
    function core_profiles(var"time"=missing, var"ids_properties"=core_profiles__ids_properties(), var"vacuum_toroidal_field"=core_profiles__vacuum_toroidal_field(), var"code"=core_profiles__code(), var"global_quantities"=core_profiles__global_quantities(), var"profiles_1d"=IDSvector(core_profiles__profiles_1d[]), _parent=WeakRef(missing))
        ids = new(var"time", var"ids_properties", var"vacuum_toroidal_field", var"code", var"global_quantities", var"profiles_1d", _parent)
        assign_expressions(ids)
        setfield!(ids.ids_properties, :_parent, WeakRef(ids))
        setfield!(ids.vacuum_toroidal_field, :_parent, WeakRef(ids))
        setfield!(ids.code, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.profiles_1d, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct dd <: IDS
    var"summary" :: Union{Missing, summary}
    var"core_sources" :: Union{Missing, core_sources}
    var"equilibrium" :: Union{Missing, equilibrium}
    var"pf_active" :: Union{Missing, pf_active}
    var"core_profiles" :: Union{Missing, core_profiles}
    var"radial_build" :: Union{Missing, radial_build}
    var"wall" :: Union{Missing, wall}
    var"dataset_description" :: Union{Missing, dataset_description}
    _parent :: WeakRef
    function dd(var"summary"=summary(), var"core_sources"=core_sources(), var"equilibrium"=equilibrium(), var"pf_active"=pf_active(), var"core_profiles"=core_profiles(), var"radial_build"=radial_build(), var"wall"=wall(), var"dataset_description"=dataset_description(), _parent=WeakRef(missing))
        ids = new(var"summary", var"core_sources", var"equilibrium", var"pf_active", var"core_profiles", var"radial_build", var"wall", var"dataset_description", _parent)
        assign_expressions(ids)
        setfield!(ids.summary, :_parent, WeakRef(ids))
        setfield!(ids.core_sources, :_parent, WeakRef(ids))
        setfield!(ids.equilibrium, :_parent, WeakRef(ids))
        setfield!(ids.pf_active, :_parent, WeakRef(ids))
        setfield!(ids.core_profiles, :_parent, WeakRef(ids))
        setfield!(ids.radial_build, :_parent, WeakRef(ids))
        setfield!(ids.wall, :_parent, WeakRef(ids))
        setfield!(ids.dataset_description, :_parent, WeakRef(ids))
        return ids
    end
end

