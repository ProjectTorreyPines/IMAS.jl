include("functionarrays.jl")

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

mutable struct wall__global_quantities__neutral___element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__global_quantities__neutral___element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__global_quantities__neutral <: IDSvectorStaticElement
    var"element" :: IDSvector{T} where {T<:wall__global_quantities__neutral___element}
    var"gas_puff" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"label" :: Union{Missing, String, Function}
    var"particle_flux_from_plasma" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particle_flux_from_wall" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"pumping_speed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"recycling_energy_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"recycling_particles_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"sputtering_chemical_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"sputtering_physical_coefficient" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"wall_inventory" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__global_quantities__neutral(var"element"=IDSvector(wall__global_quantities__neutral___element[]), var"gas_puff"=missing, var"label"=missing, var"particle_flux_from_plasma"=missing, var"particle_flux_from_wall"=missing, var"pumping_speed"=missing, var"recycling_energy_coefficient"=missing, var"recycling_particles_coefficient"=missing, var"sputtering_chemical_coefficient"=missing, var"sputtering_physical_coefficient"=missing, var"wall_inventory"=missing, _parent=WeakRef(missing))
        ids = new(var"element", var"gas_puff", var"label", var"particle_flux_from_plasma", var"particle_flux_from_wall", var"pumping_speed", var"recycling_energy_coefficient", var"recycling_particles_coefficient", var"sputtering_chemical_coefficient", var"sputtering_physical_coefficient", var"wall_inventory", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__global_quantities__electrons <: IDS
    var"gas_puff" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particle_flux_from_plasma" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particle_flux_from_wall" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"power_inner_target" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_outer_target" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pumping_speed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__global_quantities__electrons(var"gas_puff"=missing, var"particle_flux_from_plasma"=missing, var"particle_flux_from_wall"=missing, var"power_inner_target"=missing, var"power_outer_target"=missing, var"pumping_speed"=missing, _parent=WeakRef(missing))
        ids = new(var"gas_puff", var"particle_flux_from_plasma", var"particle_flux_from_wall", var"power_inner_target", var"power_outer_target", var"pumping_speed", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__global_quantities <: IDS
    var"current_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electrons" :: wall__global_quantities__electrons
    var"neutral" :: IDSvector{T} where {T<:wall__global_quantities__neutral}
    var"power_black_body" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_conducted" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_convected" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_currents" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_density_inner_target_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_density_outer_target_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_incident" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_inner_target_ion_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_neutrals" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_radiated" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_recombination_neutrals" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_recombination_plasma" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_to_cooling" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__global_quantities(var"current_tor"=missing, var"electrons"=wall__global_quantities__electrons(), var"neutral"=IDSvector(wall__global_quantities__neutral[]), var"power_black_body"=missing, var"power_conducted"=missing, var"power_convected"=missing, var"power_currents"=missing, var"power_density_inner_target_max"=missing, var"power_density_outer_target_max"=missing, var"power_incident"=missing, var"power_inner_target_ion_total"=missing, var"power_neutrals"=missing, var"power_radiated"=missing, var"power_recombination_neutrals"=missing, var"power_recombination_plasma"=missing, var"power_to_cooling"=missing, var"temperature"=missing, _parent=WeakRef(missing))
        ids = new(var"current_tor", var"electrons", var"neutral", var"power_black_body", var"power_conducted", var"power_convected", var"power_currents", var"power_density_inner_target_max", var"power_density_outer_target_max", var"power_incident", var"power_inner_target_ion_total", var"power_neutrals", var"power_radiated", var"power_recombination_neutrals", var"power_recombination_plasma", var"power_to_cooling", var"temperature", _parent)
        assign_expressions(ids)
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__first_wall_power_flux_peak <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__first_wall_power_flux_peak(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_ggd___type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary <: IDSvectorStaticElement
    var"index" :: Union{Missing, Integer, Function}
    var"neighbours" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary(var"index"=missing, var"neighbours"=missing, _parent=WeakRef(missing))
        ids = new(var"index", var"neighbours", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object <: IDSvectorStaticElement
    var"boundary" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary}
    var"geometry" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"measure" :: Union{Missing, Real, Function}
    var"nodes" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object(var"boundary"=IDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[]), var"geometry"=missing, var"measure"=missing, var"nodes"=missing, _parent=WeakRef(missing))
        ids = new(var"boundary", var"geometry", var"measure", var"nodes", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension <: IDSvectorStaticElement
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space___geometry_type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___space___geometry_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___space <: IDSvectorStaticElement
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___element___object <: IDSvectorStaticElement
    var"dimension" :: Union{Missing, Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"space" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___element___object(var"dimension"=missing, var"index"=missing, var"space"=missing, _parent=WeakRef(missing))
        ids = new(var"dimension", var"index", var"space", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___element <: IDSvectorStaticElement
    var"object" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element___object}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset___element(var"object"=IDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd___grid_subset___base <: IDSvectorStaticElement
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

mutable struct wall__description_ggd___grid_ggd___grid_subset <: IDSvectorStaticElement
    var"base" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___base}
    var"dimension" :: Union{Missing, Integer, Function}
    var"element" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element}
    var"identifier" :: wall__description_ggd___grid_ggd___grid_subset___identifier
    var"metric" :: wall__description_ggd___grid_ggd___grid_subset___metric
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd___grid_subset(var"base"=IDSvector(wall__description_ggd___grid_ggd___grid_subset___base[]), var"dimension"=missing, var"element"=IDSvector(wall__description_ggd___grid_ggd___grid_subset___element[]), var"identifier"=wall__description_ggd___grid_ggd___grid_subset___identifier(), var"metric"=wall__description_ggd___grid_ggd___grid_subset___metric(), _parent=WeakRef(missing))
        ids = new(var"base", var"dimension", var"element", var"identifier", var"metric", _parent)
        assign_expressions(ids)
        setfield!(ids.base, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.metric, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd___grid_ggd <: IDSvectorTimeElement
    var"grid_subset" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset}
    var"identifier" :: wall__description_ggd___grid_ggd___identifier
    var"space" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd___space}
    var"time" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__description_ggd___grid_ggd(var"grid_subset"=IDSvector(wall__description_ggd___grid_ggd___grid_subset[]), var"identifier"=wall__description_ggd___grid_ggd___identifier(), var"space"=IDSvector(wall__description_ggd___grid_ggd___space[]), var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"grid_subset", var"identifier", var"space", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_subset, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.space, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_ggd <: IDSvectorStaticElement
    var"grid_ggd" :: IDSvector{T} where {T<:wall__description_ggd___grid_ggd}
    var"type" :: wall__description_ggd___type
    _parent :: WeakRef
    function wall__description_ggd(var"grid_ggd"=IDSvector(wall__description_ggd___grid_ggd[]), var"type"=wall__description_ggd___type(), _parent=WeakRef(missing))
        ids = new(var"grid_ggd", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_ggd, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
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
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___element___j_tor(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit___element <: IDSvectorStaticElement
    var"j_tor" :: wall__description_2d___vessel__unit___element___j_tor
    var"name" :: Union{Missing, String, Function}
    var"outline" :: wall__description_2d___vessel__unit___element___outline
    var"resistance" :: Union{Missing, Real, Function}
    var"resistivity" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___element(var"j_tor"=wall__description_2d___vessel__unit___element___j_tor(), var"name"=missing, var"outline"=wall__description_2d___vessel__unit___element___outline(), var"resistance"=missing, var"resistivity"=missing, _parent=WeakRef(missing))
        ids = new(var"j_tor", var"name", var"outline", var"resistance", var"resistivity", _parent)
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
    var"centreline" :: wall__description_2d___vessel__unit___annular__centreline
    var"outline_inner" :: wall__description_2d___vessel__unit___annular__outline_inner
    var"outline_outer" :: wall__description_2d___vessel__unit___annular__outline_outer
    var"resistivity" :: Union{Missing, Real, Function}
    var"thickness" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit___annular(var"centreline"=wall__description_2d___vessel__unit___annular__centreline(), var"outline_inner"=wall__description_2d___vessel__unit___annular__outline_inner(), var"outline_outer"=wall__description_2d___vessel__unit___annular__outline_outer(), var"resistivity"=missing, var"thickness"=missing, _parent=WeakRef(missing))
        ids = new(var"centreline", var"outline_inner", var"outline_outer", var"resistivity", var"thickness", _parent)
        assign_expressions(ids)
        setfield!(ids.centreline, :_parent, WeakRef(ids))
        setfield!(ids.outline_inner, :_parent, WeakRef(ids))
        setfield!(ids.outline_outer, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___vessel__unit <: IDSvectorStaticElement
    var"annular" :: wall__description_2d___vessel__unit___annular
    var"element" :: IDSvector{T} where {T<:wall__description_2d___vessel__unit___element}
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__unit(var"annular"=wall__description_2d___vessel__unit___annular(), var"element"=IDSvector(wall__description_2d___vessel__unit___element[]), var"identifier"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"annular", var"element", var"identifier", var"name", _parent)
        assign_expressions(ids)
        setfield!(ids.annular, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___vessel__type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_2d___vessel__type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_2d___type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___mobile__unit___outline <: IDSvectorTimeElement
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall__description_2d___mobile__unit___outline(var"r"=missing, var"time"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"time", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct wall__description_2d___mobile__unit <: IDSvectorStaticElement
    var"closed" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    var"outline" :: IDSvector{T} where {T<:wall__description_2d___mobile__unit___outline}
    var"phi_extensions" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"resistivity" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__description_2d___mobile__unit(var"closed"=missing, var"name"=missing, var"outline"=IDSvector(wall__description_2d___mobile__unit___outline[]), var"phi_extensions"=missing, var"resistivity"=missing, _parent=WeakRef(missing))
        ids = new(var"closed", var"name", var"outline", var"phi_extensions", var"resistivity", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___mobile__type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_2d___mobile__type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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

mutable struct wall__description_2d___limiter__unit <: IDSvectorStaticElement
    var"closed" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    var"outline" :: wall__description_2d___limiter__unit___outline
    var"phi_extensions" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"resistivity" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function wall__description_2d___limiter__unit(var"closed"=missing, var"name"=missing, var"outline"=wall__description_2d___limiter__unit___outline(), var"phi_extensions"=missing, var"resistivity"=missing, _parent=WeakRef(missing))
        ids = new(var"closed", var"name", var"outline", var"phi_extensions", var"resistivity", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall__description_2d___limiter__type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function wall__description_2d___limiter__type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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

mutable struct wall__description_2d <: IDSvectorStaticElement
    var"limiter" :: wall__description_2d___limiter
    var"mobile" :: wall__description_2d___mobile
    var"type" :: wall__description_2d___type
    var"vessel" :: wall__description_2d___vessel
    _parent :: WeakRef
    function wall__description_2d(var"limiter"=wall__description_2d___limiter(), var"mobile"=wall__description_2d___mobile(), var"type"=wall__description_2d___type(), var"vessel"=wall__description_2d___vessel(), _parent=WeakRef(missing))
        ids = new(var"limiter", var"mobile", var"type", var"vessel", _parent)
        assign_expressions(ids)
        setfield!(ids.limiter, :_parent, WeakRef(ids))
        setfield!(ids.mobile, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        setfield!(ids.vessel, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct wall <: IDS
    var"description_2d" :: IDSvector{T} where {T<:wall__description_2d}
    var"description_ggd" :: IDSvector{T} where {T<:wall__description_ggd}
    var"first_wall_power_flux_peak" :: wall__first_wall_power_flux_peak
    var"first_wall_surface_area" :: Union{Missing, Real, Function}
    var"global_quantities" :: wall__global_quantities
    var"temperature_reference" :: wall__temperature_reference
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function wall(var"description_2d"=IDSvector(wall__description_2d[]), var"description_ggd"=IDSvector(wall__description_ggd[]), var"first_wall_power_flux_peak"=wall__first_wall_power_flux_peak(), var"first_wall_surface_area"=missing, var"global_quantities"=wall__global_quantities(), var"temperature_reference"=wall__temperature_reference(), var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"description_2d", var"description_ggd", var"first_wall_power_flux_peak", var"first_wall_surface_area", var"global_quantities", var"temperature_reference", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.description_2d, :_parent, WeakRef(ids))
        setfield!(ids.description_ggd, :_parent, WeakRef(ids))
        setfield!(ids.first_wall_power_flux_peak, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.temperature_reference, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___grid__space___objects_per_dimension___object___boundary <: IDSvectorStaticElement
    var"index" :: Union{Missing, Integer, Function}
    var"neighbours" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function tf__field_map___grid__space___objects_per_dimension___object___boundary(var"index"=missing, var"neighbours"=missing, _parent=WeakRef(missing))
        ids = new(var"index", var"neighbours", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__space___objects_per_dimension___object <: IDSvectorStaticElement
    var"boundary" :: IDSvector{T} where {T<:tf__field_map___grid__space___objects_per_dimension___object___boundary}
    var"geometry" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"measure" :: Union{Missing, Real, Function}
    var"nodes" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function tf__field_map___grid__space___objects_per_dimension___object(var"boundary"=IDSvector(tf__field_map___grid__space___objects_per_dimension___object___boundary[]), var"geometry"=missing, var"measure"=missing, var"nodes"=missing, _parent=WeakRef(missing))
        ids = new(var"boundary", var"geometry", var"measure", var"nodes", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___grid__space___objects_per_dimension <: IDSvectorStaticElement
    var"object" :: IDSvector{T} where {T<:tf__field_map___grid__space___objects_per_dimension___object}
    _parent :: WeakRef
    function tf__field_map___grid__space___objects_per_dimension(var"object"=IDSvector(tf__field_map___grid__space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___grid__space___identifier <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function tf__field_map___grid__space___identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__space___geometry_type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function tf__field_map___grid__space___geometry_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__space <: IDSvectorStaticElement
    var"coordinates_type" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"geometry_type" :: tf__field_map___grid__space___geometry_type
    var"identifier" :: tf__field_map___grid__space___identifier
    var"objects_per_dimension" :: IDSvector{T} where {T<:tf__field_map___grid__space___objects_per_dimension}
    _parent :: WeakRef
    function tf__field_map___grid__space(var"coordinates_type"=missing, var"geometry_type"=tf__field_map___grid__space___geometry_type(), var"identifier"=tf__field_map___grid__space___identifier(), var"objects_per_dimension"=IDSvector(tf__field_map___grid__space___objects_per_dimension[]), _parent=WeakRef(missing))
        ids = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_expressions(ids)
        setfield!(ids.geometry_type, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.objects_per_dimension, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___grid__identifier <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function tf__field_map___grid__identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__grid_subset___metric <: IDS
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___grid__grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__grid_subset___identifier <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function tf__field_map___grid__grid_subset___identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__grid_subset___element___object <: IDSvectorStaticElement
    var"dimension" :: Union{Missing, Integer, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"space" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function tf__field_map___grid__grid_subset___element___object(var"dimension"=missing, var"index"=missing, var"space"=missing, _parent=WeakRef(missing))
        ids = new(var"dimension", var"index", var"space", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__grid_subset___element <: IDSvectorStaticElement
    var"object" :: IDSvector{T} where {T<:tf__field_map___grid__grid_subset___element___object}
    _parent :: WeakRef
    function tf__field_map___grid__grid_subset___element(var"object"=IDSvector(tf__field_map___grid__grid_subset___element___object[]), _parent=WeakRef(missing))
        ids = new(var"object", _parent)
        assign_expressions(ids)
        setfield!(ids.object, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___grid__grid_subset___base <: IDSvectorStaticElement
    var"jacobian" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 3} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___grid__grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        ids = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___grid__grid_subset <: IDSvectorStaticElement
    var"base" :: IDSvector{T} where {T<:tf__field_map___grid__grid_subset___base}
    var"dimension" :: Union{Missing, Integer, Function}
    var"element" :: IDSvector{T} where {T<:tf__field_map___grid__grid_subset___element}
    var"identifier" :: tf__field_map___grid__grid_subset___identifier
    var"metric" :: tf__field_map___grid__grid_subset___metric
    _parent :: WeakRef
    function tf__field_map___grid__grid_subset(var"base"=IDSvector(tf__field_map___grid__grid_subset___base[]), var"dimension"=missing, var"element"=IDSvector(tf__field_map___grid__grid_subset___element[]), var"identifier"=tf__field_map___grid__grid_subset___identifier(), var"metric"=tf__field_map___grid__grid_subset___metric(), _parent=WeakRef(missing))
        ids = new(var"base", var"dimension", var"element", var"identifier", var"metric", _parent)
        assign_expressions(ids)
        setfield!(ids.base, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.metric, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___grid <: IDS
    var"grid_subset" :: IDSvector{T} where {T<:tf__field_map___grid__grid_subset}
    var"identifier" :: tf__field_map___grid__identifier
    var"space" :: IDSvector{T} where {T<:tf__field_map___grid__space}
    _parent :: WeakRef
    function tf__field_map___grid(var"grid_subset"=IDSvector(tf__field_map___grid__grid_subset[]), var"identifier"=tf__field_map___grid__identifier(), var"space"=IDSvector(tf__field_map___grid__space[]), _parent=WeakRef(missing))
        ids = new(var"grid_subset", var"identifier", var"space", _parent)
        assign_expressions(ids)
        setfield!(ids.grid_subset, :_parent, WeakRef(ids))
        setfield!(ids.identifier, :_parent, WeakRef(ids))
        setfield!(ids.space, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__field_map___b_field_z <: IDSvectorStaticElement
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_index" :: Union{Missing, Integer, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___b_field_z(var"coefficients"=missing, var"grid_index"=missing, var"grid_subset_index"=missing, var"values"=missing, _parent=WeakRef(missing))
        ids = new(var"coefficients", var"grid_index", var"grid_subset_index", var"values", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___b_field_tor <: IDSvectorStaticElement
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_index" :: Union{Missing, Integer, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___b_field_tor(var"coefficients"=missing, var"grid_index"=missing, var"grid_subset_index"=missing, var"values"=missing, _parent=WeakRef(missing))
        ids = new(var"coefficients", var"grid_index", var"grid_subset_index", var"values", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___b_field_r <: IDSvectorStaticElement
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_index" :: Union{Missing, Integer, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___b_field_r(var"coefficients"=missing, var"grid_index"=missing, var"grid_subset_index"=missing, var"values"=missing, _parent=WeakRef(missing))
        ids = new(var"coefficients", var"grid_index", var"grid_subset_index", var"values", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___a_field_z <: IDSvectorStaticElement
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_index" :: Union{Missing, Integer, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___a_field_z(var"coefficients"=missing, var"grid_index"=missing, var"grid_subset_index"=missing, var"values"=missing, _parent=WeakRef(missing))
        ids = new(var"coefficients", var"grid_index", var"grid_subset_index", var"values", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___a_field_tor <: IDSvectorStaticElement
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_index" :: Union{Missing, Integer, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___a_field_tor(var"coefficients"=missing, var"grid_index"=missing, var"grid_subset_index"=missing, var"values"=missing, _parent=WeakRef(missing))
        ids = new(var"coefficients", var"grid_index", var"grid_subset_index", var"values", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map___a_field_r <: IDSvectorStaticElement
    var"coefficients" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid_index" :: Union{Missing, Integer, Function}
    var"grid_subset_index" :: Union{Missing, Integer, Function}
    var"values" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__field_map___a_field_r(var"coefficients"=missing, var"grid_index"=missing, var"grid_subset_index"=missing, var"values"=missing, _parent=WeakRef(missing))
        ids = new(var"coefficients", var"grid_index", var"grid_subset_index", var"values", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__field_map <: IDSvectorTimeElement
    var"a_field_r" :: IDSvector{T} where {T<:tf__field_map___a_field_r}
    var"a_field_tor" :: IDSvector{T} where {T<:tf__field_map___a_field_tor}
    var"a_field_z" :: IDSvector{T} where {T<:tf__field_map___a_field_z}
    var"b_field_r" :: IDSvector{T} where {T<:tf__field_map___b_field_r}
    var"b_field_tor" :: IDSvector{T} where {T<:tf__field_map___b_field_tor}
    var"b_field_z" :: IDSvector{T} where {T<:tf__field_map___b_field_z}
    var"grid" :: tf__field_map___grid
    var"time" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function tf__field_map(var"a_field_r"=IDSvector(tf__field_map___a_field_r[]), var"a_field_tor"=IDSvector(tf__field_map___a_field_tor[]), var"a_field_z"=IDSvector(tf__field_map___a_field_z[]), var"b_field_r"=IDSvector(tf__field_map___b_field_r[]), var"b_field_tor"=IDSvector(tf__field_map___b_field_tor[]), var"b_field_z"=IDSvector(tf__field_map___b_field_z[]), var"grid"=tf__field_map___grid(), var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"a_field_r", var"a_field_tor", var"a_field_z", var"b_field_r", var"b_field_tor", var"b_field_z", var"grid", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.a_field_r, :_parent, WeakRef(ids))
        setfield!(ids.a_field_tor, :_parent, WeakRef(ids))
        setfield!(ids.a_field_z, :_parent, WeakRef(ids))
        setfield!(ids.b_field_r, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor, :_parent, WeakRef(ids))
        setfield!(ids.b_field_z, :_parent, WeakRef(ids))
        setfield!(ids.grid, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__delta_b_field_tor_vacuum_r <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__delta_b_field_tor_vacuum_r(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___voltage <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___voltage(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___current <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___current(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___voltage <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___voltage(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___elements__start_points <: IDS
    var"phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___elements__start_points(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"phi", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___elements__intermediate_points <: IDS
    var"phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___elements__intermediate_points(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"phi", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___elements__end_points <: IDS
    var"phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___elements__end_points(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"phi", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___elements__centres <: IDS
    var"phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___elements__centres(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"phi", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___elements <: IDS
    var"centres" :: tf__coil___conductor___elements__centres
    var"end_points" :: tf__coil___conductor___elements__end_points
    var"intermediate_points" :: tf__coil___conductor___elements__intermediate_points
    var"names" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"start_points" :: tf__coil___conductor___elements__start_points
    var"types" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    _parent :: WeakRef
    function tf__coil___conductor___elements(var"centres"=tf__coil___conductor___elements__centres(), var"end_points"=tf__coil___conductor___elements__end_points(), var"intermediate_points"=tf__coil___conductor___elements__intermediate_points(), var"names"=missing, var"start_points"=tf__coil___conductor___elements__start_points(), var"types"=missing, _parent=WeakRef(missing))
        ids = new(var"centres", var"end_points", var"intermediate_points", var"names", var"start_points", var"types", _parent)
        assign_expressions(ids)
        setfield!(ids.centres, :_parent, WeakRef(ids))
        setfield!(ids.end_points, :_parent, WeakRef(ids))
        setfield!(ids.intermediate_points, :_parent, WeakRef(ids))
        setfield!(ids.start_points, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__coil___conductor___current <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___current(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor___cross_section <: IDS
    var"delta_phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"delta_r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"delta_z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__coil___conductor___cross_section(var"delta_phi"=missing, var"delta_r"=missing, var"delta_z"=missing, _parent=WeakRef(missing))
        ids = new(var"delta_phi", var"delta_r", var"delta_z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf__coil___conductor <: IDSvectorStaticElement
    var"cross_section" :: tf__coil___conductor___cross_section
    var"current" :: tf__coil___conductor___current
    var"elements" :: tf__coil___conductor___elements
    var"resistance" :: Union{Missing, Real, Function}
    var"voltage" :: tf__coil___conductor___voltage
    _parent :: WeakRef
    function tf__coil___conductor(var"cross_section"=tf__coil___conductor___cross_section(), var"current"=tf__coil___conductor___current(), var"elements"=tf__coil___conductor___elements(), var"resistance"=missing, var"voltage"=tf__coil___conductor___voltage(), _parent=WeakRef(missing))
        ids = new(var"cross_section", var"current", var"elements", var"resistance", var"voltage", _parent)
        assign_expressions(ids)
        setfield!(ids.cross_section, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.elements, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__coil <: IDSvectorStaticElement
    var"conductor" :: IDSvector{T} where {T<:tf__coil___conductor}
    var"current" :: tf__coil___current
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"resistance" :: Union{Missing, Real, Function}
    var"turns" :: Union{Missing, Real, Function}
    var"voltage" :: tf__coil___voltage
    _parent :: WeakRef
    function tf__coil(var"conductor"=IDSvector(tf__coil___conductor[]), var"current"=tf__coil___current(), var"identifier"=missing, var"name"=missing, var"resistance"=missing, var"turns"=missing, var"voltage"=tf__coil___voltage(), _parent=WeakRef(missing))
        ids = new(var"conductor", var"current", var"identifier", var"name", var"resistance", var"turns", var"voltage", _parent)
        assign_expressions(ids)
        setfield!(ids.conductor, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct tf__b_field_tor_vacuum_r <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf__b_field_tor_vacuum_r(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct tf <: IDS
    var"b_field_tor_vacuum_r" :: tf__b_field_tor_vacuum_r
    var"coil" :: IDSvector{T} where {T<:tf__coil}
    var"coils_n" :: Union{Missing, Integer, Function}
    var"delta_b_field_tor_vacuum_r" :: tf__delta_b_field_tor_vacuum_r
    var"field_map" :: IDSvector{T} where {T<:tf__field_map}
    var"is_periodic" :: Union{Missing, Integer, Function}
    var"latency" :: Union{Missing, Real, Function}
    var"r0" :: Union{Missing, Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function tf(var"b_field_tor_vacuum_r"=tf__b_field_tor_vacuum_r(), var"coil"=IDSvector(tf__coil[]), var"coils_n"=missing, var"delta_b_field_tor_vacuum_r"=tf__delta_b_field_tor_vacuum_r(), var"field_map"=IDSvector(tf__field_map[]), var"is_periodic"=missing, var"latency"=missing, var"r0"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"b_field_tor_vacuum_r", var"coil", var"coils_n", var"delta_b_field_tor_vacuum_r", var"field_map", var"is_periodic", var"latency", var"r0", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.b_field_tor_vacuum_r, :_parent, WeakRef(ids))
        setfield!(ids.coil, :_parent, WeakRef(ids))
        setfield!(ids.delta_b_field_tor_vacuum_r, :_parent, WeakRef(ids))
        setfield!(ids.field_map, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__wall__material <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__wall__material(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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
    var"evaporation" :: summary__wall__evaporation
    var"material" :: summary__wall__material
    _parent :: WeakRef
    function summary__wall(var"evaporation"=summary__wall__evaporation(), var"material"=summary__wall__material(), _parent=WeakRef(missing))
        ids = new(var"evaporation", var"material", _parent)
        assign_expressions(ids)
        setfield!(ids.evaporation, :_parent, WeakRef(ids))
        setfield!(ids.material, :_parent, WeakRef(ids))
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
    var"argon" :: summary__volume_average__n_i__argon
    var"beryllium" :: summary__volume_average__n_i__beryllium
    var"carbon" :: summary__volume_average__n_i__carbon
    var"deuterium" :: summary__volume_average__n_i__deuterium
    var"helium_3" :: summary__volume_average__n_i__helium_3
    var"helium_4" :: summary__volume_average__n_i__helium_4
    var"hydrogen" :: summary__volume_average__n_i__hydrogen
    var"iron" :: summary__volume_average__n_i__iron
    var"krypton" :: summary__volume_average__n_i__krypton
    var"lithium" :: summary__volume_average__n_i__lithium
    var"neon" :: summary__volume_average__n_i__neon
    var"nitrogen" :: summary__volume_average__n_i__nitrogen
    var"oxygen" :: summary__volume_average__n_i__oxygen
    var"tritium" :: summary__volume_average__n_i__tritium
    var"tungsten" :: summary__volume_average__n_i__tungsten
    var"xenon" :: summary__volume_average__n_i__xenon
    _parent :: WeakRef
    function summary__volume_average__n_i(var"argon"=summary__volume_average__n_i__argon(), var"beryllium"=summary__volume_average__n_i__beryllium(), var"carbon"=summary__volume_average__n_i__carbon(), var"deuterium"=summary__volume_average__n_i__deuterium(), var"helium_3"=summary__volume_average__n_i__helium_3(), var"helium_4"=summary__volume_average__n_i__helium_4(), var"hydrogen"=summary__volume_average__n_i__hydrogen(), var"iron"=summary__volume_average__n_i__iron(), var"krypton"=summary__volume_average__n_i__krypton(), var"lithium"=summary__volume_average__n_i__lithium(), var"neon"=summary__volume_average__n_i__neon(), var"nitrogen"=summary__volume_average__n_i__nitrogen(), var"oxygen"=summary__volume_average__n_i__oxygen(), var"tritium"=summary__volume_average__n_i__tritium(), var"tungsten"=summary__volume_average__n_i__tungsten(), var"xenon"=summary__volume_average__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"dn_e_dt" :: summary__volume_average__dn_e_dt
    var"isotope_fraction_hydrogen" :: summary__volume_average__isotope_fraction_hydrogen
    var"meff_hydrogenic" :: summary__volume_average__meff_hydrogenic
    var"n_e" :: summary__volume_average__n_e
    var"n_i" :: summary__volume_average__n_i
    var"n_i_total" :: summary__volume_average__n_i_total
    var"t_e" :: summary__volume_average__t_e
    var"t_i_average" :: summary__volume_average__t_i_average
    var"zeff" :: summary__volume_average__zeff
    _parent :: WeakRef
    function summary__volume_average(var"dn_e_dt"=summary__volume_average__dn_e_dt(), var"isotope_fraction_hydrogen"=summary__volume_average__isotope_fraction_hydrogen(), var"meff_hydrogenic"=summary__volume_average__meff_hydrogenic(), var"n_e"=summary__volume_average__n_e(), var"n_i"=summary__volume_average__n_i(), var"n_i_total"=summary__volume_average__n_i_total(), var"t_e"=summary__volume_average__t_e(), var"t_i_average"=summary__volume_average__t_i_average(), var"zeff"=summary__volume_average__zeff(), _parent=WeakRef(missing))
        ids = new(var"dn_e_dt", var"isotope_fraction_hydrogen", var"meff_hydrogenic", var"n_e", var"n_i", var"n_i_total", var"t_e", var"t_i_average", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.dn_e_dt, :_parent, WeakRef(ids))
        setfield!(ids.isotope_fraction_hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.meff_hydrogenic, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__tag <: IDS
    var"comment" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__tag(var"comment"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"comment", var"name", _parent)
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
    var"heat_flux_e_decay_length" :: summary__scrape_off_layer__heat_flux_e_decay_length
    var"heat_flux_i_decay_length" :: summary__scrape_off_layer__heat_flux_i_decay_length
    var"n_e_decay_length" :: summary__scrape_off_layer__n_e_decay_length
    var"n_i_total_decay_length" :: summary__scrape_off_layer__n_i_total_decay_length
    var"power_radiated" :: summary__scrape_off_layer__power_radiated
    var"pressure_neutral" :: summary__scrape_off_layer__pressure_neutral
    var"t_e_decay_length" :: summary__scrape_off_layer__t_e_decay_length
    var"t_i_average_decay_length" :: summary__scrape_off_layer__t_i_average_decay_length
    _parent :: WeakRef
    function summary__scrape_off_layer(var"heat_flux_e_decay_length"=summary__scrape_off_layer__heat_flux_e_decay_length(), var"heat_flux_i_decay_length"=summary__scrape_off_layer__heat_flux_i_decay_length(), var"n_e_decay_length"=summary__scrape_off_layer__n_e_decay_length(), var"n_i_total_decay_length"=summary__scrape_off_layer__n_i_total_decay_length(), var"power_radiated"=summary__scrape_off_layer__power_radiated(), var"pressure_neutral"=summary__scrape_off_layer__pressure_neutral(), var"t_e_decay_length"=summary__scrape_off_layer__t_e_decay_length(), var"t_i_average_decay_length"=summary__scrape_off_layer__t_i_average_decay_length(), _parent=WeakRef(missing))
        ids = new(var"heat_flux_e_decay_length", var"heat_flux_i_decay_length", var"n_e_decay_length", var"n_i_total_decay_length", var"power_radiated", var"pressure_neutral", var"t_e_decay_length", var"t_i_average_decay_length", _parent)
        assign_expressions(ids)
        setfield!(ids.heat_flux_e_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.heat_flux_i_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.n_e_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total_decay_length, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated, :_parent, WeakRef(ids))
        setfield!(ids.pressure_neutral, :_parent, WeakRef(ids))
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
    var"current" :: summary__runaways__current
    var"particles" :: summary__runaways__particles
    _parent :: WeakRef
    function summary__runaways(var"current"=summary__runaways__current(), var"particles"=summary__runaways__particles(), _parent=WeakRef(missing))
        ids = new(var"current", var"particles", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.particles, :_parent, WeakRef(ids))
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
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position
    var"offset" :: summary__pedestal_fits__mtanh__t_e__offset
    var"pedestal_height" :: summary__pedestal_fits__mtanh__t_e__pedestal_height
    var"pedestal_position" :: summary__pedestal_fits__mtanh__t_e__pedestal_position
    var"pedestal_width" :: summary__pedestal_fits__mtanh__t_e__pedestal_width
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__t_e(var"d_dpsi_norm"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__t_e__d_dpsi_norm_max_position(), var"offset"=summary__pedestal_fits__mtanh__t_e__offset(), var"pedestal_height"=summary__pedestal_fits__mtanh__t_e__pedestal_height(), var"pedestal_position"=summary__pedestal_fits__mtanh__t_e__pedestal_position(), var"pedestal_width"=summary__pedestal_fits__mtanh__t_e__pedestal_width(), _parent=WeakRef(missing))
        ids = new(var"d_dpsi_norm", var"d_dpsi_norm_max", var"d_dpsi_norm_max_position", var"offset", var"pedestal_height", var"pedestal_position", var"pedestal_width", _parent)
        assign_expressions(ids)
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
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
    var"alpha_ratio" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio
    var"t_e_pedestal_top_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter(var"alpha_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_critical(), var"alpha_ratio"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__alpha_ratio(), var"t_e_pedestal_top_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_sauter__t_e_pedestal_top_critical(), _parent=WeakRef(missing))
        ids = new(var"alpha_critical", var"alpha_ratio", var"t_e_pedestal_top_critical", _parent)
        assign_expressions(ids)
        setfield!(ids.alpha_critical, :_parent, WeakRef(ids))
        setfield!(ids.alpha_ratio, :_parent, WeakRef(ids))
        setfield!(ids.t_e_pedestal_top_critical, :_parent, WeakRef(ids))
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
    var"alpha_ratio" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio
    var"t_e_pedestal_top_critical" :: summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__stability__bootstrap_current_hager(var"alpha_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_critical(), var"alpha_ratio"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__alpha_ratio(), var"t_e_pedestal_top_critical"=summary__pedestal_fits__mtanh__stability__bootstrap_current_hager__t_e_pedestal_top_critical(), _parent=WeakRef(missing))
        ids = new(var"alpha_critical", var"alpha_ratio", var"t_e_pedestal_top_critical", _parent)
        assign_expressions(ids)
        setfield!(ids.alpha_critical, :_parent, WeakRef(ids))
        setfield!(ids.alpha_ratio, :_parent, WeakRef(ids))
        setfield!(ids.t_e_pedestal_top_critical, :_parent, WeakRef(ids))
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
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position
    var"offset" :: summary__pedestal_fits__mtanh__pressure_electron__offset
    var"pedestal_height" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_height
    var"pedestal_position" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_position
    var"pedestal_width" :: summary__pedestal_fits__mtanh__pressure_electron__pedestal_width
    var"separatrix" :: summary__pedestal_fits__mtanh__pressure_electron__separatrix
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__pressure_electron(var"d_dpsi_norm"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__pressure_electron__d_dpsi_norm_max_position(), var"offset"=summary__pedestal_fits__mtanh__pressure_electron__offset(), var"pedestal_height"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_height(), var"pedestal_position"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_position(), var"pedestal_width"=summary__pedestal_fits__mtanh__pressure_electron__pedestal_width(), var"separatrix"=summary__pedestal_fits__mtanh__pressure_electron__separatrix(), _parent=WeakRef(missing))
        ids = new(var"d_dpsi_norm", var"d_dpsi_norm_max", var"d_dpsi_norm_max_position", var"offset", var"pedestal_height", var"pedestal_position", var"pedestal_width", var"separatrix", _parent)
        assign_expressions(ids)
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
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
    var"d_dpsi_norm" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm
    var"d_dpsi_norm_max" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position
    var"offset" :: summary__pedestal_fits__mtanh__n_e__offset
    var"pedestal_height" :: summary__pedestal_fits__mtanh__n_e__pedestal_height
    var"pedestal_position" :: summary__pedestal_fits__mtanh__n_e__pedestal_position
    var"pedestal_width" :: summary__pedestal_fits__mtanh__n_e__pedestal_width
    var"separatrix" :: summary__pedestal_fits__mtanh__n_e__separatrix
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh__n_e(var"d_dpsi_norm"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm(), var"d_dpsi_norm_max"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__mtanh__n_e__d_dpsi_norm_max_position(), var"offset"=summary__pedestal_fits__mtanh__n_e__offset(), var"pedestal_height"=summary__pedestal_fits__mtanh__n_e__pedestal_height(), var"pedestal_position"=summary__pedestal_fits__mtanh__n_e__pedestal_position(), var"pedestal_width"=summary__pedestal_fits__mtanh__n_e__pedestal_width(), var"separatrix"=summary__pedestal_fits__mtanh__n_e__separatrix(), _parent=WeakRef(missing))
        ids = new(var"d_dpsi_norm", var"d_dpsi_norm_max", var"d_dpsi_norm_max_position", var"offset", var"pedestal_height", var"pedestal_position", var"pedestal_width", var"separatrix", _parent)
        assign_expressions(ids)
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
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
    var"alpha_electron_pedestal_max" :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max
    var"alpha_electron_pedestal_max_position" :: summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position
    var"b_field_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs
    var"b_field_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs
    var"b_field_pol_pedestal_top_average" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average
    var"b_field_pol_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs
    var"b_field_pol_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs
    var"b_field_tor_pedestal_top_hfs" :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs
    var"b_field_tor_pedestal_top_lfs" :: summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs
    var"beta_pol_pedestal_top_electron_average" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average
    var"beta_pol_pedestal_top_electron_hfs" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs
    var"beta_pol_pedestal_top_electron_lfs" :: summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs
    var"coulomb_factor_pedestal_top" :: summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top
    var"energy_thermal_pedestal_electron" :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron
    var"energy_thermal_pedestal_ion" :: summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion
    var"n_e" :: summary__pedestal_fits__mtanh__n_e
    var"nustar_pedestal_top_electron" :: summary__pedestal_fits__mtanh__nustar_pedestal_top_electron
    var"parameters" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_electron" :: summary__pedestal_fits__mtanh__pressure_electron
    var"rhostar_pedestal_top_electron_hfs" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs
    var"rhostar_pedestal_top_electron_lfs" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs
    var"rhostar_pedestal_top_electron_magnetic_axis" :: summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis
    var"stability" :: summary__pedestal_fits__mtanh__stability
    var"t_e" :: summary__pedestal_fits__mtanh__t_e
    var"volume_inside_pedestal" :: summary__pedestal_fits__mtanh__volume_inside_pedestal
    _parent :: WeakRef
    function summary__pedestal_fits__mtanh(var"alpha_electron_pedestal_max"=summary__pedestal_fits__mtanh__alpha_electron_pedestal_max(), var"alpha_electron_pedestal_max_position"=summary__pedestal_fits__mtanh__alpha_electron_pedestal_max_position(), var"b_field_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_pedestal_top_hfs(), var"b_field_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_pedestal_top_lfs(), var"b_field_pol_pedestal_top_average"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_average(), var"b_field_pol_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_hfs(), var"b_field_pol_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_pol_pedestal_top_lfs(), var"b_field_tor_pedestal_top_hfs"=summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_hfs(), var"b_field_tor_pedestal_top_lfs"=summary__pedestal_fits__mtanh__b_field_tor_pedestal_top_lfs(), var"beta_pol_pedestal_top_electron_average"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_average(), var"beta_pol_pedestal_top_electron_hfs"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_hfs(), var"beta_pol_pedestal_top_electron_lfs"=summary__pedestal_fits__mtanh__beta_pol_pedestal_top_electron_lfs(), var"coulomb_factor_pedestal_top"=summary__pedestal_fits__mtanh__coulomb_factor_pedestal_top(), var"energy_thermal_pedestal_electron"=summary__pedestal_fits__mtanh__energy_thermal_pedestal_electron(), var"energy_thermal_pedestal_ion"=summary__pedestal_fits__mtanh__energy_thermal_pedestal_ion(), var"n_e"=summary__pedestal_fits__mtanh__n_e(), var"nustar_pedestal_top_electron"=summary__pedestal_fits__mtanh__nustar_pedestal_top_electron(), var"parameters"=missing, var"pressure_electron"=summary__pedestal_fits__mtanh__pressure_electron(), var"rhostar_pedestal_top_electron_hfs"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_hfs(), var"rhostar_pedestal_top_electron_lfs"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_lfs(), var"rhostar_pedestal_top_electron_magnetic_axis"=summary__pedestal_fits__mtanh__rhostar_pedestal_top_electron_magnetic_axis(), var"stability"=summary__pedestal_fits__mtanh__stability(), var"t_e"=summary__pedestal_fits__mtanh__t_e(), var"volume_inside_pedestal"=summary__pedestal_fits__mtanh__volume_inside_pedestal(), _parent=WeakRef(missing))
        ids = new(var"alpha_electron_pedestal_max", var"alpha_electron_pedestal_max_position", var"b_field_pedestal_top_hfs", var"b_field_pedestal_top_lfs", var"b_field_pol_pedestal_top_average", var"b_field_pol_pedestal_top_hfs", var"b_field_pol_pedestal_top_lfs", var"b_field_tor_pedestal_top_hfs", var"b_field_tor_pedestal_top_lfs", var"beta_pol_pedestal_top_electron_average", var"beta_pol_pedestal_top_electron_hfs", var"beta_pol_pedestal_top_electron_lfs", var"coulomb_factor_pedestal_top", var"energy_thermal_pedestal_electron", var"energy_thermal_pedestal_ion", var"n_e", var"nustar_pedestal_top_electron", var"parameters", var"pressure_electron", var"rhostar_pedestal_top_electron_hfs", var"rhostar_pedestal_top_electron_lfs", var"rhostar_pedestal_top_electron_magnetic_axis", var"stability", var"t_e", var"volume_inside_pedestal", _parent)
        assign_expressions(ids)
        setfield!(ids.alpha_electron_pedestal_max, :_parent, WeakRef(ids))
        setfield!(ids.alpha_electron_pedestal_max_position, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_average, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_average, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.coulomb_factor_pedestal_top, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_electron, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_ion, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.nustar_pedestal_top_electron, :_parent, WeakRef(ids))
        setfield!(ids.pressure_electron, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_magnetic_axis, :_parent, WeakRef(ids))
        setfield!(ids.stability, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.volume_inside_pedestal, :_parent, WeakRef(ids))
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
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__t_e__d_dpsi_norm
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__t_e__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__linear__t_e__offset
    var"pedestal_height" :: summary__pedestal_fits__linear__t_e__pedestal_height
    var"pedestal_position" :: summary__pedestal_fits__linear__t_e__pedestal_position
    var"pedestal_width" :: summary__pedestal_fits__linear__t_e__pedestal_width
    _parent :: WeakRef
    function summary__pedestal_fits__linear__t_e(var"d_dpsi_norm"=summary__pedestal_fits__linear__t_e__d_dpsi_norm(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__t_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__t_e__offset(), var"pedestal_height"=summary__pedestal_fits__linear__t_e__pedestal_height(), var"pedestal_position"=summary__pedestal_fits__linear__t_e__pedestal_position(), var"pedestal_width"=summary__pedestal_fits__linear__t_e__pedestal_width(), _parent=WeakRef(missing))
        ids = new(var"d_dpsi_norm", var"d_dpsi_norm_max", var"offset", var"pedestal_height", var"pedestal_position", var"pedestal_width", _parent)
        assign_expressions(ids)
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
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
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max
    var"d_dpsi_norm_max_position" :: summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position
    var"offset" :: summary__pedestal_fits__linear__pressure_electron__offset
    var"pedestal_height" :: summary__pedestal_fits__linear__pressure_electron__pedestal_height
    var"pedestal_position" :: summary__pedestal_fits__linear__pressure_electron__pedestal_position
    var"pedestal_width" :: summary__pedestal_fits__linear__pressure_electron__pedestal_width
    var"separatrix" :: summary__pedestal_fits__linear__pressure_electron__separatrix
    _parent :: WeakRef
    function summary__pedestal_fits__linear__pressure_electron(var"d_dpsi_norm"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max(), var"d_dpsi_norm_max_position"=summary__pedestal_fits__linear__pressure_electron__d_dpsi_norm_max_position(), var"offset"=summary__pedestal_fits__linear__pressure_electron__offset(), var"pedestal_height"=summary__pedestal_fits__linear__pressure_electron__pedestal_height(), var"pedestal_position"=summary__pedestal_fits__linear__pressure_electron__pedestal_position(), var"pedestal_width"=summary__pedestal_fits__linear__pressure_electron__pedestal_width(), var"separatrix"=summary__pedestal_fits__linear__pressure_electron__separatrix(), _parent=WeakRef(missing))
        ids = new(var"d_dpsi_norm", var"d_dpsi_norm_max", var"d_dpsi_norm_max_position", var"offset", var"pedestal_height", var"pedestal_position", var"pedestal_width", var"separatrix", _parent)
        assign_expressions(ids)
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max_position, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
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
    var"d_dpsi_norm" :: summary__pedestal_fits__linear__n_e__d_dpsi_norm
    var"d_dpsi_norm_max" :: summary__pedestal_fits__linear__n_e__d_dpsi_norm_max
    var"offset" :: summary__pedestal_fits__linear__n_e__offset
    var"pedestal_height" :: summary__pedestal_fits__linear__n_e__pedestal_height
    var"pedestal_position" :: summary__pedestal_fits__linear__n_e__pedestal_position
    var"pedestal_width" :: summary__pedestal_fits__linear__n_e__pedestal_width
    var"separatrix" :: summary__pedestal_fits__linear__n_e__separatrix
    _parent :: WeakRef
    function summary__pedestal_fits__linear__n_e(var"d_dpsi_norm"=summary__pedestal_fits__linear__n_e__d_dpsi_norm(), var"d_dpsi_norm_max"=summary__pedestal_fits__linear__n_e__d_dpsi_norm_max(), var"offset"=summary__pedestal_fits__linear__n_e__offset(), var"pedestal_height"=summary__pedestal_fits__linear__n_e__pedestal_height(), var"pedestal_position"=summary__pedestal_fits__linear__n_e__pedestal_position(), var"pedestal_width"=summary__pedestal_fits__linear__n_e__pedestal_width(), var"separatrix"=summary__pedestal_fits__linear__n_e__separatrix(), _parent=WeakRef(missing))
        ids = new(var"d_dpsi_norm", var"d_dpsi_norm_max", var"offset", var"pedestal_height", var"pedestal_position", var"pedestal_width", var"separatrix", _parent)
        assign_expressions(ids)
        setfield!(ids.d_dpsi_norm, :_parent, WeakRef(ids))
        setfield!(ids.d_dpsi_norm_max, :_parent, WeakRef(ids))
        setfield!(ids.offset, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_height, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_position, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_width, :_parent, WeakRef(ids))
        setfield!(ids.separatrix, :_parent, WeakRef(ids))
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
    var"b_field_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_pedestal_top_hfs
    var"b_field_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_pedestal_top_lfs
    var"b_field_pol_pedestal_top_average" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_average
    var"b_field_pol_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs
    var"b_field_pol_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs
    var"b_field_tor_pedestal_top_hfs" :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs
    var"b_field_tor_pedestal_top_lfs" :: summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs
    var"beta_pol_pedestal_top_electron_average" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average
    var"beta_pol_pedestal_top_electron_hfs" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs
    var"beta_pol_pedestal_top_electron_lfs" :: summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs
    var"coulomb_factor_pedestal_top" :: summary__pedestal_fits__linear__coulomb_factor_pedestal_top
    var"energy_thermal_pedestal_electron" :: summary__pedestal_fits__linear__energy_thermal_pedestal_electron
    var"energy_thermal_pedestal_ion" :: summary__pedestal_fits__linear__energy_thermal_pedestal_ion
    var"n_e" :: summary__pedestal_fits__linear__n_e
    var"nustar_pedestal_top_electron" :: summary__pedestal_fits__linear__nustar_pedestal_top_electron
    var"parameters" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_electron" :: summary__pedestal_fits__linear__pressure_electron
    var"rhostar_pedestal_top_electron_hfs" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs
    var"rhostar_pedestal_top_electron_lfs" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs
    var"rhostar_pedestal_top_electron_magnetic_axis" :: summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis
    var"t_e" :: summary__pedestal_fits__linear__t_e
    var"volume_inside_pedestal" :: summary__pedestal_fits__linear__volume_inside_pedestal
    _parent :: WeakRef
    function summary__pedestal_fits__linear(var"b_field_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_pedestal_top_hfs(), var"b_field_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_pedestal_top_lfs(), var"b_field_pol_pedestal_top_average"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_average(), var"b_field_pol_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_hfs(), var"b_field_pol_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_pol_pedestal_top_lfs(), var"b_field_tor_pedestal_top_hfs"=summary__pedestal_fits__linear__b_field_tor_pedestal_top_hfs(), var"b_field_tor_pedestal_top_lfs"=summary__pedestal_fits__linear__b_field_tor_pedestal_top_lfs(), var"beta_pol_pedestal_top_electron_average"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_average(), var"beta_pol_pedestal_top_electron_hfs"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_hfs(), var"beta_pol_pedestal_top_electron_lfs"=summary__pedestal_fits__linear__beta_pol_pedestal_top_electron_lfs(), var"coulomb_factor_pedestal_top"=summary__pedestal_fits__linear__coulomb_factor_pedestal_top(), var"energy_thermal_pedestal_electron"=summary__pedestal_fits__linear__energy_thermal_pedestal_electron(), var"energy_thermal_pedestal_ion"=summary__pedestal_fits__linear__energy_thermal_pedestal_ion(), var"n_e"=summary__pedestal_fits__linear__n_e(), var"nustar_pedestal_top_electron"=summary__pedestal_fits__linear__nustar_pedestal_top_electron(), var"parameters"=missing, var"pressure_electron"=summary__pedestal_fits__linear__pressure_electron(), var"rhostar_pedestal_top_electron_hfs"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_hfs(), var"rhostar_pedestal_top_electron_lfs"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_lfs(), var"rhostar_pedestal_top_electron_magnetic_axis"=summary__pedestal_fits__linear__rhostar_pedestal_top_electron_magnetic_axis(), var"t_e"=summary__pedestal_fits__linear__t_e(), var"volume_inside_pedestal"=summary__pedestal_fits__linear__volume_inside_pedestal(), _parent=WeakRef(missing))
        ids = new(var"b_field_pedestal_top_hfs", var"b_field_pedestal_top_lfs", var"b_field_pol_pedestal_top_average", var"b_field_pol_pedestal_top_hfs", var"b_field_pol_pedestal_top_lfs", var"b_field_tor_pedestal_top_hfs", var"b_field_tor_pedestal_top_lfs", var"beta_pol_pedestal_top_electron_average", var"beta_pol_pedestal_top_electron_hfs", var"beta_pol_pedestal_top_electron_lfs", var"coulomb_factor_pedestal_top", var"energy_thermal_pedestal_electron", var"energy_thermal_pedestal_ion", var"n_e", var"nustar_pedestal_top_electron", var"parameters", var"pressure_electron", var"rhostar_pedestal_top_electron_hfs", var"rhostar_pedestal_top_electron_lfs", var"rhostar_pedestal_top_electron_magnetic_axis", var"t_e", var"volume_inside_pedestal", _parent)
        assign_expressions(ids)
        setfield!(ids.b_field_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_average, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_pol_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_hfs, :_parent, WeakRef(ids))
        setfield!(ids.b_field_tor_pedestal_top_lfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_average, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.coulomb_factor_pedestal_top, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_electron, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal_pedestal_ion, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.nustar_pedestal_top_electron, :_parent, WeakRef(ids))
        setfield!(ids.pressure_electron, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_hfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_lfs, :_parent, WeakRef(ids))
        setfield!(ids.rhostar_pedestal_top_electron_magnetic_axis, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.volume_inside_pedestal, :_parent, WeakRef(ids))
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__midplane(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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
    var"argon" :: summary__local__separatrix__velocity_tor__argon
    var"beryllium" :: summary__local__separatrix__velocity_tor__beryllium
    var"carbon" :: summary__local__separatrix__velocity_tor__carbon
    var"deuterium" :: summary__local__separatrix__velocity_tor__deuterium
    var"helium_3" :: summary__local__separatrix__velocity_tor__helium_3
    var"helium_4" :: summary__local__separatrix__velocity_tor__helium_4
    var"hydrogen" :: summary__local__separatrix__velocity_tor__hydrogen
    var"iron" :: summary__local__separatrix__velocity_tor__iron
    var"krypton" :: summary__local__separatrix__velocity_tor__krypton
    var"lithium" :: summary__local__separatrix__velocity_tor__lithium
    var"neon" :: summary__local__separatrix__velocity_tor__neon
    var"nitrogen" :: summary__local__separatrix__velocity_tor__nitrogen
    var"oxygen" :: summary__local__separatrix__velocity_tor__oxygen
    var"tritium" :: summary__local__separatrix__velocity_tor__tritium
    var"tungsten" :: summary__local__separatrix__velocity_tor__tungsten
    var"xenon" :: summary__local__separatrix__velocity_tor__xenon
    _parent :: WeakRef
    function summary__local__separatrix__velocity_tor(var"argon"=summary__local__separatrix__velocity_tor__argon(), var"beryllium"=summary__local__separatrix__velocity_tor__beryllium(), var"carbon"=summary__local__separatrix__velocity_tor__carbon(), var"deuterium"=summary__local__separatrix__velocity_tor__deuterium(), var"helium_3"=summary__local__separatrix__velocity_tor__helium_3(), var"helium_4"=summary__local__separatrix__velocity_tor__helium_4(), var"hydrogen"=summary__local__separatrix__velocity_tor__hydrogen(), var"iron"=summary__local__separatrix__velocity_tor__iron(), var"krypton"=summary__local__separatrix__velocity_tor__krypton(), var"lithium"=summary__local__separatrix__velocity_tor__lithium(), var"neon"=summary__local__separatrix__velocity_tor__neon(), var"nitrogen"=summary__local__separatrix__velocity_tor__nitrogen(), var"oxygen"=summary__local__separatrix__velocity_tor__oxygen(), var"tritium"=summary__local__separatrix__velocity_tor__tritium(), var"tungsten"=summary__local__separatrix__velocity_tor__tungsten(), var"xenon"=summary__local__separatrix__velocity_tor__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__separatrix__position(var"psi"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"rho_tor", var"rho_tor_norm", _parent)
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
    var"argon" :: summary__local__separatrix__n_i__argon
    var"beryllium" :: summary__local__separatrix__n_i__beryllium
    var"carbon" :: summary__local__separatrix__n_i__carbon
    var"deuterium" :: summary__local__separatrix__n_i__deuterium
    var"helium_3" :: summary__local__separatrix__n_i__helium_3
    var"helium_4" :: summary__local__separatrix__n_i__helium_4
    var"hydrogen" :: summary__local__separatrix__n_i__hydrogen
    var"iron" :: summary__local__separatrix__n_i__iron
    var"krypton" :: summary__local__separatrix__n_i__krypton
    var"lithium" :: summary__local__separatrix__n_i__lithium
    var"neon" :: summary__local__separatrix__n_i__neon
    var"nitrogen" :: summary__local__separatrix__n_i__nitrogen
    var"oxygen" :: summary__local__separatrix__n_i__oxygen
    var"tritium" :: summary__local__separatrix__n_i__tritium
    var"tungsten" :: summary__local__separatrix__n_i__tungsten
    var"xenon" :: summary__local__separatrix__n_i__xenon
    _parent :: WeakRef
    function summary__local__separatrix__n_i(var"argon"=summary__local__separatrix__n_i__argon(), var"beryllium"=summary__local__separatrix__n_i__beryllium(), var"carbon"=summary__local__separatrix__n_i__carbon(), var"deuterium"=summary__local__separatrix__n_i__deuterium(), var"helium_3"=summary__local__separatrix__n_i__helium_3(), var"helium_4"=summary__local__separatrix__n_i__helium_4(), var"hydrogen"=summary__local__separatrix__n_i__hydrogen(), var"iron"=summary__local__separatrix__n_i__iron(), var"krypton"=summary__local__separatrix__n_i__krypton(), var"lithium"=summary__local__separatrix__n_i__lithium(), var"neon"=summary__local__separatrix__n_i__neon(), var"nitrogen"=summary__local__separatrix__n_i__nitrogen(), var"oxygen"=summary__local__separatrix__n_i__oxygen(), var"tritium"=summary__local__separatrix__n_i__tritium(), var"tungsten"=summary__local__separatrix__n_i__tungsten(), var"xenon"=summary__local__separatrix__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"e_field_parallel" :: summary__local__separatrix__e_field_parallel
    var"magnetic_shear" :: summary__local__separatrix__magnetic_shear
    var"momentum_tor" :: summary__local__separatrix__momentum_tor
    var"n_e" :: summary__local__separatrix__n_e
    var"n_i" :: summary__local__separatrix__n_i
    var"n_i_total" :: summary__local__separatrix__n_i_total
    var"position" :: summary__local__separatrix__position
    var"q" :: summary__local__separatrix__q
    var"t_e" :: summary__local__separatrix__t_e
    var"t_i_average" :: summary__local__separatrix__t_i_average
    var"velocity_tor" :: summary__local__separatrix__velocity_tor
    var"zeff" :: summary__local__separatrix__zeff
    _parent :: WeakRef
    function summary__local__separatrix(var"e_field_parallel"=summary__local__separatrix__e_field_parallel(), var"magnetic_shear"=summary__local__separatrix__magnetic_shear(), var"momentum_tor"=summary__local__separatrix__momentum_tor(), var"n_e"=summary__local__separatrix__n_e(), var"n_i"=summary__local__separatrix__n_i(), var"n_i_total"=summary__local__separatrix__n_i_total(), var"position"=summary__local__separatrix__position(), var"q"=summary__local__separatrix__q(), var"t_e"=summary__local__separatrix__t_e(), var"t_i_average"=summary__local__separatrix__t_i_average(), var"velocity_tor"=summary__local__separatrix__velocity_tor(), var"zeff"=summary__local__separatrix__zeff(), _parent=WeakRef(missing))
        ids = new(var"e_field_parallel", var"magnetic_shear", var"momentum_tor", var"n_e", var"n_i", var"n_i_total", var"position", var"q", var"t_e", var"t_i_average", var"velocity_tor", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
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
    var"effective_helical_ripple" :: summary__local__r_eff_norm_2_3__effective_helical_ripple
    var"iota" :: summary__local__r_eff_norm_2_3__iota
    var"plateau_factor" :: summary__local__r_eff_norm_2_3__plateau_factor
    _parent :: WeakRef
    function summary__local__r_eff_norm_2_3(var"effective_helical_ripple"=summary__local__r_eff_norm_2_3__effective_helical_ripple(), var"iota"=summary__local__r_eff_norm_2_3__iota(), var"plateau_factor"=summary__local__r_eff_norm_2_3__plateau_factor(), _parent=WeakRef(missing))
        ids = new(var"effective_helical_ripple", var"iota", var"plateau_factor", _parent)
        assign_expressions(ids)
        setfield!(ids.effective_helical_ripple, :_parent, WeakRef(ids))
        setfield!(ids.iota, :_parent, WeakRef(ids))
        setfield!(ids.plateau_factor, :_parent, WeakRef(ids))
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
    var"argon" :: summary__local__pedestal__velocity_tor__argon
    var"beryllium" :: summary__local__pedestal__velocity_tor__beryllium
    var"carbon" :: summary__local__pedestal__velocity_tor__carbon
    var"deuterium" :: summary__local__pedestal__velocity_tor__deuterium
    var"helium_3" :: summary__local__pedestal__velocity_tor__helium_3
    var"helium_4" :: summary__local__pedestal__velocity_tor__helium_4
    var"hydrogen" :: summary__local__pedestal__velocity_tor__hydrogen
    var"iron" :: summary__local__pedestal__velocity_tor__iron
    var"krypton" :: summary__local__pedestal__velocity_tor__krypton
    var"lithium" :: summary__local__pedestal__velocity_tor__lithium
    var"neon" :: summary__local__pedestal__velocity_tor__neon
    var"nitrogen" :: summary__local__pedestal__velocity_tor__nitrogen
    var"oxygen" :: summary__local__pedestal__velocity_tor__oxygen
    var"tritium" :: summary__local__pedestal__velocity_tor__tritium
    var"tungsten" :: summary__local__pedestal__velocity_tor__tungsten
    var"xenon" :: summary__local__pedestal__velocity_tor__xenon
    _parent :: WeakRef
    function summary__local__pedestal__velocity_tor(var"argon"=summary__local__pedestal__velocity_tor__argon(), var"beryllium"=summary__local__pedestal__velocity_tor__beryllium(), var"carbon"=summary__local__pedestal__velocity_tor__carbon(), var"deuterium"=summary__local__pedestal__velocity_tor__deuterium(), var"helium_3"=summary__local__pedestal__velocity_tor__helium_3(), var"helium_4"=summary__local__pedestal__velocity_tor__helium_4(), var"hydrogen"=summary__local__pedestal__velocity_tor__hydrogen(), var"iron"=summary__local__pedestal__velocity_tor__iron(), var"krypton"=summary__local__pedestal__velocity_tor__krypton(), var"lithium"=summary__local__pedestal__velocity_tor__lithium(), var"neon"=summary__local__pedestal__velocity_tor__neon(), var"nitrogen"=summary__local__pedestal__velocity_tor__nitrogen(), var"oxygen"=summary__local__pedestal__velocity_tor__oxygen(), var"tritium"=summary__local__pedestal__velocity_tor__tritium(), var"tungsten"=summary__local__pedestal__velocity_tor__tungsten(), var"xenon"=summary__local__pedestal__velocity_tor__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__pedestal__position(var"psi"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"rho_tor", var"rho_tor_norm", _parent)
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
    var"argon" :: summary__local__pedestal__n_i__argon
    var"beryllium" :: summary__local__pedestal__n_i__beryllium
    var"carbon" :: summary__local__pedestal__n_i__carbon
    var"deuterium" :: summary__local__pedestal__n_i__deuterium
    var"helium_3" :: summary__local__pedestal__n_i__helium_3
    var"helium_4" :: summary__local__pedestal__n_i__helium_4
    var"hydrogen" :: summary__local__pedestal__n_i__hydrogen
    var"iron" :: summary__local__pedestal__n_i__iron
    var"krypton" :: summary__local__pedestal__n_i__krypton
    var"lithium" :: summary__local__pedestal__n_i__lithium
    var"neon" :: summary__local__pedestal__n_i__neon
    var"nitrogen" :: summary__local__pedestal__n_i__nitrogen
    var"oxygen" :: summary__local__pedestal__n_i__oxygen
    var"tritium" :: summary__local__pedestal__n_i__tritium
    var"tungsten" :: summary__local__pedestal__n_i__tungsten
    var"xenon" :: summary__local__pedestal__n_i__xenon
    _parent :: WeakRef
    function summary__local__pedestal__n_i(var"argon"=summary__local__pedestal__n_i__argon(), var"beryllium"=summary__local__pedestal__n_i__beryllium(), var"carbon"=summary__local__pedestal__n_i__carbon(), var"deuterium"=summary__local__pedestal__n_i__deuterium(), var"helium_3"=summary__local__pedestal__n_i__helium_3(), var"helium_4"=summary__local__pedestal__n_i__helium_4(), var"hydrogen"=summary__local__pedestal__n_i__hydrogen(), var"iron"=summary__local__pedestal__n_i__iron(), var"krypton"=summary__local__pedestal__n_i__krypton(), var"lithium"=summary__local__pedestal__n_i__lithium(), var"neon"=summary__local__pedestal__n_i__neon(), var"nitrogen"=summary__local__pedestal__n_i__nitrogen(), var"oxygen"=summary__local__pedestal__n_i__oxygen(), var"tritium"=summary__local__pedestal__n_i__tritium(), var"tungsten"=summary__local__pedestal__n_i__tungsten(), var"xenon"=summary__local__pedestal__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"e_field_parallel" :: summary__local__pedestal__e_field_parallel
    var"magnetic_shear" :: summary__local__pedestal__magnetic_shear
    var"momentum_tor" :: summary__local__pedestal__momentum_tor
    var"n_e" :: summary__local__pedestal__n_e
    var"n_i" :: summary__local__pedestal__n_i
    var"n_i_total" :: summary__local__pedestal__n_i_total
    var"position" :: summary__local__pedestal__position
    var"q" :: summary__local__pedestal__q
    var"t_e" :: summary__local__pedestal__t_e
    var"t_i_average" :: summary__local__pedestal__t_i_average
    var"velocity_tor" :: summary__local__pedestal__velocity_tor
    var"zeff" :: summary__local__pedestal__zeff
    _parent :: WeakRef
    function summary__local__pedestal(var"e_field_parallel"=summary__local__pedestal__e_field_parallel(), var"magnetic_shear"=summary__local__pedestal__magnetic_shear(), var"momentum_tor"=summary__local__pedestal__momentum_tor(), var"n_e"=summary__local__pedestal__n_e(), var"n_i"=summary__local__pedestal__n_i(), var"n_i_total"=summary__local__pedestal__n_i_total(), var"position"=summary__local__pedestal__position(), var"q"=summary__local__pedestal__q(), var"t_e"=summary__local__pedestal__t_e(), var"t_i_average"=summary__local__pedestal__t_i_average(), var"velocity_tor"=summary__local__pedestal__velocity_tor(), var"zeff"=summary__local__pedestal__zeff(), _parent=WeakRef(missing))
        ids = new(var"e_field_parallel", var"magnetic_shear", var"momentum_tor", var"n_e", var"n_i", var"n_i_total", var"position", var"q", var"t_e", var"t_i_average", var"velocity_tor", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
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
    var"argon" :: summary__local__magnetic_axis__velocity_tor__argon
    var"beryllium" :: summary__local__magnetic_axis__velocity_tor__beryllium
    var"carbon" :: summary__local__magnetic_axis__velocity_tor__carbon
    var"deuterium" :: summary__local__magnetic_axis__velocity_tor__deuterium
    var"helium_3" :: summary__local__magnetic_axis__velocity_tor__helium_3
    var"helium_4" :: summary__local__magnetic_axis__velocity_tor__helium_4
    var"hydrogen" :: summary__local__magnetic_axis__velocity_tor__hydrogen
    var"iron" :: summary__local__magnetic_axis__velocity_tor__iron
    var"krypton" :: summary__local__magnetic_axis__velocity_tor__krypton
    var"lithium" :: summary__local__magnetic_axis__velocity_tor__lithium
    var"neon" :: summary__local__magnetic_axis__velocity_tor__neon
    var"nitrogen" :: summary__local__magnetic_axis__velocity_tor__nitrogen
    var"oxygen" :: summary__local__magnetic_axis__velocity_tor__oxygen
    var"tritium" :: summary__local__magnetic_axis__velocity_tor__tritium
    var"tungsten" :: summary__local__magnetic_axis__velocity_tor__tungsten
    var"xenon" :: summary__local__magnetic_axis__velocity_tor__xenon
    _parent :: WeakRef
    function summary__local__magnetic_axis__velocity_tor(var"argon"=summary__local__magnetic_axis__velocity_tor__argon(), var"beryllium"=summary__local__magnetic_axis__velocity_tor__beryllium(), var"carbon"=summary__local__magnetic_axis__velocity_tor__carbon(), var"deuterium"=summary__local__magnetic_axis__velocity_tor__deuterium(), var"helium_3"=summary__local__magnetic_axis__velocity_tor__helium_3(), var"helium_4"=summary__local__magnetic_axis__velocity_tor__helium_4(), var"hydrogen"=summary__local__magnetic_axis__velocity_tor__hydrogen(), var"iron"=summary__local__magnetic_axis__velocity_tor__iron(), var"krypton"=summary__local__magnetic_axis__velocity_tor__krypton(), var"lithium"=summary__local__magnetic_axis__velocity_tor__lithium(), var"neon"=summary__local__magnetic_axis__velocity_tor__neon(), var"nitrogen"=summary__local__magnetic_axis__velocity_tor__nitrogen(), var"oxygen"=summary__local__magnetic_axis__velocity_tor__oxygen(), var"tritium"=summary__local__magnetic_axis__velocity_tor__tritium(), var"tungsten"=summary__local__magnetic_axis__velocity_tor__tungsten(), var"xenon"=summary__local__magnetic_axis__velocity_tor__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__magnetic_axis__position(var"psi"=missing, var"r"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"r", var"rho_tor", var"rho_tor_norm", var"z", _parent)
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
    var"argon" :: summary__local__magnetic_axis__n_i__argon
    var"beryllium" :: summary__local__magnetic_axis__n_i__beryllium
    var"carbon" :: summary__local__magnetic_axis__n_i__carbon
    var"deuterium" :: summary__local__magnetic_axis__n_i__deuterium
    var"helium_3" :: summary__local__magnetic_axis__n_i__helium_3
    var"helium_4" :: summary__local__magnetic_axis__n_i__helium_4
    var"hydrogen" :: summary__local__magnetic_axis__n_i__hydrogen
    var"iron" :: summary__local__magnetic_axis__n_i__iron
    var"krypton" :: summary__local__magnetic_axis__n_i__krypton
    var"lithium" :: summary__local__magnetic_axis__n_i__lithium
    var"neon" :: summary__local__magnetic_axis__n_i__neon
    var"nitrogen" :: summary__local__magnetic_axis__n_i__nitrogen
    var"oxygen" :: summary__local__magnetic_axis__n_i__oxygen
    var"tritium" :: summary__local__magnetic_axis__n_i__tritium
    var"tungsten" :: summary__local__magnetic_axis__n_i__tungsten
    var"xenon" :: summary__local__magnetic_axis__n_i__xenon
    _parent :: WeakRef
    function summary__local__magnetic_axis__n_i(var"argon"=summary__local__magnetic_axis__n_i__argon(), var"beryllium"=summary__local__magnetic_axis__n_i__beryllium(), var"carbon"=summary__local__magnetic_axis__n_i__carbon(), var"deuterium"=summary__local__magnetic_axis__n_i__deuterium(), var"helium_3"=summary__local__magnetic_axis__n_i__helium_3(), var"helium_4"=summary__local__magnetic_axis__n_i__helium_4(), var"hydrogen"=summary__local__magnetic_axis__n_i__hydrogen(), var"iron"=summary__local__magnetic_axis__n_i__iron(), var"krypton"=summary__local__magnetic_axis__n_i__krypton(), var"lithium"=summary__local__magnetic_axis__n_i__lithium(), var"neon"=summary__local__magnetic_axis__n_i__neon(), var"nitrogen"=summary__local__magnetic_axis__n_i__nitrogen(), var"oxygen"=summary__local__magnetic_axis__n_i__oxygen(), var"tritium"=summary__local__magnetic_axis__n_i__tritium(), var"tungsten"=summary__local__magnetic_axis__n_i__tungsten(), var"xenon"=summary__local__magnetic_axis__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"b_field" :: summary__local__magnetic_axis__b_field
    var"e_field_parallel" :: summary__local__magnetic_axis__e_field_parallel
    var"magnetic_shear" :: summary__local__magnetic_axis__magnetic_shear
    var"momentum_tor" :: summary__local__magnetic_axis__momentum_tor
    var"n_e" :: summary__local__magnetic_axis__n_e
    var"n_i" :: summary__local__magnetic_axis__n_i
    var"n_i_total" :: summary__local__magnetic_axis__n_i_total
    var"position" :: summary__local__magnetic_axis__position
    var"q" :: summary__local__magnetic_axis__q
    var"t_e" :: summary__local__magnetic_axis__t_e
    var"t_i_average" :: summary__local__magnetic_axis__t_i_average
    var"velocity_tor" :: summary__local__magnetic_axis__velocity_tor
    var"zeff" :: summary__local__magnetic_axis__zeff
    _parent :: WeakRef
    function summary__local__magnetic_axis(var"b_field"=summary__local__magnetic_axis__b_field(), var"e_field_parallel"=summary__local__magnetic_axis__e_field_parallel(), var"magnetic_shear"=summary__local__magnetic_axis__magnetic_shear(), var"momentum_tor"=summary__local__magnetic_axis__momentum_tor(), var"n_e"=summary__local__magnetic_axis__n_e(), var"n_i"=summary__local__magnetic_axis__n_i(), var"n_i_total"=summary__local__magnetic_axis__n_i_total(), var"position"=summary__local__magnetic_axis__position(), var"q"=summary__local__magnetic_axis__q(), var"t_e"=summary__local__magnetic_axis__t_e(), var"t_i_average"=summary__local__magnetic_axis__t_i_average(), var"velocity_tor"=summary__local__magnetic_axis__velocity_tor(), var"zeff"=summary__local__magnetic_axis__zeff(), _parent=WeakRef(missing))
        ids = new(var"b_field", var"e_field_parallel", var"magnetic_shear", var"momentum_tor", var"n_e", var"n_i", var"n_i_total", var"position", var"q", var"t_e", var"t_i_average", var"velocity_tor", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.b_field, :_parent, WeakRef(ids))
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
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
    var"argon" :: summary__local__limiter__n_i__argon
    var"beryllium" :: summary__local__limiter__n_i__beryllium
    var"carbon" :: summary__local__limiter__n_i__carbon
    var"deuterium" :: summary__local__limiter__n_i__deuterium
    var"helium_3" :: summary__local__limiter__n_i__helium_3
    var"helium_4" :: summary__local__limiter__n_i__helium_4
    var"hydrogen" :: summary__local__limiter__n_i__hydrogen
    var"iron" :: summary__local__limiter__n_i__iron
    var"krypton" :: summary__local__limiter__n_i__krypton
    var"lithium" :: summary__local__limiter__n_i__lithium
    var"neon" :: summary__local__limiter__n_i__neon
    var"nitrogen" :: summary__local__limiter__n_i__nitrogen
    var"oxygen" :: summary__local__limiter__n_i__oxygen
    var"tritium" :: summary__local__limiter__n_i__tritium
    var"tungsten" :: summary__local__limiter__n_i__tungsten
    var"xenon" :: summary__local__limiter__n_i__xenon
    _parent :: WeakRef
    function summary__local__limiter__n_i(var"argon"=summary__local__limiter__n_i__argon(), var"beryllium"=summary__local__limiter__n_i__beryllium(), var"carbon"=summary__local__limiter__n_i__carbon(), var"deuterium"=summary__local__limiter__n_i__deuterium(), var"helium_3"=summary__local__limiter__n_i__helium_3(), var"helium_4"=summary__local__limiter__n_i__helium_4(), var"hydrogen"=summary__local__limiter__n_i__hydrogen(), var"iron"=summary__local__limiter__n_i__iron(), var"krypton"=summary__local__limiter__n_i__krypton(), var"lithium"=summary__local__limiter__n_i__lithium(), var"neon"=summary__local__limiter__n_i__neon(), var"nitrogen"=summary__local__limiter__n_i__nitrogen(), var"oxygen"=summary__local__limiter__n_i__oxygen(), var"tritium"=summary__local__limiter__n_i__tritium(), var"tungsten"=summary__local__limiter__n_i__tungsten(), var"xenon"=summary__local__limiter__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"n_e" :: summary__local__limiter__n_e
    var"n_i" :: summary__local__limiter__n_i
    var"n_i_total" :: summary__local__limiter__n_i_total
    var"name" :: summary__local__limiter__name
    var"power_flux_peak" :: summary__local__limiter__power_flux_peak
    var"t_e" :: summary__local__limiter__t_e
    var"t_i_average" :: summary__local__limiter__t_i_average
    var"zeff" :: summary__local__limiter__zeff
    _parent :: WeakRef
    function summary__local__limiter(var"flux_expansion"=summary__local__limiter__flux_expansion(), var"n_e"=summary__local__limiter__n_e(), var"n_i"=summary__local__limiter__n_i(), var"n_i_total"=summary__local__limiter__n_i_total(), var"name"=summary__local__limiter__name(), var"power_flux_peak"=summary__local__limiter__power_flux_peak(), var"t_e"=summary__local__limiter__t_e(), var"t_i_average"=summary__local__limiter__t_i_average(), var"zeff"=summary__local__limiter__zeff(), _parent=WeakRef(missing))
        ids = new(var"flux_expansion", var"n_e", var"n_i", var"n_i_total", var"name", var"power_flux_peak", var"t_e", var"t_i_average", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.flux_expansion, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.name, :_parent, WeakRef(ids))
        setfield!(ids.power_flux_peak, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
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
    var"argon" :: summary__local__itb__velocity_tor__argon
    var"beryllium" :: summary__local__itb__velocity_tor__beryllium
    var"carbon" :: summary__local__itb__velocity_tor__carbon
    var"deuterium" :: summary__local__itb__velocity_tor__deuterium
    var"helium_3" :: summary__local__itb__velocity_tor__helium_3
    var"helium_4" :: summary__local__itb__velocity_tor__helium_4
    var"hydrogen" :: summary__local__itb__velocity_tor__hydrogen
    var"iron" :: summary__local__itb__velocity_tor__iron
    var"krypton" :: summary__local__itb__velocity_tor__krypton
    var"lithium" :: summary__local__itb__velocity_tor__lithium
    var"neon" :: summary__local__itb__velocity_tor__neon
    var"nitrogen" :: summary__local__itb__velocity_tor__nitrogen
    var"oxygen" :: summary__local__itb__velocity_tor__oxygen
    var"tritium" :: summary__local__itb__velocity_tor__tritium
    var"tungsten" :: summary__local__itb__velocity_tor__tungsten
    var"xenon" :: summary__local__itb__velocity_tor__xenon
    _parent :: WeakRef
    function summary__local__itb__velocity_tor(var"argon"=summary__local__itb__velocity_tor__argon(), var"beryllium"=summary__local__itb__velocity_tor__beryllium(), var"carbon"=summary__local__itb__velocity_tor__carbon(), var"deuterium"=summary__local__itb__velocity_tor__deuterium(), var"helium_3"=summary__local__itb__velocity_tor__helium_3(), var"helium_4"=summary__local__itb__velocity_tor__helium_4(), var"hydrogen"=summary__local__itb__velocity_tor__hydrogen(), var"iron"=summary__local__itb__velocity_tor__iron(), var"krypton"=summary__local__itb__velocity_tor__krypton(), var"lithium"=summary__local__itb__velocity_tor__lithium(), var"neon"=summary__local__itb__velocity_tor__neon(), var"nitrogen"=summary__local__itb__velocity_tor__nitrogen(), var"oxygen"=summary__local__itb__velocity_tor__oxygen(), var"tritium"=summary__local__itb__velocity_tor__tritium(), var"tungsten"=summary__local__itb__velocity_tor__tungsten(), var"xenon"=summary__local__itb__velocity_tor__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function summary__local__itb__position(var"psi"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, _parent=WeakRef(missing))
        ids = new(var"psi", var"rho_tor", var"rho_tor_norm", _parent)
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
    var"argon" :: summary__local__itb__n_i__argon
    var"beryllium" :: summary__local__itb__n_i__beryllium
    var"carbon" :: summary__local__itb__n_i__carbon
    var"deuterium" :: summary__local__itb__n_i__deuterium
    var"helium_3" :: summary__local__itb__n_i__helium_3
    var"helium_4" :: summary__local__itb__n_i__helium_4
    var"hydrogen" :: summary__local__itb__n_i__hydrogen
    var"iron" :: summary__local__itb__n_i__iron
    var"krypton" :: summary__local__itb__n_i__krypton
    var"lithium" :: summary__local__itb__n_i__lithium
    var"neon" :: summary__local__itb__n_i__neon
    var"nitrogen" :: summary__local__itb__n_i__nitrogen
    var"oxygen" :: summary__local__itb__n_i__oxygen
    var"tritium" :: summary__local__itb__n_i__tritium
    var"tungsten" :: summary__local__itb__n_i__tungsten
    var"xenon" :: summary__local__itb__n_i__xenon
    _parent :: WeakRef
    function summary__local__itb__n_i(var"argon"=summary__local__itb__n_i__argon(), var"beryllium"=summary__local__itb__n_i__beryllium(), var"carbon"=summary__local__itb__n_i__carbon(), var"deuterium"=summary__local__itb__n_i__deuterium(), var"helium_3"=summary__local__itb__n_i__helium_3(), var"helium_4"=summary__local__itb__n_i__helium_4(), var"hydrogen"=summary__local__itb__n_i__hydrogen(), var"iron"=summary__local__itb__n_i__iron(), var"krypton"=summary__local__itb__n_i__krypton(), var"lithium"=summary__local__itb__n_i__lithium(), var"neon"=summary__local__itb__n_i__neon(), var"nitrogen"=summary__local__itb__n_i__nitrogen(), var"oxygen"=summary__local__itb__n_i__oxygen(), var"tritium"=summary__local__itb__n_i__tritium(), var"tungsten"=summary__local__itb__n_i__tungsten(), var"xenon"=summary__local__itb__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"e_field_parallel" :: summary__local__itb__e_field_parallel
    var"magnetic_shear" :: summary__local__itb__magnetic_shear
    var"momentum_tor" :: summary__local__itb__momentum_tor
    var"n_e" :: summary__local__itb__n_e
    var"n_i" :: summary__local__itb__n_i
    var"n_i_total" :: summary__local__itb__n_i_total
    var"position" :: summary__local__itb__position
    var"q" :: summary__local__itb__q
    var"t_e" :: summary__local__itb__t_e
    var"t_i_average" :: summary__local__itb__t_i_average
    var"velocity_tor" :: summary__local__itb__velocity_tor
    var"zeff" :: summary__local__itb__zeff
    _parent :: WeakRef
    function summary__local__itb(var"e_field_parallel"=summary__local__itb__e_field_parallel(), var"magnetic_shear"=summary__local__itb__magnetic_shear(), var"momentum_tor"=summary__local__itb__momentum_tor(), var"n_e"=summary__local__itb__n_e(), var"n_i"=summary__local__itb__n_i(), var"n_i_total"=summary__local__itb__n_i_total(), var"position"=summary__local__itb__position(), var"q"=summary__local__itb__q(), var"t_e"=summary__local__itb__t_e(), var"t_i_average"=summary__local__itb__t_i_average(), var"velocity_tor"=summary__local__itb__velocity_tor(), var"zeff"=summary__local__itb__zeff(), _parent=WeakRef(missing))
        ids = new(var"e_field_parallel", var"magnetic_shear", var"momentum_tor", var"n_e", var"n_i", var"n_i_total", var"position", var"q", var"t_e", var"t_i_average", var"velocity_tor", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.e_field_parallel, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear, :_parent, WeakRef(ids))
        setfield!(ids.momentum_tor, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.velocity_tor, :_parent, WeakRef(ids))
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
    var"argon" :: summary__local__divertor_plate___n_i__argon
    var"beryllium" :: summary__local__divertor_plate___n_i__beryllium
    var"carbon" :: summary__local__divertor_plate___n_i__carbon
    var"deuterium" :: summary__local__divertor_plate___n_i__deuterium
    var"helium_3" :: summary__local__divertor_plate___n_i__helium_3
    var"helium_4" :: summary__local__divertor_plate___n_i__helium_4
    var"hydrogen" :: summary__local__divertor_plate___n_i__hydrogen
    var"iron" :: summary__local__divertor_plate___n_i__iron
    var"krypton" :: summary__local__divertor_plate___n_i__krypton
    var"lithium" :: summary__local__divertor_plate___n_i__lithium
    var"neon" :: summary__local__divertor_plate___n_i__neon
    var"nitrogen" :: summary__local__divertor_plate___n_i__nitrogen
    var"oxygen" :: summary__local__divertor_plate___n_i__oxygen
    var"tritium" :: summary__local__divertor_plate___n_i__tritium
    var"tungsten" :: summary__local__divertor_plate___n_i__tungsten
    var"xenon" :: summary__local__divertor_plate___n_i__xenon
    _parent :: WeakRef
    function summary__local__divertor_plate___n_i(var"argon"=summary__local__divertor_plate___n_i__argon(), var"beryllium"=summary__local__divertor_plate___n_i__beryllium(), var"carbon"=summary__local__divertor_plate___n_i__carbon(), var"deuterium"=summary__local__divertor_plate___n_i__deuterium(), var"helium_3"=summary__local__divertor_plate___n_i__helium_3(), var"helium_4"=summary__local__divertor_plate___n_i__helium_4(), var"hydrogen"=summary__local__divertor_plate___n_i__hydrogen(), var"iron"=summary__local__divertor_plate___n_i__iron(), var"krypton"=summary__local__divertor_plate___n_i__krypton(), var"lithium"=summary__local__divertor_plate___n_i__lithium(), var"neon"=summary__local__divertor_plate___n_i__neon(), var"nitrogen"=summary__local__divertor_plate___n_i__nitrogen(), var"oxygen"=summary__local__divertor_plate___n_i__oxygen(), var"tritium"=summary__local__divertor_plate___n_i__tritium(), var"tungsten"=summary__local__divertor_plate___n_i__tungsten(), var"xenon"=summary__local__divertor_plate___n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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

mutable struct summary__local__divertor_plate <: IDSvectorStaticElement
    var"flux_expansion" :: summary__local__divertor_plate___flux_expansion
    var"n_e" :: summary__local__divertor_plate___n_e
    var"n_i" :: summary__local__divertor_plate___n_i
    var"n_i_total" :: summary__local__divertor_plate___n_i_total
    var"name" :: summary__local__divertor_plate___name
    var"power_flux_peak" :: summary__local__divertor_plate___power_flux_peak
    var"t_e" :: summary__local__divertor_plate___t_e
    var"t_i_average" :: summary__local__divertor_plate___t_i_average
    var"zeff" :: summary__local__divertor_plate___zeff
    _parent :: WeakRef
    function summary__local__divertor_plate(var"flux_expansion"=summary__local__divertor_plate___flux_expansion(), var"n_e"=summary__local__divertor_plate___n_e(), var"n_i"=summary__local__divertor_plate___n_i(), var"n_i_total"=summary__local__divertor_plate___n_i_total(), var"name"=summary__local__divertor_plate___name(), var"power_flux_peak"=summary__local__divertor_plate___power_flux_peak(), var"t_e"=summary__local__divertor_plate___t_e(), var"t_i_average"=summary__local__divertor_plate___t_i_average(), var"zeff"=summary__local__divertor_plate___zeff(), _parent=WeakRef(missing))
        ids = new(var"flux_expansion", var"n_e", var"n_i", var"n_i_total", var"name", var"power_flux_peak", var"t_e", var"t_i_average", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.flux_expansion, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.name, :_parent, WeakRef(ids))
        setfield!(ids.power_flux_peak, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__local <: IDS
    var"divertor_plate" :: IDSvector{T} where {T<:summary__local__divertor_plate}
    var"itb" :: summary__local__itb
    var"limiter" :: summary__local__limiter
    var"magnetic_axis" :: summary__local__magnetic_axis
    var"pedestal" :: summary__local__pedestal
    var"r_eff_norm_2_3" :: summary__local__r_eff_norm_2_3
    var"separatrix" :: summary__local__separatrix
    _parent :: WeakRef
    function summary__local(var"divertor_plate"=IDSvector(summary__local__divertor_plate[]), var"itb"=summary__local__itb(), var"limiter"=summary__local__limiter(), var"magnetic_axis"=summary__local__magnetic_axis(), var"pedestal"=summary__local__pedestal(), var"r_eff_norm_2_3"=summary__local__r_eff_norm_2_3(), var"separatrix"=summary__local__separatrix(), _parent=WeakRef(missing))
        ids = new(var"divertor_plate", var"itb", var"limiter", var"magnetic_axis", var"pedestal", var"r_eff_norm_2_3", var"separatrix", _parent)
        assign_expressions(ids)
        setfield!(ids.divertor_plate, :_parent, WeakRef(ids))
        setfield!(ids.itb, :_parent, WeakRef(ids))
        setfield!(ids.limiter, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis, :_parent, WeakRef(ids))
        setfield!(ids.pedestal, :_parent, WeakRef(ids))
        setfield!(ids.r_eff_norm_2_3, :_parent, WeakRef(ids))
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
    var"argon" :: summary__line_average__n_i__argon
    var"beryllium" :: summary__line_average__n_i__beryllium
    var"carbon" :: summary__line_average__n_i__carbon
    var"deuterium" :: summary__line_average__n_i__deuterium
    var"helium_3" :: summary__line_average__n_i__helium_3
    var"helium_4" :: summary__line_average__n_i__helium_4
    var"hydrogen" :: summary__line_average__n_i__hydrogen
    var"iron" :: summary__line_average__n_i__iron
    var"krypton" :: summary__line_average__n_i__krypton
    var"lithium" :: summary__line_average__n_i__lithium
    var"neon" :: summary__line_average__n_i__neon
    var"nitrogen" :: summary__line_average__n_i__nitrogen
    var"oxygen" :: summary__line_average__n_i__oxygen
    var"tritium" :: summary__line_average__n_i__tritium
    var"tungsten" :: summary__line_average__n_i__tungsten
    var"xenon" :: summary__line_average__n_i__xenon
    _parent :: WeakRef
    function summary__line_average__n_i(var"argon"=summary__line_average__n_i__argon(), var"beryllium"=summary__line_average__n_i__beryllium(), var"carbon"=summary__line_average__n_i__carbon(), var"deuterium"=summary__line_average__n_i__deuterium(), var"helium_3"=summary__line_average__n_i__helium_3(), var"helium_4"=summary__line_average__n_i__helium_4(), var"hydrogen"=summary__line_average__n_i__hydrogen(), var"iron"=summary__line_average__n_i__iron(), var"krypton"=summary__line_average__n_i__krypton(), var"lithium"=summary__line_average__n_i__lithium(), var"neon"=summary__line_average__n_i__neon(), var"nitrogen"=summary__line_average__n_i__nitrogen(), var"oxygen"=summary__line_average__n_i__oxygen(), var"tritium"=summary__line_average__n_i__tritium(), var"tungsten"=summary__line_average__n_i__tungsten(), var"xenon"=summary__line_average__n_i__xenon(), _parent=WeakRef(missing))
        ids = new(var"argon", var"beryllium", var"carbon", var"deuterium", var"helium_3", var"helium_4", var"hydrogen", var"iron", var"krypton", var"lithium", var"neon", var"nitrogen", var"oxygen", var"tritium", var"tungsten", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.iron, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.tungsten, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"dn_e_dt" :: summary__line_average__dn_e_dt
    var"isotope_fraction_hydrogen" :: summary__line_average__isotope_fraction_hydrogen
    var"meff_hydrogenic" :: summary__line_average__meff_hydrogenic
    var"n_e" :: summary__line_average__n_e
    var"n_i" :: summary__line_average__n_i
    var"n_i_total" :: summary__line_average__n_i_total
    var"t_e" :: summary__line_average__t_e
    var"t_i_average" :: summary__line_average__t_i_average
    var"zeff" :: summary__line_average__zeff
    _parent :: WeakRef
    function summary__line_average(var"dn_e_dt"=summary__line_average__dn_e_dt(), var"isotope_fraction_hydrogen"=summary__line_average__isotope_fraction_hydrogen(), var"meff_hydrogenic"=summary__line_average__meff_hydrogenic(), var"n_e"=summary__line_average__n_e(), var"n_i"=summary__line_average__n_i(), var"n_i_total"=summary__line_average__n_i_total(), var"t_e"=summary__line_average__t_e(), var"t_i_average"=summary__line_average__t_i_average(), var"zeff"=summary__line_average__zeff(), _parent=WeakRef(missing))
        ids = new(var"dn_e_dt", var"isotope_fraction_hydrogen", var"meff_hydrogenic", var"n_e", var"n_i", var"n_i_total", var"t_e", var"t_i_average", var"zeff", _parent)
        assign_expressions(ids)
        setfield!(ids.dn_e_dt, :_parent, WeakRef(ids))
        setfield!(ids.isotope_fraction_hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.meff_hydrogenic, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_i, :_parent, WeakRef(ids))
        setfield!(ids.n_i_total, :_parent, WeakRef(ids))
        setfield!(ids.t_e, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average, :_parent, WeakRef(ids))
        setfield!(ids.zeff, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__limiter__material <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function summary__limiter__material(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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
    var"a" :: summary__heating_current_drive__nbi___species__a
    var"label" :: summary__heating_current_drive__nbi___species__label
    var"z_n" :: summary__heating_current_drive__nbi___species__z_n
    _parent :: WeakRef
    function summary__heating_current_drive__nbi___species(var"a"=summary__heating_current_drive__nbi___species__a(), var"label"=summary__heating_current_drive__nbi___species__label(), var"z_n"=summary__heating_current_drive__nbi___species__z_n(), _parent=WeakRef(missing))
        ids = new(var"a", var"label", var"z_n", _parent)
        assign_expressions(ids)
        setfield!(ids.a, :_parent, WeakRef(ids))
        setfield!(ids.label, :_parent, WeakRef(ids))
        setfield!(ids.z_n, :_parent, WeakRef(ids))
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

mutable struct summary__heating_current_drive__nbi <: IDSvectorStaticElement
    var"angle" :: summary__heating_current_drive__nbi___angle
    var"beam_current_fraction" :: summary__heating_current_drive__nbi___beam_current_fraction
    var"beam_power_fraction" :: summary__heating_current_drive__nbi___beam_power_fraction
    var"current" :: summary__heating_current_drive__nbi___current
    var"direction" :: summary__heating_current_drive__nbi___direction
    var"energy" :: summary__heating_current_drive__nbi___energy
    var"position" :: summary__heating_current_drive__nbi___position
    var"power" :: summary__heating_current_drive__nbi___power
    var"power_launched" :: summary__heating_current_drive__nbi___power_launched
    var"species" :: summary__heating_current_drive__nbi___species
    var"tangency_radius" :: summary__heating_current_drive__nbi___tangency_radius
    _parent :: WeakRef
    function summary__heating_current_drive__nbi(var"angle"=summary__heating_current_drive__nbi___angle(), var"beam_current_fraction"=summary__heating_current_drive__nbi___beam_current_fraction(), var"beam_power_fraction"=summary__heating_current_drive__nbi___beam_power_fraction(), var"current"=summary__heating_current_drive__nbi___current(), var"direction"=summary__heating_current_drive__nbi___direction(), var"energy"=summary__heating_current_drive__nbi___energy(), var"position"=summary__heating_current_drive__nbi___position(), var"power"=summary__heating_current_drive__nbi___power(), var"power_launched"=summary__heating_current_drive__nbi___power_launched(), var"species"=summary__heating_current_drive__nbi___species(), var"tangency_radius"=summary__heating_current_drive__nbi___tangency_radius(), _parent=WeakRef(missing))
        ids = new(var"angle", var"beam_current_fraction", var"beam_power_fraction", var"current", var"direction", var"energy", var"position", var"power", var"power_launched", var"species", var"tangency_radius", _parent)
        assign_expressions(ids)
        setfield!(ids.angle, :_parent, WeakRef(ids))
        setfield!(ids.beam_current_fraction, :_parent, WeakRef(ids))
        setfield!(ids.beam_power_fraction, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.direction, :_parent, WeakRef(ids))
        setfield!(ids.energy, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
        setfield!(ids.species, :_parent, WeakRef(ids))
        setfield!(ids.tangency_radius, :_parent, WeakRef(ids))
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

mutable struct summary__heating_current_drive__lh <: IDSvectorStaticElement
    var"current" :: summary__heating_current_drive__lh___current
    var"energy_fast" :: summary__heating_current_drive__lh___energy_fast
    var"frequency" :: summary__heating_current_drive__lh___frequency
    var"n_parallel" :: summary__heating_current_drive__lh___n_parallel
    var"position" :: summary__heating_current_drive__lh___position
    var"power" :: summary__heating_current_drive__lh___power
    var"power_launched" :: summary__heating_current_drive__lh___power_launched
    _parent :: WeakRef
    function summary__heating_current_drive__lh(var"current"=summary__heating_current_drive__lh___current(), var"energy_fast"=summary__heating_current_drive__lh___energy_fast(), var"frequency"=summary__heating_current_drive__lh___frequency(), var"n_parallel"=summary__heating_current_drive__lh___n_parallel(), var"position"=summary__heating_current_drive__lh___position(), var"power"=summary__heating_current_drive__lh___power(), var"power_launched"=summary__heating_current_drive__lh___power_launched(), _parent=WeakRef(missing))
        ids = new(var"current", var"energy_fast", var"frequency", var"n_parallel", var"position", var"power", var"power_launched", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast, :_parent, WeakRef(ids))
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.n_parallel, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
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

mutable struct summary__heating_current_drive__ic <: IDSvectorStaticElement
    var"current" :: summary__heating_current_drive__ic___current
    var"e_field_plus_minus_ratio" :: summary__heating_current_drive__ic___e_field_plus_minus_ratio
    var"energy_fast" :: summary__heating_current_drive__ic___energy_fast
    var"frequency" :: summary__heating_current_drive__ic___frequency
    var"harmonic" :: summary__heating_current_drive__ic___harmonic
    var"k_perpendicular" :: summary__heating_current_drive__ic___k_perpendicular
    var"n_tor" :: summary__heating_current_drive__ic___n_tor
    var"phase" :: summary__heating_current_drive__ic___phase
    var"position" :: summary__heating_current_drive__ic___position
    var"power" :: summary__heating_current_drive__ic___power
    var"power_launched" :: summary__heating_current_drive__ic___power_launched
    _parent :: WeakRef
    function summary__heating_current_drive__ic(var"current"=summary__heating_current_drive__ic___current(), var"e_field_plus_minus_ratio"=summary__heating_current_drive__ic___e_field_plus_minus_ratio(), var"energy_fast"=summary__heating_current_drive__ic___energy_fast(), var"frequency"=summary__heating_current_drive__ic___frequency(), var"harmonic"=summary__heating_current_drive__ic___harmonic(), var"k_perpendicular"=summary__heating_current_drive__ic___k_perpendicular(), var"n_tor"=summary__heating_current_drive__ic___n_tor(), var"phase"=summary__heating_current_drive__ic___phase(), var"position"=summary__heating_current_drive__ic___position(), var"power"=summary__heating_current_drive__ic___power(), var"power_launched"=summary__heating_current_drive__ic___power_launched(), _parent=WeakRef(missing))
        ids = new(var"current", var"e_field_plus_minus_ratio", var"energy_fast", var"frequency", var"harmonic", var"k_perpendicular", var"n_tor", var"phase", var"position", var"power", var"power_launched", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.e_field_plus_minus_ratio, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast, :_parent, WeakRef(ids))
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.harmonic, :_parent, WeakRef(ids))
        setfield!(ids.k_perpendicular, :_parent, WeakRef(ids))
        setfield!(ids.n_tor, :_parent, WeakRef(ids))
        setfield!(ids.phase, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
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

mutable struct summary__heating_current_drive__ec <: IDSvectorStaticElement
    var"angle_pol" :: summary__heating_current_drive__ec___angle_pol
    var"angle_tor" :: summary__heating_current_drive__ec___angle_tor
    var"current" :: summary__heating_current_drive__ec___current
    var"energy_fast" :: summary__heating_current_drive__ec___energy_fast
    var"frequency" :: summary__heating_current_drive__ec___frequency
    var"harmonic" :: summary__heating_current_drive__ec___harmonic
    var"polarisation" :: summary__heating_current_drive__ec___polarisation
    var"position" :: summary__heating_current_drive__ec___position
    var"power" :: summary__heating_current_drive__ec___power
    var"power_launched" :: summary__heating_current_drive__ec___power_launched
    _parent :: WeakRef
    function summary__heating_current_drive__ec(var"angle_pol"=summary__heating_current_drive__ec___angle_pol(), var"angle_tor"=summary__heating_current_drive__ec___angle_tor(), var"current"=summary__heating_current_drive__ec___current(), var"energy_fast"=summary__heating_current_drive__ec___energy_fast(), var"frequency"=summary__heating_current_drive__ec___frequency(), var"harmonic"=summary__heating_current_drive__ec___harmonic(), var"polarisation"=summary__heating_current_drive__ec___polarisation(), var"position"=summary__heating_current_drive__ec___position(), var"power"=summary__heating_current_drive__ec___power(), var"power_launched"=summary__heating_current_drive__ec___power_launched(), _parent=WeakRef(missing))
        ids = new(var"angle_pol", var"angle_tor", var"current", var"energy_fast", var"frequency", var"harmonic", var"polarisation", var"position", var"power", var"power_launched", _parent)
        assign_expressions(ids)
        setfield!(ids.angle_pol, :_parent, WeakRef(ids))
        setfield!(ids.angle_tor, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast, :_parent, WeakRef(ids))
        setfield!(ids.frequency, :_parent, WeakRef(ids))
        setfield!(ids.harmonic, :_parent, WeakRef(ids))
        setfield!(ids.polarisation, :_parent, WeakRef(ids))
        setfield!(ids.position, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
        setfield!(ids.power_launched, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__heating_current_drive <: IDS
    var"ec" :: IDSvector{T} where {T<:summary__heating_current_drive__ec}
    var"ic" :: IDSvector{T} where {T<:summary__heating_current_drive__ic}
    var"lh" :: IDSvector{T} where {T<:summary__heating_current_drive__lh}
    var"nbi" :: IDSvector{T} where {T<:summary__heating_current_drive__nbi}
    var"power_additional" :: summary__heating_current_drive__power_additional
    var"power_ec" :: summary__heating_current_drive__power_ec
    var"power_ic" :: summary__heating_current_drive__power_ic
    var"power_launched_ec" :: summary__heating_current_drive__power_launched_ec
    var"power_launched_ic" :: summary__heating_current_drive__power_launched_ic
    var"power_launched_lh" :: summary__heating_current_drive__power_launched_lh
    var"power_launched_nbi" :: summary__heating_current_drive__power_launched_nbi
    var"power_launched_nbi_co_injected_ratio" :: summary__heating_current_drive__power_launched_nbi_co_injected_ratio
    var"power_lh" :: summary__heating_current_drive__power_lh
    var"power_nbi" :: summary__heating_current_drive__power_nbi
    _parent :: WeakRef
    function summary__heating_current_drive(var"ec"=IDSvector(summary__heating_current_drive__ec[]), var"ic"=IDSvector(summary__heating_current_drive__ic[]), var"lh"=IDSvector(summary__heating_current_drive__lh[]), var"nbi"=IDSvector(summary__heating_current_drive__nbi[]), var"power_additional"=summary__heating_current_drive__power_additional(), var"power_ec"=summary__heating_current_drive__power_ec(), var"power_ic"=summary__heating_current_drive__power_ic(), var"power_launched_ec"=summary__heating_current_drive__power_launched_ec(), var"power_launched_ic"=summary__heating_current_drive__power_launched_ic(), var"power_launched_lh"=summary__heating_current_drive__power_launched_lh(), var"power_launched_nbi"=summary__heating_current_drive__power_launched_nbi(), var"power_launched_nbi_co_injected_ratio"=summary__heating_current_drive__power_launched_nbi_co_injected_ratio(), var"power_lh"=summary__heating_current_drive__power_lh(), var"power_nbi"=summary__heating_current_drive__power_nbi(), _parent=WeakRef(missing))
        ids = new(var"ec", var"ic", var"lh", var"nbi", var"power_additional", var"power_ec", var"power_ic", var"power_launched_ec", var"power_launched_ic", var"power_launched_lh", var"power_launched_nbi", var"power_launched_nbi_co_injected_ratio", var"power_lh", var"power_nbi", _parent)
        assign_expressions(ids)
        setfield!(ids.ec, :_parent, WeakRef(ids))
        setfield!(ids.ic, :_parent, WeakRef(ids))
        setfield!(ids.lh, :_parent, WeakRef(ids))
        setfield!(ids.nbi, :_parent, WeakRef(ids))
        setfield!(ids.power_additional, :_parent, WeakRef(ids))
        setfield!(ids.power_ec, :_parent, WeakRef(ids))
        setfield!(ids.power_ic, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_ec, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_ic, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_lh, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_nbi, :_parent, WeakRef(ids))
        setfield!(ids.power_launched_nbi_co_injected_ratio, :_parent, WeakRef(ids))
        setfield!(ids.power_lh, :_parent, WeakRef(ids))
        setfield!(ids.power_nbi, :_parent, WeakRef(ids))
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
    var"b0" :: summary__global_quantities__b0
    var"beta_pol" :: summary__global_quantities__beta_pol
    var"beta_pol_mhd" :: summary__global_quantities__beta_pol_mhd
    var"beta_tor" :: summary__global_quantities__beta_tor
    var"beta_tor_mhd" :: summary__global_quantities__beta_tor_mhd
    var"beta_tor_norm" :: summary__global_quantities__beta_tor_norm
    var"beta_tor_norm_mhd" :: summary__global_quantities__beta_tor_norm_mhd
    var"beta_tor_thermal_norm" :: summary__global_quantities__beta_tor_thermal_norm
    var"current_alignment" :: summary__global_quantities__current_alignment
    var"current_bootstrap" :: summary__global_quantities__current_bootstrap
    var"current_non_inductive" :: summary__global_quantities__current_non_inductive
    var"current_ohm" :: summary__global_quantities__current_ohm
    var"denergy_diamagnetic_dt" :: summary__global_quantities__denergy_diamagnetic_dt
    var"denergy_thermal_dt" :: summary__global_quantities__denergy_thermal_dt
    var"energy_b_field_pol" :: summary__global_quantities__energy_b_field_pol
    var"energy_diamagnetic" :: summary__global_quantities__energy_diamagnetic
    var"energy_electrons_thermal" :: summary__global_quantities__energy_electrons_thermal
    var"energy_fast_parallel" :: summary__global_quantities__energy_fast_parallel
    var"energy_fast_perpendicular" :: summary__global_quantities__energy_fast_perpendicular
    var"energy_ion_total_thermal" :: summary__global_quantities__energy_ion_total_thermal
    var"energy_mhd" :: summary__global_quantities__energy_mhd
    var"energy_thermal" :: summary__global_quantities__energy_thermal
    var"energy_total" :: summary__global_quantities__energy_total
    var"fusion_fluence" :: summary__global_quantities__fusion_fluence
    var"fusion_gain" :: summary__global_quantities__fusion_gain
    var"greenwald_fraction" :: summary__global_quantities__greenwald_fraction
    var"h_98" :: summary__global_quantities__h_98
    var"h_mode" :: summary__global_quantities__h_mode
    var"ip" :: summary__global_quantities__ip
    var"li" :: summary__global_quantities__li
    var"li_mhd" :: summary__global_quantities__li_mhd
    var"power_bremsstrahlung" :: summary__global_quantities__power_bremsstrahlung
    var"power_line" :: summary__global_quantities__power_line
    var"power_loss" :: summary__global_quantities__power_loss
    var"power_ohm" :: summary__global_quantities__power_ohm
    var"power_radiated" :: summary__global_quantities__power_radiated
    var"power_radiated_inside_lcfs" :: summary__global_quantities__power_radiated_inside_lcfs
    var"power_radiated_outside_lcfs" :: summary__global_quantities__power_radiated_outside_lcfs
    var"power_steady" :: summary__global_quantities__power_steady
    var"power_synchrotron" :: summary__global_quantities__power_synchrotron
    var"q_95" :: summary__global_quantities__q_95
    var"r0" :: summary__global_quantities__r0
    var"ratio_tau_helium_fuel" :: summary__global_quantities__ratio_tau_helium_fuel
    var"resistance" :: summary__global_quantities__resistance
    var"tau_energy" :: summary__global_quantities__tau_energy
    var"tau_energy_98" :: summary__global_quantities__tau_energy_98
    var"tau_helium" :: summary__global_quantities__tau_helium
    var"tau_resistive" :: summary__global_quantities__tau_resistive
    var"v_loop" :: summary__global_quantities__v_loop
    var"volume" :: summary__global_quantities__volume
    _parent :: WeakRef
    function summary__global_quantities(var"b0"=summary__global_quantities__b0(), var"beta_pol"=summary__global_quantities__beta_pol(), var"beta_pol_mhd"=summary__global_quantities__beta_pol_mhd(), var"beta_tor"=summary__global_quantities__beta_tor(), var"beta_tor_mhd"=summary__global_quantities__beta_tor_mhd(), var"beta_tor_norm"=summary__global_quantities__beta_tor_norm(), var"beta_tor_norm_mhd"=summary__global_quantities__beta_tor_norm_mhd(), var"beta_tor_thermal_norm"=summary__global_quantities__beta_tor_thermal_norm(), var"current_alignment"=summary__global_quantities__current_alignment(), var"current_bootstrap"=summary__global_quantities__current_bootstrap(), var"current_non_inductive"=summary__global_quantities__current_non_inductive(), var"current_ohm"=summary__global_quantities__current_ohm(), var"denergy_diamagnetic_dt"=summary__global_quantities__denergy_diamagnetic_dt(), var"denergy_thermal_dt"=summary__global_quantities__denergy_thermal_dt(), var"energy_b_field_pol"=summary__global_quantities__energy_b_field_pol(), var"energy_diamagnetic"=summary__global_quantities__energy_diamagnetic(), var"energy_electrons_thermal"=summary__global_quantities__energy_electrons_thermal(), var"energy_fast_parallel"=summary__global_quantities__energy_fast_parallel(), var"energy_fast_perpendicular"=summary__global_quantities__energy_fast_perpendicular(), var"energy_ion_total_thermal"=summary__global_quantities__energy_ion_total_thermal(), var"energy_mhd"=summary__global_quantities__energy_mhd(), var"energy_thermal"=summary__global_quantities__energy_thermal(), var"energy_total"=summary__global_quantities__energy_total(), var"fusion_fluence"=summary__global_quantities__fusion_fluence(), var"fusion_gain"=summary__global_quantities__fusion_gain(), var"greenwald_fraction"=summary__global_quantities__greenwald_fraction(), var"h_98"=summary__global_quantities__h_98(), var"h_mode"=summary__global_quantities__h_mode(), var"ip"=summary__global_quantities__ip(), var"li"=summary__global_quantities__li(), var"li_mhd"=summary__global_quantities__li_mhd(), var"power_bremsstrahlung"=summary__global_quantities__power_bremsstrahlung(), var"power_line"=summary__global_quantities__power_line(), var"power_loss"=summary__global_quantities__power_loss(), var"power_ohm"=summary__global_quantities__power_ohm(), var"power_radiated"=summary__global_quantities__power_radiated(), var"power_radiated_inside_lcfs"=summary__global_quantities__power_radiated_inside_lcfs(), var"power_radiated_outside_lcfs"=summary__global_quantities__power_radiated_outside_lcfs(), var"power_steady"=summary__global_quantities__power_steady(), var"power_synchrotron"=summary__global_quantities__power_synchrotron(), var"q_95"=summary__global_quantities__q_95(), var"r0"=summary__global_quantities__r0(), var"ratio_tau_helium_fuel"=summary__global_quantities__ratio_tau_helium_fuel(), var"resistance"=summary__global_quantities__resistance(), var"tau_energy"=summary__global_quantities__tau_energy(), var"tau_energy_98"=summary__global_quantities__tau_energy_98(), var"tau_helium"=summary__global_quantities__tau_helium(), var"tau_resistive"=summary__global_quantities__tau_resistive(), var"v_loop"=summary__global_quantities__v_loop(), var"volume"=summary__global_quantities__volume(), _parent=WeakRef(missing))
        ids = new(var"b0", var"beta_pol", var"beta_pol_mhd", var"beta_tor", var"beta_tor_mhd", var"beta_tor_norm", var"beta_tor_norm_mhd", var"beta_tor_thermal_norm", var"current_alignment", var"current_bootstrap", var"current_non_inductive", var"current_ohm", var"denergy_diamagnetic_dt", var"denergy_thermal_dt", var"energy_b_field_pol", var"energy_diamagnetic", var"energy_electrons_thermal", var"energy_fast_parallel", var"energy_fast_perpendicular", var"energy_ion_total_thermal", var"energy_mhd", var"energy_thermal", var"energy_total", var"fusion_fluence", var"fusion_gain", var"greenwald_fraction", var"h_98", var"h_mode", var"ip", var"li", var"li_mhd", var"power_bremsstrahlung", var"power_line", var"power_loss", var"power_ohm", var"power_radiated", var"power_radiated_inside_lcfs", var"power_radiated_outside_lcfs", var"power_steady", var"power_synchrotron", var"q_95", var"r0", var"ratio_tau_helium_fuel", var"resistance", var"tau_energy", var"tau_energy_98", var"tau_helium", var"tau_resistive", var"v_loop", var"volume", _parent)
        assign_expressions(ids)
        setfield!(ids.b0, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol, :_parent, WeakRef(ids))
        setfield!(ids.beta_pol_mhd, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_mhd, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_norm, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_norm_mhd, :_parent, WeakRef(ids))
        setfield!(ids.beta_tor_thermal_norm, :_parent, WeakRef(ids))
        setfield!(ids.current_alignment, :_parent, WeakRef(ids))
        setfield!(ids.current_bootstrap, :_parent, WeakRef(ids))
        setfield!(ids.current_non_inductive, :_parent, WeakRef(ids))
        setfield!(ids.current_ohm, :_parent, WeakRef(ids))
        setfield!(ids.denergy_diamagnetic_dt, :_parent, WeakRef(ids))
        setfield!(ids.denergy_thermal_dt, :_parent, WeakRef(ids))
        setfield!(ids.energy_b_field_pol, :_parent, WeakRef(ids))
        setfield!(ids.energy_diamagnetic, :_parent, WeakRef(ids))
        setfield!(ids.energy_electrons_thermal, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast_parallel, :_parent, WeakRef(ids))
        setfield!(ids.energy_fast_perpendicular, :_parent, WeakRef(ids))
        setfield!(ids.energy_ion_total_thermal, :_parent, WeakRef(ids))
        setfield!(ids.energy_mhd, :_parent, WeakRef(ids))
        setfield!(ids.energy_thermal, :_parent, WeakRef(ids))
        setfield!(ids.energy_total, :_parent, WeakRef(ids))
        setfield!(ids.fusion_fluence, :_parent, WeakRef(ids))
        setfield!(ids.fusion_gain, :_parent, WeakRef(ids))
        setfield!(ids.greenwald_fraction, :_parent, WeakRef(ids))
        setfield!(ids.h_98, :_parent, WeakRef(ids))
        setfield!(ids.h_mode, :_parent, WeakRef(ids))
        setfield!(ids.ip, :_parent, WeakRef(ids))
        setfield!(ids.li, :_parent, WeakRef(ids))
        setfield!(ids.li_mhd, :_parent, WeakRef(ids))
        setfield!(ids.power_bremsstrahlung, :_parent, WeakRef(ids))
        setfield!(ids.power_line, :_parent, WeakRef(ids))
        setfield!(ids.power_loss, :_parent, WeakRef(ids))
        setfield!(ids.power_ohm, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated_inside_lcfs, :_parent, WeakRef(ids))
        setfield!(ids.power_radiated_outside_lcfs, :_parent, WeakRef(ids))
        setfield!(ids.power_steady, :_parent, WeakRef(ids))
        setfield!(ids.power_synchrotron, :_parent, WeakRef(ids))
        setfield!(ids.q_95, :_parent, WeakRef(ids))
        setfield!(ids.r0, :_parent, WeakRef(ids))
        setfield!(ids.ratio_tau_helium_fuel, :_parent, WeakRef(ids))
        setfield!(ids.resistance, :_parent, WeakRef(ids))
        setfield!(ids.tau_energy, :_parent, WeakRef(ids))
        setfield!(ids.tau_energy_98, :_parent, WeakRef(ids))
        setfield!(ids.tau_helium, :_parent, WeakRef(ids))
        setfield!(ids.tau_resistive, :_parent, WeakRef(ids))
        setfield!(ids.v_loop, :_parent, WeakRef(ids))
        setfield!(ids.volume, :_parent, WeakRef(ids))
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
    var"ammonia" :: summary__gas_injection_rates__ammonia
    var"ammonia_deuterated" :: summary__gas_injection_rates__ammonia_deuterated
    var"argon" :: summary__gas_injection_rates__argon
    var"beryllium" :: summary__gas_injection_rates__beryllium
    var"bottom" :: summary__gas_injection_rates__bottom
    var"carbon" :: summary__gas_injection_rates__carbon
    var"deuterium" :: summary__gas_injection_rates__deuterium
    var"ethane" :: summary__gas_injection_rates__ethane
    var"ethylene" :: summary__gas_injection_rates__ethylene
    var"helium_3" :: summary__gas_injection_rates__helium_3
    var"helium_4" :: summary__gas_injection_rates__helium_4
    var"hydrogen" :: summary__gas_injection_rates__hydrogen
    var"impurity_seeding" :: summary__gas_injection_rates__impurity_seeding
    var"krypton" :: summary__gas_injection_rates__krypton
    var"lithium" :: summary__gas_injection_rates__lithium
    var"methane" :: summary__gas_injection_rates__methane
    var"methane_carbon_13" :: summary__gas_injection_rates__methane_carbon_13
    var"methane_deuterated" :: summary__gas_injection_rates__methane_deuterated
    var"midplane" :: summary__gas_injection_rates__midplane
    var"neon" :: summary__gas_injection_rates__neon
    var"nitrogen" :: summary__gas_injection_rates__nitrogen
    var"oxygen" :: summary__gas_injection_rates__oxygen
    var"propane" :: summary__gas_injection_rates__propane
    var"silane" :: summary__gas_injection_rates__silane
    var"top" :: summary__gas_injection_rates__top
    var"total" :: summary__gas_injection_rates__total
    var"tritium" :: summary__gas_injection_rates__tritium
    var"xenon" :: summary__gas_injection_rates__xenon
    _parent :: WeakRef
    function summary__gas_injection_rates(var"ammonia"=summary__gas_injection_rates__ammonia(), var"ammonia_deuterated"=summary__gas_injection_rates__ammonia_deuterated(), var"argon"=summary__gas_injection_rates__argon(), var"beryllium"=summary__gas_injection_rates__beryllium(), var"bottom"=summary__gas_injection_rates__bottom(), var"carbon"=summary__gas_injection_rates__carbon(), var"deuterium"=summary__gas_injection_rates__deuterium(), var"ethane"=summary__gas_injection_rates__ethane(), var"ethylene"=summary__gas_injection_rates__ethylene(), var"helium_3"=summary__gas_injection_rates__helium_3(), var"helium_4"=summary__gas_injection_rates__helium_4(), var"hydrogen"=summary__gas_injection_rates__hydrogen(), var"impurity_seeding"=summary__gas_injection_rates__impurity_seeding(), var"krypton"=summary__gas_injection_rates__krypton(), var"lithium"=summary__gas_injection_rates__lithium(), var"methane"=summary__gas_injection_rates__methane(), var"methane_carbon_13"=summary__gas_injection_rates__methane_carbon_13(), var"methane_deuterated"=summary__gas_injection_rates__methane_deuterated(), var"midplane"=summary__gas_injection_rates__midplane(), var"neon"=summary__gas_injection_rates__neon(), var"nitrogen"=summary__gas_injection_rates__nitrogen(), var"oxygen"=summary__gas_injection_rates__oxygen(), var"propane"=summary__gas_injection_rates__propane(), var"silane"=summary__gas_injection_rates__silane(), var"top"=summary__gas_injection_rates__top(), var"total"=summary__gas_injection_rates__total(), var"tritium"=summary__gas_injection_rates__tritium(), var"xenon"=summary__gas_injection_rates__xenon(), _parent=WeakRef(missing))
        ids = new(var"ammonia", var"ammonia_deuterated", var"argon", var"beryllium", var"bottom", var"carbon", var"deuterium", var"ethane", var"ethylene", var"helium_3", var"helium_4", var"hydrogen", var"impurity_seeding", var"krypton", var"lithium", var"methane", var"methane_carbon_13", var"methane_deuterated", var"midplane", var"neon", var"nitrogen", var"oxygen", var"propane", var"silane", var"top", var"total", var"tritium", var"xenon", _parent)
        assign_expressions(ids)
        setfield!(ids.ammonia, :_parent, WeakRef(ids))
        setfield!(ids.ammonia_deuterated, :_parent, WeakRef(ids))
        setfield!(ids.argon, :_parent, WeakRef(ids))
        setfield!(ids.beryllium, :_parent, WeakRef(ids))
        setfield!(ids.bottom, :_parent, WeakRef(ids))
        setfield!(ids.carbon, :_parent, WeakRef(ids))
        setfield!(ids.deuterium, :_parent, WeakRef(ids))
        setfield!(ids.ethane, :_parent, WeakRef(ids))
        setfield!(ids.ethylene, :_parent, WeakRef(ids))
        setfield!(ids.helium_3, :_parent, WeakRef(ids))
        setfield!(ids.helium_4, :_parent, WeakRef(ids))
        setfield!(ids.hydrogen, :_parent, WeakRef(ids))
        setfield!(ids.impurity_seeding, :_parent, WeakRef(ids))
        setfield!(ids.krypton, :_parent, WeakRef(ids))
        setfield!(ids.lithium, :_parent, WeakRef(ids))
        setfield!(ids.methane, :_parent, WeakRef(ids))
        setfield!(ids.methane_carbon_13, :_parent, WeakRef(ids))
        setfield!(ids.methane_deuterated, :_parent, WeakRef(ids))
        setfield!(ids.midplane, :_parent, WeakRef(ids))
        setfield!(ids.neon, :_parent, WeakRef(ids))
        setfield!(ids.nitrogen, :_parent, WeakRef(ids))
        setfield!(ids.oxygen, :_parent, WeakRef(ids))
        setfield!(ids.propane, :_parent, WeakRef(ids))
        setfield!(ids.silane, :_parent, WeakRef(ids))
        setfield!(ids.top, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.tritium, :_parent, WeakRef(ids))
        setfield!(ids.xenon, :_parent, WeakRef(ids))
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
    var"beam_beam" :: summary__fusion__neutron_fluxes__tt__beam_beam
    var"beam_thermal" :: summary__fusion__neutron_fluxes__tt__beam_thermal
    var"thermal" :: summary__fusion__neutron_fluxes__tt__thermal
    var"total" :: summary__fusion__neutron_fluxes__tt__total
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__tt(var"beam_beam"=summary__fusion__neutron_fluxes__tt__beam_beam(), var"beam_thermal"=summary__fusion__neutron_fluxes__tt__beam_thermal(), var"thermal"=summary__fusion__neutron_fluxes__tt__thermal(), var"total"=summary__fusion__neutron_fluxes__tt__total(), _parent=WeakRef(missing))
        ids = new(var"beam_beam", var"beam_thermal", var"thermal", var"total", _parent)
        assign_expressions(ids)
        setfield!(ids.beam_beam, :_parent, WeakRef(ids))
        setfield!(ids.beam_thermal, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
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
    var"beam_beam" :: summary__fusion__neutron_fluxes__dt__beam_beam
    var"beam_thermal" :: summary__fusion__neutron_fluxes__dt__beam_thermal
    var"thermal" :: summary__fusion__neutron_fluxes__dt__thermal
    var"total" :: summary__fusion__neutron_fluxes__dt__total
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dt(var"beam_beam"=summary__fusion__neutron_fluxes__dt__beam_beam(), var"beam_thermal"=summary__fusion__neutron_fluxes__dt__beam_thermal(), var"thermal"=summary__fusion__neutron_fluxes__dt__thermal(), var"total"=summary__fusion__neutron_fluxes__dt__total(), _parent=WeakRef(missing))
        ids = new(var"beam_beam", var"beam_thermal", var"thermal", var"total", _parent)
        assign_expressions(ids)
        setfield!(ids.beam_beam, :_parent, WeakRef(ids))
        setfield!(ids.beam_thermal, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
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
    var"beam_beam" :: summary__fusion__neutron_fluxes__dd__beam_beam
    var"beam_thermal" :: summary__fusion__neutron_fluxes__dd__beam_thermal
    var"thermal" :: summary__fusion__neutron_fluxes__dd__thermal
    var"total" :: summary__fusion__neutron_fluxes__dd__total
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes__dd(var"beam_beam"=summary__fusion__neutron_fluxes__dd__beam_beam(), var"beam_thermal"=summary__fusion__neutron_fluxes__dd__beam_thermal(), var"thermal"=summary__fusion__neutron_fluxes__dd__thermal(), var"total"=summary__fusion__neutron_fluxes__dd__total(), _parent=WeakRef(missing))
        ids = new(var"beam_beam", var"beam_thermal", var"thermal", var"total", _parent)
        assign_expressions(ids)
        setfield!(ids.beam_beam, :_parent, WeakRef(ids))
        setfield!(ids.beam_thermal, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary__fusion__neutron_fluxes <: IDS
    var"dd" :: summary__fusion__neutron_fluxes__dd
    var"dt" :: summary__fusion__neutron_fluxes__dt
    var"thermal" :: summary__fusion__neutron_fluxes__thermal
    var"total" :: summary__fusion__neutron_fluxes__total
    var"tt" :: summary__fusion__neutron_fluxes__tt
    _parent :: WeakRef
    function summary__fusion__neutron_fluxes(var"dd"=summary__fusion__neutron_fluxes__dd(), var"dt"=summary__fusion__neutron_fluxes__dt(), var"thermal"=summary__fusion__neutron_fluxes__thermal(), var"total"=summary__fusion__neutron_fluxes__total(), var"tt"=summary__fusion__neutron_fluxes__tt(), _parent=WeakRef(missing))
        ids = new(var"dd", var"dt", var"thermal", var"total", var"tt", _parent)
        assign_expressions(ids)
        setfield!(ids.dd, :_parent, WeakRef(ids))
        setfield!(ids.dt, :_parent, WeakRef(ids))
        setfield!(ids.thermal, :_parent, WeakRef(ids))
        setfield!(ids.total, :_parent, WeakRef(ids))
        setfield!(ids.tt, :_parent, WeakRef(ids))
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
    var"current" :: summary__fusion__current
    var"neutron_fluxes" :: summary__fusion__neutron_fluxes
    var"neutron_power_total" :: summary__fusion__neutron_power_total
    var"power" :: summary__fusion__power
    _parent :: WeakRef
    function summary__fusion(var"current"=summary__fusion__current(), var"neutron_fluxes"=summary__fusion__neutron_fluxes(), var"neutron_power_total"=summary__fusion__neutron_power_total(), var"power"=summary__fusion__power(), _parent=WeakRef(missing))
        ids = new(var"current", var"neutron_fluxes", var"neutron_power_total", var"power", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.neutron_fluxes, :_parent, WeakRef(ids))
        setfield!(ids.neutron_power_total, :_parent, WeakRef(ids))
        setfield!(ids.power, :_parent, WeakRef(ids))
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
    var"mitigation_valve" :: summary__disruption__mitigation_valve
    var"time" :: summary__disruption__time
    var"time_half_ip" :: summary__disruption__time_half_ip
    var"time_radiated_power_max" :: summary__disruption__time_radiated_power_max
    var"vertical_displacement" :: summary__disruption__vertical_displacement
    _parent :: WeakRef
    function summary__disruption(var"mitigation_valve"=summary__disruption__mitigation_valve(), var"time"=summary__disruption__time(), var"time_half_ip"=summary__disruption__time_half_ip(), var"time_radiated_power_max"=summary__disruption__time_radiated_power_max(), var"vertical_displacement"=summary__disruption__vertical_displacement(), _parent=WeakRef(missing))
        ids = new(var"mitigation_valve", var"time", var"time_half_ip", var"time_radiated_power_max", var"vertical_displacement", _parent)
        assign_expressions(ids)
        setfield!(ids.mitigation_valve, :_parent, WeakRef(ids))
        setfield!(ids.time, :_parent, WeakRef(ids))
        setfield!(ids.time_half_ip, :_parent, WeakRef(ids))
        setfield!(ids.time_radiated_power_max, :_parent, WeakRef(ids))
        setfield!(ids.vertical_displacement, :_parent, WeakRef(ids))
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
    var"elongation" :: summary__boundary__elongation
    var"gap_limiter_wall" :: summary__boundary__gap_limiter_wall
    var"geometric_axis_r" :: summary__boundary__geometric_axis_r
    var"geometric_axis_z" :: summary__boundary__geometric_axis_z
    var"magnetic_axis_r" :: summary__boundary__magnetic_axis_r
    var"magnetic_axis_z" :: summary__boundary__magnetic_axis_z
    var"minor_radius" :: summary__boundary__minor_radius
    var"strike_point_configuration" :: summary__boundary__strike_point_configuration
    var"strike_point_inner_r" :: summary__boundary__strike_point_inner_r
    var"strike_point_inner_z" :: summary__boundary__strike_point_inner_z
    var"strike_point_outer_r" :: summary__boundary__strike_point_outer_r
    var"strike_point_outer_z" :: summary__boundary__strike_point_outer_z
    var"triangularity_lower" :: summary__boundary__triangularity_lower
    var"triangularity_upper" :: summary__boundary__triangularity_upper
    var"type" :: summary__boundary__type
    _parent :: WeakRef
    function summary__boundary(var"elongation"=summary__boundary__elongation(), var"gap_limiter_wall"=summary__boundary__gap_limiter_wall(), var"geometric_axis_r"=summary__boundary__geometric_axis_r(), var"geometric_axis_z"=summary__boundary__geometric_axis_z(), var"magnetic_axis_r"=summary__boundary__magnetic_axis_r(), var"magnetic_axis_z"=summary__boundary__magnetic_axis_z(), var"minor_radius"=summary__boundary__minor_radius(), var"strike_point_configuration"=summary__boundary__strike_point_configuration(), var"strike_point_inner_r"=summary__boundary__strike_point_inner_r(), var"strike_point_inner_z"=summary__boundary__strike_point_inner_z(), var"strike_point_outer_r"=summary__boundary__strike_point_outer_r(), var"strike_point_outer_z"=summary__boundary__strike_point_outer_z(), var"triangularity_lower"=summary__boundary__triangularity_lower(), var"triangularity_upper"=summary__boundary__triangularity_upper(), var"type"=summary__boundary__type(), _parent=WeakRef(missing))
        ids = new(var"elongation", var"gap_limiter_wall", var"geometric_axis_r", var"geometric_axis_z", var"magnetic_axis_r", var"magnetic_axis_z", var"minor_radius", var"strike_point_configuration", var"strike_point_inner_r", var"strike_point_inner_z", var"strike_point_outer_r", var"strike_point_outer_z", var"triangularity_lower", var"triangularity_upper", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.elongation, :_parent, WeakRef(ids))
        setfield!(ids.gap_limiter_wall, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis_r, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis_z, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis_r, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis_z, :_parent, WeakRef(ids))
        setfield!(ids.minor_radius, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_configuration, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_inner_r, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_inner_z, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_outer_r, :_parent, WeakRef(ids))
        setfield!(ids.strike_point_outer_z, :_parent, WeakRef(ids))
        setfield!(ids.triangularity_lower, :_parent, WeakRef(ids))
        setfield!(ids.triangularity_upper, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct summary <: IDS
    var"boundary" :: summary__boundary
    var"configuration" :: summary__configuration
    var"disruption" :: summary__disruption
    var"elms" :: summary__elms
    var"fusion" :: summary__fusion
    var"gas_injection_rates" :: summary__gas_injection_rates
    var"global_quantities" :: summary__global_quantities
    var"heating_current_drive" :: summary__heating_current_drive
    var"kicks" :: summary__kicks
    var"limiter" :: summary__limiter
    var"line_average" :: summary__line_average
    var"local" :: summary__local
    var"magnetic_shear_flag" :: summary__magnetic_shear_flag
    var"midplane" :: summary__midplane
    var"pedestal_fits" :: summary__pedestal_fits
    var"pellets" :: summary__pellets
    var"rmps" :: summary__rmps
    var"runaways" :: summary__runaways
    var"scrape_off_layer" :: summary__scrape_off_layer
    var"stationary_phase_flag" :: summary__stationary_phase_flag
    var"tag" :: summary__tag
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume_average" :: summary__volume_average
    var"wall" :: summary__wall
    _parent :: WeakRef
    function summary(var"boundary"=summary__boundary(), var"configuration"=summary__configuration(), var"disruption"=summary__disruption(), var"elms"=summary__elms(), var"fusion"=summary__fusion(), var"gas_injection_rates"=summary__gas_injection_rates(), var"global_quantities"=summary__global_quantities(), var"heating_current_drive"=summary__heating_current_drive(), var"kicks"=summary__kicks(), var"limiter"=summary__limiter(), var"line_average"=summary__line_average(), var"local"=summary__local(), var"magnetic_shear_flag"=summary__magnetic_shear_flag(), var"midplane"=summary__midplane(), var"pedestal_fits"=summary__pedestal_fits(), var"pellets"=summary__pellets(), var"rmps"=summary__rmps(), var"runaways"=summary__runaways(), var"scrape_off_layer"=summary__scrape_off_layer(), var"stationary_phase_flag"=summary__stationary_phase_flag(), var"tag"=summary__tag(), var"time"=missing, var"time_width"=missing, var"volume_average"=summary__volume_average(), var"wall"=summary__wall(), _parent=WeakRef(missing))
        ids = new(var"boundary", var"configuration", var"disruption", var"elms", var"fusion", var"gas_injection_rates", var"global_quantities", var"heating_current_drive", var"kicks", var"limiter", var"line_average", var"local", var"magnetic_shear_flag", var"midplane", var"pedestal_fits", var"pellets", var"rmps", var"runaways", var"scrape_off_layer", var"stationary_phase_flag", var"tag", var"time", var"time_width", var"volume_average", var"wall", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        setfield!(ids.configuration, :_parent, WeakRef(ids))
        setfield!(ids.disruption, :_parent, WeakRef(ids))
        setfield!(ids.elms, :_parent, WeakRef(ids))
        setfield!(ids.fusion, :_parent, WeakRef(ids))
        setfield!(ids.gas_injection_rates, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.heating_current_drive, :_parent, WeakRef(ids))
        setfield!(ids.kicks, :_parent, WeakRef(ids))
        setfield!(ids.limiter, :_parent, WeakRef(ids))
        setfield!(ids.line_average, :_parent, WeakRef(ids))
        setfield!(ids.local, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_shear_flag, :_parent, WeakRef(ids))
        setfield!(ids.midplane, :_parent, WeakRef(ids))
        setfield!(ids.pedestal_fits, :_parent, WeakRef(ids))
        setfield!(ids.pellets, :_parent, WeakRef(ids))
        setfield!(ids.rmps, :_parent, WeakRef(ids))
        setfield!(ids.runaways, :_parent, WeakRef(ids))
        setfield!(ids.scrape_off_layer, :_parent, WeakRef(ids))
        setfield!(ids.stationary_phase_flag, :_parent, WeakRef(ids))
        setfield!(ids.tag, :_parent, WeakRef(ids))
        setfield!(ids.volume_average, :_parent, WeakRef(ids))
        setfield!(ids.wall, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__vertical_force___force <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__vertical_force___force(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__vertical_force <: IDSvectorStaticElement
    var"combination" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"force" :: pf_active__vertical_force___force
    var"limit_max" :: Union{Missing, Real, Function}
    var"limit_min" :: Union{Missing, Real, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function pf_active__vertical_force(var"combination"=missing, var"force"=pf_active__vertical_force___force(), var"limit_max"=missing, var"limit_min"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"combination", var"force", var"limit_max", var"limit_min", var"name", _parent)
        assign_expressions(ids)
        setfield!(ids.force, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__supply___voltage <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__supply___voltage(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__supply___current <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__supply___current(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__supply <: IDSvectorStaticElement
    var"current" :: pf_active__supply___current
    var"current_limit_max" :: Union{Missing, Real, Function}
    var"current_limit_min" :: Union{Missing, Real, Function}
    var"current_limiter_gain" :: Union{Missing, Real, Function}
    var"delay" :: Union{Missing, Real, Function}
    var"energy_limit_max" :: Union{Missing, Real, Function}
    var"filter_denominator" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"filter_numerator" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"nonlinear_model" :: Union{Missing, String, Function}
    var"resistance" :: Union{Missing, Real, Function}
    var"type" :: Union{Missing, Integer, Function}
    var"voltage" :: pf_active__supply___voltage
    var"voltage_limit_max" :: Union{Missing, Real, Function}
    var"voltage_limit_min" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__supply(var"current"=pf_active__supply___current(), var"current_limit_max"=missing, var"current_limit_min"=missing, var"current_limiter_gain"=missing, var"delay"=missing, var"energy_limit_max"=missing, var"filter_denominator"=missing, var"filter_numerator"=missing, var"identifier"=missing, var"name"=missing, var"nonlinear_model"=missing, var"resistance"=missing, var"type"=missing, var"voltage"=pf_active__supply___voltage(), var"voltage_limit_max"=missing, var"voltage_limit_min"=missing, _parent=WeakRef(missing))
        ids = new(var"current", var"current_limit_max", var"current_limit_min", var"current_limiter_gain", var"delay", var"energy_limit_max", var"filter_denominator", var"filter_numerator", var"identifier", var"name", var"nonlinear_model", var"resistance", var"type", var"voltage", var"voltage_limit_max", var"voltage_limit_min", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__radial_force___force <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__radial_force___force(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__radial_force <: IDSvectorStaticElement
    var"combination" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"force" :: pf_active__radial_force___force
    var"limit_max" :: Union{Missing, Real, Function}
    var"limit_min" :: Union{Missing, Real, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function pf_active__radial_force(var"combination"=missing, var"force"=pf_active__radial_force___force(), var"limit_max"=missing, var"limit_min"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"combination", var"force", var"limit_max", var"limit_min", var"name", _parent)
        assign_expressions(ids)
        setfield!(ids.force, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__global_quantities <: IDS
    var"psi_coils_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_coils_list" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__global_quantities(var"psi_coils_average"=missing, var"psi_coils_list"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"psi_coils_average", var"psi_coils_list", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___voltage <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___voltage(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
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
    var"beta" :: Union{Missing, Real, Function}
    var"length_alpha" :: Union{Missing, Real, Function}
    var"length_beta" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__coil___element___geometry__oblique(var"alpha"=missing, var"beta"=missing, var"length_alpha"=missing, var"length_beta"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"alpha", var"beta", var"length_alpha", var"length_beta", var"r", var"z", _parent)
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
    var"geometry_type" :: Union{Missing, Integer, Function}
    var"oblique" :: pf_active__coil___element___geometry__oblique
    var"outline" :: pf_active__coil___element___geometry__outline
    var"rectangle" :: pf_active__coil___element___geometry__rectangle
    _parent :: WeakRef
    function pf_active__coil___element___geometry(var"arcs_of_circle"=pf_active__coil___element___geometry__arcs_of_circle(), var"geometry_type"=missing, var"oblique"=pf_active__coil___element___geometry__oblique(), var"outline"=pf_active__coil___element___geometry__outline(), var"rectangle"=pf_active__coil___element___geometry__rectangle(), _parent=WeakRef(missing))
        ids = new(var"arcs_of_circle", var"geometry_type", var"oblique", var"outline", var"rectangle", _parent)
        assign_expressions(ids)
        setfield!(ids.arcs_of_circle, :_parent, WeakRef(ids))
        setfield!(ids.oblique, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.rectangle, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__coil___element <: IDSvectorStaticElement
    var"area" :: Union{Missing, Real, Function}
    var"geometry" :: pf_active__coil___element___geometry
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"turns_with_sign" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function pf_active__coil___element(var"area"=missing, var"geometry"=pf_active__coil___element___geometry(), var"identifier"=missing, var"name"=missing, var"turns_with_sign"=missing, _parent=WeakRef(missing))
        ids = new(var"area", var"geometry", var"identifier", var"name", var"turns_with_sign", _parent)
        assign_expressions(ids)
        setfield!(ids.geometry, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__coil___current <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___current(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil___b_field_max_timed <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__coil___b_field_max_timed(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__coil <: IDSvectorStaticElement
    var"b_field_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_field_max_timed" :: pf_active__coil___b_field_max_timed
    var"current" :: pf_active__coil___current
    var"current_limit_max" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"element" :: IDSvector{T} where {T<:pf_active__coil___element}
    var"energy_limit_max" :: Union{Missing, Real, Function}
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"resistance" :: Union{Missing, Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"voltage" :: pf_active__coil___voltage
    _parent :: WeakRef
    function pf_active__coil(var"b_field_max"=missing, var"b_field_max_timed"=pf_active__coil___b_field_max_timed(), var"current"=pf_active__coil___current(), var"current_limit_max"=missing, var"element"=IDSvector(pf_active__coil___element[]), var"energy_limit_max"=missing, var"identifier"=missing, var"name"=missing, var"resistance"=missing, var"temperature"=missing, var"voltage"=pf_active__coil___voltage(), _parent=WeakRef(missing))
        ids = new(var"b_field_max", var"b_field_max_timed", var"current", var"current_limit_max", var"element", var"energy_limit_max", var"identifier", var"name", var"resistance", var"temperature", var"voltage", _parent)
        assign_expressions(ids)
        setfield!(ids.b_field_max_timed, :_parent, WeakRef(ids))
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active__circuit___voltage <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__circuit___voltage(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__circuit___current <: IDS
    var"data" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function pf_active__circuit___current(var"data"=missing, var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data", var"time", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct pf_active__circuit <: IDSvectorStaticElement
    var"connections" :: Union{Missing, AbstractArray{T, 2} where T<:Integer, Function}
    var"current" :: pf_active__circuit___current
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"type" :: Union{Missing, String, Function}
    var"voltage" :: pf_active__circuit___voltage
    _parent :: WeakRef
    function pf_active__circuit(var"connections"=missing, var"current"=pf_active__circuit___current(), var"identifier"=missing, var"name"=missing, var"type"=missing, var"voltage"=pf_active__circuit___voltage(), _parent=WeakRef(missing))
        ids = new(var"connections", var"current", var"identifier", var"name", var"type", var"voltage", _parent)
        assign_expressions(ids)
        setfield!(ids.current, :_parent, WeakRef(ids))
        setfield!(ids.voltage, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct pf_active <: IDS
    var"circuit" :: IDSvector{T} where {T<:pf_active__circuit}
    var"coil" :: IDSvector{T} where {T<:pf_active__coil}
    var"global_quantities" :: pf_active__global_quantities
    var"latency" :: Union{Missing, Real, Function}
    var"radial_force" :: IDSvector{T} where {T<:pf_active__radial_force}
    var"supply" :: IDSvector{T} where {T<:pf_active__supply}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vertical_force" :: IDSvector{T} where {T<:pf_active__vertical_force}
    _parent :: WeakRef
    function pf_active(var"circuit"=IDSvector(pf_active__circuit[]), var"coil"=IDSvector(pf_active__coil[]), var"global_quantities"=pf_active__global_quantities(), var"latency"=missing, var"radial_force"=IDSvector(pf_active__radial_force[]), var"supply"=IDSvector(pf_active__supply[]), var"time"=missing, var"vertical_force"=IDSvector(pf_active__vertical_force[]), _parent=WeakRef(missing))
        ids = new(var"circuit", var"coil", var"global_quantities", var"latency", var"radial_force", var"supply", var"time", var"vertical_force", _parent)
        assign_expressions(ids)
        setfield!(ids.circuit, :_parent, WeakRef(ids))
        setfield!(ids.coil, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.radial_force, :_parent, WeakRef(ids))
        setfield!(ids.supply, :_parent, WeakRef(ids))
        setfield!(ids.vertical_force, :_parent, WeakRef(ids))
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_2d___grid_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_2d___grid <: IDS
    var"dim1" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dim2" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume_element" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_2d___grid(var"dim1"=missing, var"dim2"=missing, var"volume_element"=missing, _parent=WeakRef(missing))
        ids = new(var"dim1", var"dim2", var"volume_element", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___profiles_2d <: IDSvectorStaticElement
    var"b_field_r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_field_tor" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"b_field_z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"grid" :: equilibrium__time_slice___profiles_2d___grid
    var"grid_type" :: equilibrium__time_slice___profiles_2d___grid_type
    var"j_parallel" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"j_tor" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"phi" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"psi" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"theta" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_2d(var"b_field_r"=missing, var"b_field_tor"=missing, var"b_field_z"=missing, var"grid"=equilibrium__time_slice___profiles_2d___grid(), var"grid_type"=equilibrium__time_slice___profiles_2d___grid_type(), var"j_parallel"=missing, var"j_tor"=missing, var"phi"=missing, var"psi"=missing, var"r"=missing, var"theta"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"b_field_r", var"b_field_tor", var"b_field_z", var"grid", var"grid_type", var"j_parallel", var"j_tor", var"phi", var"psi", var"r", var"theta", var"z", _parent)
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
    var"area" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_field_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_field_max" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"b_field_min" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"beta_pol" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"darea_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"darea_drho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dpressure_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dpsi_drho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dvolume_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dvolume_drho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"elongation" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"f" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"f_df_dpsi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"geometric_axis" :: equilibrium__time_slice___profiles_1d__geometric_axis
    var"gm1" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm2" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm3" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm4" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm5" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm6" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm7" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm8" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"gm9" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"magnetic_shear" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"mass_density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"phi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"q" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r_inboard" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r_outboard" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_volume_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_lower_inner" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_lower_outer" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_upper_inner" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"squareness_upper_outer" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"surface" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"trapped_fraction" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"triangularity_lower" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"triangularity_upper" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___profiles_1d(var"area"=missing, var"b_field_average"=missing, var"b_field_max"=missing, var"b_field_min"=missing, var"beta_pol"=missing, var"darea_dpsi"=missing, var"darea_drho_tor"=missing, var"dpressure_dpsi"=missing, var"dpsi_drho_tor"=missing, var"dvolume_dpsi"=missing, var"dvolume_drho_tor"=missing, var"elongation"=missing, var"f"=missing, var"f_df_dpsi"=missing, var"geometric_axis"=equilibrium__time_slice___profiles_1d__geometric_axis(), var"gm1"=missing, var"gm2"=missing, var"gm3"=missing, var"gm4"=missing, var"gm5"=missing, var"gm6"=missing, var"gm7"=missing, var"gm8"=missing, var"gm9"=missing, var"j_parallel"=missing, var"j_tor"=missing, var"magnetic_shear"=missing, var"mass_density"=missing, var"phi"=missing, var"pressure"=missing, var"psi"=missing, var"q"=missing, var"r_inboard"=missing, var"r_outboard"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, var"rho_volume_norm"=missing, var"squareness_lower_inner"=missing, var"squareness_lower_outer"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"surface"=missing, var"trapped_fraction"=missing, var"triangularity_lower"=missing, var"triangularity_upper"=missing, var"volume"=missing, _parent=WeakRef(missing))
        ids = new(var"area", var"b_field_average", var"b_field_max", var"b_field_min", var"beta_pol", var"darea_dpsi", var"darea_drho_tor", var"dpressure_dpsi", var"dpsi_drho_tor", var"dvolume_dpsi", var"dvolume_drho_tor", var"elongation", var"f", var"f_df_dpsi", var"geometric_axis", var"gm1", var"gm2", var"gm3", var"gm4", var"gm5", var"gm6", var"gm7", var"gm8", var"gm9", var"j_parallel", var"j_tor", var"magnetic_shear", var"mass_density", var"phi", var"pressure", var"psi", var"q", var"r_inboard", var"r_outboard", var"rho_tor", var"rho_tor_norm", var"rho_volume_norm", var"squareness_lower_inner", var"squareness_lower_outer", var"squareness_upper_inner", var"squareness_upper_outer", var"surface", var"trapped_fraction", var"triangularity_lower", var"triangularity_upper", var"volume", _parent)
        assign_expressions(ids)
        setfield!(ids.geometric_axis, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities__q_min <: IDS
    var"rho_tor_norm" :: Union{Missing, Real, Function}
    var"value" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities__q_min(var"rho_tor_norm"=missing, var"value"=missing, _parent=WeakRef(missing))
        ids = new(var"rho_tor_norm", var"value", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities__magnetic_axis <: IDS
    var"b_field_tor" :: Union{Missing, Real, Function}
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities__magnetic_axis(var"b_field_tor"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"b_field_tor", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities__current_centre <: IDS
    var"r" :: Union{Missing, Real, Function}
    var"velocity_z" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities__current_centre(var"r"=missing, var"velocity_z"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"velocity_z", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___global_quantities <: IDS
    var"area" :: Union{Missing, Real, Function}
    var"beta_normal" :: Union{Missing, Real, Function}
    var"beta_pol" :: Union{Missing, Real, Function}
    var"beta_tor" :: Union{Missing, Real, Function}
    var"current_centre" :: equilibrium__time_slice___global_quantities__current_centre
    var"energy_mhd" :: Union{Missing, Real, Function}
    var"ip" :: Union{Missing, Real, Function}
    var"length_pol" :: Union{Missing, Real, Function}
    var"li_3" :: Union{Missing, Real, Function}
    var"magnetic_axis" :: equilibrium__time_slice___global_quantities__magnetic_axis
    var"plasma_inductance" :: Union{Missing, Real, Function}
    var"psi_axis" :: Union{Missing, Real, Function}
    var"psi_boundary" :: Union{Missing, Real, Function}
    var"psi_external_average" :: Union{Missing, Real, Function}
    var"q_95" :: Union{Missing, Real, Function}
    var"q_axis" :: Union{Missing, Real, Function}
    var"q_min" :: equilibrium__time_slice___global_quantities__q_min
    var"surface" :: Union{Missing, Real, Function}
    var"volume" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___global_quantities(var"area"=missing, var"beta_normal"=missing, var"beta_pol"=missing, var"beta_tor"=missing, var"current_centre"=equilibrium__time_slice___global_quantities__current_centre(), var"energy_mhd"=missing, var"ip"=missing, var"length_pol"=missing, var"li_3"=missing, var"magnetic_axis"=equilibrium__time_slice___global_quantities__magnetic_axis(), var"plasma_inductance"=missing, var"psi_axis"=missing, var"psi_boundary"=missing, var"psi_external_average"=missing, var"q_95"=missing, var"q_axis"=missing, var"q_min"=equilibrium__time_slice___global_quantities__q_min(), var"surface"=missing, var"volume"=missing, _parent=WeakRef(missing))
        ids = new(var"area", var"beta_normal", var"beta_pol", var"beta_tor", var"current_centre", var"energy_mhd", var"ip", var"length_pol", var"li_3", var"magnetic_axis", var"plasma_inductance", var"psi_axis", var"psi_boundary", var"psi_external_average", var"q_95", var"q_axis", var"q_min", var"surface", var"volume", _parent)
        assign_expressions(ids)
        setfield!(ids.current_centre, :_parent, WeakRef(ids))
        setfield!(ids.magnetic_axis, :_parent, WeakRef(ids))
        setfield!(ids.q_min, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___coordinate_system__grid_type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___coordinate_system__grid_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___coordinate_system__grid <: IDS
    var"dim1" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"dim2" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume_element" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___coordinate_system__grid(var"dim1"=missing, var"dim2"=missing, var"volume_element"=missing, _parent=WeakRef(missing))
        ids = new(var"dim1", var"dim2", var"volume_element", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___coordinate_system <: IDS
    var"grid" :: equilibrium__time_slice___coordinate_system__grid
    var"grid_type" :: equilibrium__time_slice___coordinate_system__grid_type
    var"jacobian" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    var"tensor_contravariant" :: Union{Missing, AbstractArray{T, 4} where T<:Real, Function}
    var"tensor_covariant" :: Union{Missing, AbstractArray{T, 4} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 2} where T<:Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___coordinate_system(var"grid"=equilibrium__time_slice___coordinate_system__grid(), var"grid_type"=equilibrium__time_slice___coordinate_system__grid_type(), var"jacobian"=missing, var"r"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"grid", var"grid_type", var"jacobian", var"r", var"tensor_contravariant", var"tensor_covariant", var"z", _parent)
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

mutable struct equilibrium__time_slice___constraints__x_point <: IDSvectorStaticElement
    var"chi_squared_r" :: Union{Missing, Real, Function}
    var"chi_squared_z" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"position_measured" :: equilibrium__time_slice___constraints__x_point___position_measured
    var"position_reconstructed" :: equilibrium__time_slice___constraints__x_point___position_reconstructed
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__x_point(var"chi_squared_r"=missing, var"chi_squared_z"=missing, var"exact"=missing, var"position_measured"=equilibrium__time_slice___constraints__x_point___position_measured(), var"position_reconstructed"=equilibrium__time_slice___constraints__x_point___position_reconstructed(), var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared_r", var"chi_squared_z", var"exact", var"position_measured", var"position_reconstructed", var"source", var"time_measurement", var"weight", _parent)
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

mutable struct equilibrium__time_slice___constraints__strike_point <: IDSvectorStaticElement
    var"chi_squared_r" :: Union{Missing, Real, Function}
    var"chi_squared_z" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"position_measured" :: equilibrium__time_slice___constraints__strike_point___position_measured
    var"position_reconstructed" :: equilibrium__time_slice___constraints__strike_point___position_reconstructed
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__strike_point(var"chi_squared_r"=missing, var"chi_squared_z"=missing, var"exact"=missing, var"position_measured"=equilibrium__time_slice___constraints__strike_point___position_measured(), var"position_reconstructed"=equilibrium__time_slice___constraints__strike_point___position_reconstructed(), var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared_r", var"chi_squared_z", var"exact", var"position_measured", var"position_reconstructed", var"source", var"time_measurement", var"weight", _parent)
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

mutable struct equilibrium__time_slice___constraints__q <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"position" :: equilibrium__time_slice___constraints__q___position
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__q(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"position"=equilibrium__time_slice___constraints__q___position(), var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"position", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.position, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__pressure <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__pressure(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__pf_passive_current <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__pf_passive_current(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__pf_current <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__pf_current(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__n_e_line <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__n_e_line(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__n_e <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__n_e(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__mse_polarisation_angle <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__mse_polarisation_angle(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__iron_core_segment <: IDSvectorStaticElement
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
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__ip(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__flux_loop <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__flux_loop(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__faraday_angle <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__faraday_angle(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__diamagnetic_flux <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__diamagnetic_flux(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__bpol_probe <: IDSvectorStaticElement
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__bpol_probe(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints__b_field_tor_vacuum_r <: IDS
    var"chi_squared" :: Union{Missing, Real, Function}
    var"exact" :: Union{Missing, Integer, Function}
    var"measured" :: Union{Missing, Real, Function}
    var"reconstructed" :: Union{Missing, Real, Function}
    var"source" :: Union{Missing, String, Function}
    var"time_measurement" :: Union{Missing, Real, Function}
    var"weight" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints__b_field_tor_vacuum_r(var"chi_squared"=missing, var"exact"=missing, var"measured"=missing, var"reconstructed"=missing, var"source"=missing, var"time_measurement"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"exact", var"measured", var"reconstructed", var"source", var"time_measurement", var"weight", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___constraints <: IDS
    var"b_field_tor_vacuum_r" :: equilibrium__time_slice___constraints__b_field_tor_vacuum_r
    var"bpol_probe" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__bpol_probe}
    var"diamagnetic_flux" :: equilibrium__time_slice___constraints__diamagnetic_flux
    var"faraday_angle" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__faraday_angle}
    var"flux_loop" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__flux_loop}
    var"ip" :: equilibrium__time_slice___constraints__ip
    var"iron_core_segment" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__iron_core_segment}
    var"mse_polarisation_angle" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__mse_polarisation_angle}
    var"n_e" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__n_e}
    var"n_e_line" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__n_e_line}
    var"pf_current" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__pf_current}
    var"pf_passive_current" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__pf_passive_current}
    var"pressure" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__pressure}
    var"q" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__q}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__strike_point}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___constraints__x_point}
    _parent :: WeakRef
    function equilibrium__time_slice___constraints(var"b_field_tor_vacuum_r"=equilibrium__time_slice___constraints__b_field_tor_vacuum_r(), var"bpol_probe"=IDSvector(equilibrium__time_slice___constraints__bpol_probe[]), var"diamagnetic_flux"=equilibrium__time_slice___constraints__diamagnetic_flux(), var"faraday_angle"=IDSvector(equilibrium__time_slice___constraints__faraday_angle[]), var"flux_loop"=IDSvector(equilibrium__time_slice___constraints__flux_loop[]), var"ip"=equilibrium__time_slice___constraints__ip(), var"iron_core_segment"=IDSvector(equilibrium__time_slice___constraints__iron_core_segment[]), var"mse_polarisation_angle"=IDSvector(equilibrium__time_slice___constraints__mse_polarisation_angle[]), var"n_e"=IDSvector(equilibrium__time_slice___constraints__n_e[]), var"n_e_line"=IDSvector(equilibrium__time_slice___constraints__n_e_line[]), var"pf_current"=IDSvector(equilibrium__time_slice___constraints__pf_current[]), var"pf_passive_current"=IDSvector(equilibrium__time_slice___constraints__pf_passive_current[]), var"pressure"=IDSvector(equilibrium__time_slice___constraints__pressure[]), var"q"=IDSvector(equilibrium__time_slice___constraints__q[]), var"strike_point"=IDSvector(equilibrium__time_slice___constraints__strike_point[]), var"x_point"=IDSvector(equilibrium__time_slice___constraints__x_point[]), _parent=WeakRef(missing))
        ids = new(var"b_field_tor_vacuum_r", var"bpol_probe", var"diamagnetic_flux", var"faraday_angle", var"flux_loop", var"ip", var"iron_core_segment", var"mse_polarisation_angle", var"n_e", var"n_e_line", var"pf_current", var"pf_passive_current", var"pressure", var"q", var"strike_point", var"x_point", _parent)
        assign_expressions(ids)
        setfield!(ids.b_field_tor_vacuum_r, :_parent, WeakRef(ids))
        setfield!(ids.bpol_probe, :_parent, WeakRef(ids))
        setfield!(ids.diamagnetic_flux, :_parent, WeakRef(ids))
        setfield!(ids.faraday_angle, :_parent, WeakRef(ids))
        setfield!(ids.flux_loop, :_parent, WeakRef(ids))
        setfield!(ids.ip, :_parent, WeakRef(ids))
        setfield!(ids.iron_core_segment, :_parent, WeakRef(ids))
        setfield!(ids.mse_polarisation_angle, :_parent, WeakRef(ids))
        setfield!(ids.n_e, :_parent, WeakRef(ids))
        setfield!(ids.n_e_line, :_parent, WeakRef(ids))
        setfield!(ids.pf_current, :_parent, WeakRef(ids))
        setfield!(ids.pf_passive_current, :_parent, WeakRef(ids))
        setfield!(ids.pressure, :_parent, WeakRef(ids))
        setfield!(ids.q, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__x_point <: IDSvectorStaticElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_separatrix__strike_point <: IDSvectorStaticElement
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

mutable struct equilibrium__time_slice___boundary_separatrix__gap <: IDSvectorStaticElement
    var"angle" :: Union{Missing, Real, Function}
    var"identifier" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"r" :: Union{Missing, Real, Function}
    var"value" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix__gap(var"angle"=missing, var"identifier"=missing, var"name"=missing, var"r"=missing, var"value"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"angle", var"identifier", var"name", var"r", var"value", var"z", _parent)
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
    var"active_limiter_point" :: equilibrium__time_slice___boundary_separatrix__active_limiter_point
    var"closest_wall_point" :: equilibrium__time_slice___boundary_separatrix__closest_wall_point
    var"dr_dz_zero_point" :: equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point
    var"elongation" :: Union{Missing, Real, Function}
    var"elongation_lower" :: Union{Missing, Real, Function}
    var"elongation_upper" :: Union{Missing, Real, Function}
    var"gap" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__gap}
    var"geometric_axis" :: equilibrium__time_slice___boundary_separatrix__geometric_axis
    var"minor_radius" :: Union{Missing, Real, Function}
    var"outline" :: equilibrium__time_slice___boundary_separatrix__outline
    var"psi" :: Union{Missing, Real, Function}
    var"squareness_lower_inner" :: Union{Missing, Real, Function}
    var"squareness_lower_outer" :: Union{Missing, Real, Function}
    var"squareness_upper_inner" :: Union{Missing, Real, Function}
    var"squareness_upper_outer" :: Union{Missing, Real, Function}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__strike_point}
    var"triangularity" :: Union{Missing, Real, Function}
    var"triangularity_lower" :: Union{Missing, Real, Function}
    var"triangularity_upper" :: Union{Missing, Real, Function}
    var"type" :: Union{Missing, Integer, Function}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__x_point}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_separatrix(var"active_limiter_point"=equilibrium__time_slice___boundary_separatrix__active_limiter_point(), var"closest_wall_point"=equilibrium__time_slice___boundary_separatrix__closest_wall_point(), var"dr_dz_zero_point"=equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(), var"elongation"=missing, var"elongation_lower"=missing, var"elongation_upper"=missing, var"gap"=IDSvector(equilibrium__time_slice___boundary_separatrix__gap[]), var"geometric_axis"=equilibrium__time_slice___boundary_separatrix__geometric_axis(), var"minor_radius"=missing, var"outline"=equilibrium__time_slice___boundary_separatrix__outline(), var"psi"=missing, var"squareness_lower_inner"=missing, var"squareness_lower_outer"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"strike_point"=IDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[]), var"triangularity"=missing, var"triangularity_lower"=missing, var"triangularity_upper"=missing, var"type"=missing, var"x_point"=IDSvector(equilibrium__time_slice___boundary_separatrix__x_point[]), _parent=WeakRef(missing))
        ids = new(var"active_limiter_point", var"closest_wall_point", var"dr_dz_zero_point", var"elongation", var"elongation_lower", var"elongation_upper", var"gap", var"geometric_axis", var"minor_radius", var"outline", var"psi", var"squareness_lower_inner", var"squareness_lower_outer", var"squareness_upper_inner", var"squareness_upper_outer", var"strike_point", var"triangularity", var"triangularity_lower", var"triangularity_upper", var"type", var"x_point", _parent)
        assign_expressions(ids)
        setfield!(ids.active_limiter_point, :_parent, WeakRef(ids))
        setfield!(ids.closest_wall_point, :_parent, WeakRef(ids))
        setfield!(ids.dr_dz_zero_point, :_parent, WeakRef(ids))
        setfield!(ids.gap, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_secondary_separatrix__x_point <: IDSvectorStaticElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_secondary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary_secondary_separatrix__strike_point <: IDSvectorStaticElement
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
    var"distance_inner_outer" :: Union{Missing, Real, Function}
    var"outline" :: equilibrium__time_slice___boundary_secondary_separatrix__outline
    var"psi" :: Union{Missing, Real, Function}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__strike_point}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__x_point}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary_secondary_separatrix(var"distance_inner_outer"=missing, var"outline"=equilibrium__time_slice___boundary_secondary_separatrix__outline(), var"psi"=missing, var"strike_point"=IDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[]), var"x_point"=IDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[]), _parent=WeakRef(missing))
        ids = new(var"distance_inner_outer", var"outline", var"psi", var"strike_point", var"x_point", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__x_point <: IDSvectorStaticElement
    var"r" :: Union{Missing, Real, Function}
    var"z" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct equilibrium__time_slice___boundary__strike_point <: IDSvectorStaticElement
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
    var"active_limiter_point" :: equilibrium__time_slice___boundary__active_limiter_point
    var"elongation" :: Union{Missing, Real, Function}
    var"elongation_lower" :: Union{Missing, Real, Function}
    var"elongation_upper" :: Union{Missing, Real, Function}
    var"geometric_axis" :: equilibrium__time_slice___boundary__geometric_axis
    var"minor_radius" :: Union{Missing, Real, Function}
    var"outline" :: equilibrium__time_slice___boundary__outline
    var"psi" :: Union{Missing, Real, Function}
    var"psi_norm" :: Union{Missing, Real, Function}
    var"squareness_lower_inner" :: Union{Missing, Real, Function}
    var"squareness_lower_outer" :: Union{Missing, Real, Function}
    var"squareness_upper_inner" :: Union{Missing, Real, Function}
    var"squareness_upper_outer" :: Union{Missing, Real, Function}
    var"strike_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary__strike_point}
    var"triangularity" :: Union{Missing, Real, Function}
    var"triangularity_lower" :: Union{Missing, Real, Function}
    var"triangularity_upper" :: Union{Missing, Real, Function}
    var"type" :: Union{Missing, Integer, Function}
    var"x_point" :: IDSvector{T} where {T<:equilibrium__time_slice___boundary__x_point}
    _parent :: WeakRef
    function equilibrium__time_slice___boundary(var"active_limiter_point"=equilibrium__time_slice___boundary__active_limiter_point(), var"elongation"=missing, var"elongation_lower"=missing, var"elongation_upper"=missing, var"geometric_axis"=equilibrium__time_slice___boundary__geometric_axis(), var"minor_radius"=missing, var"outline"=equilibrium__time_slice___boundary__outline(), var"psi"=missing, var"psi_norm"=missing, var"squareness_lower_inner"=missing, var"squareness_lower_outer"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"strike_point"=IDSvector(equilibrium__time_slice___boundary__strike_point[]), var"triangularity"=missing, var"triangularity_lower"=missing, var"triangularity_upper"=missing, var"type"=missing, var"x_point"=IDSvector(equilibrium__time_slice___boundary__x_point[]), _parent=WeakRef(missing))
        ids = new(var"active_limiter_point", var"elongation", var"elongation_lower", var"elongation_upper", var"geometric_axis", var"minor_radius", var"outline", var"psi", var"psi_norm", var"squareness_lower_inner", var"squareness_lower_outer", var"squareness_upper_inner", var"squareness_upper_outer", var"strike_point", var"triangularity", var"triangularity_lower", var"triangularity_upper", var"type", var"x_point", _parent)
        assign_expressions(ids)
        setfield!(ids.active_limiter_point, :_parent, WeakRef(ids))
        setfield!(ids.geometric_axis, :_parent, WeakRef(ids))
        setfield!(ids.outline, :_parent, WeakRef(ids))
        setfield!(ids.strike_point, :_parent, WeakRef(ids))
        setfield!(ids.x_point, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium__time_slice <: IDSvectorTimeElement
    var"boundary" :: equilibrium__time_slice___boundary
    var"boundary_secondary_separatrix" :: equilibrium__time_slice___boundary_secondary_separatrix
    var"boundary_separatrix" :: equilibrium__time_slice___boundary_separatrix
    var"constraints" :: equilibrium__time_slice___constraints
    var"convergence" :: equilibrium__time_slice___convergence
    var"coordinate_system" :: equilibrium__time_slice___coordinate_system
    var"global_quantities" :: equilibrium__time_slice___global_quantities
    var"profiles_1d" :: equilibrium__time_slice___profiles_1d
    var"profiles_2d" :: IDSvector{T} where {T<:equilibrium__time_slice___profiles_2d}
    var"time" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function equilibrium__time_slice(var"boundary"=equilibrium__time_slice___boundary(), var"boundary_secondary_separatrix"=equilibrium__time_slice___boundary_secondary_separatrix(), var"boundary_separatrix"=equilibrium__time_slice___boundary_separatrix(), var"constraints"=equilibrium__time_slice___constraints(), var"convergence"=equilibrium__time_slice___convergence(), var"coordinate_system"=equilibrium__time_slice___coordinate_system(), var"global_quantities"=equilibrium__time_slice___global_quantities(), var"profiles_1d"=equilibrium__time_slice___profiles_1d(), var"profiles_2d"=IDSvector(equilibrium__time_slice___profiles_2d[]), var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"boundary", var"boundary_secondary_separatrix", var"boundary_separatrix", var"constraints", var"convergence", var"coordinate_system", var"global_quantities", var"profiles_1d", var"profiles_2d", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.boundary, :_parent, WeakRef(ids))
        setfield!(ids.boundary_secondary_separatrix, :_parent, WeakRef(ids))
        setfield!(ids.boundary_separatrix, :_parent, WeakRef(ids))
        setfield!(ids.constraints, :_parent, WeakRef(ids))
        setfield!(ids.convergence, :_parent, WeakRef(ids))
        setfield!(ids.coordinate_system, :_parent, WeakRef(ids))
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.profiles_1d, :_parent, WeakRef(ids))
        setfield!(ids.profiles_2d, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct equilibrium <: IDS
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_slice" :: IDSvector{T} where {T<:equilibrium__time_slice}
    var"vacuum_toroidal_field" :: equilibrium__vacuum_toroidal_field
    _parent :: WeakRef
    function equilibrium(var"time"=missing, var"time_slice"=IDSvector(equilibrium__time_slice[]), var"vacuum_toroidal_field"=equilibrium__vacuum_toroidal_field(), _parent=WeakRef(missing))
        ids = new(var"time", var"time_slice", var"vacuum_toroidal_field", _parent)
        assign_expressions(ids)
        setfield!(ids.time_slice, :_parent, WeakRef(ids))
        setfield!(ids.vacuum_toroidal_field, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct dataset_description__simulation <: IDS
    var"comment_after" :: Union{Missing, String, Function}
    var"comment_before" :: Union{Missing, String, Function}
    var"time_begin" :: Union{Missing, Real, Function}
    var"time_begun" :: Union{Missing, String, Function}
    var"time_current" :: Union{Missing, Real, Function}
    var"time_end" :: Union{Missing, Real, Function}
    var"time_ended" :: Union{Missing, String, Function}
    var"time_restart" :: Union{Missing, Real, Function}
    var"time_step" :: Union{Missing, Real, Function}
    var"workflow" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description__simulation(var"comment_after"=missing, var"comment_before"=missing, var"time_begin"=missing, var"time_begun"=missing, var"time_current"=missing, var"time_end"=missing, var"time_ended"=missing, var"time_restart"=missing, var"time_step"=missing, var"workflow"=missing, _parent=WeakRef(missing))
        ids = new(var"comment_after", var"comment_before", var"time_begin", var"time_begun", var"time_current", var"time_end", var"time_ended", var"time_restart", var"time_step", var"workflow", _parent)
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
    var"machine" :: Union{Missing, String, Function}
    var"pulse" :: Union{Missing, Integer, Function}
    var"pulse_type" :: Union{Missing, String, Function}
    var"run" :: Union{Missing, Integer, Function}
    var"user" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description__parent_entry(var"machine"=missing, var"pulse"=missing, var"pulse_type"=missing, var"run"=missing, var"user"=missing, _parent=WeakRef(missing))
        ids = new(var"machine", var"pulse", var"pulse_type", var"run", var"user", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description__data_entry <: IDS
    var"machine" :: Union{Missing, String, Function}
    var"pulse" :: Union{Missing, Integer, Function}
    var"pulse_type" :: Union{Missing, String, Function}
    var"run" :: Union{Missing, Integer, Function}
    var"user" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function dataset_description__data_entry(var"machine"=missing, var"pulse"=missing, var"pulse_type"=missing, var"run"=missing, var"user"=missing, _parent=WeakRef(missing))
        ids = new(var"machine", var"pulse", var"pulse_type", var"run", var"user", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct dataset_description <: IDS
    var"data_entry" :: dataset_description__data_entry
    var"dd_version" :: Union{Missing, String, Function}
    var"imas_version" :: Union{Missing, String, Function}
    var"parent_entry" :: dataset_description__parent_entry
    var"pulse_time_begin" :: Union{Missing, String, Function}
    var"pulse_time_begin_epoch" :: dataset_description__pulse_time_begin_epoch
    var"pulse_time_end_epoch" :: dataset_description__pulse_time_end_epoch
    var"simulation" :: dataset_description__simulation
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function dataset_description(var"data_entry"=dataset_description__data_entry(), var"dd_version"=missing, var"imas_version"=missing, var"parent_entry"=dataset_description__parent_entry(), var"pulse_time_begin"=missing, var"pulse_time_begin_epoch"=dataset_description__pulse_time_begin_epoch(), var"pulse_time_end_epoch"=dataset_description__pulse_time_end_epoch(), var"simulation"=dataset_description__simulation(), var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"data_entry", var"dd_version", var"imas_version", var"parent_entry", var"pulse_time_begin", var"pulse_time_begin_epoch", var"pulse_time_end_epoch", var"simulation", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.data_entry, :_parent, WeakRef(ids))
        setfield!(ids.parent_entry, :_parent, WeakRef(ids))
        setfield!(ids.pulse_time_begin_epoch, :_parent, WeakRef(ids))
        setfield!(ids.pulse_time_end_epoch, :_parent, WeakRef(ids))
        setfield!(ids.simulation, :_parent, WeakRef(ids))
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__source___species__type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__neutral__state__neutral_type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__source___species__neutral__state__neutral_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__neutral__state <: IDS
    var"electron_configuration" :: Union{Missing, String, Function}
    var"label" :: Union{Missing, String, Function}
    var"neutral_type" :: core_sources__source___species__neutral__state__neutral_type
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__source___species__neutral__state(var"electron_configuration"=missing, var"label"=missing, var"neutral_type"=core_sources__source___species__neutral__state__neutral_type(), var"vibrational_level"=missing, var"vibrational_mode"=missing, _parent=WeakRef(missing))
        ids = new(var"electron_configuration", var"label", var"neutral_type", var"vibrational_level", var"vibrational_mode", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___species__neutral__element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__neutral__element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__neutral <: IDS
    var"element" :: IDSvector{T} where {T<:core_sources__source___species__neutral__element}
    var"label" :: Union{Missing, String, Function}
    var"state" :: core_sources__source___species__neutral__state
    _parent :: WeakRef
    function core_sources__source___species__neutral(var"element"=IDSvector(core_sources__source___species__neutral__element[]), var"label"=missing, var"state"=core_sources__source___species__neutral__state(), _parent=WeakRef(missing))
        ids = new(var"element", var"label", var"state", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___species__ion__state <: IDS
    var"electron_configuration" :: Union{Missing, String, Function}
    var"label" :: Union{Missing, String, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"z_max" :: Union{Missing, Real, Function}
    var"z_min" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__ion__state(var"electron_configuration"=missing, var"label"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, var"z_max"=missing, var"z_min"=missing, _parent=WeakRef(missing))
        ids = new(var"electron_configuration", var"label", var"vibrational_level", var"vibrational_mode", var"z_max", var"z_min", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__ion__element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__ion__element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___species__ion <: IDS
    var"element" :: IDSvector{T} where {T<:core_sources__source___species__ion__element}
    var"label" :: Union{Missing, String, Function}
    var"state" :: core_sources__source___species__ion__state
    var"z_ion" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___species__ion(var"element"=IDSvector(core_sources__source___species__ion__element[]), var"label"=missing, var"state"=core_sources__source___species__ion__state(), var"z_ion"=missing, _parent=WeakRef(missing))
        ids = new(var"element", var"label", var"state", var"z_ion", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___species <: IDS
    var"ion" :: core_sources__source___species__ion
    var"neutral" :: core_sources__source___species__neutral
    var"type" :: core_sources__source___species__type
    _parent :: WeakRef
    function core_sources__source___species(var"ion"=core_sources__source___species__ion(), var"neutral"=core_sources__source___species__neutral(), var"type"=core_sources__source___species__type(), _parent=WeakRef(missing))
        ids = new(var"ion", var"neutral", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.ion, :_parent, WeakRef(ids))
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        setfield!(ids.type, :_parent, WeakRef(ids))
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral___state___neutral_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral___state <: IDSvectorStaticElement
    var"electron_configuration" :: Union{Missing, String, Function}
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"label" :: Union{Missing, String, Function}
    var"neutral_type" :: core_sources__source___profiles_1d___neutral___state___neutral_type
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral___state(var"electron_configuration"=missing, var"energy"=missing, var"label"=missing, var"neutral_type"=core_sources__source___profiles_1d___neutral___state___neutral_type(), var"particles"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, _parent=WeakRef(missing))
        ids = new(var"electron_configuration", var"energy", var"label", var"neutral_type", var"particles", var"vibrational_level", var"vibrational_mode", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral___element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral___element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___neutral <: IDSvectorStaticElement
    var"element" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___neutral___element}
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ion_index" :: Union{Missing, Integer, Function}
    var"label" :: Union{Missing, String, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"state" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___neutral___state}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___neutral(var"element"=IDSvector(core_sources__source___profiles_1d___neutral___element[]), var"energy"=missing, var"ion_index"=missing, var"label"=missing, var"multiple_states_flag"=missing, var"particles"=missing, var"state"=IDSvector(core_sources__source___profiles_1d___neutral___state[]), _parent=WeakRef(missing))
        ids = new(var"element", var"energy", var"ion_index", var"label", var"multiple_states_flag", var"particles", var"state", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
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

mutable struct core_sources__source___profiles_1d___ion___state <: IDSvectorStaticElement
    var"electron_configuration" :: Union{Missing, String, Function}
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"energy_decomposed" :: core_sources__source___profiles_1d___ion___state___energy_decomposed
    var"label" :: Union{Missing, String, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particles_decomposed" :: core_sources__source___profiles_1d___ion___state___particles_decomposed
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"z_max" :: Union{Missing, Real, Function}
    var"z_min" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___state(var"electron_configuration"=missing, var"energy"=missing, var"energy_decomposed"=core_sources__source___profiles_1d___ion___state___energy_decomposed(), var"label"=missing, var"particles"=missing, var"particles_decomposed"=core_sources__source___profiles_1d___ion___state___particles_decomposed(), var"vibrational_level"=missing, var"vibrational_mode"=missing, var"z_max"=missing, var"z_min"=missing, _parent=WeakRef(missing))
        ids = new(var"electron_configuration", var"energy", var"energy_decomposed", var"label", var"particles", var"particles_decomposed", var"vibrational_level", var"vibrational_mode", var"z_max", var"z_min", _parent)
        assign_expressions(ids)
        setfield!(ids.energy_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.particles_decomposed, :_parent, WeakRef(ids))
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
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal_decomposed" :: core_sources__source___profiles_1d___ion___momentum__toroidal_decomposed
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___momentum(var"diamagnetic"=missing, var"parallel"=missing, var"poloidal"=missing, var"radial"=missing, var"toroidal"=missing, var"toroidal_decomposed"=core_sources__source___profiles_1d___ion___momentum__toroidal_decomposed(), _parent=WeakRef(missing))
        ids = new(var"diamagnetic", var"parallel", var"poloidal", var"radial", var"toroidal", var"toroidal_decomposed", _parent)
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

mutable struct core_sources__source___profiles_1d___ion___element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion___element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___ion <: IDSvectorStaticElement
    var"element" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___ion___element}
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"energy_decomposed" :: core_sources__source___profiles_1d___ion___energy_decomposed
    var"label" :: Union{Missing, String, Function}
    var"momentum" :: core_sources__source___profiles_1d___ion___momentum
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"neutral_index" :: Union{Missing, Integer, Function}
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particles_decomposed" :: core_sources__source___profiles_1d___ion___particles_decomposed
    var"state" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___ion___state}
    var"z_ion" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___ion(var"element"=IDSvector(core_sources__source___profiles_1d___ion___element[]), var"energy"=missing, var"energy_decomposed"=core_sources__source___profiles_1d___ion___energy_decomposed(), var"label"=missing, var"momentum"=core_sources__source___profiles_1d___ion___momentum(), var"multiple_states_flag"=missing, var"neutral_index"=missing, var"particles"=missing, var"particles_decomposed"=core_sources__source___profiles_1d___ion___particles_decomposed(), var"state"=IDSvector(core_sources__source___profiles_1d___ion___state[]), var"z_ion"=missing, _parent=WeakRef(missing))
        ids = new(var"element", var"energy", var"energy_decomposed", var"label", var"momentum", var"multiple_states_flag", var"neutral_index", var"particles", var"particles_decomposed", var"state", var"z_ion", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.energy_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.momentum, :_parent, WeakRef(ids))
        setfield!(ids.particles_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d___grid <: IDS
    var"area" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_boundary" :: Union{Missing, Real, Function}
    var"psi_magnetic_axis" :: Union{Missing, Real, Function}
    var"rho_pol_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"surface" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___grid(var"area"=missing, var"psi"=missing, var"psi_boundary"=missing, var"psi_magnetic_axis"=missing, var"rho_pol_norm"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, var"surface"=missing, var"volume"=missing, _parent=WeakRef(missing))
        ids = new(var"area", var"psi", var"psi_boundary", var"psi_magnetic_axis", var"rho_pol_norm", var"rho_tor", var"rho_tor_norm", var"surface", var"volume", _parent)
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
    var"energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"energy_decomposed" :: core_sources__source___profiles_1d___electrons__energy_decomposed
    var"particles" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"particles_decomposed" :: core_sources__source___profiles_1d___electrons__particles_decomposed
    var"particles_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"power_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d___electrons(var"energy"=missing, var"energy_decomposed"=core_sources__source___profiles_1d___electrons__energy_decomposed(), var"particles"=missing, var"particles_decomposed"=core_sources__source___profiles_1d___electrons__particles_decomposed(), var"particles_inside"=missing, var"power_inside"=missing, _parent=WeakRef(missing))
        ids = new(var"energy", var"energy_decomposed", var"particles", var"particles_decomposed", var"particles_inside", var"power_inside", _parent)
        assign_expressions(ids)
        setfield!(ids.energy_decomposed, :_parent, WeakRef(ids))
        setfield!(ids.particles_decomposed, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___profiles_1d <: IDSvectorTimeElement
    var"conductivity_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_parallel_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electrons" :: core_sources__source___profiles_1d___electrons
    var"grid" :: core_sources__source___profiles_1d___grid
    var"ion" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___ion}
    var"j_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"momentum_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"momentum_tor_j_cross_b_field" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"neutral" :: IDSvector{T} where {T<:core_sources__source___profiles_1d___neutral}
    var"time" :: Union{Missing, Real, Function}
    var"torque_tor_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"total_ion_energy" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"total_ion_energy_decomposed" :: core_sources__source___profiles_1d___total_ion_energy_decomposed
    var"total_ion_power_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_sources__source___profiles_1d(var"conductivity_parallel"=missing, var"current_parallel_inside"=missing, var"electrons"=core_sources__source___profiles_1d___electrons(), var"grid"=core_sources__source___profiles_1d___grid(), var"ion"=IDSvector(core_sources__source___profiles_1d___ion[]), var"j_parallel"=missing, var"momentum_tor"=missing, var"momentum_tor_j_cross_b_field"=missing, var"neutral"=IDSvector(core_sources__source___profiles_1d___neutral[]), var"time"=missing, var"torque_tor_inside"=missing, var"total_ion_energy"=missing, var"total_ion_energy_decomposed"=core_sources__source___profiles_1d___total_ion_energy_decomposed(), var"total_ion_power_inside"=missing, _parent=WeakRef(missing))
        ids = new(var"conductivity_parallel", var"current_parallel_inside", var"electrons", var"grid", var"ion", var"j_parallel", var"momentum_tor", var"momentum_tor_j_cross_b_field", var"neutral", var"time", var"torque_tor_inside", var"total_ion_energy", var"total_ion_energy_decomposed", var"total_ion_power_inside", _parent)
        assign_expressions(ids)
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        setfield!(ids.grid, :_parent, WeakRef(ids))
        setfield!(ids.ion, :_parent, WeakRef(ids))
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        setfield!(ids.total_ion_energy_decomposed, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source___identifier <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_sources__source___identifier(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
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

mutable struct core_sources__source___global_quantities <: IDSvectorTimeElement
    var"current_parallel" :: Union{Missing, Real, Function}
    var"electrons" :: core_sources__source___global_quantities___electrons
    var"power" :: Union{Missing, Real, Function}
    var"time" :: Union{Missing, Real, Function}
    var"torque_tor" :: Union{Missing, Real, Function}
    var"total_ion_particles" :: Union{Missing, Real, Function}
    var"total_ion_power" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_sources__source___global_quantities(var"current_parallel"=missing, var"electrons"=core_sources__source___global_quantities___electrons(), var"power"=missing, var"time"=missing, var"torque_tor"=missing, var"total_ion_particles"=missing, var"total_ion_power"=missing, _parent=WeakRef(missing))
        ids = new(var"current_parallel", var"electrons", var"power", var"time", var"torque_tor", var"total_ion_particles", var"total_ion_power", _parent)
        assign_expressions(ids)
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_sources__source <: IDSvectorStaticElement
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

mutable struct core_sources <: IDS
    var"source" :: IDSvector{T} where {T<:core_sources__source}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vacuum_toroidal_field" :: core_sources__vacuum_toroidal_field
    _parent :: WeakRef
    function core_sources(var"source"=IDSvector(core_sources__source[]), var"time"=missing, var"vacuum_toroidal_field"=core_sources__vacuum_toroidal_field(), _parent=WeakRef(missing))
        ids = new(var"source", var"time", var"vacuum_toroidal_field", _parent)
        assign_expressions(ids)
        setfield!(ids.source, :_parent, WeakRef(ids))
        setfield!(ids.vacuum_toroidal_field, :_parent, WeakRef(ids))
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
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___zeff_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___zeff_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___t_i_average_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___t_i_average_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___state___neutral_type <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___state___neutral_type(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___state <: IDSvectorStaticElement
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"label" :: Union{Missing, String, Function}
    var"neutral_type" :: core_profiles__profiles_1d___neutral___state___neutral_type
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___state(var"density"=missing, var"density_fast"=missing, var"density_thermal"=missing, var"electron_configuration"=missing, var"label"=missing, var"neutral_type"=core_profiles__profiles_1d___neutral___state___neutral_type(), var"pressure"=missing, var"pressure_fast_parallel"=missing, var"pressure_fast_perpendicular"=missing, var"pressure_thermal"=missing, var"temperature"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, _parent=WeakRef(missing))
        ids = new(var"density", var"density_fast", var"density_thermal", var"electron_configuration", var"label", var"neutral_type", var"pressure", var"pressure_fast_parallel", var"pressure_fast_perpendicular", var"pressure_thermal", var"temperature", var"vibrational_level", var"vibrational_mode", _parent)
        assign_expressions(ids)
        setfield!(ids.neutral_type, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral___element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral___element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___neutral <: IDSvectorStaticElement
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"element" :: IDSvector{T} where {T<:core_profiles__profiles_1d___neutral___element}
    var"ion_index" :: Union{Missing, Integer, Function}
    var"label" :: Union{Missing, String, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"state" :: IDSvector{T} where {T<:core_profiles__profiles_1d___neutral___state}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___neutral(var"density"=missing, var"density_fast"=missing, var"density_thermal"=missing, var"element"=IDSvector(core_profiles__profiles_1d___neutral___element[]), var"ion_index"=missing, var"label"=missing, var"multiple_states_flag"=missing, var"pressure"=missing, var"pressure_fast_parallel"=missing, var"pressure_fast_perpendicular"=missing, var"pressure_thermal"=missing, var"state"=IDSvector(core_profiles__profiles_1d___neutral___state[]), var"temperature"=missing, _parent=WeakRef(missing))
        ids = new(var"density", var"density_fast", var"density_thermal", var"element", var"ion_index", var"label", var"multiple_states_flag", var"pressure", var"pressure_fast_parallel", var"pressure_fast_perpendicular", var"pressure_thermal", var"state", var"temperature", _parent)
        assign_expressions(ids)
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___velocity <: IDS
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___velocity(var"diamagnetic"=missing, var"parallel"=missing, var"poloidal"=missing, var"radial"=missing, var"toroidal"=missing, _parent=WeakRef(missing))
        ids = new(var"diamagnetic", var"parallel", var"poloidal", var"radial", var"toroidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___temperature_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___temperature_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state___density_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state___density_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___state <: IDSvectorStaticElement
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fit" :: core_profiles__profiles_1d___ion___state___density_fit
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"electron_configuration" :: Union{Missing, String, Function}
    var"ionisation_potential" :: Union{Missing, Real, Function}
    var"label" :: Union{Missing, String, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rotation_frequency_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vibrational_level" :: Union{Missing, Real, Function}
    var"vibrational_mode" :: Union{Missing, String, Function}
    var"z_average" :: Union{Missing, Real, Function}
    var"z_average_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_average_square_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_max" :: Union{Missing, Real, Function}
    var"z_min" :: Union{Missing, Real, Function}
    var"z_square_average" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___state(var"density"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___ion___state___density_fit(), var"density_thermal"=missing, var"electron_configuration"=missing, var"ionisation_potential"=missing, var"label"=missing, var"pressure"=missing, var"pressure_fast_parallel"=missing, var"pressure_fast_perpendicular"=missing, var"pressure_thermal"=missing, var"rotation_frequency_tor"=missing, var"temperature"=missing, var"vibrational_level"=missing, var"vibrational_mode"=missing, var"z_average"=missing, var"z_average_1d"=missing, var"z_average_square_1d"=missing, var"z_max"=missing, var"z_min"=missing, var"z_square_average"=missing, _parent=WeakRef(missing))
        ids = new(var"density", var"density_fast", var"density_fit", var"density_thermal", var"electron_configuration", var"ionisation_potential", var"label", var"pressure", var"pressure_fast_parallel", var"pressure_fast_perpendicular", var"pressure_thermal", var"rotation_frequency_tor", var"temperature", var"vibrational_level", var"vibrational_mode", var"z_average", var"z_average_1d", var"z_average_square_1d", var"z_max", var"z_min", var"z_square_average", _parent)
        assign_expressions(ids)
        setfield!(ids.density_fit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___element <: IDSvectorStaticElement
    var"a" :: Union{Missing, Real, Function}
    var"atoms_n" :: Union{Missing, Integer, Function}
    var"z_n" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___element(var"a"=missing, var"atoms_n"=missing, var"z_n"=missing, _parent=WeakRef(missing))
        ids = new(var"a", var"atoms_n", var"z_n", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion___density_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion___density_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___ion <: IDSvectorStaticElement
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fit" :: core_profiles__profiles_1d___ion___density_fit
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_validity" :: Union{Missing, Integer, Function}
    var"element" :: IDSvector{T} where {T<:core_profiles__profiles_1d___ion___element}
    var"label" :: Union{Missing, String, Function}
    var"multiple_states_flag" :: Union{Missing, Integer, Function}
    var"neutral_index" :: Union{Missing, Integer, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rotation_frequency_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"state" :: IDSvector{T} where {T<:core_profiles__profiles_1d___ion___state}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature_fit" :: core_profiles__profiles_1d___ion___temperature_fit
    var"temperature_validity" :: Union{Missing, Integer, Function}
    var"velocity" :: core_profiles__profiles_1d___ion___velocity
    var"z_ion" :: Union{Missing, Real, Function}
    var"z_ion_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_ion_square_1d" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___ion(var"density"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___ion___density_fit(), var"density_thermal"=missing, var"density_validity"=missing, var"element"=IDSvector(core_profiles__profiles_1d___ion___element[]), var"label"=missing, var"multiple_states_flag"=missing, var"neutral_index"=missing, var"pressure"=missing, var"pressure_fast_parallel"=missing, var"pressure_fast_perpendicular"=missing, var"pressure_thermal"=missing, var"rotation_frequency_tor"=missing, var"state"=IDSvector(core_profiles__profiles_1d___ion___state[]), var"temperature"=missing, var"temperature_fit"=core_profiles__profiles_1d___ion___temperature_fit(), var"temperature_validity"=missing, var"velocity"=core_profiles__profiles_1d___ion___velocity(), var"z_ion"=missing, var"z_ion_1d"=missing, var"z_ion_square_1d"=missing, _parent=WeakRef(missing))
        ids = new(var"density", var"density_fast", var"density_fit", var"density_thermal", var"density_validity", var"element", var"label", var"multiple_states_flag", var"neutral_index", var"pressure", var"pressure_fast_parallel", var"pressure_fast_perpendicular", var"pressure_thermal", var"rotation_frequency_tor", var"state", var"temperature", var"temperature_fit", var"temperature_validity", var"velocity", var"z_ion", var"z_ion_1d", var"z_ion_square_1d", _parent)
        assign_expressions(ids)
        setfield!(ids.density_fit, :_parent, WeakRef(ids))
        setfield!(ids.element, :_parent, WeakRef(ids))
        setfield!(ids.state, :_parent, WeakRef(ids))
        setfield!(ids.temperature_fit, :_parent, WeakRef(ids))
        setfield!(ids.velocity, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___grid <: IDS
    var"area" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"psi_boundary" :: Union{Missing, Real, Function}
    var"psi_magnetic_axis" :: Union{Missing, Real, Function}
    var"rho_pol_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"surface" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"volume" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___grid(var"area"=missing, var"psi"=missing, var"psi_boundary"=missing, var"psi_magnetic_axis"=missing, var"rho_pol_norm"=missing, var"rho_tor"=missing, var"rho_tor_norm"=missing, var"surface"=missing, var"volume"=missing, _parent=WeakRef(missing))
        ids = new(var"area", var"psi", var"psi_boundary", var"psi_magnetic_axis", var"rho_pol_norm", var"rho_tor", var"rho_tor_norm", var"surface", var"volume", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__temperature_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__temperature_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method <: IDS
    var"description" :: Union{Missing, String, Function}
    var"index" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(var"description"=missing, var"index"=missing, var"name"=missing, _parent=WeakRef(missing))
        ids = new(var"description", var"index", var"name", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons__density_fit <: IDS
    var"chi_squared" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"local" :: Union{Missing, AbstractArray{T, 1} where T<:Integer, AbstractRange{T} where T<:Integer, Function}
    var"measured" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parameters" :: Union{Missing, String, Function}
    var"reconstructed" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rho_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"source" :: Union{Missing, AbstractArray{T, 1} where T<:String, AbstractRange{T} where T<:String, Function}
    var"time_measurement" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method
    var"time_measurement_width" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"weight" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons__density_fit(var"chi_squared"=missing, var"local"=missing, var"measured"=missing, var"parameters"=missing, var"reconstructed"=missing, var"rho_tor_norm"=missing, var"source"=missing, var"time_measurement"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(), var"time_measurement_width"=missing, var"weight"=missing, _parent=WeakRef(missing))
        ids = new(var"chi_squared", var"local", var"measured", var"parameters", var"reconstructed", var"rho_tor_norm", var"source", var"time_measurement", var"time_measurement_slice_method", var"time_measurement_width", var"weight", _parent)
        assign_expressions(ids)
        setfield!(ids.time_measurement_slice_method, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___electrons <: IDS
    var"collisionality_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fast" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_fit" :: core_profiles__profiles_1d___electrons__density_fit
    var"density_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"density_validity" :: Union{Missing, Integer, Function}
    var"pressure" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"temperature_fit" :: core_profiles__profiles_1d___electrons__temperature_fit
    var"temperature_validity" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___electrons(var"collisionality_norm"=missing, var"density"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___electrons__density_fit(), var"density_thermal"=missing, var"density_validity"=missing, var"pressure"=missing, var"pressure_fast_parallel"=missing, var"pressure_fast_perpendicular"=missing, var"pressure_thermal"=missing, var"temperature"=missing, var"temperature_fit"=core_profiles__profiles_1d___electrons__temperature_fit(), var"temperature_validity"=missing, _parent=WeakRef(missing))
        ids = new(var"collisionality_norm", var"density", var"density_fast", var"density_fit", var"density_thermal", var"density_validity", var"pressure", var"pressure_fast_parallel", var"pressure_fast_perpendicular", var"pressure_thermal", var"temperature", var"temperature_fit", var"temperature_validity", _parent)
        assign_expressions(ids)
        setfield!(ids.density_fit, :_parent, WeakRef(ids))
        setfield!(ids.temperature_fit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__profiles_1d___e_field <: IDS
    var"diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"poloidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"radial" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"toroidal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__profiles_1d___e_field(var"diamagnetic"=missing, var"parallel"=missing, var"poloidal"=missing, var"radial"=missing, var"toroidal"=missing, _parent=WeakRef(missing))
        ids = new(var"diamagnetic", var"parallel", var"poloidal", var"radial", var"toroidal", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles__profiles_1d <: IDSvectorTimeElement
    var"conductivity_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_parallel_inside" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"e_field" :: core_profiles__profiles_1d___e_field
    var"electrons" :: core_profiles__profiles_1d___electrons
    var"grid" :: core_profiles__profiles_1d___grid
    var"ion" :: IDSvector{T} where {T<:core_profiles__profiles_1d___ion}
    var"j_bootstrap" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_non_inductive" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_ohmic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"j_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"magnetic_shear" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"momentum_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"n_i_thermal_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"n_i_total_over_n_e" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"neutral" :: IDSvector{T} where {T<:core_profiles__profiles_1d___neutral}
    var"phi_potential" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_ion_total" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_parallel" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_perpendicular" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"pressure_thermal" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"q" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"rotation_frequency_tor_sonic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_i_average" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_i_average_fit" :: core_profiles__profiles_1d___t_i_average_fit
    var"time" :: Union{Missing, Real, Function}
    var"zeff" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"zeff_fit" :: core_profiles__profiles_1d___zeff_fit
    _parent :: WeakRef
    function core_profiles__profiles_1d(var"conductivity_parallel"=missing, var"current_parallel_inside"=missing, var"e_field"=core_profiles__profiles_1d___e_field(), var"electrons"=core_profiles__profiles_1d___electrons(), var"grid"=core_profiles__profiles_1d___grid(), var"ion"=IDSvector(core_profiles__profiles_1d___ion[]), var"j_bootstrap"=missing, var"j_non_inductive"=missing, var"j_ohmic"=missing, var"j_tor"=missing, var"j_total"=missing, var"magnetic_shear"=missing, var"momentum_tor"=missing, var"n_i_thermal_total"=missing, var"n_i_total_over_n_e"=missing, var"neutral"=IDSvector(core_profiles__profiles_1d___neutral[]), var"phi_potential"=missing, var"pressure_ion_total"=missing, var"pressure_parallel"=missing, var"pressure_perpendicular"=missing, var"pressure_thermal"=missing, var"q"=missing, var"rotation_frequency_tor_sonic"=missing, var"t_i_average"=missing, var"t_i_average_fit"=core_profiles__profiles_1d___t_i_average_fit(), var"time"=missing, var"zeff"=missing, var"zeff_fit"=core_profiles__profiles_1d___zeff_fit(), _parent=WeakRef(missing))
        ids = new(var"conductivity_parallel", var"current_parallel_inside", var"e_field", var"electrons", var"grid", var"ion", var"j_bootstrap", var"j_non_inductive", var"j_ohmic", var"j_tor", var"j_total", var"magnetic_shear", var"momentum_tor", var"n_i_thermal_total", var"n_i_total_over_n_e", var"neutral", var"phi_potential", var"pressure_ion_total", var"pressure_parallel", var"pressure_perpendicular", var"pressure_thermal", var"q", var"rotation_frequency_tor_sonic", var"t_i_average", var"t_i_average_fit", var"time", var"zeff", var"zeff_fit", _parent)
        assign_expressions(ids)
        setfield!(ids.e_field, :_parent, WeakRef(ids))
        setfield!(ids.electrons, :_parent, WeakRef(ids))
        setfield!(ids.grid, :_parent, WeakRef(ids))
        setfield!(ids.ion, :_parent, WeakRef(ids))
        setfield!(ids.neutral, :_parent, WeakRef(ids))
        setfield!(ids.t_i_average_fit, :_parent, WeakRef(ids))
        setfield!(ids.zeff_fit, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct core_profiles__global_quantities <: IDS
    var"beta_pol" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"beta_tor" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"beta_tor_norm" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_bootstrap" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"current_non_inductive" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ejima" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"energy_diamagnetic" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"ip" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"li_3" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"resistive_psi_losses" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_e_peaking" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"t_i_average_peaking" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"v_loop" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z_eff_resistive" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function core_profiles__global_quantities(var"beta_pol"=missing, var"beta_tor"=missing, var"beta_tor_norm"=missing, var"current_bootstrap"=missing, var"current_non_inductive"=missing, var"ejima"=missing, var"energy_diamagnetic"=missing, var"ip"=missing, var"li_3"=missing, var"resistive_psi_losses"=missing, var"t_e_peaking"=missing, var"t_i_average_peaking"=missing, var"v_loop"=missing, var"z_eff_resistive"=missing, _parent=WeakRef(missing))
        ids = new(var"beta_pol", var"beta_tor", var"beta_tor_norm", var"current_bootstrap", var"current_non_inductive", var"ejima", var"energy_diamagnetic", var"ip", var"li_3", var"resistive_psi_losses", var"t_e_peaking", var"t_i_average_peaking", var"v_loop", var"z_eff_resistive", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct core_profiles <: IDS
    var"global_quantities" :: core_profiles__global_quantities
    var"profiles_1d" :: IDSvector{T} where {T<:core_profiles__profiles_1d}
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"vacuum_toroidal_field" :: core_profiles__vacuum_toroidal_field
    _parent :: WeakRef
    function core_profiles(var"global_quantities"=core_profiles__global_quantities(), var"profiles_1d"=IDSvector(core_profiles__profiles_1d[]), var"time"=missing, var"vacuum_toroidal_field"=core_profiles__vacuum_toroidal_field(), _parent=WeakRef(missing))
        ids = new(var"global_quantities", var"profiles_1d", var"time", var"vacuum_toroidal_field", _parent)
        assign_expressions(ids)
        setfield!(ids.global_quantities, :_parent, WeakRef(ids))
        setfield!(ids.profiles_1d, :_parent, WeakRef(ids))
        setfield!(ids.vacuum_toroidal_field, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct build__tf <: IDS
    var"coils_n" :: Union{Missing, Integer, Function}
    var"thickness" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function build__tf(var"coils_n"=missing, var"thickness"=missing, _parent=WeakRef(missing))
        ids = new(var"coils_n", var"thickness", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct build__pf_coils_rail___outline <: IDS
    var"distance" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function build__pf_coils_rail___outline(var"distance"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"distance", var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct build__pf_coils_rail <: IDSvectorStaticElement
    var"coils_cleareance" :: Union{Missing, Real, Function}
    var"coils_elements_area" :: Union{Missing, Real, Function}
    var"coils_number" :: Union{Missing, Integer, Function}
    var"name" :: Union{Missing, String, Function}
    var"outline" :: build__pf_coils_rail___outline
    _parent :: WeakRef
    function build__pf_coils_rail(var"coils_cleareance"=missing, var"coils_elements_area"=missing, var"coils_number"=missing, var"name"=missing, var"outline"=build__pf_coils_rail___outline(), _parent=WeakRef(missing))
        ids = new(var"coils_cleareance", var"coils_elements_area", var"coils_number", var"name", var"outline", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct build__oh <: IDS
    var"required_b_field" :: Union{Missing, Real, Function}
    var"required_j" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function build__oh(var"required_b_field"=missing, var"required_j"=missing, _parent=WeakRef(missing))
        ids = new(var"required_b_field", var"required_j", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct build__layer___outline <: IDS
    var"r" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"z" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function build__layer___outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        ids = new(var"r", var"z", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct build__layer <: IDSvectorStaticElement
    var"end_radius" :: Union{Missing, Real, Function}
    var"hfs" :: Union{Missing, Integer, Function}
    var"identifier" :: Union{Missing, Integer, Function}
    var"material" :: Union{Missing, String, Function}
    var"name" :: Union{Missing, String, Function}
    var"outline" :: build__layer___outline
    var"shape" :: Union{Missing, Integer, Function}
    var"shape_parameters" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    var"start_radius" :: Union{Missing, Real, Function}
    var"thickness" :: Union{Missing, Real, Function}
    var"type" :: Union{Missing, Integer, Function}
    _parent :: WeakRef
    function build__layer(var"end_radius"=missing, var"hfs"=missing, var"identifier"=missing, var"material"=missing, var"name"=missing, var"outline"=build__layer___outline(), var"shape"=missing, var"shape_parameters"=missing, var"start_radius"=missing, var"thickness"=missing, var"type"=missing, _parent=WeakRef(missing))
        ids = new(var"end_radius", var"hfs", var"identifier", var"material", var"name", var"outline", var"shape", var"shape_parameters", var"start_radius", var"thickness", var"type", _parent)
        assign_expressions(ids)
        setfield!(ids.outline, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct build__flux_swing_requirements <: IDS
    var"flattop" :: Union{Missing, Real, Function}
    var"pf" :: Union{Missing, Real, Function}
    var"rampup" :: Union{Missing, Real, Function}
    _parent :: WeakRef
    function build__flux_swing_requirements(var"flattop"=missing, var"pf"=missing, var"rampup"=missing, _parent=WeakRef(missing))
        ids = new(var"flattop", var"pf", var"rampup", _parent)
        assign_expressions(ids)
        return ids
    end
end

mutable struct build <: IDS
    var"flux_swing_requirements" :: build__flux_swing_requirements
    var"layer" :: IDSvector{T} where {T<:build__layer}
    var"oh" :: build__oh
    var"pf_coils_rail" :: IDSvector{T} where {T<:build__pf_coils_rail}
    var"tf" :: build__tf
    var"time" :: Union{Missing, AbstractArray{T, 1} where T<:Real, AbstractRange{T} where T<:Real, Function}
    _parent :: WeakRef
    function build(var"flux_swing_requirements"=build__flux_swing_requirements(), var"layer"=IDSvector(build__layer[]), var"oh"=build__oh(), var"pf_coils_rail"=IDSvector(build__pf_coils_rail[]), var"tf"=build__tf(), var"time"=missing, _parent=WeakRef(missing))
        ids = new(var"flux_swing_requirements", var"layer", var"oh", var"pf_coils_rail", var"tf", var"time", _parent)
        assign_expressions(ids)
        setfield!(ids.flux_swing_requirements, :_parent, WeakRef(ids))
        setfield!(ids.layer, :_parent, WeakRef(ids))
        setfield!(ids.oh, :_parent, WeakRef(ids))
        setfield!(ids.pf_coils_rail, :_parent, WeakRef(ids))
        setfield!(ids.tf, :_parent, WeakRef(ids))
        return ids
    end
end

mutable struct dd <: IDS
    var"build" :: Union{Missing, build}
    var"core_profiles" :: Union{Missing, core_profiles}
    var"core_sources" :: Union{Missing, core_sources}
    var"dataset_description" :: Union{Missing, dataset_description}
    var"equilibrium" :: Union{Missing, equilibrium}
    var"pf_active" :: Union{Missing, pf_active}
    var"summary" :: Union{Missing, summary}
    var"tf" :: Union{Missing, tf}
    var"wall" :: Union{Missing, wall}
    var"global_time" :: Real
    _parent :: WeakRef
    function dd(var"build"=build(), var"core_profiles"=core_profiles(), var"core_sources"=core_sources(), var"dataset_description"=dataset_description(), var"equilibrium"=equilibrium(), var"pf_active"=pf_active(), var"summary"=summary(), var"tf"=tf(), var"wall"=wall(), var"global_time"=0.0, _parent=WeakRef(missing))
        ids = new(var"build", var"core_profiles", var"core_sources", var"dataset_description", var"equilibrium", var"pf_active", var"summary", var"tf", var"wall", var"global_time", _parent)
        assign_expressions(ids)
        setfield!(ids.build, :_parent, WeakRef(ids))
        setfield!(ids.core_profiles, :_parent, WeakRef(ids))
        setfield!(ids.core_sources, :_parent, WeakRef(ids))
        setfield!(ids.dataset_description, :_parent, WeakRef(ids))
        setfield!(ids.equilibrium, :_parent, WeakRef(ids))
        setfield!(ids.pf_active, :_parent, WeakRef(ids))
        setfield!(ids.summary, :_parent, WeakRef(ids))
        setfield!(ids.tf, :_parent, WeakRef(ids))
        setfield!(ids.wall, :_parent, WeakRef(ids))
        return ids
    end
end

