include("functionarrays.jl")

convertsion_types = Union{Float64, Int64, String, AbstractFDVector{Float64}, AbstractFDVector{Int64}, AbstractFDVector{String}, Array{Float64, N} where N}

Base.@kwdef mutable struct wall__temperature_reference <: FDS
    var"data" :: Union{Missing, Float64} = missing
    var"description" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__ids_properties <: FDS
    var"provider" :: Union{Missing, String} = missing
    var"version_put" :: wall__ids_properties__version_put = wall__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Int64} = missing
    var"source" :: Union{Missing, String} = missing
    var"creation_date" :: Union{Missing, String} = missing
    var"comment" :: Union{Missing, String} = missing
    var"occurrence" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__ids_properties(var"provider"=missing, var"version_put"=wall__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        setfield!(obj.version_put, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral___element <: FDS
    var"atoms_n" :: Union{Missing, Int64} = missing
    var"z_n" :: Union{Missing, Float64} = missing
    var"multiplicity" :: Union{Missing, Float64} = missing
    var"a" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__global_quantities__neutral <: FDS
    var"label" :: Union{Missing, String} = missing
    var"sputtering_chemical_coefficient" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"gas_puff" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"recycling_particles_coefficient" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"pumping_speed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"particle_flux_from_wall" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"recycling_energy_coefficient" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"wall_inventory" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"particle_flux_from_plasma" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"sputtering_physical_coefficient" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"element" :: FDSvector{T} where {T<:wall__global_quantities__neutral___element} = FDSvector(wall__global_quantities__neutral___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities__neutral(var"label"=missing, var"sputtering_chemical_coefficient"=missing, var"gas_puff"=missing, var"recycling_particles_coefficient"=missing, var"pumping_speed"=missing, var"particle_flux_from_wall"=missing, var"recycling_energy_coefficient"=missing, var"wall_inventory"=missing, var"particle_flux_from_plasma"=missing, var"sputtering_physical_coefficient"=missing, var"element"=FDSvector(wall__global_quantities__neutral___element[]), _parent=WeakRef(missing))
        obj = new(var"label", var"sputtering_chemical_coefficient", var"gas_puff", var"recycling_particles_coefficient", var"pumping_speed", var"particle_flux_from_wall", var"recycling_energy_coefficient", var"wall_inventory", var"particle_flux_from_plasma", var"sputtering_physical_coefficient", var"element", _parent)
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__global_quantities__electrons <: FDS
    var"particle_flux_from_plasma" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gas_puff" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_outer_target" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pumping_speed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"particle_flux_from_wall" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"power_inner_target" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__global_quantities <: FDS
    var"neutral" :: FDSvector{T} where {T<:wall__global_quantities__neutral} = FDSvector(wall__global_quantities__neutral[])
    var"power_incident" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_radiated" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_inner_target_ion_total" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"temperature" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_conducted" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_convected" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"current_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"electrons" :: wall__global_quantities__electrons = wall__global_quantities__electrons()
    var"power_density_inner_target_max" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_black_body" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_recombination_neutrals" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_to_cooling" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_density_outer_target_max" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_recombination_plasma" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_currents" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"power_neutrals" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities(var"neutral"=FDSvector(wall__global_quantities__neutral[]), var"power_incident"=missing, var"power_radiated"=missing, var"power_inner_target_ion_total"=missing, var"temperature"=missing, var"power_conducted"=missing, var"power_convected"=missing, var"current_tor"=missing, var"electrons"=wall__global_quantities__electrons(), var"power_density_inner_target_max"=missing, var"power_black_body"=missing, var"power_recombination_neutrals"=missing, var"power_to_cooling"=missing, var"power_density_outer_target_max"=missing, var"power_recombination_plasma"=missing, var"power_currents"=missing, var"power_neutrals"=missing, _parent=WeakRef(missing))
        obj = new(var"neutral", var"power_incident", var"power_radiated", var"power_inner_target_ion_total", var"temperature", var"power_conducted", var"power_convected", var"current_tor", var"electrons", var"power_density_inner_target_max", var"power_black_body", var"power_recombination_neutrals", var"power_to_cooling", var"power_density_outer_target_max", var"power_recombination_plasma", var"power_currents", var"power_neutrals", _parent)
        setfield!(obj.neutral, :_parent, WeakRef(obj))
        setfield!(obj.electrons, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__first_wall_power_flux_peak <: FDS
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"data" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary <: FDS
    var"neighbours" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object <: FDS
    var"nodes" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"measure" :: Union{Missing, Float64} = missing
    var"geometry" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"boundary" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        setfield!(obj.boundary, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension(var"object"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        obj = new(var"object", _parent)
        setfield!(obj.object, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___geometry_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space <: FDS
    var"coordinates_type" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"geometry_type" :: wall__description_ggd___grid_ggd___space___geometry_type = wall__description_ggd___grid_ggd___space___geometry_type()
    var"identifier" :: wall__description_ggd___grid_ggd___space___identifier = wall__description_ggd___grid_ggd___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space(var"coordinates_type"=missing, var"geometry_type"=wall__description_ggd___grid_ggd___space___geometry_type(), var"identifier"=wall__description_ggd___grid_ggd___space___identifier(), var"objects_per_dimension"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[]), _parent=WeakRef(missing))
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        setfield!(obj.geometry_type, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        setfield!(obj.objects_per_dimension, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element___object <: FDS
    var"dimension" :: Union{Missing, Int64} = missing
    var"space" :: Union{Missing, Int64} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element <: FDS
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element___object} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___element(var"object"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[]), _parent=WeakRef(missing))
        obj = new(var"object", _parent)
        setfield!(obj.object, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___base <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___base} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[])
    var"metric" :: wall__description_ggd___grid_ggd___grid_subset___metric = wall__description_ggd___grid_ggd___grid_subset___metric()
    var"dimension" :: Union{Missing, Int64} = missing
    var"identifier" :: wall__description_ggd___grid_ggd___grid_subset___identifier = wall__description_ggd___grid_ggd___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset(var"base"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[]), var"metric"=wall__description_ggd___grid_ggd___grid_subset___metric(), var"dimension"=missing, var"identifier"=wall__description_ggd___grid_ggd___grid_subset___identifier(), var"element"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[]), _parent=WeakRef(missing))
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        setfield!(obj.base, :_parent, WeakRef(obj))
        setfield!(obj.metric, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd <: FDS
    var"time" :: Union{Missing, Float64} = missing
    var"grid_subset" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset} = FDSvector(wall__description_ggd___grid_ggd___grid_subset[])
    var"space" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space} = FDSvector(wall__description_ggd___grid_ggd___space[])
    var"identifier" :: wall__description_ggd___grid_ggd___identifier = wall__description_ggd___grid_ggd___identifier()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd(var"time"=missing, var"grid_subset"=FDSvector(wall__description_ggd___grid_ggd___grid_subset[]), var"space"=FDSvector(wall__description_ggd___grid_ggd___space[]), var"identifier"=wall__description_ggd___grid_ggd___identifier(), _parent=WeakRef(missing))
        obj = new(var"time", var"grid_subset", var"space", var"identifier", _parent)
        setfield!(obj.grid_subset, :_parent, WeakRef(obj))
        setfield!(obj.space, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___temperature <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___ggd___power_density <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_ggd___ggd <: FDS
    var"temperature" :: FDSvector{T} where {T<:wall__description_ggd___ggd___temperature} = FDSvector(wall__description_ggd___ggd___temperature[])
    var"time" :: Union{Missing, Float64} = missing
    var"power_density" :: FDSvector{T} where {T<:wall__description_ggd___ggd___power_density} = FDSvector(wall__description_ggd___ggd___power_density[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd(var"temperature"=FDSvector(wall__description_ggd___ggd___temperature[]), var"time"=missing, var"power_density"=FDSvector(wall__description_ggd___ggd___power_density[]), _parent=WeakRef(missing))
        obj = new(var"temperature", var"time", var"power_density", _parent)
        setfield!(obj.temperature, :_parent, WeakRef(obj))
        setfield!(obj.power_density, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd <: FDS
    var"grid_ggd" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd} = FDSvector(wall__description_ggd___grid_ggd[])
    var"type" :: wall__description_ggd___type = wall__description_ggd___type()
    var"ggd" :: FDSvector{T} where {T<:wall__description_ggd___ggd} = FDSvector(wall__description_ggd___ggd[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd(var"grid_ggd"=FDSvector(wall__description_ggd___grid_ggd[]), var"type"=wall__description_ggd___type(), var"ggd"=FDSvector(wall__description_ggd___ggd[]), _parent=WeakRef(missing))
        obj = new(var"grid_ggd", var"type", var"ggd", _parent)
        setfield!(obj.grid_ggd, :_parent, WeakRef(obj))
        setfield!(obj.type, :_parent, WeakRef(obj))
        setfield!(obj.ggd, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___outline <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___j_tor <: FDS
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"data" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element <: FDS
    var"name" :: Union{Missing, String} = missing
    var"j_tor" :: wall__description_2d___vessel__unit___element___j_tor = wall__description_2d___vessel__unit___element___j_tor()
    var"resistivity" :: Union{Missing, Float64} = missing
    var"outline" :: wall__description_2d___vessel__unit___element___outline = wall__description_2d___vessel__unit___element___outline()
    var"resistance" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element(var"name"=missing, var"j_tor"=wall__description_2d___vessel__unit___element___j_tor(), var"resistivity"=missing, var"outline"=wall__description_2d___vessel__unit___element___outline(), var"resistance"=missing, _parent=WeakRef(missing))
        obj = new(var"name", var"j_tor", var"resistivity", var"outline", var"resistance", _parent)
        setfield!(obj.j_tor, :_parent, WeakRef(obj))
        setfield!(obj.outline, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_outer <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_inner <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__centreline <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular <: FDS
    var"outline_inner" :: wall__description_2d___vessel__unit___annular__outline_inner = wall__description_2d___vessel__unit___annular__outline_inner()
    var"centreline" :: wall__description_2d___vessel__unit___annular__centreline = wall__description_2d___vessel__unit___annular__centreline()
    var"thickness" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"resistivity" :: Union{Missing, Float64} = missing
    var"outline_outer" :: wall__description_2d___vessel__unit___annular__outline_outer = wall__description_2d___vessel__unit___annular__outline_outer()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular(var"outline_inner"=wall__description_2d___vessel__unit___annular__outline_inner(), var"centreline"=wall__description_2d___vessel__unit___annular__centreline(), var"thickness"=missing, var"resistivity"=missing, var"outline_outer"=wall__description_2d___vessel__unit___annular__outline_outer(), _parent=WeakRef(missing))
        obj = new(var"outline_inner", var"centreline", var"thickness", var"resistivity", var"outline_outer", _parent)
        setfield!(obj.outline_inner, :_parent, WeakRef(obj))
        setfield!(obj.centreline, :_parent, WeakRef(obj))
        setfield!(obj.outline_outer, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit <: FDS
    var"name" :: Union{Missing, String} = missing
    var"annular" :: wall__description_2d___vessel__unit___annular = wall__description_2d___vessel__unit___annular()
    var"identifier" :: Union{Missing, String} = missing
    var"element" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit___element} = FDSvector(wall__description_2d___vessel__unit___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit(var"name"=missing, var"annular"=wall__description_2d___vessel__unit___annular(), var"identifier"=missing, var"element"=FDSvector(wall__description_2d___vessel__unit___element[]), _parent=WeakRef(missing))
        obj = new(var"name", var"annular", var"identifier", var"element", _parent)
        setfield!(obj.annular, :_parent, WeakRef(obj))
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___vessel <: FDS
    var"type" :: wall__description_2d___vessel__type = wall__description_2d___vessel__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit} = FDSvector(wall__description_2d___vessel__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel(var"type"=wall__description_2d___vessel__type(), var"unit"=FDSvector(wall__description_2d___vessel__unit[]), _parent=WeakRef(missing))
        obj = new(var"type", var"unit", _parent)
        setfield!(obj.type, :_parent, WeakRef(obj))
        setfield!(obj.unit, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit___outline <: FDS
    var"time" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit <: FDS
    var"name" :: Union{Missing, String} = missing
    var"resistivity" :: Union{Missing, Float64} = missing
    var"closed" :: Union{Missing, Int64} = missing
    var"phi_extensions" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"outline" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit___outline} = FDSvector(wall__description_2d___mobile__unit___outline[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"phi_extensions"=missing, var"outline"=FDSvector(wall__description_2d___mobile__unit___outline[]), _parent=WeakRef(missing))
        obj = new(var"name", var"resistivity", var"closed", var"phi_extensions", var"outline", _parent)
        setfield!(obj.outline, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___mobile <: FDS
    var"type" :: wall__description_2d___mobile__type = wall__description_2d___mobile__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit} = FDSvector(wall__description_2d___mobile__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile(var"type"=wall__description_2d___mobile__type(), var"unit"=FDSvector(wall__description_2d___mobile__unit[]), _parent=WeakRef(missing))
        obj = new(var"type", var"unit", _parent)
        setfield!(obj.type, :_parent, WeakRef(obj))
        setfield!(obj.unit, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit___outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit <: FDS
    var"name" :: Union{Missing, String} = missing
    var"resistivity" :: Union{Missing, Float64} = missing
    var"closed" :: Union{Missing, Int64} = missing
    var"outline" :: wall__description_2d___limiter__unit___outline = wall__description_2d___limiter__unit___outline()
    var"phi_extensions" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"outline"=wall__description_2d___limiter__unit___outline(), var"phi_extensions"=missing, _parent=WeakRef(missing))
        obj = new(var"name", var"resistivity", var"closed", var"outline", var"phi_extensions", _parent)
        setfield!(obj.outline, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__description_2d___limiter <: FDS
    var"type" :: wall__description_2d___limiter__type = wall__description_2d___limiter__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___limiter__unit} = FDSvector(wall__description_2d___limiter__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter(var"type"=wall__description_2d___limiter__type(), var"unit"=FDSvector(wall__description_2d___limiter__unit[]), _parent=WeakRef(missing))
        obj = new(var"type", var"unit", _parent)
        setfield!(obj.type, :_parent, WeakRef(obj))
        setfield!(obj.unit, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d <: FDS
    var"mobile" :: wall__description_2d___mobile = wall__description_2d___mobile()
    var"limiter" :: wall__description_2d___limiter = wall__description_2d___limiter()
    var"type" :: wall__description_2d___type = wall__description_2d___type()
    var"vessel" :: wall__description_2d___vessel = wall__description_2d___vessel()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d(var"mobile"=wall__description_2d___mobile(), var"limiter"=wall__description_2d___limiter(), var"type"=wall__description_2d___type(), var"vessel"=wall__description_2d___vessel(), _parent=WeakRef(missing))
        obj = new(var"mobile", var"limiter", var"type", var"vessel", _parent)
        setfield!(obj.mobile, :_parent, WeakRef(obj))
        setfield!(obj.limiter, :_parent, WeakRef(obj))
        setfield!(obj.type, :_parent, WeakRef(obj))
        setfield!(obj.vessel, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall__code__library <: FDS
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct wall__code <: FDS
    var"library" :: FDSvector{T} where {T<:wall__code__library} = FDSvector(wall__code__library[])
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"output_flag" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__code(var"library"=FDSvector(wall__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        setfield!(obj.library, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct wall <: FDS
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ids_properties" :: wall__ids_properties = wall__ids_properties()
    var"description_ggd" :: FDSvector{T} where {T<:wall__description_ggd} = FDSvector(wall__description_ggd[])
    var"description_2d" :: FDSvector{T} where {T<:wall__description_2d} = FDSvector(wall__description_2d[])
    var"first_wall_surface_area" :: Union{Missing, Float64} = missing
    var"code" :: wall__code = wall__code()
    var"global_quantities" :: wall__global_quantities = wall__global_quantities()
    var"temperature_reference" :: wall__temperature_reference = wall__temperature_reference()
    var"first_wall_power_flux_peak" :: wall__first_wall_power_flux_peak = wall__first_wall_power_flux_peak()
    _parent :: WeakRef = WeakRef(missing)
    function wall(var"time"=missing, var"ids_properties"=wall__ids_properties(), var"description_ggd"=FDSvector(wall__description_ggd[]), var"description_2d"=FDSvector(wall__description_2d[]), var"first_wall_surface_area"=missing, var"code"=wall__code(), var"global_quantities"=wall__global_quantities(), var"temperature_reference"=wall__temperature_reference(), var"first_wall_power_flux_peak"=wall__first_wall_power_flux_peak(), _parent=WeakRef(missing))
        obj = new(var"time", var"ids_properties", var"description_ggd", var"description_2d", var"first_wall_surface_area", var"code", var"global_quantities", var"temperature_reference", var"first_wall_power_flux_peak", _parent)
        setfield!(obj.ids_properties, :_parent, WeakRef(obj))
        setfield!(obj.description_ggd, :_parent, WeakRef(obj))
        setfield!(obj.description_2d, :_parent, WeakRef(obj))
        setfield!(obj.code, :_parent, WeakRef(obj))
        setfield!(obj.global_quantities, :_parent, WeakRef(obj))
        setfield!(obj.temperature_reference, :_parent, WeakRef(obj))
        setfield!(obj.first_wall_power_flux_peak, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"r0" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid <: FDS
    var"volume_element" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"dim2" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dim1" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d <: FDS
    var"psi" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"b_field_r" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"r" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"b_r" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"theta" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"b_field_z" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"j_tor" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"phi" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"z" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"b_field_tor" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"b_z" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"grid" :: equilibrium__time_slice___profiles_2d___grid = equilibrium__time_slice___profiles_2d___grid()
    var"grid_type" :: equilibrium__time_slice___profiles_2d___grid_type = equilibrium__time_slice___profiles_2d___grid_type()
    var"j_parallel" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"b_tor" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_2d(var"psi"=missing, var"b_field_r"=missing, var"r"=missing, var"b_r"=missing, var"theta"=missing, var"b_field_z"=missing, var"j_tor"=missing, var"phi"=missing, var"z"=missing, var"b_field_tor"=missing, var"b_z"=missing, var"grid"=equilibrium__time_slice___profiles_2d___grid(), var"grid_type"=equilibrium__time_slice___profiles_2d___grid_type(), var"j_parallel"=missing, var"b_tor"=missing, _parent=WeakRef(missing))
        obj = new(var"psi", var"b_field_r", var"r", var"b_r", var"theta", var"b_field_z", var"j_tor", var"phi", var"z", var"b_field_tor", var"b_z", var"grid", var"grid_type", var"j_parallel", var"b_tor", _parent)
        setfield!(obj.grid, :_parent, WeakRef(obj))
        setfield!(obj.grid_type, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d__geometric_axis <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d <: FDS
    var"b_field_max" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dvolume_drho_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm9" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dpsi_drho_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"surface" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"magnetic_shear" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"b_average" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"b_field_min" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"darea_dpsi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm3" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"squareness_upper_inner" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"squareness_lower_inner" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"elongation" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"beta_pol" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"b_field_average" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"j_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm6" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"psi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm8" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dpressure_dpsi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"triangularity_upper" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"darea_drho_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"area" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"trapped_fraction" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"volume" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dvolume_dpsi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"b_min" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"f" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"mass_density" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"r_outboard" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm4" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"phi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"squareness_lower_outer" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"triangularity_lower" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm2" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_volume_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm1" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm5" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"b_max" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"f_df_dpsi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"j_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"r_inboard" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"q" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"gm7" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"squareness_upper_outer" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"geometric_axis" :: equilibrium__time_slice___profiles_1d__geometric_axis = equilibrium__time_slice___profiles_1d__geometric_axis()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_1d(var"b_field_max"=missing, var"dvolume_drho_tor"=missing, var"gm9"=missing, var"dpsi_drho_tor"=missing, var"surface"=missing, var"rho_tor"=missing, var"magnetic_shear"=missing, var"b_average"=missing, var"b_field_min"=missing, var"darea_dpsi"=missing, var"gm3"=missing, var"squareness_upper_inner"=missing, var"squareness_lower_inner"=missing, var"rho_tor_norm"=missing, var"elongation"=missing, var"beta_pol"=missing, var"b_field_average"=missing, var"j_parallel"=missing, var"gm6"=missing, var"psi"=missing, var"gm8"=missing, var"dpressure_dpsi"=missing, var"triangularity_upper"=missing, var"darea_drho_tor"=missing, var"area"=missing, var"trapped_fraction"=missing, var"volume"=missing, var"dvolume_dpsi"=missing, var"b_min"=missing, var"f"=missing, var"mass_density"=missing, var"r_outboard"=missing, var"gm4"=missing, var"phi"=missing, var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"gm2"=missing, var"rho_volume_norm"=missing, var"gm1"=missing, var"gm5"=missing, var"b_max"=missing, var"f_df_dpsi"=missing, var"j_tor"=missing, var"r_inboard"=missing, var"q"=missing, var"gm7"=missing, var"pressure"=missing, var"squareness_upper_outer"=missing, var"geometric_axis"=equilibrium__time_slice___profiles_1d__geometric_axis(), _parent=WeakRef(missing))
        obj = new(var"b_field_max", var"dvolume_drho_tor", var"gm9", var"dpsi_drho_tor", var"surface", var"rho_tor", var"magnetic_shear", var"b_average", var"b_field_min", var"darea_dpsi", var"gm3", var"squareness_upper_inner", var"squareness_lower_inner", var"rho_tor_norm", var"elongation", var"beta_pol", var"b_field_average", var"j_parallel", var"gm6", var"psi", var"gm8", var"dpressure_dpsi", var"triangularity_upper", var"darea_drho_tor", var"area", var"trapped_fraction", var"volume", var"dvolume_dpsi", var"b_min", var"f", var"mass_density", var"r_outboard", var"gm4", var"phi", var"squareness_lower_outer", var"triangularity_lower", var"gm2", var"rho_volume_norm", var"gm1", var"gm5", var"b_max", var"f_df_dpsi", var"j_tor", var"r_inboard", var"q", var"gm7", var"pressure", var"squareness_upper_outer", var"geometric_axis", _parent)
        setfield!(obj.geometric_axis, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__q_min <: FDS
    var"value" :: Union{Missing, Float64} = missing
    var"rho_tor_norm" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__magnetic_axis <: FDS
    var"b_field_tor" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    var"b_tor" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__current_centre <: FDS
    var"velocity_z" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities <: FDS
    var"ip" :: Union{Missing, Float64} = missing
    var"li_3" :: Union{Missing, Float64} = missing
    var"beta_tor" :: Union{Missing, Float64} = missing
    var"surface" :: Union{Missing, Float64} = missing
    var"magnetic_axis" :: equilibrium__time_slice___global_quantities__magnetic_axis = equilibrium__time_slice___global_quantities__magnetic_axis()
    var"energy_mhd" :: Union{Missing, Float64} = missing
    var"psi_boundary" :: Union{Missing, Float64} = missing
    var"length_pol" :: Union{Missing, Float64} = missing
    var"area" :: Union{Missing, Float64} = missing
    var"psi_external_average" :: Union{Missing, Float64} = missing
    var"q_95" :: Union{Missing, Float64} = missing
    var"q_axis" :: Union{Missing, Float64} = missing
    var"psi_axis" :: Union{Missing, Float64} = missing
    var"w_mhd" :: Union{Missing, Float64} = missing
    var"volume" :: Union{Missing, Float64} = missing
    var"plasma_inductance" :: Union{Missing, Float64} = missing
    var"beta_pol" :: Union{Missing, Float64} = missing
    var"beta_normal" :: Union{Missing, Float64} = missing
    var"current_centre" :: equilibrium__time_slice___global_quantities__current_centre = equilibrium__time_slice___global_quantities__current_centre()
    var"q_min" :: equilibrium__time_slice___global_quantities__q_min = equilibrium__time_slice___global_quantities__q_min()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities(var"ip"=missing, var"li_3"=missing, var"beta_tor"=missing, var"surface"=missing, var"magnetic_axis"=equilibrium__time_slice___global_quantities__magnetic_axis(), var"energy_mhd"=missing, var"psi_boundary"=missing, var"length_pol"=missing, var"area"=missing, var"psi_external_average"=missing, var"q_95"=missing, var"q_axis"=missing, var"psi_axis"=missing, var"w_mhd"=missing, var"volume"=missing, var"plasma_inductance"=missing, var"beta_pol"=missing, var"beta_normal"=missing, var"current_centre"=equilibrium__time_slice___global_quantities__current_centre(), var"q_min"=equilibrium__time_slice___global_quantities__q_min(), _parent=WeakRef(missing))
        obj = new(var"ip", var"li_3", var"beta_tor", var"surface", var"magnetic_axis", var"energy_mhd", var"psi_boundary", var"length_pol", var"area", var"psi_external_average", var"q_95", var"q_axis", var"psi_axis", var"w_mhd", var"volume", var"plasma_inductance", var"beta_pol", var"beta_normal", var"current_centre", var"q_min", _parent)
        setfield!(obj.magnetic_axis, :_parent, WeakRef(obj))
        setfield!(obj.current_centre, :_parent, WeakRef(obj))
        setfield!(obj.q_min, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___z <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___theta <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___r <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___psi <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___phi <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_tor <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_parallel <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary <: FDS
    var"neighbours" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object <: FDS
    var"nodes" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"measure" :: Union{Missing, Float64} = missing
    var"geometry" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        setfield!(obj.boundary, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        obj = new(var"object", _parent)
        setfield!(obj.object, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___geometry_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space <: FDS
    var"coordinates_type" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"geometry_type" :: equilibrium__time_slice___ggd___grid__space___geometry_type = equilibrium__time_slice___ggd___grid__space___geometry_type()
    var"identifier" :: equilibrium__time_slice___ggd___grid__space___identifier = equilibrium__time_slice___ggd___grid__space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__time_slice___ggd___grid__space___geometry_type(), var"identifier"=equilibrium__time_slice___ggd___grid__space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[]), _parent=WeakRef(missing))
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        setfield!(obj.geometry_type, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        setfield!(obj.objects_per_dimension, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element___object <: FDS
    var"dimension" :: Union{Missing, Int64} = missing
    var"space" :: Union{Missing, Int64} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element___object} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[]), _parent=WeakRef(missing))
        obj = new(var"object", _parent)
        setfield!(obj.object, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___base <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___base} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[])
    var"metric" :: equilibrium__time_slice___ggd___grid__grid_subset___metric = equilibrium__time_slice___ggd___grid__grid_subset___metric()
    var"dimension" :: Union{Missing, Int64} = missing
    var"identifier" :: equilibrium__time_slice___ggd___grid__grid_subset___identifier = equilibrium__time_slice___ggd___grid__grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset(var"base"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[]), var"metric"=equilibrium__time_slice___ggd___grid__grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__time_slice___ggd___grid__grid_subset___identifier(), var"element"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[]), _parent=WeakRef(missing))
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        setfield!(obj.base, :_parent, WeakRef(obj))
        setfield!(obj.metric, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space} = FDSvector(equilibrium__time_slice___ggd___grid__space[])
    var"identifier" :: equilibrium__time_slice___ggd___grid__identifier = equilibrium__time_slice___ggd___grid__identifier()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid(var"grid_subset"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset[]), var"space"=FDSvector(equilibrium__time_slice___ggd___grid__space[]), var"identifier"=equilibrium__time_slice___ggd___grid__identifier(), _parent=WeakRef(missing))
        obj = new(var"grid_subset", var"space", var"identifier", _parent)
        setfield!(obj.grid_subset, :_parent, WeakRef(obj))
        setfield!(obj.space, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_z <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_tor <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_r <: FDS
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd <: FDS
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
        obj = new(var"b_field_z", var"psi", var"theta", var"z", var"phi", var"j_tor", var"grid", var"b_field_tor", var"b_field_r", var"r", var"j_parallel", _parent)
        setfield!(obj.b_field_z, :_parent, WeakRef(obj))
        setfield!(obj.psi, :_parent, WeakRef(obj))
        setfield!(obj.theta, :_parent, WeakRef(obj))
        setfield!(obj.z, :_parent, WeakRef(obj))
        setfield!(obj.phi, :_parent, WeakRef(obj))
        setfield!(obj.j_tor, :_parent, WeakRef(obj))
        setfield!(obj.grid, :_parent, WeakRef(obj))
        setfield!(obj.b_field_tor, :_parent, WeakRef(obj))
        setfield!(obj.b_field_r, :_parent, WeakRef(obj))
        setfield!(obj.r, :_parent, WeakRef(obj))
        setfield!(obj.j_parallel, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid <: FDS
    var"volume_element" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"dim2" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dim1" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system <: FDS
    var"jacobian" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g13_covariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g11_contravariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g13_contravariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"r" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g12_contravariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g22_contravariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"z" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g33_contravariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g22_covariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 4}} = missing
    var"g12_covariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g33_covariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"grid" :: equilibrium__time_slice___coordinate_system__grid = equilibrium__time_slice___coordinate_system__grid()
    var"g23_covariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"g11_covariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 4}} = missing
    var"g23_contravariant" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"grid_type" :: equilibrium__time_slice___coordinate_system__grid_type = equilibrium__time_slice___coordinate_system__grid_type()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___coordinate_system(var"jacobian"=missing, var"g13_covariant"=missing, var"g11_contravariant"=missing, var"g13_contravariant"=missing, var"r"=missing, var"g12_contravariant"=missing, var"g22_contravariant"=missing, var"z"=missing, var"g33_contravariant"=missing, var"g22_covariant"=missing, var"tensor_contravariant"=missing, var"g12_covariant"=missing, var"g33_covariant"=missing, var"grid"=equilibrium__time_slice___coordinate_system__grid(), var"g23_covariant"=missing, var"g11_covariant"=missing, var"tensor_covariant"=missing, var"g23_contravariant"=missing, var"grid_type"=equilibrium__time_slice___coordinate_system__grid_type(), _parent=WeakRef(missing))
        obj = new(var"jacobian", var"g13_covariant", var"g11_contravariant", var"g13_contravariant", var"r", var"g12_contravariant", var"g22_contravariant", var"z", var"g33_contravariant", var"g22_covariant", var"tensor_contravariant", var"g12_covariant", var"g33_covariant", var"grid", var"g23_covariant", var"g11_covariant", var"tensor_covariant", var"g23_contravariant", var"grid_type", _parent)
        setfield!(obj.grid, :_parent, WeakRef(obj))
        setfield!(obj.grid_type, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___convergence <: FDS
    var"iterations_n" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_reconstructed <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_measured <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point <: FDS
    var"chi_squared_z" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"position_measured" :: equilibrium__time_slice___constraints__x_point___position_measured = equilibrium__time_slice___constraints__x_point___position_measured()
    var"time_measurement" :: Union{Missing, Float64} = missing
    var"chi_squared_r" :: Union{Missing, Float64} = missing
    var"position_reconstructed" :: equilibrium__time_slice___constraints__x_point___position_reconstructed = equilibrium__time_slice___constraints__x_point___position_reconstructed()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__x_point(var"chi_squared_z"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"position_measured"=equilibrium__time_slice___constraints__x_point___position_measured(), var"time_measurement"=missing, var"chi_squared_r"=missing, var"position_reconstructed"=equilibrium__time_slice___constraints__x_point___position_reconstructed(), _parent=WeakRef(missing))
        obj = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        setfield!(obj.position_measured, :_parent, WeakRef(obj))
        setfield!(obj.position_reconstructed, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_reconstructed <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_measured <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point <: FDS
    var"chi_squared_z" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"position_measured" :: equilibrium__time_slice___constraints__strike_point___position_measured = equilibrium__time_slice___constraints__strike_point___position_measured()
    var"time_measurement" :: Union{Missing, Float64} = missing
    var"chi_squared_r" :: Union{Missing, Float64} = missing
    var"position_reconstructed" :: equilibrium__time_slice___constraints__strike_point___position_reconstructed = equilibrium__time_slice___constraints__strike_point___position_reconstructed()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__strike_point(var"chi_squared_z"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"position_measured"=equilibrium__time_slice___constraints__strike_point___position_measured(), var"time_measurement"=missing, var"chi_squared_r"=missing, var"position_reconstructed"=equilibrium__time_slice___constraints__strike_point___position_reconstructed(), _parent=WeakRef(missing))
        obj = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        setfield!(obj.position_measured, :_parent, WeakRef(obj))
        setfield!(obj.position_reconstructed, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q___position <: FDS
    var"phi" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"position" :: equilibrium__time_slice___constraints__q___position = equilibrium__time_slice___constraints__q___position()
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__q(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"position"=equilibrium__time_slice___constraints__q___position(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"position", var"time_measurement", _parent)
        setfield!(obj.position, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pressure <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_passive_current <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_current <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e_line <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__mse_polarisation_angle <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment <: FDS
    var"magnetisation_r" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r()
    var"magnetisation_z" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__iron_core_segment(var"magnetisation_r"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(), var"magnetisation_z"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(), _parent=WeakRef(missing))
        obj = new(var"magnetisation_r", var"magnetisation_z", _parent)
        setfield!(obj.magnetisation_r, :_parent, WeakRef(obj))
        setfield!(obj.magnetisation_z, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__ip <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__flux_loop <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__faraday_angle <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__diamagnetic_flux <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__bpol_probe <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__b_field_tor_vacuum_r <: FDS
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
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
        obj = new(var"faraday_angle", var"n_e", var"ip", var"n_e_line", var"pf_current", var"strike_point", var"x_point", var"iron_core_segment", var"pressure", var"diamagnetic_flux", var"pf_passive_current", var"bpol_probe", var"mse_polarisation_angle", var"q", var"b_field_tor_vacuum_r", var"flux_loop", _parent)
        setfield!(obj.faraday_angle, :_parent, WeakRef(obj))
        setfield!(obj.n_e, :_parent, WeakRef(obj))
        setfield!(obj.ip, :_parent, WeakRef(obj))
        setfield!(obj.n_e_line, :_parent, WeakRef(obj))
        setfield!(obj.pf_current, :_parent, WeakRef(obj))
        setfield!(obj.strike_point, :_parent, WeakRef(obj))
        setfield!(obj.x_point, :_parent, WeakRef(obj))
        setfield!(obj.iron_core_segment, :_parent, WeakRef(obj))
        setfield!(obj.pressure, :_parent, WeakRef(obj))
        setfield!(obj.diamagnetic_flux, :_parent, WeakRef(obj))
        setfield!(obj.pf_passive_current, :_parent, WeakRef(obj))
        setfield!(obj.bpol_probe, :_parent, WeakRef(obj))
        setfield!(obj.mse_polarisation_angle, :_parent, WeakRef(obj))
        setfield!(obj.q, :_parent, WeakRef(obj))
        setfield!(obj.b_field_tor_vacuum_r, :_parent, WeakRef(obj))
        setfield!(obj.flux_loop, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__x_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__strike_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__geometric_axis <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__gap <: FDS
    var"name" :: Union{Missing, String} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"value" :: Union{Missing, Float64} = missing
    var"identifier" :: Union{Missing, String} = missing
    var"angle" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__closest_wall_point <: FDS
    var"distance" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__active_limiter_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix <: FDS
    var"psi" :: Union{Missing, Float64} = missing
    var"elongation_lower" :: Union{Missing, Float64} = missing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_separatrix__x_point[])
    var"gap" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__gap} = FDSvector(equilibrium__time_slice___boundary_separatrix__gap[])
    var"triangularity" :: Union{Missing, Float64} = missing
    var"elongation_upper" :: Union{Missing, Float64} = missing
    var"triangularity_upper" :: Union{Missing, Float64} = missing
    var"outline" :: equilibrium__time_slice___boundary_separatrix__outline = equilibrium__time_slice___boundary_separatrix__outline()
    var"dr_dz_zero_point" :: equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point = equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point()
    var"squareness_lower_outer" :: Union{Missing, Float64} = missing
    var"triangularity_lower" :: Union{Missing, Float64} = missing
    var"minor_radius" :: Union{Missing, Float64} = missing
    var"squareness_upper_inner" :: Union{Missing, Float64} = missing
    var"squareness_upper_outer" :: Union{Missing, Float64} = missing
    var"squareness_lower_inner" :: Union{Missing, Float64} = missing
    var"geometric_axis" :: equilibrium__time_slice___boundary_separatrix__geometric_axis = equilibrium__time_slice___boundary_separatrix__geometric_axis()
    var"elongation" :: Union{Missing, Float64} = missing
    var"active_limiter_point" :: equilibrium__time_slice___boundary_separatrix__active_limiter_point = equilibrium__time_slice___boundary_separatrix__active_limiter_point()
    var"closest_wall_point" :: equilibrium__time_slice___boundary_separatrix__closest_wall_point = equilibrium__time_slice___boundary_separatrix__closest_wall_point()
    var"type" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix(var"psi"=missing, var"elongation_lower"=missing, var"strike_point"=FDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___boundary_separatrix__x_point[]), var"gap"=FDSvector(equilibrium__time_slice___boundary_separatrix__gap[]), var"triangularity"=missing, var"elongation_upper"=missing, var"triangularity_upper"=missing, var"outline"=equilibrium__time_slice___boundary_separatrix__outline(), var"dr_dz_zero_point"=equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(), var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"minor_radius"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"squareness_lower_inner"=missing, var"geometric_axis"=equilibrium__time_slice___boundary_separatrix__geometric_axis(), var"elongation"=missing, var"active_limiter_point"=equilibrium__time_slice___boundary_separatrix__active_limiter_point(), var"closest_wall_point"=equilibrium__time_slice___boundary_separatrix__closest_wall_point(), var"type"=missing, _parent=WeakRef(missing))
        obj = new(var"psi", var"elongation_lower", var"strike_point", var"x_point", var"gap", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"dr_dz_zero_point", var"squareness_lower_outer", var"triangularity_lower", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"closest_wall_point", var"type", _parent)
        setfield!(obj.strike_point, :_parent, WeakRef(obj))
        setfield!(obj.x_point, :_parent, WeakRef(obj))
        setfield!(obj.gap, :_parent, WeakRef(obj))
        setfield!(obj.outline, :_parent, WeakRef(obj))
        setfield!(obj.dr_dz_zero_point, :_parent, WeakRef(obj))
        setfield!(obj.geometric_axis, :_parent, WeakRef(obj))
        setfield!(obj.active_limiter_point, :_parent, WeakRef(obj))
        setfield!(obj.closest_wall_point, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__x_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__strike_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix <: FDS
    var"psi" :: Union{Missing, Float64} = missing
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[])
    var"distance_inner_outer" :: Union{Missing, Float64} = missing
    var"outline" :: equilibrium__time_slice___boundary_secondary_separatrix__outline = equilibrium__time_slice___boundary_secondary_separatrix__outline()
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix(var"psi"=missing, var"x_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[]), var"distance_inner_outer"=missing, var"outline"=equilibrium__time_slice___boundary_secondary_separatrix__outline(), var"strike_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[]), _parent=WeakRef(missing))
        obj = new(var"psi", var"x_point", var"distance_inner_outer", var"outline", var"strike_point", _parent)
        setfield!(obj.x_point, :_parent, WeakRef(obj))
        setfield!(obj.outline, :_parent, WeakRef(obj))
        setfield!(obj.strike_point, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__x_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__strike_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__lcfs <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__geometric_axis <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__active_limiter_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary <: FDS
    var"psi" :: Union{Missing, Float64} = missing
    var"lcfs" :: equilibrium__time_slice___boundary__lcfs = equilibrium__time_slice___boundary__lcfs()
    var"elongation_lower" :: Union{Missing, Float64} = missing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary__strike_point} = FDSvector(equilibrium__time_slice___boundary__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary__x_point} = FDSvector(equilibrium__time_slice___boundary__x_point[])
    var"triangularity" :: Union{Missing, Float64} = missing
    var"elongation_upper" :: Union{Missing, Float64} = missing
    var"triangularity_upper" :: Union{Missing, Float64} = missing
    var"outline" :: equilibrium__time_slice___boundary__outline = equilibrium__time_slice___boundary__outline()
    var"squareness_lower_outer" :: Union{Missing, Float64} = missing
    var"triangularity_lower" :: Union{Missing, Float64} = missing
    var"psi_norm" :: Union{Missing, Float64} = missing
    var"minor_radius" :: Union{Missing, Float64} = missing
    var"squareness_upper_inner" :: Union{Missing, Float64} = missing
    var"squareness_upper_outer" :: Union{Missing, Float64} = missing
    var"squareness_lower_inner" :: Union{Missing, Float64} = missing
    var"geometric_axis" :: equilibrium__time_slice___boundary__geometric_axis = equilibrium__time_slice___boundary__geometric_axis()
    var"elongation" :: Union{Missing, Float64} = missing
    var"active_limiter_point" :: equilibrium__time_slice___boundary__active_limiter_point = equilibrium__time_slice___boundary__active_limiter_point()
    var"b_flux_pol_norm" :: Union{Missing, Float64} = missing
    var"type" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary(var"psi"=missing, var"lcfs"=equilibrium__time_slice___boundary__lcfs(), var"elongation_lower"=missing, var"strike_point"=FDSvector(equilibrium__time_slice___boundary__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___boundary__x_point[]), var"triangularity"=missing, var"elongation_upper"=missing, var"triangularity_upper"=missing, var"outline"=equilibrium__time_slice___boundary__outline(), var"squareness_lower_outer"=missing, var"triangularity_lower"=missing, var"psi_norm"=missing, var"minor_radius"=missing, var"squareness_upper_inner"=missing, var"squareness_upper_outer"=missing, var"squareness_lower_inner"=missing, var"geometric_axis"=equilibrium__time_slice___boundary__geometric_axis(), var"elongation"=missing, var"active_limiter_point"=equilibrium__time_slice___boundary__active_limiter_point(), var"b_flux_pol_norm"=missing, var"type"=missing, _parent=WeakRef(missing))
        obj = new(var"psi", var"lcfs", var"elongation_lower", var"strike_point", var"x_point", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"squareness_lower_outer", var"triangularity_lower", var"psi_norm", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"b_flux_pol_norm", var"type", _parent)
        setfield!(obj.lcfs, :_parent, WeakRef(obj))
        setfield!(obj.strike_point, :_parent, WeakRef(obj))
        setfield!(obj.x_point, :_parent, WeakRef(obj))
        setfield!(obj.outline, :_parent, WeakRef(obj))
        setfield!(obj.geometric_axis, :_parent, WeakRef(obj))
        setfield!(obj.active_limiter_point, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice <: FDS
    var"time" :: Union{Missing, Float64} = missing
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
        obj = new(var"time", var"ggd", var"profiles_1d", var"boundary", var"constraints", var"global_quantities", var"convergence", var"coordinate_system", var"boundary_secondary_separatrix", var"boundary_separatrix", var"profiles_2d", _parent)
        setfield!(obj.ggd, :_parent, WeakRef(obj))
        setfield!(obj.profiles_1d, :_parent, WeakRef(obj))
        setfield!(obj.boundary, :_parent, WeakRef(obj))
        setfield!(obj.constraints, :_parent, WeakRef(obj))
        setfield!(obj.global_quantities, :_parent, WeakRef(obj))
        setfield!(obj.convergence, :_parent, WeakRef(obj))
        setfield!(obj.coordinate_system, :_parent, WeakRef(obj))
        setfield!(obj.boundary_secondary_separatrix, :_parent, WeakRef(obj))
        setfield!(obj.boundary_separatrix, :_parent, WeakRef(obj))
        setfield!(obj.profiles_2d, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__ids_properties <: FDS
    var"provider" :: Union{Missing, String} = missing
    var"version_put" :: equilibrium__ids_properties__version_put = equilibrium__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Int64} = missing
    var"source" :: Union{Missing, String} = missing
    var"creation_date" :: Union{Missing, String} = missing
    var"comment" :: Union{Missing, String} = missing
    var"occurrence" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__ids_properties(var"provider"=missing, var"version_put"=equilibrium__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        setfield!(obj.version_put, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary <: FDS
    var"neighbours" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object <: FDS
    var"nodes" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"measure" :: Union{Missing, Float64} = missing
    var"geometry" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        setfield!(obj.boundary, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension(var"object"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        obj = new(var"object", _parent)
        setfield!(obj.object, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___geometry_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space <: FDS
    var"coordinates_type" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"geometry_type" :: equilibrium__grids_ggd___grid___space___geometry_type = equilibrium__grids_ggd___grid___space___geometry_type()
    var"identifier" :: equilibrium__grids_ggd___grid___space___identifier = equilibrium__grids_ggd___grid___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__grids_ggd___grid___space___geometry_type(), var"identifier"=equilibrium__grids_ggd___grid___space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[]), _parent=WeakRef(missing))
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        setfield!(obj.geometry_type, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        setfield!(obj.objects_per_dimension, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element___object <: FDS
    var"dimension" :: Union{Missing, Int64} = missing
    var"space" :: Union{Missing, Int64} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element___object} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___element(var"object"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[]), _parent=WeakRef(missing))
        obj = new(var"object", _parent)
        setfield!(obj.object, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___base <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___base} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[])
    var"metric" :: equilibrium__grids_ggd___grid___grid_subset___metric = equilibrium__grids_ggd___grid___grid_subset___metric()
    var"dimension" :: Union{Missing, Int64} = missing
    var"identifier" :: equilibrium__grids_ggd___grid___grid_subset___identifier = equilibrium__grids_ggd___grid___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset(var"base"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[]), var"metric"=equilibrium__grids_ggd___grid___grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__grids_ggd___grid___grid_subset___identifier(), var"element"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[]), _parent=WeakRef(missing))
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        setfield!(obj.base, :_parent, WeakRef(obj))
        setfield!(obj.metric, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset} = FDSvector(equilibrium__grids_ggd___grid___grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space} = FDSvector(equilibrium__grids_ggd___grid___space[])
    var"identifier" :: equilibrium__grids_ggd___grid___identifier = equilibrium__grids_ggd___grid___identifier()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid(var"grid_subset"=FDSvector(equilibrium__grids_ggd___grid___grid_subset[]), var"space"=FDSvector(equilibrium__grids_ggd___grid___space[]), var"identifier"=equilibrium__grids_ggd___grid___identifier(), _parent=WeakRef(missing))
        obj = new(var"grid_subset", var"space", var"identifier", _parent)
        setfield!(obj.grid_subset, :_parent, WeakRef(obj))
        setfield!(obj.space, :_parent, WeakRef(obj))
        setfield!(obj.identifier, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd <: FDS
    var"time" :: Union{Missing, Float64} = missing
    var"grid" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid} = FDSvector(equilibrium__grids_ggd___grid[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd(var"time"=missing, var"grid"=FDSvector(equilibrium__grids_ggd___grid[]), _parent=WeakRef(missing))
        obj = new(var"time", var"grid", _parent)
        setfield!(obj.grid, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__code__library <: FDS
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct equilibrium__code <: FDS
    var"library" :: FDSvector{T} where {T<:equilibrium__code__library} = FDSvector(equilibrium__code__library[])
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"output_flag" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__code(var"library"=FDSvector(equilibrium__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        setfield!(obj.library, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct equilibrium <: FDS
    var"time_slice" :: FDSvector{T} where {T<:equilibrium__time_slice} = FDSvector(equilibrium__time_slice[])
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ids_properties" :: equilibrium__ids_properties = equilibrium__ids_properties()
    var"grids_ggd" :: FDSvector{T} where {T<:equilibrium__grids_ggd} = FDSvector(equilibrium__grids_ggd[])
    var"vacuum_toroidal_field" :: equilibrium__vacuum_toroidal_field = equilibrium__vacuum_toroidal_field()
    var"code" :: equilibrium__code = equilibrium__code()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium(var"time_slice"=FDSvector(equilibrium__time_slice[]), var"time"=missing, var"ids_properties"=equilibrium__ids_properties(), var"grids_ggd"=FDSvector(equilibrium__grids_ggd[]), var"vacuum_toroidal_field"=equilibrium__vacuum_toroidal_field(), var"code"=equilibrium__code(), _parent=WeakRef(missing))
        obj = new(var"time_slice", var"time", var"ids_properties", var"grids_ggd", var"vacuum_toroidal_field", var"code", _parent)
        setfield!(obj.time_slice, :_parent, WeakRef(obj))
        setfield!(obj.ids_properties, :_parent, WeakRef(obj))
        setfield!(obj.grids_ggd, :_parent, WeakRef(obj))
        setfield!(obj.vacuum_toroidal_field, :_parent, WeakRef(obj))
        setfield!(obj.code, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct dataset_description__simulation <: FDS
    var"time_ended" :: Union{Missing, String} = missing
    var"time_begun" :: Union{Missing, String} = missing
    var"time_current" :: Union{Missing, Float64} = missing
    var"time_restart" :: Union{Missing, Float64} = missing
    var"workflow" :: Union{Missing, String} = missing
    var"comment_after" :: Union{Missing, String} = missing
    var"time_begin" :: Union{Missing, Float64} = missing
    var"time_end" :: Union{Missing, Float64} = missing
    var"comment_before" :: Union{Missing, String} = missing
    var"time_step" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct dataset_description__pulse_time_end_epoch <: FDS
    var"nanoseconds" :: Union{Missing, Int64} = missing
    var"seconds" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct dataset_description__pulse_time_begin_epoch <: FDS
    var"nanoseconds" :: Union{Missing, Int64} = missing
    var"seconds" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct dataset_description__parent_entry <: FDS
    var"pulse_type" :: Union{Missing, String} = missing
    var"run" :: Union{Missing, Int64} = missing
    var"machine" :: Union{Missing, String} = missing
    var"pulse" :: Union{Missing, Int64} = missing
    var"user" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct dataset_description__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct dataset_description__ids_properties <: FDS
    var"provider" :: Union{Missing, String} = missing
    var"version_put" :: dataset_description__ids_properties__version_put = dataset_description__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Int64} = missing
    var"source" :: Union{Missing, String} = missing
    var"creation_date" :: Union{Missing, String} = missing
    var"comment" :: Union{Missing, String} = missing
    var"occurrence" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__ids_properties(var"provider"=missing, var"version_put"=dataset_description__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        setfield!(obj.version_put, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct dataset_description__data_entry <: FDS
    var"pulse_type" :: Union{Missing, String} = missing
    var"run" :: Union{Missing, Int64} = missing
    var"machine" :: Union{Missing, String} = missing
    var"pulse" :: Union{Missing, Int64} = missing
    var"user" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct dataset_description <: FDS
    var"pulse_time_begin_epoch" :: dataset_description__pulse_time_begin_epoch = dataset_description__pulse_time_begin_epoch()
    var"imas_version" :: Union{Missing, String} = missing
    var"ids_properties" :: dataset_description__ids_properties = dataset_description__ids_properties()
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dd_version" :: Union{Missing, String} = missing
    var"parent_entry" :: dataset_description__parent_entry = dataset_description__parent_entry()
    var"simulation" :: dataset_description__simulation = dataset_description__simulation()
    var"pulse_time_end_epoch" :: dataset_description__pulse_time_end_epoch = dataset_description__pulse_time_end_epoch()
    var"data_entry" :: dataset_description__data_entry = dataset_description__data_entry()
    var"pulse_time_begin" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description(var"pulse_time_begin_epoch"=dataset_description__pulse_time_begin_epoch(), var"imas_version"=missing, var"ids_properties"=dataset_description__ids_properties(), var"time"=missing, var"dd_version"=missing, var"parent_entry"=dataset_description__parent_entry(), var"simulation"=dataset_description__simulation(), var"pulse_time_end_epoch"=dataset_description__pulse_time_end_epoch(), var"data_entry"=dataset_description__data_entry(), var"pulse_time_begin"=missing, _parent=WeakRef(missing))
        obj = new(var"pulse_time_begin_epoch", var"imas_version", var"ids_properties", var"time", var"dd_version", var"parent_entry", var"simulation", var"pulse_time_end_epoch", var"data_entry", var"pulse_time_begin", _parent)
        setfield!(obj.pulse_time_begin_epoch, :_parent, WeakRef(obj))
        setfield!(obj.ids_properties, :_parent, WeakRef(obj))
        setfield!(obj.parent_entry, :_parent, WeakRef(obj))
        setfield!(obj.simulation, :_parent, WeakRef(obj))
        setfield!(obj.pulse_time_end_epoch, :_parent, WeakRef(obj))
        setfield!(obj.data_entry, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"r0" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method = core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___zeff_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method = core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___t_i_average_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___neutral_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state <: FDS
    var"label" :: Union{Missing, String} = missing
    var"vibrational_level" :: Union{Missing, Float64} = missing
    var"temperature" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"electron_configuration" :: Union{Missing, String} = missing
    var"pressure" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"vibrational_mode" :: Union{Missing, String} = missing
    var"pressure_fast_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"velocity" :: core_profiles__profiles_1d___neutral___state___velocity = core_profiles__profiles_1d___neutral___state___velocity()
    var"density" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fast" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"neutral_type" :: core_profiles__profiles_1d___neutral___state___neutral_type = core_profiles__profiles_1d___neutral___state___neutral_type()
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___state(var"label"=missing, var"vibrational_level"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"electron_configuration"=missing, var"pressure"=missing, var"density_thermal"=missing, var"vibrational_mode"=missing, var"pressure_fast_parallel"=missing, var"velocity"=core_profiles__profiles_1d___neutral___state___velocity(), var"density"=missing, var"density_fast"=missing, var"neutral_type"=core_profiles__profiles_1d___neutral___state___neutral_type(), _parent=WeakRef(missing))
        obj = new(var"label", var"vibrational_level", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"velocity", var"density", var"density_fast", var"neutral_type", _parent)
        setfield!(obj.velocity, :_parent, WeakRef(obj))
        setfield!(obj.neutral_type, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___element <: FDS
    var"atoms_n" :: Union{Missing, Int64} = missing
    var"z_n" :: Union{Missing, Float64} = missing
    var"multiplicity" :: Union{Missing, Float64} = missing
    var"a" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral <: FDS
    var"label" :: Union{Missing, String} = missing
    var"temperature" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ion_index" :: Union{Missing, Int64} = missing
    var"multiple_states_flag" :: Union{Missing, Int64} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_fast_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral___state} = FDSvector(core_profiles__profiles_1d___neutral___state[])
    var"velocity" :: core_profiles__profiles_1d___neutral___velocity = core_profiles__profiles_1d___neutral___velocity()
    var"density" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fast" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral___element} = FDSvector(core_profiles__profiles_1d___neutral___element[])
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral(var"label"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"ion_index"=missing, var"multiple_states_flag"=missing, var"pressure_fast_perpendicular"=missing, var"pressure"=missing, var"density_thermal"=missing, var"pressure_fast_parallel"=missing, var"state"=FDSvector(core_profiles__profiles_1d___neutral___state[]), var"velocity"=core_profiles__profiles_1d___neutral___velocity(), var"density"=missing, var"density_fast"=missing, var"element"=FDSvector(core_profiles__profiles_1d___neutral___element[]), _parent=WeakRef(missing))
        obj = new(var"label", var"temperature", var"pressure_thermal", var"ion_index", var"multiple_states_flag", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"pressure_fast_parallel", var"state", var"velocity", var"density", var"density_fast", var"element", _parent)
        setfield!(obj.state, :_parent, WeakRef(obj))
        setfield!(obj.velocity, :_parent, WeakRef(obj))
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___temperature_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state___density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state <: FDS
    var"label" :: Union{Missing, String} = missing
    var"vibrational_level" :: Union{Missing, Float64} = missing
    var"rotation_frequency_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"temperature" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z_min" :: Union{Missing, Float64} = missing
    var"electron_configuration" :: Union{Missing, String} = missing
    var"pressure" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"vibrational_mode" :: Union{Missing, String} = missing
    var"pressure_fast_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z_average_square_1d" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"velocity" :: core_profiles__profiles_1d___ion___state___velocity = core_profiles__profiles_1d___ion___state___velocity()
    var"z_average" :: Union{Missing, Float64} = missing
    var"z_max" :: Union{Missing, Float64} = missing
    var"z_square_average" :: Union{Missing, Float64} = missing
    var"ionisation_potential" :: Union{Missing, Float64} = missing
    var"density" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fast" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z_average_1d" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fit" :: core_profiles__profiles_1d___ion___state___density_fit = core_profiles__profiles_1d___ion___state___density_fit()
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state(var"label"=missing, var"vibrational_level"=missing, var"rotation_frequency_tor"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"z_min"=missing, var"electron_configuration"=missing, var"pressure"=missing, var"density_thermal"=missing, var"vibrational_mode"=missing, var"pressure_fast_parallel"=missing, var"z_average_square_1d"=missing, var"velocity"=core_profiles__profiles_1d___ion___state___velocity(), var"z_average"=missing, var"z_max"=missing, var"z_square_average"=missing, var"ionisation_potential"=missing, var"density"=missing, var"density_fast"=missing, var"z_average_1d"=missing, var"density_fit"=core_profiles__profiles_1d___ion___state___density_fit(), _parent=WeakRef(missing))
        obj = new(var"label", var"vibrational_level", var"rotation_frequency_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"z_min", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"z_average_square_1d", var"velocity", var"z_average", var"z_max", var"z_square_average", var"ionisation_potential", var"density", var"density_fast", var"z_average_1d", var"density_fit", _parent)
        setfield!(obj.velocity, :_parent, WeakRef(obj))
        setfield!(obj.density_fit, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___element <: FDS
    var"atoms_n" :: Union{Missing, Int64} = missing
    var"z_n" :: Union{Missing, Float64} = missing
    var"multiplicity" :: Union{Missing, Float64} = missing
    var"a" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion <: FDS
    var"label" :: Union{Missing, String} = missing
    var"rotation_frequency_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"temperature_validity" :: Union{Missing, Int64} = missing
    var"velocity_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"temperature" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z_ion_1d" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"multiple_states_flag" :: Union{Missing, Int64} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"neutral_index" :: Union{Missing, Int64} = missing
    var"pressure" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_validity" :: Union{Missing, Int64} = missing
    var"pressure_fast_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion___state} = FDSvector(core_profiles__profiles_1d___ion___state[])
    var"velocity" :: core_profiles__profiles_1d___ion___velocity = core_profiles__profiles_1d___ion___velocity()
    var"z_ion" :: Union{Missing, Float64} = missing
    var"temperature_fit" :: core_profiles__profiles_1d___ion___temperature_fit = core_profiles__profiles_1d___ion___temperature_fit()
    var"density" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"velocity_pol" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fast" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fit" :: core_profiles__profiles_1d___ion___density_fit = core_profiles__profiles_1d___ion___density_fit()
    var"z_ion_square_1d" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion___element} = FDSvector(core_profiles__profiles_1d___ion___element[])
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion(var"label"=missing, var"rotation_frequency_tor"=missing, var"temperature_validity"=missing, var"velocity_tor"=missing, var"temperature"=missing, var"z_ion_1d"=missing, var"pressure_thermal"=missing, var"multiple_states_flag"=missing, var"pressure_fast_perpendicular"=missing, var"neutral_index"=missing, var"pressure"=missing, var"density_thermal"=missing, var"density_validity"=missing, var"pressure_fast_parallel"=missing, var"state"=FDSvector(core_profiles__profiles_1d___ion___state[]), var"velocity"=core_profiles__profiles_1d___ion___velocity(), var"z_ion"=missing, var"temperature_fit"=core_profiles__profiles_1d___ion___temperature_fit(), var"density"=missing, var"velocity_pol"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___ion___density_fit(), var"z_ion_square_1d"=missing, var"element"=FDSvector(core_profiles__profiles_1d___ion___element[]), _parent=WeakRef(missing))
        obj = new(var"label", var"rotation_frequency_tor", var"temperature_validity", var"velocity_tor", var"temperature", var"z_ion_1d", var"pressure_thermal", var"multiple_states_flag", var"pressure_fast_perpendicular", var"neutral_index", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"state", var"velocity", var"z_ion", var"temperature_fit", var"density", var"velocity_pol", var"density_fast", var"density_fit", var"z_ion_square_1d", var"element", _parent)
        setfield!(obj.state, :_parent, WeakRef(obj))
        setfield!(obj.velocity, :_parent, WeakRef(obj))
        setfield!(obj.temperature_fit, :_parent, WeakRef(obj))
        setfield!(obj.density_fit, :_parent, WeakRef(obj))
        setfield!(obj.element, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___grid <: FDS
    var"psi" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"psi_boundary" :: Union{Missing, Float64} = missing
    var"volume" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"area" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_pol_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"surface" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"psi_magnetic_axis" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__temperature_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit <: FDS
    var"local" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"chi_squared" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"reconstructed" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_width" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rho_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"weight" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"source" :: Union{Missing, AbstractFDVector{String}} = missing
    var"measured" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method = core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__density_fit(var"local"=missing, var"chi_squared"=missing, var"parameters"=missing, var"reconstructed"=missing, var"time_measurement_width"=missing, var"rho_tor_norm"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(), var"time_measurement"=missing, _parent=WeakRef(missing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        setfield!(obj.time_measurement_slice_method, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons <: FDS
    var"temperature_validity" :: Union{Missing, Int64} = missing
    var"velocity_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"temperature" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_fast_perpendicular" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_validity" :: Union{Missing, Int64} = missing
    var"pressure_fast_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"velocity" :: core_profiles__profiles_1d___electrons__velocity = core_profiles__profiles_1d___electrons__velocity()
    var"temperature_fit" :: core_profiles__profiles_1d___electrons__temperature_fit = core_profiles__profiles_1d___electrons__temperature_fit()
    var"density" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"velocity_pol" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"collisionality_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fast" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"density_fit" :: core_profiles__profiles_1d___electrons__density_fit = core_profiles__profiles_1d___electrons__density_fit()
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons(var"temperature_validity"=missing, var"velocity_tor"=missing, var"temperature"=missing, var"pressure_thermal"=missing, var"pressure_fast_perpendicular"=missing, var"pressure"=missing, var"density_thermal"=missing, var"density_validity"=missing, var"pressure_fast_parallel"=missing, var"velocity"=core_profiles__profiles_1d___electrons__velocity(), var"temperature_fit"=core_profiles__profiles_1d___electrons__temperature_fit(), var"density"=missing, var"velocity_pol"=missing, var"collisionality_norm"=missing, var"density_fast"=missing, var"density_fit"=core_profiles__profiles_1d___electrons__density_fit(), _parent=WeakRef(missing))
        obj = new(var"temperature_validity", var"velocity_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"velocity", var"temperature_fit", var"density", var"velocity_pol", var"collisionality_norm", var"density_fast", var"density_fit", _parent)
        setfield!(obj.velocity, :_parent, WeakRef(obj))
        setfield!(obj.temperature_fit, :_parent, WeakRef(obj))
        setfield!(obj.density_fit, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___e_field <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__profiles_1d <: FDS
    var"pressure_ion_total" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"time" :: Union{Missing, Float64} = missing
    var"t_i_average_fit" :: core_profiles__profiles_1d___t_i_average_fit = core_profiles__profiles_1d___t_i_average_fit()
    var"neutral" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral} = FDSvector(core_profiles__profiles_1d___neutral[])
    var"n_i_thermal_total" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"magnetic_shear" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ion" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion} = FDSvector(core_profiles__profiles_1d___ion[])
    var"j_total" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"rotation_frequency_tor_sonic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"pressure_thermal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"j_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"current_parallel_inside" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"j_non_inductive" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"e_field_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"momentum_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"conductivity_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"electrons" :: core_profiles__profiles_1d___electrons = core_profiles__profiles_1d___electrons()
    var"pressure_perpendicular" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"q" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"t_i_average" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"j_ohmic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid" :: core_profiles__profiles_1d___grid = core_profiles__profiles_1d___grid()
    var"phi_potential" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"j_bootstrap" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"zeff_fit" :: core_profiles__profiles_1d___zeff_fit = core_profiles__profiles_1d___zeff_fit()
    var"pressure_parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"e_field" :: core_profiles__profiles_1d___e_field = core_profiles__profiles_1d___e_field()
    var"zeff" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"n_i_total_over_n_e" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d(var"pressure_ion_total"=missing, var"time"=missing, var"t_i_average_fit"=core_profiles__profiles_1d___t_i_average_fit(), var"neutral"=FDSvector(core_profiles__profiles_1d___neutral[]), var"n_i_thermal_total"=missing, var"magnetic_shear"=missing, var"ion"=FDSvector(core_profiles__profiles_1d___ion[]), var"j_total"=missing, var"rotation_frequency_tor_sonic"=missing, var"pressure_thermal"=missing, var"j_tor"=missing, var"current_parallel_inside"=missing, var"j_non_inductive"=missing, var"e_field_parallel"=missing, var"momentum_tor"=missing, var"conductivity_parallel"=missing, var"electrons"=core_profiles__profiles_1d___electrons(), var"pressure_perpendicular"=missing, var"q"=missing, var"t_i_average"=missing, var"j_ohmic"=missing, var"grid"=core_profiles__profiles_1d___grid(), var"phi_potential"=missing, var"j_bootstrap"=missing, var"zeff_fit"=core_profiles__profiles_1d___zeff_fit(), var"pressure_parallel"=missing, var"e_field"=core_profiles__profiles_1d___e_field(), var"zeff"=missing, var"n_i_total_over_n_e"=missing, _parent=WeakRef(missing))
        obj = new(var"pressure_ion_total", var"time", var"t_i_average_fit", var"neutral", var"n_i_thermal_total", var"magnetic_shear", var"ion", var"j_total", var"rotation_frequency_tor_sonic", var"pressure_thermal", var"j_tor", var"current_parallel_inside", var"j_non_inductive", var"e_field_parallel", var"momentum_tor", var"conductivity_parallel", var"electrons", var"pressure_perpendicular", var"q", var"t_i_average", var"j_ohmic", var"grid", var"phi_potential", var"j_bootstrap", var"zeff_fit", var"pressure_parallel", var"e_field", var"zeff", var"n_i_total_over_n_e", _parent)
        setfield!(obj.t_i_average_fit, :_parent, WeakRef(obj))
        setfield!(obj.neutral, :_parent, WeakRef(obj))
        setfield!(obj.ion, :_parent, WeakRef(obj))
        setfield!(obj.electrons, :_parent, WeakRef(obj))
        setfield!(obj.grid, :_parent, WeakRef(obj))
        setfield!(obj.zeff_fit, :_parent, WeakRef(obj))
        setfield!(obj.e_field, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__ids_properties <: FDS
    var"provider" :: Union{Missing, String} = missing
    var"version_put" :: core_profiles__ids_properties__version_put = core_profiles__ids_properties__version_put()
    var"homogeneous_time" :: Union{Missing, Int64} = missing
    var"source" :: Union{Missing, String} = missing
    var"creation_date" :: Union{Missing, String} = missing
    var"comment" :: Union{Missing, String} = missing
    var"occurrence" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__ids_properties(var"provider"=missing, var"version_put"=core_profiles__ids_properties__version_put(), var"homogeneous_time"=missing, var"source"=missing, var"creation_date"=missing, var"comment"=missing, var"occurrence"=missing, _parent=WeakRef(missing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        setfield!(obj.version_put, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__global_quantities <: FDS
    var"beta_tor_norm" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"resistive_psi_losses" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ip" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"li_3" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"t_i_average_peaking" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"t_e_peaking" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"beta_tor" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z_eff_resistive" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ejima" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"energy_diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"li" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"current_non_inductive" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"v_loop" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"beta_pol" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"current_bootstrap" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__code__library <: FDS
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
end

Base.@kwdef mutable struct core_profiles__code <: FDS
    var"library" :: FDSvector{T} where {T<:core_profiles__code__library} = FDSvector(core_profiles__code__library[])
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"output_flag" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__code(var"library"=FDSvector(core_profiles__code__library[]), var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"output_flag"=missing, var"version"=missing, _parent=WeakRef(missing))
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        setfield!(obj.library, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct core_profiles <: FDS
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"ids_properties" :: core_profiles__ids_properties = core_profiles__ids_properties()
    var"vacuum_toroidal_field" :: core_profiles__vacuum_toroidal_field = core_profiles__vacuum_toroidal_field()
    var"code" :: core_profiles__code = core_profiles__code()
    var"global_quantities" :: core_profiles__global_quantities = core_profiles__global_quantities()
    var"profiles_1d" :: FDSvector{T} where {T<:core_profiles__profiles_1d} = FDSvector(core_profiles__profiles_1d[])
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles(var"time"=missing, var"ids_properties"=core_profiles__ids_properties(), var"vacuum_toroidal_field"=core_profiles__vacuum_toroidal_field(), var"code"=core_profiles__code(), var"global_quantities"=core_profiles__global_quantities(), var"profiles_1d"=FDSvector(core_profiles__profiles_1d[]), _parent=WeakRef(missing))
        obj = new(var"time", var"ids_properties", var"vacuum_toroidal_field", var"code", var"global_quantities", var"profiles_1d", _parent)
        setfield!(obj.ids_properties, :_parent, WeakRef(obj))
        setfield!(obj.vacuum_toroidal_field, :_parent, WeakRef(obj))
        setfield!(obj.code, :_parent, WeakRef(obj))
        setfield!(obj.global_quantities, :_parent, WeakRef(obj))
        setfield!(obj.profiles_1d, :_parent, WeakRef(obj))
        return obj
    end
end

Base.@kwdef mutable struct dd <: FDS
    var"equilibrium" :: Union{Missing, equilibrium} = equilibrium()
    var"core_profiles" :: Union{Missing, core_profiles} = core_profiles()
    var"wall" :: Union{Missing, wall} = wall()
    var"dataset_description" :: Union{Missing, dataset_description} = dataset_description()
    _parent :: WeakRef = WeakRef(missing)
    function dd(var"equilibrium"=equilibrium(), var"core_profiles"=core_profiles(), var"wall"=wall(), var"dataset_description"=dataset_description(), _parent=WeakRef(missing))
        obj = new(var"equilibrium", var"core_profiles", var"wall", var"dataset_description", _parent)
        setfield!(obj.equilibrium, :_parent, WeakRef(obj))
        setfield!(obj.core_profiles, :_parent, WeakRef(obj))
        setfield!(obj.wall, :_parent, WeakRef(obj))
        setfield!(obj.dataset_description, :_parent, WeakRef(obj))
        return obj
    end
end

