abstract type FDS end

include("functionarrays.jl")

convertsion_types = Union{Float64, Int64, String, Array{Float64, N} where N, Array{Int64, N} where N, Array{String, N} where N}

Base.@kwdef mutable struct wall__temperature_reference <: FDS
    var"data" :: Union{Nothing, Float64} = nothing
    var"description" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__temperature_reference(var"data"=nothing, var"description"=nothing, _parent=WeakRef(nothing))
        obj = new(var"data", var"description", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Nothing, String} = nothing
    var"data_dictionary" :: Union{Nothing, String} = nothing
    var"access_layer" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__ids_properties__version_put(var"access_layer_language"=nothing, var"data_dictionary"=nothing, var"access_layer"=nothing, _parent=WeakRef(nothing))
        obj = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__ids_properties <: FDS
    var"provider" :: Union{Nothing, String} = nothing
    var"version_put" :: wall__ids_properties__version_put = wall__ids_properties__version_put()
    var"homogeneous_time" :: Union{Nothing, Int} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"creation_date" :: Union{Nothing, String} = nothing
    var"comment" :: Union{Nothing, String} = nothing
    var"occurrence" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__ids_properties(var"provider"=nothing, var"version_put"=wall__ids_properties__version_put(), var"homogeneous_time"=nothing, var"source"=nothing, var"creation_date"=nothing, var"comment"=nothing, var"occurrence"=nothing, _parent=WeakRef(nothing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        obj.version_put._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral___element <: FDS
    var"atoms_n" :: Union{Nothing, Int} = nothing
    var"z_n" :: Union{Nothing, Float64} = nothing
    var"multiplicity" :: Union{Nothing, Float64} = nothing
    var"a" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__global_quantities__neutral___element(var"atoms_n"=nothing, var"z_n"=nothing, var"multiplicity"=nothing, var"a"=nothing, _parent=WeakRef(nothing))
        obj = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"sputtering_chemical_coefficient" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"gas_puff" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"recycling_particles_coefficient" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"pumping_speed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"particle_flux_from_wall" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"recycling_energy_coefficient" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"wall_inventory" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"particle_flux_from_plasma" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"sputtering_physical_coefficient" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"element" :: FDSvector{T} where {T<:wall__global_quantities__neutral___element} = FDSvector(wall__global_quantities__neutral___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__global_quantities__neutral(var"label"=nothing, var"sputtering_chemical_coefficient"=nothing, var"gas_puff"=nothing, var"recycling_particles_coefficient"=nothing, var"pumping_speed"=nothing, var"particle_flux_from_wall"=nothing, var"recycling_energy_coefficient"=nothing, var"wall_inventory"=nothing, var"particle_flux_from_plasma"=nothing, var"sputtering_physical_coefficient"=nothing, var"element"=FDSvector(wall__global_quantities__neutral___element[]), _parent=WeakRef(nothing))
        obj = new(var"label", var"sputtering_chemical_coefficient", var"gas_puff", var"recycling_particles_coefficient", var"pumping_speed", var"particle_flux_from_wall", var"recycling_energy_coefficient", var"wall_inventory", var"particle_flux_from_plasma", var"sputtering_physical_coefficient", var"element", _parent)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__global_quantities__electrons <: FDS
    var"particle_flux_from_plasma" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gas_puff" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_outer_target" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pumping_speed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"particle_flux_from_wall" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"power_inner_target" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__global_quantities__electrons(var"particle_flux_from_plasma"=nothing, var"gas_puff"=nothing, var"power_outer_target"=nothing, var"pumping_speed"=nothing, var"particle_flux_from_wall"=nothing, var"power_inner_target"=nothing, _parent=WeakRef(nothing))
        obj = new(var"particle_flux_from_plasma", var"gas_puff", var"power_outer_target", var"pumping_speed", var"particle_flux_from_wall", var"power_inner_target", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__global_quantities <: FDS
    var"neutral" :: FDSvector{T} where {T<:wall__global_quantities__neutral} = FDSvector(wall__global_quantities__neutral[])
    var"power_incident" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_radiated" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_inner_target_ion_total" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_conducted" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_convected" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"current_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"electrons" :: wall__global_quantities__electrons = wall__global_quantities__electrons()
    var"power_density_inner_target_max" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_black_body" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_recombination_neutrals" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_to_cooling" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_density_outer_target_max" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_recombination_plasma" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_currents" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"power_neutrals" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__global_quantities(var"neutral"=FDSvector(wall__global_quantities__neutral[]), var"power_incident"=nothing, var"power_radiated"=nothing, var"power_inner_target_ion_total"=nothing, var"temperature"=nothing, var"power_conducted"=nothing, var"power_convected"=nothing, var"current_tor"=nothing, var"electrons"=wall__global_quantities__electrons(), var"power_density_inner_target_max"=nothing, var"power_black_body"=nothing, var"power_recombination_neutrals"=nothing, var"power_to_cooling"=nothing, var"power_density_outer_target_max"=nothing, var"power_recombination_plasma"=nothing, var"power_currents"=nothing, var"power_neutrals"=nothing, _parent=WeakRef(nothing))
        obj = new(var"neutral", var"power_incident", var"power_radiated", var"power_inner_target_ion_total", var"temperature", var"power_conducted", var"power_convected", var"current_tor", var"electrons", var"power_density_inner_target_max", var"power_black_body", var"power_recombination_neutrals", var"power_to_cooling", var"power_density_outer_target_max", var"power_recombination_plasma", var"power_currents", var"power_neutrals", _parent)
        obj.neutral._parent = WeakRef(obj)
        obj.electrons._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__first_wall_power_flux_peak <: FDS
    var"time" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"data" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__first_wall_power_flux_peak(var"time"=nothing, var"data"=nothing, _parent=WeakRef(nothing))
        obj = new(var"time", var"data", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary <: FDS
    var"neighbours" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary(var"neighbours"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"neighbours", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object <: FDS
    var"nodes" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"measure" :: Union{Nothing, Float64} = nothing
    var"geometry" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"boundary" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object(var"nodes"=nothing, var"measure"=nothing, var"geometry"=nothing, var"boundary"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(nothing))
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        obj.boundary._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension(var"object"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[]), _parent=WeakRef(nothing))
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___space___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___geometry_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___space___geometry_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space <: FDS
    var"coordinates_type" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"geometry_type" :: wall__description_ggd___grid_ggd___space___geometry_type = wall__description_ggd___grid_ggd___space___geometry_type()
    var"identifier" :: wall__description_ggd___grid_ggd___space___identifier = wall__description_ggd___grid_ggd___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___space(var"coordinates_type"=nothing, var"geometry_type"=wall__description_ggd___grid_ggd___space___geometry_type(), var"identifier"=wall__description_ggd___grid_ggd___space___identifier(), var"objects_per_dimension"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[]), _parent=WeakRef(nothing))
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        obj.geometry_type._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.objects_per_dimension._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___metric <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___grid_subset___metric(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___grid_subset___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element___object <: FDS
    var"dimension" :: Union{Nothing, Int} = nothing
    var"space" :: Union{Nothing, Int} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___grid_subset___element___object(var"dimension"=nothing, var"space"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"dimension", var"space", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element <: FDS
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element___object} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___grid_subset___element(var"object"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[]), _parent=WeakRef(nothing))
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___base <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___grid_subset___base(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___base} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[])
    var"metric" :: wall__description_ggd___grid_ggd___grid_subset___metric = wall__description_ggd___grid_ggd___grid_subset___metric()
    var"dimension" :: Union{Nothing, Int} = nothing
    var"identifier" :: wall__description_ggd___grid_ggd___grid_subset___identifier = wall__description_ggd___grid_ggd___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd___grid_subset(var"base"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[]), var"metric"=wall__description_ggd___grid_ggd___grid_subset___metric(), var"dimension"=nothing, var"identifier"=wall__description_ggd___grid_ggd___grid_subset___identifier(), var"element"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[]), _parent=WeakRef(nothing))
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        obj.base._parent = WeakRef(obj)
        obj.metric._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd <: FDS
    var"time" :: Union{Nothing, Float64} = nothing
    var"grid_subset" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset} = FDSvector(wall__description_ggd___grid_ggd___grid_subset[])
    var"space" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space} = FDSvector(wall__description_ggd___grid_ggd___space[])
    var"identifier" :: wall__description_ggd___grid_ggd___identifier = wall__description_ggd___grid_ggd___identifier()
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___grid_ggd(var"time"=nothing, var"grid_subset"=FDSvector(wall__description_ggd___grid_ggd___grid_subset[]), var"space"=FDSvector(wall__description_ggd___grid_ggd___space[]), var"identifier"=wall__description_ggd___grid_ggd___identifier(), _parent=WeakRef(nothing))
        obj = new(var"time", var"grid_subset", var"space", var"identifier", _parent)
        obj.grid_subset._parent = WeakRef(obj)
        obj.space._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___temperature <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___ggd___temperature(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___power_density <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___ggd___power_density(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd <: FDS
    var"temperature" :: FDSvector{T} where {T<:wall__description_ggd___ggd___temperature} = FDSvector(wall__description_ggd___ggd___temperature[])
    var"time" :: Union{Nothing, Float64} = nothing
    var"power_density" :: FDSvector{T} where {T<:wall__description_ggd___ggd___power_density} = FDSvector(wall__description_ggd___ggd___power_density[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd___ggd(var"temperature"=FDSvector(wall__description_ggd___ggd___temperature[]), var"time"=nothing, var"power_density"=FDSvector(wall__description_ggd___ggd___power_density[]), _parent=WeakRef(nothing))
        obj = new(var"temperature", var"time", var"power_density", _parent)
        obj.temperature._parent = WeakRef(obj)
        obj.power_density._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_ggd <: FDS
    var"grid_ggd" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd} = FDSvector(wall__description_ggd___grid_ggd[])
    var"type" :: wall__description_ggd___type = wall__description_ggd___type()
    var"ggd" :: FDSvector{T} where {T<:wall__description_ggd___ggd} = FDSvector(wall__description_ggd___ggd[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_ggd(var"grid_ggd"=FDSvector(wall__description_ggd___grid_ggd[]), var"type"=wall__description_ggd___type(), var"ggd"=FDSvector(wall__description_ggd___ggd[]), _parent=WeakRef(nothing))
        obj = new(var"grid_ggd", var"type", var"ggd", _parent)
        obj.grid_ggd._parent = WeakRef(obj)
        obj.type._parent = WeakRef(obj)
        obj.ggd._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___outline <: FDS
    var"closed" :: Union{Nothing, Int} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___element___outline(var"closed"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"closed", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___j_tor <: FDS
    var"time" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"data" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___element___j_tor(var"time"=nothing, var"data"=nothing, _parent=WeakRef(nothing))
        obj = new(var"time", var"data", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"j_tor" :: wall__description_2d___vessel__unit___element___j_tor = wall__description_2d___vessel__unit___element___j_tor()
    var"resistivity" :: Union{Nothing, Float64} = nothing
    var"outline" :: wall__description_2d___vessel__unit___element___outline = wall__description_2d___vessel__unit___element___outline()
    var"resistance" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___element(var"name"=nothing, var"j_tor"=wall__description_2d___vessel__unit___element___j_tor(), var"resistivity"=nothing, var"outline"=wall__description_2d___vessel__unit___element___outline(), var"resistance"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"j_tor", var"resistivity", var"outline", var"resistance", _parent)
        obj.j_tor._parent = WeakRef(obj)
        obj.outline._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_outer <: FDS
    var"closed" :: Union{Nothing, Int} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___annular__outline_outer(var"closed"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"closed", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_inner <: FDS
    var"closed" :: Union{Nothing, Int} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___annular__outline_inner(var"closed"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"closed", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__centreline <: FDS
    var"closed" :: Union{Nothing, Int} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___annular__centreline(var"closed"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"closed", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular <: FDS
    var"outline_inner" :: wall__description_2d___vessel__unit___annular__outline_inner = wall__description_2d___vessel__unit___annular__outline_inner()
    var"centreline" :: wall__description_2d___vessel__unit___annular__centreline = wall__description_2d___vessel__unit___annular__centreline()
    var"thickness" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"resistivity" :: Union{Nothing, Float64} = nothing
    var"outline_outer" :: wall__description_2d___vessel__unit___annular__outline_outer = wall__description_2d___vessel__unit___annular__outline_outer()
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit___annular(var"outline_inner"=wall__description_2d___vessel__unit___annular__outline_inner(), var"centreline"=wall__description_2d___vessel__unit___annular__centreline(), var"thickness"=nothing, var"resistivity"=nothing, var"outline_outer"=wall__description_2d___vessel__unit___annular__outline_outer(), _parent=WeakRef(nothing))
        obj = new(var"outline_inner", var"centreline", var"thickness", var"resistivity", var"outline_outer", _parent)
        obj.outline_inner._parent = WeakRef(obj)
        obj.centreline._parent = WeakRef(obj)
        obj.outline_outer._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"annular" :: wall__description_2d___vessel__unit___annular = wall__description_2d___vessel__unit___annular()
    var"identifier" :: Union{Nothing, String} = nothing
    var"element" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit___element} = FDSvector(wall__description_2d___vessel__unit___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__unit(var"name"=nothing, var"annular"=wall__description_2d___vessel__unit___annular(), var"identifier"=nothing, var"element"=FDSvector(wall__description_2d___vessel__unit___element[]), _parent=WeakRef(nothing))
        obj = new(var"name", var"annular", var"identifier", var"element", _parent)
        obj.annular._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel__type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel <: FDS
    var"type" :: wall__description_2d___vessel__type = wall__description_2d___vessel__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit} = FDSvector(wall__description_2d___vessel__unit[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___vessel(var"type"=wall__description_2d___vessel__type(), var"unit"=FDSvector(wall__description_2d___vessel__unit[]), _parent=WeakRef(nothing))
        obj = new(var"type", var"unit", _parent)
        obj.type._parent = WeakRef(obj)
        obj.unit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit___outline <: FDS
    var"time" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___mobile__unit___outline(var"time"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"time", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"resistivity" :: Union{Nothing, Float64} = nothing
    var"closed" :: Union{Nothing, Int} = nothing
    var"phi_extensions" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"outline" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit___outline} = FDSvector(wall__description_2d___mobile__unit___outline[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___mobile__unit(var"name"=nothing, var"resistivity"=nothing, var"closed"=nothing, var"phi_extensions"=nothing, var"outline"=FDSvector(wall__description_2d___mobile__unit___outline[]), _parent=WeakRef(nothing))
        obj = new(var"name", var"resistivity", var"closed", var"phi_extensions", var"outline", _parent)
        obj.outline._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___mobile__type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile <: FDS
    var"type" :: wall__description_2d___mobile__type = wall__description_2d___mobile__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit} = FDSvector(wall__description_2d___mobile__unit[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___mobile(var"type"=wall__description_2d___mobile__type(), var"unit"=FDSvector(wall__description_2d___mobile__unit[]), _parent=WeakRef(nothing))
        obj = new(var"type", var"unit", _parent)
        obj.type._parent = WeakRef(obj)
        obj.unit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit___outline <: FDS
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___limiter__unit___outline(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"resistivity" :: Union{Nothing, Float64} = nothing
    var"closed" :: Union{Nothing, Int} = nothing
    var"outline" :: wall__description_2d___limiter__unit___outline = wall__description_2d___limiter__unit___outline()
    var"phi_extensions" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___limiter__unit(var"name"=nothing, var"resistivity"=nothing, var"closed"=nothing, var"outline"=wall__description_2d___limiter__unit___outline(), var"phi_extensions"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"resistivity", var"closed", var"outline", var"phi_extensions", _parent)
        obj.outline._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___limiter__type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter <: FDS
    var"type" :: wall__description_2d___limiter__type = wall__description_2d___limiter__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___limiter__unit} = FDSvector(wall__description_2d___limiter__unit[])
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d___limiter(var"type"=wall__description_2d___limiter__type(), var"unit"=FDSvector(wall__description_2d___limiter__unit[]), _parent=WeakRef(nothing))
        obj = new(var"type", var"unit", _parent)
        obj.type._parent = WeakRef(obj)
        obj.unit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__description_2d <: FDS
    var"mobile" :: wall__description_2d___mobile = wall__description_2d___mobile()
    var"limiter" :: wall__description_2d___limiter = wall__description_2d___limiter()
    var"type" :: wall__description_2d___type = wall__description_2d___type()
    var"vessel" :: wall__description_2d___vessel = wall__description_2d___vessel()
    _parent :: WeakRef = WeakRef(nothing)
    function wall__description_2d(var"mobile"=wall__description_2d___mobile(), var"limiter"=wall__description_2d___limiter(), var"type"=wall__description_2d___type(), var"vessel"=wall__description_2d___vessel(), _parent=WeakRef(nothing))
        obj = new(var"mobile", var"limiter", var"type", var"vessel", _parent)
        obj.mobile._parent = WeakRef(obj)
        obj.limiter._parent = WeakRef(obj)
        obj.type._parent = WeakRef(obj)
        obj.vessel._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall__code__library <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"commit" :: Union{Nothing, String} = nothing
    var"repository" :: Union{Nothing, String} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__code__library(var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"version"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)

        return obj
    end
end

Base.@kwdef mutable struct wall__code <: FDS
    var"library" :: FDSvector{T} where {T<:wall__code__library} = FDSvector(wall__code__library[])
    var"name" :: Union{Nothing, String} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"commit" :: Union{Nothing, String} = nothing
    var"repository" :: Union{Nothing, String} = nothing
    var"output_flag" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function wall__code(var"library"=FDSvector(wall__code__library[]), var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"output_flag"=nothing, var"version"=nothing, _parent=WeakRef(nothing))
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        obj.library._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct wall <: FDS
    var"time" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ids_properties" :: wall__ids_properties = wall__ids_properties()
    var"description_ggd" :: FDSvector{T} where {T<:wall__description_ggd} = FDSvector(wall__description_ggd[])
    var"description_2d" :: FDSvector{T} where {T<:wall__description_2d} = FDSvector(wall__description_2d[])
    var"first_wall_surface_area" :: Union{Nothing, Float64} = nothing
    var"code" :: wall__code = wall__code()
    var"global_quantities" :: wall__global_quantities = wall__global_quantities()
    var"temperature_reference" :: wall__temperature_reference = wall__temperature_reference()
    var"first_wall_power_flux_peak" :: wall__first_wall_power_flux_peak = wall__first_wall_power_flux_peak()
    _parent :: WeakRef = WeakRef(nothing)
    function wall(var"time"=nothing, var"ids_properties"=wall__ids_properties(), var"description_ggd"=FDSvector(wall__description_ggd[]), var"description_2d"=FDSvector(wall__description_2d[]), var"first_wall_surface_area"=nothing, var"code"=wall__code(), var"global_quantities"=wall__global_quantities(), var"temperature_reference"=wall__temperature_reference(), var"first_wall_power_flux_peak"=wall__first_wall_power_flux_peak(), _parent=WeakRef(nothing))
        obj = new(var"time", var"ids_properties", var"description_ggd", var"description_2d", var"first_wall_surface_area", var"code", var"global_quantities", var"temperature_reference", var"first_wall_power_flux_peak", _parent)
        obj.ids_properties._parent = WeakRef(obj)
        obj.description_ggd._parent = WeakRef(obj)
        obj.description_2d._parent = WeakRef(obj)
        obj.code._parent = WeakRef(obj)
        obj.global_quantities._parent = WeakRef(obj)
        obj.temperature_reference._parent = WeakRef(obj)
        obj.first_wall_power_flux_peak._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"r0" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__vacuum_toroidal_field(var"b0"=nothing, var"r0"=nothing, _parent=WeakRef(nothing))
        obj = new(var"b0", var"r0", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___profiles_2d___grid_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid <: FDS
    var"volume_element" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"dim2" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dim1" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___profiles_2d___grid(var"volume_element"=nothing, var"dim2"=nothing, var"dim1"=nothing, _parent=WeakRef(nothing))
        obj = new(var"volume_element", var"dim2", var"dim1", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d <: FDS
    var"psi" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"b_field_r" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"b_r" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"theta" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"b_field_z" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"j_tor" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"phi" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"b_field_tor" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"b_z" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"grid" :: equilibrium__time_slice___profiles_2d___grid = equilibrium__time_slice___profiles_2d___grid()
    var"grid_type" :: equilibrium__time_slice___profiles_2d___grid_type = equilibrium__time_slice___profiles_2d___grid_type()
    var"j_parallel" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"b_tor" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___profiles_2d(var"psi"=nothing, var"b_field_r"=nothing, var"r"=nothing, var"b_r"=nothing, var"theta"=nothing, var"b_field_z"=nothing, var"j_tor"=nothing, var"phi"=nothing, var"z"=nothing, var"b_field_tor"=nothing, var"b_z"=nothing, var"grid"=equilibrium__time_slice___profiles_2d___grid(), var"grid_type"=equilibrium__time_slice___profiles_2d___grid_type(), var"j_parallel"=nothing, var"b_tor"=nothing, _parent=WeakRef(nothing))
        obj = new(var"psi", var"b_field_r", var"r", var"b_r", var"theta", var"b_field_z", var"j_tor", var"phi", var"z", var"b_field_tor", var"b_z", var"grid", var"grid_type", var"j_parallel", var"b_tor", _parent)
        obj.grid._parent = WeakRef(obj)
        obj.grid_type._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d__geometric_axis <: FDS
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___profiles_1d__geometric_axis(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d <: FDS
    var"b_field_max" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dvolume_drho_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm9" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dpsi_drho_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"surface" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"magnetic_shear" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"b_average" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"b_field_min" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"darea_dpsi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm3" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"squareness_upper_inner" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"squareness_lower_inner" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"elongation" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"beta_pol" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"b_field_average" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"j_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm6" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"psi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm8" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dpressure_dpsi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"triangularity_upper" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"darea_drho_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"area" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"trapped_fraction" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"volume" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dvolume_dpsi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"b_min" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"f" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"mass_density" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"r_outboard" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm4" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"phi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"squareness_lower_outer" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"triangularity_lower" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm2" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_volume_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm1" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm5" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"b_max" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"f_df_dpsi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"j_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"r_inboard" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"q" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"gm7" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"squareness_upper_outer" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"geometric_axis" :: equilibrium__time_slice___profiles_1d__geometric_axis = equilibrium__time_slice___profiles_1d__geometric_axis()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___profiles_1d(var"b_field_max"=nothing, var"dvolume_drho_tor"=nothing, var"gm9"=nothing, var"dpsi_drho_tor"=nothing, var"surface"=nothing, var"rho_tor"=nothing, var"magnetic_shear"=nothing, var"b_average"=nothing, var"b_field_min"=nothing, var"darea_dpsi"=nothing, var"gm3"=nothing, var"squareness_upper_inner"=nothing, var"squareness_lower_inner"=nothing, var"rho_tor_norm"=nothing, var"elongation"=nothing, var"beta_pol"=nothing, var"b_field_average"=nothing, var"j_parallel"=nothing, var"gm6"=nothing, var"psi"=nothing, var"gm8"=nothing, var"dpressure_dpsi"=nothing, var"triangularity_upper"=nothing, var"darea_drho_tor"=nothing, var"area"=nothing, var"trapped_fraction"=nothing, var"volume"=nothing, var"dvolume_dpsi"=nothing, var"b_min"=nothing, var"f"=nothing, var"mass_density"=nothing, var"r_outboard"=nothing, var"gm4"=nothing, var"phi"=nothing, var"squareness_lower_outer"=nothing, var"triangularity_lower"=nothing, var"gm2"=nothing, var"rho_volume_norm"=nothing, var"gm1"=nothing, var"gm5"=nothing, var"b_max"=nothing, var"f_df_dpsi"=nothing, var"j_tor"=nothing, var"r_inboard"=nothing, var"q"=nothing, var"gm7"=nothing, var"pressure"=nothing, var"squareness_upper_outer"=nothing, var"geometric_axis"=equilibrium__time_slice___profiles_1d__geometric_axis(), _parent=WeakRef(nothing))
        obj = new(var"b_field_max", var"dvolume_drho_tor", var"gm9", var"dpsi_drho_tor", var"surface", var"rho_tor", var"magnetic_shear", var"b_average", var"b_field_min", var"darea_dpsi", var"gm3", var"squareness_upper_inner", var"squareness_lower_inner", var"rho_tor_norm", var"elongation", var"beta_pol", var"b_field_average", var"j_parallel", var"gm6", var"psi", var"gm8", var"dpressure_dpsi", var"triangularity_upper", var"darea_drho_tor", var"area", var"trapped_fraction", var"volume", var"dvolume_dpsi", var"b_min", var"f", var"mass_density", var"r_outboard", var"gm4", var"phi", var"squareness_lower_outer", var"triangularity_lower", var"gm2", var"rho_volume_norm", var"gm1", var"gm5", var"b_max", var"f_df_dpsi", var"j_tor", var"r_inboard", var"q", var"gm7", var"pressure", var"squareness_upper_outer", var"geometric_axis", _parent)
        obj.geometric_axis._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__q_min <: FDS
    var"value" :: Union{Nothing, Float64} = nothing
    var"rho_tor_norm" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___global_quantities__q_min(var"value"=nothing, var"rho_tor_norm"=nothing, _parent=WeakRef(nothing))
        obj = new(var"value", var"rho_tor_norm", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__magnetic_axis <: FDS
    var"b_field_tor" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    var"b_tor" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___global_quantities__magnetic_axis(var"b_field_tor"=nothing, var"r"=nothing, var"z"=nothing, var"b_tor"=nothing, _parent=WeakRef(nothing))
        obj = new(var"b_field_tor", var"r", var"z", var"b_tor", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__current_centre <: FDS
    var"velocity_z" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___global_quantities__current_centre(var"velocity_z"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"velocity_z", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities <: FDS
    var"ip" :: Union{Nothing, Float64} = nothing
    var"li_3" :: Union{Nothing, Float64} = nothing
    var"beta_tor" :: Union{Nothing, Float64} = nothing
    var"surface" :: Union{Nothing, Float64} = nothing
    var"magnetic_axis" :: equilibrium__time_slice___global_quantities__magnetic_axis = equilibrium__time_slice___global_quantities__magnetic_axis()
    var"energy_mhd" :: Union{Nothing, Float64} = nothing
    var"psi_boundary" :: Union{Nothing, Float64} = nothing
    var"length_pol" :: Union{Nothing, Float64} = nothing
    var"area" :: Union{Nothing, Float64} = nothing
    var"psi_external_average" :: Union{Nothing, Float64} = nothing
    var"q_95" :: Union{Nothing, Float64} = nothing
    var"q_axis" :: Union{Nothing, Float64} = nothing
    var"psi_axis" :: Union{Nothing, Float64} = nothing
    var"w_mhd" :: Union{Nothing, Float64} = nothing
    var"volume" :: Union{Nothing, Float64} = nothing
    var"plasma_inductance" :: Union{Nothing, Float64} = nothing
    var"beta_pol" :: Union{Nothing, Float64} = nothing
    var"beta_normal" :: Union{Nothing, Float64} = nothing
    var"current_centre" :: equilibrium__time_slice___global_quantities__current_centre = equilibrium__time_slice___global_quantities__current_centre()
    var"q_min" :: equilibrium__time_slice___global_quantities__q_min = equilibrium__time_slice___global_quantities__q_min()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___global_quantities(var"ip"=nothing, var"li_3"=nothing, var"beta_tor"=nothing, var"surface"=nothing, var"magnetic_axis"=equilibrium__time_slice___global_quantities__magnetic_axis(), var"energy_mhd"=nothing, var"psi_boundary"=nothing, var"length_pol"=nothing, var"area"=nothing, var"psi_external_average"=nothing, var"q_95"=nothing, var"q_axis"=nothing, var"psi_axis"=nothing, var"w_mhd"=nothing, var"volume"=nothing, var"plasma_inductance"=nothing, var"beta_pol"=nothing, var"beta_normal"=nothing, var"current_centre"=equilibrium__time_slice___global_quantities__current_centre(), var"q_min"=equilibrium__time_slice___global_quantities__q_min(), _parent=WeakRef(nothing))
        obj = new(var"ip", var"li_3", var"beta_tor", var"surface", var"magnetic_axis", var"energy_mhd", var"psi_boundary", var"length_pol", var"area", var"psi_external_average", var"q_95", var"q_axis", var"psi_axis", var"w_mhd", var"volume", var"plasma_inductance", var"beta_pol", var"beta_normal", var"current_centre", var"q_min", _parent)
        obj.magnetic_axis._parent = WeakRef(obj)
        obj.current_centre._parent = WeakRef(obj)
        obj.q_min._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___z <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___z(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___theta <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___theta(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___r <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___r(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___psi <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___psi(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___phi <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___phi(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_tor <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___j_tor(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_parallel <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___j_parallel(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary <: FDS
    var"neighbours" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary(var"neighbours"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"neighbours", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object <: FDS
    var"nodes" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"measure" :: Union{Nothing, Float64} = nothing
    var"geometry" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object(var"nodes"=nothing, var"measure"=nothing, var"geometry"=nothing, var"boundary"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[]), _parent=WeakRef(nothing))
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        obj.boundary._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[]), _parent=WeakRef(nothing))
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__space___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___geometry_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__space___geometry_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space <: FDS
    var"coordinates_type" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"geometry_type" :: equilibrium__time_slice___ggd___grid__space___geometry_type = equilibrium__time_slice___ggd___grid__space___geometry_type()
    var"identifier" :: equilibrium__time_slice___ggd___grid__space___identifier = equilibrium__time_slice___ggd___grid__space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__space(var"coordinates_type"=nothing, var"geometry_type"=equilibrium__time_slice___ggd___grid__space___geometry_type(), var"identifier"=equilibrium__time_slice___ggd___grid__space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[]), _parent=WeakRef(nothing))
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        obj.geometry_type._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.objects_per_dimension._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___metric <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__grid_subset___metric(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__grid_subset___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element___object <: FDS
    var"dimension" :: Union{Nothing, Int} = nothing
    var"space" :: Union{Nothing, Int} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element___object(var"dimension"=nothing, var"space"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"dimension", var"space", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element___object} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[]), _parent=WeakRef(nothing))
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___base <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__grid_subset___base(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___base} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[])
    var"metric" :: equilibrium__time_slice___ggd___grid__grid_subset___metric = equilibrium__time_slice___ggd___grid__grid_subset___metric()
    var"dimension" :: Union{Nothing, Int} = nothing
    var"identifier" :: equilibrium__time_slice___ggd___grid__grid_subset___identifier = equilibrium__time_slice___ggd___grid__grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid__grid_subset(var"base"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[]), var"metric"=equilibrium__time_slice___ggd___grid__grid_subset___metric(), var"dimension"=nothing, var"identifier"=equilibrium__time_slice___ggd___grid__grid_subset___identifier(), var"element"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[]), _parent=WeakRef(nothing))
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        obj.base._parent = WeakRef(obj)
        obj.metric._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space} = FDSvector(equilibrium__time_slice___ggd___grid__space[])
    var"identifier" :: equilibrium__time_slice___ggd___grid__identifier = equilibrium__time_slice___ggd___grid__identifier()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___grid(var"grid_subset"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset[]), var"space"=FDSvector(equilibrium__time_slice___ggd___grid__space[]), var"identifier"=equilibrium__time_slice___ggd___grid__identifier(), _parent=WeakRef(nothing))
        obj = new(var"grid_subset", var"space", var"identifier", _parent)
        obj.grid_subset._parent = WeakRef(obj)
        obj.space._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_z <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___b_field_z(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_tor <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___b_field_tor(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_r <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd___b_field_r(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=WeakRef(nothing))
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
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
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___ggd(var"b_field_z"=FDSvector(equilibrium__time_slice___ggd___b_field_z[]), var"psi"=FDSvector(equilibrium__time_slice___ggd___psi[]), var"theta"=FDSvector(equilibrium__time_slice___ggd___theta[]), var"z"=FDSvector(equilibrium__time_slice___ggd___z[]), var"phi"=FDSvector(equilibrium__time_slice___ggd___phi[]), var"j_tor"=FDSvector(equilibrium__time_slice___ggd___j_tor[]), var"grid"=equilibrium__time_slice___ggd___grid(), var"b_field_tor"=FDSvector(equilibrium__time_slice___ggd___b_field_tor[]), var"b_field_r"=FDSvector(equilibrium__time_slice___ggd___b_field_r[]), var"r"=FDSvector(equilibrium__time_slice___ggd___r[]), var"j_parallel"=FDSvector(equilibrium__time_slice___ggd___j_parallel[]), _parent=WeakRef(nothing))
        obj = new(var"b_field_z", var"psi", var"theta", var"z", var"phi", var"j_tor", var"grid", var"b_field_tor", var"b_field_r", var"r", var"j_parallel", _parent)
        obj.b_field_z._parent = WeakRef(obj)
        obj.psi._parent = WeakRef(obj)
        obj.theta._parent = WeakRef(obj)
        obj.z._parent = WeakRef(obj)
        obj.phi._parent = WeakRef(obj)
        obj.j_tor._parent = WeakRef(obj)
        obj.grid._parent = WeakRef(obj)
        obj.b_field_tor._parent = WeakRef(obj)
        obj.b_field_r._parent = WeakRef(obj)
        obj.r._parent = WeakRef(obj)
        obj.j_parallel._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___coordinate_system__grid_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid <: FDS
    var"volume_element" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"dim2" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dim1" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___coordinate_system__grid(var"volume_element"=nothing, var"dim2"=nothing, var"dim1"=nothing, _parent=WeakRef(nothing))
        obj = new(var"volume_element", var"dim2", var"dim1", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g13_covariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g11_contravariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g13_contravariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"r" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g12_contravariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g22_contravariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g33_contravariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g22_covariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 4}} = nothing
    var"g12_covariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g33_covariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"grid" :: equilibrium__time_slice___coordinate_system__grid = equilibrium__time_slice___coordinate_system__grid()
    var"g23_covariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"g11_covariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 4}} = nothing
    var"g23_contravariant" :: Union{Nothing, AbstractArray{Float64, 2}} = nothing
    var"grid_type" :: equilibrium__time_slice___coordinate_system__grid_type = equilibrium__time_slice___coordinate_system__grid_type()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___coordinate_system(var"jacobian"=nothing, var"g13_covariant"=nothing, var"g11_contravariant"=nothing, var"g13_contravariant"=nothing, var"r"=nothing, var"g12_contravariant"=nothing, var"g22_contravariant"=nothing, var"z"=nothing, var"g33_contravariant"=nothing, var"g22_covariant"=nothing, var"tensor_contravariant"=nothing, var"g12_covariant"=nothing, var"g33_covariant"=nothing, var"grid"=equilibrium__time_slice___coordinate_system__grid(), var"g23_covariant"=nothing, var"g11_covariant"=nothing, var"tensor_covariant"=nothing, var"g23_contravariant"=nothing, var"grid_type"=equilibrium__time_slice___coordinate_system__grid_type(), _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"g13_covariant", var"g11_contravariant", var"g13_contravariant", var"r", var"g12_contravariant", var"g22_contravariant", var"z", var"g33_contravariant", var"g22_covariant", var"tensor_contravariant", var"g12_covariant", var"g33_covariant", var"grid", var"g23_covariant", var"g11_covariant", var"tensor_covariant", var"g23_contravariant", var"grid_type", _parent)
        obj.grid._parent = WeakRef(obj)
        obj.grid_type._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___convergence <: FDS
    var"iterations_n" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___convergence(var"iterations_n"=nothing, _parent=WeakRef(nothing))
        obj = new(var"iterations_n", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_reconstructed <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__x_point___position_reconstructed(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_measured <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__x_point___position_measured(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point <: FDS
    var"chi_squared_z" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"position_measured" :: equilibrium__time_slice___constraints__x_point___position_measured = equilibrium__time_slice___constraints__x_point___position_measured()
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    var"chi_squared_r" :: Union{Nothing, Float64} = nothing
    var"position_reconstructed" :: equilibrium__time_slice___constraints__x_point___position_reconstructed = equilibrium__time_slice___constraints__x_point___position_reconstructed()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__x_point(var"chi_squared_z"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"position_measured"=equilibrium__time_slice___constraints__x_point___position_measured(), var"time_measurement"=nothing, var"chi_squared_r"=nothing, var"position_reconstructed"=equilibrium__time_slice___constraints__x_point___position_reconstructed(), _parent=WeakRef(nothing))
        obj = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        obj.position_measured._parent = WeakRef(obj)
        obj.position_reconstructed._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_reconstructed <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__strike_point___position_reconstructed(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_measured <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__strike_point___position_measured(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point <: FDS
    var"chi_squared_z" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"position_measured" :: equilibrium__time_slice___constraints__strike_point___position_measured = equilibrium__time_slice___constraints__strike_point___position_measured()
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    var"chi_squared_r" :: Union{Nothing, Float64} = nothing
    var"position_reconstructed" :: equilibrium__time_slice___constraints__strike_point___position_reconstructed = equilibrium__time_slice___constraints__strike_point___position_reconstructed()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__strike_point(var"chi_squared_z"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"position_measured"=equilibrium__time_slice___constraints__strike_point___position_measured(), var"time_measurement"=nothing, var"chi_squared_r"=nothing, var"position_reconstructed"=equilibrium__time_slice___constraints__strike_point___position_reconstructed(), _parent=WeakRef(nothing))
        obj = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        obj.position_measured._parent = WeakRef(obj)
        obj.position_reconstructed._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q___position <: FDS
    var"phi" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__q___position(var"phi"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"phi", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"position" :: equilibrium__time_slice___constraints__q___position = equilibrium__time_slice___constraints__q___position()
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__q(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"position"=equilibrium__time_slice___constraints__q___position(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"position", var"time_measurement", _parent)
        obj.position._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pressure <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__pressure(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_passive_current <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__pf_passive_current(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_current <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__pf_current(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e_line <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__n_e_line(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__n_e(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__mse_polarisation_angle <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__mse_polarisation_angle(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment <: FDS
    var"magnetisation_r" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r()
    var"magnetisation_z" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__iron_core_segment(var"magnetisation_r"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(), var"magnetisation_z"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(), _parent=WeakRef(nothing))
        obj = new(var"magnetisation_r", var"magnetisation_z", _parent)
        obj.magnetisation_r._parent = WeakRef(obj)
        obj.magnetisation_z._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__ip <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__ip(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__flux_loop <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__flux_loop(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__faraday_angle <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__faraday_angle(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__diamagnetic_flux <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__diamagnetic_flux(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__bpol_probe <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__bpol_probe(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__b_field_tor_vacuum_r <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints__b_field_tor_vacuum_r(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
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
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___constraints(var"faraday_angle"=FDSvector(equilibrium__time_slice___constraints__faraday_angle[]), var"n_e"=FDSvector(equilibrium__time_slice___constraints__n_e[]), var"ip"=equilibrium__time_slice___constraints__ip(), var"n_e_line"=FDSvector(equilibrium__time_slice___constraints__n_e_line[]), var"pf_current"=FDSvector(equilibrium__time_slice___constraints__pf_current[]), var"strike_point"=FDSvector(equilibrium__time_slice___constraints__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___constraints__x_point[]), var"iron_core_segment"=FDSvector(equilibrium__time_slice___constraints__iron_core_segment[]), var"pressure"=FDSvector(equilibrium__time_slice___constraints__pressure[]), var"diamagnetic_flux"=equilibrium__time_slice___constraints__diamagnetic_flux(), var"pf_passive_current"=FDSvector(equilibrium__time_slice___constraints__pf_passive_current[]), var"bpol_probe"=FDSvector(equilibrium__time_slice___constraints__bpol_probe[]), var"mse_polarisation_angle"=FDSvector(equilibrium__time_slice___constraints__mse_polarisation_angle[]), var"q"=FDSvector(equilibrium__time_slice___constraints__q[]), var"b_field_tor_vacuum_r"=equilibrium__time_slice___constraints__b_field_tor_vacuum_r(), var"flux_loop"=FDSvector(equilibrium__time_slice___constraints__flux_loop[]), _parent=WeakRef(nothing))
        obj = new(var"faraday_angle", var"n_e", var"ip", var"n_e_line", var"pf_current", var"strike_point", var"x_point", var"iron_core_segment", var"pressure", var"diamagnetic_flux", var"pf_passive_current", var"bpol_probe", var"mse_polarisation_angle", var"q", var"b_field_tor_vacuum_r", var"flux_loop", _parent)
        obj.faraday_angle._parent = WeakRef(obj)
        obj.n_e._parent = WeakRef(obj)
        obj.ip._parent = WeakRef(obj)
        obj.n_e_line._parent = WeakRef(obj)
        obj.pf_current._parent = WeakRef(obj)
        obj.strike_point._parent = WeakRef(obj)
        obj.x_point._parent = WeakRef(obj)
        obj.iron_core_segment._parent = WeakRef(obj)
        obj.pressure._parent = WeakRef(obj)
        obj.diamagnetic_flux._parent = WeakRef(obj)
        obj.pf_passive_current._parent = WeakRef(obj)
        obj.bpol_probe._parent = WeakRef(obj)
        obj.mse_polarisation_angle._parent = WeakRef(obj)
        obj.q._parent = WeakRef(obj)
        obj.b_field_tor_vacuum_r._parent = WeakRef(obj)
        obj.flux_loop._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__x_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__x_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__strike_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__strike_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__outline <: FDS
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__outline(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__geometric_axis <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__geometric_axis(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__gap <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"value" :: Union{Nothing, Float64} = nothing
    var"identifier" :: Union{Nothing, String} = nothing
    var"angle" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__gap(var"name"=nothing, var"r"=nothing, var"value"=nothing, var"identifier"=nothing, var"angle"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"r", var"value", var"identifier", var"angle", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__closest_wall_point <: FDS
    var"distance" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__closest_wall_point(var"distance"=nothing, var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"distance", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__active_limiter_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix__active_limiter_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix <: FDS
    var"psi" :: Union{Nothing, Float64} = nothing
    var"elongation_lower" :: Union{Nothing, Float64} = nothing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_separatrix__x_point[])
    var"gap" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_separatrix__gap} = FDSvector(equilibrium__time_slice___boundary_separatrix__gap[])
    var"triangularity" :: Union{Nothing, Float64} = nothing
    var"elongation_upper" :: Union{Nothing, Float64} = nothing
    var"triangularity_upper" :: Union{Nothing, Float64} = nothing
    var"outline" :: equilibrium__time_slice___boundary_separatrix__outline = equilibrium__time_slice___boundary_separatrix__outline()
    var"dr_dz_zero_point" :: equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point = equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point()
    var"squareness_lower_outer" :: Union{Nothing, Float64} = nothing
    var"triangularity_lower" :: Union{Nothing, Float64} = nothing
    var"minor_radius" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_inner" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_outer" :: Union{Nothing, Float64} = nothing
    var"squareness_lower_inner" :: Union{Nothing, Float64} = nothing
    var"geometric_axis" :: equilibrium__time_slice___boundary_separatrix__geometric_axis = equilibrium__time_slice___boundary_separatrix__geometric_axis()
    var"elongation" :: Union{Nothing, Float64} = nothing
    var"active_limiter_point" :: equilibrium__time_slice___boundary_separatrix__active_limiter_point = equilibrium__time_slice___boundary_separatrix__active_limiter_point()
    var"closest_wall_point" :: equilibrium__time_slice___boundary_separatrix__closest_wall_point = equilibrium__time_slice___boundary_separatrix__closest_wall_point()
    var"type" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_separatrix(var"psi"=nothing, var"elongation_lower"=nothing, var"strike_point"=FDSvector(equilibrium__time_slice___boundary_separatrix__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___boundary_separatrix__x_point[]), var"gap"=FDSvector(equilibrium__time_slice___boundary_separatrix__gap[]), var"triangularity"=nothing, var"elongation_upper"=nothing, var"triangularity_upper"=nothing, var"outline"=equilibrium__time_slice___boundary_separatrix__outline(), var"dr_dz_zero_point"=equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(), var"squareness_lower_outer"=nothing, var"triangularity_lower"=nothing, var"minor_radius"=nothing, var"squareness_upper_inner"=nothing, var"squareness_upper_outer"=nothing, var"squareness_lower_inner"=nothing, var"geometric_axis"=equilibrium__time_slice___boundary_separatrix__geometric_axis(), var"elongation"=nothing, var"active_limiter_point"=equilibrium__time_slice___boundary_separatrix__active_limiter_point(), var"closest_wall_point"=equilibrium__time_slice___boundary_separatrix__closest_wall_point(), var"type"=nothing, _parent=WeakRef(nothing))
        obj = new(var"psi", var"elongation_lower", var"strike_point", var"x_point", var"gap", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"dr_dz_zero_point", var"squareness_lower_outer", var"triangularity_lower", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"closest_wall_point", var"type", _parent)
        obj.strike_point._parent = WeakRef(obj)
        obj.x_point._parent = WeakRef(obj)
        obj.gap._parent = WeakRef(obj)
        obj.outline._parent = WeakRef(obj)
        obj.dr_dz_zero_point._parent = WeakRef(obj)
        obj.geometric_axis._parent = WeakRef(obj)
        obj.active_limiter_point._parent = WeakRef(obj)
        obj.closest_wall_point._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__x_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_secondary_separatrix__x_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__strike_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_secondary_separatrix__strike_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__outline <: FDS
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_secondary_separatrix__outline(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix <: FDS
    var"psi" :: Union{Nothing, Float64} = nothing
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[])
    var"distance_inner_outer" :: Union{Nothing, Float64} = nothing
    var"outline" :: equilibrium__time_slice___boundary_secondary_separatrix__outline = equilibrium__time_slice___boundary_secondary_separatrix__outline()
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary_secondary_separatrix(var"psi"=nothing, var"x_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[]), var"distance_inner_outer"=nothing, var"outline"=equilibrium__time_slice___boundary_secondary_separatrix__outline(), var"strike_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[]), _parent=WeakRef(nothing))
        obj = new(var"psi", var"x_point", var"distance_inner_outer", var"outline", var"strike_point", _parent)
        obj.x_point._parent = WeakRef(obj)
        obj.outline._parent = WeakRef(obj)
        obj.strike_point._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__x_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary__x_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__strike_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary__strike_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__outline <: FDS
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary__outline(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__lcfs <: FDS
    var"r" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary__lcfs(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__geometric_axis <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary__geometric_axis(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__active_limiter_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary__active_limiter_point(var"r"=nothing, var"z"=nothing, _parent=WeakRef(nothing))
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary <: FDS
    var"psi" :: Union{Nothing, Float64} = nothing
    var"lcfs" :: equilibrium__time_slice___boundary__lcfs = equilibrium__time_slice___boundary__lcfs()
    var"elongation_lower" :: Union{Nothing, Float64} = nothing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary__strike_point} = FDSvector(equilibrium__time_slice___boundary__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary__x_point} = FDSvector(equilibrium__time_slice___boundary__x_point[])
    var"triangularity" :: Union{Nothing, Float64} = nothing
    var"elongation_upper" :: Union{Nothing, Float64} = nothing
    var"triangularity_upper" :: Union{Nothing, Float64} = nothing
    var"outline" :: equilibrium__time_slice___boundary__outline = equilibrium__time_slice___boundary__outline()
    var"squareness_lower_outer" :: Union{Nothing, Float64} = nothing
    var"triangularity_lower" :: Union{Nothing, Float64} = nothing
    var"psi_norm" :: Union{Nothing, Float64} = nothing
    var"minor_radius" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_inner" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_outer" :: Union{Nothing, Float64} = nothing
    var"squareness_lower_inner" :: Union{Nothing, Float64} = nothing
    var"geometric_axis" :: equilibrium__time_slice___boundary__geometric_axis = equilibrium__time_slice___boundary__geometric_axis()
    var"elongation" :: Union{Nothing, Float64} = nothing
    var"active_limiter_point" :: equilibrium__time_slice___boundary__active_limiter_point = equilibrium__time_slice___boundary__active_limiter_point()
    var"b_flux_pol_norm" :: Union{Nothing, Float64} = nothing
    var"type" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice___boundary(var"psi"=nothing, var"lcfs"=equilibrium__time_slice___boundary__lcfs(), var"elongation_lower"=nothing, var"strike_point"=FDSvector(equilibrium__time_slice___boundary__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice___boundary__x_point[]), var"triangularity"=nothing, var"elongation_upper"=nothing, var"triangularity_upper"=nothing, var"outline"=equilibrium__time_slice___boundary__outline(), var"squareness_lower_outer"=nothing, var"triangularity_lower"=nothing, var"psi_norm"=nothing, var"minor_radius"=nothing, var"squareness_upper_inner"=nothing, var"squareness_upper_outer"=nothing, var"squareness_lower_inner"=nothing, var"geometric_axis"=equilibrium__time_slice___boundary__geometric_axis(), var"elongation"=nothing, var"active_limiter_point"=equilibrium__time_slice___boundary__active_limiter_point(), var"b_flux_pol_norm"=nothing, var"type"=nothing, _parent=WeakRef(nothing))
        obj = new(var"psi", var"lcfs", var"elongation_lower", var"strike_point", var"x_point", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"squareness_lower_outer", var"triangularity_lower", var"psi_norm", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"b_flux_pol_norm", var"type", _parent)
        obj.lcfs._parent = WeakRef(obj)
        obj.strike_point._parent = WeakRef(obj)
        obj.x_point._parent = WeakRef(obj)
        obj.outline._parent = WeakRef(obj)
        obj.geometric_axis._parent = WeakRef(obj)
        obj.active_limiter_point._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice <: FDS
    var"time" :: Union{Nothing, Float64} = nothing
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
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__time_slice(var"time"=nothing, var"ggd"=FDSvector(equilibrium__time_slice___ggd[]), var"profiles_1d"=equilibrium__time_slice___profiles_1d(), var"boundary"=equilibrium__time_slice___boundary(), var"constraints"=equilibrium__time_slice___constraints(), var"global_quantities"=equilibrium__time_slice___global_quantities(), var"convergence"=equilibrium__time_slice___convergence(), var"coordinate_system"=equilibrium__time_slice___coordinate_system(), var"boundary_secondary_separatrix"=equilibrium__time_slice___boundary_secondary_separatrix(), var"boundary_separatrix"=equilibrium__time_slice___boundary_separatrix(), var"profiles_2d"=FDSvector(equilibrium__time_slice___profiles_2d[]), _parent=WeakRef(nothing))
        obj = new(var"time", var"ggd", var"profiles_1d", var"boundary", var"constraints", var"global_quantities", var"convergence", var"coordinate_system", var"boundary_secondary_separatrix", var"boundary_separatrix", var"profiles_2d", _parent)
        obj.ggd._parent = WeakRef(obj)
        obj.profiles_1d._parent = WeakRef(obj)
        obj.boundary._parent = WeakRef(obj)
        obj.constraints._parent = WeakRef(obj)
        obj.global_quantities._parent = WeakRef(obj)
        obj.convergence._parent = WeakRef(obj)
        obj.coordinate_system._parent = WeakRef(obj)
        obj.boundary_secondary_separatrix._parent = WeakRef(obj)
        obj.boundary_separatrix._parent = WeakRef(obj)
        obj.profiles_2d._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Nothing, String} = nothing
    var"data_dictionary" :: Union{Nothing, String} = nothing
    var"access_layer" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__ids_properties__version_put(var"access_layer_language"=nothing, var"data_dictionary"=nothing, var"access_layer"=nothing, _parent=WeakRef(nothing))
        obj = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__ids_properties <: FDS
    var"provider" :: Union{Nothing, String} = nothing
    var"version_put" :: equilibrium__ids_properties__version_put = equilibrium__ids_properties__version_put()
    var"homogeneous_time" :: Union{Nothing, Int} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"creation_date" :: Union{Nothing, String} = nothing
    var"comment" :: Union{Nothing, String} = nothing
    var"occurrence" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__ids_properties(var"provider"=nothing, var"version_put"=equilibrium__ids_properties__version_put(), var"homogeneous_time"=nothing, var"source"=nothing, var"creation_date"=nothing, var"comment"=nothing, var"occurrence"=nothing, _parent=WeakRef(nothing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        obj.version_put._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary <: FDS
    var"neighbours" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary(var"neighbours"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"neighbours", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object <: FDS
    var"nodes" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"measure" :: Union{Nothing, Float64} = nothing
    var"geometry" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object(var"nodes"=nothing, var"measure"=nothing, var"geometry"=nothing, var"boundary"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(nothing))
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        obj.boundary._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension(var"object"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[]), _parent=WeakRef(nothing))
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___space___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___geometry_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___space___geometry_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space <: FDS
    var"coordinates_type" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"geometry_type" :: equilibrium__grids_ggd___grid___space___geometry_type = equilibrium__grids_ggd___grid___space___geometry_type()
    var"identifier" :: equilibrium__grids_ggd___grid___space___identifier = equilibrium__grids_ggd___grid___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___space(var"coordinates_type"=nothing, var"geometry_type"=equilibrium__grids_ggd___grid___space___geometry_type(), var"identifier"=equilibrium__grids_ggd___grid___space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[]), _parent=WeakRef(nothing))
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        obj.geometry_type._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.objects_per_dimension._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___metric <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___grid_subset___metric(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___grid_subset___identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element___object <: FDS
    var"dimension" :: Union{Nothing, Int} = nothing
    var"space" :: Union{Nothing, Int} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___grid_subset___element___object(var"dimension"=nothing, var"space"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"dimension", var"space", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element___object} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___grid_subset___element(var"object"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[]), _parent=WeakRef(nothing))
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___base <: FDS
    var"jacobian" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, AbstractArray{Float64, 3}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___grid_subset___base(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=WeakRef(nothing))
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___base} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[])
    var"metric" :: equilibrium__grids_ggd___grid___grid_subset___metric = equilibrium__grids_ggd___grid___grid_subset___metric()
    var"dimension" :: Union{Nothing, Int} = nothing
    var"identifier" :: equilibrium__grids_ggd___grid___grid_subset___identifier = equilibrium__grids_ggd___grid___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid___grid_subset(var"base"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[]), var"metric"=equilibrium__grids_ggd___grid___grid_subset___metric(), var"dimension"=nothing, var"identifier"=equilibrium__grids_ggd___grid___grid_subset___identifier(), var"element"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[]), _parent=WeakRef(nothing))
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        obj.base._parent = WeakRef(obj)
        obj.metric._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset} = FDSvector(equilibrium__grids_ggd___grid___grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space} = FDSvector(equilibrium__grids_ggd___grid___space[])
    var"identifier" :: equilibrium__grids_ggd___grid___identifier = equilibrium__grids_ggd___grid___identifier()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd___grid(var"grid_subset"=FDSvector(equilibrium__grids_ggd___grid___grid_subset[]), var"space"=FDSvector(equilibrium__grids_ggd___grid___space[]), var"identifier"=equilibrium__grids_ggd___grid___identifier(), _parent=WeakRef(nothing))
        obj = new(var"grid_subset", var"space", var"identifier", _parent)
        obj.grid_subset._parent = WeakRef(obj)
        obj.space._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd <: FDS
    var"time" :: Union{Nothing, Float64} = nothing
    var"grid" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid} = FDSvector(equilibrium__grids_ggd___grid[])
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__grids_ggd(var"time"=nothing, var"grid"=FDSvector(equilibrium__grids_ggd___grid[]), _parent=WeakRef(nothing))
        obj = new(var"time", var"grid", _parent)
        obj.grid._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__code__library <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"commit" :: Union{Nothing, String} = nothing
    var"repository" :: Union{Nothing, String} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__code__library(var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"version"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__code <: FDS
    var"library" :: FDSvector{T} where {T<:equilibrium__code__library} = FDSvector(equilibrium__code__library[])
    var"name" :: Union{Nothing, String} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"commit" :: Union{Nothing, String} = nothing
    var"repository" :: Union{Nothing, String} = nothing
    var"output_flag" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium__code(var"library"=FDSvector(equilibrium__code__library[]), var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"output_flag"=nothing, var"version"=nothing, _parent=WeakRef(nothing))
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        obj.library._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium <: FDS
    var"time_slice" :: FDSvector{T} where {T<:equilibrium__time_slice} = FDSvector(equilibrium__time_slice[])
    var"time" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ids_properties" :: equilibrium__ids_properties = equilibrium__ids_properties()
    var"grids_ggd" :: FDSvector{T} where {T<:equilibrium__grids_ggd} = FDSvector(equilibrium__grids_ggd[])
    var"vacuum_toroidal_field" :: equilibrium__vacuum_toroidal_field = equilibrium__vacuum_toroidal_field()
    var"code" :: equilibrium__code = equilibrium__code()
    _parent :: WeakRef = WeakRef(nothing)
    function equilibrium(var"time_slice"=FDSvector(equilibrium__time_slice[]), var"time"=nothing, var"ids_properties"=equilibrium__ids_properties(), var"grids_ggd"=FDSvector(equilibrium__grids_ggd[]), var"vacuum_toroidal_field"=equilibrium__vacuum_toroidal_field(), var"code"=equilibrium__code(), _parent=WeakRef(nothing))
        obj = new(var"time_slice", var"time", var"ids_properties", var"grids_ggd", var"vacuum_toroidal_field", var"code", _parent)
        obj.time_slice._parent = WeakRef(obj)
        obj.ids_properties._parent = WeakRef(obj)
        obj.grids_ggd._parent = WeakRef(obj)
        obj.vacuum_toroidal_field._parent = WeakRef(obj)
        obj.code._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct dataset_description__simulation <: FDS
    var"time_ended" :: Union{Nothing, String} = nothing
    var"time_begun" :: Union{Nothing, String} = nothing
    var"time_current" :: Union{Nothing, Float64} = nothing
    var"time_restart" :: Union{Nothing, Float64} = nothing
    var"workflow" :: Union{Nothing, String} = nothing
    var"comment_after" :: Union{Nothing, String} = nothing
    var"time_begin" :: Union{Nothing, Float64} = nothing
    var"time_end" :: Union{Nothing, Float64} = nothing
    var"comment_before" :: Union{Nothing, String} = nothing
    var"time_step" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__simulation(var"time_ended"=nothing, var"time_begun"=nothing, var"time_current"=nothing, var"time_restart"=nothing, var"workflow"=nothing, var"comment_after"=nothing, var"time_begin"=nothing, var"time_end"=nothing, var"comment_before"=nothing, var"time_step"=nothing, _parent=WeakRef(nothing))
        obj = new(var"time_ended", var"time_begun", var"time_current", var"time_restart", var"workflow", var"comment_after", var"time_begin", var"time_end", var"comment_before", var"time_step", _parent)

        return obj
    end
end

Base.@kwdef mutable struct dataset_description__pulse_time_end_epoch <: FDS
    var"nanoseconds" :: Union{Nothing, Int} = nothing
    var"seconds" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__pulse_time_end_epoch(var"nanoseconds"=nothing, var"seconds"=nothing, _parent=WeakRef(nothing))
        obj = new(var"nanoseconds", var"seconds", _parent)

        return obj
    end
end

Base.@kwdef mutable struct dataset_description__pulse_time_begin_epoch <: FDS
    var"nanoseconds" :: Union{Nothing, Int} = nothing
    var"seconds" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__pulse_time_begin_epoch(var"nanoseconds"=nothing, var"seconds"=nothing, _parent=WeakRef(nothing))
        obj = new(var"nanoseconds", var"seconds", _parent)

        return obj
    end
end

Base.@kwdef mutable struct dataset_description__parent_entry <: FDS
    var"pulse_type" :: Union{Nothing, String} = nothing
    var"run" :: Union{Nothing, Int} = nothing
    var"machine" :: Union{Nothing, String} = nothing
    var"pulse" :: Union{Nothing, Int} = nothing
    var"user" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__parent_entry(var"pulse_type"=nothing, var"run"=nothing, var"machine"=nothing, var"pulse"=nothing, var"user"=nothing, _parent=WeakRef(nothing))
        obj = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)

        return obj
    end
end

Base.@kwdef mutable struct dataset_description__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Nothing, String} = nothing
    var"data_dictionary" :: Union{Nothing, String} = nothing
    var"access_layer" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__ids_properties__version_put(var"access_layer_language"=nothing, var"data_dictionary"=nothing, var"access_layer"=nothing, _parent=WeakRef(nothing))
        obj = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)

        return obj
    end
end

Base.@kwdef mutable struct dataset_description__ids_properties <: FDS
    var"provider" :: Union{Nothing, String} = nothing
    var"version_put" :: dataset_description__ids_properties__version_put = dataset_description__ids_properties__version_put()
    var"homogeneous_time" :: Union{Nothing, Int} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"creation_date" :: Union{Nothing, String} = nothing
    var"comment" :: Union{Nothing, String} = nothing
    var"occurrence" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__ids_properties(var"provider"=nothing, var"version_put"=dataset_description__ids_properties__version_put(), var"homogeneous_time"=nothing, var"source"=nothing, var"creation_date"=nothing, var"comment"=nothing, var"occurrence"=nothing, _parent=WeakRef(nothing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        obj.version_put._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct dataset_description__data_entry <: FDS
    var"pulse_type" :: Union{Nothing, String} = nothing
    var"run" :: Union{Nothing, Int} = nothing
    var"machine" :: Union{Nothing, String} = nothing
    var"pulse" :: Union{Nothing, Int} = nothing
    var"user" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description__data_entry(var"pulse_type"=nothing, var"run"=nothing, var"machine"=nothing, var"pulse"=nothing, var"user"=nothing, _parent=WeakRef(nothing))
        obj = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)

        return obj
    end
end

Base.@kwdef mutable struct dataset_description <: FDS
    var"pulse_time_begin_epoch" :: dataset_description__pulse_time_begin_epoch = dataset_description__pulse_time_begin_epoch()
    var"imas_version" :: Union{Nothing, String} = nothing
    var"ids_properties" :: dataset_description__ids_properties = dataset_description__ids_properties()
    var"time" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"dd_version" :: Union{Nothing, String} = nothing
    var"parent_entry" :: dataset_description__parent_entry = dataset_description__parent_entry()
    var"simulation" :: dataset_description__simulation = dataset_description__simulation()
    var"pulse_time_end_epoch" :: dataset_description__pulse_time_end_epoch = dataset_description__pulse_time_end_epoch()
    var"data_entry" :: dataset_description__data_entry = dataset_description__data_entry()
    var"pulse_time_begin" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function dataset_description(var"pulse_time_begin_epoch"=dataset_description__pulse_time_begin_epoch(), var"imas_version"=nothing, var"ids_properties"=dataset_description__ids_properties(), var"time"=nothing, var"dd_version"=nothing, var"parent_entry"=dataset_description__parent_entry(), var"simulation"=dataset_description__simulation(), var"pulse_time_end_epoch"=dataset_description__pulse_time_end_epoch(), var"data_entry"=dataset_description__data_entry(), var"pulse_time_begin"=nothing, _parent=WeakRef(nothing))
        obj = new(var"pulse_time_begin_epoch", var"imas_version", var"ids_properties", var"time", var"dd_version", var"parent_entry", var"simulation", var"pulse_time_end_epoch", var"data_entry", var"pulse_time_begin", _parent)
        obj.pulse_time_begin_epoch._parent = WeakRef(obj)
        obj.ids_properties._parent = WeakRef(obj)
        obj.parent_entry._parent = WeakRef(obj)
        obj.simulation._parent = WeakRef(obj)
        obj.pulse_time_end_epoch._parent = WeakRef(obj)
        obj.data_entry._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"r0" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__vacuum_toroidal_field(var"b0"=nothing, var"r0"=nothing, _parent=WeakRef(nothing))
        obj = new(var"b0", var"r0", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method = core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___zeff_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method = core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___t_i_average_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___velocity <: FDS
    var"parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___neutral___velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=WeakRef(nothing))
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___velocity <: FDS
    var"parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___neutral___state___velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=WeakRef(nothing))
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___neutral_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___neutral___state___neutral_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"vibrational_level" :: Union{Nothing, Float64} = nothing
    var"temperature" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"electron_configuration" :: Union{Nothing, String} = nothing
    var"pressure" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"vibrational_mode" :: Union{Nothing, String} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"velocity" :: core_profiles__profiles_1d___neutral___state___velocity = core_profiles__profiles_1d___neutral___state___velocity()
    var"density" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"neutral_type" :: core_profiles__profiles_1d___neutral___state___neutral_type = core_profiles__profiles_1d___neutral___state___neutral_type()
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___neutral___state(var"label"=nothing, var"vibrational_level"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"pressure_fast_perpendicular"=nothing, var"electron_configuration"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"vibrational_mode"=nothing, var"pressure_fast_parallel"=nothing, var"velocity"=core_profiles__profiles_1d___neutral___state___velocity(), var"density"=nothing, var"density_fast"=nothing, var"neutral_type"=core_profiles__profiles_1d___neutral___state___neutral_type(), _parent=WeakRef(nothing))
        obj = new(var"label", var"vibrational_level", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"velocity", var"density", var"density_fast", var"neutral_type", _parent)
        obj.velocity._parent = WeakRef(obj)
        obj.neutral_type._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___element <: FDS
    var"atoms_n" :: Union{Nothing, Int} = nothing
    var"z_n" :: Union{Nothing, Float64} = nothing
    var"multiplicity" :: Union{Nothing, Float64} = nothing
    var"a" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___neutral___element(var"atoms_n"=nothing, var"z_n"=nothing, var"multiplicity"=nothing, var"a"=nothing, _parent=WeakRef(nothing))
        obj = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"temperature" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ion_index" :: Union{Nothing, Int} = nothing
    var"multiple_states_flag" :: Union{Nothing, Int} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral___state} = FDSvector(core_profiles__profiles_1d___neutral___state[])
    var"velocity" :: core_profiles__profiles_1d___neutral___velocity = core_profiles__profiles_1d___neutral___velocity()
    var"density" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral___element} = FDSvector(core_profiles__profiles_1d___neutral___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___neutral(var"label"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"ion_index"=nothing, var"multiple_states_flag"=nothing, var"pressure_fast_perpendicular"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"pressure_fast_parallel"=nothing, var"state"=FDSvector(core_profiles__profiles_1d___neutral___state[]), var"velocity"=core_profiles__profiles_1d___neutral___velocity(), var"density"=nothing, var"density_fast"=nothing, var"element"=FDSvector(core_profiles__profiles_1d___neutral___element[]), _parent=WeakRef(nothing))
        obj = new(var"label", var"temperature", var"pressure_thermal", var"ion_index", var"multiple_states_flag", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"pressure_fast_parallel", var"state", var"velocity", var"density", var"density_fast", var"element", _parent)
        obj.state._parent = WeakRef(obj)
        obj.velocity._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___velocity <: FDS
    var"parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=WeakRef(nothing))
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___temperature_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___velocity <: FDS
    var"parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___state___velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=WeakRef(nothing))
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___state___density_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"vibrational_level" :: Union{Nothing, Float64} = nothing
    var"rotation_frequency_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z_min" :: Union{Nothing, Float64} = nothing
    var"electron_configuration" :: Union{Nothing, String} = nothing
    var"pressure" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"vibrational_mode" :: Union{Nothing, String} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z_average_square_1d" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"velocity" :: core_profiles__profiles_1d___ion___state___velocity = core_profiles__profiles_1d___ion___state___velocity()
    var"z_average" :: Union{Nothing, Float64} = nothing
    var"z_max" :: Union{Nothing, Float64} = nothing
    var"z_square_average" :: Union{Nothing, Float64} = nothing
    var"ionisation_potential" :: Union{Nothing, Float64} = nothing
    var"density" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z_average_1d" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fit" :: core_profiles__profiles_1d___ion___state___density_fit = core_profiles__profiles_1d___ion___state___density_fit()
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___state(var"label"=nothing, var"vibrational_level"=nothing, var"rotation_frequency_tor"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"pressure_fast_perpendicular"=nothing, var"z_min"=nothing, var"electron_configuration"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"vibrational_mode"=nothing, var"pressure_fast_parallel"=nothing, var"z_average_square_1d"=nothing, var"velocity"=core_profiles__profiles_1d___ion___state___velocity(), var"z_average"=nothing, var"z_max"=nothing, var"z_square_average"=nothing, var"ionisation_potential"=nothing, var"density"=nothing, var"density_fast"=nothing, var"z_average_1d"=nothing, var"density_fit"=core_profiles__profiles_1d___ion___state___density_fit(), _parent=WeakRef(nothing))
        obj = new(var"label", var"vibrational_level", var"rotation_frequency_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"z_min", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"z_average_square_1d", var"velocity", var"z_average", var"z_max", var"z_square_average", var"ionisation_potential", var"density", var"density_fast", var"z_average_1d", var"density_fit", _parent)
        obj.velocity._parent = WeakRef(obj)
        obj.density_fit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___element <: FDS
    var"atoms_n" :: Union{Nothing, Int} = nothing
    var"z_n" :: Union{Nothing, Float64} = nothing
    var"multiplicity" :: Union{Nothing, Float64} = nothing
    var"a" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___element(var"atoms_n"=nothing, var"z_n"=nothing, var"multiplicity"=nothing, var"a"=nothing, _parent=WeakRef(nothing))
        obj = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method = core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion___density_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"rotation_frequency_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"temperature_validity" :: Union{Nothing, Int} = nothing
    var"velocity_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z_ion_1d" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"multiple_states_flag" :: Union{Nothing, Int} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"neutral_index" :: Union{Nothing, Int} = nothing
    var"pressure" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_validity" :: Union{Nothing, Int} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion___state} = FDSvector(core_profiles__profiles_1d___ion___state[])
    var"velocity" :: core_profiles__profiles_1d___ion___velocity = core_profiles__profiles_1d___ion___velocity()
    var"z_ion" :: Union{Nothing, Float64} = nothing
    var"temperature_fit" :: core_profiles__profiles_1d___ion___temperature_fit = core_profiles__profiles_1d___ion___temperature_fit()
    var"density" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"velocity_pol" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fit" :: core_profiles__profiles_1d___ion___density_fit = core_profiles__profiles_1d___ion___density_fit()
    var"z_ion_square_1d" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion___element} = FDSvector(core_profiles__profiles_1d___ion___element[])
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___ion(var"label"=nothing, var"rotation_frequency_tor"=nothing, var"temperature_validity"=nothing, var"velocity_tor"=nothing, var"temperature"=nothing, var"z_ion_1d"=nothing, var"pressure_thermal"=nothing, var"multiple_states_flag"=nothing, var"pressure_fast_perpendicular"=nothing, var"neutral_index"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"density_validity"=nothing, var"pressure_fast_parallel"=nothing, var"state"=FDSvector(core_profiles__profiles_1d___ion___state[]), var"velocity"=core_profiles__profiles_1d___ion___velocity(), var"z_ion"=nothing, var"temperature_fit"=core_profiles__profiles_1d___ion___temperature_fit(), var"density"=nothing, var"velocity_pol"=nothing, var"density_fast"=nothing, var"density_fit"=core_profiles__profiles_1d___ion___density_fit(), var"z_ion_square_1d"=nothing, var"element"=FDSvector(core_profiles__profiles_1d___ion___element[]), _parent=WeakRef(nothing))
        obj = new(var"label", var"rotation_frequency_tor", var"temperature_validity", var"velocity_tor", var"temperature", var"z_ion_1d", var"pressure_thermal", var"multiple_states_flag", var"pressure_fast_perpendicular", var"neutral_index", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"state", var"velocity", var"z_ion", var"temperature_fit", var"density", var"velocity_pol", var"density_fast", var"density_fit", var"z_ion_square_1d", var"element", _parent)
        obj.state._parent = WeakRef(obj)
        obj.velocity._parent = WeakRef(obj)
        obj.temperature_fit._parent = WeakRef(obj)
        obj.density_fit._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___grid <: FDS
    var"psi" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"psi_boundary" :: Union{Nothing, Float64} = nothing
    var"volume" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"area" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_pol_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"surface" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"psi_magnetic_axis" :: Union{Nothing, Float64} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___grid(var"psi"=nothing, var"psi_boundary"=nothing, var"volume"=nothing, var"area"=nothing, var"rho_pol_norm"=nothing, var"rho_tor_norm"=nothing, var"surface"=nothing, var"rho_tor"=nothing, var"psi_magnetic_axis"=nothing, _parent=WeakRef(nothing))
        obj = new(var"psi", var"psi_boundary", var"volume", var"area", var"rho_pol_norm", var"rho_tor_norm", var"surface", var"rho_tor", var"psi_magnetic_axis", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__velocity <: FDS
    var"parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___electrons__velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=WeakRef(nothing))
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___electrons__temperature_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit <: FDS
    var"local" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"source" :: Union{Nothing, AbstractArray{String, 1}} = nothing
    var"measured" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method = core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___electrons__density_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=WeakRef(nothing))
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons <: FDS
    var"temperature_validity" :: Union{Nothing, Int} = nothing
    var"velocity_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_validity" :: Union{Nothing, Int} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"velocity" :: core_profiles__profiles_1d___electrons__velocity = core_profiles__profiles_1d___electrons__velocity()
    var"temperature_fit" :: core_profiles__profiles_1d___electrons__temperature_fit = core_profiles__profiles_1d___electrons__temperature_fit()
    var"density" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"velocity_pol" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"collisionality_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"density_fit" :: core_profiles__profiles_1d___electrons__density_fit = core_profiles__profiles_1d___electrons__density_fit()
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___electrons(var"temperature_validity"=nothing, var"velocity_tor"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"pressure_fast_perpendicular"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"density_validity"=nothing, var"pressure_fast_parallel"=nothing, var"velocity"=core_profiles__profiles_1d___electrons__velocity(), var"temperature_fit"=core_profiles__profiles_1d___electrons__temperature_fit(), var"density"=nothing, var"velocity_pol"=nothing, var"collisionality_norm"=nothing, var"density_fast"=nothing, var"density_fit"=core_profiles__profiles_1d___electrons__density_fit(), _parent=WeakRef(nothing))
        obj = new(var"temperature_validity", var"velocity_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"velocity", var"temperature_fit", var"density", var"velocity_pol", var"collisionality_norm", var"density_fast", var"density_fit", _parent)
        obj.velocity._parent = WeakRef(obj)
        obj.temperature_fit._parent = WeakRef(obj)
        obj.density_fit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___e_field <: FDS
    var"parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d___e_field(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=WeakRef(nothing))
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d <: FDS
    var"pressure_ion_total" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"time" :: Union{Nothing, Float64} = nothing
    var"t_i_average_fit" :: core_profiles__profiles_1d___t_i_average_fit = core_profiles__profiles_1d___t_i_average_fit()
    var"neutral" :: FDSvector{T} where {T<:core_profiles__profiles_1d___neutral} = FDSvector(core_profiles__profiles_1d___neutral[])
    var"n_i_thermal_total" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"magnetic_shear" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ion" :: FDSvector{T} where {T<:core_profiles__profiles_1d___ion} = FDSvector(core_profiles__profiles_1d___ion[])
    var"j_total" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"rotation_frequency_tor_sonic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"j_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"current_parallel_inside" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"j_non_inductive" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"e_field_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"momentum_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"conductivity_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"electrons" :: core_profiles__profiles_1d___electrons = core_profiles__profiles_1d___electrons()
    var"pressure_perpendicular" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"q" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"t_i_average" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"j_ohmic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"grid" :: core_profiles__profiles_1d___grid = core_profiles__profiles_1d___grid()
    var"phi_potential" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"j_bootstrap" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"zeff_fit" :: core_profiles__profiles_1d___zeff_fit = core_profiles__profiles_1d___zeff_fit()
    var"pressure_parallel" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"e_field" :: core_profiles__profiles_1d___e_field = core_profiles__profiles_1d___e_field()
    var"zeff" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"n_i_total_over_n_e" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__profiles_1d(var"pressure_ion_total"=nothing, var"time"=nothing, var"t_i_average_fit"=core_profiles__profiles_1d___t_i_average_fit(), var"neutral"=FDSvector(core_profiles__profiles_1d___neutral[]), var"n_i_thermal_total"=nothing, var"magnetic_shear"=nothing, var"ion"=FDSvector(core_profiles__profiles_1d___ion[]), var"j_total"=nothing, var"rotation_frequency_tor_sonic"=nothing, var"pressure_thermal"=nothing, var"j_tor"=nothing, var"current_parallel_inside"=nothing, var"j_non_inductive"=nothing, var"e_field_parallel"=nothing, var"momentum_tor"=nothing, var"conductivity_parallel"=nothing, var"electrons"=core_profiles__profiles_1d___electrons(), var"pressure_perpendicular"=nothing, var"q"=nothing, var"t_i_average"=nothing, var"j_ohmic"=nothing, var"grid"=core_profiles__profiles_1d___grid(), var"phi_potential"=nothing, var"j_bootstrap"=nothing, var"zeff_fit"=core_profiles__profiles_1d___zeff_fit(), var"pressure_parallel"=nothing, var"e_field"=core_profiles__profiles_1d___e_field(), var"zeff"=nothing, var"n_i_total_over_n_e"=nothing, _parent=WeakRef(nothing))
        obj = new(var"pressure_ion_total", var"time", var"t_i_average_fit", var"neutral", var"n_i_thermal_total", var"magnetic_shear", var"ion", var"j_total", var"rotation_frequency_tor_sonic", var"pressure_thermal", var"j_tor", var"current_parallel_inside", var"j_non_inductive", var"e_field_parallel", var"momentum_tor", var"conductivity_parallel", var"electrons", var"pressure_perpendicular", var"q", var"t_i_average", var"j_ohmic", var"grid", var"phi_potential", var"j_bootstrap", var"zeff_fit", var"pressure_parallel", var"e_field", var"zeff", var"n_i_total_over_n_e", _parent)
        obj.t_i_average_fit._parent = WeakRef(obj)
        obj.neutral._parent = WeakRef(obj)
        obj.ion._parent = WeakRef(obj)
        obj.electrons._parent = WeakRef(obj)
        obj.grid._parent = WeakRef(obj)
        obj.zeff_fit._parent = WeakRef(obj)
        obj.e_field._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Nothing, String} = nothing
    var"data_dictionary" :: Union{Nothing, String} = nothing
    var"access_layer" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__ids_properties__version_put(var"access_layer_language"=nothing, var"data_dictionary"=nothing, var"access_layer"=nothing, _parent=WeakRef(nothing))
        obj = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__ids_properties <: FDS
    var"provider" :: Union{Nothing, String} = nothing
    var"version_put" :: core_profiles__ids_properties__version_put = core_profiles__ids_properties__version_put()
    var"homogeneous_time" :: Union{Nothing, Int} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"creation_date" :: Union{Nothing, String} = nothing
    var"comment" :: Union{Nothing, String} = nothing
    var"occurrence" :: Union{Nothing, Int} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__ids_properties(var"provider"=nothing, var"version_put"=core_profiles__ids_properties__version_put(), var"homogeneous_time"=nothing, var"source"=nothing, var"creation_date"=nothing, var"comment"=nothing, var"occurrence"=nothing, _parent=WeakRef(nothing))
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        obj.version_put._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__global_quantities <: FDS
    var"beta_tor_norm" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"resistive_psi_losses" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ip" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"li_3" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"t_i_average_peaking" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"t_e_peaking" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"beta_tor" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"z_eff_resistive" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ejima" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"energy_diamagnetic" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"li" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"current_non_inductive" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"v_loop" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"beta_pol" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"current_bootstrap" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__global_quantities(var"beta_tor_norm"=nothing, var"resistive_psi_losses"=nothing, var"ip"=nothing, var"li_3"=nothing, var"t_i_average_peaking"=nothing, var"t_e_peaking"=nothing, var"beta_tor"=nothing, var"z_eff_resistive"=nothing, var"ejima"=nothing, var"energy_diamagnetic"=nothing, var"li"=nothing, var"current_non_inductive"=nothing, var"v_loop"=nothing, var"beta_pol"=nothing, var"current_bootstrap"=nothing, _parent=WeakRef(nothing))
        obj = new(var"beta_tor_norm", var"resistive_psi_losses", var"ip", var"li_3", var"t_i_average_peaking", var"t_e_peaking", var"beta_tor", var"z_eff_resistive", var"ejima", var"energy_diamagnetic", var"li", var"current_non_inductive", var"v_loop", var"beta_pol", var"current_bootstrap", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__code__library <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"commit" :: Union{Nothing, String} = nothing
    var"repository" :: Union{Nothing, String} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__code__library(var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"version"=nothing, _parent=WeakRef(nothing))
        obj = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__code <: FDS
    var"library" :: FDSvector{T} where {T<:core_profiles__code__library} = FDSvector(core_profiles__code__library[])
    var"name" :: Union{Nothing, String} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"commit" :: Union{Nothing, String} = nothing
    var"repository" :: Union{Nothing, String} = nothing
    var"output_flag" :: Union{Nothing, AbstractArray{Int, 1}} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles__code(var"library"=FDSvector(core_profiles__code__library[]), var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"output_flag"=nothing, var"version"=nothing, _parent=WeakRef(nothing))
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        obj.library._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles <: FDS
    var"time" :: Union{Nothing, AbstractArray{Float64, 1}} = nothing
    var"ids_properties" :: core_profiles__ids_properties = core_profiles__ids_properties()
    var"vacuum_toroidal_field" :: core_profiles__vacuum_toroidal_field = core_profiles__vacuum_toroidal_field()
    var"code" :: core_profiles__code = core_profiles__code()
    var"global_quantities" :: core_profiles__global_quantities = core_profiles__global_quantities()
    var"profiles_1d" :: FDSvector{T} where {T<:core_profiles__profiles_1d} = FDSvector(core_profiles__profiles_1d[])
    _parent :: WeakRef = WeakRef(nothing)
    function core_profiles(var"time"=nothing, var"ids_properties"=core_profiles__ids_properties(), var"vacuum_toroidal_field"=core_profiles__vacuum_toroidal_field(), var"code"=core_profiles__code(), var"global_quantities"=core_profiles__global_quantities(), var"profiles_1d"=FDSvector(core_profiles__profiles_1d[]), _parent=WeakRef(nothing))
        obj = new(var"time", var"ids_properties", var"vacuum_toroidal_field", var"code", var"global_quantities", var"profiles_1d", _parent)
        obj.ids_properties._parent = WeakRef(obj)
        obj.vacuum_toroidal_field._parent = WeakRef(obj)
        obj.code._parent = WeakRef(obj)
        obj.global_quantities._parent = WeakRef(obj)
        obj.profiles_1d._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct dd <: FDS
    var"equilibrium" :: Union{Nothing, equilibrium} = equilibrium()
    var"core_profiles" :: Union{Nothing, core_profiles} = core_profiles()
    var"wall" :: Union{Nothing, wall} = wall()
    var"dataset_description" :: Union{Nothing, dataset_description} = dataset_description()
    _parent :: WeakRef = WeakRef(nothing)
    function dd(var"equilibrium"=equilibrium(), var"core_profiles"=core_profiles(), var"wall"=wall(), var"dataset_description"=dataset_description(), _parent=WeakRef(nothing))
        obj = new(var"equilibrium", var"core_profiles", var"wall", var"dataset_description", _parent)
        obj.equilibrium._parent = WeakRef(obj)
        obj.core_profiles._parent = WeakRef(obj)
        obj.wall._parent = WeakRef(obj)
        obj.dataset_description._parent = WeakRef(obj)
        return obj
    end
end

