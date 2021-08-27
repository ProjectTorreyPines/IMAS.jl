include("functionarrays.jl")

conversion_types = Union{Float64, Int64, String, AbstractFDVector{Float64}, AbstractFDVector{Int64}, AbstractFDVector{String}, Array{Float64, N} where N}

Base.@kwdef mutable struct wall__temperature_reference <: FDS
    var"data" :: Union{Missing, Float64} = missing
    var"description" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__temperature_reference(var"data"=missing, var"description"=missing, _parent=WeakRef(missing))
        fds = new(var"data", var"description", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral___element <: FDSvectorElement
    var"atoms_n" :: Union{Missing, Int64} = missing
    var"z_n" :: Union{Missing, Float64} = missing
    var"multiplicity" :: Union{Missing, Float64} = missing
    var"a" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__global_quantities__neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        fds = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__global_quantities__neutral <: FDSvectorElement
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
        fds = new(var"label", var"sputtering_chemical_coefficient", var"gas_puff", var"recycling_particles_coefficient", var"pumping_speed", var"particle_flux_from_wall", var"recycling_energy_coefficient", var"wall_inventory", var"particle_flux_from_plasma", var"sputtering_physical_coefficient", var"element", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
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
    function wall__global_quantities__electrons(var"particle_flux_from_plasma"=missing, var"gas_puff"=missing, var"power_outer_target"=missing, var"pumping_speed"=missing, var"particle_flux_from_wall"=missing, var"power_inner_target"=missing, _parent=WeakRef(missing))
        fds = new(var"particle_flux_from_plasma", var"gas_puff", var"power_outer_target", var"pumping_speed", var"particle_flux_from_wall", var"power_inner_target", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"neutral", var"power_incident", var"power_radiated", var"power_inner_target_ion_total", var"temperature", var"power_conducted", var"power_convected", var"current_tor", var"electrons", var"power_density_inner_target_max", var"power_black_body", var"power_recombination_neutrals", var"power_to_cooling", var"power_density_outer_target_max", var"power_recombination_plasma", var"power_currents", var"power_neutrals", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.neutral, :_parent, WeakRef(fds))
        setfield!(fds.electrons, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__first_wall_power_flux_peak <: FDS
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"data" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__first_wall_power_flux_peak(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        fds = new(var"time", var"data", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary <: FDSvectorElement
    var"neighbours" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"neighbours", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension___object <: FDSvectorElement
    var"nodes" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"measure" :: Union{Missing, Float64} = missing
    var"geometry" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"boundary" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        fds = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___objects_per_dimension <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension___object} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___objects_per_dimension(var"object"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space___geometry_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___space <: FDSvectorElement
    var"coordinates_type" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"geometry_type" :: wall__description_ggd___grid_ggd___space___geometry_type = wall__description_ggd___grid_ggd___space___geometry_type()
    var"identifier" :: wall__description_ggd___grid_ggd___space___identifier = wall__description_ggd___grid_ggd___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space___objects_per_dimension} = FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___space(var"coordinates_type"=missing, var"geometry_type"=wall__description_ggd___grid_ggd___space___geometry_type(), var"identifier"=wall__description_ggd___grid_ggd___space___identifier(), var"objects_per_dimension"=FDSvector(wall__description_ggd___grid_ggd___space___objects_per_dimension[]), _parent=WeakRef(missing))
        fds = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.geometry_type, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.objects_per_dimension, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element___object <: FDSvectorElement
    var"dimension" :: Union{Missing, Int64} = missing
    var"space" :: Union{Missing, Int64} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"dimension", var"space", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___element <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element___object} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___element(var"object"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset___base <: FDSvectorElement
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd___grid_subset <: FDSvectorElement
    var"base" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___base} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[])
    var"metric" :: wall__description_ggd___grid_ggd___grid_subset___metric = wall__description_ggd___grid_ggd___grid_subset___metric()
    var"dimension" :: Union{Missing, Int64} = missing
    var"identifier" :: wall__description_ggd___grid_ggd___grid_subset___identifier = wall__description_ggd___grid_ggd___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset___element} = FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd___grid_subset(var"base"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___base[]), var"metric"=wall__description_ggd___grid_ggd___grid_subset___metric(), var"dimension"=missing, var"identifier"=wall__description_ggd___grid_ggd___grid_subset___identifier(), var"element"=FDSvector(wall__description_ggd___grid_ggd___grid_subset___element[]), _parent=WeakRef(missing))
        fds = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.base, :_parent, WeakRef(fds))
        setfield!(fds.metric, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___grid_ggd <: FDSvectorElement
    var"time" :: Union{Missing, Float64} = missing
    var"grid_subset" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___grid_subset} = FDSvector(wall__description_ggd___grid_ggd___grid_subset[])
    var"space" :: FDSvector{T} where {T<:wall__description_ggd___grid_ggd___space} = FDSvector(wall__description_ggd___grid_ggd___space[])
    var"identifier" :: wall__description_ggd___grid_ggd___identifier = wall__description_ggd___grid_ggd___identifier()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___grid_ggd(var"time"=missing, var"grid_subset"=FDSvector(wall__description_ggd___grid_ggd___grid_subset[]), var"space"=FDSvector(wall__description_ggd___grid_ggd___space[]), var"identifier"=wall__description_ggd___grid_ggd___identifier(), _parent=WeakRef(missing))
        fds = new(var"time", var"grid_subset", var"space", var"identifier", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.grid_subset, :_parent, WeakRef(fds))
        setfield!(fds.space, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___temperature <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd___temperature(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd___power_density <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd___power_density(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_ggd___ggd <: FDSvectorElement
    var"temperature" :: FDSvector{T} where {T<:wall__description_ggd___ggd___temperature} = FDSvector(wall__description_ggd___ggd___temperature[])
    var"time" :: Union{Missing, Float64} = missing
    var"power_density" :: FDSvector{T} where {T<:wall__description_ggd___ggd___power_density} = FDSvector(wall__description_ggd___ggd___power_density[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_ggd___ggd(var"temperature"=FDSvector(wall__description_ggd___ggd___temperature[]), var"time"=missing, var"power_density"=FDSvector(wall__description_ggd___ggd___power_density[]), _parent=WeakRef(missing))
        fds = new(var"temperature", var"time", var"power_density", _parent)
        assign_derived_quantities(fds)
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
        assign_derived_quantities(fds)
        setfield!(fds.grid_ggd, :_parent, WeakRef(fds))
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.ggd, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___outline <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element___outline(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element___j_tor <: FDS
    var"time" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"data" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element___j_tor(var"time"=missing, var"data"=missing, _parent=WeakRef(missing))
        fds = new(var"time", var"data", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___element <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"j_tor" :: wall__description_2d___vessel__unit___element___j_tor = wall__description_2d___vessel__unit___element___j_tor()
    var"resistivity" :: Union{Missing, Float64} = missing
    var"outline" :: wall__description_2d___vessel__unit___element___outline = wall__description_2d___vessel__unit___element___outline()
    var"resistance" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___element(var"name"=missing, var"j_tor"=wall__description_2d___vessel__unit___element___j_tor(), var"resistivity"=missing, var"outline"=wall__description_2d___vessel__unit___element___outline(), var"resistance"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"j_tor", var"resistivity", var"outline", var"resistance", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.j_tor, :_parent, WeakRef(fds))
        setfield!(fds.outline, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_outer <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular__outline_outer(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__outline_inner <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular__outline_inner(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular__centreline <: FDS
    var"closed" :: Union{Missing, Int64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular__centreline(var"closed"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"closed", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit___annular <: FDS
    var"outline_inner" :: wall__description_2d___vessel__unit___annular__outline_inner = wall__description_2d___vessel__unit___annular__outline_inner()
    var"centreline" :: wall__description_2d___vessel__unit___annular__centreline = wall__description_2d___vessel__unit___annular__centreline()
    var"thickness" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"resistivity" :: Union{Missing, Float64} = missing
    var"outline_outer" :: wall__description_2d___vessel__unit___annular__outline_outer = wall__description_2d___vessel__unit___annular__outline_outer()
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit___annular(var"outline_inner"=wall__description_2d___vessel__unit___annular__outline_inner(), var"centreline"=wall__description_2d___vessel__unit___annular__centreline(), var"thickness"=missing, var"resistivity"=missing, var"outline_outer"=wall__description_2d___vessel__unit___annular__outline_outer(), _parent=WeakRef(missing))
        fds = new(var"outline_inner", var"centreline", var"thickness", var"resistivity", var"outline_outer", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.outline_inner, :_parent, WeakRef(fds))
        setfield!(fds.centreline, :_parent, WeakRef(fds))
        setfield!(fds.outline_outer, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__unit <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"annular" :: wall__description_2d___vessel__unit___annular = wall__description_2d___vessel__unit___annular()
    var"identifier" :: Union{Missing, String} = missing
    var"element" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit___element} = FDSvector(wall__description_2d___vessel__unit___element[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__unit(var"name"=missing, var"annular"=wall__description_2d___vessel__unit___annular(), var"identifier"=missing, var"element"=FDSvector(wall__description_2d___vessel__unit___element[]), _parent=WeakRef(missing))
        fds = new(var"name", var"annular", var"identifier", var"element", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.annular, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel__type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___vessel <: FDS
    var"type" :: wall__description_2d___vessel__type = wall__description_2d___vessel__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___vessel__unit} = FDSvector(wall__description_2d___vessel__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___vessel(var"type"=wall__description_2d___vessel__type(), var"unit"=FDSvector(wall__description_2d___vessel__unit[]), _parent=WeakRef(missing))
        fds = new(var"type", var"unit", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.unit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit___outline <: FDSvectorElement
    var"time" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__unit___outline(var"time"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"time", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__unit <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"resistivity" :: Union{Missing, Float64} = missing
    var"closed" :: Union{Missing, Int64} = missing
    var"phi_extensions" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"outline" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit___outline} = FDSvector(wall__description_2d___mobile__unit___outline[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"phi_extensions"=missing, var"outline"=FDSvector(wall__description_2d___mobile__unit___outline[]), _parent=WeakRef(missing))
        fds = new(var"name", var"resistivity", var"closed", var"phi_extensions", var"outline", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.outline, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile__type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___mobile <: FDS
    var"type" :: wall__description_2d___mobile__type = wall__description_2d___mobile__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___mobile__unit} = FDSvector(wall__description_2d___mobile__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___mobile(var"type"=wall__description_2d___mobile__type(), var"unit"=FDSvector(wall__description_2d___mobile__unit[]), _parent=WeakRef(missing))
        fds = new(var"type", var"unit", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.unit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit___outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__unit___outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__unit <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"resistivity" :: Union{Missing, Float64} = missing
    var"closed" :: Union{Missing, Int64} = missing
    var"outline" :: wall__description_2d___limiter__unit___outline = wall__description_2d___limiter__unit___outline()
    var"phi_extensions" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__unit(var"name"=missing, var"resistivity"=missing, var"closed"=missing, var"outline"=wall__description_2d___limiter__unit___outline(), var"phi_extensions"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"resistivity", var"closed", var"outline", var"phi_extensions", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.outline, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter__type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter__type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct wall__description_2d___limiter <: FDS
    var"type" :: wall__description_2d___limiter__type = wall__description_2d___limiter__type()
    var"unit" :: FDSvector{T} where {T<:wall__description_2d___limiter__unit} = FDSvector(wall__description_2d___limiter__unit[])
    _parent :: WeakRef = WeakRef(missing)
    function wall__description_2d___limiter(var"type"=wall__description_2d___limiter__type(), var"unit"=FDSvector(wall__description_2d___limiter__unit[]), _parent=WeakRef(missing))
        fds = new(var"type", var"unit", _parent)
        assign_derived_quantities(fds)
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
        assign_derived_quantities(fds)
        setfield!(fds.mobile, :_parent, WeakRef(fds))
        setfield!(fds.limiter, :_parent, WeakRef(fds))
        setfield!(fds.type, :_parent, WeakRef(fds))
        setfield!(fds.vessel, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct wall__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function wall__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
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
        fds = new(var"time", var"ids_properties", var"description_ggd", var"description_2d", var"first_wall_surface_area", var"code", var"global_quantities", var"temperature_reference", var"first_wall_power_flux_peak", _parent)
        assign_derived_quantities(fds)
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

Base.@kwdef mutable struct equilibrium__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"r0" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        fds = new(var"b0", var"r0", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_2d___grid_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d___grid <: FDS
    var"volume_element" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"dim2" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dim1" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_2d___grid(var"volume_element"=missing, var"dim2"=missing, var"dim1"=missing, _parent=WeakRef(missing))
        fds = new(var"volume_element", var"dim2", var"dim1", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_2d <: FDSvectorElement
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
        fds = new(var"psi", var"b_field_r", var"r", var"b_r", var"theta", var"b_field_z", var"j_tor", var"phi", var"z", var"b_field_tor", var"b_z", var"grid", var"grid_type", var"j_parallel", var"b_tor", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.grid, :_parent, WeakRef(fds))
        setfield!(fds.grid_type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___profiles_1d__geometric_axis <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___profiles_1d__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"b_field_max", var"dvolume_drho_tor", var"gm9", var"dpsi_drho_tor", var"surface", var"rho_tor", var"magnetic_shear", var"b_average", var"b_field_min", var"darea_dpsi", var"gm3", var"squareness_upper_inner", var"squareness_lower_inner", var"rho_tor_norm", var"elongation", var"beta_pol", var"b_field_average", var"j_parallel", var"gm6", var"psi", var"gm8", var"dpressure_dpsi", var"triangularity_upper", var"darea_drho_tor", var"area", var"trapped_fraction", var"volume", var"dvolume_dpsi", var"b_min", var"f", var"mass_density", var"r_outboard", var"gm4", var"phi", var"squareness_lower_outer", var"triangularity_lower", var"gm2", var"rho_volume_norm", var"gm1", var"gm5", var"b_max", var"f_df_dpsi", var"j_tor", var"r_inboard", var"q", var"gm7", var"pressure", var"squareness_upper_outer", var"geometric_axis", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.geometric_axis, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__q_min <: FDS
    var"value" :: Union{Missing, Float64} = missing
    var"rho_tor_norm" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities__q_min(var"value"=missing, var"rho_tor_norm"=missing, _parent=WeakRef(missing))
        fds = new(var"value", var"rho_tor_norm", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__magnetic_axis <: FDS
    var"b_field_tor" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    var"b_tor" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities__magnetic_axis(var"b_field_tor"=missing, var"r"=missing, var"z"=missing, var"b_tor"=missing, _parent=WeakRef(missing))
        fds = new(var"b_field_tor", var"r", var"z", var"b_tor", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___global_quantities__current_centre <: FDS
    var"velocity_z" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___global_quantities__current_centre(var"velocity_z"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"velocity_z", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"ip", var"li_3", var"beta_tor", var"surface", var"magnetic_axis", var"energy_mhd", var"psi_boundary", var"length_pol", var"area", var"psi_external_average", var"q_95", var"q_axis", var"psi_axis", var"w_mhd", var"volume", var"plasma_inductance", var"beta_pol", var"beta_normal", var"current_centre", var"q_min", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.magnetic_axis, :_parent, WeakRef(fds))
        setfield!(fds.current_centre, :_parent, WeakRef(fds))
        setfield!(fds.q_min, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___z <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___z(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___theta <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___theta(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___r <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___r(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___psi <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___psi(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___phi <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___phi(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_tor <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___j_tor(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___j_parallel <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___j_parallel(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary <: FDSvectorElement
    var"neighbours" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"neighbours", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object <: FDSvectorElement
    var"nodes" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"measure" :: Union{Missing, Float64} = missing
    var"geometry" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        fds = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___objects_per_dimension <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___objects_per_dimension(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space___geometry_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__space <: FDSvectorElement
    var"coordinates_type" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"geometry_type" :: equilibrium__time_slice___ggd___grid__space___geometry_type = equilibrium__time_slice___ggd___grid__space___geometry_type()
    var"identifier" :: equilibrium__time_slice___ggd___grid__space___identifier = equilibrium__time_slice___ggd___grid__space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__space___objects_per_dimension} = FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__time_slice___ggd___grid__space___geometry_type(), var"identifier"=equilibrium__time_slice___ggd___grid__space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__time_slice___ggd___grid__space___objects_per_dimension[]), _parent=WeakRef(missing))
        fds = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.geometry_type, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.objects_per_dimension, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element___object <: FDSvectorElement
    var"dimension" :: Union{Missing, Int64} = missing
    var"space" :: Union{Missing, Int64} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"dimension", var"space", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___element <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element___object} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___element(var"object"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset___base <: FDSvectorElement
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___grid__grid_subset <: FDSvectorElement
    var"base" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___base} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[])
    var"metric" :: equilibrium__time_slice___ggd___grid__grid_subset___metric = equilibrium__time_slice___ggd___grid__grid_subset___metric()
    var"dimension" :: Union{Missing, Int64} = missing
    var"identifier" :: equilibrium__time_slice___ggd___grid__grid_subset___identifier = equilibrium__time_slice___ggd___grid__grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__time_slice___ggd___grid__grid_subset___element} = FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___grid__grid_subset(var"base"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___base[]), var"metric"=equilibrium__time_slice___ggd___grid__grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__time_slice___ggd___grid__grid_subset___identifier(), var"element"=FDSvector(equilibrium__time_slice___ggd___grid__grid_subset___element[]), _parent=WeakRef(missing))
        fds = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_derived_quantities(fds)
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
        assign_derived_quantities(fds)
        setfield!(fds.grid_subset, :_parent, WeakRef(fds))
        setfield!(fds.space, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_z <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___b_field_z(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_tor <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___b_field_tor(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___ggd___b_field_r <: FDSvectorElement
    var"grid_index" :: Union{Missing, Int64} = missing
    var"values" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"grid_subset_index" :: Union{Missing, Int64} = missing
    var"coefficients" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___ggd___b_field_r(var"grid_index"=missing, var"values"=missing, var"grid_subset_index"=missing, var"coefficients"=missing, _parent=WeakRef(missing))
        fds = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)
        assign_derived_quantities(fds)
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
        assign_derived_quantities(fds)
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
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___coordinate_system__grid_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___coordinate_system__grid <: FDS
    var"volume_element" :: Union{Missing, AbstractArray{Float64, 2}} = missing
    var"dim2" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"dim1" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___coordinate_system__grid(var"volume_element"=missing, var"dim2"=missing, var"dim1"=missing, _parent=WeakRef(missing))
        fds = new(var"volume_element", var"dim2", var"dim1", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"jacobian", var"g13_covariant", var"g11_contravariant", var"g13_contravariant", var"r", var"g12_contravariant", var"g22_contravariant", var"z", var"g33_contravariant", var"g22_covariant", var"tensor_contravariant", var"g12_covariant", var"g33_covariant", var"grid", var"g23_covariant", var"g11_covariant", var"tensor_covariant", var"g23_contravariant", var"grid_type", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.grid, :_parent, WeakRef(fds))
        setfield!(fds.grid_type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___convergence <: FDS
    var"iterations_n" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___convergence(var"iterations_n"=missing, _parent=WeakRef(missing))
        fds = new(var"iterations_n", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_reconstructed <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__x_point___position_reconstructed(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point___position_measured <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__x_point___position_measured(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__x_point <: FDSvectorElement
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
        fds = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.position_measured, :_parent, WeakRef(fds))
        setfield!(fds.position_reconstructed, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_reconstructed <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__strike_point___position_reconstructed(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point___position_measured <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__strike_point___position_measured(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__strike_point <: FDSvectorElement
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
        fds = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.position_measured, :_parent, WeakRef(fds))
        setfield!(fds.position_reconstructed, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q___position <: FDS
    var"phi" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__q___position(var"phi"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"phi", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__q <: FDSvectorElement
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
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"position", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.position, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pressure <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__pressure(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_passive_current <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__pf_passive_current(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__pf_current <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__pf_current(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e_line <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__n_e_line(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__n_e <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__n_e(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__mse_polarisation_angle <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__mse_polarisation_angle(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
    function equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__iron_core_segment <: FDSvectorElement
    var"magnetisation_r" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r()
    var"magnetisation_z" :: equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z = equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z()
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__iron_core_segment(var"magnetisation_r"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_r(), var"magnetisation_z"=equilibrium__time_slice___constraints__iron_core_segment___magnetisation_z(), _parent=WeakRef(missing))
        fds = new(var"magnetisation_r", var"magnetisation_z", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.magnetisation_r, :_parent, WeakRef(fds))
        setfield!(fds.magnetisation_z, :_parent, WeakRef(fds))
        return fds
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
    function equilibrium__time_slice___constraints__ip(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__flux_loop <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__flux_loop(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__faraday_angle <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__faraday_angle(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
    function equilibrium__time_slice___constraints__diamagnetic_flux(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___constraints__bpol_probe <: FDSvectorElement
    var"chi_squared" :: Union{Missing, Float64} = missing
    var"exact" :: Union{Missing, Int64} = missing
    var"weight" :: Union{Missing, Float64} = missing
    var"source" :: Union{Missing, String} = missing
    var"measured" :: Union{Missing, Float64} = missing
    var"reconstructed" :: Union{Missing, Float64} = missing
    var"time_measurement" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___constraints__bpol_probe(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
    function equilibrium__time_slice___constraints__b_field_tor_vacuum_r(var"chi_squared"=missing, var"exact"=missing, var"weight"=missing, var"source"=missing, var"measured"=missing, var"reconstructed"=missing, var"time_measurement"=missing, _parent=WeakRef(missing))
        fds = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)
        assign_derived_quantities(fds)
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
        assign_derived_quantities(fds)
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
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__strike_point <: FDSvectorElement
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__geometric_axis <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__gap <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"value" :: Union{Missing, Float64} = missing
    var"identifier" :: Union{Missing, String} = missing
    var"angle" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__gap(var"name"=missing, var"r"=missing, var"value"=missing, var"identifier"=missing, var"angle"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"r", var"value", var"identifier", var"angle", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__dr_dz_zero_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__closest_wall_point <: FDS
    var"distance" :: Union{Missing, Float64} = missing
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__closest_wall_point(var"distance"=missing, var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"distance", var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_separatrix__active_limiter_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_separatrix__active_limiter_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"psi", var"elongation_lower", var"strike_point", var"x_point", var"gap", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"dr_dz_zero_point", var"squareness_lower_outer", var"triangularity_lower", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"closest_wall_point", var"type", _parent)
        assign_derived_quantities(fds)
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
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__strike_point <: FDSvectorElement
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix__outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary_secondary_separatrix <: FDS
    var"psi" :: Union{Missing, Float64} = missing
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__x_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[])
    var"distance_inner_outer" :: Union{Missing, Float64} = missing
    var"outline" :: equilibrium__time_slice___boundary_secondary_separatrix__outline = equilibrium__time_slice___boundary_secondary_separatrix__outline()
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice___boundary_secondary_separatrix__strike_point} = FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary_secondary_separatrix(var"psi"=missing, var"x_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__x_point[]), var"distance_inner_outer"=missing, var"outline"=equilibrium__time_slice___boundary_secondary_separatrix__outline(), var"strike_point"=FDSvector(equilibrium__time_slice___boundary_secondary_separatrix__strike_point[]), _parent=WeakRef(missing))
        fds = new(var"psi", var"x_point", var"distance_inner_outer", var"outline", var"strike_point", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.x_point, :_parent, WeakRef(fds))
        setfield!(fds.outline, :_parent, WeakRef(fds))
        setfield!(fds.strike_point, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__x_point <: FDSvectorElement
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__x_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__strike_point <: FDSvectorElement
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__strike_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__outline <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__outline(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__lcfs <: FDS
    var"r" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"z" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__lcfs(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__geometric_axis <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__geometric_axis(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__time_slice___boundary__active_limiter_point <: FDS
    var"r" :: Union{Missing, Float64} = missing
    var"z" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__time_slice___boundary__active_limiter_point(var"r"=missing, var"z"=missing, _parent=WeakRef(missing))
        fds = new(var"r", var"z", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"psi", var"lcfs", var"elongation_lower", var"strike_point", var"x_point", var"triangularity", var"elongation_upper", var"triangularity_upper", var"outline", var"squareness_lower_outer", var"triangularity_lower", var"psi_norm", var"minor_radius", var"squareness_upper_inner", var"squareness_upper_outer", var"squareness_lower_inner", var"geometric_axis", var"elongation", var"active_limiter_point", var"b_flux_pol_norm", var"type", _parent)
        assign_derived_quantities(fds)
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
        fds = new(var"time", var"ggd", var"profiles_1d", var"boundary", var"constraints", var"global_quantities", var"convergence", var"coordinate_system", var"boundary_secondary_separatrix", var"boundary_separatrix", var"profiles_2d", _parent)
        assign_derived_quantities(fds)
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
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary <: FDSvectorElement
    var"neighbours" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary(var"neighbours"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"neighbours", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension___object <: FDSvectorElement
    var"nodes" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"measure" :: Union{Missing, Float64} = missing
    var"geometry" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension___object(var"nodes"=missing, var"measure"=missing, var"geometry"=missing, var"boundary"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object___boundary[]), _parent=WeakRef(missing))
        fds = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.boundary, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___objects_per_dimension <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension___object} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___objects_per_dimension(var"object"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space___geometry_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space___geometry_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___space <: FDSvectorElement
    var"coordinates_type" :: Union{Missing, AbstractFDVector{Int64}} = missing
    var"geometry_type" :: equilibrium__grids_ggd___grid___space___geometry_type = equilibrium__grids_ggd___grid___space___geometry_type()
    var"identifier" :: equilibrium__grids_ggd___grid___space___identifier = equilibrium__grids_ggd___grid___space___identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___space___objects_per_dimension} = FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___space(var"coordinates_type"=missing, var"geometry_type"=equilibrium__grids_ggd___grid___space___geometry_type(), var"identifier"=equilibrium__grids_ggd___grid___space___identifier(), var"objects_per_dimension"=FDSvector(equilibrium__grids_ggd___grid___space___objects_per_dimension[]), _parent=WeakRef(missing))
        fds = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.geometry_type, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        setfield!(fds.objects_per_dimension, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___metric <: FDS
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___metric(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___identifier <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___identifier(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element___object <: FDSvectorElement
    var"dimension" :: Union{Missing, Int64} = missing
    var"space" :: Union{Missing, Int64} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___element___object(var"dimension"=missing, var"space"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"dimension", var"space", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___element <: FDSvectorElement
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element___object} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___element(var"object"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element___object[]), _parent=WeakRef(missing))
        fds = new(var"object", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.object, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset___base <: FDSvectorElement
    var"jacobian" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"tensor_contravariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    var"tensor_covariant" :: Union{Missing, AbstractArray{Float64, 3}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset___base(var"jacobian"=missing, var"tensor_contravariant"=missing, var"tensor_covariant"=missing, _parent=WeakRef(missing))
        fds = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd___grid___grid_subset <: FDSvectorElement
    var"base" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___base} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[])
    var"metric" :: equilibrium__grids_ggd___grid___grid_subset___metric = equilibrium__grids_ggd___grid___grid_subset___metric()
    var"dimension" :: Union{Missing, Int64} = missing
    var"identifier" :: equilibrium__grids_ggd___grid___grid_subset___identifier = equilibrium__grids_ggd___grid___grid_subset___identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid___grid_subset___element} = FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd___grid___grid_subset(var"base"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___base[]), var"metric"=equilibrium__grids_ggd___grid___grid_subset___metric(), var"dimension"=missing, var"identifier"=equilibrium__grids_ggd___grid___grid_subset___identifier(), var"element"=FDSvector(equilibrium__grids_ggd___grid___grid_subset___element[]), _parent=WeakRef(missing))
        fds = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        assign_derived_quantities(fds)
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
        assign_derived_quantities(fds)
        setfield!(fds.grid_subset, :_parent, WeakRef(fds))
        setfield!(fds.space, :_parent, WeakRef(fds))
        setfield!(fds.identifier, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd <: FDSvectorElement
    var"time" :: Union{Missing, Float64} = missing
    var"grid" :: FDSvector{T} where {T<:equilibrium__grids_ggd___grid} = FDSvector(equilibrium__grids_ggd___grid[])
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__grids_ggd(var"time"=missing, var"grid"=FDSvector(equilibrium__grids_ggd___grid[]), _parent=WeakRef(missing))
        fds = new(var"time", var"grid", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.grid, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct equilibrium__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function equilibrium__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
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
        fds = new(var"time_slice", var"time", var"ids_properties", var"grids_ggd", var"vacuum_toroidal_field", var"code", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_slice, :_parent, WeakRef(fds))
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.grids_ggd, :_parent, WeakRef(fds))
        setfield!(fds.vacuum_toroidal_field, :_parent, WeakRef(fds))
        setfield!(fds.code, :_parent, WeakRef(fds))
        return fds
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
    function dataset_description__simulation(var"time_ended"=missing, var"time_begun"=missing, var"time_current"=missing, var"time_restart"=missing, var"workflow"=missing, var"comment_after"=missing, var"time_begin"=missing, var"time_end"=missing, var"comment_before"=missing, var"time_step"=missing, _parent=WeakRef(missing))
        fds = new(var"time_ended", var"time_begun", var"time_current", var"time_restart", var"workflow", var"comment_after", var"time_begin", var"time_end", var"comment_before", var"time_step", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__pulse_time_end_epoch <: FDS
    var"nanoseconds" :: Union{Missing, Int64} = missing
    var"seconds" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__pulse_time_end_epoch(var"nanoseconds"=missing, var"seconds"=missing, _parent=WeakRef(missing))
        fds = new(var"nanoseconds", var"seconds", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__pulse_time_begin_epoch <: FDS
    var"nanoseconds" :: Union{Missing, Int64} = missing
    var"seconds" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__pulse_time_begin_epoch(var"nanoseconds"=missing, var"seconds"=missing, _parent=WeakRef(missing))
        fds = new(var"nanoseconds", var"seconds", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__parent_entry <: FDS
    var"pulse_type" :: Union{Missing, String} = missing
    var"run" :: Union{Missing, Int64} = missing
    var"machine" :: Union{Missing, String} = missing
    var"pulse" :: Union{Missing, Int64} = missing
    var"user" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__parent_entry(var"pulse_type"=missing, var"run"=missing, var"machine"=missing, var"pulse"=missing, var"user"=missing, _parent=WeakRef(missing))
        fds = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__ids_properties__version_put <: FDS
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct dataset_description__data_entry <: FDS
    var"pulse_type" :: Union{Missing, String} = missing
    var"run" :: Union{Missing, Int64} = missing
    var"machine" :: Union{Missing, String} = missing
    var"pulse" :: Union{Missing, Int64} = missing
    var"user" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function dataset_description__data_entry(var"pulse_type"=missing, var"run"=missing, var"machine"=missing, var"pulse"=missing, var"user"=missing, _parent=WeakRef(missing))
        fds = new(var"pulse_type", var"run", var"machine", var"pulse", var"user", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"pulse_time_begin_epoch", var"imas_version", var"ids_properties", var"time", var"dd_version", var"parent_entry", var"simulation", var"pulse_time_end_epoch", var"data_entry", var"pulse_time_begin", _parent)
        assign_derived_quantities(fds)
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
    var"b0" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"r0" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__vacuum_toroidal_field(var"b0"=missing, var"r0"=missing, _parent=WeakRef(missing))
        fds = new(var"b0", var"r0", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___zeff_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___t_i_average_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___state___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state___neutral_type <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___state___neutral_type(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___state <: FDSvectorElement
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
        fds = new(var"label", var"vibrational_level", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"velocity", var"density", var"density_fast", var"neutral_type", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.neutral_type, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral___element <: FDSvectorElement
    var"atoms_n" :: Union{Missing, Int64} = missing
    var"z_n" :: Union{Missing, Float64} = missing
    var"multiplicity" :: Union{Missing, Float64} = missing
    var"a" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___neutral___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        fds = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___neutral <: FDSvectorElement
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
        fds = new(var"label", var"temperature", var"pressure_thermal", var"ion_index", var"multiple_states_flag", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"pressure_fast_parallel", var"state", var"velocity", var"density", var"density_fast", var"element", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.state, :_parent, WeakRef(fds))
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___temperature_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state___velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___state___density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___state <: FDSvectorElement
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
        fds = new(var"label", var"vibrational_level", var"rotation_frequency_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"z_min", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"z_average_square_1d", var"velocity", var"z_average", var"z_max", var"z_square_average", var"ionisation_potential", var"density", var"density_fast", var"z_average_1d", var"density_fit", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.density_fit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___element <: FDSvectorElement
    var"atoms_n" :: Union{Missing, Int64} = missing
    var"z_n" :: Union{Missing, Float64} = missing
    var"multiplicity" :: Union{Missing, Float64} = missing
    var"a" :: Union{Missing, Float64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___element(var"atoms_n"=missing, var"z_n"=missing, var"multiplicity"=missing, var"a"=missing, _parent=WeakRef(missing))
        fds = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___ion___density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___ion <: FDSvectorElement
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
        fds = new(var"label", var"rotation_frequency_tor", var"temperature_validity", var"velocity_tor", var"temperature", var"z_ion_1d", var"pressure_thermal", var"multiple_states_flag", var"pressure_fast_perpendicular", var"neutral_index", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"state", var"velocity", var"z_ion", var"temperature_fit", var"density", var"velocity_pol", var"density_fast", var"density_fit", var"z_ion_square_1d", var"element", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.state, :_parent, WeakRef(fds))
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.temperature_fit, :_parent, WeakRef(fds))
        setfield!(fds.density_fit, :_parent, WeakRef(fds))
        setfield!(fds.element, :_parent, WeakRef(fds))
        return fds
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
    function core_profiles__profiles_1d___grid(var"psi"=missing, var"psi_boundary"=missing, var"volume"=missing, var"area"=missing, var"rho_pol_norm"=missing, var"rho_tor_norm"=missing, var"surface"=missing, var"rho_tor"=missing, var"psi_magnetic_axis"=missing, _parent=WeakRef(missing))
        fds = new(var"psi", var"psi_boundary", var"volume", var"area", var"rho_pol_norm", var"rho_tor_norm", var"surface", var"rho_tor", var"psi_magnetic_axis", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__velocity <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__velocity(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__temperature_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Missing, String} = missing
    var"description" :: Union{Missing, String} = missing
    var"index" :: Union{Missing, Int64} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___electrons__density_fit__time_measurement_slice_method(var"name"=missing, var"description"=missing, var"index"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"description", var"index", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.time_measurement_slice_method, :_parent, WeakRef(fds))
        return fds
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
        fds = new(var"temperature_validity", var"velocity_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"velocity", var"temperature_fit", var"density", var"velocity_pol", var"collisionality_norm", var"density_fast", var"density_fit", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.velocity, :_parent, WeakRef(fds))
        setfield!(fds.temperature_fit, :_parent, WeakRef(fds))
        setfield!(fds.density_fit, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d___e_field <: FDS
    var"parallel" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"toroidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"diamagnetic" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"radial" :: Union{Missing, AbstractFDVector{Float64}} = missing
    var"poloidal" :: Union{Missing, AbstractFDVector{Float64}} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__profiles_1d___e_field(var"parallel"=missing, var"toroidal"=missing, var"diamagnetic"=missing, var"radial"=missing, var"poloidal"=missing, _parent=WeakRef(missing))
        fds = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d <: FDSvectorElement
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
        fds = new(var"pressure_ion_total", var"time", var"t_i_average_fit", var"neutral", var"n_i_thermal_total", var"magnetic_shear", var"ion", var"j_total", var"rotation_frequency_tor_sonic", var"pressure_thermal", var"j_tor", var"current_parallel_inside", var"j_non_inductive", var"e_field_parallel", var"momentum_tor", var"conductivity_parallel", var"electrons", var"pressure_perpendicular", var"q", var"t_i_average", var"j_ohmic", var"grid", var"phi_potential", var"j_bootstrap", var"zeff_fit", var"pressure_parallel", var"e_field", var"zeff", var"n_i_total_over_n_e", _parent)
        assign_derived_quantities(fds)
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
    var"access_layer_language" :: Union{Missing, String} = missing
    var"data_dictionary" :: Union{Missing, String} = missing
    var"access_layer" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__ids_properties__version_put(var"access_layer_language"=missing, var"data_dictionary"=missing, var"access_layer"=missing, _parent=WeakRef(missing))
        fds = new(var"access_layer_language", var"data_dictionary", var"access_layer", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.version_put, :_parent, WeakRef(fds))
        return fds
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
    function core_profiles__global_quantities(var"beta_tor_norm"=missing, var"resistive_psi_losses"=missing, var"ip"=missing, var"li_3"=missing, var"t_i_average_peaking"=missing, var"t_e_peaking"=missing, var"beta_tor"=missing, var"z_eff_resistive"=missing, var"ejima"=missing, var"energy_diamagnetic"=missing, var"li"=missing, var"current_non_inductive"=missing, var"v_loop"=missing, var"beta_pol"=missing, var"current_bootstrap"=missing, _parent=WeakRef(missing))
        fds = new(var"beta_tor_norm", var"resistive_psi_losses", var"ip", var"li_3", var"t_i_average_peaking", var"t_e_peaking", var"beta_tor", var"z_eff_resistive", var"ejima", var"energy_diamagnetic", var"li", var"current_non_inductive", var"v_loop", var"beta_pol", var"current_bootstrap", _parent)
        assign_derived_quantities(fds)
        return fds
    end
end

Base.@kwdef mutable struct core_profiles__code__library <: FDSvectorElement
    var"name" :: Union{Missing, String} = missing
    var"parameters" :: Union{Missing, String} = missing
    var"commit" :: Union{Missing, String} = missing
    var"repository" :: Union{Missing, String} = missing
    var"version" :: Union{Missing, String} = missing
    _parent :: WeakRef = WeakRef(missing)
    function core_profiles__code__library(var"name"=missing, var"parameters"=missing, var"commit"=missing, var"repository"=missing, var"version"=missing, _parent=WeakRef(missing))
        fds = new(var"name", var"parameters", var"commit", var"repository", var"version", _parent)
        assign_derived_quantities(fds)
        return fds
    end
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
        fds = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.library, :_parent, WeakRef(fds))
        return fds
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
        fds = new(var"time", var"ids_properties", var"vacuum_toroidal_field", var"code", var"global_quantities", var"profiles_1d", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.ids_properties, :_parent, WeakRef(fds))
        setfield!(fds.vacuum_toroidal_field, :_parent, WeakRef(fds))
        setfield!(fds.code, :_parent, WeakRef(fds))
        setfield!(fds.global_quantities, :_parent, WeakRef(fds))
        setfield!(fds.profiles_1d, :_parent, WeakRef(fds))
        return fds
    end
end

Base.@kwdef mutable struct dd <: FDS
    var"equilibrium" :: Union{Missing, equilibrium} = equilibrium()
    var"core_profiles" :: Union{Missing, core_profiles} = core_profiles()
    var"wall" :: Union{Missing, wall} = wall()
    var"dataset_description" :: Union{Missing, dataset_description} = dataset_description()
    _parent :: WeakRef = WeakRef(missing)
    function dd(var"equilibrium"=equilibrium(), var"core_profiles"=core_profiles(), var"wall"=wall(), var"dataset_description"=dataset_description(), _parent=WeakRef(missing))
        fds = new(var"equilibrium", var"core_profiles", var"wall", var"dataset_description", _parent)
        assign_derived_quantities(fds)
        setfield!(fds.equilibrium, :_parent, WeakRef(fds))
        setfield!(fds.core_profiles, :_parent, WeakRef(fds))
        setfield!(fds.wall, :_parent, WeakRef(fds))
        setfield!(fds.dataset_description, :_parent, WeakRef(fds))
        return fds
    end
end

