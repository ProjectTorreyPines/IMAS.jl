abstract type FDS end

include("functionarrays.jl")

supported_types = Union{Float64, Int64, String, Array{Float64, N} where N, Array{Int64, N} where N, Array{String, N} where N}

Base.@kwdef mutable struct equilibrium__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"r0" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__vacuum_toroidal_field(var"b0"=nothing, var"r0"=nothing, _parent=nothing)
        obj = new(var"b0", var"r0", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_2d__grid_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__profiles_2d__grid_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_2d__grid <: FDS
    var"volume_element" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"dim2" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"dim1" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__profiles_2d__grid(var"volume_element"=nothing, var"dim2"=nothing, var"dim1"=nothing, _parent=nothing)
        obj = new(var"volume_element", var"dim2", var"dim1", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_2d <: FDS
    var"psi" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"b_field_r" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"r" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"b_r" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"theta" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"b_field_z" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"j_tor" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"phi" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"b_field_tor" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"b_z" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"grid" :: equilibrium__time_slice__profiles_2d__grid = equilibrium__time_slice__profiles_2d__grid()
    var"grid_type" :: equilibrium__time_slice__profiles_2d__grid_type = equilibrium__time_slice__profiles_2d__grid_type()
    var"j_parallel" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"b_tor" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__profiles_2d(var"psi"=nothing, var"b_field_r"=nothing, var"r"=nothing, var"b_r"=nothing, var"theta"=nothing, var"b_field_z"=nothing, var"j_tor"=nothing, var"phi"=nothing, var"z"=nothing, var"b_field_tor"=nothing, var"b_z"=nothing, var"grid"=equilibrium__time_slice__profiles_2d__grid(), var"grid_type"=equilibrium__time_slice__profiles_2d__grid_type(), var"j_parallel"=nothing, var"b_tor"=nothing, _parent=nothing)
        obj = new(var"psi", var"b_field_r", var"r", var"b_r", var"theta", var"b_field_z", var"j_tor", var"phi", var"z", var"b_field_tor", var"b_z", var"grid", var"grid_type", var"j_parallel", var"b_tor", _parent)
        obj.grid._parent = WeakRef(obj)
        obj.grid_type._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_1d__geometric_axis <: FDS
    var"r" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__profiles_1d__geometric_axis(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__profiles_1d <: FDS
    var"b_field_max" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"dvolume_drho_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm9" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"dpsi_drho_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"surface" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"magnetic_shear" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"b_average" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"b_field_min" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"darea_dpsi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm3" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"squareness_upper_inner" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"squareness_lower_inner" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"elongation" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"beta_pol" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"b_field_average" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"j_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm6" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"psi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm8" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"dpressure_dpsi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"triangularity_upper" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"darea_drho_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"area" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"trapped_fraction" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"volume" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"dvolume_dpsi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"b_min" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"f" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"mass_density" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"r_outboard" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm4" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"phi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"squareness_lower_outer" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"triangularity_lower" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm2" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_volume_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm1" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm5" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"b_max" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"f_df_dpsi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"j_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"r_inboard" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"q" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"gm7" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"squareness_upper_outer" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"geometric_axis" :: equilibrium__time_slice__profiles_1d__geometric_axis = equilibrium__time_slice__profiles_1d__geometric_axis()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__profiles_1d(var"b_field_max"=nothing, var"dvolume_drho_tor"=nothing, var"gm9"=nothing, var"dpsi_drho_tor"=nothing, var"surface"=nothing, var"rho_tor"=nothing, var"magnetic_shear"=nothing, var"b_average"=nothing, var"b_field_min"=nothing, var"darea_dpsi"=nothing, var"gm3"=nothing, var"squareness_upper_inner"=nothing, var"squareness_lower_inner"=nothing, var"rho_tor_norm"=nothing, var"elongation"=nothing, var"beta_pol"=nothing, var"b_field_average"=nothing, var"j_parallel"=nothing, var"gm6"=nothing, var"psi"=nothing, var"gm8"=nothing, var"dpressure_dpsi"=nothing, var"triangularity_upper"=nothing, var"darea_drho_tor"=nothing, var"area"=nothing, var"trapped_fraction"=nothing, var"volume"=nothing, var"dvolume_dpsi"=nothing, var"b_min"=nothing, var"f"=nothing, var"mass_density"=nothing, var"r_outboard"=nothing, var"gm4"=nothing, var"phi"=nothing, var"squareness_lower_outer"=nothing, var"triangularity_lower"=nothing, var"gm2"=nothing, var"rho_volume_norm"=nothing, var"gm1"=nothing, var"gm5"=nothing, var"b_max"=nothing, var"f_df_dpsi"=nothing, var"j_tor"=nothing, var"r_inboard"=nothing, var"q"=nothing, var"gm7"=nothing, var"pressure"=nothing, var"squareness_upper_outer"=nothing, var"geometric_axis"=equilibrium__time_slice__profiles_1d__geometric_axis(), _parent=nothing)
        obj = new(var"b_field_max", var"dvolume_drho_tor", var"gm9", var"dpsi_drho_tor", var"surface", var"rho_tor", var"magnetic_shear", var"b_average", var"b_field_min", var"darea_dpsi", var"gm3", var"squareness_upper_inner", var"squareness_lower_inner", var"rho_tor_norm", var"elongation", var"beta_pol", var"b_field_average", var"j_parallel", var"gm6", var"psi", var"gm8", var"dpressure_dpsi", var"triangularity_upper", var"darea_drho_tor", var"area", var"trapped_fraction", var"volume", var"dvolume_dpsi", var"b_min", var"f", var"mass_density", var"r_outboard", var"gm4", var"phi", var"squareness_lower_outer", var"triangularity_lower", var"gm2", var"rho_volume_norm", var"gm1", var"gm5", var"b_max", var"f_df_dpsi", var"j_tor", var"r_inboard", var"q", var"gm7", var"pressure", var"squareness_upper_outer", var"geometric_axis", _parent)
        obj.geometric_axis._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities__q_min <: FDS
    var"value" :: Union{Nothing, Float64} = nothing
    var"rho_tor_norm" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__global_quantities__q_min(var"value"=nothing, var"rho_tor_norm"=nothing, _parent=nothing)
        obj = new(var"value", var"rho_tor_norm", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities__magnetic_axis <: FDS
    var"b_field_tor" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    var"b_tor" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__global_quantities__magnetic_axis(var"b_field_tor"=nothing, var"r"=nothing, var"z"=nothing, var"b_tor"=nothing, _parent=nothing)
        obj = new(var"b_field_tor", var"r", var"z", var"b_tor", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities__current_centre <: FDS
    var"velocity_z" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__global_quantities__current_centre(var"velocity_z"=nothing, var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"velocity_z", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__global_quantities <: FDS
    var"ip" :: Union{Nothing, Float64} = nothing
    var"li_3" :: Union{Nothing, Float64} = nothing
    var"beta_tor" :: Union{Nothing, Float64} = nothing
    var"surface" :: Union{Nothing, Float64} = nothing
    var"magnetic_axis" :: equilibrium__time_slice__global_quantities__magnetic_axis = equilibrium__time_slice__global_quantities__magnetic_axis()
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
    var"current_centre" :: equilibrium__time_slice__global_quantities__current_centre = equilibrium__time_slice__global_quantities__current_centre()
    var"q_min" :: equilibrium__time_slice__global_quantities__q_min = equilibrium__time_slice__global_quantities__q_min()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__global_quantities(var"ip"=nothing, var"li_3"=nothing, var"beta_tor"=nothing, var"surface"=nothing, var"magnetic_axis"=equilibrium__time_slice__global_quantities__magnetic_axis(), var"energy_mhd"=nothing, var"psi_boundary"=nothing, var"length_pol"=nothing, var"area"=nothing, var"psi_external_average"=nothing, var"q_95"=nothing, var"q_axis"=nothing, var"psi_axis"=nothing, var"w_mhd"=nothing, var"volume"=nothing, var"plasma_inductance"=nothing, var"beta_pol"=nothing, var"beta_normal"=nothing, var"current_centre"=equilibrium__time_slice__global_quantities__current_centre(), var"q_min"=equilibrium__time_slice__global_quantities__q_min(), _parent=nothing)
        obj = new(var"ip", var"li_3", var"beta_tor", var"surface", var"magnetic_axis", var"energy_mhd", var"psi_boundary", var"length_pol", var"area", var"psi_external_average", var"q_95", var"q_axis", var"psi_axis", var"w_mhd", var"volume", var"plasma_inductance", var"beta_pol", var"beta_normal", var"current_centre", var"q_min", _parent)
        obj.magnetic_axis._parent = WeakRef(obj)
        obj.current_centre._parent = WeakRef(obj)
        obj.q_min._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__z <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__z(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__theta <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__theta(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__r <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__r(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__psi <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__psi(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__phi <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__phi(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__j_tor <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__j_tor(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__j_parallel <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__j_parallel(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary <: FDS
    var"neighbours" :: Union{Nothing, Array{Int, 1}} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary(var"neighbours"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"neighbours", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object <: FDS
    var"nodes" :: Union{Nothing, Array{Int, 1}} = nothing
    var"measure" :: Union{Nothing, Float64} = nothing
    var"geometry" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary} = FDSvector(equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object(var"nodes"=nothing, var"measure"=nothing, var"geometry"=nothing, var"boundary"=FDSvector(equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object__boundary[]), _parent=nothing)
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        obj.boundary._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object} = FDSvector(equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__space__objects_per_dimension(var"object"=FDSvector(equilibrium__time_slice__ggd__grid__space__objects_per_dimension__object[]), _parent=nothing)
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__space__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space__geometry_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__space__geometry_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__space <: FDS
    var"coordinates_type" :: Union{Nothing, Array{Int, 1}} = nothing
    var"geometry_type" :: equilibrium__time_slice__ggd__grid__space__geometry_type = equilibrium__time_slice__ggd__grid__space__geometry_type()
    var"identifier" :: equilibrium__time_slice__ggd__grid__space__identifier = equilibrium__time_slice__ggd__grid__space__identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__space__objects_per_dimension} = FDSvector(equilibrium__time_slice__ggd__grid__space__objects_per_dimension[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__space(var"coordinates_type"=nothing, var"geometry_type"=equilibrium__time_slice__ggd__grid__space__geometry_type(), var"identifier"=equilibrium__time_slice__ggd__grid__space__identifier(), var"objects_per_dimension"=FDSvector(equilibrium__time_slice__ggd__grid__space__objects_per_dimension[]), _parent=nothing)
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        obj.geometry_type._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.objects_per_dimension._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__metric <: FDS
    var"jacobian" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__grid_subset__metric(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=nothing)
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__grid_subset__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__element__object <: FDS
    var"dimension" :: Union{Nothing, Int} = nothing
    var"space" :: Union{Nothing, Int} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__grid_subset__element__object(var"dimension"=nothing, var"space"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"dimension", var"space", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__element <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__grid_subset__element__object} = FDSvector(equilibrium__time_slice__ggd__grid__grid_subset__element__object[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__grid_subset__element(var"object"=FDSvector(equilibrium__time_slice__ggd__grid__grid_subset__element__object[]), _parent=nothing)
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset__base <: FDS
    var"jacobian" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__grid_subset__base(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=nothing)
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid__grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__grid_subset__base} = FDSvector(equilibrium__time_slice__ggd__grid__grid_subset__base[])
    var"metric" :: equilibrium__time_slice__ggd__grid__grid_subset__metric = equilibrium__time_slice__ggd__grid__grid_subset__metric()
    var"dimension" :: Union{Nothing, Int} = nothing
    var"identifier" :: equilibrium__time_slice__ggd__grid__grid_subset__identifier = equilibrium__time_slice__ggd__grid__grid_subset__identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__grid_subset__element} = FDSvector(equilibrium__time_slice__ggd__grid__grid_subset__element[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid__grid_subset(var"base"=FDSvector(equilibrium__time_slice__ggd__grid__grid_subset__base[]), var"metric"=equilibrium__time_slice__ggd__grid__grid_subset__metric(), var"dimension"=nothing, var"identifier"=equilibrium__time_slice__ggd__grid__grid_subset__identifier(), var"element"=FDSvector(equilibrium__time_slice__ggd__grid__grid_subset__element[]), _parent=nothing)
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        obj.base._parent = WeakRef(obj)
        obj.metric._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__grid_subset} = FDSvector(equilibrium__time_slice__ggd__grid__grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__grid__space} = FDSvector(equilibrium__time_slice__ggd__grid__space[])
    var"identifier" :: equilibrium__time_slice__ggd__grid__identifier = equilibrium__time_slice__ggd__grid__identifier()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__grid(var"grid_subset"=FDSvector(equilibrium__time_slice__ggd__grid__grid_subset[]), var"space"=FDSvector(equilibrium__time_slice__ggd__grid__space[]), var"identifier"=equilibrium__time_slice__ggd__grid__identifier(), _parent=nothing)
        obj = new(var"grid_subset", var"space", var"identifier", _parent)
        obj.grid_subset._parent = WeakRef(obj)
        obj.space._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__b_field_z <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__b_field_z(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__b_field_tor <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__b_field_tor(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd__b_field_r <: FDS
    var"grid_index" :: Union{Nothing, Int} = nothing
    var"values" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid_subset_index" :: Union{Nothing, Int} = nothing
    var"coefficients" :: Union{Nothing, Array{Float64, 2}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd__b_field_r(var"grid_index"=nothing, var"values"=nothing, var"grid_subset_index"=nothing, var"coefficients"=nothing, _parent=nothing)
        obj = new(var"grid_index", var"values", var"grid_subset_index", var"coefficients", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__ggd <: FDS
    var"b_field_z" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__b_field_z} = FDSvector(equilibrium__time_slice__ggd__b_field_z[])
    var"psi" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__psi} = FDSvector(equilibrium__time_slice__ggd__psi[])
    var"theta" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__theta} = FDSvector(equilibrium__time_slice__ggd__theta[])
    var"z" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__z} = FDSvector(equilibrium__time_slice__ggd__z[])
    var"phi" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__phi} = FDSvector(equilibrium__time_slice__ggd__phi[])
    var"j_tor" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__j_tor} = FDSvector(equilibrium__time_slice__ggd__j_tor[])
    var"grid" :: equilibrium__time_slice__ggd__grid = equilibrium__time_slice__ggd__grid()
    var"b_field_tor" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__b_field_tor} = FDSvector(equilibrium__time_slice__ggd__b_field_tor[])
    var"b_field_r" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__b_field_r} = FDSvector(equilibrium__time_slice__ggd__b_field_r[])
    var"r" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__r} = FDSvector(equilibrium__time_slice__ggd__r[])
    var"j_parallel" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd__j_parallel} = FDSvector(equilibrium__time_slice__ggd__j_parallel[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__ggd(var"b_field_z"=FDSvector(equilibrium__time_slice__ggd__b_field_z[]), var"psi"=FDSvector(equilibrium__time_slice__ggd__psi[]), var"theta"=FDSvector(equilibrium__time_slice__ggd__theta[]), var"z"=FDSvector(equilibrium__time_slice__ggd__z[]), var"phi"=FDSvector(equilibrium__time_slice__ggd__phi[]), var"j_tor"=FDSvector(equilibrium__time_slice__ggd__j_tor[]), var"grid"=equilibrium__time_slice__ggd__grid(), var"b_field_tor"=FDSvector(equilibrium__time_slice__ggd__b_field_tor[]), var"b_field_r"=FDSvector(equilibrium__time_slice__ggd__b_field_r[]), var"r"=FDSvector(equilibrium__time_slice__ggd__r[]), var"j_parallel"=FDSvector(equilibrium__time_slice__ggd__j_parallel[]), _parent=nothing)
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

Base.@kwdef mutable struct equilibrium__time_slice__coordinate_system__grid_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__coordinate_system__grid_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__coordinate_system__grid <: FDS
    var"volume_element" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"dim2" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"dim1" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__coordinate_system__grid(var"volume_element"=nothing, var"dim2"=nothing, var"dim1"=nothing, _parent=nothing)
        obj = new(var"volume_element", var"dim2", var"dim1", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__coordinate_system <: FDS
    var"jacobian" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g13_covariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g11_contravariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g13_contravariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"r" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g12_contravariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g22_contravariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g33_contravariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g22_covariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"tensor_contravariant" :: Union{Nothing, Array{Float64, 4}} = nothing
    var"g12_covariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g33_covariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"grid" :: equilibrium__time_slice__coordinate_system__grid = equilibrium__time_slice__coordinate_system__grid()
    var"g23_covariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"g11_covariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"tensor_covariant" :: Union{Nothing, Array{Float64, 4}} = nothing
    var"g23_contravariant" :: Union{Nothing, Array{Float64, 2}} = nothing
    var"grid_type" :: equilibrium__time_slice__coordinate_system__grid_type = equilibrium__time_slice__coordinate_system__grid_type()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__coordinate_system(var"jacobian"=nothing, var"g13_covariant"=nothing, var"g11_contravariant"=nothing, var"g13_contravariant"=nothing, var"r"=nothing, var"g12_contravariant"=nothing, var"g22_contravariant"=nothing, var"z"=nothing, var"g33_contravariant"=nothing, var"g22_covariant"=nothing, var"tensor_contravariant"=nothing, var"g12_covariant"=nothing, var"g33_covariant"=nothing, var"grid"=equilibrium__time_slice__coordinate_system__grid(), var"g23_covariant"=nothing, var"g11_covariant"=nothing, var"tensor_covariant"=nothing, var"g23_contravariant"=nothing, var"grid_type"=equilibrium__time_slice__coordinate_system__grid_type(), _parent=nothing)
        obj = new(var"jacobian", var"g13_covariant", var"g11_contravariant", var"g13_contravariant", var"r", var"g12_contravariant", var"g22_contravariant", var"z", var"g33_contravariant", var"g22_covariant", var"tensor_contravariant", var"g12_covariant", var"g33_covariant", var"grid", var"g23_covariant", var"g11_covariant", var"tensor_covariant", var"g23_contravariant", var"grid_type", _parent)
        obj.grid._parent = WeakRef(obj)
        obj.grid_type._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__convergence <: FDS
    var"iterations_n" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__convergence(var"iterations_n"=nothing, _parent=nothing)
        obj = new(var"iterations_n", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__x_point__position_reconstructed <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__x_point__position_reconstructed(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__x_point__position_measured <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__x_point__position_measured(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__x_point <: FDS
    var"chi_squared_z" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"position_measured" :: equilibrium__time_slice__constraints__x_point__position_measured = equilibrium__time_slice__constraints__x_point__position_measured()
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    var"chi_squared_r" :: Union{Nothing, Float64} = nothing
    var"position_reconstructed" :: equilibrium__time_slice__constraints__x_point__position_reconstructed = equilibrium__time_slice__constraints__x_point__position_reconstructed()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__x_point(var"chi_squared_z"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"position_measured"=equilibrium__time_slice__constraints__x_point__position_measured(), var"time_measurement"=nothing, var"chi_squared_r"=nothing, var"position_reconstructed"=equilibrium__time_slice__constraints__x_point__position_reconstructed(), _parent=nothing)
        obj = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        obj.position_measured._parent = WeakRef(obj)
        obj.position_reconstructed._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__strike_point__position_reconstructed <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__strike_point__position_reconstructed(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__strike_point__position_measured <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__strike_point__position_measured(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__strike_point <: FDS
    var"chi_squared_z" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"position_measured" :: equilibrium__time_slice__constraints__strike_point__position_measured = equilibrium__time_slice__constraints__strike_point__position_measured()
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    var"chi_squared_r" :: Union{Nothing, Float64} = nothing
    var"position_reconstructed" :: equilibrium__time_slice__constraints__strike_point__position_reconstructed = equilibrium__time_slice__constraints__strike_point__position_reconstructed()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__strike_point(var"chi_squared_z"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"position_measured"=equilibrium__time_slice__constraints__strike_point__position_measured(), var"time_measurement"=nothing, var"chi_squared_r"=nothing, var"position_reconstructed"=equilibrium__time_slice__constraints__strike_point__position_reconstructed(), _parent=nothing)
        obj = new(var"chi_squared_z", var"exact", var"weight", var"source", var"position_measured", var"time_measurement", var"chi_squared_r", var"position_reconstructed", _parent)
        obj.position_measured._parent = WeakRef(obj)
        obj.position_reconstructed._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__q__position <: FDS
    var"phi" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__q__position(var"phi"=nothing, var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"phi", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__q <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"position" :: equilibrium__time_slice__constraints__q__position = equilibrium__time_slice__constraints__q__position()
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__q(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"position"=equilibrium__time_slice__constraints__q__position(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"position", var"time_measurement", _parent)
        obj.position._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__pressure <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__pressure(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__pf_passive_current <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__pf_passive_current(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__pf_current <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__pf_current(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__n_e_line <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__n_e_line(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__n_e <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__n_e(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__mse_polarisation_angle <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__mse_polarisation_angle(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__iron_core_segment <: FDS
    var"magnetisation_r" :: equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r = equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r()
    var"magnetisation_z" :: equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z = equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__iron_core_segment(var"magnetisation_r"=equilibrium__time_slice__constraints__iron_core_segment__magnetisation_r(), var"magnetisation_z"=equilibrium__time_slice__constraints__iron_core_segment__magnetisation_z(), _parent=nothing)
        obj = new(var"magnetisation_r", var"magnetisation_z", _parent)
        obj.magnetisation_r._parent = WeakRef(obj)
        obj.magnetisation_z._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__ip <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__ip(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__flux_loop <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__flux_loop(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__faraday_angle <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__faraday_angle(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__diamagnetic_flux <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__diamagnetic_flux(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__bpol_probe <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__bpol_probe(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints__b_field_tor_vacuum_r <: FDS
    var"chi_squared" :: Union{Nothing, Float64} = nothing
    var"exact" :: Union{Nothing, Int} = nothing
    var"weight" :: Union{Nothing, Float64} = nothing
    var"source" :: Union{Nothing, String} = nothing
    var"measured" :: Union{Nothing, Float64} = nothing
    var"reconstructed" :: Union{Nothing, Float64} = nothing
    var"time_measurement" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints__b_field_tor_vacuum_r(var"chi_squared"=nothing, var"exact"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"reconstructed"=nothing, var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"chi_squared", var"exact", var"weight", var"source", var"measured", var"reconstructed", var"time_measurement", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__constraints <: FDS
    var"faraday_angle" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__faraday_angle} = FDSvector(equilibrium__time_slice__constraints__faraday_angle[])
    var"n_e" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__n_e} = FDSvector(equilibrium__time_slice__constraints__n_e[])
    var"ip" :: equilibrium__time_slice__constraints__ip = equilibrium__time_slice__constraints__ip()
    var"n_e_line" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__n_e_line} = FDSvector(equilibrium__time_slice__constraints__n_e_line[])
    var"pf_current" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__pf_current} = FDSvector(equilibrium__time_slice__constraints__pf_current[])
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__strike_point} = FDSvector(equilibrium__time_slice__constraints__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__x_point} = FDSvector(equilibrium__time_slice__constraints__x_point[])
    var"iron_core_segment" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__iron_core_segment} = FDSvector(equilibrium__time_slice__constraints__iron_core_segment[])
    var"pressure" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__pressure} = FDSvector(equilibrium__time_slice__constraints__pressure[])
    var"diamagnetic_flux" :: equilibrium__time_slice__constraints__diamagnetic_flux = equilibrium__time_slice__constraints__diamagnetic_flux()
    var"pf_passive_current" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__pf_passive_current} = FDSvector(equilibrium__time_slice__constraints__pf_passive_current[])
    var"bpol_probe" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__bpol_probe} = FDSvector(equilibrium__time_slice__constraints__bpol_probe[])
    var"mse_polarisation_angle" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__mse_polarisation_angle} = FDSvector(equilibrium__time_slice__constraints__mse_polarisation_angle[])
    var"q" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__q} = FDSvector(equilibrium__time_slice__constraints__q[])
    var"b_field_tor_vacuum_r" :: equilibrium__time_slice__constraints__b_field_tor_vacuum_r = equilibrium__time_slice__constraints__b_field_tor_vacuum_r()
    var"flux_loop" :: FDSvector{T} where {T<:equilibrium__time_slice__constraints__flux_loop} = FDSvector(equilibrium__time_slice__constraints__flux_loop[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__constraints(var"faraday_angle"=FDSvector(equilibrium__time_slice__constraints__faraday_angle[]), var"n_e"=FDSvector(equilibrium__time_slice__constraints__n_e[]), var"ip"=equilibrium__time_slice__constraints__ip(), var"n_e_line"=FDSvector(equilibrium__time_slice__constraints__n_e_line[]), var"pf_current"=FDSvector(equilibrium__time_slice__constraints__pf_current[]), var"strike_point"=FDSvector(equilibrium__time_slice__constraints__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice__constraints__x_point[]), var"iron_core_segment"=FDSvector(equilibrium__time_slice__constraints__iron_core_segment[]), var"pressure"=FDSvector(equilibrium__time_slice__constraints__pressure[]), var"diamagnetic_flux"=equilibrium__time_slice__constraints__diamagnetic_flux(), var"pf_passive_current"=FDSvector(equilibrium__time_slice__constraints__pf_passive_current[]), var"bpol_probe"=FDSvector(equilibrium__time_slice__constraints__bpol_probe[]), var"mse_polarisation_angle"=FDSvector(equilibrium__time_slice__constraints__mse_polarisation_angle[]), var"q"=FDSvector(equilibrium__time_slice__constraints__q[]), var"b_field_tor_vacuum_r"=equilibrium__time_slice__constraints__b_field_tor_vacuum_r(), var"flux_loop"=FDSvector(equilibrium__time_slice__constraints__flux_loop[]), _parent=nothing)
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

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__x_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__x_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__strike_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__strike_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__outline <: FDS
    var"r" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__outline(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__geometric_axis <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__geometric_axis(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__gap <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"value" :: Union{Nothing, Float64} = nothing
    var"identifier" :: Union{Nothing, String} = nothing
    var"angle" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__gap(var"name"=nothing, var"r"=nothing, var"value"=nothing, var"identifier"=nothing, var"angle"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"name", var"r", var"value", var"identifier", var"angle", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__closest_wall_point <: FDS
    var"distance" :: Union{Nothing, Float64} = nothing
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__closest_wall_point(var"distance"=nothing, var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"distance", var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix__active_limiter_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix__active_limiter_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_separatrix <: FDS
    var"psi" :: Union{Nothing, Float64} = nothing
    var"elongation_lower" :: Union{Nothing, Float64} = nothing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary_separatrix__strike_point} = FDSvector(equilibrium__time_slice__boundary_separatrix__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary_separatrix__x_point} = FDSvector(equilibrium__time_slice__boundary_separatrix__x_point[])
    var"gap" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary_separatrix__gap} = FDSvector(equilibrium__time_slice__boundary_separatrix__gap[])
    var"triangularity" :: Union{Nothing, Float64} = nothing
    var"elongation_upper" :: Union{Nothing, Float64} = nothing
    var"triangularity_upper" :: Union{Nothing, Float64} = nothing
    var"outline" :: equilibrium__time_slice__boundary_separatrix__outline = equilibrium__time_slice__boundary_separatrix__outline()
    var"dr_dz_zero_point" :: equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point = equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point()
    var"squareness_lower_outer" :: Union{Nothing, Float64} = nothing
    var"triangularity_lower" :: Union{Nothing, Float64} = nothing
    var"minor_radius" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_inner" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_outer" :: Union{Nothing, Float64} = nothing
    var"squareness_lower_inner" :: Union{Nothing, Float64} = nothing
    var"geometric_axis" :: equilibrium__time_slice__boundary_separatrix__geometric_axis = equilibrium__time_slice__boundary_separatrix__geometric_axis()
    var"elongation" :: Union{Nothing, Float64} = nothing
    var"active_limiter_point" :: equilibrium__time_slice__boundary_separatrix__active_limiter_point = equilibrium__time_slice__boundary_separatrix__active_limiter_point()
    var"closest_wall_point" :: equilibrium__time_slice__boundary_separatrix__closest_wall_point = equilibrium__time_slice__boundary_separatrix__closest_wall_point()
    var"type" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_separatrix(var"psi"=nothing, var"elongation_lower"=nothing, var"strike_point"=FDSvector(equilibrium__time_slice__boundary_separatrix__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice__boundary_separatrix__x_point[]), var"gap"=FDSvector(equilibrium__time_slice__boundary_separatrix__gap[]), var"triangularity"=nothing, var"elongation_upper"=nothing, var"triangularity_upper"=nothing, var"outline"=equilibrium__time_slice__boundary_separatrix__outline(), var"dr_dz_zero_point"=equilibrium__time_slice__boundary_separatrix__dr_dz_zero_point(), var"squareness_lower_outer"=nothing, var"triangularity_lower"=nothing, var"minor_radius"=nothing, var"squareness_upper_inner"=nothing, var"squareness_upper_outer"=nothing, var"squareness_lower_inner"=nothing, var"geometric_axis"=equilibrium__time_slice__boundary_separatrix__geometric_axis(), var"elongation"=nothing, var"active_limiter_point"=equilibrium__time_slice__boundary_separatrix__active_limiter_point(), var"closest_wall_point"=equilibrium__time_slice__boundary_separatrix__closest_wall_point(), var"type"=nothing, _parent=nothing)
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

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix__x_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_secondary_separatrix__x_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix__strike_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_secondary_separatrix__strike_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix__outline <: FDS
    var"r" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_secondary_separatrix__outline(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary_secondary_separatrix <: FDS
    var"psi" :: Union{Nothing, Float64} = nothing
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary_secondary_separatrix__x_point} = FDSvector(equilibrium__time_slice__boundary_secondary_separatrix__x_point[])
    var"distance_inner_outer" :: Union{Nothing, Float64} = nothing
    var"outline" :: equilibrium__time_slice__boundary_secondary_separatrix__outline = equilibrium__time_slice__boundary_secondary_separatrix__outline()
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary_secondary_separatrix__strike_point} = FDSvector(equilibrium__time_slice__boundary_secondary_separatrix__strike_point[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary_secondary_separatrix(var"psi"=nothing, var"x_point"=FDSvector(equilibrium__time_slice__boundary_secondary_separatrix__x_point[]), var"distance_inner_outer"=nothing, var"outline"=equilibrium__time_slice__boundary_secondary_separatrix__outline(), var"strike_point"=FDSvector(equilibrium__time_slice__boundary_secondary_separatrix__strike_point[]), _parent=nothing)
        obj = new(var"psi", var"x_point", var"distance_inner_outer", var"outline", var"strike_point", _parent)
        obj.x_point._parent = WeakRef(obj)
        obj.outline._parent = WeakRef(obj)
        obj.strike_point._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__x_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary__x_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__strike_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary__strike_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__outline <: FDS
    var"r" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary__outline(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__lcfs <: FDS
    var"r" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary__lcfs(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__geometric_axis <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary__geometric_axis(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary__active_limiter_point <: FDS
    var"r" :: Union{Nothing, Float64} = nothing
    var"z" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary__active_limiter_point(var"r"=nothing, var"z"=nothing, _parent=nothing)
        obj = new(var"r", var"z", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__time_slice__boundary <: FDS
    var"psi" :: Union{Nothing, Float64} = nothing
    var"lcfs" :: equilibrium__time_slice__boundary__lcfs = equilibrium__time_slice__boundary__lcfs()
    var"elongation_lower" :: Union{Nothing, Float64} = nothing
    var"strike_point" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary__strike_point} = FDSvector(equilibrium__time_slice__boundary__strike_point[])
    var"x_point" :: FDSvector{T} where {T<:equilibrium__time_slice__boundary__x_point} = FDSvector(equilibrium__time_slice__boundary__x_point[])
    var"triangularity" :: Union{Nothing, Float64} = nothing
    var"elongation_upper" :: Union{Nothing, Float64} = nothing
    var"triangularity_upper" :: Union{Nothing, Float64} = nothing
    var"outline" :: equilibrium__time_slice__boundary__outline = equilibrium__time_slice__boundary__outline()
    var"squareness_lower_outer" :: Union{Nothing, Float64} = nothing
    var"triangularity_lower" :: Union{Nothing, Float64} = nothing
    var"psi_norm" :: Union{Nothing, Float64} = nothing
    var"minor_radius" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_inner" :: Union{Nothing, Float64} = nothing
    var"squareness_upper_outer" :: Union{Nothing, Float64} = nothing
    var"squareness_lower_inner" :: Union{Nothing, Float64} = nothing
    var"geometric_axis" :: equilibrium__time_slice__boundary__geometric_axis = equilibrium__time_slice__boundary__geometric_axis()
    var"elongation" :: Union{Nothing, Float64} = nothing
    var"active_limiter_point" :: equilibrium__time_slice__boundary__active_limiter_point = equilibrium__time_slice__boundary__active_limiter_point()
    var"b_flux_pol_norm" :: Union{Nothing, Float64} = nothing
    var"type" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice__boundary(var"psi"=nothing, var"lcfs"=equilibrium__time_slice__boundary__lcfs(), var"elongation_lower"=nothing, var"strike_point"=FDSvector(equilibrium__time_slice__boundary__strike_point[]), var"x_point"=FDSvector(equilibrium__time_slice__boundary__x_point[]), var"triangularity"=nothing, var"elongation_upper"=nothing, var"triangularity_upper"=nothing, var"outline"=equilibrium__time_slice__boundary__outline(), var"squareness_lower_outer"=nothing, var"triangularity_lower"=nothing, var"psi_norm"=nothing, var"minor_radius"=nothing, var"squareness_upper_inner"=nothing, var"squareness_upper_outer"=nothing, var"squareness_lower_inner"=nothing, var"geometric_axis"=equilibrium__time_slice__boundary__geometric_axis(), var"elongation"=nothing, var"active_limiter_point"=equilibrium__time_slice__boundary__active_limiter_point(), var"b_flux_pol_norm"=nothing, var"type"=nothing, _parent=nothing)
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
    var"ggd" :: FDSvector{T} where {T<:equilibrium__time_slice__ggd} = FDSvector(equilibrium__time_slice__ggd[])
    var"profiles_1d" :: equilibrium__time_slice__profiles_1d = equilibrium__time_slice__profiles_1d()
    var"boundary" :: equilibrium__time_slice__boundary = equilibrium__time_slice__boundary()
    var"constraints" :: equilibrium__time_slice__constraints = equilibrium__time_slice__constraints()
    var"global_quantities" :: equilibrium__time_slice__global_quantities = equilibrium__time_slice__global_quantities()
    var"convergence" :: equilibrium__time_slice__convergence = equilibrium__time_slice__convergence()
    var"coordinate_system" :: equilibrium__time_slice__coordinate_system = equilibrium__time_slice__coordinate_system()
    var"boundary_secondary_separatrix" :: equilibrium__time_slice__boundary_secondary_separatrix = equilibrium__time_slice__boundary_secondary_separatrix()
    var"boundary_separatrix" :: equilibrium__time_slice__boundary_separatrix = equilibrium__time_slice__boundary_separatrix()
    var"profiles_2d" :: FDSvector{T} where {T<:equilibrium__time_slice__profiles_2d} = FDSvector(equilibrium__time_slice__profiles_2d[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__time_slice(var"time"=nothing, var"ggd"=FDSvector(equilibrium__time_slice__ggd[]), var"profiles_1d"=equilibrium__time_slice__profiles_1d(), var"boundary"=equilibrium__time_slice__boundary(), var"constraints"=equilibrium__time_slice__constraints(), var"global_quantities"=equilibrium__time_slice__global_quantities(), var"convergence"=equilibrium__time_slice__convergence(), var"coordinate_system"=equilibrium__time_slice__coordinate_system(), var"boundary_secondary_separatrix"=equilibrium__time_slice__boundary_secondary_separatrix(), var"boundary_separatrix"=equilibrium__time_slice__boundary_separatrix(), var"profiles_2d"=FDSvector(equilibrium__time_slice__profiles_2d[]), _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__ids_properties__version_put(var"access_layer_language"=nothing, var"data_dictionary"=nothing, var"access_layer"=nothing, _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__ids_properties(var"provider"=nothing, var"version_put"=equilibrium__ids_properties__version_put(), var"homogeneous_time"=nothing, var"source"=nothing, var"creation_date"=nothing, var"comment"=nothing, var"occurrence"=nothing, _parent=nothing)
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        obj.version_put._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary <: FDS
    var"neighbours" :: Union{Nothing, Array{Int, 1}} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary(var"neighbours"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"neighbours", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__objects_per_dimension__object <: FDS
    var"nodes" :: Union{Nothing, Array{Int, 1}} = nothing
    var"measure" :: Union{Nothing, Float64} = nothing
    var"geometry" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"boundary" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary} = FDSvector(equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__space__objects_per_dimension__object(var"nodes"=nothing, var"measure"=nothing, var"geometry"=nothing, var"boundary"=FDSvector(equilibrium__grids_ggd__grid__space__objects_per_dimension__object__boundary[]), _parent=nothing)
        obj = new(var"nodes", var"measure", var"geometry", var"boundary", _parent)
        obj.boundary._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__objects_per_dimension <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__space__objects_per_dimension__object} = FDSvector(equilibrium__grids_ggd__grid__space__objects_per_dimension__object[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__space__objects_per_dimension(var"object"=FDSvector(equilibrium__grids_ggd__grid__space__objects_per_dimension__object[]), _parent=nothing)
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__space__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space__geometry_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__space__geometry_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__space <: FDS
    var"coordinates_type" :: Union{Nothing, Array{Int, 1}} = nothing
    var"geometry_type" :: equilibrium__grids_ggd__grid__space__geometry_type = equilibrium__grids_ggd__grid__space__geometry_type()
    var"identifier" :: equilibrium__grids_ggd__grid__space__identifier = equilibrium__grids_ggd__grid__space__identifier()
    var"objects_per_dimension" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__space__objects_per_dimension} = FDSvector(equilibrium__grids_ggd__grid__space__objects_per_dimension[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__space(var"coordinates_type"=nothing, var"geometry_type"=equilibrium__grids_ggd__grid__space__geometry_type(), var"identifier"=equilibrium__grids_ggd__grid__space__identifier(), var"objects_per_dimension"=FDSvector(equilibrium__grids_ggd__grid__space__objects_per_dimension[]), _parent=nothing)
        obj = new(var"coordinates_type", var"geometry_type", var"identifier", var"objects_per_dimension", _parent)
        obj.geometry_type._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.objects_per_dimension._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__metric <: FDS
    var"jacobian" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__grid_subset__metric(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=nothing)
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__identifier <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__grid_subset__identifier(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__element__object <: FDS
    var"dimension" :: Union{Nothing, Int} = nothing
    var"space" :: Union{Nothing, Int} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__grid_subset__element__object(var"dimension"=nothing, var"space"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"dimension", var"space", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__element <: FDS
    var"object" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__grid_subset__element__object} = FDSvector(equilibrium__grids_ggd__grid__grid_subset__element__object[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__grid_subset__element(var"object"=FDSvector(equilibrium__grids_ggd__grid__grid_subset__element__object[]), _parent=nothing)
        obj = new(var"object", _parent)
        obj.object._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset__base <: FDS
    var"jacobian" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"tensor_contravariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    var"tensor_covariant" :: Union{Nothing, Array{Float64, 3}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__grid_subset__base(var"jacobian"=nothing, var"tensor_contravariant"=nothing, var"tensor_covariant"=nothing, _parent=nothing)
        obj = new(var"jacobian", var"tensor_contravariant", var"tensor_covariant", _parent)

        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid__grid_subset <: FDS
    var"base" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__grid_subset__base} = FDSvector(equilibrium__grids_ggd__grid__grid_subset__base[])
    var"metric" :: equilibrium__grids_ggd__grid__grid_subset__metric = equilibrium__grids_ggd__grid__grid_subset__metric()
    var"dimension" :: Union{Nothing, Int} = nothing
    var"identifier" :: equilibrium__grids_ggd__grid__grid_subset__identifier = equilibrium__grids_ggd__grid__grid_subset__identifier()
    var"element" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__grid_subset__element} = FDSvector(equilibrium__grids_ggd__grid__grid_subset__element[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid__grid_subset(var"base"=FDSvector(equilibrium__grids_ggd__grid__grid_subset__base[]), var"metric"=equilibrium__grids_ggd__grid__grid_subset__metric(), var"dimension"=nothing, var"identifier"=equilibrium__grids_ggd__grid__grid_subset__identifier(), var"element"=FDSvector(equilibrium__grids_ggd__grid__grid_subset__element[]), _parent=nothing)
        obj = new(var"base", var"metric", var"dimension", var"identifier", var"element", _parent)
        obj.base._parent = WeakRef(obj)
        obj.metric._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd__grid <: FDS
    var"grid_subset" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__grid_subset} = FDSvector(equilibrium__grids_ggd__grid__grid_subset[])
    var"space" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid__space} = FDSvector(equilibrium__grids_ggd__grid__space[])
    var"identifier" :: equilibrium__grids_ggd__grid__identifier = equilibrium__grids_ggd__grid__identifier()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd__grid(var"grid_subset"=FDSvector(equilibrium__grids_ggd__grid__grid_subset[]), var"space"=FDSvector(equilibrium__grids_ggd__grid__space[]), var"identifier"=equilibrium__grids_ggd__grid__identifier(), _parent=nothing)
        obj = new(var"grid_subset", var"space", var"identifier", _parent)
        obj.grid_subset._parent = WeakRef(obj)
        obj.space._parent = WeakRef(obj)
        obj.identifier._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium__grids_ggd <: FDS
    var"time" :: Union{Nothing, Float64} = nothing
    var"grid" :: FDSvector{T} where {T<:equilibrium__grids_ggd__grid} = FDSvector(equilibrium__grids_ggd__grid[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__grids_ggd(var"time"=nothing, var"grid"=FDSvector(equilibrium__grids_ggd__grid[]), _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__code__library(var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"version"=nothing, _parent=nothing)
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
    var"output_flag" :: Union{Nothing, Array{Int, 1}} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium__code(var"library"=FDSvector(equilibrium__code__library[]), var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"output_flag"=nothing, var"version"=nothing, _parent=nothing)
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        obj.library._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct equilibrium <: FDS
    var"time_slice" :: FDSvector{T} where {T<:equilibrium__time_slice} = FDSvector(equilibrium__time_slice[])
    var"time" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"ids_properties" :: equilibrium__ids_properties = equilibrium__ids_properties()
    var"grids_ggd" :: FDSvector{T} where {T<:equilibrium__grids_ggd} = FDSvector(equilibrium__grids_ggd[])
    var"vacuum_toroidal_field" :: equilibrium__vacuum_toroidal_field = equilibrium__vacuum_toroidal_field()
    var"code" :: equilibrium__code = equilibrium__code()
    _parent :: Union{Nothing, WeakRef} = nothing
    function equilibrium(var"time_slice"=FDSvector(equilibrium__time_slice[]), var"time"=nothing, var"ids_properties"=equilibrium__ids_properties(), var"grids_ggd"=FDSvector(equilibrium__grids_ggd[]), var"vacuum_toroidal_field"=equilibrium__vacuum_toroidal_field(), var"code"=equilibrium__code(), _parent=nothing)
        obj = new(var"time_slice", var"time", var"ids_properties", var"grids_ggd", var"vacuum_toroidal_field", var"code", _parent)
        obj.time_slice._parent = WeakRef(obj)
        obj.ids_properties._parent = WeakRef(obj)
        obj.grids_ggd._parent = WeakRef(obj)
        obj.vacuum_toroidal_field._parent = WeakRef(obj)
        obj.code._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__vacuum_toroidal_field <: FDS
    var"b0" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"r0" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__vacuum_toroidal_field(var"b0"=nothing, var"r0"=nothing, _parent=nothing)
        obj = new(var"b0", var"r0", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__zeff_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method = core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__zeff_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__zeff_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__t_i_average_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method = core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__t_i_average_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__t_i_average_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__velocity <: FDS
    var"parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__neutral__velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=nothing)
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__state__velocity <: FDS
    var"parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__neutral__state__velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=nothing)
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__state__neutral_type <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__neutral__state__neutral_type(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__state <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"vibrational_level" :: Union{Nothing, Float64} = nothing
    var"temperature" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"electron_configuration" :: Union{Nothing, String} = nothing
    var"pressure" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"vibrational_mode" :: Union{Nothing, String} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"velocity" :: core_profiles__profiles_1d__neutral__state__velocity = core_profiles__profiles_1d__neutral__state__velocity()
    var"density" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"neutral_type" :: core_profiles__profiles_1d__neutral__state__neutral_type = core_profiles__profiles_1d__neutral__state__neutral_type()
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__neutral__state(var"label"=nothing, var"vibrational_level"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"pressure_fast_perpendicular"=nothing, var"electron_configuration"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"vibrational_mode"=nothing, var"pressure_fast_parallel"=nothing, var"velocity"=core_profiles__profiles_1d__neutral__state__velocity(), var"density"=nothing, var"density_fast"=nothing, var"neutral_type"=core_profiles__profiles_1d__neutral__state__neutral_type(), _parent=nothing)
        obj = new(var"label", var"vibrational_level", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"velocity", var"density", var"density_fast", var"neutral_type", _parent)
        obj.velocity._parent = WeakRef(obj)
        obj.neutral_type._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral__element <: FDS
    var"atoms_n" :: Union{Nothing, Int} = nothing
    var"z_n" :: Union{Nothing, Float64} = nothing
    var"multiplicity" :: Union{Nothing, Float64} = nothing
    var"a" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__neutral__element(var"atoms_n"=nothing, var"z_n"=nothing, var"multiplicity"=nothing, var"a"=nothing, _parent=nothing)
        obj = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__neutral <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"temperature" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"ion_index" :: Union{Nothing, Int} = nothing
    var"multiple_states_flag" :: Union{Nothing, Int} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d__neutral__state} = FDSvector(core_profiles__profiles_1d__neutral__state[])
    var"velocity" :: core_profiles__profiles_1d__neutral__velocity = core_profiles__profiles_1d__neutral__velocity()
    var"density" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d__neutral__element} = FDSvector(core_profiles__profiles_1d__neutral__element[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__neutral(var"label"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"ion_index"=nothing, var"multiple_states_flag"=nothing, var"pressure_fast_perpendicular"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"pressure_fast_parallel"=nothing, var"state"=FDSvector(core_profiles__profiles_1d__neutral__state[]), var"velocity"=core_profiles__profiles_1d__neutral__velocity(), var"density"=nothing, var"density_fast"=nothing, var"element"=FDSvector(core_profiles__profiles_1d__neutral__element[]), _parent=nothing)
        obj = new(var"label", var"temperature", var"pressure_thermal", var"ion_index", var"multiple_states_flag", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"pressure_fast_parallel", var"state", var"velocity", var"density", var"density_fast", var"element", _parent)
        obj.state._parent = WeakRef(obj)
        obj.velocity._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__velocity <: FDS
    var"parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=nothing)
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__temperature_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__temperature_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__ion__temperature_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state__velocity <: FDS
    var"parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__state__velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=nothing)
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state__density_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method = core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__state__density_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__ion__state__density_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__state <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"vibrational_level" :: Union{Nothing, Float64} = nothing
    var"rotation_frequency_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z_min" :: Union{Nothing, Float64} = nothing
    var"electron_configuration" :: Union{Nothing, String} = nothing
    var"pressure" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"vibrational_mode" :: Union{Nothing, String} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z_average_square_1d" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"velocity" :: core_profiles__profiles_1d__ion__state__velocity = core_profiles__profiles_1d__ion__state__velocity()
    var"z_average" :: Union{Nothing, Float64} = nothing
    var"z_max" :: Union{Nothing, Float64} = nothing
    var"z_square_average" :: Union{Nothing, Float64} = nothing
    var"ionisation_potential" :: Union{Nothing, Float64} = nothing
    var"density" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z_average_1d" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fit" :: core_profiles__profiles_1d__ion__state__density_fit = core_profiles__profiles_1d__ion__state__density_fit()
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__state(var"label"=nothing, var"vibrational_level"=nothing, var"rotation_frequency_tor"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"pressure_fast_perpendicular"=nothing, var"z_min"=nothing, var"electron_configuration"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"vibrational_mode"=nothing, var"pressure_fast_parallel"=nothing, var"z_average_square_1d"=nothing, var"velocity"=core_profiles__profiles_1d__ion__state__velocity(), var"z_average"=nothing, var"z_max"=nothing, var"z_square_average"=nothing, var"ionisation_potential"=nothing, var"density"=nothing, var"density_fast"=nothing, var"z_average_1d"=nothing, var"density_fit"=core_profiles__profiles_1d__ion__state__density_fit(), _parent=nothing)
        obj = new(var"label", var"vibrational_level", var"rotation_frequency_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"z_min", var"electron_configuration", var"pressure", var"density_thermal", var"vibrational_mode", var"pressure_fast_parallel", var"z_average_square_1d", var"velocity", var"z_average", var"z_max", var"z_square_average", var"ionisation_potential", var"density", var"density_fast", var"z_average_1d", var"density_fit", _parent)
        obj.velocity._parent = WeakRef(obj)
        obj.density_fit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__element <: FDS
    var"atoms_n" :: Union{Nothing, Int} = nothing
    var"z_n" :: Union{Nothing, Float64} = nothing
    var"multiplicity" :: Union{Nothing, Float64} = nothing
    var"a" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__element(var"atoms_n"=nothing, var"z_n"=nothing, var"multiplicity"=nothing, var"a"=nothing, _parent=nothing)
        obj = new(var"atoms_n", var"z_n", var"multiplicity", var"a", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion__density_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method = core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion__density_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__ion__density_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__ion <: FDS
    var"label" :: Union{Nothing, String} = nothing
    var"rotation_frequency_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"temperature_validity" :: Union{Nothing, Int} = nothing
    var"velocity_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z_ion_1d" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"multiple_states_flag" :: Union{Nothing, Int} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"neutral_index" :: Union{Nothing, Int} = nothing
    var"pressure" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_validity" :: Union{Nothing, Int} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"state" :: FDSvector{T} where {T<:core_profiles__profiles_1d__ion__state} = FDSvector(core_profiles__profiles_1d__ion__state[])
    var"velocity" :: core_profiles__profiles_1d__ion__velocity = core_profiles__profiles_1d__ion__velocity()
    var"z_ion" :: Union{Nothing, Float64} = nothing
    var"temperature_fit" :: core_profiles__profiles_1d__ion__temperature_fit = core_profiles__profiles_1d__ion__temperature_fit()
    var"density" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"velocity_pol" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fit" :: core_profiles__profiles_1d__ion__density_fit = core_profiles__profiles_1d__ion__density_fit()
    var"z_ion_square_1d" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"element" :: FDSvector{T} where {T<:core_profiles__profiles_1d__ion__element} = FDSvector(core_profiles__profiles_1d__ion__element[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__ion(var"label"=nothing, var"rotation_frequency_tor"=nothing, var"temperature_validity"=nothing, var"velocity_tor"=nothing, var"temperature"=nothing, var"z_ion_1d"=nothing, var"pressure_thermal"=nothing, var"multiple_states_flag"=nothing, var"pressure_fast_perpendicular"=nothing, var"neutral_index"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"density_validity"=nothing, var"pressure_fast_parallel"=nothing, var"state"=FDSvector(core_profiles__profiles_1d__ion__state[]), var"velocity"=core_profiles__profiles_1d__ion__velocity(), var"z_ion"=nothing, var"temperature_fit"=core_profiles__profiles_1d__ion__temperature_fit(), var"density"=nothing, var"velocity_pol"=nothing, var"density_fast"=nothing, var"density_fit"=core_profiles__profiles_1d__ion__density_fit(), var"z_ion_square_1d"=nothing, var"element"=FDSvector(core_profiles__profiles_1d__ion__element[]), _parent=nothing)
        obj = new(var"label", var"rotation_frequency_tor", var"temperature_validity", var"velocity_tor", var"temperature", var"z_ion_1d", var"pressure_thermal", var"multiple_states_flag", var"pressure_fast_perpendicular", var"neutral_index", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"state", var"velocity", var"z_ion", var"temperature_fit", var"density", var"velocity_pol", var"density_fast", var"density_fit", var"z_ion_square_1d", var"element", _parent)
        obj.state._parent = WeakRef(obj)
        obj.velocity._parent = WeakRef(obj)
        obj.temperature_fit._parent = WeakRef(obj)
        obj.density_fit._parent = WeakRef(obj)
        obj.element._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__grid <: FDS
    var"psi" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"psi_boundary" :: Union{Nothing, Float64} = nothing
    var"volume" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"area" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_pol_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"surface" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"psi_magnetic_axis" :: Union{Nothing, Float64} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__grid(var"psi"=nothing, var"psi_boundary"=nothing, var"volume"=nothing, var"area"=nothing, var"rho_pol_norm"=nothing, var"rho_tor_norm"=nothing, var"surface"=nothing, var"rho_tor"=nothing, var"psi_magnetic_axis"=nothing, _parent=nothing)
        obj = new(var"psi", var"psi_boundary", var"volume", var"area", var"rho_pol_norm", var"rho_tor_norm", var"surface", var"rho_tor", var"psi_magnetic_axis", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__velocity <: FDS
    var"parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__electrons__velocity(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=nothing)
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__temperature_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method = core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__electrons__temperature_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__electrons__temperature_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method <: FDS
    var"name" :: Union{Nothing, String} = nothing
    var"description" :: Union{Nothing, String} = nothing
    var"index" :: Union{Nothing, Int} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method(var"name"=nothing, var"description"=nothing, var"index"=nothing, _parent=nothing)
        obj = new(var"name", var"description", var"index", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons__density_fit <: FDS
    var"local" :: Union{Nothing, Array{Int, 1}} = nothing
    var"chi_squared" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"parameters" :: Union{Nothing, String} = nothing
    var"reconstructed" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_width" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rho_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"weight" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"source" :: Union{Nothing, Array{String, 1}} = nothing
    var"measured" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time_measurement_slice_method" :: core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method = core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method()
    var"time_measurement" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__electrons__density_fit(var"local"=nothing, var"chi_squared"=nothing, var"parameters"=nothing, var"reconstructed"=nothing, var"time_measurement_width"=nothing, var"rho_tor_norm"=nothing, var"weight"=nothing, var"source"=nothing, var"measured"=nothing, var"time_measurement_slice_method"=core_profiles__profiles_1d__electrons__density_fit__time_measurement_slice_method(), var"time_measurement"=nothing, _parent=nothing)
        obj = new(var"local", var"chi_squared", var"parameters", var"reconstructed", var"time_measurement_width", var"rho_tor_norm", var"weight", var"source", var"measured", var"time_measurement_slice_method", var"time_measurement", _parent)
        obj.time_measurement_slice_method._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__electrons <: FDS
    var"temperature_validity" :: Union{Nothing, Int} = nothing
    var"velocity_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"temperature" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_fast_perpendicular" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_validity" :: Union{Nothing, Int} = nothing
    var"pressure_fast_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"velocity" :: core_profiles__profiles_1d__electrons__velocity = core_profiles__profiles_1d__electrons__velocity()
    var"temperature_fit" :: core_profiles__profiles_1d__electrons__temperature_fit = core_profiles__profiles_1d__electrons__temperature_fit()
    var"density" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"velocity_pol" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"collisionality_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fast" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"density_fit" :: core_profiles__profiles_1d__electrons__density_fit = core_profiles__profiles_1d__electrons__density_fit()
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__electrons(var"temperature_validity"=nothing, var"velocity_tor"=nothing, var"temperature"=nothing, var"pressure_thermal"=nothing, var"pressure_fast_perpendicular"=nothing, var"pressure"=nothing, var"density_thermal"=nothing, var"density_validity"=nothing, var"pressure_fast_parallel"=nothing, var"velocity"=core_profiles__profiles_1d__electrons__velocity(), var"temperature_fit"=core_profiles__profiles_1d__electrons__temperature_fit(), var"density"=nothing, var"velocity_pol"=nothing, var"collisionality_norm"=nothing, var"density_fast"=nothing, var"density_fit"=core_profiles__profiles_1d__electrons__density_fit(), _parent=nothing)
        obj = new(var"temperature_validity", var"velocity_tor", var"temperature", var"pressure_thermal", var"pressure_fast_perpendicular", var"pressure", var"density_thermal", var"density_validity", var"pressure_fast_parallel", var"velocity", var"temperature_fit", var"density", var"velocity_pol", var"collisionality_norm", var"density_fast", var"density_fit", _parent)
        obj.velocity._parent = WeakRef(obj)
        obj.temperature_fit._parent = WeakRef(obj)
        obj.density_fit._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d__e_field <: FDS
    var"parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"toroidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"radial" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"poloidal" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d__e_field(var"parallel"=nothing, var"toroidal"=nothing, var"diamagnetic"=nothing, var"radial"=nothing, var"poloidal"=nothing, _parent=nothing)
        obj = new(var"parallel", var"toroidal", var"diamagnetic", var"radial", var"poloidal", _parent)

        return obj
    end
end

Base.@kwdef mutable struct core_profiles__profiles_1d <: FDS
    var"pressure_ion_total" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"time" :: Union{Nothing, Float64} = nothing
    var"t_i_average_fit" :: core_profiles__profiles_1d__t_i_average_fit = core_profiles__profiles_1d__t_i_average_fit()
    var"neutral" :: FDSvector{T} where {T<:core_profiles__profiles_1d__neutral} = FDSvector(core_profiles__profiles_1d__neutral[])
    var"n_i_thermal_total" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"magnetic_shear" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"ion" :: FDSvector{T} where {T<:core_profiles__profiles_1d__ion} = FDSvector(core_profiles__profiles_1d__ion[])
    var"j_total" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"rotation_frequency_tor_sonic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"pressure_thermal" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"j_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"current_parallel_inside" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"j_non_inductive" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"e_field_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"momentum_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"conductivity_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"electrons" :: core_profiles__profiles_1d__electrons = core_profiles__profiles_1d__electrons()
    var"pressure_perpendicular" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"q" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"t_i_average" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"j_ohmic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"grid" :: core_profiles__profiles_1d__grid = core_profiles__profiles_1d__grid()
    var"phi_potential" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"j_bootstrap" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"zeff_fit" :: core_profiles__profiles_1d__zeff_fit = core_profiles__profiles_1d__zeff_fit()
    var"pressure_parallel" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"e_field" :: core_profiles__profiles_1d__e_field = core_profiles__profiles_1d__e_field()
    var"zeff" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"n_i_total_over_n_e" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__profiles_1d(var"pressure_ion_total"=nothing, var"time"=nothing, var"t_i_average_fit"=core_profiles__profiles_1d__t_i_average_fit(), var"neutral"=FDSvector(core_profiles__profiles_1d__neutral[]), var"n_i_thermal_total"=nothing, var"magnetic_shear"=nothing, var"ion"=FDSvector(core_profiles__profiles_1d__ion[]), var"j_total"=nothing, var"rotation_frequency_tor_sonic"=nothing, var"pressure_thermal"=nothing, var"j_tor"=nothing, var"current_parallel_inside"=nothing, var"j_non_inductive"=nothing, var"e_field_parallel"=nothing, var"momentum_tor"=nothing, var"conductivity_parallel"=nothing, var"electrons"=core_profiles__profiles_1d__electrons(), var"pressure_perpendicular"=nothing, var"q"=nothing, var"t_i_average"=nothing, var"j_ohmic"=nothing, var"grid"=core_profiles__profiles_1d__grid(), var"phi_potential"=nothing, var"j_bootstrap"=nothing, var"zeff_fit"=core_profiles__profiles_1d__zeff_fit(), var"pressure_parallel"=nothing, var"e_field"=core_profiles__profiles_1d__e_field(), var"zeff"=nothing, var"n_i_total_over_n_e"=nothing, _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__ids_properties__version_put(var"access_layer_language"=nothing, var"data_dictionary"=nothing, var"access_layer"=nothing, _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__ids_properties(var"provider"=nothing, var"version_put"=core_profiles__ids_properties__version_put(), var"homogeneous_time"=nothing, var"source"=nothing, var"creation_date"=nothing, var"comment"=nothing, var"occurrence"=nothing, _parent=nothing)
        obj = new(var"provider", var"version_put", var"homogeneous_time", var"source", var"creation_date", var"comment", var"occurrence", _parent)
        obj.version_put._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles__global_quantities <: FDS
    var"beta_tor_norm" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"resistive_psi_losses" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"ip" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"li_3" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"t_i_average_peaking" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"t_e_peaking" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"beta_tor" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"z_eff_resistive" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"ejima" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"energy_diamagnetic" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"li" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"current_non_inductive" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"v_loop" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"beta_pol" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"current_bootstrap" :: Union{Nothing, Array{Float64, 1}} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__global_quantities(var"beta_tor_norm"=nothing, var"resistive_psi_losses"=nothing, var"ip"=nothing, var"li_3"=nothing, var"t_i_average_peaking"=nothing, var"t_e_peaking"=nothing, var"beta_tor"=nothing, var"z_eff_resistive"=nothing, var"ejima"=nothing, var"energy_diamagnetic"=nothing, var"li"=nothing, var"current_non_inductive"=nothing, var"v_loop"=nothing, var"beta_pol"=nothing, var"current_bootstrap"=nothing, _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__code__library(var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"version"=nothing, _parent=nothing)
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
    var"output_flag" :: Union{Nothing, Array{Int, 1}} = nothing
    var"version" :: Union{Nothing, String} = nothing
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles__code(var"library"=FDSvector(core_profiles__code__library[]), var"name"=nothing, var"parameters"=nothing, var"commit"=nothing, var"repository"=nothing, var"output_flag"=nothing, var"version"=nothing, _parent=nothing)
        obj = new(var"library", var"name", var"parameters", var"commit", var"repository", var"output_flag", var"version", _parent)
        obj.library._parent = WeakRef(obj)
        return obj
    end
end

Base.@kwdef mutable struct core_profiles <: FDS
    var"time" :: Union{Nothing, Array{Float64, 1}} = nothing
    var"ids_properties" :: core_profiles__ids_properties = core_profiles__ids_properties()
    var"vacuum_toroidal_field" :: core_profiles__vacuum_toroidal_field = core_profiles__vacuum_toroidal_field()
    var"code" :: core_profiles__code = core_profiles__code()
    var"global_quantities" :: core_profiles__global_quantities = core_profiles__global_quantities()
    var"profiles_1d" :: FDSvector{T} where {T<:core_profiles__profiles_1d} = FDSvector(core_profiles__profiles_1d[])
    _parent :: Union{Nothing, WeakRef} = nothing
    function core_profiles(var"time"=nothing, var"ids_properties"=core_profiles__ids_properties(), var"vacuum_toroidal_field"=core_profiles__vacuum_toroidal_field(), var"code"=core_profiles__code(), var"global_quantities"=core_profiles__global_quantities(), var"profiles_1d"=FDSvector(core_profiles__profiles_1d[]), _parent=nothing)
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
    _parent :: Union{Nothing, WeakRef} = nothing
    function dd(var"equilibrium"=equilibrium(), var"core_profiles"=core_profiles(), _parent=nothing)
        obj = new(var"equilibrium", var"core_profiles", _parent)
        obj.equilibrium._parent = WeakRef(obj)
        obj.core_profiles._parent = WeakRef(obj)
        return obj
    end
end

