using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

#= ==================================== =#
# Header
#= ==================================== =#
include("functionarrays.jl")

#= ==================================== =#
# FUSE data structure
#= ==================================== =#
filenames = readdir("$omas_imas_structure_folder/$imas_version")
desired_structure = []
for filename in filenames
    if startswith(filename, "_") || ! endswith(filename, ".json")
        continue
    end
    if ! (filename in ["equilibrium.json", "core_profiles.json", "wall.json", "dataset_description.json"])
        continue
    end
    filename = replace(filename, ".json" => "")
    for item in keys(load_imasdd(filename))
        push!(desired_structure, item)
    end
end

#= ================================================= =#
#  Implementing selected IMAS DD as Julia structures
#= ================================================= =#

"""
    desired_structure::Vector{String}

Translate the IMAS `desired_structure` entries into Julia structs
"""
function imas_julia_struct(desired_structure::Vector)
    convertsion_types = Union{}
    struct_commands = String[]

    branches = Vector()
    push!(branches, "")

    # generate hierarchy of dicts
    ddict = Dict()
    for sel in sort(desired_structure)
        if endswith(sel, "error_upper") | endswith(sel, "error_lower") | endswith(sel, "error_index")
            continue
        end
        path = split(sel, ".")
        
        imasdd = load_imasdd(String(path[1]))

        h = ddict
        for (k, item) in enumerate(path)
            if (k > 1) & (k == length(path))
                if imasdd[sel]["data_type"] in ["STRUCTURE", "STRUCT_ARRAY"]
                    continue

                elseif imasdd[sel]["data_type"] == "STR_0D"
                    h[item] = ":: Union{Nothing, String} = nothing"
                    convertsion_types = Union{convertsion_types,String}
                elseif imasdd[sel]["data_type"] == "STR_1D"
                    h[item] = ":: Union{Nothing, AbstractArray{String, 1}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{String}}

                elseif imasdd[sel]["data_type"] in ["INT_0D", "INT_TYPE"]
                    h[item] = ":: Union{Nothing, Int} = nothing"
                    convertsion_types = Union{convertsion_types,Int}
                elseif imasdd[sel]["data_type"] == "INT_1D"
                    h[item] = ":: Union{Nothing, AbstractArray{Int, 1}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Int}}
                elseif imasdd[sel]["data_type"] == "INT_2D"
                    h[item] = ":: Union{Nothing, AbstractArray{Int, 2}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Int}}
                elseif imasdd[sel]["data_type"] == "INT_3D"
                    h[item] = ":: Union{Nothing, AbstractArray{Int, 3}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Int}}
                elseif imasdd[sel]["data_type"] == "INT_4D"
                    h[item] = ":: Union{Nothing, AbstractArray{Int, 4}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Int}}
                elseif imasdd[sel]["data_type"] == "INT_5D"
                    h[item] = ":: Union{Nothing, AbstractArray{Int, 5}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Int}}
                elseif imasdd[sel]["data_type"] == "INT_6D"
                    h[item] = ":: Union{Nothing, AbstractArray{Int, 6}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Int}}

                elseif imasdd[sel]["data_type"] == "FLT_0D"
                    h[item] = ":: Union{Nothing, Float64} = nothing"
                    convertsion_types = Union{convertsion_types,Float64}
                elseif imasdd[sel]["data_type"] in ["FLT_1D", "FLT_1D_TYPE"]
                    if false
                        h[item] = ":: Union{Nothing, AbstractFunctionArray{Float64,1}} = nothing"
                        convertsion_types = Union{convertsion_types,AbstractFunctionArray{Float64}}
                    else
                        h[item] = ":: Union{Nothing, AbstractArray{Float64, 1}} = nothing"
                        convertsion_types = Union{convertsion_types,Array{Float64}}
                    end
                elseif imasdd[sel]["data_type"] == "FLT_2D"
                    h[item] = ":: Union{Nothing, AbstractArray{Float64, 2}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Float64}}
                elseif imasdd[sel]["data_type"] == "FLT_3D"
                    h[item] = ":: Union{Nothing, AbstractArray{Float64, 3}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Float64}}
                elseif imasdd[sel]["data_type"] == "FLT_4D"
                    h[item] = ":: Union{Nothing, AbstractArray{Float64, 4}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Float64}}
                elseif imasdd[sel]["data_type"] == "FLT_5D"
                    h[item] = ":: Union{Nothing, AbstractArray{Float64, 5}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Float64}}
                elseif imasdd[sel]["data_type"] == "FLT_6D"
                    h[item] = ":: Union{Nothing, AbstractArray{Float64, 6}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Float64}}

                elseif imasdd[sel]["data_type"] == "CPX_0D"
                    h[item] = ":: Union{Nothing, Complex{Float64} = nothing"
                    convertsion_types = Union{convertsion_types,Complex{Float64}}
                elseif imasdd[sel]["data_type"] == "CPX_1D"
                    h[item] = ":: Union{Nothing, AbstractArray{Complex{Float64}, 1}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Complex{Float64}}}
                elseif imasdd[sel]["data_type"] == "CPX_2D"
                    h[item] = ":: Union{Nothing, AbstractArray{Complex{Float64}, 2}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Complex{Float64}}}
                elseif imasdd[sel]["data_type"] == "CPX_3D"
                    h[item] = ":: Union{Nothing, AbstractArray{Complex{Float64}, 3}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Complex{Float64}}}
                elseif imasdd[sel]["data_type"] == "CPX_4D"
                    h[item] = ":: Union{Nothing, AbstractArray{Complex{Float64}, 4}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Complex{Float64}}}
                elseif imasdd[sel]["data_type"] == "CPX_5D"
                    h[item] = ":: Union{Nothing, AbstractArray{Complex{Float64}, 5}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Complex{Float64}}}
                elseif imasdd[sel]["data_type"] == "CPX_6D"
                    h[item] = ":: Union{Nothing, AbstractArray{Complex{Float64}, 6}} = nothing"
                    convertsion_types = Union{convertsion_types,Array{Complex{Float64}}}

                else
                    throw(ArgumentError("$(sel) IMAS $(imasdd[sel]["data_type"]) has not been mapped to Julia data type"))
                end
            end
            if ! (item in keys(h))
                h[item] = Dict()
                push!(branches, path[1:k])
            end
            h = h[item]
        end
    end

    # generate Julia data structures
    is_structarray = false
    for branch in reverse(branches)
        h = ddict
        for item in branch
            h = h[item]
        end

        sep = "__"
        is_structarray = length(branch) > 0 && occursin("[", branch[end])
        struct_name_ = replace(join(branch, sep), "[:]" => "_")
        struct_name = replace(struct_name_, r"_$" => "")
        txt = String[]
        txt_parent = String[]
        inits = String[]
        for (item, info) in h
            # leaf
            if typeof(info) <: String
                push!(txt, "    var\"$(item)\" $(info)")
                push!(inits, "var\"$(item)\"=nothing")
            # branch
            else
                # top level
                if length(struct_name_) == 0
                    push!(txt, "    var\"$(item)\" :: Union{Nothing, $(item)} = $(item)()")
                    push!(txt_parent, "        obj.$(item)._parent = WeakRef(obj)")
                    push!(inits, "var\"$(item)\"=$(item)()")
                # arrays of structs
                elseif occursin("[:]", item)
                    item_ = replace(item, "[:]" => "_")
                    item = replace(item, "[:]" => "")
                    item_=item
                    push!(txt, "    var\"$(item)\" :: FDSvector{T} where {T<:$(struct_name_)$(sep)$(item_)} = FDSvector($(struct_name_)$(sep)$(item_)[])")
                    push!(txt_parent, "        obj.$(item)._parent = WeakRef(obj)")
                    push!(inits, "var\"$(item)\"=FDSvector($(struct_name_)$(sep)$(item_)[])")
                # structs
                else
                    push!(txt, "    var\"$(item)\" :: $(struct_name_)$(sep)$(item) = $(struct_name_)$(sep)$(item)()")
                    push!(txt_parent, "        obj.$(item)._parent = WeakRef(obj)")
                    push!(inits, "var\"$(item)\"=$(struct_name_)$(sep)$(item)()")
                end
            end
        end

        txt = join(txt, "\n")
        txt_parent = join(txt_parent, "\n")
        inits = join(inits, ", ")
        if length(struct_name) == 0
            struct_name = "dd"
        end
        txt_parent = """
    _parent :: WeakRef = WeakRef(nothing)
    function $(struct_name)($(inits), _parent=WeakRef(nothing))
        obj = new($(join(map(x -> split(x, "=")[1], split(inits, ", ")), ", ")), _parent)
$(txt_parent)
        return obj
    end
"""

        txt = """
Base.@kwdef mutable struct $(struct_name) <: FDS
$(txt)
$(rstrip(txt_parent))
end
"""

        push!(struct_commands, txt)

#       uncomment to debug data structure
#        println(txt)
#        eval(Meta.parse(txt))
    end

    return struct_commands, convertsion_types
end

struct_commands, convertsion_types = imas_julia_struct(desired_structure)

# Parse the Julia structs to make sure there are no issues
using ProgressMeter
ProgressMeter.ijulia_behavior(:clear)
p = Progress(length(struct_commands); desc="Compile IMAS structs", showspeed=true)
for txt in struct_commands
    ProgressMeter.next!(p)
    try
        eval(Meta.parse(txt))
    catch e
        println(txt)
        throw(e)
    end
end

open("$(dirname(@__FILE__))/dd.jl","w") do io
    println(io, "abstract type FDS end\n")
    println(io, "include(\"functionarrays.jl\")\n")
    println(io, "convertsion_types = $(string(convertsion_types))\n")
    println(io, join(struct_commands, "\n"))
end
