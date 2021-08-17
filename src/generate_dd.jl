using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

const imas_version = "3_33_0"
ENV["OMAS_ROOT"] = "/Users/meneghini/Coding/atom/omas"
const omas_imas_structure_folder = joinpath(ENV["OMAS_ROOT"], "omas", "imas_structures")
run(`sh -c "cp -rf $(omas_imas_structure_folder)/$(imas_version)/*.json $(dirname(dirname(@__FILE__)))/data_structures"`)

#= ==================================== =#
# Header
#= ==================================== =#
include("functionarrays.jl")

#= ==================================== =#
# FUSE data structure
#= ==================================== =#
filenames = readdir(joinpath(dirname(dirname(@__FILE__)), "data_structures"))
desired_structure = []
for filename in filenames
    if startswith(filename, "_") || ! endswith(filename, ".json")
        continue
    end
    if ! (filename in ["equilibrium.json", "core_profiles.json", "wall.json", "dataset_description.json"])
        continue
    end
    filename = replace(filename, ".json" => "")
    for item in keys(imas_load_dd(filename))
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
        
        imasdd = imas_load_dd(path[1])

        h = ddict
        for (k, item) in enumerate(path)
            if (k > 1) & (k == length(path))
                if imasdd[sel]["data_type"] in ["STRUCTURE", "STRUCT_ARRAY"]
                    continue
                end

                if imasdd[sel]["data_type"] == "INT_TYPE"
                    imasdd[sel]["data_type"] = "INT_0D"
                elseif imasdd[sel]["data_type"] == "FLT_1D_TYPE"
                    imasdd[sel]["data_type"] = "FLT_1D"
                end

                (tp, dim) = split(imasdd[sel]["data_type"], "_")
                dim = parse(Int, replace(dim, "D" => ""))

                if (tp == "STR") & (dim == 0)
                    h[item] = ":: Union{Missing, String} = missing"
                    convertsion_types = Union{convertsion_types,String}
                elseif tp == "STR"
                    h[item] = ":: Union{Missing, AbstractArray{String, $dim}} = missing"
                    convertsion_types = Union{convertsion_types,Array{String}}

                elseif (tp == "INT") & (dim == 0)
                    h[item] = ":: Union{Missing, Int} = missing"
                    convertsion_types = Union{convertsion_types,Int}
                elseif tp == "INT"
                    h[item] = ":: Union{Missing, AbstractArray{Int, $dim}} = missing"
                    convertsion_types = Union{convertsion_types,Array{Int}}

                elseif (tp == "FLT") & (dim == 0)
                    h[item] = ":: Union{Missing, Float64} = missing"
                    convertsion_types = Union{convertsion_types,Float64}
                elseif tp == "FLT"
                    if false
                        h[item] = ":: Union{Missing, AbstractFunctionArray{Float64,$dim}} = missing"
                        convertsion_types = Union{convertsion_types,AbstractFunctionArray{Float64}}
                    else
                        h[item] = ":: Union{Missing, AbstractArray{Float64, $dim}} = missing"
                        convertsion_types = Union{convertsion_types,Array{Float64}}
                    end

                elseif (tp == "CPX") & (dim == 0)
                    h[item] = ":: Union{Missing, Complex{Float64} = missing"
                    convertsion_types = Union{convertsion_types,Complex{Float64}}
                elseif tp == "CPX"
                    h[item] = ":: Union{Missing, AbstractArray{Complex{Float64}, $dim}} = missing"
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
                push!(inits, "var\"$(item)\"=missing")
            # branch
            else
                # top level
                if length(struct_name_) == 0
                    push!(txt, "    var\"$(item)\" :: Union{Missing, $(item)} = $(item)()")
                    push!(txt_parent, "        obj.$(item)._parent = WeakRef(obj)")
                    push!(inits, "var\"$(item)\"=$(item)()")
                # arrays of structs
                elseif occursin("[:]", item)
                    item_ = replace(item, "[:]" => "_")
                    item = replace(item, "[:]" => "")
                    item_ = item
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
    _parent :: WeakRef = WeakRef(missing)
    function $(struct_name)($(inits), _parent=WeakRef(missing))
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
