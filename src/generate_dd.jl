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
desired_structure = String[]
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

const type_translator = Dict("STR" => String,
                             "INT" => Int64,
                             "FLT" => Float64,
                             "CPX" => Complex{Float64})

"""
    imas_julia_struct(desired_structure::Vector{String})

Translate the IMAS `desired_structure` entries into Julia structs
Note that `desired_structure` are fully qualified IMAS locations
"""
function imas_julia_struct(desired_structure::Vector{String})
    convertsion_types = Union{}
    struct_commands = String[]

    branches = Vector()
    push!(branches, "")

    # generate hierarchy of dicts
    ddict = Dict()
    for sel in sort(desired_structure)
        # ignore error fields
        if endswith(sel, "error_upper") | endswith(sel, "error_lower") | endswith(sel, "error_index")
            continue
        end

        # split IMAS path in 
        path = split(sel, ".")
        
        # load imas data structure
        imasdd = imas_load_dd(path[1])

        h = ddict
        for (k, item) in enumerate(path)
            if (k > 1) & (k == length(path))
                if imasdd[sel]["data_type"] in ["STRUCTURE", "STRUCT_ARRAY"]
                    continue
                end

                # fix some issues in IMAS DD type definitions
                if imasdd[sel]["data_type"] == "INT_TYPE"
                    imasdd[sel]["data_type"] = "INT_0D"
                elseif imasdd[sel]["data_type"] == "FLT_1D_TYPE"
                    imasdd[sel]["data_type"] = "FLT_1D"
                end

                # find data type and dimension
                (tp, dim) = split(imasdd[sel]["data_type"], "_")
                dim = parse(Int, replace(dim, "D" => ""))

                # translate from IMAS data type to Julia types
                if tp in keys(type_translator)
                    if dim == 0
                        h[item] = ":: Union{Missing, $(type_translator[tp])} = missing"
                        convertsion_types = Union{convertsion_types,type_translator[tp]}
                    elseif dim == 1
                        h[item] = ":: Union{Missing, AbstractFDVector{$(type_translator[tp])}} = missing"
                        convertsion_types = Union{convertsion_types,AbstractFDVector{type_translator[tp]}}
                    else
                        h[item] = ":: Union{Missing, AbstractArray{$(type_translator[tp]), $dim}} = missing"
                        convertsion_types = Union{convertsion_types,Array{type_translator[tp]}}
                    end
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

    # generate source code of Julia structures
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
                    push!(inits, "var\"$(item)\"=$(item)()")
                # arrays of structs
                elseif occursin("[:]", item)
                    item_ = replace(item, "[:]" => "_")
                    item = replace(item, "[:]" => "")
                    item_ = item
                    push!(txt, "    var\"$(item)\" :: FDSvector{T} where {T<:$(struct_name_)$(sep)$(item_)} = FDSvector($(struct_name_)$(sep)$(item_)[])")
                    push!(inits, "var\"$(item)\"=FDSvector($(struct_name_)$(sep)$(item_)[])")
                # structs
                else
                    push!(txt, "    var\"$(item)\" :: $(struct_name_)$(sep)$(item) = $(struct_name_)$(sep)$(item)()")
                    push!(inits, "var\"$(item)\"=$(struct_name_)$(sep)$(item)()")
                end
                push!(txt_parent, "        setfield!(obj.$(item), :_parent, WeakRef(obj))")
            end
        end

        txt = join(txt, "\n")
        txt_parent = join(txt_parent, "\n")
        inits = join(inits, ", ")
        if length(struct_name) == 0
            struct_name = "dd"
        end
        if length(txt_parent) > 0
            txt_parent = """
    _parent :: WeakRef = WeakRef(missing)
    function $(struct_name)($(inits), _parent=WeakRef(missing))
        obj = new($(join(map(x -> split(x, "=")[1], split(inits, ", ")), ", ")), _parent)
$(txt_parent)
        return obj
    end
"""
        else
            txt_parent = """
    _parent :: WeakRef = WeakRef(missing)
"""
        end
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
    println(io, "include(\"functionarrays.jl\")\n")
    println(io, "convertsion_types = $(string(convertsion_types))\n")
    println(io, join(struct_commands, "\n"))
end
