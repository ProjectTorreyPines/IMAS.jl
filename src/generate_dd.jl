using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using ProgressMeter
ProgressMeter.ijulia_behavior(:clear)

const imas_version = "3_33_0"
if !("OMAS_ROOT" in keys(ENV))
    error("Environmental variable OMAS_ROOT is not set")
end
const omas_imas_structure_folder = joinpath(ENV["OMAS_ROOT"], "omas", "imas_structures")
run(`sh -c "rm -rf $(dirname(dirname(@__FILE__)))/data_structures"`)
run(`sh -c "mkdir $(dirname(dirname(@__FILE__)))/data_structures"`)
run(`sh -c "cp -rf $(omas_imas_structure_folder)/$(imas_version)/*.json $(dirname(dirname(@__FILE__)))/data_structures"`)

#= ==================================== =#
# Header
#= ==================================== =#
include("functionarrays.jl")
assign_expressions = x -> x

#= ==================================== =#
# IMAS data structure
#= ==================================== =#
ids_names = imas_dd_ids_names()
ids_names = ["core_profiles", "core_sources", "dataset_description", "equilibrium", "pf_active", "summary", "tf", "wall","nbi", "ec_launchers", "lh_antennas", "build"]

p = Progress(length(ids_names); desc = "Parse JSON structs ", showspeed = true)
desired_structure = String[]
for ids_name in ids_names
    ProgressMeter.next!(p)
    for item in sort(collect(keys(imas_dd_ids(ids_name))))
        push!(desired_structure, item)
    end
end

#= ================================================= =#
#  Implementing selected IMAS DD as Julia structures
#= ================================================= =#

const type_translator = Dict("STR" => String, "INT" => Integer, "FLT" => Real, "CPX" => Complex{Real})

"""
    imas_julia_struct(desired_structure::Vector{String})

Translate the IMAS `desired_structure` entries into Julia structs
Note that `desired_structure` are fully qualified IMAS locations
"""
function imas_julia_struct(desired_structure::Vector{String})
    struct_commands = String[]

    branches = Vector()
    timedep_structures = Vector()
    push!(branches, "")

    # generate hierarchy of dicts
    p = Progress(length(desired_structure); desc = "Generate hierarchy of dicts ", showspeed = true)
    ddict = Dict()
    for sel in sort(desired_structure)
        ProgressMeter.next!(p)

        # shed some weight (1): ignore error fields
        if endswith(sel, "error_upper") | endswith(sel, "error_lower") | endswith(sel, "error_index")
            continue
        end

        # split IMAS path in 
        path = split(sel, ".")

        # shed some weight (2): no ggd, code, ids_properties
        if ("ggd" in path) || ("ggd[:]" in path) || ("grids_ggd" in path) || ("grids_ggd[:]" in path) || ("code" in path) || ("ids_properties" in path)
            continue
        end

        # load imas data structure
        ddids = imas_dd_ids(path[1])

        # shed some weight (3): no obsolescent
        if ("lifecycle_status" in keys(ddids[sel])) && (ddids[sel]["lifecycle_status"] == "obsolescent")
            continue
        end

        # detect if this is a time dependendent array of structure
        if ("full_path" in keys(ddids[sel])) && (endswith(ddids[sel]["full_path"], "(itime)"))
            push!(timedep_structures, sel * "[:]")
        end

        h = ddict
        for (k, item) in enumerate(path)
            if (k > 1) & (k == length(path))
                if ddids[sel]["data_type"] in ["STRUCTURE", "STRUCT_ARRAY"]
                    continue
                end

                # find data type and dimension
                (tp, dim) = split(ddids[sel]["data_type"], "_")
                dim = parse(Int, replace(dim, "D" => ""))

                # translate from IMAS data type to Julia types
                if tp in keys(type_translator)
                    if dim == 0
                        h[item] = ":: Union{Missing, $(type_translator[tp]), Function}"
                    elseif dim == 1
                        h[item] = ":: Union{Missing, AbstractArray{T, $dim} where T<:$(type_translator[tp]), AbstractRange{T} where T<:$(type_translator[tp]), Function}"
                    else
                        h[item] = ":: Union{Missing, AbstractArray{T, $dim} where T<:$(type_translator[tp]), Function}"
                    end
                else
                    throw(ArgumentError("$(sel) IMAS $(ddids[sel]["data_type"]) has not been mapped to Julia data type"))
                end
            end
            if !(item in keys(h))
                h[item] = Dict()
                push!(branches, path[1:k])
            end
            h = h[item]
        end
    end

    # generate source code of Julia structures
    is_structarray = false
    p = Progress(length(branches); desc = "Generate source code of Julia structures ", showspeed = true)
    for branch in reverse(branches)
        ProgressMeter.next!(p)

        h = ddict
        for item in branch
            h = h[item]
        end

        is_structarray = length(branch) > 0 && occursin("[", branch[end])
        if is_structarray
            if join(branch, ".") in timedep_structures
                struct_type = "IDSvectorTimeElement"
            else
                struct_type = "IDSvectorStaticElement"
            end
        else
            struct_type = "IDS"
        end

        sep = "__"
        struct_name_ = replace(join(branch, sep), "[:]" => "_")
        struct_name = replace(struct_name_, r"_$" => "")
        txt = String[]
        txt_parent = String[]
        inits = String[]
        for item in sort(collect(keys(h)))
            info = h[item]
            # leaf
            if typeof(info) <: String
                push!(txt, "    var\"$(item)\" $(info)")
                push!(inits, "var\"$(item)\"=missing")
                # branch
            else
                # top level
                if length(struct_name_) == 0
                    push!(txt, "    var\"$(item)\" :: Union{Missing, $(item)}")
                    push!(inits, "var\"$(item)\"=$(item)()")
                    # arrays of structs
                elseif occursin("[:]", item)
                    item_ = replace(item, "[:]" => "_")
                    item = replace(item, "[:]" => "")
                    item_ = item
                    push!(txt, "    var\"$(item)\" :: IDSvector{T} where {T<:$(struct_name_)$(sep)$(item_)}")
                    push!(inits, "var\"$(item)\"=IDSvector($(struct_name_)$(sep)$(item_)[])")
                    # structs
                else
                    push!(txt, "    var\"$(item)\" :: $(struct_name_)$(sep)$(item)")
                    push!(inits, "var\"$(item)\"=$(struct_name_)$(sep)$(item)()")
                end
                push!(txt_parent, "        setfield!(ids.$(item), :_parent, WeakRef(ids))")
            end
        end

        txt_parent = join(txt_parent, "\n")
        if length(txt_parent) > 0
            txt_parent = "\n$(txt_parent)"
        end
        if length(struct_name) == 0
            struct_name = "dd"
            global_time_varname = "global_time"
            push!(txt, "    var\"$(global_time_varname)\" :: Real")
            push!(inits, "var\"$(global_time_varname)\"=0.0")
        end
        inits = join(inits, ", ")
        txt_parent = """
    _parent :: WeakRef
    function $(struct_name)($(inits), _parent=WeakRef(missing))
        ids = new($(join(map(x -> split(x, "=")[1], split(inits, ", ")), ", ")), _parent)
        assign_expressions(ids)$(txt_parent)
        return ids
    end
"""

        txt = join(txt, "\n")
        txt = """
mutable struct $(struct_name) <: $(struct_type)
$(txt)
$(rstrip(txt_parent))
end
"""

        push!(struct_commands, txt)

        # #uncomment to debug data structure
        # println(txt)
        # eval(Meta.parse(txt))
    end

    return struct_commands
end

struct_commands = imas_julia_struct(desired_structure)

# Parse the Julia structs to make sure there are no issues
p = Progress(length(struct_commands); desc = "Compile IMAS.jl structs ", showspeed = true)
for txt in struct_commands
    ProgressMeter.next!(p)
    try
        eval(Meta.parse(txt))
    catch e
        println(txt)
        throw(e)
    end
end

open("$(dirname(@__FILE__))/dd.jl", "w") do io
    println(io, "include(\"functionarrays.jl\")\n")
    println(io, join(struct_commands, "\n"))
end
