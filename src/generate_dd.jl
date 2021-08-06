using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

const imas_version = "3_33_0"
const omas_imas_structure_folder = "/Users/meneghini/Coding/atom/omas/omas/imas_structures"

import JSON
using Memoize

"""
Function used to read the IMAS data structures in the OMAS JSON format
"""
@memoize function load_imasdd(ids;imas_version=imas_version)
    JSON.parsefile("$omas_imas_structure_folder/$imas_version/$ids.json")  # parse and transform data
end

#= ==================================== =#
# FUSE data structure
#= ==================================== =#

desired_structure = split(chomp("""
core_profiles.profiles_1d[:].grid.rho_tor_norm
core_profiles.profiles_1d[:].electrons.density
core_profiles.profiles_1d[:].electrons.density_fast
equilibrium.time_slice[:].global_quantities.ip
"""))

filenames = readdir("$omas_imas_structure_folder/$imas_version")
desired_structure = []
for filename in filenames
    if startswith(filename, "_") || ! endswith(filename, ".json")
        continue
    end
    filename = replace(filename, ".json" => "")
    for item in keys(load_imasdd(filename))
        push!(desired_structure, item)
    end
end

#= ==================================== =#
# Implementing selected IMAS DD as Julia structures
#= ==================================== =#

abstract type FDS end

"""
Function used to translate the desired IMAS data structures into Julia structs
"""
function imas_julia_struct(desired_structure)
    struct_commands = []

    branches = Vector()
    push!(branches, "")

    # generate hierarchy of dicts
    ddict = Dict()
    for sel in desired_structure
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
                elseif imasdd[sel]["data_type"] == "STR_1D"
                    h[item] = ":: Union{Nothing, Array{String, 1}} = nothing"

                elseif imasdd[sel]["data_type"] in ["INT_0D", "INT_TYPE"]
                    h[item] = ":: Union{Nothing, Int32} = nothing"
                elseif imasdd[sel]["data_type"] == "INT_1D"
                    h[item] = ":: Union{Nothing, Array{Int32, 1}} = nothing"
                elseif imasdd[sel]["data_type"] == "INT_2D"
                    h[item] = ":: Union{Nothing, Array{Int32, 2}} = nothing"
                elseif imasdd[sel]["data_type"] == "INT_3D"
                    h[item] = ":: Union{Nothing, Array{Int32, 3}} = nothing"
                elseif imasdd[sel]["data_type"] == "INT_4D"
                    h[item] = ":: Union{Nothing, Array{Int32, 4}} = nothing"
                elseif imasdd[sel]["data_type"] == "INT_5D"
                    h[item] = ":: Union{Nothing, Array{Int32, 5}} = nothing"
                elseif imasdd[sel]["data_type"] == "INT_6D"
                    h[item] = ":: Union{Nothing, Array{Int32, 6}} = nothing"

                elseif imasdd[sel]["data_type"] == "FLT_0D"
                    h[item] = ":: Union{Nothing, Float64} = nothing"
                elseif imasdd[sel]["data_type"] in ["FLT_1D", "FLT_1D_TYPE"]
                    h[item] = ":: Union{Nothing, Array{Float64, 1}} = nothing"
                elseif imasdd[sel]["data_type"] == "FLT_2D"
                    h[item] = ":: Union{Nothing, Array{Float64, 2}} = nothing"
                elseif imasdd[sel]["data_type"] == "FLT_3D"
                    h[item] = ":: Union{Nothing, Array{Float64, 3}} = nothing"
                elseif imasdd[sel]["data_type"] == "FLT_4D"
                    h[item] = ":: Union{Nothing, Array{Float64, 4}} = nothing"
                elseif imasdd[sel]["data_type"] == "FLT_5D"
                    h[item] = ":: Union{Nothing, Array{Float64, 5}} = nothing"
                elseif imasdd[sel]["data_type"] == "FLT_6D"
                    h[item] = ":: Union{Nothing, Array{Float64, 6}} = nothing"

                elseif imasdd[sel]["data_type"] == "CPX_0D"
                    h[item] = ":: Union{Nothing, Complex{Float64} = nothing"
                elseif imasdd[sel]["data_type"] == "CPX_1D"
                    h[item] = ":: Union{Nothing, Array{Complex{Float64}, 1}} = nothing"
                elseif imasdd[sel]["data_type"] == "CPX_2D"
                    h[item] = ":: Union{Nothing, Array{Complex{Float64}, 2}} = nothing"
                elseif imasdd[sel]["data_type"] == "CPX_3D"
                    h[item] = ":: Union{Nothing, Array{Complex{Float64}, 3}} = nothing"
                elseif imasdd[sel]["data_type"] == "CPX_4D"
                    h[item] = ":: Union{Nothing, Array{Complex{Float64}, 4}} = nothing"
                elseif imasdd[sel]["data_type"] == "CPX_5D"
                    h[item] = ":: Union{Nothing, Array{Complex{Float64}, 5}} = nothing"
                elseif imasdd[sel]["data_type"] == "CPX_6D"
                    h[item] = ":: Union{Nothing, Array{Complex{Float64}, 6}} = nothing"

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
        is_structarray = length(branch)>0 && occursin("[",branch[end])
        struct_name = replace(join(branch, sep), "[:]" => "")
        txt = []
        for (item, info) in h
            # leaf
            if typeof(info) <: String
                push!(txt, "    var\"$(item)\" $(info)")
            # branch
            else
                # top level
                if length(struct_name) == 0
                    push!(txt, "    var\"$(item)\" :: Union{Nothing, $(item)} = $(item)()")
                    #push!(txt, "    var\"$(item)\" :: Union{Nothing, $(item)} = nothing")
                # arrays of structs
                elseif occursin("[:]", item)
                    item = replace(item, "[:]" => "")
                    push!(txt, "    var\"$(item)\" :: Vector{$(struct_name)$(sep)$(item)} = [$(struct_name)$(sep)$(item)()]")
                # structs
                else
                    push!(txt, "    var\"$(item)\" :: $(struct_name)$(sep)$(item) = $(struct_name)$(sep)$(item)()")
                end
            end
        end
        
        txt = join(txt, "\n")
        if length(struct_name) == 0
            txt = """
Base.@kwdef mutable struct dd <: FDS
$(txt)
end
"""
        elseif is_structarray
            txt = """
Base.@kwdef mutable struct $(struct_name) <: FDS
$(txt)
end
"""
        else
            txt = """
Base.@kwdef mutable struct $(struct_name) <: FDS
$(txt)
end
"""
        end
        push!(struct_commands, txt)
    end

    return struct_commands
end


const struct_commands = imas_julia_struct(desired_structure)

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
    println(io, join(struct_commands, "\n"))
end
