__precompile__()

module FUSE

const FUSE_verbose = false
const imas_version = "3_32_1"
const omas_imas_structure_folder = "/Users/meneghini/Coding/atom/omas/omas/imas_structures"

#= ==================================== =#
# FUSE data structure
#= ==================================== =#

desired_structure = split(chomp("""
core_profiles.profiles_1d[:].electrons.density
core_profiles.profiles_1d[:].electrons.density_fast
equilibrium.time_slice[:].global_quantities.ip
core_profiles.profiles_1d[:].grid.rho_tor_norm
"""))

#= ==================================== =#
# Implementing selected IMAS DD as Julia structures
#= ==================================== =#
using Memoize
using StructArrays:StructArray

using ProgressMeter
ProgressMeter.ijulia_behavior(:clear)

import JSON

@memoize function load_imasdd(ids;imas_version=imas_version)
    JSON.parsefile("$omas_imas_structure_folder/$imas_version/$ids.json")  # parse and transform data
end

function imas_julia_struct(desired_structure)
    struct_commands = []

    branches = Vector()
    push!(branches, "")

    # generate hierarchy of dicts
    p = Progress(length(desired_structure); showspeed=true)
    ddict = Dict()
    for sel in desired_structure
        # FUSE_verbose && println(sel)
        ProgressMeter.next!(p)# ; showvalues = [(:state,"hierarchy: $(sel)")])
        if endswith(sel, "error_upper") | endswith(sel, "error_lower") | endswith(sel, "error_index")
            continue
        end
        if (occursin("local", sel)) | (occursin("module", sel)) | (occursin("function", sel))
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
                    h[item] = ":: String = \"\""
                elseif imasdd[sel]["data_type"] == "STR_1D"
                    h[item] = ":: Array{String, 1} = String[]"

                elseif imasdd[sel]["data_type"] in ["INT_0D", "INT_TYPE"]
                    h[item] = ":: Int32 = 0"
                elseif imasdd[sel]["data_type"] == "INT_1D"
                    h[item] = ":: Array{Int32, 1} = zeros(Int32,(0))"
                elseif imasdd[sel]["data_type"] == "INT_2D"
                    h[item] = ":: Array{Int32, 2} = zeros(Int32,(0, 0))"
                elseif imasdd[sel]["data_type"] == "INT_3D"
                    h[item] = ":: Array{Int32, 3} = zeros(Int32,(0, 0, 0))"

                elseif imasdd[sel]["data_type"] == "FLT_0D"
                    h[item] = ":: Float64 = 0.0"
                elseif imasdd[sel]["data_type"] in ["FLT_1D", "FLT_1D_TYPE"]
                    h[item] = ":: Array{Float64, 1} = zeros(Float64,(0))"
                elseif imasdd[sel]["data_type"] == "FLT_2D"
                    h[item] = ":: Array{Float64, 2} = zeros(Float64,(0,0))"
                elseif imasdd[sel]["data_type"] == "FLT_3D"
                    h[item] = ":: Array{Float64, 3} = zeros(Float64,(0,0,0))"
                elseif imasdd[sel]["data_type"] == "FLT_4D"
                    h[item] = ":: Array{Float64, 4} = zeros(Float64,(0,0,0,0))"
                elseif imasdd[sel]["data_type"] == "FLT_5D"
                    h[item] = ":: Array{Float64, 5} = zeros(Float64,(0,0,0,0,0))"
                elseif imasdd[sel]["data_type"] == "FLT_6D"
                    h[item] = ":: Array{Float64, 6} = zeros(Float64,(0,0,0,0,0,0))"

                elseif imasdd[sel]["data_type"] == "CPX_0D"
                    h[item] = ":: Complex{Float64} = 0.0+0.0im"
                elseif imasdd[sel]["data_type"] == "CPX_1D"
                    h[item] = ":: Array{Complex{Float64}, 1} = zeros(Complex{Float64},(0))"
                elseif imasdd[sel]["data_type"] == "CPX_2D"
                    h[item] = ":: Array{Complex{Float64}, 2} = zeros(Complex{Float64},(0,0))"
                elseif imasdd[sel]["data_type"] == "CPX_3D"
                    h[item] = ":: Array{Complex{Float64}, 3} = zeros(Complex{Float64},(0,0,0))"
                elseif imasdd[sel]["data_type"] == "CPX_4D"
                    h[item] = ":: Array{Complex{Float64}, 4} = zeros(Complex{Float64},(0,0,0,0))"
                elseif imasdd[sel]["data_type"] == "CPX_5D"
                    h[item] = ":: Array{Complex{Float64}, 5} = zeros(Complex{Float64},(0,0,0,0,0))"
                elseif imasdd[sel]["data_type"] == "CPX_6D"
                    h[item] = ":: Array{Complex{Float64}, 6} = zeros(Complex{Float64},(0,0,0,0,0,0))"

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
    # FUSE_verbose && JSON.print(ddict,4)

    # generate Julia data structures
    p = Progress(length(branches); showspeed=true)
    for branch in reverse(branches)
        ProgressMeter.next!(p)# ; showvalues = [(:branch,branch)])
        h = ddict
        for item in branch
            h = h[item]
        end

        sep = "__"
        struct_name = replace(join(branch, sep), "[:]" => "")
        txt = []
        for (item, info) in h
            if typeof(info) <: String
                push!(txt, "    $(item) $(info)")
            else
                if length(struct_name) == 0
                    push!(txt, "    $(item) :: $(item) = $(item)()")
                elseif occursin("[:]", item)
                    item = replace(item, "[:]" => "")
                    push!(txt, "    $(item) :: StructArray{$(struct_name)$(sep)$(item)} = StructArray($(struct_name)$(sep)$(item)() for k in 1:1)")
                else
                    push!(txt, "    $(item) :: $(struct_name)$(sep)$(item) = $(struct_name)$(sep)$(item)()")
                end
            end
        end
        
        txt = join(txt, "\n")
        if length(struct_name) == 0
            txt = """
Base.@kwdef mutable struct dd
$(txt)
end
"""
        else
            txt = """
Base.@kwdef mutable struct $(struct_name)
$(txt)
end
"""
        end
        FUSE_verbose && println(txt)
        push!(struct_commands, txt)
    end

    return struct_commands
end

desired_structure = []
for filename in readdir("$omas_imas_structure_folder/$imas_version")
    if startswith(filename, "_") || ! endswith(filename, ".json")
        continue
    end
    filename = replace(filename, ".json" => "")
    for item in keys(load_imasdd(filename))
        push!(desired_structure, item)
    end
end

struct_commands = imas_julia_struct(desired_structure)

for txt in struct_commands
    try
        eval(Meta.parse(txt))
    catch e
        println(txt)
        throw(e)
    end
end

"""
Resize StructArray to contain n elements.
If n is smaller than the current collection length, the first n elements will be retained.
If n is larger, the new elements are not guaranteed to be initialized.
"""
function resize!(collection::StructArray, n::Int)
    if n > length(collection)
        for k in length(collection):n - 1
            push!(collection, eltype(collection)())
            # println("push $(length(collection))")
        end
    elseif n < length(collection)
        for k in n:length(collection) - 1
            pop!(collection)
            # println("pop $(length(collection))")
        end
    end
    return collection
end

end # module
