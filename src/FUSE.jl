__precompile__()

module FUSE

include("dd.jl")

"""
    dict2fuse(h, f, ;path=[], verbose=false)

Populate FUSE structures based on data contained in dictionary

"""
function dict2fuse(h, f::FDS=dd() ;path=[], verbose=false)
    if verbose
        level = length(path)
    else
        level = -1
    end
    for (k, v) in h
        # Struct
        if typeof(v) <: Dict
            if level >= 0 println((" "^(level + 1)) * string(k)) end
            if typeof(f) <: FDS
                ff = eval(Meta.parse("FUSE.$k"))()
            else
                ff = getfield(f, Symbol(k))
            end
            dict2fuse(v, ff; path=vcat(path, [k]), verbose=verbose)

        # Array of struct
        elseif (typeof(v) <: Array) && (length(v) > 0) && (typeof(v[1]) <: Dict)
            ff = getfield(f, Symbol(k))
            if level >= 0 println((" "^(level + 1)) * string(k)) end
            for i in 1:length(v)
                if level >= 0 println((" "^(level + 2)) * string(i)) end
                dict2fuse(v[i], ff[i]; path = vcat(path, [k,"[$i]"]), verbose=verbose)
            end

        # Leaf
        else
            if level >= 0 println((" "^(level + 1)) * string(k)) end
            target_type = typeof(getfield(f, Symbol(k)))
            if target_type <: Array
                if ndims(target_type) == 2
                    v = hcat(v...)'
                end
            end
            setfield!(f, Symbol(k), convert(target_type, v))
        end
    end
end

end # module
