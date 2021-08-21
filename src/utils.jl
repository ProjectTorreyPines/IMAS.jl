"""
    common_base_string(s1::String, s2::String)::Vector{String}

given two strings it returns a tuple of 3 strings that is the common initial part, and then the remaining parts
"""
function common_base_string(s1::String, s2::String)::Tuple{String, String, String}
    index = nothing
    for k in 1:min(length(s1), length(s2))
        sub = SubString(s2, 1, k)
        if startswith(s1, sub)
            index = k
        end
    end
    if index === nothing
        return "", s1, s2
    else
        return string(SubString(s1, 1, index)), string(SubString(s1, index + 1, length(s1))), string(SubString(s2, index + 1, length(s2)))
    end
end

function iscallable(f)
    return !isempty(methods(f))
end