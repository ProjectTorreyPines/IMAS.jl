document[Symbol("Physics diagnostics")] = Symbol[]

"""
    extract_subsystem_from_channel_name(name::AbstractString)

Extract subsystem identifier from channel name.
Returns the subsystem string, or the full name if no pattern is recognized.

Handles various naming conventions:

  - Position-based: "TS_core_r+0_38" → "TS_core"
  - Position-based: "TS_tangential_r+1_25" → "TS_tangential"
  - Simple numbered: "TS_CORE_01" → "TS_CORE"
  - Generic pattern: "DIVERTOR_05" → "DIVERTOR"
"""
function extract_subsystem_from_channel_name(name::AbstractString)
    # First try pattern with position info
    ts_pattern = r"^(TS_[A-Za-z_]+?)_r[+-][0-9]+_[0-9]+$"
    m = match(ts_pattern, name)
    if m !== nothing
        return m.captures[1]
    end

    # Try general pattern: SUBSYSTEM_extras_NUMBER
    general_pattern = r"^([A-Za-z_]+?)(?:_[A-Za-z0-9+-]+)*_[0-9]+$"
    m = match(general_pattern, name)
    if m !== nothing
        subsys = m.captures[1]
        return rstrip(subsys, '_')
    end

    # Fallback: treat the entire name as subsystem
    return string(name)
end

@compat public extract_subsystem_from_channel_name
push!(document[Symbol("Physics diagnostics")], :extract_subsystem_from_channel_name)