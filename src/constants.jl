import PhysicalConstants.CODATA2018 as PCs

"""
Named tuple with physics constants:

    μ_0 = 1.25663706212e-6
    c = 2.99792458e8
    ϵ_0 = 8.8541878128e-12
    k_B = 1.380649e-23
    e = 1.602176634e-19
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    m_n = 1.67492749804e-27
    atm = 101325.0
    m_u = 1.6605390666e-27
"""
const constants = (
    μ_0 = convert(Float64, PCs.μ_0).val,
    c = convert(Float64, PCs.c_0).val,
    ϵ_0 = convert(Float64, PCs.ε_0).val,
    k_B = convert(Float64, PCs.k_B).val,
    e = convert(Float64, PCs.e).val,
    m_e = convert(Float64, PCs.m_e).val,
    m_p = convert(Float64, PCs.m_p).val,
    m_n = convert(Float64, PCs.m_n).val,
    atm = convert(Float64, PCs.atm).val,
    m_u = convert(Float64, PCs.m_u).val
)

function index_2_name(ids::Union{IDS, IDSvector}, field::Symbol)
    return index_2_name(ids, Val(field))
end

function index_2_name(ids::IDSvector{core_transport__model}, field::Val{:identifier})
    return Dict(0 => :unspecified,
    1 => :combined, 2 => :transport_solver,
    3 => :background, 4 => :database, 5 => :neoclassical,
    6 => :anomalous, 19 => :mhd, 20 => :ntm,
    21 => :sawteeth, 22 => :elm_continuous, 23 => :elm_resolved,
    24 => :pedestal, 25 => :unknown)
end

function name_2_index(ids::Union{IDS, IDSvector}, field::Any)
    dict = index_2_name(ids,field)
    return Dict(v=>k for (k,v) in dict)
end
