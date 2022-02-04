import PhysicalConstants.CODATA2018 as PCs

struct PhysicsConstants
    μ_0 :: Float64
    c_0 :: Float64
    ϵ_0 :: Float64
    k_B :: Float64
    e :: Float64
    m_e :: Float64
    m_p :: Float64
    m_n :: Float64
    atm :: Float64
    m_u :: Float64
end

Constants = PhysicsConstants(
    float(Float64, PCs.μ_0).val,
    float(Float64, PCs.c_0).val,
    float(Float64, PCs.ε_0).val,
    float(Float64, PCs.k_B).val,
    float(Float64, PCs.e).val,
    float(Float64, PCs.m_e).val,
    float(Float64, PCs.m_p).val,
    float(Float64, PCs.m_n).val,
    float(Float64, PCs.atm).val,
    float(Float64, PCs.m_u).val
    )
