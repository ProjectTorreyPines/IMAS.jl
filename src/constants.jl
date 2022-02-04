import PhysicalConstants.CODATA2018 as PCs

@kwdef struct PhysicsConstants
    μ_0 :: Float64 = float(Float64, PCs.μ_0).val
    c_0 :: Float64 = float(Float64, PCs.c_0).val
    ϵ_0 :: Float64 = float(Float64, PCs.ε_0).val
    k_B :: Float64 = float(Float64, PCs.k_B).val
    e :: Float64 = float(Float64, PCs.e).val
    m_e :: Float64 =  float(Float64, PCs.m_e).val
    m_p :: Float64 =  float(Float64, PCs.m_p).val
    m_n :: Float64 =  float(Float64, PCs.m_n).val
    atm :: Float64 = float(Float64, PCs.atm).val
    m_u :: Float64 =  float(Float64, PCs.m_u).val
end

Constants = PhysicsConstants()
