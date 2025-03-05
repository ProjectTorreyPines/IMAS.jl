document[Symbol("Physics constants")] = Symbol[]

import PhysicalConstants.CODATA2018 as PCs

"""
Named tuple with physics constants in mks (and eV):

    μ_0 = 1.25663706212e-6 [N A^-2]  # Vacuum permeability
    c = 2.99792458e8 [m s^-1]        # Speed of light in vacuum
    ϵ_0 = 8.8541878128e-12 [F m^-1]  # Vacuum permittivity
    k_B = 1.380649e-23 [J K^-1]      # Boltzmann constant
    e = 1.602176634e-19 [C]          # Elementary charge
    m_e = 9.1093837015e-31 [kg]      # Electron mass
    m_p = 1.67262192369e-27 [kg]     # Proton mass
    m_n = 1.67492749804e-27 [kg]     # Neutron mass
    m_d = 3.3435837768e-27 [kg]      # Deuteron mass
    atm = 101325.0 [Pa]              # Standard atmosphere
    m_u = 1.6605390666e-27 [kg]      # Atomic mass constant
    avog = 6.02214076e23 [mol^-1]    # Avogadro constant
    E_α = 3.518e6 [eV]               # Alpha particle energy
    E_n = 14.072e6 [eV]              # Neutron energy
"""
const mks = (
    μ_0=float(Float64, PCs.μ_0).val,
    c=float(Float64, PCs.c_0).val,
    ϵ_0=float(Float64, PCs.ε_0).val,
    k_B=float(Float64, PCs.k_B).val,
    e=float(Float64, PCs.e).val,
    m_e=float(Float64, PCs.m_e).val,
    m_p=float(Float64, PCs.m_p).val,
    m_n=float(Float64, PCs.m_n).val,
    m_d=3.3435837768e-27,
    atm=float(Float64, PCs.atm).val,
    m_u=float(Float64, PCs.m_u).val,
    avog=float(Float64, PCs.AvogadroConstant).val,
    E_α=3.518e6,
    E_n=14.072e6
)

@compat public mks
push!(document[Symbol("Physics constants")], :mks)

"""
Named tuple with physics constants used for interfacing with CGS codes:

    e=4.8032e-10 # stacoul
    k=1.6022e-12 # erg/eV
    c=2.9979e10 # cm/s
    me=9.1094e-28 # g
    mp=1.6726e-24 # g
    md=3.3435837768e-24 # g
    T_to_Gauss=1e4
    Erg_to_J=1e-7
    m_to_cm=1e2
    m²_to_cm²=1E4
    m³_to_cm³=1e6
"""
const cgs = (
    e=4.8032e-10, # stacoul
    k=1.6022e-12, # erg/eV
    c=2.9979e10, # cm/s
    me=9.1094e-28, # g
    mp=1.6726e-24, # g
    md=3.3435837768e-24, # g
    T_to_Gauss=1e4,
    Erg_to_J=1e-7,
    m_to_cm=1e2,
    m²_to_cm²=1E4,
    m³_to_cm³=1e6,
)

@compat public cgs
push!(document[Symbol("Physics constants")], :cgs)
