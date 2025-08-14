using IMAS
using IMAS: cgs, mks

"""
    c_s(cp1d::IMAS.core_profiles__profiles_1d)

Sounds speed in [cm/s]
"""
function c_s(cp1d::IMAS.core_profiles__profiles_1d)
    Te = cp1d.electrons.temperature
    return c_s.(Te)
end
function c_s(Te::Real)
    return sqrt(cgs.k * Te / cgs.md)
end

"""
    rho_s(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

sound gyro radius in [cm]
"""
function rho_s(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    eqt1d = eqt.profiles_1d
    rho_cp = cp1d.grid.rho_tor_norm
    Te = cp1d.electrons.temperature
    bu_itp = IMAS.interp1d(eqt1d.rho_tor_norm, bunit(eqt1d))
    return rho_s.(rho_cp, Te, Ref(bu_itp))
end

function rho_s(rho::Real, Te::Real, bu_itp)
    return c_s(Te) / (cgs.e * abs(bu_itp(rho)) * cgs.T_to_Gauss) * (cgs.md * cgs.c)
end

"""
    r_min_core_profiles(eqt::IMAS.equilibrium__time_slice, rho_tor_norm::AbstractVector)

Geometric minor radius in [cm] evaluated
"""
function r_min_core_profiles(eqt1d::IMAS.equilibrium__time_slice___profiles_1d, rho_tor_norm::AbstractVector)
    return IMAS.interp1d(eqt1d.rho_tor_norm, cgs.m_to_cm * 0.5 * (eqt1d.r_outboard - eqt1d.r_inboard)).(rho_tor_norm)
end

"""
    bunit(eqt1d::IMAS.equilibrium__time_slice___profiles_1d)

Calculate bunit from equilibrium
"""
function bunit(eqt1d::IMAS.equilibrium__time_slice___profiles_1d)
    rmin = 0.5 .* (eqt1d.r_outboard .- eqt1d.r_inboard)
    phi = eqt1d.phi
    bunit = similar(phi)
    IMAS.gradient!(bunit, 2π .* rmin, phi)
    bunit ./= rmin
    return bunit
end

function bunit(eqt::IMAS.equilibrium__time_slice)
    return bunit(eqt.profiles_1d)
end