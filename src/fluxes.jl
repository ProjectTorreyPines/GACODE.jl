using IMAS
using IMAS: cgs, mks
using IMASutils: argmin_abs
using Printf

"""
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
    PARTICLE_FLUX_e::T
    PARTICLE_FLUX_i::Vector{T}
    STRESS_TOR_i::T

Structure used to store fluxes information from GA code in Gyrobohm units
"""
struct FluxSolution{T<:Real}
    ENERGY_FLUX_e::T
    ENERGY_FLUX_i::T
    PARTICLE_FLUX_e::T
    PARTICLE_FLUX_i::Vector{T}
    STRESS_TOR_i::T
end

function Base.show(io::IO, ::MIME"text/plain", sol::FluxSolution)
    txt = """
    Qe =  $(sol.ENERGY_FLUX_e)
    Qi =  $(sol.ENERGY_FLUX_i)
    Γe =  $(sol.PARTICLE_FLUX_e)
    Γi = [$(join( map(string, sol.PARTICLE_FLUX_i),", "))]
    Πi =  $(sol.STRESS_TOR_i)
    """
    return print(io, txt)
end

"""
    gyrobohm_energy_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Gyrobohm energy flux
"""
function gyrobohm_energy_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature
    rhos = rho_s(cp1d, eqt)
    a = eqt.boundary.minor_radius
    return gyrobohm_energy_flux.(ne, Te, rhos, a)
end

function gyrobohm_energy_flux(ne::Real, Te::Real, rhos::Real, a::Real)
    return ne / cgs.m³_to_cm³ * cgs.k * Te * c_s(Te) *
           (rhos / (a * cgs.m_to_cm))^2 * cgs.Erg_to_J * cgs.m²_to_cm²
end

"""
    gyrobohm_particle_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Gyrobohm particle flux
"""
function gyrobohm_particle_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature
    rhos = rho_s(cp1d, eqt)
    a = eqt.boundary.minor_radius
    return gyrobohm_particle_flux.(ne, Te, rhos, a)
end

function gyrobohm_particle_flux(ne::Real, Te::Real, rhos::Real, a::Real)
    return gyrobohm_energy_flux(ne, Te, rhos, a) / (mks.e * Te)
end

"""
    gyrobohm_momentum_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)

Gyrobohm momentum flux
"""
function gyrobohm_momentum_flux(cp1d::IMAS.core_profiles__profiles_1d, eqt::IMAS.equilibrium__time_slice)
    ne = cp1d.electrons.density_thermal
    Te = cp1d.electrons.temperature
    rhos = rho_s(cp1d, eqt)
    a = eqt.boundary.minor_radius
    return gyrobohm_momentum_flux.(ne, Te, rhos, a)
end

function gyrobohm_momentum_flux(ne::Real, Te::Real, rhos::Real, a::Real)
    return ne / cgs.m³_to_cm³ * cgs.k * Te *
           rhos^2 / (a * cgs.m_to_cm) * cgs.Erg_to_J * cgs.m²_to_cm²
end

"""
    volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)

Correction to account for transformation from Miller r grid in GA code equilibrium to Psi grid in FUSE equilibrium
"""
function volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)
    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    return IMAS.gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
end

"""
    flux_gacode_to_imas(
        flux_types::Tuple{Vararg{Symbol}},
        flux_solutions::Vector{FluxSolution{T}},
        m1d::IMAS.core_transport__model___profiles_1d,
        eqt::IMAS.equilibrium__time_slice,
        cp1d::IMAS.core_profiles__profiles_1d
    )

Normalizes specified transport fluxes output by GA code via gyrobohm normalization and Miller volume correction
"""
function flux_gacode_to_imas(
    flux_types::Tuple{Vararg{Symbol}},
    flux_solutions::Vector{FluxSolution{T}},
    m1d::IMAS.core_transport__model___profiles_1d{T},
    eqt::IMAS.equilibrium__time_slice{T},
    cp1d::IMAS.core_profiles__profiles_1d{T}) where {T<:Real}

    rho_cp = cp1d.grid.rho_tor_norm
    rho_eq_idxs = [argmin_abs(eqt.profiles_1d.rho_tor_norm, rho) for rho in m1d.grid_flux.rho_tor_norm]
    rho_cp_idxs = [argmin_abs(rho_cp, rho) for rho in m1d.grid_flux.rho_tor_norm]

    @views vprime_miller = volume_prime_miller_correction(eqt)[rho_eq_idxs]
    @views ne = cp1d.electrons.density_thermal[rho_cp_idxs]
    @views Te = cp1d.electrons.temperature[rho_cp_idxs]
    bu_itp = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, bunit(eqt.profiles_1d))
    @views rhos = rho_s.(rho_cp[rho_cp_idxs], Te, Ref(bu_itp))
    a = eqt.boundary.minor_radius

    if :ion_energy_flux in flux_types
        m1d.total_ion_energy.flux = @. gyrobohm_energy_flux(ne, Te, rhos, a) * (f.ENERGY_FLUX_i for f in flux_solutions) * vprime_miller
    end

    if :electron_energy_flux in flux_types
        m1d.electrons.energy.flux = @. gyrobohm_energy_flux(ne, Te, rhos, a) * (f.ENERGY_FLUX_e for f in flux_solutions) * vprime_miller
    end

    if :electron_particle_flux in flux_types
        m1d.electrons.particles.flux = @. gyrobohm_particle_flux(ne, Te, rhos, a) * (f.PARTICLE_FLUX_e for f in flux_solutions) * vprime_miller
    end

    if :momentum_flux in flux_types
        m1d.momentum_tor.flux = @. gyrobohm_momentum_flux(ne, Te, rhos, a) * (f.STRESS_TOR_i for f in flux_solutions) * vprime_miller
    end

    if :ion_particle_flux in flux_types
        for (kk, ion) in enumerate(cp1d.ion)
            ion = resize!(m1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label; wipe=false)
            ion.particles.flux = gyrobohm_particle_flux.(ne, Te, rhos, a) .* (pick_ion_flux(f.PARTICLE_FLUX_i, kk) for f in flux_solutions) .* vprime_miller
        end
    end
end

"""
    pick_ion_flux(ion_fluxes::AbstractVector{T}, kk::Int) where {T<:Real}

Select which ion flux to take
"""
function pick_ion_flux(ion_fluxes::AbstractVector{T}, kk::Int) where {T<:Real}
    if isempty(ion_fluxes)
        return zero(T)
    elseif kk <= length(ion_fluxes)
        return ion_fluxes[kk]
    else
        return ion_fluxes[end]
    end
end

