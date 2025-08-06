module GACODE

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

export FluxSolution

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

export c_s

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

export rho_s

"""
    r_min_core_profiles(eqt::IMAS.equilibrium__time_slice, rho_tor_norm::AbstractVector)

Geometric minor radius in [cm] evaluated
"""
function r_min_core_profiles(eqt1d::IMAS.equilibrium__time_slice___profiles_1d, rho_tor_norm::AbstractVector)
    return IMAS.interp1d(eqt1d.rho_tor_norm, cgs.m_to_cm * 0.5 * (eqt1d.r_outboard - eqt1d.r_inboard)).(rho_tor_norm)
end

export r_min_core_profiles

##### Gyrobohm normalizations from gacode
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

export gyrobohm_energy_flux

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

export gyrobohm_particle_flux

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

export gyrobohm_momentum_flux

"""
    volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)

Correction to account for transformation from Miller r grid in GA code equilibrium to Psi grid in FUSE equilibrium
"""
function volume_prime_miller_correction(eqt::IMAS.equilibrium__time_slice)
    a_minor = (eqt.profiles_1d.r_outboard .- eqt.profiles_1d.r_inboard) ./ 2.0
    return IMAS.gradient(a_minor, eqt.profiles_1d.volume) ./ eqt.profiles_1d.surface
end

export volume_prime_miller_correction

"""
    flux_gacode_to_imas(
    	flux_types::Tuple{Vararg{Symbol}},
    	flux_solutions::Vector{<:GACODE.FluxSolution},
    	m1d::IMAS.core_transport__model___profiles_1d,
    	eqt::IMAS.equilibrium__time_slice,
    	cp1d::IMAS.core_profiles__profiles_1d
    )

Normalizes specified transport fluxes output by GA code via gyrobohm normalization and Miller volume correction
"""
function flux_gacode_to_imas(
    flux_types::Tuple{Vararg{Symbol}},
    flux_solutions::Vector{<:GACODE.FluxSolution},
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
            ion = resize!(m1d.ion, "element[1].a" => ion.element[1].z_n, "element[1].z_n" => ion.element[1].z_n, "label" => ion.label)
            ion.particles.flux = gyrobohm_particle_flux.(ne, Te, rhos, a) .* (pick_ion_flux(f.PARTICLE_FLUX_i, kk) for f in flux_solutions) .* vprime_miller
        end
    end
end

export flux_gacode_to_imas

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



# Define the mutable struct for GACODE input
Base.@kwdef mutable struct InputGACODE
    # Header information
    header_lines::Union{Vector{String},Missing} = missing

    # Global parameters
    nion::Union{Int,Missing} = missing
    nexp::Union{Int,Missing} = missing
    shot::Union{Int,Missing} = missing
    time::Union{Int,Missing} = missing
    bcentr::Union{Real,Missing} = missing
    current::Union{Real,Missing} = missing
    rcentr::Union{Real,Missing} = missing
    torfluxa::Union{Real,Missing} = missing

    # Ion species information
    name::Union{Vector{AbstractString},Missing} = missing
    type::Union{Vector{AbstractString},Missing} = missing
    mass::Union{Vector{Real},Missing} = missing
    z::Union{Vector{Real},Missing} = missing

    # Electron profiles
    ne::Union{Vector{Real},Missing} = missing
    te::Union{Vector{Real},Missing} = missing

    # Ion profiles
    ni::Union{Matrix{Real},Missing} = missing
    ti::Union{Matrix{Real},Missing} = missing
    vpol::Union{Matrix{Real},Missing} = missing
    vtor::Union{Matrix{Real},Missing} = missing

    # Radial coordinate and profiles
    rho::Union{Vector{Real},Missing} = missing
    polflux::Union{Vector{Real},Missing} = missing
    rmin::Union{Vector{Real},Missing} = missing
    rmaj::Union{Vector{Real},Missing} = missing
    kappa::Union{Vector{Real},Missing} = missing
    delta::Union{Vector{Real},Missing} = missing
    zeta::Union{Vector{Real},Missing} = missing
    zmag::Union{Vector{Real},Missing} = missing
    q::Union{Vector{Real},Missing} = missing
    fpol::Union{Vector{Real},Missing} = missing

    # Shapes!
    shape_cos0::Union{Vector{Real},Missing} = missing
    shape_cos1::Union{Vector{Real},Missing} = missing
    shape_cos2::Union{Vector{Real},Missing} = missing
    shape_cos3::Union{Vector{Real},Missing} = missing
    shape_cos4::Union{Vector{Real},Missing} = missing
    shape_cos5::Union{Vector{Real},Missing} = missing
    shape_cos6::Union{Vector{Real},Missing} = missing
    shape_sin3::Union{Vector{Real},Missing} = missing
    shape_sin4::Union{Vector{Real},Missing} = missing
    shape_sin5::Union{Vector{Real},Missing} = missing
    shape_sin6::Union{Vector{Real},Missing} = missing

    # Currents
    johm::Union{Vector{Real},Missing} = missing
    jbs::Union{Vector{Real},Missing} = missing
    jrf::Union{Vector{Real},Missing} = missing
    jnb::Union{Vector{Real},Missing} = missing
    jbstor::Union{Vector{Real},Missing} = missing
    sigmapar::Union{Vector{Real},Missing} = missing


    # heating   
    qbeame::Union{Vector{Real},Missing} = missing
    qbeami::Union{Vector{Real},Missing} = missing
    qfuse::Union{Vector{Real},Missing} = missing
    qfusi::Union{Vector{Real},Missing} = missing
    qohme::Union{Vector{Real},Missing} = missing
    qpar_beam::Union{Vector{Real},Missing} = missing
    qmom::Union{Vector{Real},Missing} = missing
    qrfe::Union{Vector{Real},Missing} = missing
    qrfi::Union{Vector{Real},Missing} = missing
    qline::Union{Vector{Real},Missing} = missing
    qbrem::Union{Vector{Real},Missing} = missing
    qsync::Union{Vector{Real},Missing} = missing
    qei::Union{Vector{Real},Missing} = missing
    qione::Union{Vector{Real},Missing} = missing
    qioni::Union{Vector{Real},Missing} = missing
    qcxi::Union{Vector{Real},Missing} = missing
    qpar_wall::Union{Vector{Real},Missing} = missing

    # Other profiles
    ptot::Union{Vector{Real},Missing} = missing
    z_eff::Union{Vector{Real},Missing} = missing
    w0::Union{Vector{Real},Missing} = missing

end


function save(input_gacode::InputGACODE, filename::String)

    nexp = input_gacode.nexp
    nion = input_gacode.nion

    # Write header
    open(filename, "w") do io

        println(io, "#  *original : null")
        println(io, "# *statefile : null")
        println(io, "#     *gfile : null")
        println(io, "#   *cerfile : null")
        println(io, "#      *vgen : null")
        println(io, "#     *tgyro : null")

        println(io, "#")
        # Write data

        expro_writei(io, input_gacode.nexp, "nexp")
        expro_writei(io, input_gacode.nion, "nion")
        expro_writei(io, input_gacode.shot, "shot")
        expro_writei(io, input_gacode.time, "time")

        println(io, "# " * "name")
        println(io, join([strip(input_gacode.name[i]) for i in 1:nion], " "))

        println(io, "# " * "type") #[label, Z, A, type]
        println(io, join([strip(input_gacode.type[i]) for i in 1:nion], " "))

        println(io, "# " * "masse")
        @printf(io, "%.7e\n", IMAS.mks.m_e / IMAS.mks.m_p)

        println(io, "# " * "mass")
        println(io, join([@sprintf("%.7e", input_gacode.mass[i]) for i in 1:nion], " "))

        println(io, "# " * "ze")
        @printf(io, "%.7e\n", -1.0)

        println(io, "# " * "z")
        println(io, join([@sprintf("%.7e", input_gacode.z[i]) for i in 1:nion], " "))
        # Write vector/array data, skipping objects that are 0.0

        expro_writes(io, input_gacode.torfluxa, "torfluxa", "Wb/radian")

        expro_writes(io, input_gacode.rcentr, "rcentr", "m")
        expro_writes(io, input_gacode.bcentr, "bcentr", "T")
        expro_writes(io, input_gacode.current, "current", "MA")
        expro_writev(io, input_gacode.rho, nexp, "rho", "-")
        expro_writev(io, input_gacode.rmin, nexp, "rmin", "m")
        expro_writev(io, input_gacode.polflux, nexp, "polflux", "Wb/radian")
        expro_writev(io, input_gacode.q, nexp, "q", "-")
        expro_writev(io, input_gacode.w0, nexp, "w0", "rad/s")
        expro_writev(io, input_gacode.rmaj, nexp, "rmaj", "m")
        expro_writev(io, input_gacode.zmag, nexp, "zmag", "m")
        expro_writev(io, input_gacode.kappa, nexp, "kappa", "-")
        expro_writev(io, input_gacode.delta, nexp, "delta", "-")
        expro_writev(io, input_gacode.zeta, nexp, "zeta", "-")
        expro_writev(io, input_gacode.shape_cos0, nexp, "shape_cos0", "-")
        expro_writev(io, input_gacode.shape_cos1, nexp, "shape_cos1", "-")
        expro_writev(io, input_gacode.shape_cos2, nexp, "shape_cos2", "-")
        expro_writev(io, input_gacode.shape_cos3, nexp, "shape_cos3", "-")
        expro_writev(io, input_gacode.shape_cos4, nexp, "shape_cos4", "-")
        expro_writev(io, input_gacode.shape_cos5, nexp, "shape_cos5", "-")
        expro_writev(io, input_gacode.shape_cos6, nexp, "shape_cos6", "-")
        expro_writev(io, input_gacode.shape_sin3, nexp, "shape_sin3", "-")
        expro_writev(io, input_gacode.shape_sin4, nexp, "shape_sin4", "-")
        expro_writev(io, input_gacode.shape_sin5, nexp, "shape_sin5", "-")
        expro_writev(io, input_gacode.shape_sin6, nexp, "shape_sin6", "-")
        expro_writev(io, input_gacode.ne, nexp, "ne", "10^19/m^3")
        expro_writea(io, input_gacode.ni, nion, nexp, "ni", "10^19/m^3")
        expro_writev(io, input_gacode.te, nexp, "te", "keV")
        expro_writea(io, input_gacode.ti, nion, nexp, "ti", "keV")
        expro_writev(io, input_gacode.ptot, nexp, "ptot", "Pa")
        expro_writev(io, input_gacode.fpol, nexp, "fpol", "T-m")
        expro_writev(io, input_gacode.johm, nexp, "johm", "MA/m^2")
        expro_writev(io, input_gacode.jbs, nexp, "jbs", "MA/m^2")
        expro_writev(io, input_gacode.jrf, nexp, "jrf", "MA/m^2")
        expro_writev(io, input_gacode.jnb, nexp, "jnb", "MA/m^2")
        expro_writev(io, input_gacode.jbstor, nexp, "jbstor", "MA/m^2")
        expro_writev(io, input_gacode.sigmapar, nexp, "sigmapar", "MSiemens/m")
        expro_writev(io, input_gacode.z_eff, nexp, "z_eff", "-")
        expro_writea(io, input_gacode.vpol, nion, nexp, "vpol", "m/s")
        expro_writea(io, input_gacode.vtor, nion, nexp, "vtor", "m/s")
        expro_writev(io, input_gacode.qohme, nexp, "qohme", "MW/m^3")
        expro_writev(io, input_gacode.qbeame, nexp, "qbeame", "MW/m^3")
        expro_writev(io, input_gacode.qbeami, nexp, "qbeami", "MW/m^3")
        expro_writev(io, input_gacode.qrfe, nexp, "qrfe", "MW/m^3")
        expro_writev(io, input_gacode.qrfi, nexp, "qrfi", "MW/m^3")
        expro_writev(io, input_gacode.qfuse, nexp, "qfuse", "MW/m^3")
        expro_writev(io, input_gacode.qfusi, nexp, "qfusi", "MW/m^3")
        expro_writev(io, input_gacode.qbrem, nexp, "qbrem", "MW/m^3")
        expro_writev(io, input_gacode.qsync, nexp, "qsync", "MW/m^3")
        expro_writev(io, input_gacode.qline, nexp, "qline", "MW/m^3")
        expro_writev(io, input_gacode.qei, nexp, "qei", "MW/m^3")
        expro_writev(io, input_gacode.qione, nexp, "qione", "MW/m^3")
        expro_writev(io, input_gacode.qioni, nexp, "qioni", "MW/m^3")
        expro_writev(io, input_gacode.qcxi, nexp, "qcxi", "MW/m^3")
        expro_writev(io, input_gacode.qpar_beam, nexp, "qpar_beam", "1/m^3/s")
        expro_writev(io, input_gacode.qpar_wall, nexp, "qpar_wall", "1/m^3/s")
        return expro_writev(io, input_gacode.qmom, nexp, "qmom", "N/m^2")
    end
end

function expro_writes(io::IO, x::Real, xs1::String, xs2::String)
    """
    Write scalar value with identifier and units if non-zero
    """
    if abs(x) > 1e-16
        println(io, "# " * xs1 * " | " * xs2)
        @printf(io, "%.7e\n", x)
    end
end

function expro_writei(io::IO, i::Integer, xs1::String)
    """
    Write integer value with identifier if positive
    """

    if i > 0
        println(io, "# " * xs1)
        println(io, i)
    end
end

function expro_writev(io::IO, x::Union{Vector{<:Real},Missing}, n::Integer, xs1::String, xs2::String)
    """
    Write vector values with indices if non-zero
    """
    if ~ismissing(x) && sum(abs.(x)) > 1e-16
        println(io, "# " * xs1 * " | " * xs2)
        for i in 1:n
            @printf(io, "%3d %.7e\n", i, x[i])
        end
    end
end

function expro_writea(io::IO, x::Union{Matrix{<:Real},Missing}, m::Integer, n::Integer, xs1::String, xs2::String)
    """
    Write array values (matrix) with row indices if non-zero
    """
    if ~ismissing(x) && sum(abs.(x)) > 1e-16
        println(io, "# " * xs1 * " | " * xs2)
        for i in 1:n
            @printf(io, "%3d ", i)
            for j in 1:m
                @printf(io, "%.7e  ", x[j, i])
            end
            println(io)
        end
    end
end

function update_hcd(dd, input_gacode)
    ngrid = length(input_gacode.rho)
    energy_scale = 1e-6
    particle_scale = 1.0
    for item in ["qbeame", "qbeami", "qpar_beam", "qmom", "qrfe", "qrfi", "qline", "qbrem", "qsync", "qei", "qpar_wall", "qohme", "qfuse", "qfusi"]
        setproperty!(input_gacode, Symbol(item), zeros(ngrid))
    end

    for source in findall(:nbi, dd.core_sources.source)
        input_gacode.qbeame += source.profiles_1d[1].electrons.energy * energy_scale
        input_gacode.qbeami += source.profiles_1d[1].total_ion_energy * energy_scale
        input_gacode.qpar_beam += source.profiles_1d[1].electrons.particles * particle_scale
        input_gacode.qmom += source.profiles_1d[1].momentum_tor
    end

    # Add sawtooth to beam since is the only source with all channels
    for source in findall(:sawteeth, dd.core_sources.source)
        input_gacode.qbeame += source.profiles_1d[1].electrons.energy * energy_scale
        input_gacode.qbeami += source.profiles_1d[1].total_ion_energy * energy_scale
        input_gacode.qpar_beam += source.profiles_1d[1].electrons.particles * particle_scale
        input_gacode.qmom += source.profiles_1d[1].momentum_tor
    end

    for sid in [:ic, :ec, :lh]
        for source in findall(sid, dd.core_sources.source)
            input_gacode.qrfe += source.profiles_1d[1].electrons.energy * energy_scale
            input_gacode.qrfi += source.profiles_1d[1].total_ion_energy * energy_scale
        end
    end

    for source in findall(:line_radiation, dd.core_sources.source)
        input_gacode.qline -= source.profiles_1d[1].electrons.energy * energy_scale
    end

    for source in findall(:fusion, dd.core_sources.source)
        input_gacode.qfuse += source.profiles_1d[1].electrons.energy * energy_scale
        input_gacode.qfusi += source.profiles_1d[1].total_ion_energy * energy_scale
    end

    for source in findall(:bremsstrahlung, dd.core_sources.source)
        input_gacode.qbrem -= source.profiles_1d[1].electrons.energy * energy_scale
    end

    for source in findall(:synchrotron_radiation, dd.core_sources.source)
        input_gacode.qsync -= source.profiles_1d[1].electrons.energy * energy_scale
    end

    for source in findall(:collisional_equipartition, dd.core_sources.source)
        input_gacode.qei -= source.profiles_1d[1].electrons.energy * energy_scale
    end

    for source in findall(:gas_puff, dd.core_sources.source)
        input_gacode.qpar_wall += source.profiles_1d[1].electrons.particles * particle_scale
    end

    for source in findall(:ohmic, dd.core_sources.source)
        input_gacode.qohme += source.profiles_1d[1].electrons.energy * energy_scale
    end
end


"""
    process_ion_species(ion, k, cp1d)

Handles the ion species and splits DT into D and T, this assumes 50/50 DT!
"""
function process_ion_species(ion::IMAS.core_profiles__profiles_1d___ion, k::Int, cp1d::IMAS.core_profiles__profiles_1d)
    """
    Process a single ion species and return a vector of species dictionaries.
    Handles DT splitting into separate D and T species.
    """
    species_list = []
    
    A = ion.element[1].a
    Z = ion.element[1].z_n
    label = ion.label
    
    # Check if this is a DT species that needs to be split
    if label == "DT"
        # Create D species
        d_species = Dict(
            :label => "D",
            :z => 1,
            :mass => 2.014,  # Deuterium atomic mass
            :density_thermal => ion.density_thermal * 0.5,  # Split 50/50
            :density_fast => ion.density_fast * 0.5,
            :temperature => cp1d.ion[k].temperature,
            :pressure_fast_perpendicular => ion.pressure_fast_perpendicular * 0.5,
            :pressure_fast_parallel => ion.pressure_fast_parallel * 0.5
        )
        push!(species_list, d_species)
        
        # Create T species
        t_species = Dict(
            :label => "T",
            :z => 1,
            :mass => 3.016,  # Tritium atomic mass
            :density_thermal => ion.density_thermal * 0.5,  # Split 50/50
            :density_fast => ion.density_fast * 0.5,
            :temperature => cp1d.ion[k].temperature,
            :pressure_fast_perpendicular => ion.pressure_fast_perpendicular * 0.5,
            :pressure_fast_parallel => ion.pressure_fast_parallel * 0.5
        )
        push!(species_list, t_species)
    else
        # Regular species - no splitting needed
        regular_species = Dict(
            :label => label,
            :z => Z,
            :mass => A,
            :density_thermal => ion.density_thermal,
            :density_fast => ion.density_fast,
            :temperature => cp1d.ion[k].temperature,
            :pressure_fast_perpendicular => ion.pressure_fast_perpendicular,
            :pressure_fast_parallel => ion.pressure_fast_parallel
        )
        push!(species_list, regular_species)
    end
    
    return species_list
end

"""
    InputGACODE(dd::IMAS.dd)

Creates an input_gacode from dd
"""
function InputGACODE(dd::IMAS.dd)
    input_gacode = InputGACODE()
    cocosio = 2  # GACODE uses COCOS 2

    rho = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d

    # Set basic parameters
    input_gacode.bcentr = -1.0 * (@ddtime dd.equilibrium.vacuum_toroidal_field.b0)
    input_gacode.current = -1.0 * eqt.global_quantities.ip / 1e6
    input_gacode.rcentr = dd.equilibrium.vacuum_toroidal_field.r0
    input_gacode.rho = rho
    input_gacode.nexp = length(rho)

    if ismissing(dd.dataset_description.data_entry, :pulse)
        input_gacode.shot = 1
    else
        input_gacode.shot = dd.dataset_description.data_entry.pulse
    end

    if iszero(dd.global_time)
        input_gacode.time = 1
    else
        input_gacode.time = Int64(floor(1e3 * dd.global_time))
    end

    # Set geometric profiles
    rho_eq = eqt1d.rho_tor_norm
    input_gacode.polflux = IMAS.interp1d(rho_eq, eqt1d.psi .- eqt1d.psi[1]).(rho)
    input_gacode.rmin = GACODE.r_min_core_profiles(eqt1d, rho) / IMAS.cgs.m_to_cm
    input_gacode.rmaj = IMAS.interp1d(rho_eq, 0.5 .* (eqt1d.r_outboard .+ eqt1d.r_inboard)).(rho)
    input_gacode.kappa = IMAS.interp1d(rho_eq, eqt1d.elongation).(rho)
    input_gacode.delta = IMAS.interp1d(rho_eq, 0.5 .* (eqt1d.triangularity_lower .+ eqt1d.triangularity_upper)).(rho)
    input_gacode.zeta =
        IMAS.interp1d(rho_eq, 0.25 .* (
            eqt1d.squareness_upper_outer .+ eqt1d.squareness_upper_inner .+
            eqt1d.squareness_lower_outer .+ eqt1d.squareness_lower_inner
        )).(rho)
    input_gacode.zmag = IMAS.interp1d(rho_eq, eqt1d.geometric_axis.z).(rho)
    input_gacode.q = IMAS.interp1d(rho_eq, eqt1d.q).(rho)
    input_gacode.torfluxa = -1.0 .* eqt1d.phi[end] ./ 2π
    
    # Set electron profiles
    input_gacode.ne = cp1d.electrons.density / 1e19
    input_gacode.te = cp1d.electrons.temperature / 1e3

    # Initialize total pressure with electron contribution
    input_gacode.ptot = @. cp1d.electrons.density_thermal * cp1d.electrons.temperature * IMAS.mks.e

    # Process all ion species and create the expanded species list
    all_species = []
    for (k, ion) in enumerate(cp1d.ion)
        species_list = process_ion_species(ion, k, cp1d)
        append!(all_species, species_list)
    end

    # Count thermal and fast species
    nion = 0
    for species in all_species
        if sum(abs.(species[:density_thermal])) > 0
            nion += 1
        end
        if sum(abs.(species[:density_fast])) > 0
            nion += 1
        end
    end

    input_gacode.nion = nion

    # Initialize arrays
    input_gacode.ni = zeros(nion, input_gacode.nexp)
    input_gacode.ti = zeros(nion, input_gacode.nexp)
    input_gacode.vtor = zeros(nion, input_gacode.nexp)
    input_gacode.vpol = zeros(nion, input_gacode.nexp)

    input_gacode.type = Vector{String}(undef, nion)
    input_gacode.name = Vector{String}(undef, nion)
    input_gacode.z = zeros(nion)
    input_gacode.mass = zeros(nion)

    # Populate the arrays with processed species
    i = 0
    for species in all_species
        # Handle thermal ions
        if sum(abs.(species[:density_thermal])) > 0
            i += 1
            input_gacode.type[i] = "thermal"
            input_gacode.name[i] = species[:label]
            input_gacode.z[i] = species[:z]
            input_gacode.mass[i] = species[:mass]

            input_gacode.ni[i, :] = species[:density_thermal] / 1e19
            input_gacode.ti[i, :] = species[:temperature] / 1e3
            input_gacode.vtor[i, :] = 0.0 * species[:density_thermal]
            input_gacode.vpol[i, :] = 0.0 * species[:density_thermal]
            input_gacode.ptot += species[:density_thermal] .* species[:temperature] .* IMAS.mks.e
        end

        # Handle fast ions
        if sum(abs.(species[:density_fast])) > 0
            i += 1
            input_gacode.type[i] = "fast"
            input_gacode.name[i] = species[:label]
            input_gacode.z[i] = species[:z]
            input_gacode.mass[i] = species[:mass]

            Ti_fast = (((2 .* species[:pressure_fast_perpendicular] .+ species[:pressure_fast_parallel]) ./ species[:density_fast]) ./ IMAS.mks.e ./ 1e3)
            ni_fast = species[:density_fast] / 1e19

            # Set finite temperature when density ~0 for TGYRO splines
            navg = mean(species[:density_fast])
            Ti_fast[species[:density_fast].<1e-6*navg] .= mean(Ti_fast)

            input_gacode.ni[i, :] = ni_fast
            input_gacode.ti[i, :] = Ti_fast
            input_gacode.vtor[i, :] = 0.0 * species[:density_thermal]
            input_gacode.vpol[i, :] = 0.0 * species[:density_thermal]

            input_gacode.ptot += 2 .* species[:pressure_fast_perpendicular] .+ species[:pressure_fast_parallel]
        end
    end

    input_gacode.z_eff = cp1d.zeff
    input_gacode.w0 = cp1d.rotation_frequency_tor_sonic

    update_hcd(dd, input_gacode)

    return input_gacode
end

function load(filename::String)
    input_gacode = InputGACODE()
    lines = readlines(filename)
    function get_varname(line)
        if startswith(line, "# ")
            return split(line, " ")[2]
        else
            return nothing
        end
    end

    for field_name in fieldnames(InputGACODE)

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Vector{AbstractString}}
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            setproperty!(input_gacode, field_name, collect(split(lines[iline+1], " ")))
        end

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Int}
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            setproperty!(input_gacode, field_name, parse(Int, lines[iline+1]))
        end

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Real}
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            setproperty!(input_gacode, field_name, parse(Float64, lines[iline+1]))
        end

        nexp = input_gacode.nexp
        nion = input_gacode.nion
        if field_name == :mass || field_name == :z #nion length vectors
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            setproperty!(input_gacode, field_name, parse.(Float64, split(lines[iline+1])))
        elseif fieldtype(typeof(input_gacode), field_name) == Union{Missing,Vector{Real}} #nexp length vectors
            tmp_array = zeros(nexp)
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)

            if ~isnothing(iline)
                for (i, line) in enumerate(lines[(iline+1):(iline+nexp)])
                    tmp_array[i] = parse(Float64, split(lines[iline+i])[2])
                end
                setproperty!(input_gacode, field_name, tmp_array)
            end
        end

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Matrix{Real}}
            tmp_array = zeros(nion, nexp)
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            if ~isnothing(iline)
                for (i, line) in enumerate(lines[(iline+1):(iline+nexp)])
                    tmp_array[:, i] = parse.(Float64, split(lines[iline+i])[2:nion+1])
                end
                setproperty!(input_gacode, field_name, tmp_array)
            end
        end
    end
    return input_gacode
end



const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
