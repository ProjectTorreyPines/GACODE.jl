module GACODE

using CoordinateConventions: transform_cocos
using IMAS.IMASdd: @cocos2
using IMAS
using IMAS: cgs, mks
using IMASutils: argmin_abs
using Printf
using Statistics: mean

# Include the separate modules
include("physics.jl")
include("fluxes.jl")
include("inputgacode.jl")

# Export all functions from physics.jl
export c_s, rho_s, r_min_core_profiles, bunit

# Export all functions and types from fluxes.jl
export FluxSolution, gyrobohm_energy_flux, gyrobohm_particle_flux, gyrobohm_momentum_flux
export volume_prime_miller_correction, flux_gacode_to_imas, pick_ion_flux, sources_to_gyrobohm

# Export functions and types from inputgacode.jl
export InputGACODE, save, load, dd!

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end