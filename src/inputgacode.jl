using IMAS
using IMAS: cgs, mks
using IMASutils: argmin_abs
using Printf

Base.@kwdef mutable struct InputGACODE
    # Header information
    header_lines::Union{Vector{String},Missing} = missing

    # Global parameters
    nion::Union{Int,Missing} = missing
    nexp::Union{Int,Missing} = missing
    shot::Union{Int,Missing} = missing
    time::Union{Int,Missing} = missing
    bcentr::Union{Float64,Missing} = missing
    current::Union{Float64,Missing} = missing
    rcentr::Union{Float64,Missing} = missing
    torfluxa::Union{Float64,Missing} = missing

    # Ion species information
    name::Union{Vector{String},Missing} = missing
    type::Union{Vector{String},Missing} = missing
    mass::Union{Vector{Float64},Missing} = missing
    z::Union{Vector{Float64},Missing} = missing

    # Electron profiles
    ne::Union{Vector{Float64},Missing} = missing
    te::Union{Vector{Float64},Missing} = missing

    # Ion profiles
    ni::Union{Matrix{Float64},Missing} = missing
    ti::Union{Matrix{Float64},Missing} = missing
    vpol::Union{Matrix{Float64},Missing} = missing
    vtor::Union{Matrix{Float64},Missing} = missing

    # Radial coordinate and profiles
    rho::Union{Vector{Float64},Missing} = missing
    polflux::Union{Vector{Float64},Missing} = missing
    rmin::Union{Vector{Float64},Missing} = missing
    rmaj::Union{Vector{Float64},Missing} = missing
    kappa::Union{Vector{Float64},Missing} = missing
    delta::Union{Vector{Float64},Missing} = missing
    zeta::Union{Vector{Float64},Missing} = missing
    zmag::Union{Vector{Float64},Missing} = missing
    q::Union{Vector{Float64},Missing} = missing
    fpol::Union{Vector{Float64},Missing} = missing

    # Shapes!
    shape_cos0::Union{Vector{Float64},Missing} = missing
    shape_cos1::Union{Vector{Float64},Missing} = missing
    shape_cos2::Union{Vector{Float64},Missing} = missing
    shape_cos3::Union{Vector{Float64},Missing} = missing
    shape_cos4::Union{Vector{Float64},Missing} = missing
    shape_cos5::Union{Vector{Float64},Missing} = missing
    shape_cos6::Union{Vector{Float64},Missing} = missing
    shape_sin3::Union{Vector{Float64},Missing} = missing
    shape_sin4::Union{Vector{Float64},Missing} = missing
    shape_sin5::Union{Vector{Float64},Missing} = missing
    shape_sin6::Union{Vector{Float64},Missing} = missing

    # Currents
    johm::Union{Vector{Float64},Missing} = missing
    jbs::Union{Vector{Float64},Missing} = missing
    jrf::Union{Vector{Float64},Missing} = missing
    jnb::Union{Vector{Float64},Missing} = missing
    jbstor::Union{Vector{Float64},Missing} = missing
    sigmapar::Union{Vector{Float64},Missing} = missing

    # heating
    qbeame::Union{Vector{Float64},Missing} = missing
    qbeami::Union{Vector{Float64},Missing} = missing
    qfuse::Union{Vector{Float64},Missing} = missing
    qfusi::Union{Vector{Float64},Missing} = missing
    qohme::Union{Vector{Float64},Missing} = missing
    qpar_beam::Union{Vector{Float64},Missing} = missing
    qmom::Union{Vector{Float64},Missing} = missing
    qrfe::Union{Vector{Float64},Missing} = missing
    qrfi::Union{Vector{Float64},Missing} = missing
    qline::Union{Vector{Float64},Missing} = missing
    qbrem::Union{Vector{Float64},Missing} = missing
    qsync::Union{Vector{Float64},Missing} = missing
    qei::Union{Vector{Float64},Missing} = missing
    qione::Union{Vector{Float64},Missing} = missing
    qioni::Union{Vector{Float64},Missing} = missing
    qcxi::Union{Vector{Float64},Missing} = missing
    qpar_wall::Union{Vector{Float64},Missing} = missing

    # Other profiles
    ptot::Union{Vector{Float64},Missing} = missing
    z_eff::Union{Vector{Float64},Missing} = missing
    w0::Union{Vector{Float64},Missing} = missing
end

canonical_ion_type(token::AbstractString) = begin
    t = lowercase(strip(strip(String(token)), ['[', ']']))
    t in ("therm", "thermal") && return "therm"
    t == "fast" && return "fast"
    return t
end

format_ion_type(token::AbstractString) = "[$(canonical_ion_type(token))]"

# ======================== #
# save InputGACODE to file #
# ======================== #

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
        println(io, join([format_ion_type(input_gacode.type[i]) for i in 1:nion], " "))

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

"""
    expro_writes(io::IO, x::Real, xs1::String, xs2::String)

Write scalar value with identifier and units if non-zero
"""
function expro_writes(io::IO, x::Real, xs1::String, xs2::String)
    if abs(x) > 1e-16
        println(io, "# " * xs1 * " | " * xs2)
        @printf(io, "%.7e\n", x)
    end
end

"""
    expro_writei(io::IO, i::Integer, xs1::String)

Write integer value with identifier if positive
"""
function expro_writei(io::IO, i::Integer, xs1::String)

    if i > 0
        println(io, "# " * xs1)
        println(io, i)
    end
end

"""
    expro_writev(io::IO, x::Union{Vector{<:Real},Missing}, n::Integer, xs1::String, xs2::String)

Write vector values with indices if non-zero
"""
function expro_writev(io::IO, x::Union{Vector{<:Real},Missing}, n::Integer, xs1::String, xs2::String)
    if ~ismissing(x) && sum(abs.(x)) > 1e-16
        println(io, "# " * xs1 * " | " * xs2)
        for i in 1:n
            @printf(io, "%3d %.7e\n", i, x[i])
        end
    end
end

"""
    expro_writea(io::IO, x::Union{Matrix{<:Real},Missing}, m::Integer, n::Integer, xs1::String, xs2::String)

Write array values (matrix) with row indices if non-zero
"""
function expro_writea(io::IO, x::Union{Matrix{<:Real},Missing}, m::Integer, n::Integer, xs1::String, xs2::String)
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

# ========================== #
# load InputGACODE from file #
# ========================== #

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

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Vector{String}}
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            if !isnothing(iline)
                vals = collect(split(lines[iline+1]))
                field_name == :type && (vals = canonical_ion_type.(vals))
                setproperty!(input_gacode, field_name, vals)
            end
        end

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Int}
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            isnothing(iline) || setproperty!(input_gacode, field_name, parse(Int, lines[iline+1]))
        end

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Float64}
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            isnothing(iline) || setproperty!(input_gacode, field_name, parse(Float64, lines[iline+1]))
        end

        nexp = input_gacode.nexp
        nion = input_gacode.nion
        if field_name == :mass || field_name == :z #nion length vectors
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)
            setproperty!(input_gacode, field_name, parse.(Float64, split(lines[iline+1])))
        elseif fieldtype(typeof(input_gacode), field_name) == Union{Missing,Vector{Float64}} #nexp length vectors
            tmp_array = zeros(nexp)
            iline = findfirst(line -> get_varname(line) == String(field_name), lines)

            if ~isnothing(iline)
                for (i, line) in enumerate(lines[(iline+1):(iline+nexp)])
                    tmp_array[i] = parse(Float64, split(lines[iline+i])[2])
                end
                setproperty!(input_gacode, field_name, tmp_array)
            end
        end

        if fieldtype(typeof(input_gacode), field_name) == Union{Missing,Matrix{Float64}}
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

    # Handle legacy file-format alias: "qpar" is the old name for "qpar_beam"
    if ismissing(input_gacode.qpar_beam)
        nexp = input_gacode.nexp
        iline = findfirst(line -> get_varname(line) == "qpar", lines)
        if !isnothing(iline)
            tmp_array = zeros(nexp)
            for i in 1:nexp
                tmp_array[i] = parse(Float64, split(lines[iline+i])[2])
            end
            input_gacode.qpar_beam = tmp_array
        end
    end

    return input_gacode
end

# =================== #
# InputGACODE from dd #
# =================== #

"""
    InputGACODE(dd::IMAS.dd)

Creates an input_gacode from dd
"""
function InputGACODE(dd::IMAS.dd)
    input_gacode = InputGACODE()

    rho = dd.core_profiles.profiles_1d[].grid.rho_tor_norm
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    eqt1d = eqt.profiles_1d
    vtf = dd.equilibrium.vacuum_toroidal_field
    _b0_fac = transform_cocos(IMAS.internal_cocos, 2)["TOR"]

    # Set basic parameters
    input_gacode.bcentr = _b0_fac * (@ddtime vtf.b0)
    input_gacode.current = @cocos2(eqt.global_quantities.ip) / 1e6
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
    psi2 = @cocos2(eqt1d.psi)
    input_gacode.polflux = IMAS.interp1d(rho_eq, psi2 .- psi2[1]).(rho)
    input_gacode.rmin = r_min_core_profiles(eqt1d, rho) / IMAS.cgs.m_to_cm
    input_gacode.rmaj = IMAS.interp1d(rho_eq, 0.5 .* (eqt1d.r_outboard .+ eqt1d.r_inboard)).(rho)
    input_gacode.kappa = IMAS.interp1d(rho_eq, eqt1d.elongation).(rho)
    input_gacode.delta = IMAS.interp1d(rho_eq, 0.5 .* (eqt1d.triangularity_lower .+ eqt1d.triangularity_upper)).(rho)
    input_gacode.zeta =
        IMAS.interp1d(rho_eq, 0.25 .* (
            eqt1d.squareness_upper_outer .+ eqt1d.squareness_upper_inner .+
            eqt1d.squareness_lower_outer .+ eqt1d.squareness_lower_inner
        )).(rho)
    input_gacode.zmag = IMAS.interp1d(rho_eq, eqt1d.geometric_axis.z).(rho)
    input_gacode.q = IMAS.interp1d(rho_eq, @cocos2(eqt1d.q)).(rho)
    input_gacode.torfluxa = @cocos2(eqt1d.phi)[end] / 2π

    # fpol
    if !ismissing(eqt1d, :f)
        input_gacode.fpol = IMAS.interp1d(rho_eq, @cocos2(eqt1d.f)).(rho)
    end

    # Shape harmonics from profiles_2d (grid_type 57, TEQUILA convention)
    for eq2d in eqt.profiles_2d
        !ismissing(eq2d.grid_type, :index) && eq2d.grid_type.index == 57 || continue
        rho2 = eq2d.grid.dim1
        mat  = eq2d.psi
        nm   = (size(mat, 2) - 5) ÷ 2
        itp(col) = IMAS.interp1d(rho2, mat[:, col]).(rho)
        nm >= 1 && (input_gacode.shape_cos0 = itp(5))
        nm >= 1 && (input_gacode.shape_cos1 = itp(6))
        nm >= 2 && (input_gacode.shape_cos2 = itp(7))
        nm >= 3 && (input_gacode.shape_cos3 = itp(8))
        nm >= 4 && (input_gacode.shape_cos4 = itp(9))
        nm >= 5 && (input_gacode.shape_cos5 = itp(10))
        nm >= 6 && (input_gacode.shape_cos6 = itp(11))
        nm >= 3 && (input_gacode.shape_sin3 = itp(5 + nm + 3))
        nm >= 4 && (input_gacode.shape_sin4 = itp(5 + nm + 4))
        nm >= 5 && (input_gacode.shape_sin5 = itp(5 + nm + 5))
        nm >= 6 && (input_gacode.shape_sin6 = itp(5 + nm + 6))
        break
    end

    # Set electron profiles
    input_gacode.ne = cp1d.electrons.density / 1e19
    input_gacode.te = cp1d.electrons.temperature / 1e3

    # Process all ion species - now returns single species per ion
    all_species = [ion_as_dict(ion) for ion in cp1d.ion]

    # Count thermal and fast species separately
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
            input_gacode.type[i] = "therm"
            input_gacode.name[i] = species[:label]
            input_gacode.z[i] = species[:z]
            input_gacode.mass[i] = species[:mass]

            input_gacode.ni[i, :] = species[:density_thermal] / 1e19
            input_gacode.ti[i, :] = species[:temperature] / 1e3
            input_gacode.vtor[i, :] = 0.0 * species[:density_thermal]
            input_gacode.vpol[i, :] = 0.0 * species[:density_thermal]
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

        end
    end

    # Use directly stored total pressure when available; fall back to species reconstruction.
    if IMAS.hasdata(cp1d, :pressure)
        input_gacode.ptot = cp1d.pressure
    else
        input_gacode.ptot = @. cp1d.electrons.density_thermal * cp1d.electrons.temperature * IMAS.mks.e
        for species in all_species
            if sum(abs.(species[:density_thermal])) > 0
                input_gacode.ptot .+= species[:density_thermal] .* species[:temperature] .* IMAS.mks.e
            end
            if sum(abs.(species[:density_fast])) > 0
                input_gacode.ptot .+= 2 .* species[:pressure_fast_perpendicular] .+ species[:pressure_fast_parallel]
            end
        end
    end

    # Second pass: read vtor from ion.rotation_frequency_tor (preferred) or ion.velocity.toroidal
    rmaj = input_gacode.rmaj
    i = 0
    has_any_vpol = false
    for (k_ion, species) in enumerate(all_species)
        ion = cp1d.ion[k_ion]
        if sum(abs.(species[:density_thermal])) > 0
            i += 1
            if IMAS.hasdata(ion, :rotation_frequency_tor)
                input_gacode.vtor[i, :] = ion.rotation_frequency_tor .* rmaj
            elseif !ismissing(ion.velocity, :toroidal)
                input_gacode.vtor[i, :] = ion.velocity.toroidal
            end
            if !ismissing(ion.velocity, :poloidal)
                input_gacode.vpol[i, :] = ion.velocity.poloidal
                has_any_vpol = true
            end
        end
        if sum(abs.(species[:density_fast])) > 0
            i += 1  # fast-ion vtor stays zero (initialized to zeros)
        end
    end
    has_any_vpol || (input_gacode.vpol = missing)

    input_gacode.z_eff = cp1d.zeff
    input_gacode.w0 = cp1d.rotation_frequency_tor_sonic

    update_hcd(dd, input_gacode)

    # Current profiles — use hasdata to detect directly stored values (vs dynamic expressions)
    if IMAS.hasdata(cp1d, :j_ohmic)
        input_gacode.johm = cp1d.j_ohmic ./ 1e6
    end
    if IMAS.hasdata(cp1d, :j_bootstrap)
        input_gacode.jbs = cp1d.j_bootstrap ./ 1e6
    end
    # Prefer the dedicated toroidal bootstrap channel, with a guarded fallback:
    # dynamic expression may fail when equilibrium metric coefficients are unavailable.
    if hasproperty(cp1d, :j_bootstrap_tor)
        try
            jbstor_val = cp1d.j_bootstrap_tor
            if !ismissing(jbstor_val)
                input_gacode.jbstor = jbstor_val ./ 1e6
            end
        catch
            # Leave jbstor as missing when conversion data are unavailable.
        end
    end
    if IMAS.hasdata(cp1d, :conductivity_parallel)
        input_gacode.sigmapar = cp1d.conductivity_parallel ./ 1e6
    end

    # RF/NB current densities summed from core_sources
    jrf_val = zeros(length(rho))
    has_jrf = false
    for sid in [:ic, :ec, :lh]
        for src in findall(sid, dd.core_sources.source)
            isempty(src.profiles_1d) && continue
            ismissing(src.profiles_1d[1], :j_parallel) && continue
            has_jrf = true
            jrf_val .+= src.profiles_1d[1].j_parallel ./ 1e6
        end
    end
    has_jrf && (input_gacode.jrf = jrf_val)

    jnb_val = zeros(length(rho))
    has_jnb = false
    for src in findall(:nbi, dd.core_sources.source)
        isempty(src.profiles_1d) && continue
        ismissing(src.profiles_1d[1], :j_parallel) && continue
        has_jnb = true
        jnb_val .+= src.profiles_1d[1].j_parallel ./ 1e6
    end
    has_jnb && (input_gacode.jnb = jnb_val)

    return input_gacode
end

function update_hcd(dd, input_gacode)
    energy_scale = 1e-6
    particle_scale = 1.0
    for item in ["qbeame", "qbeami", "qpar_beam", "qmom", "qrfe", "qrfi", "qline", "qbrem", "qsync", "qei", "qpar_wall", "qohme", "qfuse", "qfusi", "qione", "qioni", "qcxi"]
        setproperty!(input_gacode, Symbol(item), missing)
    end

    function accumulate!(field::Symbol, values)
        vals = collect(Float64.(values))
        current = getfield(input_gacode, field)
        if ismissing(current)
            setproperty!(input_gacode, field, vals)
        else
            current .+= vals
        end
    end

    for source in findall(:nbi, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qbeame, p1d.electrons.energy .* energy_scale)
        ismissing(p1d, :total_ion_energy) || accumulate!(:qbeami, p1d.total_ion_energy .* energy_scale)
        ismissing(p1d.electrons, :particles) || accumulate!(:qpar_beam, p1d.electrons.particles .* particle_scale)
        ismissing(p1d, :momentum_tor) || accumulate!(:qmom, p1d.momentum_tor)
    end

    # Add sawtooth to beam since is the only source with all channels
    for source in findall(:sawteeth, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qbeame, p1d.electrons.energy .* energy_scale)
        ismissing(p1d, :total_ion_energy) || accumulate!(:qbeami, p1d.total_ion_energy .* energy_scale)
        ismissing(p1d.electrons, :particles) || accumulate!(:qpar_beam, p1d.electrons.particles .* particle_scale)
        ismissing(p1d, :momentum_tor) || accumulate!(:qmom, p1d.momentum_tor)
    end

    for sid in [:ic, :ec, :lh]
        for source in findall(sid, dd.core_sources.source)
            isempty(source.profiles_1d) && continue
            p1d = source.profiles_1d[1]
            ismissing(p1d.electrons, :energy) || accumulate!(:qrfe, p1d.electrons.energy .* energy_scale)
            ismissing(p1d, :total_ion_energy) || accumulate!(:qrfi, p1d.total_ion_energy .* energy_scale)
        end
    end

    for source in findall(:line_radiation, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qline, .-p1d.electrons.energy .* energy_scale)
    end

    for source in findall(:fusion, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qfuse, p1d.electrons.energy .* energy_scale)
        ismissing(p1d, :total_ion_energy) || accumulate!(:qfusi, p1d.total_ion_energy .* energy_scale)
    end

    for source in findall(:bremsstrahlung, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qbrem, .-p1d.electrons.energy .* energy_scale)
    end

    for source in findall(:synchrotron_radiation, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qsync, .-p1d.electrons.energy .* energy_scale)
    end

    for source in findall(:collisional_equipartition, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        if !ismissing(p1d.electrons, :energy)
            accumulate!(:qei, .-p1d.electrons.energy .* energy_scale)
        elseif !ismissing(p1d, :total_ion_energy)
            accumulate!(:qei, p1d.total_ion_energy .* energy_scale)
        end
    end

    for source in findall(:gas_puff, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :particles) || accumulate!(:qpar_wall, p1d.electrons.particles .* particle_scale)
    end

    for source in findall(:ohmic, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qohme, p1d.electrons.energy .* energy_scale)
    end

    for source in findall(:recombination, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d.electrons, :energy) || accumulate!(:qione, p1d.electrons.energy .* energy_scale)
        ismissing(p1d, :total_ion_energy) || accumulate!(:qioni, p1d.total_ion_energy .* energy_scale)
    end

    for source in findall(:charge_exchange, dd.core_sources.source)
        isempty(source.profiles_1d) && continue
        p1d = source.profiles_1d[1]
        ismissing(p1d, :total_ion_energy) && continue
        accumulate!(:qcxi, p1d.total_ion_energy .* energy_scale)
    end
end

"""    dd(input_gacode::InputGACODE)
Create an IMAS `dd` from an `InputGACODE` structure. Inverse of `InputGACODE(dd::IMAS.dd)`.
"""
function IMAS.dd(input_gacode::InputGACODE)
    dd = IMAS.dd()
    dd!(dd, input_gacode)
    return dd
end

"""
    dd!(dd::IMAS.dd, ig::InputGACODE)

Populate an IMAS `dd` from an `InputGACODE` structure. Inverse of `InputGACODE(dd::IMAS.dd)`.
"""
function dd!(dd::IMAS.dd, ig::InputGACODE)
    # Global time & dataset description
    t = ismissing(ig.time) ? 0.0 : Float64(ig.time) * 1e-3
    dd.global_time = t
    ismissing(ig.shot) || (dd.dataset_description.data_entry.pulse = ig.shot)

    # Vacuum toroidal field
    vtf = dd.equilibrium.vacuum_toroidal_field
    vtf.r0 = ig.rcentr
    _b0_fac = transform_cocos(IMAS.internal_cocos, 2)["TOR"]   # ±1, self-inverse for TOR
    @ddtime(vtf.b0 = _b0_fac * ig.bcentr)

    # Equilibrium time slice
    eqt = resize!(dd.equilibrium.time_slice)
    eqt.time = t
    @cocos2(eqt.global_quantities.ip = ig.current * 1e6)   # MA→A, COCOS 2→internal
    eqt.global_quantities.magnetic_axis.r = ig.rmaj[1]
    eqt.global_quantities.magnetic_axis.z = ig.zmag[1]

    eqt1d = eqt.profiles_1d
    # ig.rmin is already in meters: r_min_core_profiles returns cm, then /m_to_cm gives m
    rmin_m = ig.rmin

    # psi is the coordinate for rho_tor_norm — must be set first
    @cocos2(eqt1d.psi = ig.polflux)
    eqt1d.rho_tor_norm        = ig.rho
    @cocos2(eqt1d.phi = ig.torfluxa .* 2π .* ig.rho .^ 2)
    @cocos2(eqt1d.q   = ig.q)
    ismissing(ig.fpol) || @cocos2(eqt1d.f = ig.fpol)
    eqt1d.r_outboard          = ig.rmaj .+ rmin_m
    eqt1d.r_inboard           = ig.rmaj .- rmin_m
    eqt1d.elongation          = ig.kappa
    eqt1d.triangularity_upper = ig.delta
    eqt1d.triangularity_lower = ig.delta
    eqt1d.geometric_axis.z    = ig.zmag
    if !ismissing(ig.zeta)
        eqt1d.squareness_upper_outer = ig.zeta   # inverse of averaging 4 components
        eqt1d.squareness_upper_inner = ig.zeta
        eqt1d.squareness_lower_outer = ig.zeta
        eqt1d.squareness_lower_inner = ig.zeta
    end

    # Core profiles — electrons
    cp1d = resize!(dd.core_profiles.profiles_1d)
    cp1d.time = t
    cp1d.grid.rho_tor_norm         = ig.rho
    cp1d.electrons.density         = ig.ne .* 1e19   # 10¹⁹m⁻³ → m⁻³
    cp1d.electrons.density_thermal = ig.ne .* 1e19
    cp1d.electrons.temperature     = ig.te .* 1e3    # keV → eV
    ismissing(ig.z_eff) || (cp1d.zeff = ig.z_eff)
    ismissing(ig.w0)    || (cp1d.rotation_frequency_tor_sonic = ig.w0)

    # Core profiles — ions: two-pass merge so thermal+fast rows of the same species
    # share one IMAS ion entry (matching OMFIT to_omas behavior).
    # Pass 1: build species_map (stripped name → IMAS ion index, unique species only)
    species_map = Dict{String,Int}()
    ion_count = 0
    for k in 1:ig.nion
        name = strip(ig.name[k])
        if !haskey(species_map, name)
            ion_count += 1
            species_map[name] = ion_count
        end
    end
    resize!(cp1d.ion, ion_count)

    # Pass 2: fill ion entries (thermal and fast rows of same species share one entry)
    for k in 1:ig.nion
        name = strip(ig.name[k])
        idx = species_map[name]
        ion = cp1d.ion[idx]

        # Always set element/label (safe to repeat for identical-named species)
        resize!(ion.element, 1)
        ion.element[1].z_n = ig.z[k]
        ion.element[1].a   = ig.mass[k]
        ion.label          = ig.name[k]

        ni = ig.ni[k, :] .* 1e19   # 10¹⁹m⁻³ → m⁻³
        ti = ig.ti[k, :] .* 1e3    # keV → eV

        is_fast = !ismissing(ig.type) && lowercase(strip(ig.type[k], ['[', ']'])) == "fast"
        if is_fast
            # Reconstruct isotropic fast-ion pressure: Ti = (2p⊥+p∥)/(n·e·1e3) → p = n·Ti_eV·e/3
            p_fast = @. ni * ti * mks.e / 3.0
            ion.density_fast                = ni
            ion.pressure_fast_perpendicular = p_fast
            ion.pressure_fast_parallel      = p_fast
            # Do not set temperature for fast rows (avoids overwriting thermal temperature)
        else
            ion.density         = ni
            ion.density_thermal = ni
            ion.temperature     = ti
            ismissing(ig.vtor) || (ion.rotation_frequency_tor = ig.vtor[k, :] ./ ig.rmaj)
            ismissing(ig.vpol) || (ion.velocity.poloidal = ig.vpol[k, :])
        end
    end

    # Pressure, current profiles, conductivity — direct storage overrides dynamic expressions
    ismissing(ig.ptot) || (cp1d.pressure = ig.ptot)
    ismissing(ig.johm)     || (cp1d.j_ohmic              = ig.johm     .* 1e6)
    ismissing(ig.jbs)      || (cp1d.j_bootstrap           = ig.jbs      .* 1e6)
    if hasproperty(cp1d, :j_bootstrap_tor)
        ismissing(ig.jbstor) || (setproperty!(cp1d, :j_bootstrap_tor, ig.jbstor .* 1e6))
    end
    ismissing(ig.sigmapar) || (cp1d.conductivity_parallel = ig.sigmapar .* 1e6)

    # Core sources — inverse of update_hcd
    # Sign rule: fields using += in update_hcd are stored positive; -= fields are stored negative.
    E = 1e6   # MW/m³ → W/m³
    function make_source!(name; qe=nothing, qi=nothing, qpart=nothing, qmom=nothing, jp=nothing)
        (isnothing(qe) && isnothing(qi) && isnothing(qpart) && isnothing(qmom) && isnothing(jp)) && return
        src = resize!(dd.core_sources.source, Symbol(name); wipe=false)
        p = resize!(src.profiles_1d)
        p.time = t
        p.grid.rho_tor_norm = ig.rho
        isnothing(qe)    || (p.electrons.energy    = qe)
        isnothing(qi)    || (p.total_ion_energy     = qi)
        isnothing(qpart) || (p.electrons.particles  = qpart)
        isnothing(qmom)  || (p.momentum_tor         = qmom)
        isnothing(jp)    || (p.j_parallel           = jp)
        return src
    end

    # NBI (index 2): qbeame, qbeami, qpar_beam, qmom all use += in update_hcd
    make_source!("nbi";
        qe=ismissing(ig.qbeame) ? nothing : ig.qbeame .* E,
        qi=ismissing(ig.qbeami) ? nothing : ig.qbeami .* E,
        qpart=ismissing(ig.qpar_beam) ? nothing : ig.qpar_beam,
        qmom=ismissing(ig.qmom) ? nothing : ig.qmom,
        jp=ismissing(ig.jnb) ? nothing : ig.jnb .* 1e6)

    # EC (index 3): electron RF heating and RF current drive
    make_source!("ec";
        qe=ismissing(ig.qrfe) ? nothing : ig.qrfe .* E,
        jp=ismissing(ig.jrf) ? nothing : ig.jrf .* 1e6)

    # IC (index 5): ion RF heating
    make_source!("ic";
        qi=ismissing(ig.qrfi) ? nothing : ig.qrfi .* E)

    # Fusion (index 6): qfuse and qfusi use +=
    make_source!("fusion";
        qe=ismissing(ig.qfuse) ? nothing : ig.qfuse .* E,
        qi=ismissing(ig.qfusi) ? nothing : ig.qfusi .* E)

    # Ohmic (index 7): qohme uses +=
    make_source!("ohmic"; qe=ismissing(ig.qohme) ? nothing : ig.qohme .* E)

    # Bremsstrahlung (index 8): qbrem uses -= → store as negative
    make_source!("bremsstrahlung"; qe=ismissing(ig.qbrem) ? nothing : .-ig.qbrem .* E)

    # Synchrotron (index 9): qsync uses -= → store as negative
    make_source!("synchrotron_radiation"; qe=ismissing(ig.qsync) ? nothing : .-ig.qsync .* E)

    # Line radiation (index 10): qline uses -= → store as negative
    make_source!("line_radiation"; qe=ismissing(ig.qline) ? nothing : .-ig.qline .* E)

    # Collisional equipartition (index 11): qei uses -= for electrons, += for ions
    make_source!("collisional_equipartition";
        qe=ismissing(ig.qei) ? nothing : .-ig.qei .* E,
        qi=ismissing(ig.qei) ? nothing : ig.qei .* E)

    # Gas puff (index 108): qpar_wall uses +=
    make_source!("gas_puff"; qpart=ismissing(ig.qpar_wall) ? nothing : ig.qpar_wall)

    # Recombination (index 602): qione (electron energy) and qioni (ion energy)
    make_source!("recombination";
        qe=ismissing(ig.qione) ? nothing : ig.qione .* E,
        qi=ismissing(ig.qioni) ? nothing : ig.qioni .* E)

    # Charge exchange (index 305): qcxi (ion energy loss)
    make_source!("charge_exchange";
        qi=ismissing(ig.qcxi) ? nothing : ig.qcxi .* E)

    # Store MXH shape harmonics in profiles_2d (grid_type 57, TEQUILA convention)
    # Columns: [rmaj, zmag, rmin/rmaj, kappa, cos0..nm, asin(δ), -ζ, sin3..nm]
    shape_modes = Int[]
    for m in 0:6
        !ismissing(getfield(ig, Symbol("shape_cos$m"))) && push!(shape_modes, m)
    end
    for m in 3:6
        !ismissing(getfield(ig, Symbol("shape_sin$m"))) && push!(shape_modes, m)
    end
    if !isempty(shape_modes)
        modes = maximum(shape_modes)
        resize!(eqt.profiles_2d, 1)
        eq2d = eqt.profiles_2d[1]
        eq2d.grid.dim1 = ig.rho
        eq2d.grid.dim2 = vcat(zeros(5), Float64.(1:modes), Float64.(-(1:modes)))
        eq2d.grid_type.index = 57
        n = length(ig.rho)
        get_shape(f) = ismissing(getfield(ig, f)) ? zeros(n) : getfield(ig, f)
        mat = zeros(n, 2 * modes + 5)
        mat[:, 1]  = ig.rmaj
        mat[:, 2]  = ig.zmag
        mat[:, 3]  = ig.rmin ./ ig.rmaj
        mat[:, 4]  = ig.kappa
        for m in 0:modes
            mat[:, 5 + m] = get_shape(Symbol("shape_cos$m"))
        end
        mat[:, 5 + modes + 1] = ismissing(ig.delta) ? zeros(n) : asin.(clamp.(ig.delta, -1.0, 1.0))
        mat[:, 5 + modes + 2] = ismissing(ig.zeta)  ? zeros(n) : .-ig.zeta
        for m in 3:modes
            mat[:, 5 + modes + m] = get_shape(Symbol("shape_sin$m"))
        end
        eq2d.psi = mat
    end

    return dd
end

"""
    ion_as_dict(ion::IMAS.core_profiles__profiles_1d___ion)

Process a single ion species and return species dictionary. Performs verbatim conversion without modifying species.
"""
function ion_as_dict(ion::IMAS.core_profiles__profiles_1d___ion)
    return Dict(
        :label => ion.label,
        :z => ion.element[1].z_n,
        :mass => ion.element[1].a,
        :density_thermal => ion.density_thermal,
        :density_fast => ion.density_fast,
        :temperature => ion.temperature,
        :pressure_fast_perpendicular => ion.pressure_fast_perpendicular,
        :pressure_fast_parallel => ion.pressure_fast_parallel
    )
end
