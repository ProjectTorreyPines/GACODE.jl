using GACODE
using GACODE.IMAS
using Test

const SAMPLE = joinpath(@__DIR__, "input.gacode")
const ALL_FIELDS = Set(fieldnames(InputGACODE))
const ALLOWED_IGNORES = Set([:time])
const HAS_JBSTOR_FIELD = let dd = IMAS.dd()
    cp1d = resize!(dd.core_profiles.profiles_1d)
    hasproperty(cp1d, :j_bootstrap_tor)
end

is_fast_t(t) = lowercase(strip(string(t), ['[', ']', ' '])) == "fast"

function ion_permutation(ig1, ig2)
    p = zeros(Int, ig1.nion)
    used = falses(ig2.nion)
    for i in 1:ig1.nion
        n1 = strip(ig1.name[i])
        f1 = is_fast_t(ig1.type[i])
        for j in 1:ig2.nion
            used[j] && continue
            if strip(ig2.name[j]) == n1 && is_fast_t(ig2.type[j]) == f1
                p[i] = j
                used[j] = true
                break
            end
        end
    end
    return p
end

function assert_values_equal(v1, v2)
    if v1 isa AbstractArray
        if eltype(v1) <: Number
            @test v1 ≈ v2
        else
            @test v1 == v2
        end
    elseif v1 isa Number
        @test v1 ≈ v2
    else
        @test v1 == v2
    end
end

function assert_inputgacode_equal(ig1, ig2; check_missing=true, ignore_fields=Symbol[])
    ignore_set = Set(ignore_fields)
    @test ignore_set ⊆ ALLOWED_IGNORES
    @test isempty(setdiff(ignore_set, ALL_FIELDS))

    perm = ion_permutation(ig1, ig2)
    @test all(!iszero, perm)

    ion_vector_fields = Set((:name, :type, :mass, :z))
    ion_matrix_fields = Set((:ni, :ti, :vpol, :vtor))
    checked_fields = Set{Symbol}()

    for field in fieldnames(InputGACODE)
        field in ignore_set && continue
        push!(checked_fields, field)
        v1 = getfield(ig1, field)
        v2 = getfield(ig2, field)

        if check_missing
            @test ismissing(v1) == ismissing(v2)
        end
        (ismissing(v1) || ismissing(v2)) && continue

        if field in ion_vector_fields
            v2 = v2[perm]
        elseif field in ion_matrix_fields
            v2 = v2[perm, :]
        end
        assert_values_equal(v1, v2)
    end

    @test union(checked_fields, ignore_set) == ALL_FIELDS
end

function synthetic_input_gacode(; nexp=11)
    rho = collect(range(0.0, 1.0; length=nexp))

    ig = GACODE.InputGACODE()
    ig.nexp = nexp
    ig.nion = 3
    ig.shot = 123456
    ig.time = 789
    ig.bcentr = 2.2
    ig.current = -1.35
    ig.rcentr = 1.68
    ig.torfluxa = 0.71

    ig.name = ["D", "C", "N"]
    ig.type = ["therm", "therm", "therm"]
    ig.mass = [2.014, 12.011, 14.007]
    ig.z = [1.0, 6.0, 7.0]

    ig.rho = rho
    ig.polflux = -0.32 .* rho .^ 2
    ig.rmin = 0.58 .* rho
    ig.rmaj = 1.68 .+ 0.08 .* rho
    ig.kappa = 1.55 .+ 0.22 .* rho
    ig.delta = 0.20 .* rho
    ig.zeta = 0.05 .* rho
    ig.zmag = 0.01 .* (1 .- rho)
    ig.q = 1.2 .+ 3.5 .* rho .^ 2
    ig.fpol = 2.1 .+ 0.12 .* rho

    ig.shape_cos0 = 0.01 .+ 0.02 .* rho
    ig.shape_cos1 = -0.04 .* rho
    ig.shape_cos2 = 0.03 .* rho .^ 2
    ig.shape_cos3 = -0.02 .* rho .^ 3
    ig.shape_cos4 = 0.015 .* rho
    ig.shape_cos5 = -0.012 .* rho .^ 2
    ig.shape_cos6 = 0.010 .* rho .^ 3
    ig.shape_sin3 = 0.018 .* rho
    ig.shape_sin4 = -0.016 .* rho .^ 2
    ig.shape_sin5 = 0.014 .* rho .^ 3
    ig.shape_sin6 = -0.010 .* rho .^ 4

    ig.ne = 5.0 .* (1 .- 0.7 .* rho .^ 2) .+ 0.1
    ig.te = 7.0 .* (1 .- 0.8 .* rho .^ 2) .+ 0.2
    ig.z_eff = 2.0 .+ 0.3 .* rho
    ig.w0 = 4.0e4 .* (1 .- 0.5 .* rho)

    ig.ni = zeros(ig.nion, nexp)
    ig.ti = zeros(ig.nion, nexp)
    ig.vtor = zeros(ig.nion, nexp)
    ig.vpol = zeros(ig.nion, nexp)
    for k in 1:ig.nion
        ig.ni[k, :] = (2.6 - 0.4 * k) .* (1 .- 0.65 .* rho .^ 2) .+ 0.05
        ig.ti[k, :] = (6.0 - 0.8 * k) .* (1 .- 0.75 .* rho .^ 2) .+ 0.1
        ig.vtor[k, :] = (1.0e5 + 2.0e4 * k) .* (1 .- 0.6 .* rho)
        ig.vpol[k, :] = (4.0e4 + 8.0e3 * k) .* (1 .- 0.5 .* rho .^ 2)
    end

    ig.ptot = 8.0e4 .* (1 .- 0.75 .* rho .^ 2) .+ 2.0e3

    prof(scale, off=0.0) = scale .* (1 .- 0.72 .* rho .^ 2) .+ off

    ig.johm = prof(0.6)
    ig.jbs = prof(0.4)
    ig.jrf = prof(0.2)
    ig.jnb = prof(0.3)
    ig.jbstor = prof(0.35)
    ig.sigmapar = prof(1.4, 0.2)

    ig.qohme = prof(0.04)
    ig.qbeame = prof(0.08)
    ig.qbeami = prof(0.16)
    ig.qrfe = prof(0.05)
    ig.qrfi = prof(0.03)
    ig.qfuse = prof(0.07)
    ig.qfusi = prof(0.06)
    ig.qbrem = prof(0.01)
    ig.qsync = prof(3e-4)
    ig.qline = prof(0.004)
    ig.qei = @. 0.03 * (1 - 1.2 * rho^2)
    ig.qione = @. -0.02 * (1 - 0.5 * rho^2)
    ig.qioni = @. 0.01 * (1 - 0.4 * rho^2)
    ig.qcxi = @. -0.015 * (1 - 0.6 * rho^2)
    ig.qpar_beam = @. 1.3e19 * (1 - 0.65 * rho^2) + 1.0e17
    ig.qpar_wall = @. 0.8e19 * (1 - 0.50 * rho^2) + 2.0e17
    ig.qmom = @. -0.15 * (1 - 0.55 * rho^2)

    return ig
end

@testset "GACODE full-field assertion policy" begin
    @test ALLOWED_IGNORES == Set([:time])
    @test HAS_JBSTOR_FIELD
end

@testset "GACODE round-trip: file -> dd -> InputGACODE" begin
    ig1 = GACODE.load(SAMPLE)
    missing_non_time = [f for f in fieldnames(InputGACODE) if f ∉ ALLOWED_IGNORES && ismissing(getfield(ig1, f))]
    @test !isempty(missing_non_time)

    dd = IMAS.dd()
    GACODE.dd!(dd, ig1)
    ig2 = GACODE.InputGACODE(dd)
    assert_inputgacode_equal(ig1, ig2; ignore_fields=collect(ALLOWED_IGNORES))
end

@testset "GACODE round-trip: synthetic full-field coverage" begin
    ig1 = synthetic_input_gacode()
    non_missing_fields = [f for f in fieldnames(InputGACODE) if f != :header_lines]
    @test all(!ismissing(getfield(ig1, f)) for f in non_missing_fields)

    dd = IMAS.dd(ig1)
    ig2 = GACODE.InputGACODE(dd)

    assert_inputgacode_equal(ig1, ig2)
end

@testset "GACODE round-trip: file-level save/load after conversion" begin
    ig1 = GACODE.load(SAMPLE)
    dd = IMAS.dd(ig1)
    ig2 = GACODE.InputGACODE(dd)

    tmp = tempname()
    GACODE.save(ig2, tmp)
    ig3 = GACODE.load(tmp)

    assert_inputgacode_equal(ig2, ig3)
end

@testset "GACODE source guards: partial recombination/charge-exchange channels" begin
    ig = synthetic_input_gacode()
    ig.qioni = missing
    ig.qcxi = missing

    dd = IMAS.dd(ig)

    # Create a charge_exchange source with empty profiles_1d to exercise isempty guards.
    resize!(dd.core_sources.source, :charge_exchange; wipe=false)

    ig2 = @test_nowarn GACODE.InputGACODE(dd)
    @test ismissing(ig2.qioni)
    @test ismissing(ig2.qcxi)
end

@testset "GACODE jbstor guard when expression cannot evaluate" begin
    ig = synthetic_input_gacode()
    ig.jbstor = missing

    dd = IMAS.dd(ig)

    ig2 = @test_nowarn GACODE.InputGACODE(dd)
    @test ismissing(ig2.jbstor) || length(ig2.jbstor) == ig.nexp
end
