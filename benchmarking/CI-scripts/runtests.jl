# Adapter from https://github.com/sumiya11/Groebner.jl/blob/master/benchmark/CI-scripts/runtests.jl
# GPL license

# Test for performance regressions.
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse, GitHubActions, GitHub, Random, Logging
using Test, TestSetExtensions, InteractiveUtils, PrettyTables
using Base.Threads, Statistics

const MAX_DEVIATION = 0.2
const IGNORE_SMALL = 1e-3
const SAMPLES = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 1

const dir_master = (@__DIR__) * "/run-on-master"
const dir_nightly = (@__DIR__) * "/run-on-nightly"

function runbench()
    @info "Start benchmarking.."
    @info "Using $(nthreads()) threads"
    @info "Using $SAMPLES samples"

    for i in 1:SAMPLES
        # Run benchmarks on master
        @info "Benchmarking StructuralIdentifiability.jl, master, running $dir_master"
        @time run(
            `$(Base.julia_cmd()) --startup-file=no --threads=$(nthreads()) --project=$dir_master $dir_master/run_benchmarks.jl $i`,
            wait=true
        )

        # Run benchmarks on nightly
        @info "Benchmarking StructuralIdentifiability.jl, nightly, running $dir_nightly"
        @time run(
            `$(Base.julia_cmd()) --startup-file=no --threads=$(nthreads()) --project=$dir_nightly $dir_nightly/run_benchmarks.jl $i`,
            wait=true
        )
    end
end

# Adapted from https://github.com/MakieOrg/Makie.jl/blob/v0.21.0/metrics/ttfp/run-benchmark.jl.
# License is MIT.
function best_unit(m)
    if m < 1e3
        return 1, "ns"
    elseif m < 1e6
        return 1e3, "μs"
    elseif m < 1e9
        return 1e6, "ms"
    else
        return 1e9, "s"
    end
end

function load_data()
    results = []
    for i in 1:SAMPLES
        results_master_i = nothing
        results_nightly_i = nothing
        try
            results_master_file = open(dir_master * "/results_$i", "r")
            @info "Reading $results_master_file"
            results_master_i = readlines(results_master_file)
            close(results_master_file)
        catch e
            @warn "Error when reading the file with results"
        end
        try
            results_nightly_file = open(dir_nightly * "/results_$i", "r")
            @info "Reading $results_nightly_file"
            results_nightly_i = readlines(results_nightly_file)
            close(results_nightly_file)
        catch e
            @warn "Error when reading the file with results, sample $i"
        end
        @assert !(results_master_i === nothing) && !(results_nightly_i === nothing)
        @assert length(results_master_i) == length(results_nightly_i)
        @assert !isempty(results_master_i)
        @assert results_master_i[1] == "master" && results_nightly_i[1] == "nightly"
        results_master_i = results_master_i[2:end]
        results_nightly_i = results_nightly_i[2:end]
        push!(results, (master=results_master_i, nightly=results_nightly_i))
    end
    results
end

function clean_data(results)
    nrecords = length(results[1][1])
    results_problems = Vector{Any}(undef, nrecords)
    results_types = Vector{Any}(undef, nrecords)
    results_master = [[] for _ in 1:nrecords]
    results_nightly = [[] for _ in 1:nrecords]
    for i in 1:SAMPLES
        for j in 1:nrecords
            master = results[i].master[j]
            nightly = results[i].nightly[j]
            problem_name_master, type, times_master = split(master, ":")
            problem_name_nightly, type, times_nightly = split(nightly, ":")
            @assert problem_name_master == problem_name_nightly
            times_master = map(
                x -> parse(Float64, String(strip(x, ['[', ']', ' ', '\t']))),
                split(times_master, ",")
            )
            times_nightly = map(
                x -> parse(Float64, String(strip(x, ['[', ']', ' ', '\t']))),
                split(times_nightly, ",")
            )
            append!(results_master[j], times_master)
            append!(results_nightly[j], times_nightly)
            results_problems[j] = join(split(problem_name_master, ","), " ")
            results_types[j] = type
        end
    end
    results_problems, results_types, results_master, results_nightly
end

# Compare results
function compare()
    results = load_data()
    results_problems, results_types, results_master, results_nightly = clean_data(results)
    table = Matrix{Any}(undef, length(results_master), 4)
    fail = false
    tolerance = 0.02
    for (i, (master, nightly)) in enumerate(zip(results_master, results_nightly))
        if results_types[i] == "time"
            master = 1e9 .* master
            nightly = 1e9 .* nightly
            f, unit = best_unit(maximum(master))
            m1 = round(mean(master) / f, digits=2)
            d1 = round(std(master) / f, digits=2)
            label_master = "$m1 ± $d1 $unit"
            m2 = round(mean(nightly) / f, digits=2)
            d2 = round(std(nightly) / f, digits=2)
            label_nightly = "$m2 ± $d2 $unit"
            indicator = if mean(master) < 1e9 * IGNORE_SMALL
                0, "insignificant"
            elseif (1 + MAX_DEVIATION) * m1 < m2
                fail = true
                2, "worse❌"
            elseif m1 > (1 + MAX_DEVIATION) * m2
                1, "better✅"
            else
                0, "don't care"
            end
        elseif results_types[i] == "allocs"
            label_master = mean(master)
            label_nightly = mean(nightly)
            indicator = if label_master < (1 - tolerance) * label_nightly
                fail = true
                2, "worse❌"
            elseif label_master > (1 + tolerance) * label_nightly
                1, "better✅"
            else
                0, "don't care"
            end
        else
            error("Beda!")
        end
        table[i, 1] = results_problems[i]
        table[i, 2] = label_master
        table[i, 3] = label_nightly
        table[i, 4] = indicator[2]
    end
    fail, table
end

function post(fail, table)
    comment_header = """
    ## Running times benchmark

    Note, that these numbers may fluctuate on the CI servers, so take them with a grain of salt.

    """
    io = IOBuffer()
    println(io, comment_header)
    if fail
        println(io, "Potential regressions detected❌")
    else
        println(io, "No regressions detected✅")
    end
    table_header = ["Problem", "Master", "This commit", "Result"]
    pretty_table(io, table, header=table_header, alignment=[:l, :r, :r, :r])
    comment_str = String(take!(io))
    println(comment_str)
end

function main()
    runbench()
    fail, table = compare()
    post(fail, table)
    versioninfo(verbose=true)
    @testset "Benchmarks" begin
        @test !fail
        if fail
            exit(1)
        end
    end
end

main()
