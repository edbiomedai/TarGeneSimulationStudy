using JLD2
using DataFrames
using TMLE
using TMLECLI
using OrderedCollections
using CairoMakie
using Makie
using Arrow
using Simulations

include("src/wrangling.jl")
include("src/plots.jl")

function main()
    # Null Simulation
    ## Files
    null_dataset_file = joinpath("null_simulation", "results", "datasets", "all_genotypes.data.arrow")
    null_results_file = joinpath("null_simulation", "results", "null_simulation_results.hdf5")
    ## Data Loading
    null_dataset = Arrow.Table(null_dataset_file) |> DataFrame
    null_results = jldopen(io -> io["results"], null_results_file)
    ## Plots
    add_estimator_info_cols!(null_results)
    plot_coverage_per_estimator(null_results)
    constrained_estimands_coverages = get_constrained_estimands_coverages(null_results, null_dataset)
    plot_coverage_by_positivity_across_estimators(constrained_estimands_coverages; sample_size=500_000)
    plot_coverage_by_positivity_across_estimators(constrained_estimands_coverages; sample_size=50_000)

    # Realistic Simulation
    ## Files
    realistic_dataset_file = joinpath("null_simulation", "results", "datasets", "all_genotypes.data.arrow")
    realistic_results_file = joinpath("null_simulation", "results", "null_simulation_results.hdf5")
    ## Data Loading
    realistic_dataset = Arrow.Table(realistic_dataset_file) |> DataFrame
    realistic_results = jldopen(io -> io["results"], realistic_results_file)
    ## Plots
    add_estimator_info_cols!(realistic_results)
    plot_coverage_per_estimator(realistic_results)
    constrained_estimands_coverages = get_constrained_estimands_coverages(realistic_results, realistic_dataset)
    plot_coverage_by_positivity_across_estimators(constrained_estimands_coverages; sample_size=500_000)
    plot_coverage_by_positivity_across_estimators(constrained_estimands_coverages; sample_size=50_000)

end