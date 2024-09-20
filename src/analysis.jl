using JLD2
using DataFrames
using TMLE
using TMLECLI
using OrderedCollections
using CairoMakie
using Makie
using Arrow
using Simulations
using TargeneCore

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
    fig = plot_coverage_per_estimator(null_results)
    save("coverage_per_estimator_null.png", fig)
    null_constrained_estimands_coverages = get_constrained_estimands_coverages(null_results, null_dataset)
    fig = plot_coverage_by_positivity_across_estimators(null_constrained_estimands_coverages; sample_size=500_000)
    save("coverage_by_positivity_null_500000.png", fig)
    fig = plot_coverage_by_positivity_across_estimators(null_constrained_estimands_coverages; sample_size=50_000)
    save("coverage_by_positivity_null_50000.png", fig)

    # Realistic Simulation
    ## Files
    realistic_dataset_file = joinpath("realistic_simulation", "results", "realistic_simulation_inputs", "ga_sim_input.data.arrow")
    realistic_results_file = joinpath("realistic_simulation", "results", "realistic_simulation_results.hdf5")
    density_estimates_prefix = joinpath("realistic_simulation", "results", "density_estimates", "ga_sim_input.conditional_density_estimate")
    ## Data Loading
    realistic_dataset = Arrow.Table(realistic_dataset_file) |> DataFrame
    realistic_results = jldopen(io -> io["results"], realistic_results_file)
    ## Plots
    fig = plot_loss_relative_difference(density_estimates_prefix)
    save("density_estimates.png", fig)
    add_estimator_info_cols!(realistic_results)
    fig = plot_coverage_per_estimator(realistic_results)
    save("coverage_per_estimator_realistic.png", fig)
    realistic_constrained_estimands_coverages = get_constrained_estimands_coverages(realistic_results, realistic_dataset)
    fig = plot_coverage_by_positivity_across_estimators(realistic_constrained_estimands_coverages; sample_size=500_000)
    save("coverage_by_positivity_realistic_500000.png", fig)
    fig = plot_coverage_by_positivity_across_estimators(realistic_constrained_estimands_coverages; sample_size=50_000)
    save("coverage_by_positivity_realistic_50000.png", fig)
end