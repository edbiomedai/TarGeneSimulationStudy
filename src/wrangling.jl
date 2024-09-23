
add_power_col!(results) = 
results.POWER = [1 - Simulations.mean_coverage(Ψ̂s, zero(TMLE.estimate(first(Ψ̂s)))) for Ψ̂s in results.ESTIMATES]

add_estimand_type_col!(results) =     
    results.ESTIMAND_TYPE = [infer_estimand_type(Ψ) for Ψ ∈ results.ESTIMAND]

add_model!(results) = results.MODEL = [split(estimator, "_")[end] for estimator in string.(results.ESTIMATOR)]

add_cv_info!(results) = results.CV_INFO = [startswith(estimator, "CV") ? "CV" : "Canonical" for estimator in string.(results.ESTIMATOR)]

infer_estimand_type(Ψ) = typeof(Ψ)

function infer_estimand_type(Ψ::TMLE.JointEstimand)
    infered_type = infer_estimand_type(first(Ψ.args))
    @assert all(infer_estimand_type(arg) == infered_type for arg ∈ Ψ.args)
    return infered_type
end

add_estimator_type!(results) = results.ESTIMATOR_TYPE = [occursin("OSE", estimator) ? "OSE" : "wTMLE" for estimator in string.(results.ESTIMATOR)]

function add_estimator_info_cols!(results)
    add_model!(results)
    add_cv_info!(results)
    add_estimator_type!(results)
end

function satisfies_positivity(Ψ, freq_table; positivity_threshold=0.01)
    for treatment_value ∈ keys(TMLE.indicator_fns(Ψ))
        if freq_table[treatment_value] < positivity_threshold
            return false
        end
    end
    return true
end

function get_frequency_table(columns, dataset)
    return DataFrames.combine(groupby(dataset, collect(columns), skipmissing=true), proprow)
end

function get_treatment_frequency_map(estimands, dataset)
    treatments_frequencies = Dict()
    for treatment in unique(Simulations.get_treatments(Ψ) for Ψ in estimands)
        treatment_frequencies = Dict()
        for row in eachrow(get_frequency_table(treatment, dataset))
            vals = values(row)
            treatment_frequencies[vals[1:end-1]] = vals[end]
        end
        treatments_frequencies[treatment] = treatment_frequencies
    end
    return treatments_frequencies
end

function positivity_respecting_estimands_map(estimands, treatments_frequencies;positivity_threshold=0.01)
    estimand_map = Dict()
    for Ψ ∈ estimands
        filtered_Ψᵢs_and_indices = filter(collect(enumerate(Ψ.args))) do (i, Ψᵢ)
            treatments = Simulations.get_treatments(Ψᵢ)
            freq_table = treatments_frequencies[treatments]
            satisfies_positivity(Ψᵢ, freq_table; positivity_threshold=positivity_threshold)
        end
        if isempty(filtered_Ψᵢs_and_indices)
            @info(string("Dropping a complete estimand"))
            continue
        end
        estimand_map[Ψ] = (
            indices = first.(filtered_Ψᵢs_and_indices), 
            new_estimand=JointEstimand(args=Tuple(last.(filtered_Ψᵢs_and_indices)))
        )
    end
    return estimand_map
end

function positivity_constrained(results, treatments_frequencies; positivity_threshold=0.01)
    unique_estimands = unique(results.ESTIMAND)

    # Creates a map to store new estimands as well as indices or remaining sub estimands for each original estimand
    estimand_map = positivity_respecting_estimands_map(unique_estimands, treatments_frequencies;
        positivity_threshold=positivity_threshold
    )
    # Build filtered results
    constrained_estimates = Any[]
    constrained_estimands = Any[]
    constrained_true_effects = Any[]
    constrained_indices = Int[]
    for (index, (Ψ, Ψ̂s, Ψ₀)) ∈ enumerate(zip(results.ESTIMAND, results.ESTIMATES, results.TRUE_EFFECT))
        if haskey(estimand_map, Ψ)
            indices, new_estimand = estimand_map[Ψ]
            newΨ̂s = map(Ψ̂s) do Ψ̂
                new_estimates = collect(Ψ̂.estimates)[indices]
                new_cov = Ψ̂.cov[indices, indices]
                n = Ψ̂.n
                TMLE.JointEstimate(;estimand=new_estimand, estimates=new_estimates, cov=new_cov, n=n)
            end
            push!(constrained_estimates, newΨ̂s)
            push!(constrained_estimands, new_estimand)
            push!(constrained_true_effects, Ψ₀[indices])
            push!(constrained_indices, index)
        end
    end
    constrained_results = DataFrame(
        ESTIMAND  = constrained_estimands,
        ESTIMATES = constrained_estimates,
        ESTIMATOR = results.ESTIMATOR[constrained_indices],
        MODEL     = results.MODEL[constrained_indices],
        CV_INFO   = results.CV_INFO[constrained_indices],
        ESTIMATOR_TYPE = results.ESTIMATOR_TYPE[constrained_indices],
        SAMPLE_SIZE = results.SAMPLE_SIZE[constrained_indices],
        TRUE_EFFECT = constrained_true_effects
    )

    Simulations.add_mean_coverage_col!(constrained_results)
    return constrained_results
end

function get_constrained_estimands_coverages(results, dataset)
    # Compute the frequency table for each treatment in the estimands
    unique_estimands = unique(results.ESTIMAND)
    treatments_frequencies = get_treatment_frequency_map(unique_estimands, dataset)
    # Constrain estimands to respect positivity and recompute statistics
    mean_coverages = []
    for positivity_threshold in (0, 0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.02, 0.03, 0.05)
        constrained_results = positivity_constrained(results, treatments_frequencies; positivity_threshold=positivity_threshold)
        mean_coverages_ = combine(groupby(constrained_results, [:SAMPLE_SIZE, :ESTIMATOR]), :MEAN_COVERAGE => mean)
        n_unique_estimands = sum(length(Ψ.args) for Ψ ∈ unique(constrained_results.ESTIMAND))
        mean_coverages_.POSITIVITY_CONSTRAINT .= positivity_threshold
        mean_coverages_.N_UNIQUE_ESTIMANDS .= n_unique_estimands
        push!(mean_coverages, mean_coverages_)
    end
    mean_coverages = vcat(mean_coverages...)
    add_estimator_info_cols!(mean_coverages)
    mean_coverages.ESTIMANDS_RATIO = mean_coverages.N_UNIQUE_ESTIMANDS ./ maximum(mean_coverages.N_UNIQUE_ESTIMANDS)
    mean_coverages.CONFINT = [confint(OneSampleTTest(x, x*(1-x), 500)) for x in mean_coverages.MEAN_COVERAGE_mean]
    return mean_coverages
end

function load_density_results(density_estimates_prefix)
    density_results = DataFrame()
    for (distribution_id, file) ∈ enumerate(TargeneCore.files_matching_prefix(density_estimates_prefix))
        jldopen(file) do io
            outcome = string(io["outcome"])
            losses = io["metrics"]
            estimators = [string(typeof(x)) for x in io["estimators"]]
            @assert estimators[1] == "SieveNeuralNetworkEstimator"
            relative_train_improvement = -100(losses[1].train_loss - losses[2].train_loss) / losses[2].train_loss
            relative_test_improvement = -100(losses[1].test_loss - losses[2].test_loss) / losses[2].test_loss
            push!(density_results, (OUTCOME_ID=distribution_id, OUTCOME=outcome, RELATIVE_IMPROVEMENT=relative_train_improvement, TYPE=1))
            push!(density_results, (OUTCOME_ID=distribution_id, OUTCOME=outcome, RELATIVE_IMPROVEMENT=relative_test_improvement, TYPE=2))
        end
    end
    return density_results
end