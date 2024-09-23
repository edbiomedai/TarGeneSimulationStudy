function get_colormap()
    return Dict(
        "GLMNET" => :orange,
        "XGBOOST" => :blue
        )
end

function add_labels(fig)
    Label(fig[0, 1], "CV", fontsize=30, font=:bold, tellwidth=false, tellheight=true)
    Label(fig[0, 2], "Canonical", fontsize=30, font=:bold, tellwidth=false, tellheight=true)
    Label(fig[1, 0], "OSE", rotation = pi/2, fontsize=30, font=:bold, tellwidth=true, tellheight=false)
    Label(fig[2, 0], "wTMLE", rotation = pi/2, fontsize=30, font=:bold, tellwidth=true, tellheight=false)
end

function plot_coverage_per_estimator(results)
    color_map = get_colormap()
    sample_sizes = sort(unique(results.SAMPLE_SIZE))
    yticks = (1:size(sample_sizes, 1), string.(sample_sizes))
    fig = Figure(size=(1000, 800))
    for (cv_id, (cv_key, cv_group)) ∈ enumerate(pairs(groupby(results, :CV_INFO, sort=true)))
        for (estimator_type_id, (estimator_type_key, estimator_type_group)) ∈ enumerate(pairs(groupby(cv_group, :ESTIMATOR_TYPE, sort=true)))
            xlabel = estimator_type_id == 2 ? "Coverage" : ""
            ylabel = cv_id == 1 ? "Sample Size" : ""
            ax = Axis(fig[estimator_type_id, cv_id], 
                xlabel=xlabel, 
                ylabel=ylabel, 
                yticks=yticks,
                yticklabelsvisible=cv_id==1, 
                yticksvisible=cv_id==1
            )
            vlines!(ax, 0.95, color=:green, linewidth=3)
            for (ss_id, (ss_key, ss_group)) ∈ enumerate(pairs(groupby(estimator_type_group, :SAMPLE_SIZE, sort=true)))
                for (model_key, model_group) ∈ pairs(groupby(ss_group, :MODEL))
                    hist!(ax, model_group.MEAN_COVERAGE, 
                        color=(color_map[model_key.MODEL], 0.5), 
                        offset=ss_id, 
                        normalization=:probability,
                        label=model_key.MODEL
                    )
                end
            end
        end
    end
    # Legend
    fig[end+1, :] = Legend(fig,
        current_axis(),
        "Model",
        orientation = :horizontal, 
        tellwidth = true, 
        tellheight = true,
        unique=true
        )
    add_labels(fig)
    return fig
end

function plot_coverage_by_positivity_across_estimators(constrained_estimands_coverages; sample_size=500_000, xcrop=nothing)
    constrained_estimands_coverages_ss = filter(:SAMPLE_SIZE => ==(sample_size), constrained_estimands_coverages)
    ymin = minimum(constrained_estimands_coverages.MEAN_COVERAGE_mean)
    fig = Figure(size=(1000, 800))
    color_map = get_colormap()
    for (cv_id, (cv_key, cv_group)) ∈ enumerate(pairs(groupby(constrained_estimands_coverages_ss, :CV_INFO, sort=true)))
        for (estimator_type_id, (estimator_type_key, estimator_type_group)) ∈ enumerate(pairs(groupby(cv_group, :ESTIMATOR_TYPE, sort=true)))
            ylabel = cv_id == 1 ? "Coverage" : ""
            xlabel = estimator_type_id == 2 ? "Positivity Threshold" : ""
            ax = Axis(fig[estimator_type_id, cv_id], 
                titlesize=20,
                xlabel=xlabel, 
                ylabel=ylabel,
                yticklabelsvisible=cv_id==1, 
                yticksvisible=cv_id==1,
                limits=((nothing, xcrop), (ymin, nothing))
            )
            hlines!(ax, 0.95, color=:green, linewidth=3)
            model_groups = pairs(groupby(estimator_type_group, :MODEL, sort=true))
            for (model_key, model_group) ∈ model_groups
                scatterlines!(ax,
                    markersize = 12,
                    model_group.POSITIVITY_CONSTRAINT,
                    model_group.MEAN_COVERAGE_mean,
                    color=(color_map[model_key.MODEL], 1.),
                    label=model_key.MODEL
                )
                band!(ax, 
                    model_group.POSITIVITY_CONSTRAINT, 
                    first.(model_group.CONFINT), 
                    last.(model_group.CONFINT),
                    color=(color_map[model_key.MODEL], 0.2),
                )
            end
            ylabel = cv_id == 2 ? "Estimands %" : ""
            ax2 = Axis(fig[estimator_type_id, cv_id], 
                ylabel=ylabel, 
                yaxisposition = :right,
                yticklabelsvisible=cv_id==2, 
                yticksvisible=cv_id==2,
                limits=((nothing, xcrop), nothing)
                )
            hidespines!(ax2)
            hidexdecorations!(ax2)
            model_group = first(model_groups)[2]
            lines!(ax2, 
                model_group.POSITIVITY_CONSTRAINT,
                model_group.ESTIMANDS_RATIO,
                linestyle=:dash,
                color=:black
            )
        end
    end
    fig[end+1, :] = Legend(fig,
        fig.content[1],
        "Model",
        orientation = :horizontal, 
        tellwidth = true, 
        tellheight = true,
        unique=true
        )
    add_labels(fig)
    return fig
end

function plot_loss_relative_difference(density_estimates_prefix)
    # Retrieve Results
    density_results = load_density_results(density_estimates_prefix)
    sort!(density_results, order(:RELATIVE_IMPROVEMENT, by=x->abs(x)))
    xs_dict = Dict(outcome => index for (index, outcome) ∈ enumerate(unique(density_results.OUTCOME)))
    density_results.OUTCOME_ID = [xs_dict[outcome] for outcome ∈ density_results.OUTCOME]
    # Plot
    colors = Makie.wong_colors()
    fig = Figure(size=(1000, 1000))
    xs_to_labels = sort(unique(DataFrames.select(density_results, [:OUTCOME_ID, :OUTCOME])), :OUTCOME_ID)
    ax = Axis(fig[1, 1], 
        title="Relative Loss Improvement of SNNE over GLM", 
        ylabel="Density", 
        xlabel="%", 
        yticks=(xs_to_labels.OUTCOME_ID, xs_to_labels.OUTCOME)
    )
    barplot!(ax,
        density_results.OUTCOME_ID, 
        density_results.RELATIVE_IMPROVEMENT,
        dodge = density_results.TYPE,
        color = colors[density_results.TYPE],
        direction=:x,
    )
    vlines!(0, color=:black)
    elements = [PolyElement(polycolor = colors[1]), PolyElement(polycolor = colors[2])]
    Legend(fig[2,:], elements, ["Train", "Validation"], "Set",
        orientation = :horizontal, tellwidth = false, tellheight = true)
    return fig
end

function plot_power_cv_canonical_comparison(results;
    model = "XGBOOST",
    estimator_type = "wTMLE",
    sample_size = 500_000)

    results_subset = filter(
        x -> x.MODEL == model && x.ESTIMATOR_TYPE == estimator_type && x.SAMPLE_SIZE == sample_size, 
        results
    )
    cv_results = filter(:CV_INFO => ==("CV"), results_subset)
    canonical_results = filter(:CV_INFO => ==("Canonical"), results_subset)
    joined = innerjoin(
        select(cv_results, 
            :ESTIMAND_TYPE => ByRow(x -> string(x)) => :ESTIMAND_TYPE, 
            :TRUE_EFFECT, 
            :TRUE_EFFECT => ByRow(x -> minimum(abs.(x))) => :INF_EFFECT_NORM,
            :TRUE_EFFECT => ByRow(x -> sum(x.^2)) => :EFFECT_SQUARED_NORM, 
            :ESTIMAND, 
            :ESTIMATES => :CV_ESTIMATES,
            :POWER => :CV_POWER,
            :MEAN_VARIANCE => :CV_MEAN_VARIANCE,
            [:MEAN_VARIANCE, :TRUE_EFFECT] => ByRow((x, y) -> sqrt(sum(y.^2) / x)) => :CV_REL_STD,
            ),
        select(canonical_results, 
            :ESTIMAND, 
            :POWER => :CANONICAL_POWER,
            :MEAN_VARIANCE => :CANONICAL_MEAN_VARIANCE,
            [:MEAN_VARIANCE, :TRUE_EFFECT] => ByRow((x, y) -> sqrt(sum(y.^2) / x)) => :CANONICAL_REL_STD,
            ),
        on=:ESTIMAND
    )
    sort!(joined, [:ESTIMAND_TYPE, :CV_POWER])

    fig = Figure()
    border = findfirst(x -> x == "TMLE.StatisticalATE", joined.ESTIMAND_TYPE) - 0.5
    ax = Axis(fig[1, 1], ylabel="Power", xticklabelsvisible=false, xticksvisible=false, xgridvisible=false, ygridvisible=false)
    scatter!(ax, joined.CV_POWER, color=(:blue, 0.5), label="CV")
    scatter!(ax, joined.CANONICAL_POWER, color=(:orange, 0.5), label="Canonical")
    vlines!(ax, border, color=:black, linestyle=:dash)
    ax = Axis(fig[2, 1], xlabel="Estimands", ylabel=L"\frac{||\Psi_0||_2}{\sqrt{tr(\Sigma)}}", xticklabelsvisible=false, xticksvisible=false, xgridvisible=false, ygridvisible=false)
    scatter!(ax, joined.CV_REL_STD, color=(:blue, 0.5), label="CV")
    scatter!(ax, joined.CANONICAL_REL_STD, color=(:orange, 0.5), label="Canonical")
    vlines!(ax, border, color=:black, linestyle=:dash)
    fig[3, :] = Legend(fig, ax, orientation=:horizontal, tellheoght=true, tellwidth=false)
    return fig
end

function power_plot(results; 
    estimator_type = "wTMLE", 
    sample_size_labels = Dict(
        500_000 => "Sample Size\n500 000",
        50_000 => "Sample Size\n50 000"
    )
    )
    results_subset = filter(
        x -> x.ESTIMATOR_TYPE == estimator_type, 
        results
    )
    fig = Figure(size=(1000, 800))
    for (ax_row, (sample_size_key, sample_size_group)) in enumerate(pairs(groupby(results_subset, :SAMPLE_SIZE, sort=true)))
        for (ax_col, (model_key, model_group)) in enumerate(pairs(groupby(sample_size_group, :MODEL, sort=true)))
            ylabel = ax_col == 1 ? "Power" : ""
            ax = Axis(fig[ax_row, ax_col], 
                ylabel=ylabel, 
                limits=(nothing, (0, 1.05)),
                xticklabelsvisible=false, 
                xticksvisible=false, 
                xgridvisible=false, 
                ygridvisible=false
            )
            cv_results = filter(:CV_INFO => ==("CV"), model_group)
            canonical_results = filter(:CV_INFO => ==("Canonical"), model_group)
            joined = innerjoin(
                select(cv_results, 
                    :ESTIMAND_TYPE => ByRow(x -> string(x)) => :ESTIMAND_TYPE, 
                    :ESTIMAND, 
                    :POWER => :CV_POWER,
                    ),
                select(canonical_results, 
                    :ESTIMAND, 
                    :POWER => :CANONICAL_POWER,
                    ),
                on=:ESTIMAND
            )
            sort!(joined, [:ESTIMAND_TYPE, :CV_POWER])
            border = findfirst(x -> x == "TMLE.StatisticalATE", joined.ESTIMAND_TYPE) - 0.5
            scatter!(ax, joined.CV_POWER, color=(:blue, 0.5), label="CV", markersize=12)
            scatter!(ax, joined.CANONICAL_POWER, color=(:orange, 0.5), label="Canonical", markersize=12)
            vlines!(ax, border, color=:black, linestyle=:dash)
            ax_row == 1 && Label(fig[0, ax_col], string(model_key.MODEL), tellwidth=false)
            ax_col == 1 && Label(fig[ax_row, 0], sample_size_labels[sample_size_key.SAMPLE_SIZE], tellheight=false, rotation=π/2)
        end
    end
    
    fig[3, :] = Legend(fig, fig.current_axis[], orientation=:horizontal, tellheight=true, tellwidth=false)
    return fig
end

