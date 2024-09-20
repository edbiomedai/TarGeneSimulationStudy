function get_colormap()
    return Dict(
        "GLMNET" => :blue,
        "XGBOOST" => :red
        )
end

function add_labels(fig)
    Label(fig[0, 1], "CV", fontsize=20, font=:bold, tellwidth=false, tellheight=true)
    Label(fig[0, 2], "Canonical", fontsize=20, font=:bold, tellwidth=false, tellheight=true)
    Label(fig[1, 0], "OSE", rotation = pi/2, fontsize=20, font=:bold, tellwidth=true, tellheight=false)
    Label(fig[2, 0], "wTMLE", rotation = pi/2, fontsize=20, font=:bold, tellwidth=true, tellheight=false)
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
                yticks=yticks
            )
            vlines!(ax, 0.95, color=:black)
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
    # Title
    Label(fig[0, :, Top()], "Coverage by Estimator", fontsize = 24, font=:bold)
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

function plot_coverage_by_positivity_across_estimators(constrained_estimands_coverages; sample_size=500_000)
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
                limits=(nothing, (ymin, nothing))
            )
            hlines!(ax, 0.95, color=:black)
            model_groups = pairs(groupby(estimator_type_group, :MODEL, sort=true))
            for (model_key, model_group) ∈ model_groups
                scatterlines!(ax,
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
            ax2 = Axis(fig[estimator_type_id, cv_id], ylabel=ylabel, yaxisposition = :right)
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
    fig = Figure(size=(1000, 800))
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