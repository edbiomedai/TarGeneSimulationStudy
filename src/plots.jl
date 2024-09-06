function colors_and_color_map(models)
    color_map = OrderedDict(model => index for (index, model) ∈ enumerate(sort(unique(models))))
    colors = Makie.Colors.distinguishable_colors(maximum(values(color_map)), [Makie.RGB(1,1,1), Makie.RGB(0,0,0)], dropseed=true)
    return colors, color_map
end

function add_labels(fig)
    Label(fig[0, 1], "CV", fontsize=20, font=:bold, tellwidth=false, tellheight=true)
    Label(fig[0, 2], "Canonical", fontsize=20, font=:bold, tellwidth=false, tellheight=true)
    Label(fig[1, 0], "OSE", rotation = pi/2, fontsize=20, font=:bold, tellwidth=true, tellheight=false)
    Label(fig[2, 0], "wTMLE", rotation = pi/2, fontsize=20, font=:bold, tellwidth=true, tellheight=false)
end

function plot_coverage_per_estimator(results)
    colors, color_map = colors_and_color_map(results.MODEL)
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
                        color=(colors[color_map[model_key.MODEL]], 0.5), 
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
    fig = Figure(size=(1000, 800))
    colors, color_map = colors_and_color_map(constrained_estimands_coverages.MODEL)
    for (cv_id, (cv_key, cv_group)) ∈ enumerate(pairs(groupby(constrained_estimands_coverages_ss, :CV_INFO, sort=true)))
        for (estimator_type_id, (estimator_type_key, estimator_type_group)) ∈ enumerate(pairs(groupby(cv_group, :ESTIMATOR_TYPE, sort=true)))
            ylabel = cv_id == 1 ? "Coverage" : ""
            xlabel = estimator_type_id == 2 ? "Positivity Threshold" : ""
            ax = Axis(fig[estimator_type_id, cv_id], 
                titlesize=20,
                xlabel=xlabel, 
                ylabel=ylabel
            )
            hlines!(ax, 0.95, color=:black)
            model_groups = pairs(groupby(estimator_type_group, :MODEL, sort=true))
            for (model_key, model_group) ∈ model_groups
                scatterlines!(ax,
                    model_group.POSITIVITY_CONSTRAINT,
                    model_group.MEAN_COVERAGE_mean,
                    color=(colors[color_map[model_key.MODEL]], 0.5),
                    label=model_key.MODEL
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