#!/usr/bin/env julia
# perf_comp.jl (perform comparison)
# Unlike param_sweep.jl, which sweeps through a range of parameter values, this
#  script runs the sim multiple times for a small set of different settings, and
#  produces comparison statistics and plots.

using Revise
using RCall
@rlibrary DescTools
using DataFrames
using CSV
using Bootstrap
include("sim.jl")
include("setup_params.jl")


num_trials=40
original_seed=params[:random_seed]

settings_sets = [
    Dict(:label => "empty", :whichGraph => "empty"),
    Dict(:label => "λ=0", :whichGraph => "erdos_renyi", :λ => 0.0),
#   Dict(:label => "λ=¼", :whichGraph => "erdos_renyi", :λ => 0.25),
    Dict(:label => "λ=½", :whichGraph => "erdos_renyi", :λ => 0.5),
#   Dict(:label => "λ=¾", :whichGraph => "erdos_renyi", :λ => 0.75),
    Dict(:label => "λ=1", :whichGraph => "erdos_renyi", :λ => 1.0),
#   Dict(:label => "λ=1¼", :whichGraph => "erdos_renyi", :λ => 1.25),
    Dict(:label => "λ=1½", :whichGraph => "erdos_renyi", :λ => 1.5),
#   Dict(:label => "λ=1¾", :whichGraph => "erdos_renyi", :λ => 1.75),
    Dict(:label => "λ=2", :whichGraph => "erdos_renyi", :λ => 2.0),
    Dict(:label => "λ=3", :whichGraph => "erdos_renyi", :λ => 3.0),
    Dict(:label => "λ=4", :whichGraph => "erdos_renyi", :λ => 4.0),
    Dict(:label => "λ=5", :whichGraph => "erdos_renyi", :λ => 5.0),
    Dict(:label => "λ=6", :whichGraph => "erdos_renyi", :λ => 6.0),
    Dict(:label => "complete", :whichGraph => "complete"),
]

function perform_comparison(settings_sets)

    prc("Starting comparison...\n")

    trials_df = DataFrame(
        label=String[],
        seed=Int[],
        gini_pre=Float64[],
        gini_post=Float64[],
        num_comps=Float64[],
        avg_proto_size_pre=Float64[],
        avg_proto_size_post=Float64[],
        num_agents_pre=Float64[],
        num_agents_post=Float64[],
        num_protos_pre=Float64[],
        num_protos_post=Float64[],
    )

    for settings in settings_sets

        # refresh default parameters before overriding for next settings
        include("setup_params.jl")
        @assert all_parameters_legit(settings)
        merge!(params, Dict(settings))
        params[:make_anims] = false  # never true for a comparison
        params[:make_sim_plots] = false  # never true for a comparison

        prc("Running settings $(settings)...\n")
        for j = 1:num_trials

            prc(".")
            if j % 10 == 0 prc("$(j)\n") end

            # Actually run the simulation!
            results=specnet()
            overall_results = results[:overall_results]
            starvation_results = results[:starvation_results]

            push!(trials_df, (
                settings[:label],
                params[:random_seed],
                overall_results[:gini],
                starvation_results[:gini],
                overall_results[:num_comps],
                overall_results[:average_proto_size],
                starvation_results[:average_proto_size],
                overall_results[:num_living_agents],
                starvation_results[:num_living_agents],
                overall_results[:num_protos],
                starvation_results[:num_protos],
            ))

            params[:random_seed]+=1
        end
        params[:random_seed]=original_seed
    end

    for var in setdiff(names(trials_df), [:label,:seed])
        plot_df = DataFrame(
            label=String[],
            value=Float64[],
            high=Float64[],
            low=Float64[]
        )
        for label ∈  unique(trials_df[:label])
            # Compute averages and confidence intervals
            bs = bootstrap(mean, trials_df[trials_df[:label] .== label, var],
                BasicSampling(params[:num_boot_samples]))
            ci = confint(bs, BasicConfInt(.95))[1]
            push!(plot_df, (label, ci[1], ci[2], ci[3]))
            draw_plot(plot_df, var)
global t = plot_df
        end
    end
end

# draw_plot() -- create, write to file, and return a comparison boxplot.
function draw_plot(plot_df, var_name)
    p = plot(plot_df,
        x=:label,
        y=:value,
        ymax=:high,
        ymin=:low,
        Geom.errorbar, Geom.point,
        Guide.xlabel("Graph type"),
        Guide.ylabel(string(var_name)),
        Theme(default_color="navy")
    );
    draw(PNG("$(tempdir())/$(var_name).png"), p)
end

perform_comparison(settings_sets)
