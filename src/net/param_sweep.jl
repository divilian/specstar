#!/usr/bin/env julia
using Revise
using RCall
@rlibrary DescTools
using DataFrames
using CSV
using Bootstrap
include("sim.jl")
include("setup_params.jl")


param_to_sweep=:salary   #parameter to iterate over *any parameter*
                     #for graph sweep param_to_sweep should not be exclusive to one graph type e.g. SF_prob
start_value=10      #value to begin sweep
end_value=200         #value to end sweep
num_values=50        #number of distinct values to run
trials_per_value=4   #for each distinct value, number of independent sims to run
graph_sweep=false    #run the sweep once for each graph type
original_seed=params[:random_seed]

components=[[],[]]
global comp_df=DataFrame(size_largest_comp=Int[],num_comps=Int[],sim_tag=Int[])


function param_sweeper(graph_name; additional_params...)

    @assert all_parameters_legit(additional_params)
    merge!(params, Dict(additional_params))

    println("Starting sweep..")
    print("Sweeping for: $(param_to_sweep)")
    counter=start_value

    global agent_line_df = DataFrame(
        replace_this=Float64[],
        seed=Int[],
        agent=String[],
        sugar=Float64[],
        proto_id=Int[],
        sim_tag=Int[]
    )
    names!(agent_line_df,
        prepend!(names(agent_line_df)[2:end], [param_to_sweep]))

    global trial_line_df=DataFrame(
        replace_this=Float64[],
        seed=Int[],
        gini=Float64[],
        sim_tag=Int[]
    )
    names!(trial_line_df,
        prepend!(names(trial_line_df)[2:end], [param_to_sweep]))

    params[:make_anims] = false  # We would never want this true for a sweep
    params[:make_sim_plots] = false  # We would never want this true for a sweep

    for i = 1:num_values
        for j = 1:trials_per_value
            #setting the random seed and adding it to the DataFrame of the final gini of each simulation
            Random.seed!(params[:random_seed])

            if typeof(params[param_to_sweep])==Int
                param_counter=convert(Int,floor(counter))
             else
                param_counter=counter
            end

            params[param_to_sweep]=param_counter

            # Actually run the simulation!
            results=specnet()
            agent_results = results[:agent_results]
            overall_results = results[:overall_results]

            push!(comp_df,
                (overall_results[:size_largest_comp],
                 overall_results[:num_comps],
                 (i*num_values+j)))

            insertcols!(agent_results, 1, param_to_sweep => repeat(counter:counter,nrow(agent_results)))
            insertcols!(agent_results, 2, :seed => repeat(params[:random_seed]:params[:random_seed],nrow(agent_results)))
            insertcols!(agent_results, 3, :sim_tag => repeat((i*num_values+j):(i*num_values+j),nrow(agent_results)))

            agent_line_df=[agent_line_df;agent_results]

            #increment the random seed to vary the results of simulations with the same params
            params[:random_seed]+=1
        end
        params[:random_seed]=original_seed
        counter += (end_value-start_value)/num_values
    end
    sim_tag=0
    rm("$(tempdir())/$(graph_name)_agent_results.csv", force=true)
    rm("$(tempdir())/$(graph_name)_simulation_results.csv", force=true)
    rm("$(tempdir())/$(graph_name)GiniSweepPlot.png", force=true)
    rm("$(tempdir())/$(graph_name)_wealth_heatmap.png", force=true)
    rm("$(tempdir())/$(graph_name)Component_GiniSweepPlots.png", force=true)
    rm("$(tempdir())/$(graph_name)ComponentSweepPlot.png", force=true)


    #this file contains all info with one line per agent in a given run of siml.jl
    CSV.write("$(tempdir())/$(graph_name)_agent_results.csv",agent_line_df)


    #once the main dataframe is made, plots may be drawn with data from agent_line_df
    #plotting and construction of dataframe is in separate for loop for optimal organiztion

    #dataframe containing only values to be plotted

    global plot_df=DataFrame(
        replace_this=Float64[],
        gini=Float64[],
        gini_lowCI=Float64[],
        gini_highCI=Float64[],
        number_components=Float64[],
        number_components_lowCI=Float64[],
        number_components_highCI=Float64[],
        size_largest_component=Float64[],
        size_largest_component_lowCI=Float64[],
        size_largest_component_highCI=Float64[],
    )
    names!(plot_df, prepend!(names(plot_df)[2:end], [param_to_sweep]))


    #change values of this dataframe to create other plots

    #loop to populate dataframe
    counter=start_value
    #weak solution to acccessing seed value, should be reworked
    mark_seed_value=original_seed

    for j = 1:num_values

        curr_param_value_ginis = []
        global curr_param_value_sizes = []
        global curr_param_value_components=[]
        for i=1:trials_per_value

            sim_tag=(j*num_values+i)

            # When not given a conf.level parameter, the R function Gini() from
            #   DescTools returns a single value: the Gini coefficient. We're
            #   not using DescTools bootstrapping to estimate the CI here,
            #   because we want a CI from the set of Gini coefficients from our
            #   trials_per_value runs, not the CI of a single run.
            current_sim_gini = convert(Float64, Gini(
                agent_line_df[agent_line_df.sim_tag.==sim_tag,:sugar]))
            push!(curr_param_value_ginis, current_sim_gini)
            push!(curr_param_value_sizes,comp_df[comp_df[:sim_tag].==sim_tag,:size_largest_comp])
            push!(curr_param_value_components,comp_df[comp_df[:sim_tag].==sim_tag,:num_comps])
            #adding results to the df
            push!(trial_line_df,(counter,mark_seed_value,current_sim_gini,
                sim_tag))

            #averaging and pushing data for wealth histogram
            mark_seed_value+=1
        end
        mark_seed_value=original_seed

        # Compute the average Gini, with CI, for this set of param values.
        bs = bootstrap(mean, curr_param_value_ginis,
            BasicSampling(params[:num_boot_samples]))
        ciGinis = confint(bs, BasicConfInt(.95))[1]

        #Compute the average size of the largest component, with a CI for current params

        bs1 = bootstrap(mean, curr_param_value_sizes, BasicSampling(params[:num_boot_samples]))
        ciSizes = confint(bs1, BasicConfInt(.95))[1]

        #Compute the average size of the largest component, with a CI for current params

        bs2 = bootstrap(mean, curr_param_value_components, BasicSampling(params[:num_boot_samples]))
        ciNumbers = confint(bs2, BasicConfInt(.95))[1]

        push!(plot_df, (counter,ciGinis[1], ciGinis[2], ciGinis[3],
                                ciNumbers[1], ciNumbers[2], ciNumbers[3],
                                ciSizes[1], ciSizes[2], ciSizes[3]))
        counter+=((end_value-start_value)/num_values)


    end
    #this file contains (currently) only the resulting Gini index from each simulation
    #trial_line_df=hcat(trial_line_df,comp_df)
    trial_line_df=join(trial_line_df, comp_df, on = :sim_tag)

    CSV.write("$(tempdir())/$(graph_name)_simulation_results.csv",trial_line_df)

    #drawing plots
    println("Creating $(param_to_sweep) Gini plot...")
    plotLG=plot(plot_df,
        layer(
            x=param_to_sweep, y=:gini_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:gini,
            Geom.line,
            Theme(default_color=colorant"navy", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:gini_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, ymin=:gini_lowCI, ymax=:gini_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"lightblue")
        ),
        Theme(background_color=colorant"white"),
        Guide.xlabel(string(param_to_sweep)), Guide.ylabel("Gini Index"))

    plotComponents=plot(plot_df,
        layer(
            x=param_to_sweep, y=:number_components_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgreen")
        ),
        layer(
            x=param_to_sweep, y=:number_components,
            Geom.line,
            Theme(default_color=colorant"green", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:number_components_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgreen")
        ),
        layer(
            x=param_to_sweep, ymin=:number_components_lowCI, ymax=:number_components_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"yellow")
        ),
        layer(
            x=param_to_sweep, y=:size_largest_component_lowCI,
            Geom.line,
            Theme(default_color=colorant"orange")
        ),
        layer(
            x=param_to_sweep, y=:size_largest_component,
            Geom.line,
            Theme(default_color=colorant"red", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:size_largest_component_highCI,
            Geom.line,
            Theme(default_color=colorant"orange")
        ),
        layer(
            x=param_to_sweep, ymin=:size_largest_component_lowCI, ymax=:size_largest_component_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"pink",key_position=:top)
        ),
        Guide.xlabel(string(param_to_sweep)),
        Guide.ylabel("Components", orientation=:vertical),
        Guide.xticks(ticks=:auto, label=true, orientation=:horizontal),
        Guide.manual_color_key("Legend",
            ["Agents in Largest Component      .", "Number of Components"],
            ["red", "green"]),
        style(background_color=colorant"white",key_position=:bottom))

    compGiniPlot=vstack(plotLG,plotComponents)

    draw(PNG("$(tempdir())/$(graph_name)GiniSweepPlot.png"), plotLG)
    draw(PNG("$(tempdir())/$(graph_name)Component_GiniSweepPlots.png",
        4inch, 6inch), compGiniPlot)
    draw(PNG("$(tempdir())/$(graph_name)ComponentSweepPlot.png"), plotComponents)
    println("Creating $(param_to_sweep) agent heatmap...")
    wealth_heatmap=plot(x=agent_line_df.sugar,y=agent_line_df[param_to_sweep],
        Geom.histogram2d,
        Theme(background_color=colorant"white"),
        Guide.ylabel(string(param_to_sweep)),
        Guide.xlabel("Agent Wealth"))
    draw(PNG("$(tempdir())/$(graph_name)_wealth_heatmap.png"), wealth_heatmap)

    return Dict(:agent_line_df => agent_line_df,
        :trial_line_df => trial_line_df,
        :plot_df => plot_df)
end

#runs a sweep for a given parameter once for each graph type,
#saving the dataframes and plots to multiple files

if graph_sweep
    sweep_results = Dict()
    graph_types=["erdos_renyi","scale_free","small_world","complete","empty"]
    println("sweeping for graph type")
    for graph_type in graph_types
        params[:whichGraph]=graph_type
        sweep_results[graph_type] = param_sweeper(graph_type)
    end
else
    sweep_results = param_sweeper(params[:whichGraph])
end
