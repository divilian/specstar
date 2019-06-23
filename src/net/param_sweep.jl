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
global social_connectivity_df=DataFrame(size_largest_comp=Int[],num_comps=Int[],average_proto_size=Float64[], num_protos=Int[],sim_tag=Int[])


function param_sweeper(graph_name; additional_params...)

    @assert all_parameters_legit(additional_params)
    merge!(params, Dict(additional_params))

    println("Starting sweep..")
    print("Sweeping for: $(param_to_sweep)")
    counter=start_value

    global agent_line_df = DataFrame(
        replace_this=Float64[],
        seed=Int[],
        sim_tag=Int[],
        agent=String[],
        sugar=Float64[],
        proto_id=Int[],
    )
    names!(agent_line_df,
        prepend!(names(agent_line_df)[2:end], [param_to_sweep]))

    global trial_line_df=DataFrame(
        replace_this=Float64[],
        seed=Int[],
        sim_tag=Int[],
        gini=Float64[],
    )
    names!(trial_line_df,
        prepend!(names(trial_line_df)[2:end], [param_to_sweep]))

    global iter_line_df=DataFrame(
        replace_this=Float64[],
        seed=Int[],
        sim_tag=Int[],
        iter=Int[],
        stage=Int[],
    )
    names!(iter_line_df,
        prepend!(names(iter_line_df)[2:end], [param_to_sweep]))

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
            iter_results = results[:iter_results]
            overall_results = results[:overall_results]

            push!(social_connectivity_df,
                (overall_results[:size_largest_comp],
                 overall_results[:num_comps],
				 overall_results[:average_proto_size],
	             overall_results[:num_protos],

                 (i*num_values+j)))

            add_sim_info!(agent_results, param_counter, counter, i*num_values+j)
            iter_results = iter_results[[:iter,:stage]]  # only need these for now
            add_sim_info!(iter_results, param_counter, counter, i*num_values+j)

            agent_line_df=[agent_line_df;agent_results]
            iter_line_df=[iter_line_df;iter_results]

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
    rm("$(tempdir())/$(graph_name)ProtoPropertiesSweep.png", force=true)

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
        time_to_stage2=Float64[],
        time_to_stage2_lowCI=Float64[],
        time_to_stage2_highCI=Float64[],
        time_to_stage3=Float64[],
        time_to_stage3_lowCI=Float64[],
        time_to_stage3_highCI=Float64[],
		average_proto_size=Float64[],
        average_proto_size_lowCI=Float64[],
        average_proto_size_highCI=Float64[],
		num_protos=Float64[],
        num_protos_lowCI=Float64[],
        num_protos_highCI=Float64[]
    )
    names!(plot_df, prepend!(names(plot_df)[2:end], [param_to_sweep]))


    #change values of this dataframe to create other plots

    #loop to populate dataframe
    counter=start_value
    #weak solution to acccessing seed value, should be reworked
    mark_seed_value=original_seed

    for j = 1:num_values

        curr_param_value_ginis = []
        curr_param_value_sizes = []
        curr_param_value_components = []
		curr_param_value_proto_size=[]
		curr_param_value_num_protos=[]
        curr_param_value_p2s=[]
        curr_param_value_p3s=[]
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
            push!(curr_param_value_sizes,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:size_largest_comp])
            push!(curr_param_value_components,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_comps])
            push!(curr_param_value_proto_size,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:average_proto_size])
			push!(curr_param_value_num_protos,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_protos])
			push!(curr_param_value_p2s,first_iter_of_stage(iter_line_df, 2, sim_tag))
            push!(curr_param_value_p3s,first_iter_of_stage(iter_line_df, 3, sim_tag))
            #adding results to the df
            push!(trial_line_df,(counter,mark_seed_value,sim_tag,
                current_sim_gini))

            #averaging and pushing data for wealth histogram
            mark_seed_value+=1
        end
        mark_seed_value=original_seed

        # Compute the average Gini, with CI, for this set of param values.
        bs = bootstrap(mean, curr_param_value_ginis,
            BasicSampling(params[:num_boot_samples]))
        ciGinis = confint(bs, BasicConfInt(.95))[1]

        #Compute the average size of the largest component, with a CI for current params
        bs = bootstrap(mean, curr_param_value_sizes, BasicSampling(params[:num_boot_samples]))
        ciSizes = confint(bs, BasicConfInt(.95))[1]
        #CI for average proto size for current params
        bs = bootstrap(mean, curr_param_value_proto_size, BasicSampling(params[:num_boot_samples]))
        ciProtoSizes = confint(bs, BasicConfInt(.95))[1]
        #CI for number of protos for current params
        bs = bootstrap(mean, curr_param_value_num_protos, BasicSampling(params[:num_boot_samples]))
        ciNumProtos = confint(bs, BasicConfInt(.95))[1]
        #Compute the average number of components, with a CI for current params
        bs = bootstrap(mean, curr_param_value_components, BasicSampling(params[:num_boot_samples]))
        ciNumbers = confint(bs, BasicConfInt(.95))[1]

        #Compute the average time to stage 2, with a CI for current params
        bs = bootstrap(mean, curr_param_value_p2s, BasicSampling(params[:num_boot_samples]))
        p2Times = confint(bs, BasicConfInt(.95))[1]

        #Compute the average time to stage 3, with a CI for current params
        bs = bootstrap(mean, curr_param_value_p3s, BasicSampling(params[:num_boot_samples]))
        p3Times = confint(bs, BasicConfInt(.95))[1]

        push!(plot_df, (counter,ciGinis[1], ciGinis[2], ciGinis[3],
                                ciNumbers[1], ciNumbers[2], ciNumbers[3],
                                ciSizes[1], ciSizes[2], ciSizes[3],
                                p2Times[1], p2Times[2], p2Times[3],
                                p3Times[1], p3Times[2], p3Times[3],
								ciProtoSizes[1], ciProtoSizes[2], ciProtoSizes[3],
								ciNumProtos[1], ciNumProtos[2], ciNumProtos[3]))
        counter+=((end_value-start_value)/num_values)


    end
    #this file contains (currently) only the resulting Gini index from each simulation
    #trial_line_df=hcat(trial_line_df,social_connectivity_df)
    trial_line_df=join(trial_line_df, social_connectivity_df, on = :sim_tag)

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
		

    println("Creating $(param_to_sweep) components plot...")
    plotComponents=plot(plot_df,
        layer(
            x=param_to_sweep, y=:number_components_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:number_components,
            Geom.line,
            Theme(default_color=colorant"green", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:number_components_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, ymin=:number_components_lowCI, ymax=:number_components_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"lightgreen")
        ),
        layer(
            x=param_to_sweep, y=:size_largest_component_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:size_largest_component,
            Geom.line,
            Theme(default_color=colorant"red", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:size_largest_component_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
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

    println("Creating $(param_to_sweep) time-to-stop plot...")
    to_stage_df = plot_df[plot_df[:time_to_stage3] .> 0, :]
    plotTTS=plot(to_stage_df,
        layer(
            x=param_to_sweep, y=:time_to_stage2_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:time_to_stage2,
            Geom.line,
            Theme(default_color=colorant"blue", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:time_to_stage2_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, ymin=:time_to_stage2_lowCI, ymax=:time_to_stage2_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"lightblue",key_position=:top)
        ),
        layer(
            x=param_to_sweep, y=:time_to_stage3_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:time_to_stage3,
            Geom.line,
            Theme(default_color=colorant"red", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:time_to_stage3_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, ymin=:time_to_stage3_lowCI, ymax=:time_to_stage3_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"pink",key_position=:top)
        ),
        Guide.xlabel(string(param_to_sweep)),
        Guide.ylabel("Time to reach stage", orientation=:vertical),
        Guide.xticks(ticks=:auto, label=true, orientation=:horizontal),
        Guide.manual_color_key("Legend",
            ["Stage 2      .", "Stage 3"],
            ["blue", "red"]),
        style(background_color=colorant"white",key_position=:bottom))
    plotProtos=plot(plot_df,
        layer(
            x=param_to_sweep, y=:num_protos_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:num_protos,
            Geom.line,
            Theme(default_color=colorant"darkgreen", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:num_protos_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, ymin=:num_protos_lowCI, ymax=:num_protos_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"lightgreen",key_position=:top)
        ),
        layer(
            x=param_to_sweep, y=:average_proto_size_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, y=:average_proto_size,
            Geom.line,
            Theme(default_color=colorant"brown", line_width=.5mm)
        ),
        layer(
            x=param_to_sweep, y=:average_proto_size_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=param_to_sweep, ymin=:average_proto_size_lowCI, ymax=:average_proto_size_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"orange",key_position=:top)
        ),
        Guide.xlabel(string(param_to_sweep)),
        Guide.ylabel("Protos", orientation=:vertical),
        Guide.xticks(ticks=:auto, label=true, orientation=:horizontal),
        Guide.manual_color_key("Legend",
            [ "Average Size of Proto      .","Number of Protos"],
            ["brown", "green"]),
        style(background_color=colorant"white",key_position=:bottom))
		
		
    tallPlot=vstack(plotLG,plotComponents,plotProtos,plotTTS)

    draw(PNG("$(tempdir())/$(graph_name)GiniSweepPlot.png"), plotLG)
    draw(PNG("$(tempdir())/$(graph_name)ComponentSweepPlot.png"), plotComponents)
    draw(PNG("$(tempdir())/$(graph_name)TimeToStages.png"), plotTTS)
	draw(PNG("$(tempdir())/$(graph_name)ProtoPropertiesSweep.png"), plotProtos)

    draw(PNG("$(tempdir())/$(graph_name)TallPlot.png",
        5inch, 9inch), tallPlot)
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

function add_sim_info!(df, param_counter, counter, sim_tag)
    insertcols!(df, 1, param_to_sweep => repeat(counter:counter,nrow(df)))
    insertcols!(df, 2, :seed => repeat(params[:random_seed]:params[:random_seed],nrow(df)))
    insertcols!(df, 3, :sim_tag => repeat(sim_tag:sim_tag,nrow(df)))
end

function first_iter_of_stage(df, stage_num, sim_tag)
    iters = df[(df[:sim_tag].==sim_tag).&(df[:stage].==stage_num), :iter]
    length(iters) > 0 ? minimum(iters) : 0
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
