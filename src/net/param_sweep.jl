#!/usr/bin/env julia
using Revise
using RCall
@rlibrary ineq
using DataFrames
using CSV
include("sim.jl")
include("setup_params.jl")


param_to_sweep=:λ             #parameter to iterate over *any parameter*
                              #for graph sweep param_to_sweep should not be exclusive to one graph type e.g. SF_prob
start_value=1                 #value to begin sweep
end_value=4                   #value to end sweep
num_values=3                  #number of distinct values to run
trials_per_value=4            #for each distinct value, number of independent sims to run
graph_sweep=false             #run the sweep once for each graph type
original_seed=params[:random_seed]

components=[[],[]]
global comp_df=convert(DataFrame,components')
rename!(comp_df,Dict(:x1 => :Agents_per_component, :x2 => :Number_of_components))

function param_sweeper(graph_name)
    println("Starting sweep..")
    print("Sweeping for: $(param_to_sweep)")
    counter=start_value

    agent_line_df = DataFrame(
        replace_this=Float64[],
        seed=Int64[],
        agent=String[],
        sugar=Float64[],
        proto_id=Int[])
    names!(agent_line_df,[param_to_sweep,:seed,:agent,:sugar,:proto_id])

    trial_line_df=DataFrame(
        replace_this=Float64[],
        seed=Int64[],
        gini=Float64[])
    names!(trial_line_df,[param_to_sweep,:seed,:gini])

    params[:make_anims] = false  # We would never want this true for a sweep

    for i = 1:num_values
        for j = 1:trials_per_value
            #setting the random seed and adding it to the DataFrame of the final gini of each simulation
            Random.seed!(params[:random_seed])
            
            if typeof(params[param_to_sweep])==Int64
                param_counter=convert(Int64,floor(counter))
            else
                param_counter=counter
            end

            params[param_to_sweep]=param_counter

            # Actually run the simulation!
            results=specnet()

			#finding the breakdown of graph components and pushing that to component data frame
			component_vertices=connected_components(graph)
	        single_sim_components=[]
            for comp in component_vertices
	            push!(single_sim_components,length(comp))
	        end
			push!(comp_df,(single_sim_components',[length(component_vertices)]'))
            
			
          

            insertcols!(results, 1, param_to_sweep => repeat(counter:counter,nrow(results)))
            insertcols!(results, 2, :seed => repeat(params[:random_seed]:params[:random_seed],nrow(results)))

            agent_line_df=[agent_line_df;results]

            #increment the random seed to vary the results of simulations with the same params
            params[:random_seed]+=1
        end
        params[:random_seed]=original_seed
        counter += (end_value-start_value)/num_values
    end
    rm("$(tempdir())/$(graph_name)_agent_results.csv", force=true)
    rm("$(tempdir())/$(graph_name)_simulation_results.csv", force=true)
    rm("$(tempdir())/$(graph_name)ParameterSweepPlot.png", force=true)
    #this file contains all info with one line per agent in a given run of siml.jl
    CSV.write("$(tempdir())/$(graph_name)_agent_results.csv",agent_line_df)
    

    #once the main dataframe is made, plots may be drawn with data from agent_line_df
    #plotting and construction of dataframe is in separate for loop for optimal organiztion

    #dataframe containing only values to be plotted

    plot_df=DataFrame(gini=Float64[],param_to_sweep=Float64[])


    #change values of this dataframe to create other plots
	
    #loop to populate dataframe
    counter=start_value
    #weak solution to acccessing seed value, should be reworked
    mark_seed_value=original_seed
    for j = 1:num_values
        total_gini=0
        for i=1:trials_per_value
            current_sim_gini=convert(Float64,

                ineq((agent_line_df[agent_line_df.seed.==mark_seed_value,:sugar]), "Gini"))
            #adding to total gini to be averaged for plot (plot_df)
            total_gini+=current_sim_gini
            #adding the gini sim results to the df
            push!(trial_line_df,(counter,mark_seed_value,current_sim_gini))
            mark_seed_value+=1

        end
        mark_seed_value=original_seed
        #change this code to populate above dataframe with alt data
        push!(plot_df, (total_gini/trials_per_value,counter))
        counter+=((end_value-start_value)/num_values)
    end
    #this file contains (currently) only the resulting Gini index from each simulation

	#one line per run of sim.jl
    global matching_comp_df = comp_df[setdiff(1:end, 1), :]
    iter_line_df=hcat(trial_line_df,matching_comp_df)
    CSV.write("$(tempdir())/$(graph_name)_simulation_results.csv",trial_line_df)

    #drawing plot
    println("Creating $(param_to_sweep) plot...")
    plotLG=plot(x=plot_df.param_to_sweep,y=plot_df.gini, Geom.point, Geom.line,
        Guide.xlabel(string(param_to_sweep)), Guide.ylabel("Gini Index"))
    draw(PNG("$(tempdir())/$(graph_name)ParameterSweepPlot.png"), plotLG)
end

#runs a sweep for a given parameter once for each graph type,
#saving the dataframes and plots to multiple files

if graph_sweep
    graphs=["erdos_renyi","scale_free","small_world","complete","empty"]
    println("sweeping for graph type")
    for graph in graphs
        params[:whichGraph]=graph
        param_sweeper(string(graph))
    end
else
    param_sweeper(params[:whichGraph])
end
