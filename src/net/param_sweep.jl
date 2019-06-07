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
graph_sweep=false             #run the sweep three times, once for each graph type
num_steps=3
iter_per_step=4
original_seed=params[:random_seed]


function param_sweeper(graph_name)
    println("Starting sweep..")
    print("Sweeping for: $(param_to_sweep)")
    counter=start_value
    agent_line_df=DataFrame(agent=String[],sugar=Float64[], proto_id=Int[], counter_value=Float64[],iter_num_sweep=Int64[])
	iter_line_df=DataFrame(gini=Float64[], temp_random_seed=Int64[])
    params[:make_anims] = false  # We would never want this true for a sweep

    for i= 1:num_steps
        for j = 1:iter_per_step
        #setting the random seed and adding it to the DataFrame of the final gini of each simulation
	    Random.seed!(params[:random_seed])
		
		  if typeof(params[param_to_sweep])==Int64

                param_counter=convert(Int64,floor(counter))


            else
                param_counter=counter
            end

            params[param_to_sweep]=param_counter
            results=specnet()
            insert!(results,4,repeat(counter:counter,nrow(results)),:counter_value)
            insert!(results,5,repeat((i*iter_per_step + j):(i*iter_per_step + j),nrow(results)),:iter_num_sweep)
            agent_line_df=[agent_line_df;results]
		#increment the random seed to vary the results of simulations with the same params
		params[:random_seed]+=1
        end
        params[:random_seed]=original_seed
		counter+=(( end_value-  start_value)/num_steps)
    end
    rm("$(tempdir())/$(graph_name)_agent_results.csv", force=true)
	rm("$(tempdir())/$(graph_name)_simulation_results.csv", force=true)
    rm("$(tempdir())/$(graph_name)ParameterSweepPlot.png", force=true)
    #this file contains all info with one line per agent in a given run of siml.jl
    CSV.write("$(tempdir())/$(graph_name)_agent_results.csv",agent_line_df)
    

    #once the main dataframe is made, plots may be drawn with data from agent_line_df
    #plotting and construction of dataframe is in separate for loop for optimal organiztion

    #dataframe containing only values to be plotted
    plot_df=DataFrame(gini=Float64[],counter_value=Float64[])

    #change values of this dataframe to create other plots

    #loop to populate dataframe
    counter=start_value
    #weak solution to acccessing seed value, should be reworked
	mark_seed_value=original_seed
	for j = 1:num_steps
        total_gini=0
        for i=1:iter_per_step
            current_sim_gini=convert(Float64,
                ineq((agent_line_df[agent_line_df.iter_num_sweep.==(j*iter_per_step+i),
                    :sugar]), "Gini"))
			#adding to total gini to be averaged for plot (plot_df)
			total_gini+=current_sim_gini
			#adding the gini sim results to the df
			push!(iter_line_df,(current_sim_gini,mark_seed_value))
			mark_seed_value+=1
        end
		mark_seed_value=original_seed
        #change this code to populate above dataframe with alt data
        push!(plot_df, (total_gini/iter_per_step,counter))
        counter+=((end_value-start_value)/num_steps)
   end
    #this file contains (currently) only the resulting Gini index from each simulation
	#one line per run of sim.jl
	CSV.write("$(tempdir())/$(graph_name)_simulation_results.csv",iter_line_df)
    #drawing plot
    println("Creating $(param_to_sweep) plot...")
    plotLG=plot(x=plot_df.counter_value,y=plot_df.gini, Geom.point, Geom.line,
        Guide.xlabel(string(param_to_sweep)), Guide.ylabel("Gini Index"))
    draw(PNG("$(tempdir())/$(graph_name)ParameterSweepPlot.png"), plotLG)
end

#runs a sweep for a given parameter three times: one for each graph type
#saving the dataframes and plots to six different files

if graph_sweep
    graphs=["erdos_renyi","scale_free", "small_world"]
println("sweeping for graph type")
    for k in 1:3
        params[:whichGraph]=graphs[k]
        param_sweeper(string(graphs[k]))

    end
else
    param_sweeper(params[:whichGraph])
end
