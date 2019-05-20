#!/usr/bin/env julia
using Revise
using DataFrames
using CSV
include("sim.jl")
include("setup_params.jl")


param_to_sweep=:salary_range    #parameter to iterate over
start_value=10                  #value to begin sweep
end_value=210                   #value to end sweep
graph_sweep=false               #run sweep for each graph type
num_steps=100



function param_sweeper(graph_name)
    counter=start_value
    large_df=DataFrame(agent=String[],sugar=Float64[], proto_id=Int[], counter_value=Float64[])

    params[:make_anims] = false  # We would never want this true for a sweep 

    for i= 1:num_steps

        if typeof(params[param_to_sweep])==Int64

            param_counter=convert(Int64,floor(counter))


        else
            param_counter=counter
        end

        params[param_to_sweep]=param_counter
        results=specnet()
        insert!(results,4,repeat(counter:counter,nrow(results)),:counter_value)
        large_df=[large_df;results]
        counter+=(( end_value-  start_value)/num_steps)

    end
    rm("$(tempdir())/$(graph_name)ParameterSweepDataFrame.csv", force=true)
    rm("$(tempdir())/$(graph_name)ParameterSweepPlot.png", force=true)
    CSV.write("$(tempdir())/$(graph_name)ParameterSweep.csv",large_df)

    #once the main dataframe is made, plots may be drawn with data from large_df
    #plotting and construction of dataframe is in separate for loop for optimal organiztion

    #dataframe containing only values to be plotted
    plot_df=DataFrame(gini=Float64[],counter_value=Float64[])  #change values of this dataframe
                                                               #to create other plots
    #loop to populate dataframe

    counter=start_value
    for j = 1:num_steps
        counter+=((end_value-start_value)/num_steps)

        #change this code to populate above dataframe with alt data
        push!(plot_df, (ineq((large_df[large_df.counter_value.==counter, :sugar]), "Gini"),counter))
    end

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

    for k in 1:3
        params[:whichGraph]=graphs[k]
        param_sweeper(string(graphs[k]))

    end
else
    param_sweeper(params[:whichGraph])
end

