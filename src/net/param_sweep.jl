#!/usr/bin/env julia
using Revise
using RCall
@rlibrary DescTools
using DataFrames
using CSV
using Bootstrap
include("sim.jl")
include("setup_params.jl")


param_to_sweep=:λ   #parameter to iterate over *any parameter*
                     #for graph sweep param_to_sweep should not be exclusive to one graph type e.g. SF_prob
start_value=0.06     #value to begin sweep
end_value=6         #value to end sweep
num_values=100        #number of distinct values to run
trials_per_value=14   #for each distinct value, number of independent sims to run
graph_sweep=false    #run the sweep once for each graph type
original_seed=params[:random_seed]
repeat_param=:N
repeat_start_value=10
repeat_end_value=15
run_repeat_sweep=false

components=[[],[]]
global social_connectivity_df=DataFrame(
    size_largest_comp=Int[],
    num_comps=Int[],
    average_proto_size_pre=Float64[],
    average_proto_size_post=Float64[],
    num_agents_in_proto=Int[],
    num_living_agents_pre=Float64[],
    num_living_agents_post=Float64[],
    num_protos_pre=Float64[],
    num_protos_post=Float64[],
    avg_lifespan_isos=Float64[],
    avg_lifespan_nonisos=Float64[],
    sim_tag=Int[])


function param_sweeper(; additional_params...)

    @assert all_parameters_legit(additional_params)
    merge!(params, Dict(additional_params))

    prc("Starting sweep...\n")
    prc("Sweeping for: $(param_to_sweep)\n")
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
            starvation_results = results[:starvation_results]
            avg_lifespans = results[:avg_lifespans]

            push!(social_connectivity_df,
                (overall_results[:size_largest_comp],
                 overall_results[:num_comps],
                 overall_results[:average_proto_size],
                 starvation_results[:average_proto_size],
                 overall_results[:num_agents_in_proto],
                 overall_results[:num_living_agents],     # shouldn't be here
                 starvation_results[:num_living_agents],  # shouldn't be here
                 overall_results[:num_protos],     # shouldn't be here
                 starvation_results[:num_protos],  # shouldn't be here
                 avg_lifespans[1],     # shouldn't be here
                 avg_lifespans[2],     # shouldn't be here

                 (i*trials_per_value+j)))

            add_sim_info!(agent_results, param_counter, counter, i*trials_per_value+j)
            iter_results = iter_results[[:iter,:stage]]  # only need these for now
            add_sim_info!(iter_results, param_counter, counter, i*trials_per_value+j)

            agent_line_df=[agent_line_df;agent_results]
            iter_line_df=[iter_line_df;iter_results]

            #increment the random seed to vary the results of simulations with the same params
            params[:random_seed]+=1
        end
        params[:random_seed]=original_seed
        counter += (end_value-start_value)/num_values
    end
    sim_tag=0

    rm("$(tempdir())/agent_results.csv", force=true)
    rm("$(tempdir())/simulation_results.csv", force=true)
    rm("$(tempdir())/GiniSweepPlot.png", force=true)
    rm("$(tempdir())/wealth_heatmap.png", force=true)
    rm("$(tempdir())/Component_GiniSweepPlots.png", force=true)
    rm("$(tempdir())/ComponentSweepPlot.png", force=true)
    rm("$(tempdir())/ProtoPropertiesSweep.png", force=true)
    rm("$(tempdir())/repeat_sweep_plot.png", force=true)


    #this file contains all info with one line per agent in a given run of siml.jl
    CSV.write("$(tempdir())/agent_results.csv",agent_line_df)


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
        average_proto_size_pre=Float64[],
        average_proto_size_pre_lowCI=Float64[],
        average_proto_size_pre_highCI=Float64[],
        average_proto_size_post=Float64[],
        average_proto_size_post_lowCI=Float64[],
        average_proto_size_post_highCI=Float64[],
        num_agents_in_proto=Float64[],
        num_agents_in_proto_lowCI=Float64[],
        num_agents_in_proto_highCI=Float64[],
        num_living_agents_pre=Float64[],
        num_living_agents_pre_lowCI=Float64[],
        num_living_agents_pre_highCI=Float64[],
        num_living_agents_post=Float64[],
        num_living_agents_post_lowCI=Float64[],
        num_living_agents_post_highCI=Float64[],
        num_protos_pre=Float64[],
        num_protos_pre_lowCI=Float64[],
        num_protos_pre_highCI=Float64[],
        num_protos_post=Float64[],
        num_protos_post_lowCI=Float64[],
        num_protos_post_highCI=Float64[],
        avg_lifespan_isos=Float64[],
        avg_lifespan_isos_lowCI=Float64[],
        avg_lifespan_isos_highCI=Float64[],
        avg_lifespan_nonisos=Float64[],
        avg_lifespan_nonisos_lowCI=Float64[],
        avg_lifespan_nonisos_highCI=Float64[],
    )
    names!(plot_df, prepend!(names(plot_df)[2:end], [param_to_sweep]))


    #change values of this dataframe to create other plots

    #loop to populate dataframe
    counter=start_value
    #weak solution to acccessing seed value, should be reworked
    mark_seed_value=original_seed

    for j = 1:num_values

        # "curr" means "for all trials of the CURRent parameter value."
        curr_ginis = []
        curr_sizes = []
        curr_components = []
        curr_proto_size_pre=[]
        curr_proto_size_post=[]

        curr_num_agents_in_proto=[]
        curr_s2s=[]
        curr_s3s=[]
        curr_num_living_agents_pre=[]
        curr_num_living_agents_post=[]
        curr_num_protos_pre=[]
        curr_num_protos_post=[]
        curr_avg_lifespan_isos=[]
        curr_avg_lifespan_nonisos=[]
        for i=1:trials_per_value

            sim_tag=(j*trials_per_value+i)

            # When not given a conf.level parameter, the R function Gini() from
            #   DescTools returns a single value: the Gini coefficient. We're
            #   not using DescTools bootstrapping to estimate the CI here,
            #   because we want a CI from the set of Gini coefficients from our
            #   trials_per_value runs, not the CI of a single run.
            current_sim_gini = convert(Float64, Gini(
                agent_line_df[agent_line_df.sim_tag.==sim_tag,:sugar]))
            push!(curr_ginis, current_sim_gini)
            push!(curr_sizes,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:size_largest_comp][1])
            push!(curr_components,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_comps][1])
            push!(curr_proto_size_pre,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:average_proto_size_pre][1])
            push!(curr_proto_size_post,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:average_proto_size_post][1])
            push!(curr_num_agents_in_proto,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_agents_in_proto][1])
            push!(curr_s2s,first_iter_of_stage(iter_line_df, 2, sim_tag))
            push!(curr_s3s,first_iter_of_stage(iter_line_df, 3, sim_tag))
            push!(curr_num_living_agents_pre,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_living_agents_pre][1])
            push!(curr_num_living_agents_post,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_living_agents_post][1])
            push!(curr_num_protos_pre,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_protos_pre][1])
            push!(curr_num_protos_post,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:num_protos_post][1])
            push!(curr_avg_lifespan_isos,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:avg_lifespan_isos][1])
            push!(curr_avg_lifespan_nonisos,social_connectivity_df[social_connectivity_df[:sim_tag].==sim_tag,:avg_lifespan_nonisos][1])
            #adding results to the df
            push!(trial_line_df,(counter,mark_seed_value,sim_tag,
                current_sim_gini))

            mark_seed_value+=1
        end
        mark_seed_value=original_seed

        # Compute the average Gini, with CI, for this set of param values.
        if sum([!isnan(g) for g ∈  curr_ginis]) < 3
            ciGinis = [NaN, NaN, NaN]
        else
            curr_ginis = [ g for g ∈  curr_ginis if !isnan(g) ]
            bs = bootstrap(mean, curr_ginis,
                BasicSampling(params[:num_boot_samples]))
            ciGinis = confint(bs, BasicConfInt(.95))[1]
        end
        filter!(x->x≠0,curr_proto_size_pre)
        filter!(x->x≠0,curr_proto_size_post)
        if length(curr_proto_size_pre)==0
            push!(curr_proto_size_pre,0)
        end
        if length(curr_proto_size_post)==0
            push!(curr_proto_size_post,0)
        end
        #Compute the average size of the largest component, with a CI for current params
        bs = bootstrap(mean, curr_sizes, BasicSampling(params[:num_boot_samples]))
        ciSizes = confint(bs, BasicConfInt(.95))[1]
        #CI for average proto size for current params
        bs = bootstrap(mean, curr_proto_size_pre, BasicSampling(params[:num_boot_samples]))
        ciProtoSizesPre = confint(bs, BasicConfInt(.95))[1]
        bs = bootstrap(mean, curr_proto_size_post, BasicSampling(params[:num_boot_samples]))
        ciProtoSizesPost = confint(bs, BasicConfInt(.95))[1]
        #CI for number of agents in proto for current params
        bs = bootstrap(mean, curr_num_agents_in_proto, BasicSampling(params[:num_boot_samples]))
        ciNumAginP = confint(bs, BasicConfInt(.95))[1]
        #Compute the average number of components, with a CI for current params
        bs = bootstrap(mean, curr_components, BasicSampling(params[:num_boot_samples]))
        ciNumbers = confint(bs, BasicConfInt(.95))[1]

        #Compute the average time to stage 2, with a CI for current params
        bs = bootstrap(mean, curr_s2s, BasicSampling(params[:num_boot_samples]))
        ciS2Times = confint(bs, BasicConfInt(.95))[1]

        #Compute the average time to stage 3, with a CI for current params
        bs = bootstrap(mean, curr_s3s, BasicSampling(params[:num_boot_samples]))
        ciS3Times = confint(bs, BasicConfInt(.95))[1]

        #CI for average number of living agents, pre-starvation
        bs = bootstrap(mean, curr_num_living_agents_pre,
            BasicSampling(params[:num_boot_samples]))
        ciNumLivingPre = confint(bs, BasicConfInt(.95))[1]
        #CI for average number of living agents, post-starvation
        bs = bootstrap(mean, curr_num_living_agents_post,
            BasicSampling(params[:num_boot_samples]))
        ciNumLivingPost = confint(bs, BasicConfInt(.95))[1]

        #CI for average number of protos, pre-starvation
        bs = bootstrap(mean, curr_num_protos_pre, BasicSampling(params[:num_boot_samples]))
        ciNumProtosPre = confint(bs, BasicConfInt(.95))[1]

        #CI for average number of protos, pre-starvation
        bs = bootstrap(mean, curr_num_protos_post, BasicSampling(params[:num_boot_samples]))
        ciNumProtosPost = confint(bs, BasicConfInt(.95))[1]

        #CIs for average lifespan
        bs = bootstrap(mean, curr_avg_lifespan_isos, BasicSampling(params[:num_boot_samples]))
        ciAvgLifespanIsos = confint(bs, BasicConfInt(.95))[1]
        bs = bootstrap(mean, curr_avg_lifespan_nonisos, BasicSampling(params[:num_boot_samples]))
        ciAvgLifespanNonisos = confint(bs, BasicConfInt(.95))[1]

        push!(plot_df, (counter,ciGinis[1], ciGinis[2], ciGinis[3],
                                ciNumbers[1], ciNumbers[2], ciNumbers[3],
                                ciSizes[1], ciSizes[2], ciSizes[3],
                                ciS2Times[1], ciS2Times[2], ciS2Times[3],
                                ciS3Times[1], ciS3Times[2], ciS3Times[3],
                                ciProtoSizesPre[1], ciProtoSizesPre[2], ciProtoSizesPre[3],
                                ciProtoSizesPost[1], ciProtoSizesPost[2], ciProtoSizesPost[3],
                                ciNumAginP[1], ciNumAginP[2], ciNumAginP[3],
                                ciNumLivingPre[1], ciNumLivingPre[2], ciNumLivingPre[3],
                                ciNumLivingPost[1], ciNumLivingPost[2], ciNumLivingPost[3],
                                ciNumProtosPre[1], ciNumProtosPre[2], ciNumProtosPre[3],
                                ciNumProtosPost[1], ciNumProtosPost[2], ciNumProtosPost[3],
                                ciAvgLifespanIsos[1], ciAvgLifespanIsos[2], ciAvgLifespanIsos[3],
                                ciAvgLifespanNonisos[1], ciAvgLifespanNonisos[2], ciAvgLifespanNonisos[3],
        ))
        counter+=((end_value-start_value)/num_values)


    end
    #this file contains (currently) only the resulting Gini index from each simulation
    trial_line_df=join(trial_line_df, social_connectivity_df, on = :sim_tag)


    CSV.write("$(tempdir())/simulation_results.csv",trial_line_df)

    global repeat_sweep_layer=layer(x=plot_df[param_to_sweep],y=plot_df[:gini],Geom.line,
                                    Theme(line_width=1mm,
                                    default_color=RGB((params[repeat_param]-repeat_start_value)/(repeat_end_value-repeat_start_value),
                                    (params[repeat_param]-repeat_start_value)/(repeat_end_value-repeat_start_value),
                                    (params[repeat_param]-repeat_start_value)/2*(repeat_end_value-repeat_start_value)) ))

    ginip = draw_plot(plot_df, param_to_sweep, Dict("gini"=>"navy"),
        y_label="Gini")

    compp = draw_plot(plot_df, param_to_sweep,
        Dict("number_components"=>"green", "size_largest_component" => "red"),
        y_label="Components")

    protop = draw_plot(plot_df, param_to_sweep,
        Dict(
            "num_protos_pre"=>"blue",
            "num_protos_post"=>"orange",
            "average_proto_size_pre" => "green",
            "average_proto_size_post" => "red",
        ),
        y_label="Protos",
        )

    numagentspp = draw_plot(plot_df, param_to_sweep,
        Dict("num_agents_in_proto" => "green"),
            
        Dict("num_agents_in_proto"=> :solid
            ),
        y_label="Number agents in proto",
        extra=[Guide.annotation(compose(context(), Compose.text(minimum(plot_df[param_to_sweep]), params[:N], "N=$(params[:N])", hleft, vtop))),
         Coord.Cartesian(ymin=0, ymax=params[:N]),
         layer(yintercept=[params[:N]], Geom.hline(style=:dot, color=colorant"navy"))[1]])

    numlivingp = draw_plot(plot_df, param_to_sweep,
        Dict("num_living_agents_pre" => "green",
            "num_living_agents_post" => "red",
            "num_protos_pre" => "blue",
            "num_protos_post" => "orange"),
        Dict("num_living_agents_pre"=> :solid,
            "num_living_agents_post" => :solid,
            "num_protos_pre" => :dot,
            "num_protos_post" => :dot),
        y_label="Number of living agents/protos",
        extra=[Guide.annotation(compose(context(), Compose.text(minimum(plot_df[param_to_sweep]), params[:N], "N=$(params[:N])", hleft, vtop))),
         Coord.Cartesian(ymin=0, ymax=params[:N]),
         layer(yintercept=[params[:N]], Geom.hline(style=:dot, color=colorant"navy"))[1]])

    # zeros are meaningless in time-to-stage plots
    to_stage_df = plot_df[plot_df[:time_to_stage3] .> 0, :]
    ttsp = draw_plot(to_stage_df, param_to_sweep,
        Dict("time_to_stage2"=>"blue", "time_to_stage3" => "red"),
        y_label="Time to reach stage")

    tallPlot=vstack(ginip,compp,ttsp,protop)
    draw(PNG("$(tempdir())/tallPlot.png", 5inch, 9inch), tallPlot)

    lifespanp = draw_plot(plot_df, param_to_sweep,
        Dict("avg_lifespan_isos"=>"orange", "avg_lifespan_nonisos" => "navy"),
        extra=[
            Coord.Cartesian(ymin=0),
            Guide.title("λ=$(params[:λ])"),   # specifically for figure 4
            Guide.xlabel("White noise (σ²)"),
            Guide.manual_color_key("", ["Isolates    .", "Non-isolates"],
                ["orange","navy"])
        ],
        y_label="Mean lifespan (iterations)")

    wealth_heatmap=plot(x=agent_line_df.sugar,y=agent_line_df[param_to_sweep],
        Geom.histogram2d,
        Theme(background_color=colorant"white"),
        Guide.ylabel(string(param_to_sweep)),
        Guide.xlabel("Agent Wealth"))
    draw(PNG("$(tempdir())/wealth_heatmap.png"), wealth_heatmap)

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

function repeat_sweep()
    params[repeat_param]=repeat_start_value
    param_sweeper()
    plot_repeat=plot(repeat_sweep_layer,Guide.XLabel("$(param_to_sweep)"),
         Guide.YLabel("Gini"),
         Guide.Title("Repeat sweep across $(repeat_param)"))
    for i=(repeat_start_value+1):repeat_end_value
        params[repeat_param]=i
        param_sweeper()
        append!(plot_repeat.layers,repeat_sweep_layer)
    end

    draw(PNG("$(tempdir())/repeat_sweep_plot.png"), plot_repeat)
end

# draw_plot() -- create, write to file, and return a parameter sweep line plot.
# Parameters:
#  plot_df -- the master plotting data frame, created in param_sweeper().
#  param_to_sweep -- the name of the main independent variable.
#  vars_colors -- maps names of dependent variable(s) (which, if plot_CIs
#    is true, must also be present with _lowCI and _highCI suffixes) to
#    line colors.
#  y_label -- text to plot vertically on the left
#  extra -- an optional array of other layers or geoms to add to plot.
#  plot_CIs -- if true, include confidence interval bands (in lighter color)
function draw_plot(plot_df, param_to_sweep, vars_colors=Dict{String,String},
        vars_linestyles=Dict{String,Symbol}();
        y_label=nothing, extra=[], plot_CIs=true)

    prd("Drawing: $(vars_colors)\n")
    plot_df = consolidate(plot_df, param_to_sweep, vars_colors, plot_CIs)
    if nrow(plot_df) == 0
        prc("WARNING: Cannot produce plot for $(param_to_sweep).")
        return
    end
    layers = Layer[]
    for (var, color) in vars_colors
        append!(layers, layer(plot_df,
            x=param_to_sweep, y=var,
            Geom.line,
            Theme(default_color=color,
                line_style=[get(vars_linestyles, var, :solid)],
                line_width=.5mm)
        ))
        if plot_CIs
            append!(layers, layer(plot_df,
                x=param_to_sweep, y=var*"_lowCI",
                Geom.line,
                Theme(default_color=colorant"lightgray")
            ))
            append!(layers, layer(plot_df,
                x=param_to_sweep, y=var*"_highCI",
                Geom.line,
                Theme(default_color=colorant"lightgray")
            ))
            append!(layers, layer(plot_df,
                x=param_to_sweep, ymin=var*"_lowCI", ymax=var*"_highCI",
                Geom.ribbon,
                Theme(default_color=lighter_shade_of[color])
            ))
        end

    end

    p = plot(plot_df, layers,
        Theme(background_color=colorant"white"),
        Guide.xlabel(string(param_to_sweep)),
        Guide.ylabel(y_label, orientation=:vertical),
        Guide.xticks(ticks=:auto, label=true, orientation=:horizontal),
        style(background_color=colorant"white",key_position=:bottom)
    )
    if length(vars_colors) > 1
        vars = sort(collect(keys(vars_colors)))
        colors = [ vars_colors[x] for x in vars ]
        vars[1:end-1] = [ v*"   ." for v in vars[1:end-1] ]
        # TODO: add legend that combines both color and linestyle, where appropriate.
        push!(p, Guide.manual_color_key(nothing, vars, colors))
    end
    [ push!(p, e) for e in extra ]

    draw(PNG("$(tempdir())/$(join(keys(vars_colors),"_")).png"), p)
    return p
end
lighter_shade_of = Dict(
    "navy" => "lightblue",
    "blue" => "lightblue",
    "orange" => "pink",
    "red" => "pink",
    "green" => "lightgreen",
    "brown" => "orange",
    "black" => "grey",
    "grey" => "lightgrey",
)

# Keep only columns we will use for this plot, and remove rows with missing
#   data so it is not (misleadingly) plotted.
function consolidate(df, param_to_sweep, vars_colors, plot_CIs)
    the_cols = collect(keys(vars_colors))
    if plot_CIs
        extended_cols = vcat(the_cols,
            [string(q)*"_lowCI" for q in the_cols])
        extended_cols = vcat(extended_cols,
            [string(q)*"_highCI" for q in the_cols])
        the_cols = extended_cols
    end
    the_cols = vcat(the_cols, string(param_to_sweep))
    df = df[[q for q in names(df) if string(q) in the_cols]]
    for col in names(df)
        df[col] = map(x->isnan(x) ? missing : x, df[col])
    end
    df[completecases(df),:]
end

#runs a sweep for a given parameter once for each graph type,
#saving the dataframes and plots to multiple files

if graph_sweep
    sweep_results = Dict()
    graph_types=["erdos_renyi","scale_free","small_world","complete","empty"]
    prc("sweeping for graph type")
    for graph_type in graph_types
        params[:whichGraph]=graph_type
        sweep_results[graph_type] = param_sweeper()
    end
elseif run_repeat_sweep
    repeat_sweep()
else
    sweep_results = param_sweeper()
end
