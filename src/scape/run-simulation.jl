include("../star/Protos.jl")
include("ScapeAgents.jl")
include("Sugarscape.jl")
include("max-num-generator.jl")
include("../star/misc.jl")

using Statistics
using Random
using Distributions
using CSV
using DataFrames
using RCall

function specscape()
    """
    Run the SPECscape simulation once, for the set of parameters in the global
    variable "params".
    Returns a single row, consisting of all of the params + gini values
    of sugar across all the time periods.
    """
    prd("SPECscape simulation parameters:")
    for (param, val) in params
        prd("   $(param) = $(val)")
    end

    metabol_distrib =  DiscreteUniform(params[:metab_range_tpl][1], params[:metab_range_tpl][2])
    vision_distrib = DiscreteUniform(params[:vision_range_tpl][1], params[:vision_range_tpl][2])
    suglvl_distrib = DiscreteUniform(params[:suglvl_range_tpl][1], params[:suglvl_range_tpl][2])
    arr_agent_ginis = zeros(params[:num_iter])
    ## the following is a hack because creating an empty array of Array{Proto, 1}
    ## and adding Proto objects via push! is resulting in errors.
    ## So to add type-checking on arr_protos, we're going to initialize it with
    ## a dummy Proto object
    # arr_protos = Array{Proto, 1}
    arr_protos = [Proto(-1, -1, false, ["-"], [Transaction(-1, -1, "", "-")])]

    dict_objs = set_up_environment()
    sugscape_obj = dict_objs["sugscape_obj"]
    arr_agents = dict_objs["arr_agents"]

    for period in 1:params[:num_iter]
        for ind in shuffle(1:length(arr_agents))
            locate_move_feed!(arr_agents[ind], sugscape_obj, arr_agents, arr_protos, period)
        end
        regenerate_sugar!(sugscape_obj)
        perform_birth_inbound_outbound!(arr_agents, sugscape_obj, params[:birth_rate],
                                        params[:inbound_rate], params[:outbound_rate],
                                        vision_distrib, metabol_distrib,
                                        suglvl_distrib)
        form_possible_protos!(arr_agents, sugscape_obj,
                                    arr_protos, period)

        arr_agents = life_check!(arr_agents)
        @assert all([aggobj.a.alive for aggobj in arr_agents])
        # println("HERHEREHERE")
        # readline()
        update_occupied_status!(arr_agents, sugscape_obj)
        update_proto_statuses!(arr_protos, period)
        arr_agent_ginis[period] = compute_Gini(arr_agents)
        # println("No. of protos: ", length(arr_protos))
        # println("No. of agents: ", string(length(arr_agents)))
        # println("Finished time-step: ", string(period), "\n\n")
        ## println("Enter enter")
        ## readline()
    end## end of params[:num_iter] for loop
    return(arr_agent_ginis)
end ## end specscape()


# Take the global params variable, and return a dict with an appropriate
# sugarscape object and array of agents.
function set_up_environment()
    """
    Arguments used (from global params):
    :scape_side
    :scape_carry_cap
    :scape_growth_rate
    :pop_density
    :metab_range_tpl
    :vision_range_tpl
    :suglvl_range_tpl

    Returns: dictionary {sugscape object =>, arr_agents => }
    """
    ## Generate an empty sugarscape
    sugscape_obj = generate_sugarscape(
        params[:scape_side],
        params[:scape_growth_rate],
        params[:scape_carry_cap],
        3)
    stats = get_sugarscape_stats(sugscape_obj)

    no_agents = Int(ceil(params[:pop_density] * params[:scape_side]^2))

    metabol_distrib =  DiscreteUniform(
        params[:metab_range_tpl][1], params[:metab_range_tpl][2])
    vision_distrib = DiscreteUniform(
        params[:vision_range_tpl][1], params[:vision_range_tpl][2])
    suglvl_distrib = DiscreteUniform(
        params[:suglvl_range_tpl][1], params[:suglvl_range_tpl][2])


    # metabol_distrib = Uniform(metab_range_tpl[1], metab_range_tpl[2]);
    # vision_distrib = Uniform(vision_range_tpl[1], vision_range_tpl[2]);
    # suglvl_distrib = Uniform(suglvl_range_tpl[1], suglvl_range_tpl[2]);

    arr_poss_locations = sample([
        (x,y) for x in 1:params[:scape_side], y in 1:params[:scape_side]],
                                no_agents, replace=false)

    arr_agents = [ScapeAgent(AGENT_IDS[agg_id],
                        rand(metabol_distrib),
                        rand(suglvl_distrib),
                        true,
                        -1,
                        arr_poss_locations[agg_id][1],
                        arr_poss_locations[agg_id][2],
                        rand(vision_distrib))
                  for agg_id in 1:no_agents]

    ## mark as occupied the cells in sugarscape corresponding to the agents' locs
    for loc in arr_poss_locations
        sugscape_obj[loc[1], loc[2]].occupied = true
    end
    # println("Created a sugarscape of size: ",
    #         string(size(sugscape_obj)[1] * size(sugscape_obj)[2]))
    # println("Created ", string(length(arr_agents)), " agents.")
    return(Dict("sugscape_obj" => sugscape_obj,
                "arr_agents" => arr_agents))
end ## end of set_up_environment()


function run_all_param_ranges_for(givenseed)
    Random.seed!(givenseed)
    params_df = CSV.read("parameter-ranges-testing.csv")
    outfile_name = "outputs-new.csv"
    time_periods = 100

    temp_out = DataFrame(zeros(nrow(params_df), time_periods))
    names!(temp_out, Symbol.(["prd_"*string(i) for i in 1:time_periods]))

    out_df = DataFrame()
    for colname in names(params_df)
        out_df[Symbol(colname)] = params_df[Symbol(colname)]
    end

    for colname in names(temp_out)
        out_df[Symbol(colname)] = temp_out[Symbol(colname)]
    end

    for rownum in 1:nrow(params_df)

        ## println(get_sugarscape_stats(sugscape_obj))
        ## println("\n\n")
        # plot_sugar_concentrations!(sugscape_obj)

        # Create the global dict of parameters for this run of the sim.
        global params = Dict{Symbol,Any}(
            :scape_side => params_df[rownum, :Side],
            :scape_carry_cap => params_df[rownum, :Capacity],
            :scape_growth_rate => params_df[rownum, :RegRate],
            :metab_range_tpl => (1, params_df[rownum, :MtblRate]),
            :vision_range_tpl => (1, params_df[rownum, :VsnRng]),
            :suglvl_range_tpl => (1, params_df[rownum, :InitSgLvl]),
            :pop_density => params_df[rownum, :Adensity],
            :birth_rate => params_df[rownum, :Birthrate],
            :inbound_rate => params_df[rownum, :InbndRt],
            :outbound_rate => params_df[rownum, :OtbndRt],
            :proto_threshold => params_df[rownum, :Threshold],
            :num_iter => time_periods,
        )

        ## next, actually run the simulation for this set of parameters:
        ## move the agents, have them consume sugar, reduce the sugar in
        ## sugscape cells, regrow the sugar....and collect the array of gini
        ## coeffs
        arr_agent_ginis = specscape()

        # for colname in names(params_df)
        #     out_df[rownum, Symbol(colname)] = params_df[rownum, Symbol(colname)]
        # end

        for colnum in ncol(params_df)+1 : ncol(out_df)
            out_df[rownum, colnum] = arr_agent_ginis[colnum - ncol(params_df)]
        end


        ## create a row
        println("Finished combination $rownum")
        # println("Here's the out_df")
        # println(out_df)
        # readline()
    end #end iterate over param rows

    return(out_df)
end ## run_all_param_ranges_for

# run_all_param_ranges_for() |> CSV.write("output.csv")
arr_seeds = [10, 80085, 4545, 4543543535, 87787765, 63542, 34983, 596895, 2152, 434];
for seednum in arr_seeds
    println("Beginning simulation run for seed: ", string(seednum))
    fname = "outputfile-" * string(seednum) * ".csv"
    run_all_param_ranges_for(seednum) |> CSV.write(fname)
    println("Completed run for seed: ", string(seednum))
end
