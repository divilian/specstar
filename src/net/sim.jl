using RCall
@rlibrary DescTools
using Gadfly
using LightGraphs
using GraphPlot, Compose
using ColorSchemes, Colors
using Random
using Distributions
using Cairo, Fontconfig
using DataFrames
using Glob

include("NetAgents.jl")
include("../star/misc.jl")
include("setup_params.jl")

mutable struct SimState
    graph::SimpleGraph
    AN::Dict{NetAgent,Any}
    arr_protos::Array{Proto,1}
end

# Run the SPECnet simulation once, for the set of parameters in the global
# variable "params".
# Any keyword arguments provided will override settings in "params". For
# instance:
#     specnet()
#     specnet(max_iters=20)
#     specnet(max_iters=20, N=100)
function specnet(;additional_params...)

    @assert all_parameters_legit(additional_params)
    merge!(params, Dict(additional_params))

    pri("SPECnet simulation parameters:")
    for param in sort(collect(keys(params)), by=x->lowercase(string(x)))
        pri("   $(param) = $(params[param])")
    end

    global locs_x, locs_y

    Random.seed!(params[:random_seed])


    # Note on representation:
    #
    # An "agent" is identified by a string (see AGENT_IDS variable in
    # StarAgent.jl.) An agent's letter *always stays the same*, even when
    # agents die.
    #
    # The word "node," on the other hand, refers to an integer in the range
    # 1:N_t, which is the current number of living agents at time t. The node
    # numbers will always occupy the entire contiguous sequence 1:N_t (with no
    # "holes").
    #
    # Example: suppose N_0 is 5, so there are 5 agents initially. At the start
    # of the sim, the agents are A, B, C, D, E, and the nodes are 1, 2, 3, 4,
    # 5. If agents B and D die, then the living agents would be A, C, E, and
    # the corresponding nodes would be 1, 2, 3.
    #
    # The global variable AN is an (always currently maintained) mapping from
    # agent numbers to node numbers. To find the graph node that corresponds to
    # a particular agent a, do this:
    #
    #     AN[a]
    #
    # To find the agent that corresponds to a particular node n, do this:
    #
    #     [ a for a in keys(AN) if AN[a] == n ][1]
    #
    # (Yes, this reverse lookup is not efficient. If this turns out to be a
    # problem, will implement two dictionaries, one in each direction, and keep
    # them ruthlessly in sync.)


    println("Run SPECnet...")


    # The initial social network.
    global graph = choose_graph()

    arr_dead_protos=[]

    suglvl_distrib = DiscreteUniform(1, params[:init_sg_lvl])
    white_noise_distrib = Normal(0, params[:white_noise_intensity])

    global AN = Dict{NetAgent,Any}(
        NetAgent(
                AGENT_IDS[k],
                params[:metabolic_rate],  # constant for now
                rand(suglvl_distrib),
                true,-1)

            =>k for k in 1:params[:N])

    ## the following is a hack; see comment in ../scape/run-simulation.jl
    arr_protos = [Proto(-1, -1, false, ["-"], [Transaction(-1, -1, "", "-")])]


    # (Erase old images.)
    save_dir = pwd()
    cd("$(tempdir())")
    rm.(glob("graph*.png"))
    rm.(glob("graph*.svg"))
    rm.(glob("wealth*.png"))
    rm.(glob("wealth*.svg"))
    rm.(glob("GiniPlot.png"))
    rm.(glob("final_wealth_histogram.png"))
    cd(save_dir)

    locs_x, locs_y = nothing, nothing


    println("Iterations:")


    # The number of stage 3 iterations run so far.
    local starvation_timer = 0

    # The actual total number of iterations run. Assume pessimistically that it
    # will be terminated by the iteration limit; if the stopping condition is
    # encountered before that time, the loop will break out and this variable
    # will be set to a lower number at that time.
    local total_iters = params[:max_iters]

    # Each row will have three values: the low CI, the value, and the high CI.
    local ginis=Array{Float16}(undef,params[:max_iters],3)
    local stages=[]
    local nums_agents=[]

    global overall_results
    for iter in 1:params[:max_iters]

        global graph, locs_x, locs_y

        push!(stages, get_stage(SimState(graph, AN, arr_protos)))
        @assert stages[end] ∈  [1,2,3]
        if stages[end] == 1
            print("-")
        elseif stages[end] == 2
            print("+")
        elseif stages[end] == 3
            print("#")
        end
        if iter % 10 == 0 println(iter) end

        if stages[end] == 3
            if starvation_timer == 0
                # Collect stats on the model right before the starvation period.
                overall_results = collect_stats(SimState(graph, AN, arr_protos))
            end
            starvation_timer += 1
            if starvation_timer == params[:starvation_period]
                total_iters = iter - 1
                pop!(stages)    # (Remove what we prematurely added.)
                break   # End simulation after starvation.
            end
        end

        push!(nums_agents, length([ a for a ∈  keys(AN) if a.a.alive ]))

        if params[:make_anims]
            if locs_x == nothing
                locs_x, locs_y = spring_layout(graph)
            else
                locs_x, locs_y = spring_layout(graph, locs_x, locs_y)
            end
        end


        # Plot the graph at the *start* of each iteration (i.e., before any
        # activity has occurred) rather than the end.
        if params[:make_anims]
            plot_iteration_graphs(iter)
        end

        form_possible_protos!(collect(keys(AN)), graph, arr_protos, iter)

        [ assert_no_dead_neighbors(graph, ag) for ag ∈  keys(AN) ]

        # TODO: make agents graph neighbors who were not previously graph
        # neighbors but who are now in a common proto.


        # Payday!
        for ag in keys(AN)
            if starvation_timer == 0
                ag.a.sugar_level += params[:salary]
            end
            ag.a.sugar_level -= ag.a.metabolic_rate
            ag.a.sugar_level += rand(white_noise_distrib)
        end

        #Makes the agents deposit any excess sugar into their proto
        for ag in keys(AN)
            if ag.a.proto_id ≠ -1
                current = fetch_specific_proto_obj(arr_protos, ag.a.proto_id)
                if ag.a.sugar_level > params[:proto_threshold]
                    deposit_amt = ag.a.sugar_level - params[:proto_threshold]
                    transaction = Transaction(deposit_amt,
                                               iter, "deposit",
                                               ag.a.agent_id)
                    ag.a.sugar_level -= deposit_amt
                    current.balance += deposit_amt
                    push!(current.ledger_transactions, transaction)
                end
            end
        end

        dying_agents =
            [ ag for ag in keys(AN)
                if ag.a.sugar_level < 0 && ag.a.alive ]
        for dying_agent in dying_agents
            try
                in_the_hole = dying_agent.a.sugar_level
                withdraw_from_proto!(dying_agent, arr_protos, iter)
                prd("Agent $(dying_agent) got some money! " *
                    "(had $(in_the_hole), " *
                    "now has $(dying_agent.a.sugar_level), " *
                    "proto still has " *
                    "$(arr_protos[dying_agent.a.proto_id].balance))")
            catch exc
                if isa(exc, NotEnoughSugarException)
                    prd("Agent $(dying_agent) died!! " *
                        "(needed $(-dying_agent.a.sugar_level), only had " *
                        "$(arr_protos[dying_agent.a.proto_id].balance) in proto)")
                else
                    prd("Agent $(dying_agent) died!! (no proto)")
                end
                kill_agent(dying_agent)
            end
        end

        wealthArray=[]
        empty!(wealthArray)
        for ag in keys(AN)
            if ag.a.sugar_level>=-400
                push!(wealthArray,ag.a.sugar_level)
            end
        end

        kill_sugarless_protos(SimState(graph, AN, arr_protos), arr_dead_protos)

        # The R function Gini() from DescTools returns a vector of three values
        #   if conf.level is not NA: (1) the estimate, (2) the lower bound of
        #   the CI, and (3) the upper bound.
        # If we're going to end up making a single-sim plot, use Gini's CI's
        #   (which use bootstrapping, and are therefore time-consuming) to get
        #   all three. Otherwise, just get (1).
        if params[:make_sim_plots]
            if length(wealthArray) > 0
                rGini=Gini(wealthArray; Symbol("conf.level")=>.95,
                    :R=>params[:num_boot_samples])
                ginis[iter,:] = [ convert(Float16,r) for r in rGini ]
            else
                ginis[iter,:] = [ NaN, NaN, NaN ]
            end
        else
            rGini=Gini(wealthArray)
            ginis[iter,1] = convert(Float16,rGini)
        end

    end   # End main simulation for loop

    terminated_early = false
    if total_iters == params[:max_iters]
        pri("\n*WARNING* sim terminated early at iteration $(total_iters)!\n")
        terminated_early = true
    end

    # Collect results in DataFrame.
    agent_results = DataFrame(
        agent = [ ag.a.agent_id for ag in keys(AN) ],
        sugar = [ ag.a.sugar_level for ag in keys(AN) ],
        proto_id = [ ag.a.proto_id for ag in keys(AN) ]
    )
    iter_results = DataFrame(
        iter = 1:length(stages),
        gini_lowCI = ginis[1:length(stages),2],
        gini = ginis[1:length(stages),1],
        gini_highCI = ginis[1:length(stages),3],
        stage = stages,
        num_agents = nums_agents
    )

    # Collect stats on the model after the starvation period.
    starvation_results = collect_stats(SimState(graph, AN, arr_protos))
    if !isdefined(Main, :overall_results)
        # Not totally sure what to do in this WARNING case.
        overall_results = copy(starvation_results)
    end

    if params[:make_sim_plots]
        plot_gini_livingfrac_over_time(iter_results,
            [:proto_threshold, :salary, :white_noise_intensity])
        plot_final_wealth_hist(SimState(graph, AN, arr_protos))
    end

    if params[:make_anims]
        println("Building wealth animation (be unbelievably patient)...")
        run(`convert -delay $(params[:animation_delay]) $(joinpath(tempdir(),"wealth"))"*".svg $(joinpath(tempdir(),"wealth.gif"))`)
        println("Building graph animation (be mind-bogglingly patient)...")
        run(`convert -delay $(params[:animation_delay]) $(joinpath(tempdir(),"graph"))"*".svg $(joinpath(tempdir(),"graph.gif"))`)
    end
    println("\n...ending SPECnet.")

    return Dict(:agent_results => sort(agent_results, :agent),
        :iter_results => iter_results,
        :overall_results => overall_results,
        :starvation_results => starvation_results,
        :terminated_early => terminated_early)
end


################################ functions ################################

# Remove all protos that have zero (or less) sugar balance, adding them to the
#   arr_dead_protos variable and setting all its members back to a -1 ("no
#   proto") proto_id.
function kill_sugarless_protos(sim_state::SimState, arr_dead_protos)
    index=length(sim_state.arr_protos)
    while index>0
        if(sim_state.arr_protos[index].proto_id ≠ -1 &&
                sim_state.arr_protos[index].balance<=0)
            pri("Killing proto $(sim_state.arr_protos[index].proto_id)")
            for id in sim_state.arr_protos[index].arr_member_ids
                for ag in keys(AN)
                    if ag.a.agent_id==id
                        ag.a.proto_id=-1
                    end
                end
            end
            push!(arr_dead_protos,sim_state.arr_protos[index])
            deleteat!(sim_state.arr_protos,index)
        end
        index-=1
    end
end

# Mark the agent "dead" whose agent number is passed. This involves
# surgically removing it from the graph, adjusting the agent-to-node mappings,
# and deleting it from the list of last-frame's plot coordinates.
function kill_agent(dying_agent)
    global graph, AN, locs_x, locs_y
    dying_node = AN[dying_agent]
    if locs_x ≠ nothing
        deleteat!(locs_x, dying_node)
        deleteat!(locs_y, dying_node)
    end

    neighbor_nodes = neighbors(graph, dying_node)
    while length(neighbor_nodes) > 0
        rem_edge!(graph, dying_node, neighbor_nodes[1])
        neighbor_nodes = neighbors(graph, dying_node)
    end
    AN = Dict{NetAgent,Any}(a=>
        (a==dying_agent ? Nothing :
            (AN[a] == nv(graph) ? AN[dying_agent] : AN[a]))
        for (a,n) in AN)
    dying_agent.a.alive = false
    # TODO: the rem_vertices!() function that Simon Schoelly wrote returns
    # a map of old-to-new node numbers, and might be safer.
    rem_vertex!(graph, dying_node)
    pop!(AN,dying_agent)
end

# Return true if the agent is a member of any proto institution.
in_proto(agent) = agent.a.proto_id > 0

# Return true if the agent is a rich enough to joing a proto, and available
# to do so.
function eligible_for_proto(agent)
    return wealth_eligible(agent) && !in_proto(agent)
end

#return true if agent is wealthy enough to join proto
function wealth_eligible(agent)
    return agent.a.sugar_level > params[:proto_threshold]
end

function assert_no_dead_neighbors(graph, agent)
    node = AN[agent]
    nodes_to_agents = rev_dict(AN)
    if !all([ nodes_to_agents[n].a.alive for n in neighbors(graph, node)])
        error("Dead neighbor of $(agent)!\n" *
            "neighbors: $(nodes_to_agents[neighbors(graph, node)].a.agent_id)")
    end
end

let possible_colors = Random.shuffle(ColorSchemes.rainbow.colors)
    global compute_colors
    function compute_colors()
        agents_in_node_order = [
            [ a for a in keys(AN)
                        if AN[a] == node ][1] for node in 1:nv(graph) ]
        return [ in_proto(ag) ? possible_colors[ag.a.proto_id] :
            (!ag.a.alive ? colorant"pink" : colorant"lightgrey")
                        for ag in agents_in_node_order ]
    end
end

function choose_graph()
    if params[:whichGraph]=="erdos_renyi"
        ER_prob=params[:λ]/params[:N]
        graph = LightGraphs.SimpleGraphs.erdos_renyi(
            params[:N], ER_prob)

    elseif params[:whichGraph]=="scale_free"
        graph = LightGraphs.SimpleGraphs.static_scale_free(
            params[:N], params[:SF_edges], params[:SF_degree])

    elseif params[:whichGraph]=="small_world"
        graph = LightGraphs.SimpleGraphs.watts_strogatz(
            params[:N], params[:k], params[:β])

    elseif params[:whichGraph]=="complete"
        graph = LightGraphs.SimpleGraphs.CompleteGraph(params[:N])

    elseif params[:whichGraph]=="empty"
        graph = LightGraphs.SimpleGraphs.SimpleGraph(params[:N], 0)

    else
        throw(DomainError("Invalid graph type \"$(params[:whichGraph])\""))
    end
    return graph
end

rev_dict(d) = Dict(y=>x for (x,y) in d)

function plot_final_wealth_hist(sim_state::SimState)
    final_wealths=[]
    final_proto_wealths_alive=[]
    final_proto_wealths=[]
    [push!(final_wealths,ag.a.sugar_level) for ag in keys(AN) ]
    for p in sim_state.arr_protos
        if(p.balance!=-1)
            push!(final_proto_wealths_alive, p.balance)
            push!(final_proto_wealths, p.balance)
        else
            push!(final_proto_wealths, p.balance)
        end
    end

    final_wealthp = plot(
        x=final_wealths,
        Geom.histogram(density=true, bincount=20),
        Theme(background_color=colorant"white"),
        Guide.xlabel("Agent Wealth"),
        Guide.ylabel("Density of agents"))
    final_proto_wealthp = plot(
        x=final_proto_wealths,
        Geom.histogram(density=true, bincount=20),
        Theme(background_color=colorant"white"),
        Guide.xlabel("Proto Wealth"),
        Guide.ylabel("Density of protos"))
    final_proto_wealth_alivep = plot(
        x=final_proto_wealths_alive,
        Geom.histogram(density=true, bincount=20),
        Theme(background_color=colorant"white"),
        Guide.xlabel("Proto Wealth"),
        Guide.ylabel("Density of protos"))
    draw(PNG("$(tempdir())/final_agent_wealth_histogram.png"),final_wealthp)
    draw(PNG("$(tempdir())/final_proto_wealth_histogram_alive.png"),final_proto_wealth_alivep)
    draw(PNG("$(tempdir())/final_proto_wealth_histogram.png"),final_proto_wealthp)

end

# Plot the Gini coefficient, and the fraction of agents still alive, over time.
# The callouts parameter should be a list of symbols which should appear in the
# title's plot.
function plot_gini_livingfrac_over_time(iter_results, callouts)

    # Remove rows with NaNs.
    for col in names(iter_results)
        iter_results[col] = map(x->isnan(x) ? missing : x, iter_results[col])
    end
    iter_results = iter_results[completecases(iter_results),:]

    stage_names = Dict(1=>"1", 2=>"2", 3=>"3")
    iter_results.stage = [ stage_names[s] for s in iter_results.stage ]
    iter_results.num_agents /= maximum(iter_results.num_agents)
    giniPlot=plot(iter_results,
        layer(
            x=:iter, y=:gini_lowCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=:iter, y=:gini,
            Geom.line, Geom.point,
            color=:stage,
        ),
        layer(
            x=:iter, y=:gini_highCI,
            Geom.line,
            Theme(default_color=colorant"lightgray")
        ),
        layer(
            x=:iter, y=:gini, ymin=:gini_lowCI, ymax=:gini_highCI,
            Geom.ribbon,
            Theme(default_color=colorant"lightblue")
        ),
        layer(
            x=:iter, y=:num_agents,
            Geom.line,
            Theme(default_color=colorant"black")
        ),
        Guide.xlabel("Iteration"), Guide.ylabel("Gini / fraction living"),
        Theme(background_color=colorant"white"),
        Scale.color_discrete_manual("blue","green","red",
            levels=["1","2","3"]),
        Guide.title(join([ "$(c)=$(params[c])\n" for c in callouts ])))
    draw(PNG("$(tempdir())/GiniPlot.png"), giniPlot)
end

function plot_iteration_graphs(iter)

    # Plot graph for this iteration.
    colors = compute_colors()

    remember_layout = x -> spring_layout(x, locs_x, locs_y)

    labels_to_plot = map(node->
            [ ag.a.agent_id for ag in keys(AN) if AN[ag] == node ][1],
        1:length(AN))
    wealths_to_plot = map(node->
            [ ag.a.sugar_level for ag in keys(AN)
                                    if AN[ag] == node ][1],
        1:length(AN))
    graphp = gplot(graph,
        layout=remember_layout,
        nodelabel=labels_to_plot,
        NODESIZE=.08,
        nodesize=ifelse.(wealths_to_plot .> 0,
                         wealths_to_plot*4,
                         maximum(wealths_to_plot)*2),
        nodestrokec=colorant"grey",
        nodestrokelw=.5,
        nodefillc=colors)
    draw(PNG("$(tempdir())/graph$(lpad(string(iter),3,'0')).png"),
        graphp)

    # Plot wealth histogram for this iteration.
    in_proto_wealths=[]
    not_in_proto_wealths=[]
    [  in_proto(ag) ?
            push!(in_proto_wealths,ag.a.sugar_level) :
            push!(not_in_proto_wealths,ag.a.sugar_level)
        for ag in keys(AN) ]

    wealthp = plot(
        layer(x=in_proto_wealths,
            Geom.histogram(density=true, bincount=20),
            Theme(default_color=Colors.RGBA(255,165,0, 0.6))),

        layer(x=not_in_proto_wealths,
            Geom.histogram(density=true, bincount=20),
            Theme(default_color=Colors.RGBA(255,0,255,0.6))),
            Guide.xlabel("Wealth"),
            Guide.ylabel("Density of agents"),
            Guide.title("Wealth distribution at iteration $(iter)"),
            # Hard to know what to set the max value to.
            Scale.x_continuous(minvalue=0,
                maxvalue=params[:init_sg_lvl]*
                    params[:max_iters]/10),

        Theme(background_color=colorant"white"),
    )

    draw(PNG("$(tempdir())/wealth$(lpad(string(iter),3,'0')).png"),
        wealthp)

    colors = compute_colors()

    remember_layout = x -> spring_layout(x, locs_x, locs_y)

    labels_to_plot = map(
        node->[ag.a.agent_id for ag in keys(AN) if AN[ag]==node][1],
        1:length(AN))
    wealths_to_plot = map(
        node->[ag.a.sugar_level for ag in keys(AN) if AN[ag]==node][1],
        1:length(AN))
    graphp = gplot(graph,
        layout=remember_layout,
        nodelabel=labels_to_plot,
        NODESIZE=.08,
        nodesize=ifelse.(wealths_to_plot .> 0,
                         wealths_to_plot*4,
                         maximum(wealths_to_plot)*2),
        nodestrokec=colorant"grey",
        nodestrokelw=.5,
        nodefillc=colors)
    draw(PNG("$(tempdir())/graph$(lpad(string(iter),3,'0')).png"),
                                                            graphp)

    #iteration label for svg files
    run(`mogrify -format svg -gravity South -pointsize 15 -annotate 0 "Iteration $(iter) of up to $(params[:max_iters])"  $(joinpath(tempdir(),"graph"))$(lpad(string(iter),3,'0')).png`)
    run(`mogrify -format svg $(joinpath(tempdir(),"wealth"))$(lpad(string(iter),3,'0')).png`)

end


# Return an integer ∈  {1,2,3} for the stage that the simulation is currently in:
#   S1: No agent has yet created a proto.
#   S2: Agents are creating and joining protos.
#   S3: All agents are either (a) graph isolates, (b) dead, or (c) in a proto.
function get_stage(sim_state::SimState)
    if length(sim_state.arr_protos) == 1   # See arr_protos hack
        return 1
    else
        for agent in keys(AN)
            if length(neighbors(graph,AN[agent])) > 0  &&
                    agent.a.alive  &&
                    agent.a.proto_id == -1
                return 2
            end
        end
        return 3
    end
end

function collect_stats(s::SimState)
    total_proto_sz = sum([ length(p.arr_member_ids) for p ∈  s.arr_protos ])
    component_vertices = connected_components(s.graph)
    proto_average_size=total_proto_sz/length(s.arr_protos)
    num_agents_in_proto = sum([ ag.a.proto_id ≠ -1 for ag ∈  keys(AN) ])
    return Dict(:size_largest_comp => nv(s.graph) == 0 ? 0 :
            findmax(length.(component_vertices))[1][1],
        :num_comps => nv(s.graph) == 0 ? 0 :
            length(component_vertices),
        :average_proto_size => proto_average_size,
        :num_protos => length(s.arr_protos),
        :num_agents_in_proto => num_agents_in_proto,
        :num_living_agents => length([ a for a ∈  keys(s.AN) if a.a.alive ]),
   )
end

function all_parameters_legit(additional_params)
    if !all([ word in keys(params) for word in keys(additional_params) ])
        for word in keys(additional_params)
            if word ∉  keys(params)
                println("No such SPECnet parameter \"$(word)\".")
            end
        end
        return false
    end
    return true
end
