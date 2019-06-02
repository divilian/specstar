using RCall
@rlibrary ineq
using Gadfly
using LightGraphs
using GraphPlot, Compose
using ColorSchemes, Colors
using Random
using Distributions
using Cairo, Fontconfig
using DataFrames

include("NetAgents.jl")
include("../star/misc.jl")

# Run the SPECnet simulation once, for the set of parameters in the global
# variable "params".
function specnet()
    pri("SPECnet simulation parameters:")
    for (param, val) in params
        pri("   $(param) = $(val)")
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

    while !is_connected(graph)
        global graph = choose_graph()

        pri("Not connected; regenerating...")
        graph = choose_graph()
    end

    suglvl_distrib = DiscreteUniform(1, params[:init_sg_lvl])
    metabolic_distribution = DiscreteUniform(1, params[:metabolic_rate])

    global AN = Dict{NetAgent,Any}(
        NetAgent(
                AGENT_IDS[k],
                rand(metabolic_distribution),
                rand(suglvl_distrib),
				
                true,-1)
				
            =>k for k in 1:params[:N])

    ## the following is a hack; see comment in ../scape/run-simulation.jl
    arr_protos = [Proto(-1, -1, false, ["-"], [Transaction(-1, -1, "", "-")])]


    # (Erase old images.)
    rm("$(tempdir())/graph"*".png", force=true)
    rm("$(tempdir())/graph"*".svg", force=true)
    rm("$(tempdir())/wealth"*".png", force=true)
    rm("$(tempdir())/wealth"*".svg", force=true)
    rm("$(tempdir())/GiniPlot.png", force=true)
    locs_x, locs_y = nothing, nothing


    println("Iterations:")

    for iter in 1:params[:num_iters]

        global graph, locs_x, locs_y

        if verbosity == 1
            if iter % 10 == 0 println(iter) else print(".") end
        else
            println(" --- Iteration $(iter) of $(params[:num_iters]) ---")
        end

        if locs_x == nothing
            locs_x, locs_y = spring_layout(graph)
        else
            locs_x, locs_y = spring_layout(graph, locs_x, locs_y)
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
								#used to be -0.5 temp solution          
		    ag.a.sugar_level += (rand(Float16)) * params[:salary_range]
		    ag.a.sugar_level -= ag.a.metabolic_rate;

        end

        dying_agents =
            [ ag for ag in keys(AN)
                if ag.a.sugar_level < 0 && ag.a.alive ]
        for dying_agent in dying_agents
            try
                in_the_hole = ag.a.sugar_level
                withdraw_from_proto!(dying_agent, arr_protos, iter)
                prd("Agent $(dying_agent) got some money! " *
                    "(had $(in_the_hole), " *
                    "now has $(dying_agent.a.sugar_level), " *
                    "proto still has " *
                    "$(arr_protos[dying_agent.a.proto_id].balance))")
            catch exc
                if isa(exc, NotEnoughSugarException)
                    prd("Agent $(dying_agent) died!! " *
                        "(needed -$(dying_agent.a.sugar_level), only had " *
                        "$(arr_protos[dying_agent.a.proto_id].balance) in proto)")
                else
                    prd("Agent $(dying_agent) died!! (no proto)")
                end
                kill_agent(dying_agent)
            end
        end

        #adding current gini index to ginis array
        wealthArray=[]
        empty!(wealthArray)
        for ag in keys(AN)
			if ag.a.sugar_level>=-400
                push!(wealthArray,ag.a.sugar_level)
            end
        end
        cGini=ineq(wealthArray, type="Gini")

        gIndex=convert(Float16,cGini)
        push!(ginis,gIndex)

    end   # End main simulation for loop

    # Collect results in DataFrame.
    results = DataFrame(
        agent = [ ag.a.agent_id for ag in keys(AN) ],
        sugar = [ ag.a.sugar_level for ag in keys(AN) ],
        proto_id = [ ag.a.proto_id for ag in keys(AN) ]
    )

    if params[:make_anims]
        #drawing the gini index plot
        giniPlot=plot(x=1:params[:num_iters],y=ginis, Geom.point, Geom.line,
            Guide.xlabel("Iteration"), Guide.ylabel("Gini Index"))
        draw(PNG("$(tempdir())/GiniPlot.png"), giniPlot)

        println("Building wealth animation (be unbelievably patient)...")
        run(`convert -delay $(params[:animation_delay]) $(joinpath(tempdir(),"wealth"))"*".svg $(joinpath(tempdir(),"wealth.gif"))`)
        println("Building graph animation (be mind-bogglingly patient)...")
        run(`convert -delay $(params[:animation_delay]) $(joinpath(tempdir(),"graph"))"*".svg $(joinpath(tempdir(),"graph.gif"))`)
    end
 
    println("...ending SPECnet.")

    return sort(results, :agent)
end


################################ functions ################################

# Mark the agent "dead" whose agent number is passed. This involves
# surgically removing it from the graph, adjusting the agent-to-node mappings,
# and deleting it from the list of last-frame's plot coordinates.
function kill_agent(dying_agent)
    global graph, AN, locs_x, locs_y
    dying_node = AN[dying_agent]
    deleteat!(locs_x, dying_node)
    deleteat!(locs_y, dying_node)

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
    end

    if params[:whichGraph]=="scale_free"
        graph = LightGraphs.SimpleGraphs.static_scale_free(
            params[:N], params[:SF_edges], params[:SF_degree])
    end

    if params[:whichGraph]=="small_world"
        graph = LightGraphs.SimpleGraphs.watts_strogatz(
            params[:N], params[:λ], params[:SW_prob])
    end
    return graph
end

ginis=[]
rev_dict(d) = Dict(y=>x for (x,y) in d)


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
                    params[:num_iters]/10)
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
    run(`mogrify -format svg -gravity South -pointsize 15 -annotate 0 "Iteration $(iter) of $(params[:num_iters])"  $(joinpath(tempdir(),"graph"))$(lpad(string(iter),3,'0')).png`)
    run(`mogrify -format svg $(joinpath(tempdir(),"wealth"))$(lpad(string(iter),3,'0')).png`)

end
