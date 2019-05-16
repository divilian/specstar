using RCall
@rlibrary ineq
using Gadfly
using LightGraphs
using GraphPlot, Compose
using ColorSchemes, Colors
using Random
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
    global AN = Dict{NetAgent,Any}(
        NetAgent(
                AGENT_IDS[k],
                0,
                rand(Float16) * params[:max_starting_wealth],true,0) 
            =>k for k in 1:params[:N])


    # (Erase old images.)
    rm("$(tempdir())/graph"*".png", force=true)
    rm("$(tempdir())/graph"*".svg", force=true)
    rm("$(tempdir())/wealth"*".png", force=true)
    rm("$(tempdir())/wealth"*".svg", force=true)
    rm("$(tempdir())/GiniPlot.png", force=true)
    locs_x, locs_y = nothing, nothing

    println("Iterations:")

    for iter in 1:params[:num_iter]

        if iter % 10 == 0 println(iter) else print(".") end

        global graph, locs_x, locs_y

        if locs_x == nothing
            locs_x, locs_y = spring_layout(graph)
        else
            locs_x, locs_y = spring_layout(graph, locs_x, locs_y)
        end

        agent1 = rand(keys(AN))
        assert_no_dead_neighbors(graph, agent1)
        if rand(Float16) < params[:openness]  ||
                length(neighbors(graph, AN[agent1])) == 0
            # Choose from the graph at large.
            agent2 = rand(filter(x->x!=agent1,keys(AN)))
            prd("$(agent1) encounters at-large $(agent1)")
        else
            # Choose from a neighbor.
            node2 = rand(neighbors(graph,AN[agent1]))
            nodes_to_agents = rev_dict(AN)
            agent2 = nodes_to_agents[node2]
            prd("$(agent1) encounters neighbor $(agent2)")
        end
        assert_no_dead_neighbors(graph, agent2)

        if eligible_for_proto(agent1) && eligible_for_proto(agent2)
            form_proto(agent1, agent2)
            # Since they're forming a proto, they also become socially
            # connected (if they weren't already.)
            if !has_edge(graph, AN[agent1], AN[agent2])
                add_edge!(graph, AN[agent1], AN[agent1])
            end
        end

        if wealth_eligible(agent1) && wealth_eligible(agent2)
            if in_proto(agent1)&&!in_proto(agent2)
                agent2.a.proto_id=agent1.a.proto_id
            end

            if in_proto(agent2)&&!in_proto(agent1)
                agent1.a.proto_id=agent2.a.proto_id
            end
        end

        # Plot graph for this iteration.
        if params[:make_anims]
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

            # Plot wealth histogram for this iteration.
            wealthp = plot(x=[ ag.a.sugar_level for ag in keys(AN) ],
                Geom.histogram(density=true, bincount=20),
                Guide.xlabel("Wealth"),
                Guide.ylabel("Density of agents"),
                Guide.title("Wealth distribution at iteration $(iter)"),
                # Hard to know what to set the max value to.
                Scale.x_continuous(minvalue=0, maxvalue=
                    params[:max_starting_wealth]*params[:num_iter]/10), theme)

            draw(PNG("$(tempdir())/wealth$(lpad(string(iter),3,'0')).png"),
                                                                    wealthp)

            #iteration label for svg files
            run(`mogrify -format svg -gravity South -pointsize 15 -annotate 0 "Iteration $(iter) of $(params[:num_iter])"  $(joinpath(tempdir(),"graph"))$(lpad(string(iter),3,'0')).png`)
            run(`mogrify -format svg $(joinpath(tempdir(),"wealth"))$(lpad(string(iter),3,'0')).png`)
        end

        # Payday!
        for ag in keys(AN)
            ag.a.sugar_level += (rand(Float16) - .5) * params[:salary_range]
            if ag.a.proto_id > 0
                ag.a.sugar_level += (rand(Float16)*10)
            end
        end

        dying_agents =
            [ ag for ag in keys(AN) 
                if ag.a.sugar_level < 0 && ag.a.alive ]
        for dying_agent in dying_agents
            prd("Agent $(dying_agent) died!")
            kill_agent(dying_agent)
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
        wealth = [ ag.a.sugar_level for ag in keys(AN) ],
        proto_id = [ ag.a.proto_id for ag in keys(AN) ]
    )

    if params[:make_anims]
        #drawing the gini index plot
        giniPlot=plot(x=1:params[:num_iter],y=ginis, Geom.point, Geom.line,
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

# Form a new proto between two agents.
function form_proto(agent1, agent2)
    next_proto_id = maximum([ ag.a.proto_id for ag in keys(AN) ]) + 1
    agent1.a.proto_id = next_proto_id
    agent2.a.proto_id = next_proto_id
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
        graph = LightGraphs.SimpleGraphs.erdos_renyi(
            params[:N], params[:ER_prob])
    end

    if params[:whichGraph]=="scale_free"
        graph = LightGraphs.SimpleGraphs.static_scale_free(
            params[:N], params[:SF_edges], params[:SF_degree])
    end

    if params[:whichGraph]=="small_world"
        graph = LightGraphs.SimpleGraphs.watts_strogatz(
            params[:N], params[:SW_degree], params[:SW_prob])
    end
    return graph
end

ginis=[]
rev_dict(d) = Dict(y=>x for (x,y) in d)
###########################################################################


