include("../star/StarAgents.jl")
include("../star/Protos.jl")

using LightGraphs

mutable struct NetAgent
    a::StarAgent

    function NetAgent(
        agent_id::String,
        metabolic_rate::Int64,
        sugar_level,
        alive::Bool,
        proto_id::Int64)
        starAgent = StarAgent(agent_id, metabolic_rate, sugar_level, alive,
            proto_id)
        return new(starAgent)
    end
end

function Base.hash(na::NetAgent, h::UInt)
    hash(na.a) + h
end

function Base.isequal(na1::NetAgent, na2::NetAgent)
    na1.a == na2.a
end

function Base.show(io::IO, na::NetAgent)
    show(io, na.a)
end


function fetch_eligible_neighbors(agobj::NetAgent,
    arr_agents::Array{NetAgent,1},
    specnet_graph::SimpleGraph{Int64})

    """
    For a given agent, identify some agents that it "may be interested in,"
    and whose sugar_level > proto_threshold, and return them.

    The phrase "may be interested in" refers to the openness parameter. Our
    selection algorithm is as follows:
        1. Identify all eligible graph-neighbors of this agent.
        2. For each of them, with probability openness, replace that
            graph-neighbor with an agent drawn from the graph at
            large, if such an eligible one exists.
    """

    nodes_to_agents = rev_dict(AN)

    neighbor_agents = [ nodes_to_agents[n]
        for n ∈  neighbors(graph, AN[agobj])
        if nodes_to_agents[n].a.sugar_level ≥ params[:proto_threshold] ]

    other_eligible_agents = [ ag for ag ∈  keys(AN)
        if AN[ag] ∉  neighbors(graph, AN[ag])  &&
        ag.a.sugar_level ≥ params[:proto_threshold] ]

    agents_to_return = neighbor_agents

    for i ∈  1:length(agents_to_return)
        if length(other_eligible_agents) == 0
            break
        end
        if rand(Float16) < params[:openness]
            replacement_idx = rand(1:length(other_eligible_agents))
            agents_to_return[i] = other_eligible_agents[replacement_idx]
            deleteat!(other_eligible_agents, replacement_idx)
            prd("Agent $(agobj.a.agent_id) considers $(agents_to_return[i].a.agent_id) at-large")
        else
            prd("Agent $(agobj.a.agent_id) considers $(agents_to_return[i].a.agent_id) from neighbors")
        end
    end

    prd("Returning: $(agobj.a.agent_id) -- $(agents_to_return)")
    return shuffle(agents_to_return)
end ## end fetch_eligible_neighbors

