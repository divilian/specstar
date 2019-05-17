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

function Base.show(io::IO, na::NetAgent)
    show(io, na.a)
end


function fetch_eligible_neighbors(agobj, arr_agents,
    specnet_graph::SimpleGraph{Int64})

    """
    For a given agent, identify agents that have a connection to it and
    whose sugar_level > proto_threshold and return them.
    """

    return shuffle([n for n in neighbors(specnet_graph, AN[agobj])
        if n.a.sugar_level > params[:proto_threshold]])
end ## end fetch_eligible_neighbors

