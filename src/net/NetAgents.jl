include("../star/StarAgents.jl")
include("../star/Protos.jl")

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
