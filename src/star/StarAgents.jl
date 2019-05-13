
module StarAgents

export StarAgent

mutable struct StarAgent
    agent_id::Int64
    metabolic_rate::Int64
    sugar_level::Float64
    alive::Bool
    proto_id::Int64
end

end
