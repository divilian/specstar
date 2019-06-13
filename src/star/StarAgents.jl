
mutable struct StarAgent
    agent_id::String
    metabolic_rate::Int64
    sugar_level::Float64
    alive::Bool
    proto_id::Int64
end

function Base.hash(sa::StarAgent, h::UInt)
    hash(sa.agent_id) + h
end

function Base.isequal(sa1::StarAgent, sa2::StarAgent)
    sa1.agent_id == sa2.agent_id
end


function Base.show(io::IO, sa::StarAgent)
    show(io, "Agent($(sa.agent_id))")
end

# For agent_ids.
uppers = collect('A':'Z')
lowers = collect('a':'z')
AGENT_IDS = [ string(x) for x in uppers ]
[ push!(AGENT_IDS, string(x)) for x in lowers ]
[ push!(AGENT_IDS, x) for x in [ "$a$b" for a in 'A':'Z' for b in 'A':'Z' ] ]
[ push!(AGENT_IDS, x) for x in [ "$a$b" for a in 'a':'z' for b in 'a':'z' ] ]

