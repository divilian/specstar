
using Random

export Proto, Transaction, form_possible_protos!, update_proto_statuses!,
    fetch_specific_proto_obj

mutable struct Transaction
    transaction_amount::Float64
    period::Int32
    type::String ## deposit or withdrawal or closure
    agent_id::String    
end

mutable struct Proto
    proto_id::Int64
    sugar_level::Float64
    alive::Bool
    arr_member_ids::Array{String, 1}
    ledger_transactions::Array{Transaction, 1}
end

function fetch_specific_proto_obj(arr_protos, proto_id)
    """
    Return a proto object from the given array of proto objects
    that has the given proto id.
    """
    return [probj for probj in arr_protos 
            if probj.proto_id == proto_id][1] 
end

function form_possible_protos!(arr_agents, agent_environment, arr_protos,
        timeperiod)

    """
    Called during each time period to enact association among agents to form
    protos. The following conditions exist:
    (a) none of the focal agent's neighbors are part of a proto
    (b) one of the focal agent's neighbors is part of a proto
    (c) two or more of the focal agent's neighbors are part of a proto.
    
    In the case of (a), a randomly chosen neighbor and the focal agent join to
    form a proto.
    In the case of (b), the focal agent joins the proto of the neighbor.
    In the case of (c), the focal agent joins the proto of a 
    randomly-chosen neighbor.
    
    Modifies focal and neighbor agents' proto_id field, if proto-
    formation or joining occurs.
    Reduces the focal agent and neighbor agents' sugarlevels by an amount 
    = sugarlevel - params[:threshold].

    Modifies the array of protos by adding, as needed, new proto struct.
    
    Returns an array of Proto structures to the caller.
    """
    # println("At the beginning form_possible_protos! the array of protos is:")
    # ## println(arr_protos)
    # println("Length of arr_protos is: ", string(length(arr_protos)))
    
    threshold = params[:proto_threshold]
    start_proto_id = begin
        if length(arr_protos) == 1 && arr_protos[1].proto_id == -1
            1
        else
            maximum([protoobj.proto_id for protoobj in arr_protos]) + 1
        end
    end
    ## println("The start_proto_id is: ", string(start_proto_id))
    for agobj in arr_agents
        ## update their proto_ids. 
        ## println("Current object:", string(agobj.a.agent_id), 
                # " its excess sugarlevel: ", agobj.a.sugar_level - threshold
        ## println("Enter enter")
        ## readline()
        if agobj.a.proto_id == -1 && agobj.a.sugar_level > threshold
            ## fetch neighbors
            ## println("Checking to see if ", string(agobj.a.agent_id), " can join an", " proto")
            arr_neighbors = fetch_eligible_neighbors(agobj, arr_agents,
                agent_environment)
            ## println("Here here here!")
            if !isempty(arr_neighbors)

                partner_agent = arr_neighbors[1]

                if partner_agent.a.proto_id == -1 ## (a)
                    ## println("No neighbor found with an existing membership, ",
                            # "so creating a new proto with agent: ", 
                            # string(partner_agent.a.agent_id))
                    agobj.a.proto_id = start_proto_id
                    agobj.a.sugar_level = agobj.a.sugar_level - threshold
                    transaction1 = Transaction(agobj.a.sugar_level - threshold,
                                               timeperiod, "deposit", 
                                               agobj.a.agent_id)

                    partner_agent.a.proto_id = start_proto_id
                    partner_agent.a.sugar_level = partner_agent.a.sugar_level - threshold
                    transaction2 = Transaction(partner_agent.a.sugar_level - threshold,
                                               timeperiod, "deposit", 
                                               partner_agent.a.agent_id)

                    prototn_obj = Proto(start_proto_id, 
                                       transaction1.transaction_amount +
                                       transaction2.transaction_amount,
                                       true,
                                       [agobj.a.agent_id, partner_agent.a.agent_id],
                                       [transaction1, transaction2])

                    push!(arr_protos, prototn_obj)
                    ## println("Added ", string(agobj.a.agent_id), " to a new ",
                            # "proto with id:", start_proto_id)

                else ## (b) and (c)
                    ## println("Inside the (b) and (c) conditional branch")
                    protot_obj = fetch_specific_proto_obj(arr_protos, 
                                                        partner_agent.a.proto_id)
                    @assert protot_obj.proto_id > 0
                    agobj.a.proto_id = partner_agent.a.proto_id
                    agobj.a.sugar_level = agobj.a.sugar_level - threshold
                    transaction1 = Transaction(agobj.a.sugar_level - threshold,
                                               timeperiod, "deposit", 
                                               agobj.a.agent_id)
                    push!(protot_obj.ledger_transactions, transaction1)
                    
                end ## (a), (b), (c)
            end ## !isempty(arr_neighbors)
        end ## agobj.a.sugar_level > threshold & agobj.a.proto_id > 0
    end ## for agobj    
end ## end of form_possible_protos

function update_proto_statuses!(arr_protos, timeperiod)
    """
    Goes through the array of protos and updates the status
    of each proto: if all of the member agents are dead,
    the proto's alive status is set to false, and a closing transaction
    is added, with a transaction amount of zero, date = current period, and
    agent_id = -999.
    """
    for proto_obj in arr_protos
        if proto_obj.sugar_level <= 0
            proto_obj.alive = false
            proto_obj.arr_member_ids = []
            push!(proto_obj.ledger_transactions, 
                  Transaction(0, timeperiod, "closure", "-"))
        end
    end
end
