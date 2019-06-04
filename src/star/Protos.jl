
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
    balance::Float64
    alive::Bool
    arr_member_ids::Array{String, 1}
    ledger_transactions::Array{Transaction, 1}
end

struct NotEnoughSugarException <: Exception
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
    
    threshold = params[:proto_threshold]
    start_proto_id = begin
        if length(arr_protos) == 1 && arr_protos[1].proto_id == -1
            1
        else
            maximum([protoobj.proto_id for protoobj in arr_protos]) + 1
        end
    end
    next_proto_id = start_proto_id

    for agobj in arr_agents

        if agobj.a.proto_id == -1 && agobj.a.sugar_level > threshold

            arr_neighbors = fetch_eligible_neighbors(agobj, arr_agents,
                agent_environment)

            if !isempty(arr_neighbors)

                partner_agent = arr_neighbors[1]

                if partner_agent.a.proto_id == -1 ## (a)
                    agobj.a.proto_id = next_proto_id
                    deposit_amt = agobj.a.sugar_level - threshold
                    transaction1 = Transaction(deposit_amt,
                                               timeperiod, "deposit", 
                                               agobj.a.agent_id)
                    agobj.a.sugar_level -= deposit_amt

                    partner_agent.a.proto_id = next_proto_id
                    deposit_amt = partner_agent.a.sugar_level - threshold
                    transaction2 = Transaction(deposit_amt,
                                               timeperiod, "deposit", 
                                               partner_agent.a.agent_id)
                    partner_agent.a.sugar_level -= deposit_amt

                    prototn_obj = Proto(next_proto_id, 
                                       transaction1.transaction_amount +
                                       transaction2.transaction_amount,
                                       true,
                                       [agobj.a.agent_id, partner_agent.a.agent_id],
                                       [transaction1, transaction2])
                    next_proto_id += 1
                    push!(arr_protos, prototn_obj)

                else ## (b) and (c)
                    protot_obj = fetch_specific_proto_obj(arr_protos, 
                                                        partner_agent.a.proto_id)
                    @assert protot_obj.proto_id > 0
                    agobj.a.proto_id = partner_agent.a.proto_id
                    deposit_amt = agobj.a.sugar_level - threshold
                    transaction1 = Transaction(deposit_amt,
                                               timeperiod, "deposit", 
                                               agobj.a.agent_id)
                    agobj.a.sugar_level -= deposit_amt
                    protot_obj.balance += deposit_amt
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
        if proto_obj.balance <= 0
            proto_obj.alive = false
            proto_obj.arr_member_ids = []
            push!(proto_obj.ledger_transactions, 
                  Transaction(0, timeperiod, "closure", "-"))
        end
    end
end

# Attempt to withdraw sugar from a desperately poor agent's proto. If enough
# funds are available for it to continue to survive this iteration, make that
# withdrawal and return silently. Otherwise, leave its proto alone and throw a
# NotEnoughSugarException.
function withdraw_from_proto!(agobj, arr_protos, timeperiod)
    probj = fetch_specific_proto_obj(arr_protos,
                                agobj.a.proto_id)
    @assert agobj.a.proto_id > 0

    needed_amount = agobj.a.metabolic_rate - agobj.a.sugar_level
    if agobj.a.proto_id > 0 && (probj.balance > needed_amount)
        agobj.a.sugar_level = 0
        probj.balance -= needed_amount
        push!(probj.ledger_transactions,
              Transaction(needed_amount, timeperiod, "withdrawal",
                          agobj.a.agent_id))
    else
        throw(NotEnoughSugarException())
    end
end
