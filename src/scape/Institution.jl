using Random

mutable struct Transaction
    transaction_amount::Float64
    period::Int32
    type::String ## deposit or withdrawal or closure
    agent_id::Int64    
end

mutable struct Institution
    institution_id::Int64
    sugar_level::Float64
    alive::Bool
    arr_member_ids::Array{Int64, 1}
    ledger_transactions::Array{Transaction, 1}
end

function fetch_specific_inst_obj(arr_institutions, inst_id)
    """
    Return an institution object from the given array of institution objects
    that has the given institution id.
    """
    return [inobj for inobj in arr_institutions 
            if inobj.institution_id == inst_id][1] 
end

function fetch_neighbors(agobj, arr_agents, sugscape_obj, threshold)
    """
    For a given agent, identify agents that are in the n-s-e-w cells
    whose sugar_level > threshold and return them.
    """

    neighbor_locs = [(agobj.location_x - 1, agobj.location_y), 
                     (agobj.location_x + 1, agobj.location_y),
                     (agobj.location_x, agobj.location_y - 1), 
                     (agobj.location_x, agobj.location_y + 1)]

    neighbor_locs = [loc for loc in neighbor_locs if 
                     loc[1] >= 1 & loc[1] <= size(sugscape_obj)[1] &
                     loc[2] >= 1 & loc[2] <= size(sugscape_obj)[2]]

    neighbor_agents = [nagobj for nagobj in arr_agents 
                       if nagobj.a.sugar_level > threshold &&
                       (nagobj.location_x, nagobj.location_y) 
                       in neighbor_locs]
    ## println("Entered fetch_neighbors")
    ## readline()
    try
        @assert size(neighbor_agents)[1] in 0:4
    catch except
        println("********************************************************************" )
        println("")
        println("No. of neighbor agents: ", length(neighbor_agents))
        println("********************************************************************" )
    end 
    shuffle!(neighbor_agents)
    return(neighbor_agents)
end ## end fetch_neighbors

function form_possible_institutions!(arr_agents, threshold, sugscape_obj, 
                                     arr_institutions, timeperiod)
    """
    Called during each time period to enact association among agents to form
    institutions. The following conditions exist:
    (a) none of the focal agent's neighbors are part of an institution
    (b) one of the focal agent's neighbors is part of an institution
    (c) two or more of the focal agent's neighbors are part of an institution.
    
    In the case of (a), a randomly chosen neighbor and the focal agent join to
    form an institution.
    In the case of (b), the focal agent joins the institution of the neighbor.
    In the case of (c), the focal agent joins the institution of a 
    randomly-chosen neighbor.
    
    Modifies focal and neighbor agents' institution_id field, if institution-
    formation or joining occurs.
    Reduces the focal agent and neighbor agents' sugarlevels by an amount 
    = sugarlevel - threshold.

    Modifies the array of institutions by adding, as needed, new institution struct.
    
    Returns an array of Institution structures to the caller.
    """
    # println("At the beginning form_possible_institutions! the array of institutions is:")
    # ## println(arr_institutions)
    # println("Length of arr_institutions is: ", string(length(arr_institutions)))
    
    start_inst_id = begin
        if length(arr_institutions) == 1 && arr_institutions[1].institution_id == -1
            1
        else
            maximum([instobj.institution_id for instobj in arr_institutions]) + 1
        end
    end
    ## println("The start_inst_id is: ", string(start_inst_id))
    for agobj in arr_agents
        ## update their institution_ids. 
        ## println("Current object:", string(agobj.a.agent_id), 
                # " its excess sugarlevel: ", agobj.a.sugar_level - threshold)
        ## println("Enter enter")
        ## readline()
        if agobj.a.institution_id == -1 && agobj.a.sugar_level > threshold 
            ## fetch neighbors
            ## println("Checking to see if ", string(agobj.a.agent_id), " can join an", " institution")
            arr_neighbors = fetch_neighbors(agobj, arr_agents, sugscape_obj, threshold)
            ## println("Here here here!")
            if !isempty(arr_neighbors)

                partner_agent = arr_neighbors[1]

                if partner_agent.a.institution_id == -1 ## (a)
                    ## println("No neighbor found with an existing membership, ",
                            # "so creating a new institution with agent: ", 
                            # string(partner_agent.a.agent_id))
                    agobj.a.institution_id = start_inst_id
                    agobj.a.sugar_level = agobj.a.sugar_level - threshold
                    transaction1 = Transaction(agobj.a.sugar_level - threshold,
                                               timeperiod, "deposit", 
                                               agobj.a.agent_id)

                    partner_agent.a.institution_id = start_inst_id
                    partner_agent.a.sugar_level = partner_agent.a.sugar_level - threshold
                    transaction2 = Transaction(partner_agent.a.sugar_level - threshold,
                                               timeperiod, "deposit", 
                                               partner_agent.a.agent_id)

                    insttn_obj = Institution(start_inst_id, 
                                             transaction1.transaction_amount +
                                             transaction2.transaction_amount,
                                             true,
                                             [agobj.a.agent_id, partner_agent.a.agent_id],
                                             [transaction1, transaction2])

                    push!(arr_institutions, insttn_obj)
                    ## println("Added ", string(agobj.a.agent_id), " to a new ",
                            # "institution with id:", start_inst_id)

                else ## (b) and (c)
                    ## println("Inside the (b) and (c) conditional branch")
                    instt_obj = fetch_specific_inst_obj(arr_institutions, 
                                                        partner_agent.a.institution_id)
                    @assert instt_obj.institution_id > 0
                    agobj.a.institution_id = partner_agent.a.institution_id
                    agobj.a.sugar_level = agobj.a.sugar_level - threshold
                    transaction1 = Transaction(agobj.a.sugar_level - threshold,
                                               timeperiod, "deposit", 
                                               agobj.a.agent_id)
                    push!(instt_obj.ledger_transactions, transaction1)
                    
                end ## (a), (b), (c)
            end ## !isempty(arr_neighbors)
        end ## agobj.a.sugar_level > threshold & agobj.a.institution_id > 0
    end ## for agobj    
end ## end of form_possible_institutions

function update_institution_statuses!(arr_institutions, timeperiod)
    """
    Goes through the array of institutions and updates the status
    of each institution: if all of the member agents are dead,
    the institution's alive status is set to false, and a closing transaction
    is added, with a transaction amount of zero, date = current period, and
    agent_id = -999.
    """
    for inst_obj in arr_institutions
        if inst_obj.sugar_level <= 0
            inst_obj.alive = false
            inst_obj.arr_member_ids = []
            push!(inst_obj.ledger_transactions, 
                  Transaction(0, timeperiod, "closure", -999))
        end
    end
end
