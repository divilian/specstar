using Random
using RCall

mutable struct Agent
    agent_id::Int64
    location_x::Int64
    location_y::Int64
    vision::Int64
    metabolic_rate::Int64
    sugar_level::Float64
    alive::Bool
    institution_id::Int64
end

function fetch_best_location(agent_obj, sugscape_obj)
    """
    Returns a tuple representing the location of the sugarscape cell that is
    not occupied and has the highest sugarlevel in the Von-Neumann neighborhood.
    Returns null if no such cell can be found.
    """
    low_row = max(1, agent_obj.location_x - agent_obj.vision)
    low_col = max(1, agent_obj.location_y - agent_obj.vision)

    hi_row = min(size(sugscape_obj, 1), agent_obj.location_x + agent_obj.vision)
    hi_col = min(size(sugscape_obj, 2), agent_obj.location_y + agent_obj.vision)

    # try
    #     @assert all(map(x -> x > 1, [low_row, low_col, hi_row, hi_col]))
    # catch
    #     println("--------------------------")
    #     println("low_row:", string(low_row), " hi_row:", string(hi_row),
    #           " low_col:", string(low_col), " hi_col:", string(hi_col))
    #     println("---------------------------")
    # end
    poss_cells = [sugscape_obj[x, y] for x in low_row:hi_row,
                  y in low_col:hi_col
                  if (!sugscape_obj[x, y].occupied &
                      (sugscape_obj[x, y].sugar_level > 0))]
    # println("Here are the potential cells")
    # x = readline()
    # println([(cellobj.location_x, cellobj.location_y, cellobj.sugar_level)
    #          for cellobj in poss_cells])
    if length(poss_cells) > 0
        a_suglevels = [cellobj.sugar_level for cellobj in poss_cells]
        a_cells =[cellobj for cellobj in poss_cells
                  if cellobj.sugar_level == maximum(a_suglevels)]
        return((a_cells[1].location_x, a_cells[1].location_y))
    else
        return(nothing)
    end    
end ## end fetch_best_location

function locate_move_feed!(agent_obj, sugscape_obj)
    """
    For a given agent, performs the feeding operation first. If sugar is not 
    available within the agent or in conjunction with current location, moves
    to a new location, if available. If not available, dies. 
    
    TODO: in a future version, will have a clock to keep of number of time
    periods that have elapsed in starvation mode, and performs death when
    a threshold has passed.
    """
    # println("Performing locate-move-feed on agent:", string(agent_obj.agent_id), "")
    
    if(agent_obj.alive)
        if agent_obj.sugar_level >= agent_obj.metabolic_rate
            agent_obj.sugar_level = agent_obj.sugar_level - agent_obj.metabolic_rate
            # println("Agent ", string(agent_obj.agent_id), " just drew from its self ",
            #       "sugar reserve!")
            ## x = readline() 
        elseif sugscape_obj[agent_obj.location_x,
                            agent_obj.location_y].sugar_level +
                                agent_obj.sugar_level >= agent_obj.metabolic_rate

            agent_obj.sugar_level = sugscape_obj[agent_obj.location_x,
                                                 agent_obj.location_y].sugar_level +
                                                     agent_obj.sugar_level - agent_obj.metabolic_rate
            # println("Agent ", string(agent_obj.agent_id), " loaded up at its current location!")
            sugscape_obj[agent_obj.location_x, agent_obj.location_y].sugar_level -=
                sugscape_obj[agent_obj.location_x, agent_obj.location_y].sugar_level +
                agent_obj.sugar_level - agent_obj.metabolic_rate
            ## x = readline()
        else ## need to move
            ## identify best location
            new_location = fetch_best_location(agent_obj, sugscape_obj)
            if isnothing(new_location)
                ## no food available at current location and no new source
                ## of food available, so set alive status to false
                agent_obj.alive = false
                agent_obj.location_x, agent_obj.location_y = -1, -1
                # println("Agent ", string(agent_obj.agent_id), " starved to death!")
                ## x = readline()
            else
                ## move to and load from the new cell location
                sugscape_obj[agent_obj.location_x,
                             agent_obj.location_y].occupied = false
                # println("Agent ", string(agent_obj.agent_id, " moving from ",
                #                        string(agent_obj.location_x), ",",
                #                        string(agent_obj.location_y), ",", " to "))
                
                agent_obj.location_x, agent_obj.location_y = new_location[1], new_location[2]
                
                sugscape_obj[agent_obj.location_x, agent_obj.location_y].sugar_level -=
                    sugscape_obj[agent_obj.location_x, agent_obj.location_y].sugar_level +
                    agent_obj.sugar_level - agent_obj.metabolic_rate
                
                sugscape_obj[agent_obj.location_x,
                             agent_obj.location_y].occupied = true
                # println(string(agent_obj.location_x), ",",
                #         string(agent_obj.location_y), "")
                ## x = readline()
            end
            ## consume 
        end 
        ## move
    else ## agent is dead
        # println("Tried to animate a dead agent: ", string(agent_obj.agent_id))
        ## x = readline()
    end 
end ## locate_move_feed!()

function life_check!(arr_agents)
    for ag_obj in arr_agents
        if ag_obj.sugar_level < 0
            ag_obj.alive = false
            # println("Agent ", string(ag_obj), " has died!")
            ## x = readline()
        end
    end
end

function compute_Gini(arr_agents)
    arr_suglevels = [agobj.sugar_level for agobj in
                     arr_agents]
    R"library(ineq)"
    gini = R"ineq($arr_suglevels, type='Gini')"
    # println(gini)
    return(gini)    
end
