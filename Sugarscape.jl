using Random
using Plots
# using VegaLite
using DataFrames
using Statistics
## using Match

mutable struct Sugarcell
    location_x::Int64
    location_y::Int64
    growth_rate::Float64
    carrying_capacity::Float64    
    sugar_level::Float64
    occupied::Bool
end

function generate_sugarscape(side, growth_rate, carrying_cap, scenario=1)
    """
    Create an array of cells, representing, together the sugarscape object.
    Each cell's sugar grows back at rate given by growth_rate.
    The max. value of sugar in a cell is specified by its carrying capacity -
    constant across all cells.
    Each cell's initial sugarlevel is initialized to 
    0) a random value, 
    1) a central region with maximal sugar value, and sugarlevel falling off
       uniformly as an inverse power of 0.3 with distance.
    2) two regions of maximal sugar value, located in the north-east and south-
       west quadrants
    3) four regions of maximal sugar value, located roughly in the middle of 
       each of the four quadrants
    """
    array_suglevels = zeros(side, side)
    if scenario == 0
        array_suglevels = [rand()* carrying_cap for x in 1:side, y in 1:side]
        # return([Sugarcell(x, y, growth_rate, carrying_cap, ceil(rand()*5), false)
        #         for x in 1:side, y in 1:side])
    elseif scenario == 1
        central_x, central_y = side/2, side/2
        array_suglevels = [carrying_cap/(sqrt((x - central_x)^2 + (y - central_y)^2)^0.3)
                           for x in 1:side, y in 1:side]
        # return([Sugarcell(x, y, growth_rate, carrying_cap, array_suglevels[x,y], false)
        #         for x in 1:side, y in 1:side])
    elseif scenario == 2
        ne_row = side/3
        ne_col = side * 2/3
        sw_row = side * 2/3
        sw_col = side/3
        for x in 1:side
            for y in 1:side
                dist_sw = sqrt((x - sw_row)^2 + (y - sw_col)^2)
                dist_ne = sqrt((x - ne_row)^2 + (y - ne_col)^2)
                array_suglevels[x, y] = carrying_cap/(min(dist_ne, dist_sw)^0.3)
                end
            end
        # return([Sugarcell(x, y, growth_rate, carrying_cap, array_suglevels[x,y], false)
        #         for x in 1:side, y in 1:side])
    elseif scenario == 3
        ne_row = side/4
        ne_col = side * 3/4
        
        nw_row = side/4
        nw_col = side/4
        
        se_row = side * 3/4
        se_col = side * 3/4

        sw_row = side * 3/4
        sw_col = side/4
        
        array_suglevels = zeros(side, side)
        for x in 1:side
            for y in 1:side
                dist_sw = sqrt((x - sw_row)^2 + (y - sw_col)^2)
                dist_ne = sqrt((x - ne_row)^2 + (y - ne_col)^2)
                dist_se = sqrt((x - se_row)^2 + (y - se_col)^2)
                dist_nw = sqrt((x - nw_row)^2 + (y - nw_col)^2)
                
                array_suglevels[x, y] = carrying_cap/(min(dist_ne, dist_nw,
                                                         dist_se, dist_sw)^0.3)
            end
        end ## end of outer for
    end ## end of scenarios
    return([Sugarcell(x, y, growth_rate, carrying_cap, array_suglevels[x,y], false)
            for x in 1:side, y in 1:side])
end ## end of function generate_sugarscape

# function get_sugscape_suglevels(sugscape_obj)
#     return([sugarobj.sugar_level for sugarobj in sugscape_obj])
# end ## get_sugscape_suglevels()

function get_sugarscape_stats(sugscape_obj)
    """
    Returns min., max., mean, median, and sd. of sugar levels.
    """
    # arr_suglevels = get_sugscape_suglevels(sugscape_obj)
    arr_suglevels = [sugarobj.sugar_level for sugarobj in sugscape_obj]
    return(minimum(arr_suglevels), maximum(arr_suglevels), mean(arr_suglevels),
           median(arr_suglevels), std(arr_suglevels))
end ## get_sugarscape_stats()


function regenerate_sugar!(sugscape_obj)
    """
    Modifies the cell objects in place by regenerating sugar in each cell by an
    amount specified by growth_date, with an upper bound of carrying_capacity.
    Does not return a new object, since the cell objects are modified in place.
    """
    print("Inside regenerate_sugar!()\n\n")
    for cellobj in sugscape_obj
        old_val = cellobj.sugar_level
        new_val = 
        cellobj.sugar_level = max(0, min(cellobj.sugar_level*(1 + cellobj.growth_rate),
                                         cellobj.carrying_capacity)) 
        # println(string("Changing sugar level from ", old_val, " to ",
        #                cellobj.sugar_level, "\n"))
    end
    return
end ## regenerate_sugar!()

function plot_sugar_concentrations!(sugscape_obj)
    """
    Generates a heatmap of the sugarlevels of all the cells
    """
    print("Inside plot_sugar_concentrations\n")
    gr()
    nrows, ncols = size(sugscape_obj)
    xs = [string(i) for i = 1:nrows]
    ys = [string(i) for i = 1:ncols]
    z = [sugscape_obj[row, col].sugar_level for row in 1:nrows, col in 1:ncols] 
    heatmap(xs, ys, z, aspect_ratio=1)
end ## plot_sugar_concentrations()
