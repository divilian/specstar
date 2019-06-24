
## Input parameters
params = Dict{Symbol,Any}(
    :N => 30,                       # number of agents
    :max_iters => 50,               # maximum number of iterations the simulation
                                    #   runs before termination (will terminate
                                    #   earlier if stopping condition reached)
    :starvation_period => 10,       # the number of iterations to starve all agents
                                    #   after the simulation has reached stage 3
                                    #   (all live non-isolate agents in protos)
    :openness => float(0.0),        # 0 <=> openness <=> 1
                                    #   (0: always choose from neighbor,
                                    #    1: always choose from entire city)
    :init_sg_lvl => 100,            # each agent starts with wealth
                                    #   ~U(1, init_sg_lvl)
    :metabolic_rate => 5,           # each iteration, each agent consumes exactly
                                    #   this much sugar
    :salary => 10,                  # each iteration, each agent receives exactly
                                    #   this much sugar
    :white_noise_intensity =>       # each iteration, an agent receives (or loses)
        float(1.0),                 #   ~N(0, white_noise_intensity) extra sugar
    :proto_threshold => 50,         # each agent in an encounter must have
                                    #   > wealth than this to form a proto
    :make_anims => false,           # create animations of results?
    :make_sim_plots => true,        # create plots (for individual simulations)?
    :num_boot_samples => 1000,      # for single sims, the number of bootstrap
                                    #   samples used in computing CI for single
                                    #   Gini
                                    # for param sweeps, the number of bootstrap
                                    #   samples used in computing CI for all
                                    #   Ginis in runs with same parameters
    :animation_delay => 20,         # milliseconds between animation frames
    :random_seed => 1234,           # random number generator starting seed

    :whichGraph => "erdos_renyi",   # the name of graph to generate e.g.
                                    #   - "erdos_renyi"
                                    #   - "scale_free"
                                    #   - "small_world"
                                    #   - "complete"
                                    #   - "empty"

    :λ => float(1),                 # ER: expected number if edges per node

    :SF_edges => 40,                # SF: number of edges
    :SF_degree => 2,                # SF: exponent of expected power law degree
                                    #   distribution

    :β => float(0.2),               # SW: prob. of rewiring
    :k => 2,                        # SW: degree of nodes
)
