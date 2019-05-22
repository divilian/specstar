
## Input parameters
params = Dict{Symbol,Any}(
    :N => 30,                       # number of agents
    :num_iters => 50,               # number of iterations the simulation runs
    :openness => 0.0,               # 0 <=> openness <=> 1
                                    #   (0: always choose from neighbor,
                                    #    1: always choose from entire city)
    :init_sg_lvl => 100,            # each agent starts with wealth
                                    #   U~(1, init_sg_lvl)
    :salary_range => 10,            # each iteration, each agent receives/loses
                                    #   U~(-salary_range, salary_range) wealth
    :proto_threshold => 50,         # each agent in an encounter must have
                                    #   > wealth than this to form a proto
    :make_anims => false,           # create animations of results?
    :animation_delay => 20,         # milliseconds between animation frames
    :random_seed => 1234,           # random number generator starting seed

    :whichGraph=>"erdos_renyi",     # the name of graph to generate e.g.
                                    #   - "erdos_renyi"
                                    #   - "scale_free"
                                    #   - "small_world"

    :ER_prob=>0.2,                  # probability of edge between two vertices

    :SF_edges=>40,                  # number of edges in SF model
    :SF_degree=>2,                  # exponent of expected power law degree
                                    #   distribution

    :SW_degree=>4,                  # degree (even integer) number of neighbors
    :SW_prob=>0.2,                  # probability of rewiring (between 0/1)
)
