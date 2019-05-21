#!/usr/bin/env julia
using Revise

include("sim.jl")
include("setup_params.jl")

results = specnet()
