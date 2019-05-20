#!/usr/bin/env julia
using Revise
using RCall 
@rlibrary ineq

include("sim.jl")
include("setup_params.jl")

results = specnet()
