
# Confirm/deny the hypothesis that life expectancy (mean and stdev) of isolates
#   and non-isolates are *not* different, across a range of parameter settings.

include("sim.jl")
include("analysis.jl")

λs = [ .1, .5, 1.0, 2.0, 6.0 ]
wn_intensities = [ 1, 10, 25 ]
salaries = [ 10, 40 ]

comparison_results = DataFrame(
    λ=Float16[],
    wn_intensity=Float16[],
    salary=Float16[],
    p_value=Float16[],
)

params[:make_anims] = false
params[:make_sim_plots] = false

for λ in λs
    for wn_intensity in wn_intensities
        for salary in salaries
            prc("Running λ=$(λ), wn=$(wn_intensity), salary=$(salary)...\n")
            results = specnet(λ=λ,
                white_noise_intensity=wn_intensity,
                salary=salary)
            mwut = perform_analysis(results,
                "l$(λ)_wn$(wn_intensity)_s$(salary)")
            if !isnothing(mwut)
                push!(comparison_results,
                    (λ, wn_intensity, salary, pvalue(mwut)))
            else
                push!(comparison_results,
                    (λ, wn_intensity, salary, NaN))
            end
        end
    end
end
