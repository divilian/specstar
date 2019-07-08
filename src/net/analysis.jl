
using HypothesisTests

# Assumption: sim.jl has been executed, and its return value saved in a
#   variable called "r".

# Perform Wilcoxon rank-sum test (a.k.a. "Mann-Whitney U Test") to confirm or
# deny the hypothesis that isolates die sooner/later than non-isolates on
# average.
lh = r[:life_history]
lifespans = by(lh, [:agent, :original_isolate], lifespan = :iter => maximum)

lifespanp = plot(lifespans, x=:original_isolate, y=:lifespan, Geom.boxplot,
    color=:original_isolate,
    Scale.y_continuous(minvalue=0))
draw(PNG("$(tempdir())/lifespanComparison.png"), lifespanp)

ilifespans = lifespans[lifespans[:original_isolate].==true,:lifespan]
nlifespans = lifespans[lifespans[:original_isolate].==false,:lifespan]
srt = MannWhitneyUTest(ilifespans, nlifespans)
