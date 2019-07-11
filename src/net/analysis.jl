
using HypothesisTests

# Assumption: sim.jl has been executed, and its return value saved in a
#   variable called "r".                                                                               

lh = r[:life_history]
lifespans = by(lh, [:agent, :stage3_isolate], lifespan = :iter => maximum)

if nrow(lifespans[lifespans[:stage3_isolate],:]) ≤ 1  ||
    nrow(lifespans[.!lifespans[:stage3_isolate],:]) ≤ 1
    prc("WARNING: could not perform analysis on $(filename_components).\n")
    prc("(There were $(sum(lifespans[:stage3_isolate])) isolates and " *
        "$(sum(.!lifespans[:stage3_isolate])) non-isolates.)\n")
    error(1)
end

lifespanp = plot(lifespans, x=:stage3_isolate, y=:lifespan, Geom.boxplot,
    color=:stage3_isolate,
    Scale.color_discrete_manual("navy","orange",
        levels=[false, true]),
    style(background_color=colorant"white",key_position=:none),
    Scale.y_continuous(minvalue=0))
draw(PNG("$(tempdir())/lifespanComparison.png"), lifespanp)

ihistogramp = plot(lifespans[lifespans[:stage3_isolate],:],
    x=:lifespan, Geom.histogram(bincount=30, density=true), Coord.cartesian(xmax=maximum(lifespans[:lifespan])+10),
    style(background_color=colorant"white", default_color=colorant"orange"),
    Guide.xlabel("Lifespan of isolates"),
    Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))

nhistogramp = plot(lifespans[.!lifespans[:stage3_isolate],:],
    x=:lifespan, Geom.histogram(bincount=30, density=true), Coord.cartesian(xmax=maximum(lifespans[:lifespan])+10),
    style(background_color=colorant"white", default_color=colorant"navy"),
    Guide.xlabel("Lifespan of non-isolates"),
    Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))

draw(PNG("$(tempdir())/lifespanHist.png", 4inch, 6inch), vstack(ihistogramp, nhistogramp))

densityp = plot(lifespans,
    color=:stage3_isolate, x=:lifespan,
    Geom.density(bandwidth=3),
    Scale.color_discrete_manual("navy","orange",
        levels=[false, true]),
    style(background_color=colorant"white"),
    Guide.xlabel("KDE of lifespan"),
    Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))
draw(PNG("$(tempdir())/lifespanKDE.png"), densityp)

# Perform Wilcoxon rank-sum test (a.k.a. "Mann-Whitney U Test") to confirm or
# deny the hypothesis that isolates die sooner/later than non-isolates on
# average.
ilifespans = lifespans[lifespans[:stage3_isolate].==true,:lifespan]
nlifespans = lifespans[lifespans[:stage3_isolate].==false,:lifespan]
if length(ilifespans) > 0 && length(nlifespans) > 0
    srt = MannWhitneyUTest(ilifespans, nlifespans)
end
