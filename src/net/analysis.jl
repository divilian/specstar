
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
    Scale.color_discrete_manual("navy","orange",
        levels=[false, true]),
    style(background_color=colorant"white",key_position=:none),
    Scale.y_continuous(minvalue=0))
draw(PNG("$(tempdir())/lifespanComparison.png"), lifespanp)

ihistogramp = plot(lifespans[lifespans[:original_isolate],:],
    x=:lifespan, Geom.histogram(bincount=30, density=true), Coord.cartesian(xmax=maximum(lifespans[:lifespan])+10),
    style(background_color=colorant"white", default_color=colorant"orange"),
    Guide.xlabel("Lifespan of isolates"),
    Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))

nhistogramp = plot(lifespans[.!lifespans[:original_isolate],:],
    x=:lifespan, Geom.histogram(bincount=30, density=true), Coord.cartesian(xmax=maximum(lifespans[:lifespan])+10),
    style(background_color=colorant"white", default_color=colorant"navy"),
    Guide.xlabel("Lifespan of non-isolates"),
    Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))

draw(PNG("$(tempdir())/lifespanHist.png", 4inch, 6inch), vstack(ihistogramp, nhistogramp))

densityp = plot(lifespans,
    color=:original_isolate, x=:lifespan,
    Geom.density(bandwidth=3),
    Scale.color_discrete_manual("navy","orange",
        levels=[false, true]),
    style(background_color=colorant"white"),
    Guide.xlabel("KDE of lifespan"),
    Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))
draw(PNG("$(tempdir())/lifespanKDE.png"), densityp)

ilifespans = lifespans[lifespans[:original_isolate].==true,:lifespan]
nlifespans = lifespans[lifespans[:original_isolate].==false,:lifespan]
srt = MannWhitneyUTest(ilifespans, nlifespans)
