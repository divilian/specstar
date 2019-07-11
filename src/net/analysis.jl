
using HypothesisTests

# This function returns either the result of a Mann-Whitney test comparing
#   lifespan of isolates with non-isolates, or "nothing" if no such comparison
#   could be made (for example, if there weren't sufficient numbers of both
#   isolates and non-isolates.) As a side effect, produces a number of plots in
#   /tmp comparing the two.
# The "results" argument should be the return value from specnet().
#   "id_string" is an arbitrary string used in plot filenames and titles.
function perform_analysis(results, id_string="")

    lh = results[:life_history]
    lifespans = by(lh, [:agent, :stage3_isolate], lifespan = :iter => maximum)

    if nrow(lifespans[lifespans[:stage3_isolate],:]) ≤ 1  ||
        nrow(lifespans[.!lifespans[:stage3_isolate],:]) ≤ 1
        prc("WARNING: could not perform analysis on $(id_string).\n")
        prc("(There were $(sum(lifespans[:stage3_isolate])) isolates and " *
            "$(sum(.!lifespans[:stage3_isolate])) non-isolates.)\n")
        return
    end

    lifespanp = plot(lifespans, x=:stage3_isolate, y=:lifespan, Geom.boxplot,
        color=:stage3_isolate,
        Scale.color_discrete_manual("navy","orange",
            levels=[false, true]),
        style(background_color=colorant"white",key_position=:none),
        Guide.title(id_string),
        Scale.y_continuous(minvalue=0))
    draw(PNG("$(tempdir())/lifespanComparison_$(id_string).png"), lifespanp)

    ihistogramp = plot(lifespans[lifespans[:stage3_isolate],:],
        x=:lifespan, Geom.histogram(bincount=30, density=true), Coord.cartesian(xmax=maximum(lifespans[:lifespan])+10),
        style(background_color=colorant"white", default_color=colorant"orange"),
        Guide.xlabel("Lifespan of isolates"),
        Guide.title(id_string),
        Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))

    nhistogramp = plot(lifespans[.!lifespans[:stage3_isolate],:],
        x=:lifespan, Geom.histogram(bincount=30, density=true), Coord.cartesian(xmax=maximum(lifespans[:lifespan])+10),
        style(background_color=colorant"white", default_color=colorant"navy"),
        Guide.xlabel("Lifespan of non-isolates"),
        Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))

    draw(PNG("$(tempdir())/lifespanHist_$(id_string).png", 4inch, 6inch), vstack(ihistogramp, nhistogramp))

    densityp = plot(lifespans,
        color=:stage3_isolate, x=:lifespan,
        Geom.density(bandwidth=3),
        Scale.color_discrete_manual("navy","orange",
            levels=[false, true]),
        style(background_color=colorant"white"),
        Guide.xlabel("KDE of lifespan"),
        Guide.title(id_string),
        Scale.x_continuous(minvalue=0, maxvalue=maximum(lifespans[:lifespan])))
    draw(PNG("$(tempdir())/lifespanKDE_$(id_string).png"), densityp)

    # Perform Wilcoxon rank-sum test (a.k.a. "Mann-Whitney U Test") to confirm or
    # deny the hypothesis that isolates die sooner/later than non-isolates on
    # average.
    ilifespans = lifespans[lifespans[:stage3_isolate].==true,:lifespan]
    nlifespans = lifespans[lifespans[:stage3_isolate].==false,:lifespan]
    if length(ilifespans) > 0 && length(nlifespans) > 0
        mwut = MannWhitneyUTest(ilifespans, nlifespans)
        return mwut
    end
end
