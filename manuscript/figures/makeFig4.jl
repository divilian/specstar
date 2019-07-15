
# Assumptions: plot_df_wn5.csv and plot_df_wn20.csv are in the current
# directory, giving the output of a parameter sweep of Gini vs. λ with white
# noise of 5, and of 20, respectively.

using DataFrames
using CSV
using Cairo, Fontconfig
using Gadfly

pdf5 = CSV.read("plot_df_wn5.csv")
pdf20 = CSV.read("plot_df_wn20.csv")
pdf5[:wn] = "5"
pdf20[:wn] = "20"
pdf = [pdf5;pdf20]

l = layer(pdf5,
    x=:λ, y=:gini, Geom.line,
    Theme(default_color="blue",
        line_style=[:solid],
        line_width=.5mm))
l2 = layer(pdf20,
    x=:λ, y=:gini, Geom.line,
    Theme(default_color="red",
        line_style=[:solid],
        line_width=.5mm))
lhi = layer(pdf5,
    x=:λ, y=:gini_highCI, Geom.line,
    Theme(default_color=colorant"lightgray"))
llo = layer(pdf5,
    x=:λ, y=:gini_lowCI, Geom.line,
    Theme(default_color=colorant"lightgray"))
lhi2 = layer(pdf20,
    x=:λ, y=:gini_highCI, Geom.line,
    Theme(default_color=colorant"lightgray"))
llo2 = layer(pdf20,
    x=:λ, y=:gini_lowCI, Geom.line,
    Theme(default_color=colorant"lightgray"))
lfill = layer(pdf5,
    x=:λ, ymin="gini_lowCI", ymax="gini_highCI",
    Theme(default_color=colorant"lightblue"),
    Geom.ribbon)
lfill2 = layer(pdf20,
    x=:λ, ymin="gini_lowCI", ymax="gini_highCI",
    Theme(default_color=colorant"pink"),
    Geom.ribbon)
layers = Layer[]
append!(layers, l)
append!(layers, l2)
append!(layers, lhi)
append!(layers, llo)
append!(layers, lhi2)
append!(layers, llo2)
append!(layers, lfill)
append!(layers, lfill2)

p = plot(pdf, layers,
    Scale.color_discrete_manual("blue","red"),
    Guide.xlabel("λ"),
    Guide.ylabel("Gini (of effective wealth)"),
    Guide.manual_color_key("White noise (σ²)", ["5     .","20"], ["blue","red"]),
    style(background_color=colorant"white",key_position=:bottom))

draw(PNG("giniVsLambda.png"), p)
