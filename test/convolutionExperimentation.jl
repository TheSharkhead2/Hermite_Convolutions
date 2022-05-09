include("../src/HermiteConvolutions.jl")

using .HermiteConvolutions

using Plots, PlotThemes, QuadGK

# testConv(x) = generic_convolution_operator(x->exp(-x^2), x->sin(x))(x)
# # testConv(x) = (quadgk(t->exp(-t^2)*sin(x-t), -Inf, Inf))[1]

# println(typeof(testConv))

# theme(:dark)

# display(plot([testConv, x->exp(-x^2), x->sin(x)], 0, 2))

# readline()
