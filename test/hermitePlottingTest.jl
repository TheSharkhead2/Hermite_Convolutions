include("../src/HermiteConvolutions.jl")

using .HermiteConvolutions

using Plots, PlotThemes

theme(:dark)

h_0(x) = 1 # initial polynomial

# generate the next polynomials 
h_1 = next_hermite_polynomial(h_0)
h_2 = next_hermite_polynomial(h_1)
h_3 = next_hermite_polynomial(h_2)

println(typeof(h_0))

testGeneration = hermite_generation(3, Vector{Function}([h_0]))

display(plot([h_0, h_1, h_2, h_3, testGeneration...], 0, 1))

readline()