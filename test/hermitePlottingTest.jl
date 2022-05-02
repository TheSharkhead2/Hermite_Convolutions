include("../src/HermiteConvolutions.jl")

using .HermiteConvolutions

using Plots, PlotThemes

theme(:dark)

h_0(x) = 1/sqrt(sqrt(2*π)) # initial polynomial

# generate the next polynomials 
polynomials = hermite_generation(3, Vector{Function}([h_0]))

display(plot(polynomials, 0, 2))

for polynomial in polynomials 
    println(norm(polynomial))
end # for

println(inner_product(polynomials[1], polynomials[2]))
println(inner_product(polynomials[1], polynomials[3]))
println(inner_product(polynomials[2], polynomials[3]))
println(inner_product(polynomials[1], polynomials[4]))

readline()