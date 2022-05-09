include("../src/HermiteConvolutions.jl")

using .HermiteConvolutions

using Plots, PlotThemes, QuadGK

theme(:dark)

# generate the next polynomials 
polynomials = normalized_hermite_generation(4)

functions = []

for i ∈ polynomials
    for j ∈ polynomials
        push!(functions, t->quadgk(s->i(s)*j(t-s), -Inf, Inf)[1])
    end # for
end # for

println(functions[1](3))

display(plot(functions, 0, 2))

readline()