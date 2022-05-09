include("../src/HermiteConvolutions.jl")

using .HermiteConvolutions

using Plots, PlotThemes

theme(:dark)

# generate the next polynomials 
polynomials = normalized_hermite_generation(5)
polynomialsVec = norm_hermite_generation(5)
polyVecFunc = vec_to_func.(polynomialsVec)

println(hermite_generation(5, [[1.0]]))

display(plot([polynomials..., polyVecFunc...], 0, 2))

# for polynomial in polynomials 
#     println(norm(polynomial))
# end # for

# println(inner_product(polynomials[1], polynomials[2]))
# println(inner_product(polynomials[1], polynomials[3]))
# println(inner_product(polynomials[2], polynomials[3]))
# println(inner_product(polynomials[1], polynomials[4]))

readline()