include("../src/HermiteConvolutions.jl")

using .HermiteConvolutions

using Plots, PlotThemes

theme(:dark)

n = 5 # size of all vectors and whatnot 

# define functions 
f(x) = sin(x)
g(x) = exp(-x^2)

# approximate in hermite basis 
hermiteF = hermite_fit(f, n)
hermiteG = hermite_fit(g, n)

basisChange = basis_change(n) # get basis change operator from hermite basis to standard basis 

# apply basis change 
standardF = basisChange*hermiteF 
standardG = basisChange*hermiteG

convOperatorF = convolution_operator(standardF, length(standardG)) # get convolution operator for f (at right length for g)

convolvedG = convOperatorF*standardG # convolve g with f 

inverse = left_inverse(convOperatorF) # get left inverse of convolution operator

# see if inverse worked 
display(plot([vec_to_func(inv(basisChange)*inverse*convolvedG), g], 0, 2))

readline()