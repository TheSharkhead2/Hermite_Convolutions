module HermiteConvolutions

# define addition on functions and scalar multiplication. idea from: https://www.reddit.com/r/Julia/comments/ibcdo6/sum_of_functions/
import Base.+, Base.-, Base.*

+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)
*(a::ComplexF64, f::Function) = (x...) -> a*f(x...)
*(f::Function, g::Function) = (x...) -> f(x...) * g(x...)

using ForwardDiff
using QuadGK

include("vectorSpace.jl")

export inner_product

include("hermitePolynomials.jl")

export next_hermite_polynomial, hermite_generation

end # module HermiteConvolutions
