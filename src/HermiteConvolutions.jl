module HermiteConvolutions

# define addition on functions and scalar multiplication. idea from: https://www.reddit.com/r/Julia/comments/ibcdo6/sum_of_functions/
import Base.+, Base.-, Base.*

+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)
*(a::Float64, f::Function) = (x...) -> a*f(x...)
*(f::Function, g::Function) = (x...) -> f(x...) * g(x...)

using ForwardDiff
using QuadGK

include("vectorSpace.jl")

export inner_product, norm

include("hermitePolynomials.jl")

export next_hermite_polynomial, hermite_generation, normalized_hermite_generation

include("convolution.jl")

export generic_convolution_operator

end # module HermiteConvolutions
