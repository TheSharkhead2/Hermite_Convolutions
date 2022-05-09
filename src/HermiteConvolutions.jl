module HermiteConvolutions

# define addition on functions and scalar multiplication. idea from: https://www.reddit.com/r/Julia/comments/ibcdo6/sum_of_functions/
import Base.+, Base.-, Base.*

+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)
*(a::Float64, f::Function) = (x...) -> a*f(x...)
*(f::Function, g::Function) = (x...) -> f(x...) * g(x...)

using ForwardDiff
using QuadGK
using LinearAlgebra

include("vectorSpace.jl")

export inner_product, norm

include("hermitePolynomials.jl")

export next_hermite_polynomial, hermite_generation, normalized_hermite_generation, vec_to_func, next_hermite, hermite_generation, norm_hermite_generation, basis_change, hermite_fit

include("convolution.jl")

export generic_convolution_operator, convolution_operator, left_inverse

end # module HermiteConvolutions
