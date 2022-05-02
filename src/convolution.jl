"""
Performs generic convolution (without use of Hermite basis)
on normal functions

"""
function generic_convolution_operator(f::Function, g::Function)
    x -> (quadgk(t -> f(t)*g(x-t), -Inf, Inf))[1]
end # function generic_convolution_operator
