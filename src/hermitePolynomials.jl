"""
Generate the k+1 Hermite polynomial

"""
function next_hermite_polynomial(h_k::Function)
    (x -> x)*h_k - (x -> ForwardDiff.derivative(h_k, x)) # hₖ₊₁ = x*hₖ - hₖ'
end # function hermite_generation

"""
Generate n Hermite polynomials (up to k=n) starting with h₀ 
Input initial h₀ as only item in polynomials vector

"""
function hermite_generation(n::Int, polynomials::Vector{Function})
    if n == 0 # if down to just first element, then return the vector 
        return polynomials 
    else 
        push!(polynomials, next_hermite_polynomial(polynomials[length(polynomials)])) # add the next polynomial to the vector
        return hermite_generation(n-1, polynomials) # recurse
    end # if 
end # function hermite_generation

"""
Generate n Hermite polynomials that are normalized (using hermite_generation method )

"""
function normalized_hermite_generation(n::Int; h_0::Function=(x->1))
    nonscaled = hermite_generation(n, Vector{Function}([h_0]))

    [1/norm(polynomial) * polynomial for polynomial in nonscaled] # normalize each polynomial
end # function normalized_hermite_generation

"""
Fit a set of n Hermite polynomials to a function f 

"""
function hermite_fit(f::Function, n::Int; h_0::Function=(x->1))
    polynomials = normalized_hermite_generation(n; h_0=h_0) # generate n hermite polynomials 
    [inner_product(polynomial, f) for polynomial in polynomials] # evaluate the inner product of each polynomial with f
end # function hermite_fit

"""
Takes vector representing hermite polynomial and converts it to a Julia function

"""
function vec_to_func(h::Vector{Float64})
    f(x) = 0 # blank function to start with 

    xterms = [x->x^i for i ∈ 0:length(h)-1] # create xⁱ for each term in vector (essentially standard basis vectors)

    for (a, x) ∈ zip(h, xterms) # for each term in vector
        f = f + a * x # add term to function
    end # for 

    f # return function
end # function vec_to_func

"""
Using vectors instead of functions to represent Hermite polynomials, 
generate the next polynomial based on the last one

"""
function next_hermite(h_k::Vector{Float64})
    rightshift = [0] # need to account for multiplying by x, so shift vector to right. This 0 is that shifting factor
    append!(rightshift, h_k) # append the last polynomial to the rightshift vector to shift it 

    leftshift = copy(h_k) # copy vector before left shifting (derivative)
    popat!(leftshift, 1) # remove the first element from the vector, completing the left shift (remove the x⁰ term)
    multipliers = [x for x ∈ 1:length(leftshift)] # create a vector of multipliers for the left shift (essentially the exponent that is taken down with power rule) 
    leftshift .*= multipliers # multiply the leftshift vector by the multipliers vector 

    append!(leftshift, [0, 0]) # append the 0s to the leftshift vector to match length of rightshift (1 for leftshift loss and 1 for right shift gain)

    h_k1 = rightshift - leftshift # hₖ₊₁ = hₖ - hₖ'
end # function next_hermite

"""
Using the vector representation of the polynomials (standard basis), generate n 
hermite polynomials 

"""
function hermite_generation(n::Int, polynomials::Vector{Vector{Float64}})
    if n <= 0 # when the bottom of the recursion is reached, return the vector
        return polynomials 
    else 
        push!(polynomials, next_hermite(polynomials[length(polynomials)])) # add the next polynomial to the vector
        return hermite_generation(n-1, polynomials) # recurse
    end # if
end # function hermite_generation

"""
Using the vector representation of the hermite polynomials (standard basis), generate 
n normalized hermite polynomials

"""
function norm_hermite_generation(n::Int; h_0::Vector{Float64}=[1.0])
    nonscaled = hermite_generation(n, [h_0]) # generate n hermite polynomials

    [1/norm(vec_to_func(polynomial)) * polynomial for polynomial ∈ nonscaled] # normalize each polynomial (need to convert to function first to calculate norm)
end # function norm_hermite_polynomials