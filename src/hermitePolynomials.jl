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