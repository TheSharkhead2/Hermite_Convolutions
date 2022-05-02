function inner_product(f::Function, g::Function)
    integrand = f*g*(x->exp((-x^2)/2))

    integral, err = quadgk(integrand, -Inf, Inf)

    integral 
end # function inner_product

function norm(f::Function)
    sqrt(inner_product(f, f))
end # function norm