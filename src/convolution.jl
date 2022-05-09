"""
Performs generic convolution (without use of Hermite basis)
on normal functions

"""
function generic_convolution_operator(f::Function, g::Function)
    x -> (quadgk(t -> f(t)*g(x-t), -Inf, Inf))[1]
end # function generic_convolution_operator

"""
Gets Toeplitz matrix representing convolution operator defined with f 
to act on a vector of length dim

"""
function convolution_operator(f::Vector{Float64}, dim::Int)
    toeplitzMatrix = zeros(length(f)+dim-1, dim) # create empty matrix of right size

    for i ∈ 1:dim # for each column 
        column = zeros(i-1)
        if size(toeplitzMatrix)[1] - (i-1 + length(f)) >= 0 # to prevent negative vector dimensions, if the column is over dimension, don't add vector of negative size  
            append!(column, f, zeros(size(toeplitzMatrix)[1] - (i-1 + length(f)))) # add vector of zeros to the end of the column
        else 
            append!(column, f)
            column = column[1:size(toeplitzMatrix)[1]] # truncate column to fit toeplitz matrix
        end # if

        toeplitzMatrix[1:size(toeplitzMatrix)[1], i] = column # set column of toeplitz matrix
    end # for 

    toeplitzMatrix # return toeplitz matrix
end # function convolution

"""
Gets "left inverse" of mxn matrix A where m>n and dim range(A) = n

"""
function left_inverse(A::Matrix{Float64})
    inv(transpose(A)*A)*transpose(A) # return left inverse
end # function left_inverse