

function generateCost(n::Int64)
    c = zeros(n,n)
    sumC = 0
    for i in 1:n
        for j in 1:n
            c[i,j] = rand([i for i in 2:10])
            sumC += c[i,j]
        end
    end
    Amin = div(60*n*n, 100)
    Amax = div(70*n*n, 100)
    B = div(70*sumC, 100)*10
    lambda = 20
    return c, Amin, Amax, B, lambda
end