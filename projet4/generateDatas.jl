

function buildtij(n::Int64, m::Int64)
    t = zeros(n+2, m+2)
    r = [i for i in 60:99]
    for i in 2:n+1
        for j in 2:m+1
            t[i,j] = rand(r)
        end
    end
    return t
end

function instance2()
    t = [ 
        0 0 0 0 0 0 0 ;
        0 10 10 10 1 10 0 ;
        0 10 10 1 1 10 0 ;
        0 10 10 1 10 10 0 ;
        0 1 10 10 10 10 0 ;
        0 1 10 10 10 10 0 ;
        0 0 0 0 0 0 0
    ]
    return t
end