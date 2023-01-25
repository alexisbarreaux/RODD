using Random

function generateDatas(n::Int64, k::Int64, e::Int64)
    p = zeros(k,n+2,n+2)
    for l in 1:k
        for i in 2:n+1
            for j in 2:n+1
                cpt = rand()
                if l <= e && cpt <= 0.1
                    p[l,i,j] = round(rand()*0.5, digits=1)
                    while p[l,i,j]==1.0
                        p[l,i,j] = round(rand()*0.5, digits=1)
                    end
                elseif l > e && cpt <= 0.2
                    p[l,i,j] = round(rand()*0.5, digits=2)
                    while p[l,i,j]==1.0
                        p[l,i,j] = round(rand()*0.5, digits=2)
                    end
                end
            end
        end
    end
    return p
end


function test(p)
    for l in 1:size(p,1)
        for i in 2:size(p,2)-1
            for j in 2:size(p,3)-1
                if p[l,i,j]!=0
                    println("p[",l,",",i,",",j,"]=",p[l,i,j])
                end
            end
        end
    end
end


function generateCost(n::Int64)
    a = zeros(n+2,n+2)
    for i in 1:n+2
        for j in 1:n+2
            a[i,j] = rand([i for i in 2:9])
        end
    end
    return a
end