

function parser(inst::String, n::Int64=10, k::Int64=6)
    f    = open(inst)
    readline(f)
    readline(f)
    p = zeros(k,n+2,n+2)
    cpt = 1
    line = readline(f)
    while cpt <= k
        line = split(line, " ")
        if line[1] == "\t" || line[1] == ""
            cpt +=1
            line = readline(f)
            continue
        end
        p[cpt, parse(Int,line[2]), parse(Int,line[3])] = parse(Float64,line[4])
        line = readline(f)
    end
    return p, n 
end


function cost()
    a = [ 
        6 6 6 4 4 4 4 8 8 8 ;
        6 6 6 4 4 4 4 8 8 8 ;
        6 6 6 4 4 4 4 8 8 8 ;
        5 5 5 3 3 3 3 7 7 7 ;
        5 5 5 3 3 3 3 7 7 7 ;
        5 5 5 3 3 3 3 7 7 7 ;
        5 5 5 3 3 3 3 7 7 7 ;
        4 4 4 6 6 6 6 5 5 5 ;
        4 4 4 6 6 6 6 5 5 5 ;
        4 4 4 6 6 6 6 5 5 5 
    ]
    return a
end



