

function parser(inst::String="data/DivGenetique.txt")
    f = open(inst)
    for i in 1:5
        readline(f)
    end
    M = zeros(8,5,2)
    line = readline(f)
    while true
        if line==""
            line = readline(f)
        end
        if line==""
            break
        end
        line = split(line, " ")
        println(line)
        M[parse(Int64,line[1]),parse(Int64,line[3]), parse(Int64,line[5])] += 0.5
        line = readline(f)
    end
    return 8, 5, 2, M
end
    


    