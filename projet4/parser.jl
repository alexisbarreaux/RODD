

function parser(inst::String="projet4/data/ExplForet.txt")
    f = open(inst)
    line = readline(f)
    line = split(line, " ")
    size = parse(Int64,line[length(line)-1])
    t = zeros(size+2, size+2)
    for i in 1:size
        line = readline(f)
        line = split(line, " ")
        for j in 1:size
            t[i+1,j+1] = parse(Int64, line[j+1])
        end
    end
    return t
end
