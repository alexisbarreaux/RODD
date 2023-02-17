using DataFrames
using CSV

DIR_PATH = @__DIR__

RESULTS_DIR_PATH = DIR_PATH * "\\results"

"""
include("./projet4/main.jl")
buildResultsDf()
"""

include("./projet4_P1.jl")
include("./projet4_P2.jl")
include("./generateDatas.jl")


function buildResultsDf( resultFile::String="res")::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    results = DataFrame( gridEdge=Int64[],
    time1=Float64[], time2=Float64[], nodes1=Int64[], nodes2=Int64[], value1=Float64[], value2=Float64[])
    
    for n in 10:20:150
        println("Running with " * string(n) * " size.")
        t = buildtij(n,n)
        solveTime1, nodes1, value1 = P1_gen(t,n)
        solveTime2, nodes2, value2 = P2_gen(t,n)
        push!(results, [n solveTime1 solveTime2 nodes1 nodes2 value1 value2])
    end

    # Run
    CSV.write(filePath, results, delim=";")

    return 
end