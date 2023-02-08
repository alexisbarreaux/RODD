using DataFrames
using CSV

DIR_PATH = @__DIR__

RESULTS_DIR_PATH = DIR_PATH * "\\results"

"""
include("./projet3/main.jl")
"""

include("./projet3_generated.jl")


function buildResultsDf( resultFile::String="test")::Nothing
    # Loading
    filePath =RESULTS_DIR_PATH * "\\" * resultFile * ".csv"
    results = DataFrame( nbIndividuals=Int64[],
     time=Float64[], nodes=Int64[], value=Float64[], bound=Float64[])

     nbGenes = 100
     nbAlleles = 3
     
        for nbIndividuals in 2:10:60
        println("Running with " * string(nbIndividuals) * " individuals.")
        nbMale = div(nbIndividuals, 2)
        solveTime, nodes, realValue, bound = geneticsGenerated(nbIndividuals, nbMale, nbGenes, nbAlleles)
        push!(results, [nbIndividuals solveTime nodes realValue bound])
    end

    # Run
    CSV.write(filePath, results, delim=";")

    return 
end