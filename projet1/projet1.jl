using JuMP
using CPLEX

include("./constants.jl")


"""
include("./projet1/projet1.jl")
reserveSolve("taxe_grille_2x3.txt")
"""
function neighbors(i::Int64, j::Int64)::Vector{Tuple{Int64, Int64}}
    return [(i+1, j), (i+1, j+1), (i+1, j-1), (i, j-1), (i,j+1), (i-1,j), (i-1, j+1), (i-1, j-1)]
end

function returnAlphas(case::Int64)::Vector{Float64}
    if case == 1:
        return [0.5 for _ in 1:6]
    elseif case == 2
        return vcat([0.9 for _ in 1:3],[0.5 for _ in 1:3])
    elseif case == 3
        return vcat([0.5 for _ in 1:3],[0.9 for _ in 1:3])
    elseif case == 4
        return vcat([0.8 for _ in 1:3],[0.6 for _ in 1:3])
    else
        println("Unknown case")
    end
end

function reserveSolve(case::Int64, silent::Bool=true)::Any
    """
    """
    # Directly load data file
    # TODO parsing
    include(PROJET_1_PATH * "\\" * inputFile)
    n = 10
    alphas = returnAlphas(case)

    reserve = [(i,j) for i in 2:n+1 for j in 2:n+1]
    border = vcat([(i,j) for i in 1:n+2 for j in [1, n+2]], [(i,j) for i in [1,n+2] for j in 2:n+1])
    #totalReserve = [(i,j) for i in 1:n+2 for j in 1:n+2]

    dangerSpecies = [1,2,3]
    communeSpecies = [4,5,6]

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ### Variables
    # Variable for sites
    @variable(model, x[i in 1:n+2, j in 1:n+2], Bin)
    # Variables for central zones
    @variable(model, c[i in 1:n+2, j in 1:n+2], Bin)

    # Objective : sum on all sites if they are used
    @objective(model, Max, sum(a[i,j] * x[i,j] for (i,j) in reserve))

    ### Constraints
    # Links between x and c
    @constraint(model, [(i,j) in reserve], x[i,j] >= c[i,j])
    @constraint(model, [(i,j) in reserve], 8*c[i,j] <= sum(x[l,k] for (l,k) in neighbors(i,j)))
    # x and c must be empty on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)
    @constraint(model, [(i,j) in border], c[i,j] == 0)
    # Constraint on probabilities for species
    @constraint(model, [k in dangerSpecies], sum(c[i,j]*log(1-p[k,i,j]) for (i,j) in reserve) <= log(1- alphas[k]))
    @constraint(model, [k in communeSpecies], sum(x[i,j]*log(1-p[k,i,j]) for (i,j) in reserve) <= log(1- alphas[k]))

end