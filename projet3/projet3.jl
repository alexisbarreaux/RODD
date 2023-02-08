using JuMP
using CPLEX

include("./parser.jl")

"""
include("./projet3/projet3.jl")
genetics("./projet3/data/data.jl")
"""

function createThetas(R::Int64=50, theta_1::Float64=1e-3)
    return [theta_1^((R - r)/(R - 1)) for r in 1:R]
end

function genetics(R::Int64=50, inputFile::String="./projet3/data/data.txt", showResult::Bool= false, silent::Bool=true)::Any
    """
    P : number of individuals
    nbMale : number of males
    """
    thetas = createThetas(R)
    #include("./projet3/data/" * inputFile * ".jl")
    nbIndividuals, nbMale, nbGenes, nbAlleles, proportion=  parser(inputFile)

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ##### Variables #####
    @variable(model, x[i in 1:nbIndividuals] >= 0, Int)
    @variable(model, pi[i in 1:nbGenes, j in 1:nbAlleles] >= 0)
    @variable(model, t[i in 1:nbGenes, j in 1:nbAlleles] >= 0)

    ##### Objective #####
    @objective(model, Min, sum(pi[i,j] for i in 1:nbGenes for j in 1:nbAlleles))

    ##### Constraints #####
    # Valid children numbers constraint
    @constraint(model, sum(x[i] for i in 1:nbMale) == sum(x[i] for i in nbMale+1:nbIndividuals))
    # Constraint on the number of children per individual
    @constraint(model, [i in 1:nbIndividuals], x[i] <= 3)
    #@constraint(model, [i in 1:nbIndividuals], x[i] == 2)

    # Number of individuals stays constant
    @constraint(model, sum(x[i] for i in 1:nbIndividuals) == 2* nbIndividuals)
    
    # Constraints on pi
    @constraint(model, [i in 1:nbGenes, j in 1:nbAlleles], pi[i,j] >= t[i,j] - sum(x[p] for p in 1:nbIndividuals if proportion[p,i,j]==1))

    # Linearized log constraints
    for r in 1:R 
        @constraint(model, [i in 1:nbGenes, j in 1:nbAlleles], log(thetas[r]) + (1/thetas[r])*(t[i,j] - thetas[r])>=  
                        sum(x[p]*log(1-proportion[p,i,j]) for p in 1:nbIndividuals if proportion[p,i,j] != 1))
    end

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound && isOptimal
        value = JuMP.objective_value(model)
        x_val = JuMP.value.(x)
        x_val = [round(x_val[i]) for i in 1:length(x_val)]

        solveTime = round(JuMP.solve_time(model), digits= 5)
        gap = JuMP.relative_gap(model)
        bound = JuMP.objective_bound(model)
        println("Relaxation bound is " * string(value) * " in " * string(solveTime) * " s and " * string(JuMP.node_count(model))* " nodes.")
        println()
        for i in 1:nbIndividuals
            println("Parent " * string(i) * " has " * string(x_val[i]) * " children.")
        end

        println()
        println()
        realValue = 0
        for i in 1:nbGenes
            for j in 1:nbAlleles
                alleleDisappearance = prod((1 - proportion[p,i,j])^x_val[p] for p in 1:nbIndividuals)
                println("Disappearance probability for gene " * string(i) *" allele " * string(j) *" is " * string(alleleDisappearance))
                realValue += alleleDisappearance
            end
        end
        println()
        println("Real value is " * string(realValue))

        return
    else
        println("Not feasible!")
        return
    end
end
