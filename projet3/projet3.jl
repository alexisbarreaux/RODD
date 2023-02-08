using JuMP
using CPLEX


"""
include("./projet3/projet3.jl")
genetics("./projet3/data/data.jl")
"""

function createThetas(R::Int64=50, theta_1::Float64=1e-3)
    return [theta_1^((h - r)/(h - 1)) for r in 1:R]
end

function genetics(inputFile::String="data", showResult::Bool= false, silent::Bool=true)::Any
    """
    P : number of individuals
    nbMale : number of males
    """
    R = 50
    thetas = createThetas(R)
    include("./projet3/data/" * inputFile * ".jl")
    numberOfSites = n*n

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ##### Variables #####
    @variable(model, x[i in 1:nbIndividuals] >= 0, Int)
    @variable(model, pi[i in 1:nbGenes, j in 1:nbAlleles] >= 0)
    @variable(model, y[i in 1:nbGenes, j in 1:nbAlleles], Bin)


    ##### Objective #####
    @objective(model, Min, sum(pi[i,j] for i in 1:nbGenes for j in 1:nbAlleles))

    ##### Constraints #####
    # Valid children numbers constraint
    @constraint(model, sum(x[i] for i in 1:nbMale) == sum(x[i] for i in nbMale+1:nbIndividuals))
    # Constraint on the number of children per individual
    @constraint(model, [i in 1:nbIndividuals], x[i] <= 3)
    # Number of individuals stays constant
    @constraint(model, sum(x[i] for i in 1:nbIndividuals) == nbIndividuals)

    # Linearization constraints
    @constraint(model, [i in 1:nbGenes, j in 1:nbAlleles], y[i,j] >= sum(x[p] for p in 1:nbIndividuals for i in 1:nbGenes 
                                for j in 1:nbIndividuals if proportion[p,i,j]==1) / nbIndividuals)
    # Linearized log constraints
    for r in 1:R 
        @constraint(model, [i in 1:nbGenes, j in 1:nbAlleles], log(theta[r]) + (1/theta[r])*(pi[i,j] - y[i,j] - theta[r])>=  
                        sum(x[p]*log(1-proportion[p,i,j])))
    end

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound
        value = JuMP.objective_value(model)
        x_val = JuMP.value.(x)
        solveTime = round(JuMP.solve_time(model), digits= 5)
        gap = JuMP.relative_gap(model)
        bound = JuMP.objective_bound(model)
        println("Current value is " * string(value))
        return
    else
        println("Not feasible or not optimal!!!")
        return
    end
end
