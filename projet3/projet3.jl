using JuMP
using CPLEX


"""
include("./projet3/projet3.jl")
genetics("./projet3/data/data.jl")
"""

function genetics(inputFile::String="data", showResult::Bool= false, silent::Bool=true)::Any
    """
    P : number of individuals
    nbMale : number of males
    """
    include("./projet3/data/" * inputFile * ".jl")
    numberOfSites = n*n

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ##### Variables #####
    @variable(model, x[i in 1:nbIndividuals] >= 0, Int)
    @variable(model, pi[i in 1:nbGenes, j in 1:nbAlleles], Bin)

    ##### Objective #####
    @objective(model, Min, sum(pi[i,j] for i in 1:nbGenes for j in 1:nbAlleles))

    ##### Constraints #####
    # Valid children numbers constraint
    @constraint(model, sum(x[i] for i in 1:nbMale) == sum(x[i] for i in nbMale+1:nbIndividuals))
    # Constraint on the number of children per individual
    @constraint(model, [i in 1:nbIndividuals], x[i] <= 3)

    # Linearization constraints
    @constraint(model, [i in 1:nbIndividuals], x[i] <= 3)

    return

    currentLambda = lambda
    
    done = false
    iteration = 0
    while true
        iteration += 1
        println()
        println("Iteration " * string(iteration))
        println("Current lambda is " * string(currentLambda))
        ##### Objective #####
        # Solve current state
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
        else
            println("Not feasible or not optimal!!!")
            return
        end

        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)

        fractValue = sum(y_val[i,j]*distance(i,j,n) for i in 1:numberOfSites for j in 1:numberOfSites if i!=j) / sum(x_val[i] for i in 1:numberOfSites)

        if value <= 0 - 1e-6
            currentLambda = fractValue
            continue
        else
            for i in 1:n
                println()
                for j in 1:n
                    if x_val[10*(i-1) + j] >= 1 - 1e-6
                        print("C ")
                    else
                        print("_ ")
                    end
                end
            end
            println()
            println("Optimal value is " * string(fractValue))
            println("Nodes " * string(JuMP.node_count(model)))
            println("Number of sites " * string(round(sum(x_val))))
            println("Sum of y " *  string(sum(y_val)))
            println("time : ", round(time()-start, digits=2),"sec")

            """
            println()
            println()
            for i in 1:numberOfSites
                for j in 1:numberOfSites
                    if y_val[i,j] > 1 - 1e-6
                        println(i," ", j," ", y_val[i,j], " ", distance(i,j))
                    end
                end
            end
            """

            
            return
        end
    end
end
