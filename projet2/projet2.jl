using JuMP
using CPLEX



"""
include("./projet2/projet2.jl")
include("./projet2/data/data.jl")
dinkelbach("./projet2/data/data.jl")
"""
 
function distance(i::Int64,j::Int64)::Float64
    i_line, i_col = div(i -1,10) + 1, (i-1)%10 + 1
    j_line, j_col = div(j - 1,10) + 1, (j-1)%10 + 1

    return sqrt((i_line - j_line)^2 + (i_col - j_col)^2)
end

function dinkelbach(inputFile::String="data", showResult::Bool= false, silent::Bool=true)::Any
    """
    n : grid edge size
    """
    
    include("./projet2/data/" * inputFile * ".jl")
    numberOfSites = n*n

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ##### Variables #####
    @variable(model, x[i in 1:numberOfSites], Bin)
    @variable(model, y[i in 1:numberOfSites, j in 1:numberOfSites], Bin)
    
    ##### Constraints #####
    # Area constraints
    @constraint(model, sum(x[i] for i in 1:numberOfSites) <= Amax)
    @constraint(model, sum(x[i] for i in 1:numberOfSites) >= Amin)
    # Costs constraint
    @constraint(model, sum(x[i]*c[i]*10 for i in 1:numberOfSites) <= B)
    # Constraint for y
    for i in 1:numberOfSites
        @constraint(model, sum(y[i,j] for j in 1:numberOfSites) == x[i])
        @constraint(model, y[i,i] == 0)
    end
    for j in 1:numberOfSites
        for i in 1:numberOfSites
            # Redundant
            #@constraint(model, y[i,j] <= x[i])
            @constraint(model, y[i,j] <= x[j])
        end
    end

    currentLambda = lambda
    
    done = false
    iteration = 0
    while true
        iteration += 1
        println()
        println("Iteration " * string(iteration))
        println("Current lambda is " * string(currentLambda))
        ##### Objective #####
        @objective(model, Min, sum(y[i,j]*distance(i,j) for i in 1:numberOfSites for j in 1:numberOfSites if i!=j)
         - currentLambda *sum(x[i] for i in 1:numberOfSites))
        # Solve current state
        optimize!(model)

        feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
        isOptimal = termination_status(model) == MOI.OPTIMAL
        if feasibleSolutionFound && isOptimal
            value = JuMP.objective_value(model)
            println("Current value is " * string(value))
        else
            println("Not feasible or not optimal!!!")
            return
        end

        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)

        fractValue = sum(y_val[i,j]*distance(i,j) for i in 1:numberOfSites for j in 1:numberOfSites if i!=j) / sum(x_val[i] for i in 1:numberOfSites)

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
