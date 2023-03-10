using JuMP
using CPLEX

include("./parser.jl")
include("./generateDatas.jl")

"""
include("./projet4/projet4_P1.jl")
P1()
"""

function neighbors(i::Int64, j::Int64)::Vector{Tuple{Int64, Int64}}
    return [(i+1, j), (i, j-1), (i,j+1), (i-1,j)]
end


function P1(instance::Int64=1, showResult::Bool= false, silent::Bool=true)::Any
    """
    """
    if instance==1
        m = n = 10
        w1 = 1
        w2 = 5
        l = 3
        g = 1.26157
        t = parser()
    else
        m = n = 5
        w1 = 2
        w2 = 1
        l = 3
        g = 1.26157
        t = instance2()
    end

    border = vcat([(i,j) for i in 1:m+2 for j in [1, n+2]], [(i,j) for i in [1,m+2] for j in 2:n+1])

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ##### Variables #####
    @variable(model, x[i in 1:m+2, j in 1:n+2], Bin)
    @variable(model, d[i in 2:m+1, j in 2:n+1] >= 0)

    ##### Objective #####
    @objective(model, Max, sum(w1*t[i,j]*(1 - x[i,j]) for i in 2:m+1, j in 2:n+1) + w2*g*l*sum(4*x[i,j] - d[i,j] for i in 2:m+1, j in 2:n+1))
    
    ##### Constraints #####
    # Constraint on d
    @constraint(model, [i in 2:m+1, j in 2:n+1], d[i,j] >= sum(x[k,l] for (k,l) in neighbors(i,j)) - 4 * (1 - x[i,j]) )
    # x is zero on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)
    #@constraint(model, sum(x[i,j] for i in 2:m+1, j in 2:n+1) >= 60)

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound && isOptimal
        value = round(JuMP.objective_value(model), digits=2)
        x_val = JuMP.value.(x)
        d_val = JuMP.value.(d)
        
        solveTime = round(JuMP.solve_time(model), digits= 5)
        gap = JuMP.relative_gap(model)
        bound = JuMP.objective_bound(model)

        for i in 2:n+1
            println()
            for j in 2:n+1
                if x_val[i,j] >= 1 - 1e-6
                    print("C ")
                else
                    print("_ ")
                end
            end
        end

        println()
        println("Nombre de parcelles non coup??es " * string(round(sum(x_val))))
        println("Value is " * string(value) * " in " * string(solveTime) * " s and " * string(JuMP.node_count(model))* " nodes.")


        return
    else
        println("Not feasible!")
        return
    end
end

function P1_gen(t,n::Int64, showResult::Bool= false, silent::Bool=true)::Any
    """
    """
    m = n
    w1 = 1
    w2 = 5
    l = 3
    g = 1.26157

    border = vcat([(i,j) for i in 1:m+2 for j in [1, n+2]], [(i,j) for i in [1,m+2] for j in 2:n+1])

    # Creating the model
    model = Model(CPLEX.Optimizer)
    if silent
        set_silent(model)
    end

    ##### Variables #####
    @variable(model, x[i in 1:m+2, j in 1:n+2], Bin)
    @variable(model, d[i in 2:m+1, j in 2:n+1] >= 0)

    ##### Objective #####
    @objective(model, Max, sum(w1*t[i,j]*(1 - x[i,j]) for i in 2:m+1, j in 2:n+1) + w2*g*l*sum(4*x[i,j] - d[i,j] for i in 2:m+1, j in 2:n+1))
    
    ##### Constraints #####
    # Constraint on d
    @constraint(model, [i in 2:m+1, j in 2:n+1], d[i,j] >= sum(x[k,l] for (k,l) in neighbors(i,j)) - 4 * (1 - x[i,j]) )
    # x is zero on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)
    #@constraint(model, sum(x[i,j] for i in 2:m+1, j in 2:n+1) >= 60)

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound && isOptimal
        value = round(JuMP.objective_value(model), digits=2)
        x_val = JuMP.value.(x)
        d_val = JuMP.value.(d)
        
        solveTime = round(JuMP.solve_time(model), digits= 5)
        gap = JuMP.relative_gap(model)
        bound = JuMP.objective_bound(model)
        nodes= JuMP.node_count(model)
        """
        for i in 2:n+1
            println()
            for j in 2:n+1
                if x_val[i,j] >= 1 - 1e-6
                    print("C ")
                else
                    print("_ ")
                end
            end
        end

        println()
        println("Nombre de parcelles non coup??es " * string(round(sum(x_val))))
        println("Value is " * string(value) * " in " * string(solveTime) * " s and " * string(nodes)* " nodes.")
        """

        return solveTime, nodes, value
    else
        println("Not feasible!")
        return
    end
end
