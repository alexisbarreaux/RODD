using JuMP
using CPLEX

include("./parser.jl")
include("./generateDatas.jl")

"""
include("./projet4/projet4_P2.jl")
P2_bis()
"""

function belowRightNeighbors(i::Int64, j::Int64)::Vector{Tuple{Int64, Int64}}
    return [(i+1, j), (i,j+1)]
end

function neighbors(i::Int64, j::Int64, m::Int64, n::Int64)::Vector{Tuple{Int64, Int64}}
    res = []
    if i > 1
        res = vcat(res, [(i-1, j)])
    end
    if j > 1
        res = vcat(res, [(i, j-1)])
    end
    if i < m + 2
        res = vcat(res, [(i+1, j)])
    end
    if j < n +2
        res = vcat(res, [(i, j + 1)])
    end
    return res
end

function P2(instance::Int64=1, showResult::Bool= false, silent::Bool=true)::Any
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
    @variable(model, x[i in 1:m+2, j in 1:n+2] >=0)
    @variable(model, y[i in 1:m+2, j in 1:n+2, (k,l) in neighbors(i,j,m,n)] >= 0)
    
    ##### Objective #####
    @objective(model, Max, sum(w1*t[i,j]*(1 - x[i,j]) for i in 2:m+1, j in 2:n+1) + w2*g*l*
    sum(x[i,j] + x[k,l] - 2*y[i,j,(k,l)] for i in 1:m+1, j in 1:n+1, (k,l) in neighbors(i,j,m,n)))
    ##### Constraints #####
    # Constraints on y
    @constraint(model, [i in 1:m+1, j in 1:n+1, (k,l) in neighbors(i,j,m,n)], y[i,j,(k,l)] >= x[i,j] + x[k,l] - 1)
    # x is zero on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)
    @constraint(model, [i in 1:m+2, j in 1:n+2], x[i,j] <= 1)
    
    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound && isOptimal
        value = round(JuMP.objective_value(model), digits=2)
        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)
        
        solveTime = round(JuMP.solve_time(model), digits= 5)
        
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
        println("Nombre de parcelles non coupées " * string(round(sum(x_val))))
        println("Value is " * string(value) * " in " * string(solveTime) * " s and " * string(JuMP.node_count(model))* " nodes.")
        
        
        return
    else
        println("Not feasible!")
        return
    end
end


function P2_bis(instance::Int64=1, showResult::Bool= false, silent::Bool=true)::Any
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
    @variable(model, 0 <= x[i in 1:m+2, j in 1:n+2] <= 1)
    @variable(model, y[i in 1:m+2, j in 1:n+2, (k,l) in belowRightNeighbors(i,j)] >=0)

    ##### Objective #####
    leftTerm  = w1*sum(t[i,j]*(1 - x[i,j]) for i=2:m+1, j=2:n+1)
    righTerm = 0
    for i=1:m+1
        for j=1:n+1
            for (k,l) in belowRightNeighbors(i,j)
                righTerm += x[i,j] + x[k,l] - 2*y[i,j,(k,l)]
            end
        end
    end
    righTerm = w2*g*l*righTerm

    @objective(model, Max, leftTerm + righTerm)
    ##### Constraints #####
    # Constraints on y
    @constraint(model, [i in 1:m+1, j in 1:n+1, (k,l) in belowRightNeighbors(i,j)], y[i,j,(k,l)] >= x[i,j] + x[k,l] - 1)
    # x is zero on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL
    if feasibleSolutionFound && isOptimal
        value = round(JuMP.objective_value(model), digits=2)
        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)
        
        solveTime = round(JuMP.solve_time(model), digits= 5)
        #gap = JuMP.relative_gap(model)
        #bound = JuMP.objective_bound(model)
        
        for i in 1:m+2
            println()
            for j in 1:n+2
                if x_val[i,j] >= 1 - 1e-6
                    print("C ")
                else
                    print("_ ")
                end
            end
        end
        println()

        println(sum(y_val))
        println("Nombre de parcelles non coupées " * string(round(sum(x_val))))
        println("Value is " * string(value) * " in " * string(solveTime) * " s and " * string(JuMP.node_count(model))* " nodes.")
        
        
        return
    else
        println("Not feasible!")
        return
    end
end