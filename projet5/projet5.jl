using JuMP
using CPLEX

"""
include("./projet5/projet5.jl")
rollingSolve()
"""

function buildGraph(;Emax::Int64=3)
    T=12
    Obj = Array{Float64}(undef, 0)
    C = Array{Float64}(undef, 0)
    d = reportD
    #d = [rand([i for i in 20:70]) for j in 1:T]
    for i in 1:12
        o,c = rollingSolve(R=i, Emax=Emax)
        append!(Obj, o)
        append!(C,c)
    end
    println("Objs : ", Obj)
    println("C :", C)
end

commonD = [31, 62, 70, 25, 36, 20, 52, 32, 59, 43, 51, 29]
reportD= [69, 25, 26, 54, 39, 66, 67, 48, 36, 41, 70, 67]
function rollingSolve(d=reportD; R::Int64=2, Emax::Int64=3, display::Bool=false)
    T=12 #horizon de temps
    M=4 #nombre de modes
    f = [10, 30, 60, 90] #cout d'approvisionnement de chaque mode
    e = [8, 6, 4, 2] #émission carbone de chaque mode
    h = 1 
    p = 0

    model = Model(CPLEX.Optimizer)
    set_silent(model)

    ### Variables
    @variable(model, x[t in 1:T, m in 1:M]>= 0)
    @variable(model, s[t in 0:T]>= 0)
    @variable(model, y[t in 1:T, m in 1:M], Bin)

    ### Constraints
    @constraint(model, s[0]==0)
    @constraint(model, s[T]==0)
    @constraint(model, [t in 1:T], sum(x[t,m] for m in 1:M) - s[t] + s[t-1] == d[t])
    @constraint(model, [t in 1:T, m in 1:M], x[t,m] <= (sum(d[t2] for t2 in t:T)*y[t,m]))
    # Périodique
    #@constraint(model, [t in 1:T], sum((e[m] - Emax)  *x[t,m] for m in 1:M) <= 0)
    # Glissant
    @constraint(model, [t in R:T], sum(sum((e[m] - Emax) *x[t2,m] for m in 1:M) for t2 in (t-R+1):t) <= 0)
    # Cumulatif
    #@constraint(model, [t in 1:T], sum(sum((e[m] - Emax)*x[t2,m]  for m in 1:M) for t2 in 1:t) <= 0)
    # Global
    #@constraint(model, sum(sum((e[m] - Emax)*x[t2,m]  for m in 1:M) for t2 in 1:T) <= 0)

    @objective(model, Min, sum(sum(p*x[t,m] + f[m]*y[t,m] for t in 1:T) for m in 1:M) + sum(h*s[t] for t in 1:T))

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL


    if feasibleSolutionFound && isOptimal
        value = round(JuMP.objective_value(model), digits=2)
        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)
        s_val = JuMP.value.(s)

        C = sum(sum((x_val[t,m]*e[m]) for m in 1:M) for t in 1:T)/ sum(sum((x_val[t,m]) for m in 1:M) for t in 1:T)

        solveTime = round(JuMP.solve_time(model), digits= 5)
        nodes = JuMP.node_count(model)
        bound = JuMP.objective_bound(model)
        if display
            println("Optimal value is " * string(value))
            println("Nodes " * string(nodes))
            println("time : ", round(solveTime, digits=2),"sec")
            println("X : ", x_val)
            println("Y : ", y_val)
            println("S : ", s_val)
        end
        return value, C
    else
        println("Problem is not feasible !!!")
        return
    end
end
