using JuMP
using CPLEX

"""
include("./projet5/projet5.jl")
rollingSolve()
"""

function buildGraph(d::Union{Vector, Nothing}=reportD;Emax::Int64=3, T::Int64=12, modes::Int64=4)
    Obj = Array{Float64}(undef, 0)
    C = Array{Float64}(undef, 0)
    
    if d == nothing
        d = [rand([i for i in 20:70]) for j in 1:T]
    end
    for i in 1:T
        if modes == 4
            o,c = rollingSolve(d, R=i, Emax=Emax, T=T)
        elseif modes == 6
            o,c = rollingSolve6Modes(d, R=i, Emax=Emax, T=T)
        elseif modes == 8
            o,c = rollingSolve8Modes(d, R=i, Emax=Emax, T=T)
        elseif modes == 10
            o,c = rollingSolve10Modes(d, R=i, Emax=Emax, T=T)
        end
        append!(Obj, o)
        append!(C,c)
    end
    println("Objs : ", Obj)
    println("C :", C)
end

commonD = [31, 62, 70, 25, 36, 20, 52, 32, 59, 43, 51, 29]
reportD= [69, 25, 26, 54, 39, 66, 67, 48, 36, 41, 70, 67]
eighteen_D = [59, 28, 23, 25, 53, 20, 27, 31, 29, 30, 64, 24, 65, 23, 42, 51, 60, 27]
twentyfour_D = [63, 69, 40, 59, 35, 20, 32, 57, 58, 27, 59, 69, 27, 37, 20, 28, 37, 33, 60, 61, 70, 60, 32, 54]
thirty_D = [23, 45, 61, 31, 66, 29, 40, 49, 36, 52, 33, 30, 56, 42, 63, 63, 36, 59, 41, 32, 52, 28, 53, 55, 69, 66, 37, 24, 57, 59]

function rollingSolve(d=reportD; R::Int64=2, Emax::Int64=3, T::Int64=12, M::Int64=4, f::Vector{Int64}=[10, 30, 60, 90], e::Vector{Int64}=[8, 6, 4, 2], display::Bool=false)
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

function rollingSolve6Modes(d=reportD; R::Int64=2, Emax::Int64=3, T::Int64=12, display::Bool=false)
    e6 = [8, 7, 6, 4, 3, 2]
    f6 = [10, 20, 30, 60, 70, 90]
    rollingSolve(d, R=R, Emax=Emax, T=T, M=6, f=f6, e=e6, display=display)
end

function rollingSolve8Modes(d=reportD; R::Int64=2, Emax::Int64=3, T::Int64=12, display::Bool=false)
    e8 = [8, 7, 6, 5, 4, 3, 2, 1]
    f8 = [10, 20, 30, 45, 60, 70, 90, 120]
    rollingSolve(d, R=R, Emax=Emax, T=T, M=8, f=f8, e=e8, display=display)
end

function rollingSolve10Modes(d=reportD; R::Int64=2, Emax::Int64=3, T::Int64=12, display::Bool=false)
    e10 = [12, 8, 7, 6, 5, 4, 3, 2, 1, 0]
    f10 = [5, 10, 20, 30, 45, 60, 70, 90, 120, 150]
    rollingSolve(d, R=R, Emax=Emax, T=T, M=10, f=f10, e=e10, display=display)
end