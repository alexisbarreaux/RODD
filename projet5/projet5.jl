using JuMP
using CPLEX

function rollingSolve()
    T=12 #horizon de temps
    M=4 #nombre de modes
    Emax = 3 #émission carbone maximum à chaque période
    d = [rand([i for i in 20:70]) for j in 1:T] #demande de chaque période
    f = [10, 30, 60, 90] #cout d'approvisionnement de chaque mode
    e = [8, 6, 4, 2] #émission carbone de chaque mode
    h = 1 
    p = 0
    R = 2 #intervalle du rolling

    model = Model(CPLEX.Optimizer)
    set_silent(model)

    ### Variables
    @variable(model, x[t in 1:T, m in 1:M]>= 0)
    @variable(model, s[t in 1:T]>= 0)
    @variable(model, y[t in 1:T, m in 1:M], Bin)

    ### Constraints 
    @constraint(model, sum(x[1,m] - s[1]  for m in 1:M) == d[1])
    @constraint(model, [t in 2:T], sum(x[t,m] - s[t] + s[t-1] for m in 1:M) == d[t])
    @constraint(model, [t in 1:T, m in 1:M], x[t,m] <= (sum(d[t2] for t2 in t:T)*y[t,m]))
    @constraint(model, [t in R:T], sum(sum(e[m] - Emax  for m in 1:M)*x[t2,m] for t2 in (t-R+1):t) <= 0)

    @objective(model, Min, sum(sum(p*x[t,m] + f[m]*y[t,m] for t in 1:T) for m in 1:M) + sum(h*s[t] for t in 1:T))

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL


    if feasibleSolutionFound && isOptimal
        value = round(JuMP.objective_value(model), digits=2)
        x_val = JuMP.value.(x)
        y_val = JuMP.value.(y)
        s_val = JuMP.value.(s)

        solveTime = round(JuMP.solve_time(model), digits= 5)
        nodes = JuMP.node_count(model)
        bound = JuMP.objective_bound(model)

        println("Optimal value is " * string(value))
        println("Nodes " * string(nodes))
        println("time : ", round(solveTime, digits=2),"sec")
        println("X : ", x_val)
        println("Y : ", y_val)
        println("S : ", s_val)
    else
        println("Problem is not feasible !!!")
        return
    end
end
