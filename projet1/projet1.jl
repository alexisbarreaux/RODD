using JuMP
using CPLEX

include("./constants.jl")
include("./parser.jl")
include("./generateDatas.jl")


"""
include("./projet1/projet1.jl")
reserveSolve(1)
randomReserveSolve(1, 20,4)
"""
function neighbors(i::Int64, j::Int64)::Vector{Tuple{Int64, Int64}}
    return [(i+1, j), (i+1, j+1), (i+1, j-1), (i, j-1), (i,j+1), (i-1,j), (i-1, j+1), (i-1, j-1)]
end

function returnAlphas(case::Int64,k::Int64=6, e::Int64=3)::Vector{Float64}
    if case == 1
        return [0.5 for _ in 1:k]
    elseif case == 2
        return vcat([0.9 for _ in 1:e],[0.5 for _ in 1:(k-e)])
    elseif case == 3
        return vcat([0.5 for _ in 1:e],[0.9 for _ in 1:(k-e)])
    elseif case == 4
        return vcat([0.8 for _ in 1:e],[0.6 for _ in 1:(k-e)])
    elseif case == 5
        return vcat([0.9 for _ in 1:e],[0.9 for _ in 1:(k-e)])
    else
        println("Unknown case")
    end
end

function reserveSolve(case::Int64, instance::String= "probabilities.txt",silent::Bool=true)::Any
    """
    """
    # Directly load data file
    alphas = returnAlphas(case)
    a = cost()
    p, n = parser(PROJET_1_PATH * "\\" * instance)

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
    @objective(model, Min, sum(a[i,j] * x[i,j] for (i,j) in reserve))

    ### Constraints
    # Links between x and c
    @constraint(model, [(i,j) in reserve], x[i,j] >= c[i,j])
    @constraint(model, [(i,j) in reserve], 8*c[i,j] <= sum(x[l,k] for (l,k) in neighbors(i,j)))
    @constraint(model, [(i,j) in reserve],  c[i,j] >= sum(x[l,k] for (l,k) in neighbors(i,j)) - 7)
    # x and c must be empty on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)
    @constraint(model, [(i,j) in border], c[i,j] == 0)
    # Constraint on probabilities for species
    @constraint(model, [k in dangerSpecies], sum(c[i,j]*log(1-p[k,i,j]) for (i,j) in reserve) <= log(1- alphas[k]))
    @constraint(model, [k in communeSpecies], sum(x[i,j]*log(1-p[k,i,j]) for (i,j) in reserve) <= log(1- alphas[k]))

    """
    f = open("model.lp", "w")
    print(f, model)
    close(f)
    """

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL

    if feasibleSolutionFound
        x_val = JuMP.value.(x)
        c_val = JuMP.value.(c)
        dangerSpecies = [1,2,3]
        communeSpecies = [4,5,6]
        for i in 2:n+1
            println()
            for j in 2:n+1
                if c_val[i,j] == 1
                    print("C ")
                elseif x_val[i,j] == 1
                    print("T ")
                else
                    print("_ ")
                end
            end
        end
        println()

        for k in dangerSpecies
            println("Specie " * string(k) * " survival rate " * string(round(1 - prod([1- p[k,i,j]*c_val[i,j] for (i,j) in reserve]),digits=4)))
        end
        for k in communeSpecies
            println("Specie " * string(k) * " survival rate " * string(round(1 - prod([1- p[k,i,j]*x_val[i,j] for (i,j) in reserve]),digits=4)))
        end

        # Result returning
        print("Nodes " * string(JuMP.node_count(model)) * " ")
        print("Time " * string(round(JuMP.solve_time(model), digits= 5)) * " ")
        print("Objective " * string(JuMP.objective_value(model)) * " ")
        println()
        return JuMP.objective_value(model)
    else
        println("Problem is not feasible !!!")
        return
    end
end

function randomReserveSolve(case::Int64, n::Int64, k::Int64, silent::Bool=true)::Any
    """
    """
    # Directly load data file
    e = Int64(round(k/2))
    alphas = returnAlphas(case, k, e)
    a = generateCost(n)
    p = generateDatas(n,k,e)

    reserve = [(i,j) for i in 2:n+1 for j in 2:n+1]
    border = vcat([(i,j) for i in 1:n+2 for j in [1, n+2]], [(i,j) for i in [1,n+2] for j in 2:n+1])
    #totalReserve = [(i,j) for i in 1:n+2 for j in 1:n+2]

    dangerSpecies = [i for i in 1:e]
    communeSpecies = [i for i in e+1:k]

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
    @objective(model, Min, sum(a[i,j] * x[i,j] for (i,j) in reserve))

    ### Constraints
    # Links between x and c
    @constraint(model, [(i,j) in reserve], x[i,j] >= c[i,j])
    @constraint(model, [(i,j) in reserve], 8*c[i,j] <= sum(x[l,k] for (l,k) in neighbors(i,j)))
    @constraint(model, [(i,j) in reserve],  c[i,j] >= sum(x[l,k] for (l,k) in neighbors(i,j)) - 7)
    # x and c must be empty on border
    @constraint(model, [(i,j) in border], x[i,j] == 0)
    @constraint(model, [(i,j) in border], c[i,j] == 0)
    # Constraint on probabilities for species
    @constraint(model, [k in dangerSpecies], sum(c[i,j]*log(1-p[k,i,j]) for (i,j) in reserve) <= log(1- alphas[k]))
    @constraint(model, [k in communeSpecies], sum(x[i,j]*log(1-p[k,i,j]) for (i,j) in reserve) <= log(1- alphas[k]))

    """
    f = open("model.lp", "w")
    print(f, model)
    close(f)
    """

    optimize!(model)
    feasibleSolutionFound = primal_status(model) == MOI.FEASIBLE_POINT
    isOptimal = termination_status(model) == MOI.OPTIMAL

    if feasibleSolutionFound
        x_val = JuMP.value.(x)
        c_val = JuMP.value.(c)
        for i in 2:n+1
            println()
            for j in 2:n+1
                if c_val[i,j] >= 1 - 1e-6
                    print("C ")
                elseif x_val[i,j] >= 1 - 1e-6
                    print("T ")
                else
                    print("_ ")
                end
            end
        end
        println()

        for k in dangerSpecies
            println("Specie " * string(k) * " survival rate " * string(round(1 - prod([1- p[k,i,j]*c_val[i,j] for (i,j) in reserve]),digits=4)))
        end
        for k in communeSpecies
            println("Specie " * string(k) * " survival rate " * string(round(1 - prod([1- p[k,i,j]*x_val[i,j] for (i,j) in reserve]),digits=4)))
        end

        # Result returning
        print("Nodes " * string(JuMP.node_count(model)) * " ")
        print("Time " * string(round(JuMP.solve_time(model), digits= 5)) * " ")
        print("Objective " * string(JuMP.objective_value(model)) * " ")
        println()
        return JuMP.objective_value(model)
    else
        println("Problem is not feasible !!!")
        return
    end
end
