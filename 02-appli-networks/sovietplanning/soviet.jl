#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Soviet Planning

=#

# Modules used


using JuMP, Gurobi, DataFrames, CSV

# Getting the data

path = pwd()
filename = "distances.csv"
data = CSV.read(string(path, "/", filename), weakrefstrings = false)

# getting nodes names

colnames = Array(CSV.read(string(path, /, filename), rows = 1, datarow = 1)[1, 2:11])
colnames = reshape(colnames, length(colnames), 1)

# Organizing the data

nsources = 68
ndests = 10

dists = Array(data[1:68, 2:11]) # data in matrix form
p = Array(data[1:68, 12])
q = Array(data[69, 2:11])

nonzeros = find(!ismissing, dists) # finding the non zeros entries

nbNodes = nsources+ndests
nbArcs = length(nonzeros)

rows = mod.(nonzeros-1,nsources) + 1
cols = div.(nonzeros-1,nsources) + 1
costs = dists[nonzeros] # cost vector
arcs = [rows,cols+nsources,costs] # arcs vector

n = vcat(-p, q')

nameNodes = vcat(Array(data[1:nsources, 1]), colnames)

# construct node-incidence matrix :

Nabla = sparse(1:nbArcs, arcs[1], -1, nbArcs, nbNodes) + sparse(1:nbArcs, arcs[2], 1, nbArcs, nbNodes)

# setting and solving the model

m = Model(solver=GurobiSolver(Presolve=0))
@variable(m, x[1:nbArcs])
@objective(m, Min, costs'*x)
@constraint(m, Nabla'*x .== n)
@constraint(m, x .>= 0.0)
status = solve(m)

π = getvalue(x)
total_distance = getobjectivevalue(m)

println()
println("******************")
println("Optimal Solution")
println()
println("Optimal Distance = $total_distance")
println()
println("******************")
