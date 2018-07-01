#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Optimal Diet Problem

=#

# Packages use

using JuMP, Gurobi, DataFrames, CSV

# Getting the data

path = pwd()
filename = "StiglerData1939.csv"
data = CSV.read(string(path, /, filename); null = "NA")

nbCommodities = count(typeof(data[:Commodity][i]) == String for i in 1:nrow(data))-1

names = data[:Commodity][1:nbCommodities]

mat = Array(data[1:nrow(data), 3:13])
mat = collect(Missings.replace(mat, 0.0))

# Calling Gurobi

N = mat[1:nbCommodities,3:11]'
d = mat[nbCommodities+1,3:11]
c = fill(1.0, nbCommodities)

m = Model(solver=GurobiSolver(Presolve=0))
@variable(m, x[1:nbCommodities])
@objective(m, Min, c'*x)
@constraint(m, N*x .>= d)
@constraint(m, x .>= 0.0)
status = solve(m)

q = getvalue(x)
obj = getobjectivevalue(m)

q_yearly = q*365

# display optimal solution

toKeep = find(!iszero, q_yearly)
foods = q_yearly[toKeep]
foodnames = names[toKeep]
diet = DataFrame(food = foodnames, yearly_cost = foods)
totc = q_yearly'*c

toKeepStigler = [1,15,46,52,69]
foods_Stigler = [13.33, 3.84,4.11,1.85,16.80]
foodnames_Stigler = names[toKeepStigler]
diet_Stigler = DataFrame(food = foodnames_Stigler, yearly_cost = foods_Stigler)
totcS = sum(foods_Stigler)

println("******************")
println("Optimal Solution")
println()
println(diet)
println()
println("Total cost (optimal) = $(round(totc, 2))")
println("******************")
println("Comparing to Stigler's Solution")
println()
println(diet_Stigler)
println()
println("Total cost (Stigler) = $(round(totcS, 2))")
println("******************")
