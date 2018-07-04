#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Optimal Diet Problem

=#

# Packages use

using JuMP, Gurobi, DataFrames, CSV

# Obtaining the data

path = pwd()
file1 = "affinitymatrix.csv"
file2 = "Xvals.csv"
file3 = "Yvals.csv"

nbcar = 10

A = Array(CSV.read(string(path, /, file1))[1:nbcar,2:(nbcar+1)])
Xvals = Array(CSV.read(string(path, /, file2))[1:nbcar])
Yvals= Array(CSV.read(string(path, /, file3))[1:nbcar])

sdX = mapslices(std, Xvals, 1)
mX = mapslices(mean, Xvals, 1)
sdY = mapslices(std, Yvals, 1)
mY = mapslices(mean, Yvals, 1)

Xvals_σ = (Xvals .- mX) ./ sdX
Yvals_σ = (Yvals .- mY) ./ sdY

nobs = size(Xvals_σ)[1]

Phi = Xvals_σ * A * Yvals_σ'

c = vec(Phi)

p = fill(1,nobs)
q = fill(1,nobs)

N = size(Phi)[1]
M = size(Phi)[2]

A1 = kron(fill(1, M)', sparse(1:N, 1:N, 1))
A2 = kron(sparse(1:M, 1:M, 1), fill(1, N)')

B = vcat(A1, A2)

d = vcat(p, q)

# Solving the model

m = Model(solver=GurobiSolver(Presolve=0))
@variable(m, x[1:N, 1:M])
@objective(m, Max, x'*Phi)
@objective(m, x .> 0)
@constraint(m, x*fill(1, M) = p)
@constraint(m, fill(1, N)'*x = q')
status = solve(m)

if status == :Optimal
    π = getvalue(x)
else
    warn("Optimization problem with Gurobi. status = $status")


pi = getvalue(x)
val = getobjectivevalue(m)



println("Value of the problem (Gurobi) = $total_distance")
println()
println()

print(paste0("Value of the problem (Gurobi) =",val))
print(u[1:10])
print(v[1:10])

result   = gurobi ( list(A=A,obj=c,modelsense="max",rhs=d,sense="="), params=list(OutputFlag=0) )
if (result$status=="OPTIMAL") {
  pi = matrix(result$x,nrow=N)
  u = result$pi[1:N]
  v = result$pi[(N+1):(N+M)]
  val = result$objval
} else {stop("optimization problem with Gurobi.") }
