#################################################
###             Quantile Methods              ###
#################################################
##########       Alfred Galichon       ##########
#################################################
# References:
# R. Koenker, G. Bassett (1978): Regression Quantiles. Econometrica
#
# R. Koenker (2005): Quantile Regression. Cambridge.
#
# G. Carlier, V. Chernozhukov, A. Galichon (2016):
# "Vector Quantile Regression: an optimal transport approach". Annals of Statistics.



# Data

using Gurobi, CSV, DataFrames, Plots, QuantileRegression, MathProgBase, JuMP

#Pkg.clone("https://github.com/pkofod/QuantileRegression.jl")

path = pwd()
file = "EngelDataset.csv"

engeldataset = CSV.read(string(path, /, file))

t = 0.5
n = size(engeldataset)[1]
ncoldata = size(engeldataset)[2]

X0 = engeldataset[:, ncoldata]

Y = engeldataset[: ,1]

scatter(X0, Y, m=(0.5, [:cross], 2), smooth = true, legend=:none)

# Using the quantile regression function

data = DataFrame(food_expenditure = Y, total_income = X0)

ResultQR = qreg(@formula(food_expenditure ~ total_income), data, t)

# Do it ourselves

X = [fill(1, nrow(data)) X0]
k = ncol(data)

obj = vcat( fill(t, n), fill(1-t, n), fill(0, k) )
A = hcat(sparse(1:n,1:n,1), - sparse(1:n, 1:n, 1), X)
rhs = Y
lb = vcat(fill(0, 2*n), fill(-Inf, k))
ub = Inf

sol = linprog(obj, A, '=', rhs, lb, ub, GurobiSolver(Presolve=0))

beta = sol.sol[2*n+1:2*n+k]

println("################")
println()
println("Results of quantile regression (integrated reg):")
println(ResultQR)
println()
println("################")
println()
println("Results of quantile regression (manual reg):")
println(beta)
println()
println("################")

# Vector Quantile Regression

function VQRTp(X, Y, U, mu, nu)
    n = size(Y)[1]
    d = 1
    r = size(X)[2]
    m = size(U)[1]

    if (n != size(X)[1]) | (d != 1)
        warning("Wrong dimensions.")
    end

    xbar = nu' * X
    c = -vec((U * Y'))

    A1 = kron(sparse(1:n, 1:n, true), fill(1, m)')
    A2 = kron(X', sparse(1:m, 1:m, true))

    f1 = nu
    f2 = vec(mu * xbar)

    e = fill(1, m*n)

    A = vcat(A1, A2)
    f = vcat(f1, f2)

    π_init = vec(mu * nu')

    # Solving
    M = Model(solver = GurobiSolver(Presolve=0))
    @variable(M, x[i=1:2585], start=π_init[i])
    @objective(M, Min, c'*x)
    @constraint(M, x .<= e)
    @constraint(M, cons, A*x .== f )
    status = solve(M)

    if status == :Optimal
        pivec = getvalue(x)
        Lvec = getdual(cons)
    else
        warning("Optimization problem with Gurobi.")
    end

  #############################################

  π = reshape(pivec, m, :)
  L1 = Lvec[1:n]
  L2vec = Lvec[n+1:n+m*r]
  L2 = reshape(L2vec, 11, :)

  ψ = -L1'
  b = -L2
  val = U' * π * Y

  println(π)
  println(ψ)
  println(val)

  return Dict("π" => π, "ψ" => ψ, "b" => b, "val" => val)

end

function ComputeBeta1D(mu, b)
    m = size(mu)[1]
    D = diagm(fill(1, m)) + diagm(fill(-1, m-1), -1)
    beta = diagm(1./mu) * D * b
    return beta
end

# Application

nu	= fill(1/n,n)
step = 0.1
U = round.(collect(linspace(0, 1, 11)), 1)
m = length(U)
mu = fill(1/m, m)
d = 1

sols = VQRTp(X, Y, U, mu, nu)

π = sols["π"]
val = sols["val"]
ψ = sols["ψ"]
b = sols["b"]

betasVQR = ComputeBeta1D(mu, b)

thebetaVQR = betasVQR[6,:]
