#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Bus Maintenance

=#

# Packages used

using Gurobi, JuMP

# Data generating process

nbX = 10
nbT = 40
nbY = 2 # choice set: 1 = run as usual; 2 = overhaul

IdX = sparse(1:nbX, 1:nbX, 1.0)
LX = sparse([nbX; collect(1:(nbX-1))], 1:nbX, 1.0)
RX = sparse(1:nbX, fill(1.0, nbX), 1.0, nbX, nbX)

P = kron([1, 0], 0.75 * IdX .+ 0.25 * LX ) + kron([0, 1], RX)

IdT = sparse(1:nbT, 1:nbT, 1.0)
NT = sparse(1:(nbT-1), 2:nbT, 1.0, nbT, nbT)
A = kron( kron(IdT, fill(1.0,nbY)), IdX ) - kron(NT, P)

overhaulCost = 8e3

function mainCost(x)
    x * 5e2
end

beta = 0.9

n1_x = fill(1.0,nbX)

b_xt = [n1_x; fill(0.0,nbX*(nbT-1))]
u_xy = [-mainCost(collect(1:(nbX-1))) ; fill(-overhaulCost,nbX+1)]

u_xyt = vec(kron(beta .^(1:nbT),u_xy))

# Solving the problem with Gurobi

env = Gurobi.Env()

setparam!(env, "Presolve", 0) # set presolve to 0
model = Gurobi.Model(env, "opt", :minimize)

add_cvars!(model, b_xt, -Inf, Inf)  # x: x >= 45 # y: y >= 5
update_model!(model)
add_constrs!(model, A, '>', u_xyt) # 50 x + 24 y <= 2400
update_model!(model)

println(model)
optimize(model)

sol = get_solution(model)

U_x_t_gurobi = reshape(sol, (nbX, nbT))

# Solving the problem with backward induction

U_x_t = fill(0.0,nbX,nbT)

X = reshape(u_xyt, (nbX, nbY, nbT))

contVals = maximum(X[:,:,nbT], 2)

U_x_t[:,nbT] = contVals

for t in (nbT-1):-1:1

  myopic = X[:,:,t]
  Econtvals = reshape(P * contVals, (nbX, nbY))

  contVals = maximum(myopic + Econtvals, 2)
  U_x_t[:,t] = contVals

end

# Comparing solutions

println()
println("##################")
println()
println("Solution (Gurobi):")
println(U_x_t_gurobi[:,1])
println()
println("Solution (Backward induction):")
println(U_x_t[:,1])
println()
println("##################")
println()
