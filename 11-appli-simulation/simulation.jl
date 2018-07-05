#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Probit Simulation

=#

using JuMP, Gurobi, CSV

nbDraws = Int(1e4)
U_y = [1.6, 3.2, 1.1, 0.0]
nbY = length(U_y)

ρ = 0.5

Covar = ρ * fill(1, (nbY, nbY)) + (1 - ρ) * diagm(fill(1,nbY))

E = eig(Covar)
V = E[1]
Q = E[2]

SqrtCovar = Q * diagm(sqrt.(V)) * Q'

ϵ_iy = reshape(randn(nbDraws*nbY), (:,4)) * SqrtCovar

u_iy = (ϵ_iy' .+ U_y)'

u_i = maximum(u_iy, 2)
s_y = fill(1/nbDraws, 4)

for i in 1:4
    s_y[i] = s_y[i] * length( find( (u_iy .- u_i)[:,i] .== 0))
end

A1 = kron(fill(true, nbY)', sparse(1:nbDraws, 1:nbDraws, true))
A2 = kron(sparse(1:nbY, 1:nbY, true), fill(true, nbDraws)')

A = vcat(A1, A2)

rhs = vcat(fill(1/nbDraws, nbDraws), s_y)

obj = vec(ϵ_iy)

m = Model(solver=GurobiSolver(Presolve=0))
@variable(m, x[1:nbY*nbDraws])
@objective(m, Max, obj' * x )
@constraint(m, cons, A * x .== rhs)
status = solve(m)

π = getdual(cons)

Uhat_y = - π[(1+nbDraws):(nbY+nbDraws)] + π[nby+nbDraws]

println("U_y (true vs recovered): ")
println(U_y)
println(Uhat_y)
