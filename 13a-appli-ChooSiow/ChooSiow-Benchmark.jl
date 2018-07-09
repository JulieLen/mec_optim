#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Choo-Siow benchmark

=#

# Packages use

# Pkg.add("RData")
# Pkg.add("NLopt")
# Pkg.add("SpecialFunctions")

using Gurobi, RData, NLopt, SpecialFunctions


path = pwd()
file1 = "nSinglesv4.RData"
file2 = "nMarrv4.RData"
file3 = "nAvailv4.RData"

obj1 = load(string(path, "/ChooSiowData/", file1))
obj2 = load(string(path, "/ChooSiowData/", file2))
obj3 = load(string(path, "/ChooSiowData/", file3))

nbCateg = 25

nSingles = obj1["nSingles70n"]
marr = obj2["marr70nN"]
nAvail = obj3["avail70n"]

muhatx0 = nSingles[1:nbCateg, 1]
muhat0y = nSingles[1:nbCateg, 2]

muhatxt = marr[1][1:nbCateg, 1:nbCateg]
then = nAvail[1][1:nbCateg, 1]
them = nAvail[1][1:nbCateg, 2]

nbIndiv = sum(then) + sum(them)

then = then / nbIndiv
them = them / nbIndiv

nbX = length(then)
nbY = length(them)

Xs = collect(1:nbCateg) + 15
Ys = collect(1:nbCateg) + 15




theϕ = - abs.( kron(Xs, fill(1, nbX)') - kron(Xs, fill(1, nbX)')' ) / 20

# test

theU = theϕ / 2
Phi = theϕ
n = then
m = them

#= Bilan tests

eval_f = OK



=#

# edgeGradient

function edgeGradient(Phi, n, m, xtol_rel = 1e-8, ftol_rel = 1e-15)
    nbX = length(n)
    nbY = length(m)

    function eval_f(theU)
        theU = reshape(theU, (nbX, nbY))
        theV = Phi - theU
        denomG = 1 + sum(exp.(theU), 2)
        denomH = 1 + sum(exp.(theV), 1)
        valG = sum(n .* log.(denomG))
        valH = sum(m .* log.(denomH'))
        gradG = exp.(theU) .* (n ./ denomG)
        gradH = ( exp.(theV)' .* (m ./ denomH') )'

        ret = Dict("obj" => valG + valH, "gradient" => vec(gradG - gradH))

        return ret
    end



    resopt = nloptr(x0 = U_init, eval_f = eval_f,
                    opt = list("algorithm" = "NLOPT_LD_LBFGS",
                               "xtol_rel"= xtol_rel,
                               "ftol_rel"= ftol_rel))
    Usol = matrix(resopt$solution,nbX,nbY)


    mu = exp.(Usol) .* (n ./ (1 + sum(exp.(Usol), 2)))
    mux0 = n - sum(mu, 2)
    mu0y = m - sum(mu, 1)

    val = sum(mu .* Phi) - 2 * sum( mu ./ sqrt(n * m) ) - sum(mux0 .* log(mux0/n)) - sum(mu0y .* log(mu0y / m))

    return Dict("mu" => mu, "mux0" => mux0, "mu0y" => mu0y, "val" => val, iter = #resopt$iterations# )

end


#edge Newton

function edgeNewton(Phi, n, m, xtol = 1e-5)
    nbX = length(n)
    nbY = length(m)

    function Z(theU)
        theU = reshape(theU, (nbX, nbY))
        theV = Phi - theU
        denomG = 1 + sum(exp.(theU), 2)
        denomH = 1 + sum(exp.(theV), 1)
        gradG = exp.(theU) .* (n ./ denomG)
        gradH = ( exp.(theV)' .* (m ./ denomH') )'

        return vec(gradG - gradH)
    end

    function JZ(theU)
        hessG = fill(0, (nbX*nbY, nbX*nbY))

        for x in 1:nbx, y in 1:nbY, yprime in 1:nbY
            hessG[x + nbX*(y-1), x + nbX*(yprime-1)] = ifelse(y == yprime, gradG[x,y] * (1-gradG[x,y] / n[x]), -gradG[x,y] * gradG[x,yprime] / n[x] )
            hessH[(x-1)*nbY + y, (x-1)*nbY + yprime] = ifelse(y == yprime, gradH[x,y] * (1-gradH[x,y] / n[x]), -gradH[x,y] * gradH[x,yprime] / n[x] )
        end

        return hessG + hessH
    end

    U_init = Phi / 2

    #=sol = nleqslv(x = U_init,
                  fn = Z, jac = JZ,
                  method = "Broyden", # "Newton"
                  control = list(xtol=xtol))


    Usol = matrix(sol$x,nbX,nbY)=#

    mu = exp.(Usol) .* n ./ (1 + sum(exp.(Usol), 2))
    mux0 = n - sum(mu, 2)
    mu0y = m - sum(mu, 1)

    val = sum(mu .* Phi) - 2 * sum(mu .* log.(mu ./ sqrt(n * m)) ) - sum(mux0 .* log(mux0/n)) - sum(mu0y .* log(mu0y / m))


    return Dict("mu" => mu, "mux0" => mux0, "mu0y" => mu0y, "val" => val, iter = #resopt$iterations# )

end

# simulatedLinprogr

function simulatedLinprogr(Phi, n, m, nbDraws = Int(1e3))
nbX = length(n)
nbY = length(m)
nbI = nbX * nbDraws
nbJ = nbY * nbDraws

ϵ_iy = reshape( digamma(1) - log.(-log.(rand(nbI*nbY))), (nbI,nbY))
ϵ0_i = vec(digamma(1) - log.(-log.(rand(nbI))))

I_ix = fill(0, (nbI, nbX))
for x in 1:nbX
    I_ix[Int(nbDraws * (x-1) + 1):Int(nbDraws * x), x] = 1
end

η_xj = reshape( digamma(1) - log.(-log.(rand(nbX*nbJ))), (nbX,nbJ))
η0_i = vec(digamma(1) - log.(-log.(rand(nbI))))

I_yj = fill(0, (nbY, nbJ))
for y in 1:nbY
    I_yj[y,Int(nbDraws * (y-1) + 1):Int(nbDraws * y)] = 1
end

ni = vec(I_ix * n) / nbDraws
mj = vec(m' * I_yj) / nbDraws

@time begin
    kron(fill(1, nbY), sparse(1:nbI, 1:nbI, 1))
end

@time begin
    kron(fill(1, nbY), sparse(1:nbI, 1:nbI, true))
end

@time begin
    kron(fill(true, nbY), sparse(1:nbI, 1:nbI, true))
end


A11 = kron(fill(1, nbY), sparse(1:nbI, 1:nbI, true))
A12 = sparse(1:nbI*nbY, 1:nbJ, 0)
A13 = kron(sparse(1:nbY, 1:nbY, -1), I_ix)

A1 = vcat(A11, A12, A13)

A21 = sparse(1:nbX*nbJ, 1:nbI, 0)
A22 = kron(sparse(1:nbJ, 1:nbJ, 1), fill(1, nbX))
A23 = kron(I_yj', sparse(1:nbX, 1:nbX, 1))

A2 = vcat(A11, A12, A13)

A = hcat(A1, A2)

lb = vcat(ϵ0_i, η0_j, fill(-Inf, nbX*nbY))
rhs = vcat(ϵ_iy, eta_xj + Phi % Y_yj)
obj = vcat[ni, mj, fill(0, nbX*nbY)]


sol = linprog(obj, A, '=', rhs, lb, GurobiSolver(Presolve=0))

muiy =
mu = I_ix' * muiy

val = sum()

muiy = matrix(result$pi[1:(nbI*nbY)],nrow=nbI)
mu = t(I_ix) %*% muiy
val = sum(ni*result$x[1:nbI]) + sum(mj*result$x[(nbI+1):(nbI+nbJ)])

mux0 = n - sum(mu, 2)
mu0y = m - sum(mu, 1)'

return Dict("mu" => mu, "mux0" => mux0, "mu0y" => mu0y, "val" => val, iter = Nullable(1)))



mux0 = n - apply(mu,1,sum)
mu0y = m - apply(mu,2,sum)
return(list(mu = mu, mux0 = mux0, mu0y = mu0y, val = val, iter = NA))




    #
