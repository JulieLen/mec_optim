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

using Gurobi, RData, NLopt, SpecialFunctions, MathProgBase, JuMP, NLsolve


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


Φ = - abs.( kron(Xs, fill(1, nbX)') - kron(Xs, fill(1, nbX)')' ) / 20

# test


Phi = Φ
n = then
m = them

#= Bilan tests

eval_f = OK



=#

# edgeGradient

function edgeGradient(Phi, n, m, xtol_rel = 1e-8, ftol_rel = 1e-15)
    nbX = length(n)
    nbY = length(m)

    function eval_f2(theU)
        theU = reshape(theU, (25, 25))
        theV = Phi - theU
        denomG = 1 + sum(exp.(theU), 2)
        denomH = 1 + sum(exp.(theV), 1)
        valG = sum(n .* log.(denomG))
        valH = sum(m .* log.(denomH'))
        gradG = exp.(theU) .* (n ./ denomG)
        gradH = ( exp.(theV)' .* (m ./ denomH') )'

        return Dict("obj" => valG + valH, "gradient" => vec(gradG - gradH))
    end

    U_init = vec(Phi / 2)

    opt = Opt(:LD_LBFGS,nbX*nbY)
    min_objective!(opt, eval_f2)
    xtol_rel!(opt, xtol_rel)
    ftol_rel!(opt, ftol_rel)

    (optf,optx,ret) = NLopt.optimize(opt, U_init)



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

    theU = reshape(theU, (nbX, nbY))
    theV = Phi - theU
    denomG = 1 + sum(exp.(theU), 2)
    denomH = 1 + sum(exp.(theV), 1)
    gradG = exp.(theU) .* (n ./ denomG)
    gradH = ( exp.(theV)' .* (m ./ denomH') )'

    hessG = hessH = fill(0.0, (nbX*nbY, nbX*nbY))

    for x in 1:nbX, y in 1:nbY, yprime in 1:nbY
        hessG[(x + nbX*(y-1)), (x + nbX*(yprime-1))] = ifelse(y == yprime, gradG[x,y] * (1-gradG[x,y] / n[x]), -gradG[x,y] * gradG[x,yprime] / n[x] )
        hessH[(x-1)*nbY + y, (x-1)*nbY + yprime] = ifelse(y == yprime, gradH[x,y] * (1-gradH[x,y] / n[x]), -gradH[x,y] * gradH[x,yprime] / n[x] )
    end

        return hessG + hessH
end

U_init = Phi / 2

sol = nlsolve(Z, JZ, U_init, xtol = xtol, method = :newton;inplace = false)

Usol = sol.zero

mu = exp.(Usol) .* n ./ (1 + sum(exp.(Usol), 2))
mux0 = n - sum(mu, 2)
mu0y = m - sum(mu, 1)'

val = sum(mu .* Phi) - 2 * sum(mu .* log.(mu ./ sqrt.(n * m'))) - sum(mux0 .* log.(mux0 ./ n)) - sum(mu0y.*log.(mu0y./m))


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

end

# nodal Gradient

function nodalGradient(Phi,n,m, xtol_rel = 1e-8 ,ftol_rel=1e-15)
    K = exp.(Phi / 2)
    tK = K'
    nbX = length(n)
    nbY = length(m)
    function eval_f(ab, n, nbX, nbY)
        a = ab[1:nbX]
        b = ab[(1+nbX):(nbX+nbY)]
        A = exp.(-a / 2)
        B = exp.(-b / 2)
        val = sum(n*a) + sum(m*b) + 2 * A * K * B + sum(A^2) + sum(B^2)
        grada = n - A .* (K * B) - A^2
        gradb = n - B .* (tK * A) - B^2
        grad = vcat(grada, gradb)
        return Dict("objective" => val, "gradient" = grad)
    end

    ab_init = -vcat(log.(n/2) + log.(m/2)) #

    resopt = nlopt()#

    absol = resopt.sol #
    a = absol[1:nbX]
    b = absol[(1+nbX):(nbX+nbY)]
    A = exp.(-a / 2)
    B = exp.(-b / 2)
    mu = vec(A) * (vec(B)*tK)'
    mux0 = n - sum(mu, 2)
    mu0y = m - sum(mu, 1)

    val = sum(mu * Phi) - 2 * sum(mu .* log.(mu ./ sqrt(n))) - sum(mux0 * log.(mux0 ./ n)) - sum(mu0y*log.(mu0y./m))

    return Dict("mu" => mu, "mux0" => mux0, "mu0y" => mu0y, "val" => val, iter = ###)

end

# nodal Newton

function nodalNewton(Phi, n, m, xtol = 1e-8)
    K = exp.(Phi / 2)
    tK = K'
    nbX = length(n)
    nbY = length(m)

    function Z(ab)
        a = ab[1:nbX]
        b = ab[(1 + nbX):(nbX + nbY)]
        A = exp.(-a /2)
        B = exp.(-b /2)
        A2 = A.^2
        B2 = B.^2
        sumx = A * K .* B
        sumy = B * tK .* A
        grada = n - sumx - A2
        gradb = m - sumy - B2
        grad = vcat(grada, gradb)
        return grad
    end

    function JZ(ab)
        J11 = diagm(0.5 * sumx + A2)
        J22 = diagm(0.5 * sumy + B2)
        J12 = 0.5 * A + (B * tK)'
        J21 = J12'
        J = vcat(hcat(J11, J12), hcat(J21, J22))
        return J
    end

    ab_init = vcat(log.(n/2), log.(m/2))

    sol = 1

    absol = 2

    a = absol[1:nbX]
    b = absol[(1+nbX):(1+nbY)]
    A = exp.(-a/2)
    B = exp.(-b/2)
    mu = vec(A) * (vec(B)*tK)'
    mux0 = n - sum(mu, 2)
    mu0y = m - sum(mu, 1)

    val = sum(mu * Phi) - 2 * sum(mu .* log.(mu ./ sqrt(n))) - sum(mux0 * log.(mux0 ./ n)) - sum(mu0y*log.(mu0y./m))

    return Dict("mu" => mu, "mux0" => mux0, "mu0y" => mu0y, "val" => val, iter = 1)

end

# IPFP

function ipfp(Phi, n, m, tol = 1e-6)

    K = exp.(Phi/2)
    tK = K'
    B = sqrt.(m)
    cont = true
    iter = 0
    while cont
        iter = iter + 1
        KBover2 = K * B / 2
        A = sqrt.(n .+ KBover2.^2) - KBover2
        tKAover2 = tK * A / 2
        B = sqrt.(m + tKAover2.^2) - tKAover2

        discrepancy = maximum( abs.(A .* (K * B + A) - n ) ./ n)

        if discrepancy < tol
            cont = false
        end
    end

    mu = A .* (B .* tK)'
    mux0 = n - sum(mu, 2)
    mu0y = m - sum(mu, 1)'
    val = sum(mu .* Phi) - 2 * sum(mu .* log.(mu ./ sqrt.(n * m'))) - sum(mux0 .* log.(mux0 ./ n)) - sum(mu0y.*log.(mu0y./m))

    return Dict("mu" => mu, "mux0" => mux0, "mu0y" => mu0y, "val" => val, "iter" => iter)

end


# Results



function printStats(n, m, mu, Phi, lambda)
    avgAbsDiff = - sum( mu * Phi) / sum(mu) # average absolute age difference between matched partners
    fractionMarried = 2 * sum(mu) / (sum(n) + sum(m)) # fraction of married individuals
    println("Average absolute age difference between matched partners = $avgAbsDiff")
    println("Fraction of married individuals = $fractionMarried")
end

λ = 1

res_edgeGradient = edgeGradient(λ * Φ, then, them)
res_edgeNewton = edgeNewton(λ * Φ, then, them)
res_nodalGradient = nodalGradient(λ * Φ, then, them)
res_nodalNewton = nodalNewton(λ * Φ, then, them)
res_ipfp = ipfp(λ * Φ, then, them)
res_simulatedLinprogr = simulatedLinprogr(λ * Φ, then, them)

printStats(then, them, res_ipfp["mu"], Φ, λ)

println()
prinln("###################")
println()
println("Edge gradient = $(res_edgeGradient["val"])")
println()
println("Edge Newton = $(res_edgeNewton["val"])")
println()
println("Nodal gradient = $(res_nodalGradient["val"])")
println()
println("Nodal Newton = $(res_nodalNewton["val"])")
println()
println("IPFP = $(res_ipfp["val"])")
println()
println("Linear progr = $(res_simulatedLinprogr["val"])")
