


using CSV, Gurobi, JuMP

syntheticData = true
doGurobi = false
doIPFP1 = true
doIPFP2 = true

tol=1e-9
maxiter = 100
sigma = 0.1 #0.001 # note: 0.1 to 0.001


if syntheticData
  #seed = 777
  nbX=10
  nbY=8
  #set.seed(seed)
  Phi=rand(nbX, nbY)
  p=fill(1/nbX,nbX)
  q=fill(1/nbY,nbY)
else
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
  p = fill(1,nobs)
  q = fill(1,nobs)

  nbX = length(p)
  nbY = length(q)
end

if doGurobi

  A1 = kron(fill(1, nbY)', sparse(1:nbX, 1:nbX, 1))
  A2 = kron(sparse(1:nbY, 1:nbY, 1), fill(1, nbX)')

  B = vcat(A1, A2)

  d = vcat(p, q)

  # Solving the model

  @time begin
    m = Model(solver=GurobiSolver(Presolve=0))
    @variable(m, x[1:size(c)[1]])
    @objective(m, Max, c'*x)
    @constraint(m, x*B' .== d)
    @constraint(m, x .>= 0.0)
    status = solve(m)
  end

  nrow = min(8,nbX)
  ncol = min(8,nbY)

  if status == :Optimal
        π = getvalue(x)
        u_gurobi = π[1:nbX]
        v_gurobi = π[(nbX+1):(nbX+nbY)]
        val_gurobi = getobjectivevalue(m)
  else
        warn("Optimization problem with Gurobi. status = $status")
  end

  println("Value of the problem (Gurobi) = $total_distance")
  println()
  println(u_gurobi[1:nrow]-u_gurobi[nrow])
  println(v_gurobi[1:ncol]+u_gurobi[nrow])
  println("***********************")
end

# computation of the  regularized problem with the  IPFP 1

if doIPFP1
  @time begin
    cont = true
    iter = 0
    K = exp.(Phi/sigma)
    B = fill(1, nbY)

    A = p ./ (K * B)
    KA = A' * K

    while cont
    iter = iter + 1
    err = maximum(abs.(KA' .* B ./ q - 1))
      if (err < tol) | (iter >= maxiter)
        cont = false
      end
    B = q ./ KA'
    A = p ./ (K * B)
    KA = A' * K
    end

    u = - sigma * log.(A)
    v = - sigma * log.(B)
  end

  π = K .* A .* B'
  val = sum(π .* Phi) - sigma * sum(π .* log.(π))

  if iter >= maxiter
    println("Maximum number of iterations reached in IPFP1.")
  else
    println("IPFP1 converged in $iter steps")
    println("Value of the problem (IPFP1) = $val")
    println("Sum(π.*Phi) (IPFP1) = $(sum(π.*Phi))")
    println()
    println(u-u[nbX])
    println()
    println(v+u[nbX])
  end

  println("####################")

end

# computation of the  regularized problem with the  IPFP 1(bis)

if true
  @time begin
    iter = 0
    cont = true
    v = fill(0,nbY)
    mu = - sigma * log.(p)
    nu = - sigma * log.(q)

    while cont
      iter = iter + 1

      u = mu + sigma * log.( sum( exp.( (Phi .- v') / sigma ), 2))
      KA = sum( exp.( (Phi .- u) / sigma ), 1)

      err = maximum( abs.( KA' .* exp.(-v / sigma) ./ q - 1))

      if (err < tol) | (iter >= maxiter)
        cont = false
      end

      v = nu + (sigma * log.(KA))'

    end

    π = exp.( ( Phi .- u .- v' ) / sigma )
    val = sum(π .* Phi) - sigma * sum( (π .* log.(π))[find(!isnull, π)] )

  end
  println()
  if iter >= maxiter
    println("Maximum number of iterations reached in IPFP 1(bis).")
  else
    println("IPFP 1(bis) converged in $iter steps")
    println("Value of the problem (IPFP 1(bis)) = $val")
    println("Sum(π.*Phi) (IPFP 1(bis)) = $(sum(π.*Phi))")
    println(u-u[nbX])
    println(v+u[nbX])
  end

  println("########################")

end

println()

# computation of the  regularized problem with the  IPFP 2

if doIPFP2
  @time begin
    iter = 0
    cont = true
    v = fill(0,nbY)
    mu = - sigma * log.(p)
    nu = - sigma * log.(q)
    uprec = -Inf

    while cont

      iter = iter + 1

      vstar = maximum(Phi .- v', 2)

      u = mu .+ vstar .+ sigma * log.( sum(exp.( (Phi .- v' .- vstar) / sigma) ,2) )

      err = maximum( abs.( u - uprec) )

      uprec = u

      ustar = maximum(Phi .- u, 1)

      v = nu .+ ustar' .+ sigma * log.( sum( exp.( (Phi .- u .- ustar) / sigma ), 1)')

      if (err < tol) | (iter >= maxiter)
        cont = false
      end

    end

    π = exp.( ( Phi .- u .- v') / sigma)
    val = sum( π .* Phi) - sigma * sum( π .* log.(π))
  end

  println()
  if iter >= maxiter
    println("Maximum number of iterations reached in IPFP 2.")
  else
    println("IPFP 2 converged in $iter steps")
    println("Value of the problem (IPFP 2 = $val")
    println("Sum(π.*Phi) (IPFP 2) = $(sum(π.*Phi))")
    println()
    println(u-u[nbX])
    println()
    println(v+u[nbX])
  end

  println("########################")

end
