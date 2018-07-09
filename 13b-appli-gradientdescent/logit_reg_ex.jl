#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Example of logit regression

=#

#Pkg.add("Optim")

# Packages used

using Distributions, Optim, RCall

# Code

srand(2018)

function sigm(x)
    1 ./ (1 + exp.(-x))
end

n_dim = 5
n_samp = 1000

X = randn((n_samp, n_dim))

θ_0 = rand(n_dim)

μ = sigm( X * θ_0)

y = []

for i in 1:n_samp
    push!(y, rand(Binomial(1, μ[i])))
end

function F(θ)

    μ = sigm( X * θ)
    val = y .* log.(μ) + (1-y) .* log.(1 - μ)

    return -sum(val)

end


function Grad!(A, θ)
    μ = sigm( X * θ)
    A[1:5] = X' * (μ - y)
end


x0 = fill(1.0, n_dim)

# without gradient : optimize(b -> fn(b, y, X), x0, BFGS())

sol = optimize(F, Grad!, x0, BFGS())


println("#### Solution #######")
println()
println(sol)
println()
println()

# Testing that the result is the same in R and Julia

reval("n_samp = 1000
n_dim = 5
g <- function(x){
return( 1 / (1 + exp(-x)))
}
gendata <- function(n_dim, n_samp){
	set.seed(2018)
	X <- matrix(rnorm(n_dim*n_samp),ncol=n_dim)
	theta <- matrix(runif(n_dim,1,4))
	mu = g(X%*%theta)
	y <- numeric(n_samp)
	for (i in 1:n_samp){y[i] <- rbinom(1,1,mu[i])}
	return(list(X, y))
}
a = gendata(n_dim, n_samp)");

x_R = rcopy(reval(("a[[1]]")))
y_R = rcopy(reval("a[[2]]"))


opt = optimize(b -> fn(b, y_R, x_R), ones(n_dim), BFGS())

println("#### Testing R-Julia matching #######")
println()
println("Value of the minimum - R result= 232.3052")
println()
println("Value of the minimum - Julia=  Optim.minimum(opt)")
println()
println("###################")
