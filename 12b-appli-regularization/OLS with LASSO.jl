#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        OLS with LASSO

=#

# Packages use

# Pkg.add("GLMNEt")

using Plots, GLMNet, GLM

# Code

nObs = 1000   # num of observations
nFeature = 100
nonZero  = 10 # num of non-zero components

x = reshape( (2*(rand(nObs*nFeature)) - 1), (nObs,nFeature))

beta = vcat(1:nonZero, fill(0, nFeature-nonZero))

eps = 0.5 * randn(nObs)

y = x * beta + eps

scatter(y, m=(0.5, [:cross], 2), legend=:none)

FIT = glmnet(x, y, alpha = 1, intercept = false)

# plot(FIT)


cvfit = glmnetcv(x, y, alpha = 1, intercept = false)

# cvfit.lambda

# coef(cvfit, s = "lambda.1se")

OLS = lm(x[:,1:10], y)

println(OLS)
