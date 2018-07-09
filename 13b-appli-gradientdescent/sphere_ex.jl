#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Example of sphere

=#

# Packages used

using Optim

function fn(x)
    return sum(x.^2)
end

function g!(G, x)
    G[1:n] = 2 * x
end

n = 10

x0 = fill(2.0, n)

sol = optimize(fn, g!, x0, BFGS())


println("#### Solution #######")
println()
println(sol)
println()
println()
