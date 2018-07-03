#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Possitive Assortative Matching

=#

# Packages use

using Distributions, QuadGK

# Setting useful functions

function cdf_P(x)
    cdf(Uniform(), x)
end

function quant_P(x)
    quantile(Uniform(), x)
end


function cdf_Q(x)
    cdf(Normal(), x)
end

function quant_Q(x)
    quantile(Normal(), x)
end

function Φ(x, y)
    x*y
end

function dΦ_dx(x, y)
    y
end

function dΦ_dy(x, y)
    x
end

function Tx(x)
    quant_Q(cdf_P(x))
end

function Tinvy(y)
    quant_P(cdf_Q(y))
end

function ux(x)
    quadgk(x -> dΦ_dx(x, Tx(x)), 0, x)[1]
end

function vy(x)
    quadgk(x -> dΦ_dy(Tinvy(x), x), Tx(0), x)[1] - Φ(0,Tx(0))
end
