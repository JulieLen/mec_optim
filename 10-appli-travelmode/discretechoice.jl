#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Discrete Choice

=#

# Packages use

using CSV, DataFrames

# Data shaping

path = pwd()
filename = "travelmodedata.csv"
travelmodedata = CSV.read(string(path, /, filename))

nobs = size(travelmodedata)[1]
choice = fill(0, nobs)
for i in 1:nobs
    if travelmodedata[:choice][i] == "yes"
        choice[i] = 1
    end
end

mode = fill(1, nobs)

for i in 1:nobs
    if travelmodedata[:mode][i] == "train"
        mode[i] = 2
    elseif travelmodedata[:mode][i] == "bus"
        mode[i] = 3
    elseif travelmodedata[:mode][i] == "car"
        mode[i] = 4
    end
end

travelmodedata[:choice] = choice
travelmodedata[:mode] = mode

nind = nobs / 4

ncols = size(travelmodedata)[2]

travelmodedataset = reshape(Array(travelmodedata), (:, 210, 9))

choices = travelmodedataset[:,:,3]

# Market Shares

s = mean(choices, 2)

println("Market Shares:")
println(DataFrame(names = ["air","train","bus","car"], market_share = s[1:4]))
println()

# Logit model

Ulogit = log.(s[1:4]/s[4])

println("Systematic Utilies: ")
println(Ulogit)
println()

# Nested logit model: two nests = (car, nocar)

lambda = [0.5 0.5]

Unocar = lambda[1] * log.(s[1:3]) + (1 - lambda[1]) * log.(sum(s[1:3]))
Ucar = lambda[2] * log(s[4]) + (1 - lambda[2]) * log(sum(s[4]))
Unested = vcat(Unocar, Ucar) - Ucar

println("Systematic utilities (nested logit):")
println(Unested)
println()

println("Choice probabilities within nocar nest (predicted vs observed):")
println( exp.(Unested[1:3] / lambda[1]) / sum(exp.( Unested[1:3] / lambda[1])))
println( s[1:3] ./ sum( s[1:3] ) )
println()

println("Choice probabilities of car nest (predicted vs observed):")
println( 1 / ( sum( exp.( Unested[1:3] / lambda[1] ))^lambda[1] + 1))
println(s[4])
println()
