#=
    Math + Econ + Code
        Julie Lenoir
        Alfred Galichon

        Rust

=#

# Packages used

using CSV

# Code

thepath = pwd()

n = 90 # number of discretization points
omax = 450000 # odometer value

fileArr = ["g870.asc", "rt50.asc" ,"t8h203.asc",  "a530875.asc", "a530874.asc","a452374.asc", "a530872.asc", "a452372.asc"]
nbRowsArr = [36,60,81,128,137,137,137,137]
nbColsArr  = [15, 4,48,37,12,10,18,18]

toInclude = 1:8
fileArr = fileArr[toInclude]
nbRowsArr = nbRowsArr[toInclude]
nbColsArr = nbColsArr[toInclude]

#= note: it looks like Rust does not include the d309 model,
which is why he has 162 buses in his dataset, not 166 =#

nbBuses = sum(nbColsArr)
nbMonths = maximum(nbRowsArr) - 11 # number of months in the period
data = fill(0, (nbMonths,nbBuses,3))
inmo = 12 # initial month
inye = 74 # initial year

# A452372.ASC 137x18 matrix for GMC model A4523 buses, model year 1972
# A530872.ASC 137x18 matrix for GMC model A5308 buses, model year 1972
# A530875.ASC 128x37 matrix for GMC model A5308 buses, model year 1975
# G870.ASC     36x15 matrix for Grumman model 870 buses
# T8H203.ASC   81x48 matrix for GMC model T8H203 buses
# A452374.ASC 137x10 matrix for GMC model A4523 buses, model year 1974
# A530874.ASC 137x12 matrix for GMC model A5308 buses, model year 1974
# D309.ASC     110x4 matrix for Davidson model 309 buses
# RT50.ASC     60x4  matrix for Chance model RT50 buses

curbus = 0
output = fill(0,(nbBuses,nbMonths,3))
outputdiscr = fill(0,(nbBuses,nbMonths,3))
transitions = fill(0, (n,n))
themax = 0

# TEST

busType = 1

nbRows = nbRowsArr[busType]
nbCols = nbColsArr[busType]

# TEST

for busType in 1:length(fileArr)

  thefile = fileArr[busType]
  nbRows = nbRowsArr[busType]
  nbCols = nbColsArr[busType]
  tmpdata = CSV.read(string(thepath,"/datafiles/", thefile), datarow=1)

  if size(tmpdata)[1] != nbRows*nbCols
      warning("Unexpected size")
  end

  tmpdata = reshape(Array(tmpdata), (nbRows, nbCols))

  println("Group = $busType")
  println("Nb at least one = $(length(find(tmpdata[6,:] .!= 0)))")
  println("Nb no repl = $(length(find(tmpdata[6,:] .== 0)))")

  for busId in 1:nbCols
      curbus = curbus + 1

      mo1stRepl = tmpdata[4,busId]
      ye1stRepl = tmpdata[5,busId]
      odo1stRep = tmpdata[6,busId]

      mo2ndRepl = tmpdata[7,busId]
      ye2ndRepl = tmpdata[8,busId]
      odo2ndRep = tmpdata[9,busId]

      moDataBegins = tmpdata[10,busId]
      yeDataBegins = tmpdata[11,busId]

      odoReadings = tmpdata[12:nbRows,busId]

      wasreplacedonce = ifelse.((tmpdata[12:nbRows,busId] .>= odo1stRep) .& (odo1stRep.>0) , 1, 0)
      wasreplacedtwice = ifelse.((tmpdata[12:nbRows,busId] .>= odo2ndRep) .& (odo2ndRep.>0) , 1, 0)
      howmanytimesreplaced = wasreplacedonce + wasreplacedtwice

      correctedmileage = tmpdata[12:nbRows,busId] + howmanytimesreplaced .* (howmanytimesreplaced-2) .* odo1stRep - 0.5 * howmanytimesreplaced .* (howmanytimesreplaced-1) .* odo1stRep


      output[curbus,1:(nbRows-12),1] =  howmanytimesreplaced[2:(nbRows-11)] - howmanytimesreplaced[1:(nbRows-12)]
      output[curbus,1:(nbRows-12),2] =  correctedmileage[1:(nbRows-12)]
      output[curbus,1:(nbRows-12),3] = tmpdata[13:nbRows,busId] - tmpdata[12:(nbRows-1),busId]

      outputdiscr[curbus,:,1] = output[curbus,:,1]
      outputdiscr[curbus,:,2:3] = ceil.( n * output[curbus,:,2:3] / omax )

      for t in 1:(nbRows - 13)
          i = outputdiscr[curbus,t,2]
          j = outputdiscr[curbus,t+1,2]
          transitions[i,j]  = transitions[i,j] + 1
      end
  end
end
