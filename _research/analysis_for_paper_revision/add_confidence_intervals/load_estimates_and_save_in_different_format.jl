
#region import and models
using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using DataFrames
using PyPlot
pygui(true)
using JLD2
using JLD
using BSON
using ScoreDrivenERGM
using Statistics
using Logging
using SharedArrays
#endregion



loadFolder = "samDgpFiltSD_conf"

@time loadedDf = collect_results( datadir("sims", loadFolder)) 
a=1
#loadedDf["modelTag"] = string.(loadedDf["model"]) 


begin

for loadedRow in eachrow(loadedDf) 

    dict =DrWatson.tostringdict(loadedRow)# merge(, coverDict)

    loadPath = loadedRow.path

    saveName = loadPath[findlast("\\", loadPath)[end]+1:end-3] * "jld2"

    timeSave = @elapsed DrWatson.save( datadir("sims", loadFolder * "_2", saveName), dict)

end

end


