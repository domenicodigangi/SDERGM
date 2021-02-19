"""
Load and plot simulations results on confidence bands' coverages
"""

#region import and models
using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using DataFrames
using PyPlot
pygui(true)
using JLD2
using ScoreDrivenERGM
using Statistics
using Logging
using SharedArrays
#endregion



df = collect_results!( datadir("sims", "sampleDgpFilterSD_est")) 
df["modelTag"] = string.(df["model"]) 


function average_coverages(res; limitSample=nothing, quantilesVals = [[0.975, 0.025]])
    
    N = res.N
    T = res.T
    nErgmPar = DynNets.number_ergm_par(res.model)
    nBands = length(quantilesVals)

    if isnothing(limitSample)
        nSample = res.nSample
    else 
        nSample = limitSample
    end     

    avgCover = SharedArray(zeros(nErgmPar, 2, nSample))
    allmvSDUnParEstCovWhite = SharedArray(zeros(3*nErgmPar, 3*nErgmPar,nSample))
    allConfBandsPar = SharedArray(zeros(nErgmPar, T, nBands, 2, nSample))
    allConfBandsFiltPar = SharedArray(zeros(nErgmPar, T, nBands, 2, nSample))
    constInds = SharedArray{Bool}((2, nSample))
    errInds = SharedArray{Bool}((2, nSample))

    errInds .= true

    count = SharedArray(ones(1))
    
    Threads.@threads for n=1:nSample
        
        Logging.@info("Estimating Conf Bands  N = $N , T=$T iter n $(count[1])")
        count[1] += 1
        try
            ~, allConfBandsFiltPar[:,:,:,:,n], allConfBandsPar[:,:,:,:,n], errFlag, allmvSDUnParEstCovWhite[:,:,n], distribFilteredSD = DynNets.conf_bands_given_SD_estimates(res.model, res.allObsT[n], N, res.allvEstSdResPar[:,n], quantilesVals;  parUncMethod = "WHITE-MLE")

            coverFiltParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,:,n],   allConfBandsFiltPar[:,:,:,:,n])

            coverParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,:,n],   allConfBandsPar[:,:,:,:,n])

            constInds[:, n]  .= any(res.allvEstSdResPar[3:3:end,n].<=0.02, dims=1)
            for indPar in 1:2
                avgCover[indPar, 1, n] =  mean(coverFiltParUnc[indPar,:,1]) 

                avgCover[indPar, 2, n] =  mean(coverParUnc[indPar,:,1]) 
            end
                errInds[:, n]  .= errFlag
        catch
            Logging.@warn("Error in estimating confidence bands $T, $N, $n")
        end

    end
    return sdata(avgCover), sdata(constInds), sdata(errInds), sdata(allConfBandsFiltPar), sdata(allConfBandsPar), sdata(allmvSDUnParEstCovWhite)
end

begin

for estResRow in eachrow(df)

    avgCover, constInds, errInds, allConfBandsFiltPar, allConfBandsPar, allmvSDUnParEstCovWhite = average_coverages(estResRow)
                    
    coverDict =  (;avgCover, constInds, errInds, allConfBandsFiltPar,  allConfBandsPar, allmvSDUnParEstCovWhite) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
    
    coverDict = merge(DrWatson.tostringdict(estResRow), coverDict)

    loadPath = estResRow.path

    saveName = loadPath[findlast("\\", loadPath)[end]+1:end]

    timeSave = @elapsed @tagsave( datadir("sims", "sampleDgpFilterSD_est_confBands", saveName), coverDict)

    Logging.@info("time save = $timeSave, errRatio = $(mean(errInds)) ")

end

end


