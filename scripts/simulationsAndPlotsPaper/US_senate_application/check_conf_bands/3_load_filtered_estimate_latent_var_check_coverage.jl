#########
#Created Date: Wednesday April 28th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday June 3rd 2021 5:18:32 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:  Load estimates of latent AR dynamics with model and dgp parameters (variance of latent parameters) equal to that estimated in the application to US congress data
#-----
#-----
########


using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
using DataFrames
using PyPlot
pygui(true)
using JLD2
using ScoreDrivenERGM
using ScoreDrivenERGM.DynNets
using Statistics
using Logging
using SharedArrays
using LinearAlgebra
using Logging


@elapsed df = collect_results( datadir("sims", "dgpAR_Fil_SS")) 
df["modelTag"] = string.(df["model"]) 



# begin

# limitSample = 10
limitSample = nothing
winsorProp=0
listParUncMethod = ["WHITE-MLE"]

@elapsed subDf = df[df.T .==73,: ]


@elapsed estResRow = eachrow(subDf)[1]
m = listParUncMethod[1]

for estResRow in eachrow(subDf)

    processFlag = false
    if  (estResRow.T in [73])  &  (estResRow.model.scoreScalingType == "FISH_D") & (estResRow.dgpSettings.type == "AR")
            processFlag = true
    end


    if processFlag
        for m in listParUncMethod

            if winsorProp!=0
                m = "$m$winsorProp"
            end

            timeConf = @elapsed avgCover, constInds, errInds, allConfBandsFiltPar, allConfBandsPar, allmvSDUnParEstCovWhite = average_coverages(estResRow, m; limitSample = limitSample, winsorProp=winsorProp)
                            
            coverDict =  (;m, avgCover, constInds, errInds, allConfBandsFiltPar,  allConfBandsPar, allmvSDUnParEstCovWhite) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
            
            coverDict = merge(DrWatson.tostringdict(estResRow), coverDict)
            delete!(coverDict, "path")
            delete!(coverDict, "modelTag")
            
            coverDict["dgp"] = coverDict["dgpSettings"]
            delete!(coverDict, "dgpSettings")
            coverDict["S"] = coverDict["nSample"]
            delete!(coverDict, "nSample")
            

            # loadPath = estResRow.path

            saveName = replace.( savename(coverDict, "jld2"; allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.SdErgm) ), r"[\"]" => "") #loadPath[findlast("\\", loadPath)[end]+1:end]

            timeSave = save( datadir("sims", "dgp&FIl_conf", "cover", saveName), coverDict)

            Logging.@info(" errRatio = $(mean(errInds)) , time = $timeConf")

        end
    end
end

# end



avgCover[:,1,:]
errInds

plot(allConfBandsFiltPar[2,:,1,:,3])
plot(estResRow.allfVecT_filt[2,:,1])
estResRow.allvEstSdResPar[end,:]
names(estResRow)

avgCover[2,1,.!errInds[2,:]]
constInds

@elapsed res = estResRow


N = res.N
T = res.T
nErgmPar = DynNets.number_ergm_par(res.model)
quantilesVals = [ [0.975, 0.025]]
nBands = length(quantilesVals)

if isnothing(limitSample)
    nSample = res.nSample
elseif res.nSample >= limitSample 
    nSample = limitSample
else 
    error("sample size?")
end     

avgCover = SharedArray(zeros(nErgmPar, 2, nSample))
allmvSDUnParEstCovWhite = SharedArray(zeros(3*nErgmPar, 3*nErgmPar,nSample))
allConfBandsPar = SharedArray(zeros(nErgmPar, T, nBands, 2, nSample))
allConfBandsFiltPar = SharedArray(zeros(nErgmPar, T, nBands, 2, nSample))
constInds = SharedArray{Bool}((2, nSample))
errInds = SharedArray{Bool}((2, nSample))

errInds .= true

count = SharedArray(ones(1))

mvSDUnParEstCov = zeros(3,3)
import ScoreDrivenERGM.DynNets:unrestrict_all_par, conf_bands_given_SD_estimates, conf_bands_coverage


Threads.@threads for n=1:nSample
    
    Logging.@info("Estimating Conf Bands  N = $N , T=$T, $(DynNets.name(res.model)), $(res.dgpSettings) iter n $(count[1]), ")
    count[1] += 1

    vEstSdUnPar = unrestrict_all_par(res.model, res.allvEstSdResPar[:,n])

    ~, allConfBandsFiltPar[:,:,:,:,n], allConfBandsPar[:,:,:,:,n], errFlag, allmvSDUnParEstCovWhite[:,:,n], distribFilteredSD = DynNets.conf_bands_given_SD_estimates(res.model, N, res.allObsT[n], vEstSdUnPar, res.allftot_0[:,n], quantilesVals;  parUncMethod = m, mvSDUnParEstCov=mvSDUnParEstCov, winsorProp=winsorProp )

    coverFiltParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,:,n],   allConfBandsFiltPar[:,:,:,:,n])

    coverParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,:,n],   allConfBandsPar[:,:,:,:,n])

    constInds[:, n]  .= any(res.allvEstSdResPar[3:3:end,n].<=0.02, dims=1)
    for indPar in 1:2
        avgCover[indPar, 1, n] =  mean(coverFiltParUnc[indPar,:,1]) 

        avgCover[indPar, 2, n] =  mean(coverParUnc[indPar,:,1]) 
    end
        errInds[:, n]  .= errFlag

    # catch
    #     Logging.@warn("Error in estimating confidence bands $T, $N, $n")
    # end

end

