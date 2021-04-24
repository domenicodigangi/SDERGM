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
using LinearAlgebra
#endregion



@elapsed dfEst = collect_results( datadir("sims", "dgp&FIl_est")) 
dfEst["modelTag"] = string.(dfEst["model"]) 

df = dfEst

# @elapsed dfConf = collect_results( datadir("sims", "dgp&FIl_conf")) 
# dfConf["modelTag"] = string.(dfConf["model"]) 

# @elapsed df = antijoin(dfEst, dfConf, on = [:allParDgpT, :modelTag])


function average_coverages(res, m; limitSample=nothing, quantilesVals = [ [0.975, 0.025]], winsorProp=0)
    
    N = res.N
    T = res.T
    nErgmPar = DynNets.number_ergm_par(res.model)
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
    
    # if m == "PB-MVN"
    #     if res.dgpSettings.type == "SD"
    #         vEstSdUnParBootDist = mapslices(x -> DynNets.unrestrict_all_par(res.model, res.model.indTvPar, x), res.allvEstSdResPar, dims=1)

    #         covParBoot = cov(Utilities.drop_bad_un_estimates(vEstSdUnParBootDist)')

    #         covParBoot, minEigenVal = Utilities.make_pos_def(covParBoot)

    #         mvSDUnParEstCov = Symmetric(covParBoot)

    #     end
    # end

    Threads.@threads for n=1:nSample
        
        Logging.@info("Estimating Conf Bands  N = $N , T=$T, $(DynNets.name(res.model)), $(res.dgpSettings) iter n $(count[1]), ")
        count[1] += 1

        vEstSdUnPar = unrestrict_all_par(res.model, res.model.indTvPar, res.allvEstSdResPar[:,n])

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
    return sdata(avgCover), sdata(constInds), sdata(errInds), sdata(allConfBandsFiltPar), sdata(allConfBandsPar), sdata(allmvSDUnParEstCovWhite)
end
begin

# limitSample = 10
limitSample = nothing
# listParUncMethod = ["PB-MVN"] #"WHITE-MLE"
winsorProp=0.2
listParUncMethod = ["WHITE-MLE"]
#listParUncMethod = ["PB-MVN"] #"WHITE-MLE"
# m = listParUncMethod[1]
# estResRow = collect(eachrow(df))[1]
for estResRow in eachrow(df)

    processFlag = false
    if (estResRow.N in [50, 100, 500]) & (estResRow.T in [300, 600])  &  (estResRow.model.scoreScalingType == "FISH_D") & (estResRow.dgpSettings.type == "AR")

        if  estResRow.dgpSettings.opt.sigma ==  ScoreDrivenERGM.DynNets.list_example_dgp_settings(DynNets.reference_model(estResRow.model)).dgpSetARlowlow.opt.sigma 
            processFlag = true
        end
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

            timeSave = save( datadir("sims", "dgp&FIl_conf", saveName), coverDict)

            Logging.@info(" errRatio = $(mean(errInds)) , time = $timeConf")

        end
    end
end

end





# T = 100
# N = 100
# nSample = 100

# modelTag = string(DynNets.SdErgmDirBin0Rec0_pmle("HESS_D"))

# dgpSetting = DynNets.list_example_dgp_settings(DynNets.SdErgmDirBin0Rec0_mle()).dgpSetSDhigh


# res = filter([:modelTag, :T, :N, :dgpSettings, :nSample] => (m,t,n, d, s) -> all((m==modelTag, t==T, n==N, d == dgpSetting, s == nSample)), df)[1,:]

# res.allvEstSdResPar
# vEstSdUnParBootDist = mapslices(x -> DynNets.unrestrict_all_par(res.model, res.model.indTvPar, x), res.allvEstSdResPar, dims=1)

# # @show res.allvEstSdResPar
# covParBoot = cov(Utilities.drop_bad_un_estimates(vEstSdUnParBootDist)')

# covParBoot, minEigenVal = Utilities.make_pos_def(covParBoot)

# mvSDUnParEstCov = Symmetric(covParBoot)
