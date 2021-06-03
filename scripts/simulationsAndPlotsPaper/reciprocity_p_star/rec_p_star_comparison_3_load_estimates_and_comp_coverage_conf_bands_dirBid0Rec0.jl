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
using ScoreDrivenERGM.DynNets
using Statistics
using Logging
using SharedArrays
using LinearAlgebra
using ProjUtilities

#endregion



@elapsed dfEst = collect_results( datadir("sims", "dgp&FIl_est")) 
dfEst["modelTag"] = string.(dfEst["model"]) 

@elapsed df = dfEst

# @elapsed dfConf = collect_results( datadir("sims", "dgp&FIl_conf")) 
# dfConf["modelTag"] = string.(dfConf["model"]) 

# @elapsed df = antijoin(dfEst, dfConf, on = [:allParDgpT, :modelTag])

df.dgpSettings

begin

# limitSample = 10
limitSample = nothing
# listParUncMethod = ["PB-MVN"] #"WHITE-MLE"
winsorProp=0
listParUncMethod = ["WHITE-MLE"]
#listParUncMethod = ["PB-MVN"] #"WHITE-MLE"
# m = listParUncMethod[1]
# estResRow = collect(eachrow(df))[1]
for estResRow in eachrow(df)

    processFlag = false
    if (estResRow.N in [100]) & (estResRow.T in [300, 3000])  &  (estResRow.model.scoreScalingType == "FISH_D") & (estResRow.dgpSettings.type == "AR")

        if  estResRow.dgpSettings.opt.B[1] == 1
            processFlag = true
            if haskey(estResRow.model.options, "integrated")
                if estResRow.model.options["integrated"]
                end
            end
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
# vEstSdUnParBootDist = mapslices(x -> DynNets.unrestrict_all_par(res.model, x), res.allvEstSdResPar, dims=1)

# # @show res.allvEstSdResPar
# covParBoot = cov(Utilities.drop_bad_un_estimates(vEstSdUnParBootDist)')

# covParBoot, minEigenVal = Utilities.make_pos_def(covParBoot)

# mvSDUnParEstCov = Symmetric(covParBoot)
