"""
Simulations to estimate coverage of confidence bands with Blasques methods and Buccheri's enhancement
"""


using ScoreDrivenExponentialRandomGraphs

import ScoreDrivenExponentialRandomGraphs:StaticNets, DynNets

import ScoreDrivenExponentialRandomGraphs.DynNets:GasNetModel,GasNetModelDirBin0Rec0, sample_dgp, statsFromMat, array2VecGasPar, unrestrict_all_par, conf_bands_par_uncertainty, avg_grad_and_hess_obj_SD_filter_time_seq, conf_bands_par_uncertainty, number_ergm_par, filter_and_conf_bands, conf_bands_coverage, estimate
using ScoreDrivenExponentialRandomGraphs.Utilities

using PyPlot
pygui(true)


using ForwardDiff
using StatsBase
using LinearAlgebra
using Distributions
using Statistics



model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
indTvPar = trues(2)





T=200
quantilesVals = [0.975, 0.95, 0.05, 0.025]
begin
N=150
parDgpT = DynNets.dgp_misspecified(model_mle, "sin", N, T;  minValAlpha = 0.1, maxValAlpha = 0.2, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=true)
# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
A_T_dgp = sample_dgp(model_mle, parDgpT,N)
res_mle = filter_and_conf_bands(model_mle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

res_pmle = filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

end

model = model_mle
allCover, allvEstSdResPar, allfVecT_filt, allConfBands, fractErr = conf_bands_coverage(model, parDgpT, N; nSampleCoverage=2)

model = model_pmle
allCover, allvEstSdResPar, allfVecT_filt, allConfBands, fractErr = conf_bands_coverage(model, parDgpT, N; nSampleCoverage=2)

coverage = mean(allCover)

#endregion

#region tests

T=200
nVals = [25, 50, 75, 100, 150, 200]
models = [model_mle, model_pmle]
nNVals = length(nVals)
nModels = length(models)
nSampleCoverage=10

allCoverVarN = Array{typeof(allCover),2}(undef, nNVals, nModels)
allvEstSdResParVarN = Array{typeof(allvEstSdResPar),2}(undef, nNVals, nModels)
allfVecT_filtVarN = Array{typeof(allfVecT_filt),2}(undef, nNVals, nModels)
allConfBandsVarN = Array{typeof(allConfBands),2}(undef, nNVals, nModels)
fractErrVarN = Array{typeof(fractErr),2}(undef, nNVals, nModels)

T=200
quantilesVals = [0.975, 0.95, 0.05, 0.025]

for (indN, N) in Iterators.enumerate(nVals) 
    parDgpT = DynNets.dgp_misspecified(model, "sin", N, T;  minValAlpha = 0.1, maxValAlpha = 0.2, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=true)

    for (indM, model) in Iterators.enumerate(models)

        filter_and_conf_bands(model, sample_dgp(model, parDgpT,N), quantilesVals; plotFlag =true, parDgpT = parDgpT)
        
        @elapsed allCover, allvEstSdResPar, allfVecT_filt, allConfBands, fractErr = conf_bands_coverage(model, parDgpT, N; nSampleCoverage=nSampleCoverage)

        allCoverVarN[indN, indM] = allCover
        allvEstSdResParVarN[indN, indM] = allvEstSdResPar
        allfVecT_filtVarN[indN, indM] = allfVecT_filt
        allConfBandsVarN[indN, indM] = allConfBands
        fractErrVarN[indN, indM] = fractErr
        
    end
end
#endregion 


for (indM, model) in Iterators.enumerate(models)
    cover95 = [mean(allCover[1,:,:]) for allCover in allCoverVarN[:,indM]]
    cover90 = [mean(allCover[2,:,:]) for allCover in allCoverVarN[:,indM]]

    plot(nVals, cover95)
end
