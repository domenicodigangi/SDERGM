"""
Simulations to estimate coverage of confidence bands with Blasques methods and Buccheri's enhancement
"""

#region import and models
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

using JLD

model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
indTvPar = trues(2)

#endregion

#region quick checks 
begin
T=00
N=50
quantilesVals = [0.975, 0.95, 0.05, 0.025]
parDgpT = DynNets.dgp_misspecified(model_mle, "AR", N, T;  minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=true, B =0.99, sigma = 0.005)
# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
A_T_dgp = sample_dgp(model_mle, parDgpT,N)
res_mle = filter_and_conf_bands(model_mle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

# res_pmle = filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

end


#endregion

#region coverage simulations
model = model_pmle
N=50
T=200
quantilesVals = [0.975, 0.95, 0.05, 0.025]
parDgpT = DynNets.dgp_misspecified(model_mle, "sin", N, T;  minValAlpha = 0.1, maxValAlpha = 0.2, nCycles=1.5, phaseAlpha = 0.5π, phaseshift = 0.1, plotFlag=false)

allCover, allvEstSdResPar, allfVecT_filt, allConfBands, fractErr = conf_bands_coverage(model, parDgpT, N; nSampleCoverage=4)


nVals = [ 50, 100, 200, 300]
tVals = [100, 200, 300]
models = [model_mle, model_pmle]
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)
nSampleCoverage=50

allCoverVarN = Array{typeof(allCover),3}(undef, nNVals, nTVals, nModels)
allvEstSdResParVarN = Array{typeof(allvEstSdResPar),3}(undef, nNVals, nTVals, nModels)
allfVecT_filtVarN = Array{typeof(allfVecT_filt),3}(undef, nNVals, nTVals, nModels)
allConfBandsVarN = Array{typeof(allConfBands),3}(undef, nNVals, nTVals, nModels)
fractErrVarN = Array{typeof(fractErr),3}(undef, nNVals, nTVals, nModels)

parDgpTvarN = Array{Array{Float64,2},2}(undef, nNVals, nTVals)

for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        parDgpT = DynNets.dgp_misspecified(model, "sin", N, T;  minValAlpha = 0.1, maxValAlpha = 0.2, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false)

        parDgpTvarN[indN, indT] = parDgpT
        for (indM, model) in Iterators.enumerate(models)

            filter_and_conf_bands(model, sample_dgp(model, parDgpT,N), quantilesVals; plotFlag =true, parDgpT = parDgpT)
            
            @elapsed allCover, allvEstSdResPar, allfVecT_filt, allConfBands, fractErr = conf_bands_coverage(model, parDgpT, N; nSampleCoverage=nSampleCoverage)

            allCoverVarN[indN, indT, indM] = allCover
            allvEstSdResParVarN[indN, indT, indM] = allvEstSdResPar
            allfVecT_filtVarN[indN, indT, indM] = allfVecT_filt
            allConfBandsVarN[indN, indT, indM] = allConfBands
            fractErrVarN[indN, indT, indM] = fractErr
            # for k=1:nSampleCoverage
            # for b in 1:nBands
            #     for p in 1:nErgmPar 
            #         for t in 1:T
            #             ub = allConfBandsVarN[indN, indT, indM][p, t, b, k] 
            #             lb = allConfBandsVarN[indN, indT, indM][p, t, end-b+1, k] 
            #             ub<lb ? error("wrong bands ordering") : ()
            #             isCovered = lb <= parDgpT[p, t] <= ub 
            #             allCoverVarN[indN, indT, indM][p, t, b, k] = isCovered
            #         end
            #     end
            # end
            # end
             
        end
    end
end
nErgmPar=2

@save("./data/confBands$(nVals)_$(tVals).jld", allCoverVarN, allfVecT_filtVarN,allConfBandsVarN, fractErrVarN, parDgpTvarN)


cover90 =zeros(2,nNVals, nTVals,2)
cover95 =zeros(2,nNVals, nTVals,2)
for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, model) in Iterators.enumerate(models)
            for indPar in 1:2
                cover95[indPar, indN, indT, indM] = mean(allCoverVarN[indN, indT, indM][indPar,:,1,:]) 
                cover90[indPar, indN, indT, indM] = mean(allCoverVarN[indN, indT, indM][indPar,:,2,:]) 
            end
        end
    end
end

begin
figure()
indT = 3
plot(nVals, (cover95[1,:,indT,1]), "--*b")
plot(nVals, (cover95[2,:,indT,1]),"-*b")
plot(nVals, (cover95[1,:,indT,1].+cover95[2,:,indT,1])./2, ":*b")
plot(nVals, (cover95[1,:,indT,2]), "--*r")
plot(nVals, (cover95[2,:,indT,2]),"-*r")
plot(nVals, (cover95[1,:,indT,2].+cover95[2,:,indT,2])./2, ":*r")
title("Coverage Conf Bands 95% Dgp Sin, T=$(tVals[indT])")
legend(reduce(vcat,["$(mod)  $par" for par in ["θ" "η" "θ and η"], mod in ["MLE", "PMLE"]]  ))

end

#endregion
