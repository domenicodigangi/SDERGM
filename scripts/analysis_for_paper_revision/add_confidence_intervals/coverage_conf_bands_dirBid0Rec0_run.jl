"""
Simulations to estimate coverage of confidence bands with Blasques methods and Buccheri's enhancement
"""

#region import and models
using ScoreDrivenExponentialRandomGraphs

import ScoreDrivenExponentialRandomGraphs:StaticNets, DynNets

import ScoreDrivenExponentialRandomGraphs.DynNets:GasNetModel,GasNetModelDirBin0Rec0, sample_dgp, statsFromMat, array2VecGasPar, unrestrict_all_par, conf_bands_par_uncertainty, avg_grad_and_hess_obj_SD_filter_time_seq, conf_bands_par_uncertainty, number_ergm_par, filter_and_conf_bands, conf_bands_coverage, estimate, mle_distrib_filtered_par, plot_filtered_and_conf_bands
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



# NEXT TO DO:

# LANCIARE STIMA COVERAGE DOPO  CONTROLLO E PER DUE dgp

# Perché la banda di Blasques viene asimmetrica? Dovuto alla non linearitá in B credo

# Posso provare a riapplicare Blasques e Buccheri ai filtri riscalati
# Posso provare a ricostruire la parametric uncertainty usando parametric bootstrap sui parametri del filtro SD


#region coverage simulations
begin 
model = model_pmle
N=50
T=200
dgpType =  "AR"
dgpOptions = (minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.98, sigma = 0.001)
quantilesVals = [[0.975, 0.025]]

parDgpT = DynNets.dgp_misspecified(model_mle, dgpType, N, T;  dgpOptions...)

allCoverBuccheri, allCoverBlasques, allvEstSdResPar, allfVecT_filt, allConfBandsBuccheri,allConfBandsBlasques, allErrFlags = conf_bands_coverage(model, parDgpT, N, 2, quantilesVals)


nVals = [ 50, 100, 200, 300]
tVals = [100, 200, 300]
models = [model_mle, model_pmle]
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)
nSampleCoverage=50

allCoverBuccheriVarN = Array{typeof(allCoverBuccheri),3}(undef, nNVals, nTVals, nModels)
allCoverBlasquesVarN = Array{typeof(allCoverBlasques),3}(undef, nNVals, nTVals, nModels)
allvEstSdResParVarN = Array{typeof(allvEstSdResPar),3}(undef, nNVals, nTVals, nModels)
allfVecT_filtVarN = Array{typeof(allfVecT_filt),3}(undef, nNVals, nTVals, nModels)
allConfBandsBuccheriVarN = Array{typeof(allConfBandsBuccheri),3}(undef, nNVals, nTVals, nModels)
allConfBandsBlasquesVarN = Array{typeof(allConfBandsBlasques),3}(undef, nNVals, nTVals, nModels)
fractErrVarN = Array{typeof(allErrFlags),3}(undef, nNVals, nTVals, nModels)

parDgpTvarN = Array{Array{Float64,2},2}(undef, nNVals, nTVals)

for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        
        parDgpT = DynNets.dgp_misspecified(model_mle, dgpType, N, T;  dgpOptions...)

        parDgpTvarN[indN, indT] = parDgpT
        for (indM, model) in Iterators.enumerate(models)

            filter_and_conf_bands(model, sample_dgp(model, parDgpT,N), quantilesVals; plotFlag =true, parDgpT = parDgpT)
            
            @elapsed allCoverBuccheri, allCoverBlasques, allvEstSdResPar, allfVecT_filt, allConfBandsBuccheri, allConfBandsBlasques, allErrFlags = conf_bands_coverage(model, parDgpT, N,  nSampleCoverage, quantilesVals)

            allCoverBuccheriVarN[indN, indT, indM] = allCoverBuccheri
            allCoverBlasquesVarN[indN, indT, indM] = allCoverBlasques
            allvEstSdResParVarN[indN, indT, indM] = allvEstSdResPar
            allfVecT_filtVarN[indN, indT, indM] = allfVecT_filt
            allConfBandsBuccheriVarN[indN, indT, indM] = allConfBandsBuccheri
            allConfBandsBlasquesVarN[indN, indT, indM] = allConfBandsBlasques
            fractErrVarN[indN, indT, indM] = allErrFlags
             
        end
    end
end
nErgmPar=2

@save("./data/confBands_$(dgpType)_B_$(dgpOptions.B)_sig_$(dgpOptions.sigma)_(nVals)_$(tVals)_nSample_$nSampleCoverage.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN)
end




begin 
using JLD
nSampleCoverage=55
dgpType = "AR"
models = [model_mle, model_pmle]


nVals = [ 50, 100, 200, 300]
tVals = [100, 200, 300]
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)

@load("./data/confBands_$(dgpType)_B_$(dgpOptions.B)_sig_$(dgpOptions.sigma)_(nVals)_$(tVals)_nSample_$nSampleCoverage.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN)

avgCover =zeros(2,nNVals, nTVals,nModels, 2)
for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, model) in Iterators.enumerate(models)
            for indPar in 1:2
                avgCover[indPar, indN, indT, indM, 1] = mean(allCoverBuccheriVarN[indN, indT, indM][indPar,:,1,:,:]) 
                avgCover[indPar, indN, indT, indM, 2] = mean(allCoverBlasquesVarN[indN, indT, indM][indPar,:,1,:]) 
            end
        end
    end
end


end

begin
figure()
indT = 3
plot(nVals, (avgCover[1,:,indT,1,2]), "--*b")
plot(nVals, (avgCover[2,:,indT,1,2]),"-*b")
plot(nVals, (avgCover[1,:,indT,2,2]), "--*r")
plot(nVals, (avgCover[2,:,indT,2,2]),"-*r")
title("Coverage Conf Bands 95% Dgp $dgpType, T=$(tVals[indT])")
legend(reduce(vcat,["$(mod)  $par" for par in ["θ" "η" ], mod in ["MLE", "PMLE"]]  ))

end

#endregion
