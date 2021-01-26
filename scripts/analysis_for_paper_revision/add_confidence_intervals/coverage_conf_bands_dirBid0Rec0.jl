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
# Filtro e conf bands con blasques per un seno con T lungo
# SE COVERAGES ANCORA NON BUONE, AGGIUNGERE RISCALAMENTO DELLO SCORE IN ENTRAMBI I MODELLI. PER PMLE POSSO USARE LO STIMATORE DI WHITE AD UN TEMPO, O EWMA SU TUTTI I TEMPI PRECEDENTI 


#region filter checks 
begin
T=300
N=100
quantilesVals = [[0.975, 0.025]]
parDgpT = DynNets.dgp_misspecified(model_mle, "AR", N, T;  minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.98, sigma = 0.001)

# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
A_T_dgp = sample_dgp(model_mle, parDgpT,N)


res_mle = filter_and_conf_bands(model_mle, A_T_dgp, quantilesVals; plotFlag =true, parDgpT = parDgpT)

res_pmle = filter_and_conf_bands(model_pmle, A_T_dgp, quantilesVals; parDgpT = parDgpT, plotFlag=true)


end
#endregion

res_mle[8]
res_pmle[8]


res_mle[8]
res_pmle[8]

#region distribution of filtered parameters checks 
begin
T=100
N=50
model = model_pmle
quantilesVals = [[0.975, 0.025] ]

#parDgpT = DynNets.dgp_misspecified(model_mle, "SIN", N, T;  minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.990, sigma = 0.005)
# quick visual checks
# DynNets.sample_est_mle_pmle(model_mle, parDgpT, N, 1; plotFlag = true)
#A_T_dgp = sample_dgp(model_mle, parDgpT,N)
obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

estSdResPar, conv_flag, UM_mple, ftot_0 = estimate(model, obsT;indTvPar=indTvPar, indTargPar=falses(2))

vEstSdResPar = array2VecGasPar(model, estSdResPar, indTvPar)

fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter(model,  vEstSdResPar, indTvPar; obsT = obsT, ftot_0 = ftot_0)


distribFilteredSD, filtCovHatSample, errFlag = mle_distrib_filtered_par(model, obsT, indTvPar, ftot_0, vEstSdResPar)



confBandsParFilt, confBandsParGauss = DynNets.conf_bands_buccheri(model, obsT, indTvPar, fVecT_filt, distribFilteredSD, filtCovHatSample, quantilesVals)

confBandsParSimul = DynNets.conf_bands_par_uncertainty_blasques(model, obsT, fVecT_filt, distribFilteredSD, quantilesVals)




plot_filtered_and_conf_bands(model, N, fVecT_filt, confBandsParFilt;  parDgpT=parDgpT, confBands2=confBandsParSimul)



deltas = distribFilteredSD[:,1,end-10,] 

fig1, ax = subplots(2,1)
for p in 1:2
    x = 1:T
  
    ax[p].hist(deltas)
  
    ax[p].grid()  
end

end


#endregion

#region coverage simulations
begin 
model = model_pmle
N=50
T=200
dgpType =  "AR"
dgpOptions = (minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.95, sigma = 0.001)
quantilesVals = [[0.975, 0.025]]

parDgpT = DynNets.dgp_misspecified(model_mle, dgpType, N, T;  dgpOptions...)

allCoverBuccheri, allCoverBlasques, allvEstSdResPar, allfVecT_filt, allConfBandsBuccheri,allConfBandsBlasques, allErrFlags = conf_bands_coverage(model, parDgpT, N, 2, quantilesVals)


nVals = [ 50, 100, 200, 300]
tVals = [100, 200, 300]
models = [model_mle, model_pmle]
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)
nSampleCoverage=2

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

@save("./data/confBands_$dgpType_$(nVals)_$(tVals)_nSample_$nSampleCoverage.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN)


avgCover =zeros(2,nNVals, nTVals,length(models), 2)
for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, model) in Iterators.enumerate(models)
            for indPar in 1:2
                avgCover[indPar, indN, indT, indM, 1] = mean(allConfBandsBuccheriVarN[indN, indT, indM][indPar,:,1,:]) 
                avgCover[indPar, indN, indT, indM, 2] = mean(allConfBandsBlasquesVarN[indN, indT, indM][indPar,:,1,:]) 
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

end

#endregion
