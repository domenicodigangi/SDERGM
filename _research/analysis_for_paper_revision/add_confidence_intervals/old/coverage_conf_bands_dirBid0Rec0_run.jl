"""
Simulations to estimate coverage of confidence bands with Blasques methods and Buccheri's enhancement
"""

#region import and models

begin 
using ScoreDrivenERGM

import ScoreDrivenERGM:StaticNets, DynNets

import ScoreDrivenERGM.DynNets:SdErgm,SdErgmDirBin0Rec0, sample_mats_sequence, stats_from_mat, array2VecGasPar, unrestrict_all_par, conf_bands_par_uncertainty, avg_grad_and_hess_obj_SD_filter_time_seq, conf_bands_par_uncertainty, number_ergm_par, estimate_filter_and_conf_bands, conf_bands_coverage, estimate, mle_distrib_filtered_par, plot_filtered_and_conf_bands
using ScoreDrivenERGM.Utilities

using PyPlot
pygui(true)


using ForwardDiff
using StatsBase
using LinearAlgebra
using Distributions
using Statistics

using JLD2

model_mle = DynNets.SdErgmDirBin0Rec0_mle()
model_pmle = DynNets.SdErgmDirBin0Rec0_pmle()
indTvPar = trues(2)

end
#endregion





#region coverage simulations
begin 
model = model_mle
N=50
T=200
dgpType =  "SD"

nVals = [ 50, 100, 300]
tVals = [100, 300]
models = [model_mle]
nSampleCoverage=50


dgpOptions = (minAlpha = 0.2, maxAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false, B =0.98, sigma = 0.0005, A = 0.3)
quantilesVals = [[0.975, 0.025]]

parDgpT = DynNets.sample_time_var_par_from_dgp(model_mle, dgpType, N, T;  dgpOptions...)

@elapsed allCoverBuccheri, allCoverBlasques, allvEstSdResPar, allfVecT_filt, allConfBandsBuccheri,allConfBandsBlasques, allErrFlags = conf_bands_coverage(model, dgpType, dgpOptions, T, N, 2, quantilesVals)

nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)

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
        

        parDgpTvarN[indN, indT] = parDgpT
        for (indM, model) in Iterators.enumerate(models)

            @show dgpType
            @show model

            parDgpT = DynNets.sample_time_var_par_from_dgp(model_mle, dgpType, N, T;  dgpOptions...)

            estimate_filter_and_conf_bands(model, sample_mats_sequence(model, parDgpT,N), quantilesVals; plotFlag =false, parDgpT = parDgpT)
            
            @elapsed allCoverBuccheri, allCoverBlasques, allvEstSdResPar, allfVecT_filt, allConfBandsBuccheri, allConfBandsBlasques, allErrFlags = conf_bands_coverage(model, dgpType, dgpOptions, T, N,  nSampleCoverage, quantilesVals)

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

# @save("./data/confBands_$(dgpType)_B_$(dgpOptions.B)_A_$(dgpOptions.A)_sig_$(dgpOptions.sigma)_(nVals)_$(tVals)_nSample_$(nSampleCoverage)_mle.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN)
end

figure()
plot(allfVecT_filtVarN[4,2,1][1,:,:])



begin 
using JLD2
nSampleCoverage=50
dgpType = "SD"
models = [model_mle]


nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)

@load("$(datadir())\\old_pre_drWatson_git_ignored\\confBands_$(dgpType)_B_$(dgpOptions.B)_A_$(dgpOptions.A)_sig_$(dgpOptions.sigma)_(nVals)_$(tVals)_nSample_$(nSampleCoverage)_mle.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN)

avgCover =zeros(2,nNVals, nTVals,nModels, 2, nSampleCoverage)
for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, model) in Iterators.enumerate(models)
            for indPar in 1:2
                for n=1:nSampleCoverage
                    avgCover[indPar, indN, indT, indM, 1, n] = mean(allCoverBuccheriVarN[indN, indT, indM][indPar,:,1,n]) 
                    avgCover[indPar, indN, indT, indM, 2, n] = mean(allCoverBlasquesVarN[indN, indT, indM][indPar,:,1,n]) 
                end
            end
        end
    end
end



indM =1
indB =2

nominalLevel = 0.95
parNames = ["θ", "η"]
BandNames = ["Parameters + Filtering Uncertainty", "Parameters Uncertainty"]

fig, ax1 = plt.subplots(2, length(tVals),figsize=(12, 6), sharey =true)
fig.canvas.set_window_title("Confidence Bands' Coverages $(BandNames[indB])")
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
fig.suptitle("Confidence Bands' Coverages $(BandNames[indB]) DGP = $(dgpType)")

for (indT, T) in Iterators.enumerate(tVals) 
    for indPar in 1:2
        
        data = avgCover[indPar,:,indT,indM,indB,:]'
        bp = ax1[indPar, indT].boxplot(data, notch=0, sym="+", vert=1, whis=1.5)


        ax1[indPar, indT].yaxis.grid(true, linestyle="-", which="major", color="lightgrey", alpha=0.5)

        # Hide these grid behind plot objects
        xlims = ax1[indPar, indT].get_xlim()
        ax1[indPar, indT].hlines(nominalLevel, xlims[1], xlims[2], linestyle=":" , colors = "r")
        ax1[indPar, indT].set_axisbelow(true)
        ax1[indPar, indT].set_title("T = $T")
        ax1[indPar, indT].set_xlabel("Network Size")
        ax1[indPar, indT].set_ylabel("Coverages for $(parNames[indPar])")
        ax1[indPar, indT].set_xticklabels(nVals, rotation=45, fontsize=8)
    end
end
tight_layout()

end


