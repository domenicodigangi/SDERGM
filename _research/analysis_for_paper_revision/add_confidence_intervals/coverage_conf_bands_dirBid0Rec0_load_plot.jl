"""
Load and plot simulations results on confidence bands' coverages
"""

#region import and models
using DrWatson
using PyPlot
pygui(true)
using JLD
#endregion


# #region load and plot coverage simulations


nVals = [300]
tVals = [300]

nSampleCoverage=8

dgpSettings = (dgpType = "AR" , dgpOptions = (minValAlpha = 0.2, maxValAlpha = 0.3, plotFlag=false, B =0.98, sigma = 0.0005))

dgpSettings = (dgpType = "SIN" , dgpOptions = (minValAlpha = 0.2, maxValAlpha = 0.3, nCycles=1.5, phaseAlpha = 0.1π, phaseshift = 0.1, plotFlag=false))

dgpSettings = (dgpType = "SD" , dgpOptions = (plotFlag=false, B =0.98, A = 0.3))

df = collect_results!( datadir("sims", "cofBandsCover")) 



for N ∈ nVals, T ∈ tVals

          
    resDict = @strdict res simulSettings 
    load
    @tagsave( , resDict)


end


begin 

nSampleCoverage=50
dgpType = "SD"
models = [model_pmle]


nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(models)

@load("$(datadir())\\old_pre_drWatson_git_ignored\\confBands_$(dgpType)_B_$(dgpOptions.B)_A_$(dgpOptions.A)_sig_$(dgpOptions.sigma)_$(nVals)_$(tVals)_nSample_$(nSampleCoverage)_pmle.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN, allvEstSdResParVarN)

avgCover =zeros(2,nNVals, nTVals,nModels, 2, nSampleCoverage)

for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, model) in Iterators.enumerate(models)

            allCoverFiltPar, allCoverPar, allvEstSdResPar, allfVecT_filt, allParDgpT, allConfBandsFiltPar,allConfBandsPar, allErrFlags

                    @load("$(datadir())\\old_pre_drWatson_git_ignored\\confBands_$(dgpType)_B_$(dgpOptions.B)_A_$(dgpOptions.A)_sig_$(dgpOptions.sigma)_$(nVals)_$(tVals)_nSample_$(nSampleCoverage)_pmle.jld", allCoverBuccheriVarN, allCoverBlasquesVarN, allfVecT_filtVarN, allConfBandsBuccheriVarN, allConfBandsBlasquesVarN, fractErrVarN, parDgpTvarN, allvEstSdResParVarN)

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


