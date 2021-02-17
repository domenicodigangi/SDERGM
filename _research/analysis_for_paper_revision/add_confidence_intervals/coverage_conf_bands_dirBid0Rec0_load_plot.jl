"""
Load and plot simulations results on confidence bands' coverages
"""

#region import and models
using DrWatson
using DataFrames
using TableView
using PyPlot
pygui(true)
using JLD
using ScoreDrivenERGM
using Statistics
using Logging
#endregion


# #region load and plot coverage simulations



df = collect_results!( datadir("sims", "sampleDgpFilterSD_est")) 
df["modelTag"] = string.(df["model"]) 



begin
tVals = [100, 300, 600]
nVals = [100, 200, 300]
modelTags = unique(df["modelTag"])
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(modelTags)
nSampleCoverage = length(unique(df["nSampleCoverage"]))>1 ? missing : (unique(df["nSampleCoverage"])[1] ) 

names(df)
mean(df["allErrFlags"][3])

dgpType = "SD"

avgCover =zeros(2,nNVals, nTVals, nModels, 2, nSampleCoverage)
constInds =falses(2,nNVals, nTVals, nModels, nSampleCoverage)


for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, model) in Iterators.enumerate(modelTags)

            res = filter([:modelTag, :T, :N, :dgpSettings] => (m,t,n,d) -> all((m==model, t==T, n==N, d.type == dgpType)), df)

            size(res)[1] > 1 ? error() : ()            

            @show N, T, mean(res.allErrFlags[1])

            

        fVecT_filt, confBandsFiltPar, confBandsPar, errFlag, mvSDUnParEstCov, distribFilteredSD = conf_bands_given_SD_estimates(model, obsT, vEstSdResPar, quantilesVals; indTvPar = indTvPar, parDgpT=parDgpT, plotFlag=plotFlag, parUncMethod = parUncMethod)

        coverFiltParUnc = conf_bands_coverage(parDgpT,  confBandsFiltPar)

        coverParUnc = conf_bands_coverage(parDgpT,  confBandsPar)

            constInds[:, indN, indT, indM, :]  .= any(res["allvEstSdResPar"][1][3:3:end, :].<=0.02, dims=1)
            

            for indPar in 1:2
                for n=1:nSampleCoverage

                    avgCover[indPar, indN, indT, indM, 1, n] =  mean(res["allCoverFiltPar"][1][indPar,:,1,n]) 

                    avgCover[indPar, indN, indT, indM, 2, n] = mean(res["allCoverPar"][1][indPar,:,1,n]) 
                end
            end
        end
    end
end

indM =2
indB =1

nominalLevel = 0.95
parNames = ["θ", "η"]
BandNames = ["Parameters + Filtering Uncertainty", "Parameters Uncertainty"]

fig, ax1 = plt.subplots(2, length(tVals),figsize=(12, 6), sharey =true)
fig.canvas.set_window_title("Confidence Bands' Coverages $(BandNames[indB])")
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
fig.suptitle("Confidence Bands' Coverages $(BandNames[indB]) DGP = $(dgpType), filter = $(modelTags[indM])")

for (indT, T) in Iterators.enumerate(tVals) 
    for indPar in 1:2
        
        data = [c[.!ind] for (c, ind) in Iterators.zip(eachrow(avgCover[indPar,:,indT,indM,indB,:]), eachrow(constInds[indPar, :, indT, indM, :]))]

        bp = ax1[indPar, indT].boxplot(data, notch=0, sym="+", vert=1, whis=1.5, showfliers =true)


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



dfMleSD = filter( :modelTag=>isequal(modelTags[1]), df)

allVResEstPar = dfMle[(getfield.(dfMle.dgpSettings,:type).=="SD") .& (dfMle.N .== 100 ).& (dfMle.T .== 100),:]["allvEstSdResPar"][1]

mean(allVResEstPar[3:3:end, :].<=0.02, dims=2)


