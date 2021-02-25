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

#endregion


# #region load and plot coverage simulations

df = collect_results( datadir("sims", "samDgpFiltSD_est_conf")) 
a=1

begin
tVals = [100, 300]
nVals = [100, 200, 300]
# modelTags = ["GasNetModelDirBin0Rec0_mle(Bool[1, 1], scal = HESS)"]# unique(df["modelTag"])
modelTags =unique(df["modelTag"])
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(modelTags)
nSample = 100 #length(unique(df["nSample"]))>1 ? missing : (unique(df["nSample"])[1] ) 
limitSample =nSample

dgpType = "AR"

allAvgCover =zeros(2,nNVals, nTVals, nModels, 2, nSample)
allConstInds = falses(2,nNVals, nTVals, nModels, nSample)
allErrInds = trues(2,nNVals, nTVals, nModels, nSample)


for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, modelTag) in Iterators.enumerate(modelTags)

            res = filter([:modelTag, :T, :N, :dgpSettings, :nSample] => (m,t,n, d, s) -> all((m==modelTag, t==T, n==N, d.type == dgpType, s == nSample)), df)[1,:]


            allAvgCover[:, indN, indT, indM, :, 1:limitSample] = res.avgCover 
            allErrInds[:, indN, indT, indM, 1:limitSample] = res.errInds 
            allConstInds[:, indN, indT, indM, 1:limitSample] = res.constInds
        end
    end
end

end


begin

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
    for indPar = 1:2
        
        data = [c[.!ind] for (c, indC, indE) in Iterators.zip(eachrow(allAvgCover[indPar,:,indT,indM,indB,:]), eachrow(allConstInds[indPar, :, indT, indM, :]), eachrow(allConstInds[indPar, :, indT, indM, :]))]
        # data = [c for (c, ind) in Iterators.zip(eachrow(allAvgCover[indPar,:,indT,indM,indB,:]), eachrow(allConstInds[indPar, :, indT, indM, :]))]

        bp = ax1[indPar, indT].boxplot(data, notch=0, sym="+", vert=1, whis=1.5, showfliers =true)


        ax1[indPar, indT].yaxis.grid(true, linestyle="-", which="major", color="lightgrey", alpha=0.5)

        # Hide these grid behind plot objects
        xlims = ax1[indPar, indT].get_xlim()
        ax1[indPar, indT].hlines(nominalLevel, xlims[1], xlims[2], linestyle=":" , colors = "r")
        ax1[indPar, indT].set_ylim([0.5, 1])
        ax1[indPar, indT].set_axisbelow(true)
        ax1[indPar, indT].set_title("T = $T")
        ax1[indPar, indT].set_xlabel("Network Size")
        ax1[indPar, indT].set_ylabel("Coverages for $(parNames[indPar])")
        ax1[indPar, indT].set_xticklabels(nVals, rotation=45, fontsize=8)
    end
end
tight_layout()

end
