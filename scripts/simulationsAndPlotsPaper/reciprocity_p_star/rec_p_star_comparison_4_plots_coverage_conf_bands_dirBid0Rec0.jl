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
using SharedArrays
using Statistics
using ScoreDrivenERGM
using TableView
using Blink
viewtab(df) = body!(Window(), showtable(df))
#endregion


# #region load and plot coverage simulations

@elapsed df = collect_results( datadir("sims", "dgp&FIl_conf")) 
df["modelTag"] = string.(df["model"]) 


df = df[getfield.(df.model, :scoreScalingType) .== "FISH_D",:] 
df = df[getfield.(df.dgp, :type) .== "AR",:] 
viewtab(df[[:modelTag, :N, :T]])


df.confBndCoverFiltPar = [zeros(1, 50) for i = 1: nrow(df)]            
df.confBndCoverPar = [zeros(1, 50) for i = 1: nrow(df)]            
df.avgConfBndCoverFiltPar = zeros(nrow(df))            
df.avgConfBndCoverPar = zeros(nrow(df))            


obsShift = 1
for res in eachrow(df)
    T = res.T
    limitSample = sum(res.allConfBandsFiltPar[1,1,1,1,:] .!= 0 )
    nBands = size(res.allConfBandsFiltPar)[3]
    avgCoverFiltPar = zeros(nBands, limitSample)
    avgCoverPar = zeros(nBands, limitSample)

    constInds = any(res.allvEstSdResPar[3:3:6, 1:limitSample] .< 0.00005, dims=1)

    # Threads.@threads 
    for n = 1:limitSample
        coverFiltParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,1:T-obsShift,n],   res.allConfBandsFiltPar[:,1+obsShift:end,:,:,n])

        coverParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,1:T-obsShift,n],   res.allConfBandsPar[:,1+obsShift:end,:,:,n])

        avgCoverFiltPar[1, n] =  mean(coverFiltParUnc[:,:,1]) 
        avgCoverPar[1, n] =  mean(coverFiltParUnc[:,:,1]) 

    end

    res.avgConfBndCoverFiltPar = mean(avgCoverFiltPar[.!constInds])            
    res.avgConfBndCoverPar = mean(avgCoverPar[.!constInds])        

    res.confBndCoverFiltPar = avgCoverFiltPar            
    res.confBndCoverPar = avgCoverPar        
end    

viewtab(df[df.N.==100,:][[:modelTag, :m,  :dgp, :N, :T, :avgConfBndCoverFiltPar, :confBndCoverFiltPar]])


#region obtain box plots to compare subsets of the simulations run
begin 

parUncMethod = "WHITE-MLE" 
parUncMethod = "NPB-MVN" #
limitSample =50
indB = 1
tVals = [300]
nVals = [100  ]
model = DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType="FISH_D")
model = DynNets.SdErgmDirBin0Rec0_pmle(scoreScalingType="FISH_D")
dgpSetting = DynNets.list_example_dgp_settings(DynNets.SdErgmDirBin0Rec0_mle()).dgpSetARlowlow


modelTags = [DynNets.name(model)] # ["SdErgmDirBin0Rec0_mle(Bool[1, 1], scal = HESS_D)"]# unique(df["modelTag"])
# modelTags =unique(df["modelTag"])
nNVals = length(nVals)
nTVals = length(tVals)
nModels = length(modelTags)
nSample = 50 #length(unique(df["nSample"]))>1 ? missing : (unique(df["nSample"])[1] ) 
nErgmPar = 2



indQuant = 1
nBands = 2
allMeanDgpStd =zeros(nNVals, nTVals, nModels)
allMeanDgpStdConf =zeros(nNVals, nTVals, nModels)
allAvgCover =zeros(nErgmPar ,nNVals, nTVals, nModels, nBands, limitSample)
allConstInds = falses(nErgmPar, nNVals, nTVals, nModels, limitSample)
allErrInds = trues(nErgmPar, nNVals, nTVals, nModels, limitSample)

avgCover = SharedArray(zeros(nErgmPar, nBands, limitSample))
avgCoverJoint = SharedArray(zeros(nErgmPar, nBands, limitSample))
indM=1
obsShift = 1



for (indT, T) in Iterators.enumerate(tVals) 
    for (indN, N) in Iterators.enumerate(nVals) 
        for (indM, modelTag) in Iterators.enumerate(modelTags)

            res = filter([:modelTag, :T, :N, :dgp, :S, :m] => (mtag,t,n, d, s, m) -> all((mtag==modelTag, t==T, n==N, d == dgpSetting, s == nSample, m==parUncMethod)), df)

            nrow(res) != 1 ? error("$N, $T,  $(size(res))") : res = res[1,:]
         

            # Threads.@threads 
            for n = 1:limitSample
                coverFiltParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,1:T-obsShift,n],   res.allConfBandsFiltPar[:,1+obsShift:end,:,:,n])

                coverParUnc = DynNets.conf_bands_coverage(res.allParDgpT[:,1:T-obsShift,n],   res.allConfBandsPar[:,1+obsShift:end,:,:,n])

                for indPar = 1:2
                    avgCover[indPar, 1, n] =  mean(coverFiltParUnc[indPar,:,1]) 

                    avgCover[indPar, 2, n] =  mean(coverParUnc[indPar,:,1]) 
                end
                
            end
            
            variability = res.allfVecT_filt |>  x -> dropdims(std(x, dims=2), dims=2) |> x-> replace(e -> (isnan(e) ? 0 : e), x ) |> x ->((mean = mean(x) , conf =  1.96.*std(x) ))


            allMeanDgpStd[indN, indT, indM] = round.(variability.mean, sigdigits=2)  
            allMeanDgpStdConf[indN, indT, indM] = round.(variability.conf, sigdigits=2)  
            allAvgCover[:, indN, indT, indM, :, 1:limitSample] = avgCover 
            allErrInds[:, indN, indT, indM, 1:limitSample] = res.errInds[:, 1:limitSample]
            allConstInds[:, indN, indT, indM, 1:limitSample] .= any(res.allvEstSdResPar[3:3:6, 1:limitSample] .< 0.00005, dims=1)

        end
    end
end


nominalLevel = 0.95
parNames = ["θ", "η", "mean θ η"]
BandNames = ["Parameters + Filtering Uncertainty", "Parameters Uncertainty"]

fig, ax1 = plt.subplots(3, length(tVals),figsize=(12, 6), sharey =true)
fig.canvas.set_window_title("Confidence Bands' Coverages $(BandNames[indB])")
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
fig.suptitle("Confidence Bands' Coverages $(BandNames[indB]) DGP = $(dgpSetting.type),\n filter = $(modelTags[indM]), Cov-Estimate : $(parUncMethod)")


for (indT, T) in Iterators.enumerate(tVals) 
    for indPar = 1:3
        
        if indPar == 3
            data = [c[(.!indC).&(.!indE)] for (c, indC, indE) in Iterators.zip(eachrow(dropdims(mean(allAvgCover[:,:,indT,indM,indB,:], dims=1), dims=1)), eachrow(allConstInds[1, :, indT, indM, :]), eachrow(allErrInds[1, :, indT, indM, :]))]
        else
            data = [c[(.!indC).&(.!indE)] for (c, indC, indE) in Iterators.zip(eachrow(allAvgCover[indPar,:,indT,indM,indB,:]), eachrow(allConstInds[indPar, :, indT, indM, :]), eachrow(allErrInds[indPar, :, indT, indM, :]))]
        end
        # data = [c for (c, ind) in Iterators.zip(eachrow(allAvgCover[indPar,:,indT,indM,indB,:]), eachrow(allConstInds[indPar, :, indT, indM, :]))]

        bp = ax1[indPar, indT].boxplot(data, notch=0, sym="+", vert=1, whis=1.5, showfliers =true, showmeans=true)


        ax1[indPar, indT].yaxis.grid(true, linestyle="-", which="major", color="lightgrey", alpha=0.5)

        # Hide these grid behind plot objects
        xlims = ax1[indPar, indT].get_xlim()
        ax1[indPar, indT].hlines(nominalLevel, xlims[1], xlims[2], linestyle=":" , colors = "r")
        ax1[indPar, indT].set_ylim([0.70, 1])
        ax1[indPar, indT].set_axisbelow(true)
        ax1[indPar, indT].set_title("T = $T")
        # ax1[indPar, indT].set_xlabel("Network Size")
        ax1[indPar, indT].set_ylabel("$(parNames[indPar])")
        if indPar ==3
            xlab = ["N = $n \n <std($(parNames[indPar]) )> =  $(allMeanDgpStd[indN, indT, indM])  " for (indN,n) in Iterators.enumerate(nVals) ] #\n ($(allMeanDgpStdConf[indPar, indN, indT, indM])
            # ax1[indPar, indT].set_xticklabels(xlab, rotation=0, fontsize=8)
        end
    end
end
tight_layout()

end
#endregion