"""
Simulate various Dgps for dirBin0Rec0 various T N and filter
"""


#region import and models

using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using JLD2
using Distributed
using SharedArrays
using ScoreDrivenERGM
using Logging

begin
nWorkers = 8
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 
    Pkg.instantiate() 
    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:GasNetModel,GasNetModelDirBin0Rec0, simulate_and_estimate_parallel
    using ScoreDrivenERGM.Utilities

    model_mle = DynNets.GasNetModelDirBin0Rec0_mle(scoreScalingType="FISH_D")

    end
end
#endregion



dgpSetARlowlow, dgpSetARlow, dgpSetARmed, dgpSetARhigh, dgpSetSIN, dgpSetSDlow, dgpSetSD, dgpSetSDhigh = ScoreDrivenERGM.DynNets.list_example_dgp_settings(model_mle)

sigmaVals = round.((10).^(range(log10(0.005), log10(0.1), length=6)), sigdigits=1)
dgpList = [deepcopy(dgpSetARlowlow) for i in 1:length(sigmaVals)]
setindex!.(getfield.(getfield.(dgpList, :opt),:sigma), sigmaVals) 

c= Dict{String, Any}()
c["model"] =[DynNets.GasNetModelDirBin0Rec0_mle(scoreScalingType="FISH_D")] 
c["T"] = [500]
c["N"] = [100]
c["dgpSettings"] = dgpList
c["nSample"] = 50

list = sort(dict_list(c), by=x->(string(x["model"])))

for d in list

    timeSim = @elapsed allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0,  allfVecT_filt_SS =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    modelTag = string(d["model"])

    res1 =  (;modelTag, allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0, allfVecT_filt_SS) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)

    saveName = replace.( savename(d, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.GasNetModel) ), r"[\"]" => "")

    timeSave = @elapsed save( datadir("sims", "dgpAR_Fil_SS", saveName), estDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")

    # res2 =  (;allAT ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
    # matsDict = merge(res2, d)
    # @tagsave( datadir("sims", "sampleDgpFilterSD_sampled_mats", saveName), matsDict)

end



# #region load and plot estimates

using Statistics
using DataFrames
using PyPlot
import ScoreDrivenERGM.Utilities.AReg

@elapsed df = collect_results( datadir("sims", "dgpAR_Fil_SS")) 
df["modelTag"] = string.(df["model"]) 


df = transform(df, :allfVecT_filt => ByRow(x-> dropdims(std(x, dims=2), dims=2)) => :std_SD )

df = transform(df, :allfVecT_filt_SS => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(rolling_mean(Float64.(x)), 1)[2] ), y, dims=(2) ), dims=2)) => :std_SS )

df = transform(df, :allfVecT_filt => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_SD )

df = transform(df, :allParDgpT => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_DGP_est )

df = transform(df, :dgpSettings => ByRow(x-> getfield(getfield(x,:opt), :sigma)[1]) => :sigma )

names(df)
qVals = unique(getfield.(getfield.(df.dgpSettings,:opt), :sigma))

df = sort(df, :sigma)

rolling_mean(x::Array{Float64}, size = 3) = vcat(fill(mean(x[i:i+size]), size), [mean(x[i:i+size]) for i =1:length(x)-size])

rolling_mean(randn(100))


df.allvEstSdResPar[6][3p, :]

begin
fig, ax = plt.subplots(nrows=2, ncols=length(qVals), figsize=(9, 4))

# generate some random test data
for i=1:nrow(df)
    for p=1:2
        all_data = [df.std_SS[i][p,:] df.std_SD[i][p,:] df.std_DGP_est[i][p,:]]
        # plot violin plot
        ax[p, i].boxplot(all_data,
                        showmeans=true)
        xlims = ax[p, i].get_xlim()
        ax[p, i].hlines(df.sigma[i], xlims[1], xlims[2], linestyle=":" , colors = "r")
    end
end
   
end

# plot box plot
ax[2].boxplot(all_data)
ax[2].set_title("Box plot")

# adding horizontal grid lines
for ax in ax
    ax.yaxis.grid(true)
    ax.set_xticks([y + 1 for y in 1:length(all_data)])
    ax.set_xlabel("Four separate samples")
    ax.set_ylabel("Observed values")
end
# add x-tick labels
# plt.setp(ax, xticks=[y + 1 for y in 1:length(all_data)],
        #  xticklabels=["x1", "x2", "x3", "x4"])
plt.show()
end

end


nominalLevel = 0.95
parNames = ["θ", "η", "mean θ η"]
BandNames = ["Parameters + Filtering Uncertainty", "Parameters Uncertainty"]

fig, ax1 = plt.subplots(3, length(tVals),figsize=(12, 6), sharey =true)
fig.canvas.set_window_title("Confidence Bands" Coverages $(BandNames[indB])")
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
fig.suptitle("Confidence Bands" Coverages $(BandNames[indB]) DGP = $(dgpSetting.type),\n filter = $(modelTags[indM]), Cov-Estimate : $(parUncMethod)")


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