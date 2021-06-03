


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
    import ScoreDrivenERGM.DynNets:SdErgm,SdErgmDirBin0Rec0, simulate_and_estimate_parallel
    using ScoreDrivenERGM.Utilities

    model_mle = DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType="FISH_D")

    end
end

#endregion


model_mle = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [true, true], scoreScalingType="FISH_D")

#region simulations
dgpSetARlowlow, dgpSetARlow, dgpSetARmed, dgpSetARhigh, dgpSetSIN, dgpSetSDlow, dgpSetSD, dgpSetSDhigh = ScoreDrivenERGM.DynNets.list_example_dgp_settings(model_mle)

sigmaVals = round.((10).^(range(log10(0.005), log10(0.1), length=6)), sigdigits=1)
dgpList = [deepcopy(dgpSetARlowlow) for i in 1:length(sigmaVals)]
setindex!.(getfield.(getfield.(dgpList, :opt),:sigma), sigmaVals) 

c= Dict{String, Any}()
c["model"] =[DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType="FISH_D")] 
c["T"] = [500]
c["N"] = [500]
c["dgpSettings"] = dgpList
c["nSample"] = 50

list = sort(dict_list(c), by=x->(string(x["model"])))

for d in list

    timeSim = @elapsed allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0,  allfVecT_filt_SS =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    modelTag = string(d["model"])

    res1 =  (;modelTag, allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0, allfVecT_filt_SS) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)

    saveName = replace.( savename(d, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.SdErgm) ), r"[\"]" => "")

    timeSave = @elapsed save( datadir("sims", "dgpAR_Fil_SS", saveName), estDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")

    # res2 =  (;allAT ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
    # matsDict = merge(res2, d)
    # @tagsave( datadir("sims", "sampleDgpFilterSD_sampled_mats", saveName), matsDict)

end

#endregion

# #region load and plot estimates
begin
using Statistics
using DataFrames
using PyPlot
import ScoreDrivenERGM.Utilities.AReg

rolling_mean(x::Array{Float64}, size = 3) = vcat(fill(mean(x[i:i+size]), size), [mean(x[i:i+size]) for i =1:length(x)-size])


@elapsed df = collect_results( datadir("sims", "dgpAR_Fil_SS")) 
df["modelTag"] = string.(df["model"]) 

N = 100
modTag = string(model_mle)

df = df[df.N .== N, :]
df = transform(df, :allfVecT_filt => ByRow(x-> dropdims(std(x, dims=2), dims=2)) => :std_SD )

df = transform(df, :allfVecT_filt_SS => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_SS )

df = transform(df, :allfVecT_filt => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_SD )

df = transform(df, :allParDgpT => ByRow(y-> dropdims(mapslices( (x-> AReg.fitARp(Float64.(x), 1)[2] ), y, dims=(2) ), dims=2)) => :std_DGP_est )

df = transform(df, :dgpSettings => ByRow(x-> getfield(getfield(x,:opt), :sigma)[1]) => :sigma )

end



begin
df = sort(df, :sigma)
qVals = unique(getfield.(getfield.(df.dgpSettings,:opt), :sigma))
fig, ax = plt.subplots(nrows=2, ncols=length(qVals), figsize=(12, 8), sharey = "row")

labNames = ["theta", "eta"]
for i=1:nrow(df)
    for p=1:2
        all_data = ([df.std_SS[i][p,:] df.std_SD[i][p,:] df.std_DGP_est[i][p,:]] .-df.sigma[i])#./df.sigma[i]
        # plot violin plot
        ax[p, i].boxplot(all_data,
                        showmeans=true, showfliers=false)
        xlims = ax[p, i].get_xlim()
        ax[p, i].hlines(0, xlims[1], xlims[2], linestyle=":" , colors = "r")

        ax[p,i].set_xticklabels(["SS", "SD", "DGP"])
        ax[p,i].grid()
        ax[p,1].set_ylabel(labNames[p])
        ax[1,i].set_title("var dgp = $(qVals[i][1]) ")
    end
end
   
fig.suptitle("Estimate Variance of Latent AR1 in reciprocity p* model N = $N")
plt.tight_layout()

plt.savefig("C:\\Users\\digan\\Dropbox\\Dynamic_Networks\\repos\\ScoreDrivenExponentialRandomGraphs\\_research\\analysis_for_paper_revision\\plots_for_revision\\estimate_var_latent_AR1_N_$(N)_$modTag")
end
#endregion

