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
nWorkers = 15
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

    indTvPar = trues(2)
    end
end
#endregion



dgpSetARlowlow, dgpSetARlow, dgpSetARmed, dgpSetARhigh, dgpSetSIN, dgpSetSDlow, dgpSetSD, dgpSetSDhigh = ScoreDrivenERGM.DynNets.list_example_dgp_settings(model_mle)



c= Dict{String, Any}()
c["model"] =[DynNets.GasNetModelDirBin0Rec0_mle(scoreScalingType="FISH_D"), DynNets.GasNetModelDirBin0Rec0_pmle(scoreScalingType="FISH_D")] 
c["T"] = [100, 300, 600]
c["N"] = [50, 100, 500]
c["dgpSettings"] = [dgpSetARlowlow]
c["nSample"] = 50

list = sort(dict_list(c), by=x->(string(x["model"])))

for d in list

    timeSim = @elapsed allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0, allfVecT_filt_SS =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    modelTag = string(d["model"])

    res1 =  (;modelTag, allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)

    saveName = replace.( savename(d, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.GasNetModel) ), r"[\"]" => "")

    timeSave = @elapsed save( datadir("sims", "samDgpFiltSD_est", saveName), estDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")

    # res2 =  (;allAT ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
    # matsDict = merge(res2, d)
    # @tagsave( datadir("sims", "sampleDgpFilterSD_sampled_mats", saveName), matsDict)

end

