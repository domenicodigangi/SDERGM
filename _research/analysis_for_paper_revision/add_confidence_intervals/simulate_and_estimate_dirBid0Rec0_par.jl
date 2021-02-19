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
nWorkers = 11
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 
    Pkg.instantiate() 
    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:GasNetModel,GasNetModelDirBin0Rec0, simulate_and_estimate_parallel
    using ScoreDrivenERGM.Utilities

    model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
    model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
    indTvPar = trues(2)
    end
end
#endregion


dgpSetAR, ~, dgpSetSD = ScoreDrivenERGM.DynNets.list_example_dgp_settings_for_paper(model_mle)


c= Dict{String, Any}()
c["model"] =[model_mle, model_pmle] 
c["T"] = [100, 300, 600]
c["N"] = [100, 200, 300]
c["dgpSettings"] = [dgpSetSD, dgpSetAR]
c["nSample"] = 100

list = sort(dict_list(c), by=x->(string(x["model"])))

for d in list

    timeSim = @elapsed allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0 =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    modelTag = string(d["model"])

    res1 =  (;modelTag, allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)
    

    saveName = replace.( savename(d, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.GasNetModel) ), r"[\"]" => "")

    timeSave = @elapsed @tagsave( datadir("sims", "sampleDgpFilterSD_est", saveName), estDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")

    # res2 =  (;allAT ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
    # matsDict = merge(res2, d)
    # @tagsave( datadir("sims", "sampleDgpFilterSD_sampled_mats", saveName), matsDict)

end

