"""
Simulate various Dgps for dirBin0Rec0 various T N and filter
"""


#region import and models

using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using JLD
using Distributed
using SharedArrays
using ScoreDrivenERGM
using Logging

begin
nWorkers = 12
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
c["T"] = [100, 300, 600]
c["N"] = [100, 200, 300]
c["dgpSettings"] = [dgpSetSD, dgpSetAR]
c["model"] =[model_pmle] 
c["nSample"] = 180

list = dict_list(c)


for d in dict_list(c)

    timeSim = @elapsed allAT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    res1 =  (;allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)
    
    timeObs = @elapsed estDict["allObsT"] = pmap(AT -> seq_of_obs_from_seq_of_mats(model, AT), eachslice(allAT, dims=4))

    saveName = replace.( savename(d, "jld";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.GasNetModel) ), r"[\"]" => "")

    timeSave = @elapsed @tagsave( datadir("sims", "sampleDgpFilterSD_est", saveName), estDict)

    Logging.@info("Time sim = $timeSim , time obs = $timeObs ,  time save = $timeSave ")

    # res2 =  (;allAT ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict
    # matsDict = merge(res2, d)
    # @tagsave( datadir("sims", "sampleDgpFilterSD_sampled_mats", saveName), matsDict)

end


