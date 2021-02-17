
using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using JLD
using Distributed
using SharedArrays
using ScoreDrivenERGM
using Logging


model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
indTvPar = trues(2)



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



dgpSetAR, ~, dgpSetSD = ScoreDrivenERGM.DynNets.list_example_dgp_settings_for_paper(model_mle)



c= Dict{String, Any}()
c["T"] = [100]
c["N"] = [100]
c["dgpSettings"] = [dgpSetSD, dgpSetAR]
c["model"] =[model_pmle] 
c["nSample"] = 60

d = dict_list(c)[1]
for d in dict_list(c)
    @show d

    saveName = replace.( savename(d, "jld";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.GasNetModel) ), r"[\"]" => "")

    loadPath = "$(datadir())\\sims\\sampleDgpFilterSD_est\\"
    estDict = wload(loadPath * saveName)

    loadPath = "$(datadir())\\sims\\sampleDgpFilterSD_sampled_mats\\"
    matsDict = wload(loadPath * saveName)
    


    timeObs = @elapsed estDict["allObsT"] = pmap(AT -> seq_of_obs_from_seq_of_mats(model, AT), eachslice(allAT, dims=4))


    timeSave = @elapsed wsave( datadir("sims", "sampleDgpFilterSD_est", saveName), estDict)
    
    Logging.@info("time obs = $timeObs ,  time save = $timeSave ")

end





