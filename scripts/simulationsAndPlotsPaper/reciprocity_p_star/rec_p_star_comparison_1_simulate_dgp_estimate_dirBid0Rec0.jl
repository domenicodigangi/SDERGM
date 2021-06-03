"""
Simulate various Dgps for dirBin0Rec0 various T N and filter
"""


#region import modules and init workers

using Pkg
Pkg.activate(".") 
Pkg.instantiate() 
using DrWatson
using JLD2
using Distributed
using SharedArrays
using ScoreDrivenERGM
using Logging
import ScoreDrivenERGM:StaticNets, DynNets
import ScoreDrivenERGM.DynNets:SdErgm,SdErgmDirBin0Rec0, simulate_and_estimate_parallel
using ScoreDrivenERGM.Utilities

begin
nWorkers = 10
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

    indTvPar = trues(2)
    end
end
#endregion



dgpSetARlowlow, dgpSetARlow, dgpSetARmed, dgpSetARhigh, dgpSetSIN, dgpSetSDlow, dgpSetSD, dgpSetSDhigh = ScoreDrivenERGM.DynNets.list_example_dgp_settings(model_mle)


# define dictionary with all the different values for each setting of the simulation. The product of all settings will be executed using DrWatson.jl functionalities
c= Dict{String, Any}()
c["model"] =[DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType="FISH_D"), DynNets.SdErgmDirBin0Rec0_pmle(scoreScalingType="FISH_D")] 
c["model"][1].options["integrated"] = true
c["model"][2].options["integrated"] = true


c["T"] = [300, 3000]
c["N"] = [100]

dgpSetARlowlow.opt.B[1] = 1
c["dgpSettings"] = [dgpSetARlowlow]
c["nSample"] = 50

list = sort(dict_list(c), by=x->(string(x["model"])))

for d in list

    timeSim = @elapsed allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0, allfVecT_filt_SS =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    modelTag = string(d["model"])

    res1 =  (;modelTag, allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)

    saveName = replace.( savename(d, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.SdErgm) ), r"[\"]" => "")
   

    timeSave = @elapsed save( datadir("sims", "dgp&Fil_est", saveName), estDict)
   
    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")


end

