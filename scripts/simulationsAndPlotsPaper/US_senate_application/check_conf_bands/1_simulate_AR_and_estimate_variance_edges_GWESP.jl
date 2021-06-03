#########
#Created Date: Wednesday April 28th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Tuesday June 1st 2021 7:55:51 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########





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
nWorkers = 10
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 

    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:SdErgm,SdErgmDirBin0Rec0, simulate_and_estimate_parallel
    using ScoreDrivenERGM.Utilities


    end
end

#endregion

model = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [true, true], scoreScalingType="FISH_D")


dgpSettings = (type = "AR", opt = (B =[1], sigma = [0, 0.18],  minVals = [-3.0, 0], maxVals = [-0.1, 2] ))


sigmaVals = [0.005, 0.1, 0.2]
dgpList = [deepcopy(dgpSettings) for i in 1:length(sigmaVals)]
[dgpList[ind].opt.sigma[2] = sigma for (ind, sigma) in enumerate(sigmaVals)]
 


#region simulations
c= Dict{String, Any}()
c["model"] =[ DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [false, true], scoreScalingType="FISH_D")] 
c["T"] = [73]
c["N"] = [100]
c["dgpSettings"] = dgpList
c["nSample"] = 50

list = sort(dict_list(c), by=x->(string(x["model"])))

for d in list
    nSample = d["nSample"]
    T = d["N"]
    timeSim = @elapsed allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0,  allfVecT_filt_SS =  ScoreDrivenERGM.DynNets.simulate_and_estimate_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSample"];)
                
    modelTag = string(d["model"])

    res1 =  (;modelTag, allObsT, allvEstSdResPar, allfVecT_filt, allParDgpT, allConvFlag, allftot_0, allfVecT_filt_SS) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    estDict = merge(res1, d)

    saveName = replace.( savename(d;allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.SdErgm) ), r"[\"]" => "")

    estDict["saveName"] = saveName

    timeSave = @elapsed save( datadir("sims", "dgpAR_Fil_SS", string(hash(saveName))) *".jld2", estDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")

end

 