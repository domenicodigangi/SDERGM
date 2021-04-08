"""
Simulations to estimate coverage of confidence bands with Blasques and Buccheri's methods 
To Do : 
- profile R interface to see if simulations for PMLE can be faster
- 
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

begin
nWorkers = 12
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 
    Pkg.instantiate() 
    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:GasNetModel,GasNetModelDirBin0Rec0, sample_mats_sequence, stats_from_mat, array2VecGasPar, unrestrict_all_par, number_ergm_par, estimate_filter_and_conf_bands, conf_bands_coverage_parallel, estimate, plot_filtered_and_conf_bands
    using ScoreDrivenERGM.Utilities

    model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
    model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
    indTvPar = trues(2)

    end

end
#endregion


# #region coverage simulations
dgpSetAR, ~, dgpSetSD = ScoreDrivenERGM.DynNets.list_example_dgp_settings(model_mle)


c= Dict{String, Any}()
c["T"] = [100, 300, 600]
c["N"] = [100, 200, 300]
c["dgpSettings"] = [dgpSetSD, dgpSetAR]
c["model"] =[ model_pmle] 
c["nSampleCoverage"] = 60

list = dict_list(c)


for d in dict_list(c)

    allCoverFiltPar, allCoverPar, allvEstSdResPar, allfVecT_filt, allParDgpT, allConfBandsFiltPar,allConfBandsPar, allErrFlags =  ScoreDrivenERGM.DynNets.conf_bands_coverage_parallel(d["model"], d["dgpSettings"], d["T"], d["N"],  d["nSampleCoverage"]; quantilesVals = [[0.975, 0.025]])
                
    res =  (;allCoverFiltPar, allCoverPar, allvEstSdResPar, allfVecT_filt, allParDgpT, allConfBandsFiltPar,allConfBandsPar, allErrFlags ) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    #res = load(datadir("sims", "cofBandsCover", saveName))

    resDict = merge(res, d)
    
    saveName = replace.( savename(d, "jld";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.GasNetModel) ), r"[\"]" => "")

   @tagsave( datadir("sims", "cofBandsCover", saveName), resDict)

end


