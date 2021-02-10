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
using JLD
using Distributed
using SharedArrays
using ScoreDrivenERGM

begin
nWorkers = 8
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 
    Pkg.instantiate() 
    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:GasNetModel,GasNetModelDirBin0Rec0, sample_dgp, statsFromMat, array2VecGasPar, unrestrict_all_par, number_ergm_par, filter_and_conf_bands, conf_bands_coverage_parallel, estimate, plot_filtered_and_conf_bands
    using ScoreDrivenERGM.Utilities

    model_mle = DynNets.GasNetModelDirBin0Rec0_mle()
    model_pmle = DynNets.GasNetModelDirBin0Rec0_pmle()
    indTvPar = trues(2)

    end

end
#endregion


# #region coverage simulations
dgpSetAR, dgpSetSIN, dgpSetSD = ScoreDrivenERGM.DynNets.list_example_dgp_settings_for_paper(model_mle)

nVals = [100, 300]
tVals = [100, 200, 300]
models = [model_mle, model_pmle] 
nSampleCoverage= 160




for N ∈ nVals, T ∈ tVals, dgpSettings ∈ [dgpSetAR, dgpSetSD], model ∈ models 

    allCoverFiltPar, allCoverPar, allvEstSdResPar, allfVecT_filt, allParDgpT, allConfBandsFiltPar,allConfBandsPar, allErrFlags =  ScoreDrivenERGM.DynNets.conf_bands_coverage_parallel(model, dgpSettings, T, N,  nSampleCoverage; quantilesVals = [[0.975, 0.025]])
                
    res = (;allCoverFiltPar, allCoverPar, allvEstSdResPar, allfVecT_filt, allParDgpT, allConfBandsFiltPar,allConfBandsPar, allErrFlags )

    simulSettings = (;dgpSettings..., N, T, nSampleCoverage, model = ScoreDrivenERGM.DynNets.name(model) )

        
    resDict = @strdict res simulSettings 
    @tagsave( datadir("sims", "cofBandsCover", savename(simulSettings, "jld")), resDict)

end

