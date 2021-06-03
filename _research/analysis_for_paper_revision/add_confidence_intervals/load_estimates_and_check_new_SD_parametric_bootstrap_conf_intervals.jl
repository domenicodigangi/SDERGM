#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\add_confidence_intervals\load_estimates_and_check_parametric_bootstrap_conf_intervals.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\add_confidence_intervals
#Created Date: Saturday April 24th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Tuesday April 27th 2021 10:34:44 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########




#region import and models
#%%
using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"

using ProjUtilities
using Logging
using DataFrames
using ScoreDrivenERGM
import ScoreDrivenERGM:DynNets, Utilities, StaticNets, Scalings
using PyPlot
pygui(true)
#endregion


@elapsed df = collect_results( datadir("sims", "dgp&FIl_est"))
df["modelTag"] = string.(df["model"]) 

T = 300
N = 50
nSample = 50
indQuant = 1
modelTag = string(DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType = "FISH_D"))

dgpSetting = DynNets.list_example_dgp_settings(DynNets.SdErgmDirBin0Rec0_mle()).dgpSetARlowlow

# # visual inspection of filter, dgp and conf bands
res = filter([:modelTag, :T, :N, :dgpSettings, :nSample, ] => (mTag,t,n, d, s) -> all((mTag==modelTag, t==T, n==N, d == dgpSetting, s == nSample)), df)[1,:]

using Statistics

model = res.model
quantilesVals = [0.975, 0.025]

n = 15
parDgp_T = res.allParDgpT[:,:,n]
parSdFilt_T = res.allfVecT_filt[:,:, n]
estSdResPar = res.allvEstSdResPar[:,n]
obsT = res.allObsT[n]
parSSFilt_T =  DynNets.estimate_single_snap_sequence(model, obsT)

# coverFlags = @sync @distributed  hcat
nQuant = length(quantilesVals)
confBandsT = zeros(model.staticModel.nErgmPar, T, 1, nQuant)
confBandsT_SS = zeros(model.staticModel.nErgmPar, T, 1, nQuant)
for t = 1:T
    parDgp_t = parDgp_T[:,t]
    parFiltSS_t = parSSFilt_T[:,t]
    parFiltSD_t = res.allfVecT_filt[:,t,n]
    confBandsT[:,t ,1, :] = mapslices((x -> quantile(x, quantilesVals)), StaticNets.get_par_boot_ergm_distrib(model.staticModel, parFiltSD_t, N; nSample = 100), dims=2)
    confBandsT_SS[:,t ,1, :] = mapslices((x -> quantile(x, quantilesVals)), StaticNets.get_par_boot_ergm_distrib(model.staticModel, parFiltSS_t, N; nSample = 100), dims=2)

end

DynNets.plot_filtered_and_conf_bands(model, N, parSdFilt_T, confBandsT; parDgpTIn=parDgp_T)
DynNets.plot_filtered_and_conf_bands(model, N, parSSFilt_T, confBandsT_SS; parDgpTIn=parDgp_T)


begin
obsT = res.allObsT[n]
ftot_0 = res.allftot_0[:,n]
modelIntegrated = deepcopy(model)
modelIntegrated.options["integrated"] = true

~, vEstSdResParIntegrated, fVecT_filt, ~, ~, conv_flag, ftot_0 = DynNets.estimate_and_filter(modelIntegrated, N, obsT; indTvPar = model.indTvPar, show_trace=true)

# vEstSdResParIntegrated[6] = 0.02
testIntFilt_T, _, sVecT, invScalT = DynNets.score_driven_filter( model, N, obsT, vEstSdResParIntegrated, model.indTvPar; ftot_0 = ftot_0 )
DynNets.plot_filtered(model, N, testIntFilt_T; parDgpTIn=parDgp_T)


end


invScalT
sVecT
#%% 
## Sequence for rescaled Par Bootstrap confidence intervals for SDERGM
nParBootSample = 100
ftot_0 = res.allftot_0[:,n]

parBootSampMatsSdFilt_T_S = StaticNets.sample_ergm_sequence(model.staticModel, N, parSdFilt_T, nParBootSample)

#Estimate SD on data and filter
parBootSdFiltDistr = zeros(model.staticModel.nErgmPar, T, nParBootSample)

for pbs=1:nParBootSample
    @info "$pbs"
    obsT = DynNets.seq_of_obs_from_seq_of_mats(model, parBootSampMatsSdFilt_T_S[:,:,:,pbs] )


    # ~, vEstSdResPar, fVecT_filt, ~, ~, conv_flag, ftot_0 = DynNets.estimate_and_filter(model, N, obsT; indTvPar = model.indTvPar)
    # parBootSdFiltDistr[:,:, pbs] = fVecT_filt

    parBootSdFiltDistr[:,:, pbs] = DynNets.score_driven_filter( model, N, obsT, estSdResPar, model.indTvPar; ftot_0 = ftot_0 )[1]

end


nQuant = length(quantilesVals)
confBandsSdParBoot = mapslices((x -> quantile(x, quantilesVals)), parBootSdFiltDistr, dims=3)
confBandsSdParBoot = reshape(confBandsSdParBoot, (model.staticModel.nErgmPar, T, 1, nQuant))
DynNets.plot_filtered_and_conf_bands(model, N, parSdFilt_T, confBandsSdParBoot; parDgpTIn=parDgp_T)



# Obtain Distribution of observations from the filtered path, sampling each ergm_t using the filtered parameters nSample times



# filter each 


