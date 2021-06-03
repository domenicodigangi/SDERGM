#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\scripts\simulationsAndPlotsPaper\reciprocity_p_star\simulate_AR_and_estimate_variance_edges_GWESP.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\scripts\simulationsAndPlotsPaper\reciprocity_p_star
#Created Date: Wednesday April 28th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Sunday May 30th 2021 12:04:11 pm
#Modified By:  Domenico Di Gangi
#-----
#Description: test sequence of simulations in the current folder for a single realization
#-----
########



using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"

using ProjUtilities
using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets, DynNets, Utilities


model = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [true, true], scoreScalingType="FISH_D")


N = 100
T = 110

import ScoreDrivenERGM.DynNets: reference_model, seq_of_obs_from_seq_of_mats, estimate_and_filter
import ScoreDrivenERGM.Utilities.AReg

dgpSettings = (type = "AR", opt = (B =[0.95], sigma = [0, 0.2],  minVals = [-3.0, 0], maxVals = [-0.1, 2] ))

parDgpT = DynNets.sample_time_var_par_from_dgp(reference_model(model), dgpSettings.type, N, T;  dgpSettings.opt..., maxAttempts=100000)

@elapsed A_T_dgp = StaticNets.sample_ergm_sequence(reference_model(model).staticModel, N, parDgpT, 1)[:,:,:,1]

obsT = seq_of_obs_from_seq_of_mats(model, A_T_dgp)

filtModel = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [false, true], scoreScalingType="FISH_D")

@elapsed ~, vEstSdResPar, fVecT_filt, ~, ~, conv_flag, ftot_0 = estimate_and_filter(filtModel, N, obsT;  show_trace = true)

using PyPlot

plot(fVecT_filt[2,:])
plot(parDgpT[2,:])


mapslices( x-> AReg.fitARp(Float64.(x), 1)[2], fVecT_filt, dims=2 )


fVecT_filt
#endregion
