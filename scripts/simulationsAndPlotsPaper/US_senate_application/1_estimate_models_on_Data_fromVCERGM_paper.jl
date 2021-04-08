
# sample sequences of ergms with different parameters' values from R package ergm
# and test the PseudoLikelihoodScoreDrivenERGM filter
using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
DrWatson.greet()

using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets,DynNets, Utilities, ErgmRcall
using LinearAlgebra
using Statistics
using StatsBase
using DataStructures
using JLD2
using StatsBase,CSV, RCall


model = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges +  gwesp(decay = 0.25,fixed = TRUE)", true), indTvPar = [false, true], scoreScalingType="FISH_D")

# load cleaned data
r_data_path = datadir("US_congr_covoting", "clean_data") * "\\Rollcall_VCERGM.RData"
R"load($r_data_path)"
R"obsMat_T = Rollcall$network "
@rget obsMat_T
A_T = [Int8.(m) for m in obsMat_T[2:end]]
obsT = DynNets.seq_of_obs_from_seq_of_mats(model, A_T)

#Single snapshot MLE estimation using MCMC in ergm R package can take a very long time. Load previously estimated data 
# estParSS_T = [ErgmRcall.get_one_mle(m , model.staticModel.ergmTermsString) for m in obsMat_T]
JLD2.@load( datadir("US_congr_covoting", "estimates") * "\\SS_est_MCMC.jld2",estParSS_T)#, obsMat_T, changeStats_T, stats_T)


ENV["JULIA_DEBUG"] = ScoreDrivenERGM

N= [size(A)[1] for A in A_T]

res_est = DynNets.estimate_and_filter(model, N, obsT; show_trace = true)
fig, ax = DynNets.plot_filtered(model, N, res[3])


res_conf = DynNets.conf_bands_given_SD_estimates(model, N, obsT, DynNets.unrestrict_all_par(model, model.indTvPar, res.vEstSdResParAll), res.ftot_0, [[0.975, 0.025]]; indTvPar = model.indTvPar, offset=0, plotFlag=true, parUncMethod = "WHITE-MLE")

DynNets.plot_filtered(model, N, estParSS_T; lineType = ".", lineColor = "k", fig = res_conf.fig, ax=res_conf.ax, gridFlag=false)


estParStatic = ErgmRcall.get_one_mple(obsT, model.staticModel.ergmTermsString)


JLD2.@save( datadir("US_congr_covoting", "estimates")*"\\SD_est_$(DynNets.name(model))", res_conf, res_est, estParStatic)


