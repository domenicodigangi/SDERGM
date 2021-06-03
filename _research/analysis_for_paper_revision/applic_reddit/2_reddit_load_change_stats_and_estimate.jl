#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_\wiki_load_change_stats_and_estimate.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_
#Created Date: Friday April 23rd 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday June 3rd 2021 5:18:21 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########


#%%
using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
DrWatson.greet()

using Distributed

begin
using TableView
using Blink
viewtab(df) = body!(Window(), showtable(df))
using Logging
using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets,DynNets, Utilities, ErgmRcall
using LinearAlgebra
using CSV 
using DataFrames
using Dates
using PyPlot
pygui(true)

end
#%%
#load data
@elapsed dfEst = collect_results( datadir("reddit_hyperlinks", "ch_stats_present")) 
viewtab(dfEst[Not(:ch_stats)])
begin
row = dfEst[dfEst.f .== "w",:][1,:]
T = length(row.ch_stats)
model = row.model
t0_train = 1
tend_train = T
N = row.n_nodes_T[t0_train:tend_train]
obsT = row.ch_stats[t0_train:tend_train]

ENV["JULIA_DEBUG"] = ScoreDrivenERGM


res_est = DynNets.estimate_and_filter(model, N, obsT; show_trace = true)
fig, ax = DynNets.plot_filtered(model, N, res_est[3])

# res_conf = DynNets.conf_bands_given_SD_estimates(model, N, obsT, DynNets.unrestrict_all_par(model,  res_est.vEstSdResParAll), res_est.ftot_0, [[0.975, 0.025]];  offset=0, plotFlag=true, parUncMethod = "WHITE-MLE", winsorProp = 0.25)
# res_conf.ax[1].set_title("college msgs $(model.staticModel.ergmTermsString)")

end
estSS_T = DynNets.estimate_single_snap_sequence(model, obsT)
fig, ax = DynNets.plot_filtered(model, N, estSS_T, ax=ax, lineType = ".", lineColor="r")

using Statistics
var(res_est.fVecT_filt[:, 100:end], dims = 2)

DynNets.plot_filtered(model, N, estParSS_T; lineType = ".", lineColor = "k", fig = res_conf.fig, ax=res_conf.ax, gridFlag=false)

DynNets.plot_filtered(model, row.n_nodes_T, estSS_T; lineType = ".", lineColor = "k")

DynNets.static_estimate(model, obsT)


vEstSdResPar = res_est.vEstSdResParAll
fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter(model, N, obsT,  vEstSdResPar, model.indTvPar;ftot_0 = res_est.ftot_0)