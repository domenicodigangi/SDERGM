#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_\wiki_load_change_stats_and_estimate.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_
#Created Date: Friday April 23rd 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday June 3rd 2021 5:18:23 pm
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
begin
data_path = "./data/collegeMsg/raw_data/"

df = CSV.read(data_path * "collegeMsg_daily_edges.csv", DataFrame)

df.date = Dates.Date.(df.date)

dates = unique(df.date)

# get sequence of edge lists at a given frequency
edge_list_T = [Matrix(gdf[[:source, :target]]) for gdf in groupby(df,[:date])]
edge_list_pres_T = [Utilities.edge_list_pres(l) for l in edge_list_T]

N = length(unique(reshape(Matrix(df[[:source,:target]]), :)))
n_links_T = [size(e)[1] for e in edge_list_T]
n_nodes_T = [ length(unique(reshape(e, :))) for e in edge_list_T]

end

@elapsed dfEst = collect_results( datadir("collegeMsg", "ch_stats_all")) 
# viewtab(dfEst[Not(:ch_stats)])
begin
row = dfEst[1,:]
model = row.model
T = row.Tf
t0_train = 10
tend_train = T

obsT = row.ch_stats[t0_train:tend_train]

ENV["JULIA_DEBUG"] = ScoreDrivenERGM

res_est = DynNets.estimate_and_filter(model, N, obsT; show_trace = true)
fig, ax = DynNets.plot_filtered(model, N, res_est[3])

# res_conf = DynNets.conf_bands_given_SD_estimates(model, N, obsT, DynNets.unrestrict_all_par(model, res_est.vEstSdResParAll), res_est.ftot_0, [[0.975, 0.025]]; , offset=0, plotFlag=true, parUncMethod = "WHITE-MLE", xval = dates[t0_train:tend_train], winsorProp = 0.25)
# res_conf.ax[1].set_title("college msgs $(model.staticModel.ergmTermsString)")

estSS_T = DynNets.estimate_single_snap_sequence(model, obsT)
fig, ax = DynNets.plot_filtered(model, N, estSS_T, ax=ax, lineType = ".", lineColor="r")
end

using Statistics
var(res_est.fVecT_filt[:, 100:end], dims = 2)

DynNets.plot_filtered(model, N, estParSS_T; lineType = ".", lineColor = "k", fig = res_conf.fig, ax=res_conf.ax, gridFlag=false)

dates[100]