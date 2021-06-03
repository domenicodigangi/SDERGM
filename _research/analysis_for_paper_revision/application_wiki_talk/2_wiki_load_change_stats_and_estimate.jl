#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_\wiki_load_change_stats_and_estimate.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_
#Created Date: Friday April 23rd 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday June 3rd 2021 5:18:20 pm
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
data_path = "./data/wiki_tk/raw_data/"

df = CSV.read(data_path * "wiki_talk_daily_edges.csv", DataFrame)

df.date = Dates.Date.(df.date)

dates = unique(df.date)

# get sequence of edge lists at a given frequency
edge_list_T = [Matrix(gdf[[:source, :target]]) for gdf in groupby(df,[:date])]
edge_list_pres_T = [Utilities.edge_list_pres(l) for l in edge_list_T]

n_links_T = [size(e)[1] for e in edge_list_T]
n_nodes_T = [ length(unique(reshape(e, :))) for e in edge_list_T]

end

# t_ref = argmax(n_nodes_T)
# n_links_sub = [100, 500, 1000, 2000, 3000, 5000, 10000]
# time = [@time DynNets.stats_from_mat(model1, edge_list_pres(edge_list_pres_T[t_ref][1:n_links, :])) for n_links in n_links_sub]
# plot(n_links_sub, time)
# time = [@time DynNets.stats_from_mat(model_rec_p_star, edge_list_pres(edge_list_pres_T[t_ref][1:n_links, :])) for n_links in n_links_sub]


@elapsed dfEst = collect_results( datadir("wiki_tk", "ch_stats_present")) 
viewtab(dfEst[Not(:ch_stats)])

row = dfEst[dfEst.ergmString .==  "edges + gwidegree(decay = 1.5, fixed = TRUE, cutoff=10) + gwodegree(decay = 1.5, fixed = TRUE, cutoff=10)",:][1,:]


row = dfEst[dfEst.ergmString .==  "edges + mutual",:][1,:]

row = dfEst[dfEst.ergmString .==  "edges + gwesp(decay = 0.25, fixed = TRUE, cutoff=10)",:][1,:]
begin
row = dfEst[dfEst.ergmString .==  "edges + mutual + gwidegree(decay = 1.5, fixed = TRUE, cutoff=10) + gwodegree(decay = 1.5, fixed = TRUE, cutoff=10)",:][1,:]

model = row.model
T = row.Tf
t0_train = 1
tend_train = 365
N = n_nodes_T[t0_train:tend_train]
obsT = row.ch_stats[t0_train:tend_train]

ENV["JULIA_DEBUG"] = ScoreDrivenERGM

res_est = DynNets.estimate_and_filter(model, N, obsT; show_trace = true)
fig, ax = DynNets.plot_filtered(model, N, res_est[3])
end

# res_conf = DynNets.conf_bands_given_SD_estimates(model, N, obsT, DynNets.unrestrict_all_par(model, res_est.vEstSdResParAll), res_est.ftot_0, [[0.975, 0.025]];  offset=0, plotFlag=true, parUncMethod = "WHITE-MLE", xval = dates[t0_train:tend_train], winsorProp = 0)
# res_conf.ax[1].set_title("wiki talk $(model.staticModel.ergmTermsString)")

est_SS = DynNets.estimate_single_snap_sequence(model, obsT)
# fig, ax = DynNets.plot_filtered(model, N, est_SS, ax=res_conf.ax, lineType = ".", lineColor="r")

dates[365]
using Statistics
var(res_est.fVecT_filt[:, 1:365], dims = 2)

