#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_\wiki_load_change_stats_and_estimate.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_
#Created Date: Friday April 23rd 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Friday April 23rd 2021 3:39:44 pm
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

model_edge_gwd(decay_par) = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + gwidegree(decay = $decay_par, fixed = TRUE, cutoff=10) + gwodegree(decay = $decay_par, fixed = TRUE, cutoff=10)", true), indTvPar = [true, true, true], scoreScalingType="FISH_D")

model_edge_gwesp = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + gwesp(decay = 0.25, fixed = TRUE, cutoff=10)", true), indTvPar = [true, true], scoreScalingType="FISH_D")

model_edge_mutual_gwd(decay_par) = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + mutual + gwidegree(decay = $decay_par, fixed = TRUE, cutoff=10) + gwodegree(decay = $decay_par, fixed = TRUE, cutoff=10)", true), indTvPar = [true, true, true, true], scoreScalingType="FISH_D")

model_rec_p_star = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + mutual ", true), indTvPar = [true, true], scoreScalingType="FISH_D")
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


@elapsed dfEst = collect_results( datadir("wiki_tk", "ch_stats")) 
viewtab(dfEst[Not(:ch_stats)])

row = dfEst[dfEst.ergmString .==  "edges + mutual + gwidegree(decay = 1.5, fixed = TRUE, cutoff=10) + gwodegree(decay = 1.5, fixed = TRUE, cutoff=10)",:][1,:]

row = dfEst[dfEst.ergmString .==  "edges + gwidegree(decay = 1.5, fixed = TRUE, cutoff=10) + gwodegree(decay = 1.5, fixed = TRUE, cutoff=10)",:][1,:]

row = dfEst[dfEst.ergmString .==  "edges + gwesp(decay = 0.25, fixed = TRUE, cutoff=10)",:][1,:]
model = row.model
T = row.Tf
t0_train = 1
tend_train = 400
N = n_nodes_T[t0_train:tend_train]
obsT = row.ch_stats[t0_train:tend_train]

ENV["JULIA_DEBUG"] = ScoreDrivenERGM

res_est = DynNets.estimate_and_filter(model, N, obsT; show_trace = true)

fig, ax = DynNets.plot_filtered(model, N, res_est[3])

DynNets.number_ergm_par(model)

est_SS = DynNets.estimate_single_snap_sequence(model, obsT)
fig, ax = DynNets.plot_filtered(model, N, est_SS)


res_conf = DynNets.conf_bands_given_SD_estimates(model, N, obsT, DynNets.unrestrict_all_par(model, model.indTvPar, res_est.vEstSdResParAll), res_est.ftot_0, [[0.975, 0.025]]; indTvPar = model.indTvPar, offset=0, plotFlag=true, parUncMethod = "WHITE-MLE", xval = dates[t0_train:tend_train])

using Statistics
var(res_est.fVecT_filt, dims = 2)
DynNets.plot_filtered(model, N, estParSS_T; lineType = ".", lineColor = "k", fig = res_conf.fig, ax=res_conf.ax, gridFlag=false)

