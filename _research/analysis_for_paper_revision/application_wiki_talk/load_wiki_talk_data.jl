#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_\load_wiki_talk_data.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\new_application_
#Created Date: Wednesday April 21st 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Saturday April 24th 2021 4:26:42 pm
#Modified By:  Domenico Di Gangi
#-----
#Description:
#-----
########


#%%
using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"
DrWatson.greet()

using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets,DynNets, Utilities, ErgmRcall
using LinearAlgebra
using CSV 
using DataFrames
using Dates
using PyPlot
pygui(true)

#%%



data_path = "./data/wiki_tk/raw_data/"

df = CSV.read(data_path * "wiki_talk_daily_edges.csv", DataFrame)

df.date = Dates.Date.(df.date)

dates = unique(df.date)

# get sequence of edge lists

edge_list_T = [Matrix(gdf[[:source, :target]]) for gdf in groupby(df,[:date])]
edge_list_pres_T = [edge_list_pres(l) for l in edge_list_T]


n_links_T = [size(e)[1] for e in edge_list_T]
n_nodes_T = [ length(unique(reshape(e, :))) for e in edge_list_T]


# t_ref = argmax(n_nodes_T)
# n_links_sub = [100, 500, 1000, 2000, 3000, 5000, 10000]
# time = [@time DynNets.stats_from_mat(model1, edge_list_pres(edge_list_pres_T[t_ref][1:n_links, :])) for n_links in n_links_sub]
# plot(n_links_sub, time)
# time = [@time DynNets.stats_from_mat(model_rec_p_star, edge_list_pres(edge_list_pres_T[t_ref][1:n_links, :])) for n_links in n_links_sub]


using Distributed
begin
nWorkers = 10
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 
    Pkg.instantiate() 
    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:SdErgm,SdErgmDirBin0Rec0, simulate_and_estimate_parallel
    using ScoreDrivenERGM.Utilities

    gwd_decay = 2
    model_gwd = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + mutual + gwidegree(decay = $gwd_decay,fixed = TRUE, cutoff=10) + gwodegree(decay = $gwd_decay,fixed = TRUE, cutoff=10)", true), indTvPar = [true, true], scoreScalingType="FISH_D")
    model_rec_p_star = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + mutual ", true), indTvPar = [true, true], scoreScalingType="FISH_D")

    end
end


T = length(n_nodes_T)

t0_ch_stats = 1
tend_ch_stats = T
ch_stats_rec_p_star = @sync @distributed vcat for edl in edge_list_pres_T[t0_ch_stats:tend_ch_stats]
    [DynNets.stats_from_mat(model_rec_p_star, edge_list_pres(edl))]
end

model = model_rec_p_star
t0_train = 500
tend_train = T-50
N = n_nodes_T[t0_train:tend_train]
obsT = ch_stats_rec_p_star[t0_train:tend_train]


res_est = DynNets.estimate_and_filter(model, N, obsT; show_trace = true)
fig, ax = DynNets.plot_filtered(model, N, res_est[3])


res_conf = DynNets.conf_bands_given_SD_estimates(model, N, obsT, DynNets.unrestrict_all_par(model, model.indTvPar, res_est.vEstSdResParAll), res_est.ftot_0, [[0.975, 0.025]]; indTvPar = model.indTvPar, offset=0, plotFlag=true, parUncMethod = "WHITE-MLE", xval = dates[t0_train:tend_train])

using Statistics
var(res_est.fVecT_filt, dims = 2)
DynNets.plot_filtered(model, N, estParSS_T; lineType = ".", lineColor = "k", fig = res_conf.fig, ax=res_conf.ax, gridFlag=false)
