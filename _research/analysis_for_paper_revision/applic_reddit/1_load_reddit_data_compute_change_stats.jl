#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\applic_reddit\1_load_reddit_data_compute_change_stats.jl
#Created Date: Tuesday May 4th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Friday May 7th 2021 6:09:42 pm
#Modified By:  Domenico Di Gangi
#-----
#Description: Load reddit edges and save change statistics for a set of models
#-----
########




#%%

using DrWatson
@quickactivate "ScoreDrivenExponentialRandomGraphs"

begin
using Distributed
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

    end
end

#%%
#load data


begin
data_path = "./data/reddit_hyperlinks/raw_data/"
df = CSV.read(data_path * "reddit_hyperlinks_edges.csv", DataFrame)
df.date = Dates.Date.(df.date)
df.week = week.(df.date)
df.year = year.(df.date)
end
 
begin
c= Dict{String, Any}()
c["model"] =[DynNets.model_edge_mutual_gwesp] 
c["f"] = ["w"]
list = sort(dict_list(c), by=x->(string(x["model"])))
d = list[1]
only_present_nodes = true
for d in list
    if d["f"] == "w"
        groupCols = [:year, :week]
    elseif d["f"] == "d"
        groupCols = [:date]
    end
    # get sequence of edge lists at a given frequency
    edge_list_T = [Matrix(gdf[[:source, :target]]) for gdf in groupby(df,groupCols)]
    edge_list_pres_T = [Utilities.edge_list_pres(l) for l in edge_list_T]

    n_links_T = [size(e)[1] for e in edge_list_T]
    n_nodes_T = [ length(unique(reshape(e, :))) for e in edge_list_T]
    dates = unique(df.date)
    figure()
    plot(n_links_T)
    plot(n_nodes_T)

    if only_present_nodes
        timeSim = @elapsed ch_stats = @sync @distributed vcat for edl in edge_list_pres_T
            [DynNets.stats_from_mat(d["model"], Utilities.edge_list_pres(edl))]
        end
    else
        timeSim = @elapsed ch_stats = @sync @distributed vcat for edl in edge_list_T
            [DynNets.stats_from_mat(d["model"], edl)]
        end
    end

                
    modelTag =  string(d["model"])
    ergmString = d["model"].staticModel.ergmTermsString

    res1 =  (;modelTag, ergmString, ch_stats, n_links_T, n_nodes_T) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    saveDict = merge(res1, d)

    d_save = deepcopy(d)
    d_save["ter"] = replace.( ergmString,  ", fixed = TRUE, cutoff=10" => "")

    saveName = replace.( savename(d_save, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.SdErgm) ), r"[\"]" => "")

    only_present_nodes ? nodes_subset =  "present" : nodes_subset = "all"

    timeSave = @elapsed save( datadir("reddit_hyperlinks", "ch_stats_$nodes_subset", saveName), saveDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")


end

end



