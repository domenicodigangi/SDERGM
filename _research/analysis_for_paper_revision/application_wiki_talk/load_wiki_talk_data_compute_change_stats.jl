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
using Distributed
begin
using Logging
using ScoreDrivenERGM
import ScoreDrivenERGM:StaticNets,DynNets, Utilities, ErgmRcall
using LinearAlgebra
using CSV 
using DataFrames
using Dates
using PyPlot
pygui(true)

nWorkers = 15
addprocs(nWorkers - nprocs())
@sync @everywhere begin 
    using Pkg
    Pkg.activate(".") 
    Pkg.instantiate() 
    using ScoreDrivenERGM
    import ScoreDrivenERGM:StaticNets, DynNets
    import ScoreDrivenERGM.DynNets:SdErgm,SdErgmDirBin0Rec0, simulate_and_estimate_parallel
    using ScoreDrivenERGM.Utilities

    model_edge_gwd(decay_par) = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + gwidegree(decay = $decay_par, fixed = TRUE, cutoff=10) + gwodegree(decay = $decay_par, fixed = TRUE, cutoff=10)", true), indTvPar = [true, true, true], scoreScalingType="FISH_D")
    
    model_edge_gwesp = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + gwesp(decay = 0.25, fixed = TRUE, cutoff=10)", true), indTvPar = [true, true], scoreScalingType="FISH_D")

    model_edge_mutual_gwd(decay_par) = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + mutual + gwidegree(decay = $decay_par, fixed = TRUE, cutoff=10) + gwodegree(decay = $decay_par, fixed = TRUE, cutoff=10)", true), indTvPar = [true, true, true, true], scoreScalingType="FISH_D")

    model_rec_p_star = DynNets.SdErgmPml(staticModel = StaticNets.NetModeErgmPml("edges + mutual ", true), indTvPar = [true, true], scoreScalingType="FISH_D")



    """ return an edge list with nodes indexed starting from 1 to the number of nodes present in the list """
    function edge_list_pres(list)
        pres = sort(unique(reshape(list, :)))
        list_pres =  [findfirst(isequal(n), pres) for n in list]
        sorted_list_pres = sortslices(list_pres, dims=1)
        return sorted_list_pres
    end

    end
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
edge_list_pres_T = [edge_list_pres(l) for l in edge_list_T]

n_links_T = [size(e)[1] for e in edge_list_T]
n_nodes_T = [ length(unique(reshape(e, :))) for e in edge_list_T]

end

# t_ref = argmax(n_nodes_T)
# n_links_sub = [100, 500, 1000, 2000, 3000, 5000, 10000]
# time = [@time DynNets.stats_from_mat(model1, edge_list_pres(edge_list_pres_T[t_ref][1:n_links, :])) for n_links in n_links_sub]
# plot(n_links_sub, time)
# time = [@time DynNets.stats_from_mat(model_rec_p_star, edge_list_pres(edge_list_pres_T[t_ref][1:n_links, :])) for n_links in n_links_sub]

begin
T = length(n_nodes_T)
c= Dict{String, Any}()
c["model"] =[model_rec_p_star, model_edge_gwd.([0, 0.5, 0.8, 1.2, 1.5, 2])..., model_edge_mutual_gwd.([0, 0.5, 0.8, 1.2, 1.5, 2])..., model_edge_gwesp] 
c["Ti"] = [1]
c["Tf"] = [T]
c["f"] = ["d"]
list = sort(dict_list(c), by=x->(string(x["model"])))
d = list[1]
for d in list

    timeSim = @elapsed ch_stats = @sync @distributed vcat for edl in edge_list_pres_T[d["Ti"]:d["Tf"]]
        [DynNets.stats_from_mat(d["model"], edge_list_pres(edl))]
    end

                
    modelTag =  string(d["model"])
    ergmString = d["model"].staticModel.ergmTermsString

    res1 =  (;modelTag, ergmString, ch_stats) |>DrWatson.ntuple2dict |> DrWatson.tostringdict

    saveDict = merge(res1, d)

    d_save = deepcopy(d)
    d_save["ter"] = replace.( ergmString,  ", fixed = TRUE, cutoff=10" => "")

    saveName = replace.( savename(d_save, "jld2";allowedtypes = (Real, String, Symbol, NamedTuple, Tuple, ScoreDrivenERGM.DynNets.SdErgm) ), r"[\"]" => "")

    timeSave = @elapsed save( datadir("wiki_tk", "ch_stats", saveName), saveDict)

    Logging.@info("Time sim = $timeSim ,  time save = $timeSave ")


end

end



