#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\add_confidence_intervals\load_estimates_and_check_parametric_bootstrap_conf_intervals.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\add_confidence_intervals
#Created Date: Saturday April 24th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Saturday April 24th 2021 7:16:20 pm
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

using DataFrames
using ScoreDrivenERGM
import ScoreDrivenERGM:DynNets, Utilities, StaticNets, Scalings
using PyPlot
pygui(true)
#endregion

using Distributed

nWorkers = 10
addprocs(nWorkers - nprocs())
@everywhere begin
    import ScoreDrivenERGM:DynNets, Utilities, StaticNets, Scalings
end


## check coverage of confidence intervals in static ergm ~ edges + mutual
model = StaticNets.ErgmDirBin0Rec0()
@everywhere begin
N=500
Ldgp = Scalings.n_pox_dir_links(N)/4
Rdgp = Scalings.n_pox_pairs(N)/4-100
parDgp = StaticNets.estimate(model, Ldgp, Rdgp, N)
quantilesVals = [0.025, 0.975]
end

coverFlags = StaticNets.get_coverage_conf_int_par_boot_ergm(model, N, parDgp, quantilesVals, nSampleDgp = 100, nSampleBands=500)

@show mean(coverFlags, dims=2)

function get_coverage_conf_int_par_boot_ergm_parallel(model, N, parDgp, quantilesVals; nSampleDgp=100, nSampleBands = 100)

    # sample from dgp
    matsSample = StaticNets.sample_ergm(model, N, parDgp, nSampleDgp)

    estimates = reduce(hcat,[StaticNets.estimate(model, A) for A in matsSample])

    coverFlags =falses(size(estimates))

    coverFlags = @sync @distributed  hcat for n in 1:nSampleDgp
        confInt = StaticNets.get_conf_int_ergm_par_boot(model, N, estimates[:, n], quantilesVals; nSample=nSampleBands)

        [StaticNets.is_between(parDgp[p], confInt[p,:]) for p in 1:model.nErgmPar]
    end

    return coverFlags
end

coverFlags = get_coverage_conf_int_par_boot_ergm_parallel(model, N, parDgp, quantilesVals)

@show mean(coverFlags, dims=2)


@elapsed df = collect_results( datadir("sims", "dgp&FIl_est"))

df["modelTag"] = string.(df["model"]) 
fieldnames(typeof(df.model[1]))



using JLD2
allData = [load(datadir("sims", "dgp&FIl_conf") *"//" * file) for file in readdir(datadir("sims", "dgp&FIl_conf"))]


for f in  readdir(datadir("sims", "dgp&FIl_conf"))
    data = load(datadir("sims", "dgp&FIl_conf") *"//" * f)
    isPseudo = (s -> contains(s,"pmle") )(string(data["model"]))
        @show isPseudo
    if isPseudo
        model = DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType="FISH_D")
    else
        model = DynNets.SdErgmDirBin0Rec0_pmle(scoreScalingType="FISH_D")
    end
    data["model"] = model

    save( datadir("sims", "dgp&Fil_conf", f), data)

end

allData

DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType="FISH_D")

df.model[1].staticModel


T = 100
N = 100
nSample = 50
indQuant = 1
modelTag = string(DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType = "FISH_D"))

dgpSetting = DynNets.list_example_dgp_settings(DynNets.SdErgmDirBin0Rec0_mle()).dgpSetARlowlow

# # visual inspection of filter, dgp and conf bands
res = filter([:modelTag, :T, :N, :dgpSettings, :nSample, ] => (mTag,t,n, d, s) -> all((mTag==modelTag, t==T, n==N, d == dgpSetting, s == nSample)), df)

nrow(res) != 1 ? error("$N, $T,  $(size(res))") : res = res[1,:]
    for nPlot = 1:10
        if !any(res.errInds[:,nPlot])
            DynNets.plot_filtered_and_conf_bands(res.model, N, res.allfVecT_filt[:,1:end,nPlot], res.allConfBandsFiltPar[:,1:end,:,:,nPlot] ; parDgpTIn=res.allParDgpT[:,1:T,nPlot], confBands2In = res.allConfBandsPar[:,1:end,:,:,nPlot] , offset = 1, indBand = indQuant)
        end
    end

using Statistics









