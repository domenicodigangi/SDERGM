#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\add_confidence_intervals\load_estimates_and_check_parametric_bootstrap_conf_intervals.jl
#Project: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\add_confidence_intervals
#Created Date: Saturday April 24th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Tuesday April 27th 2021 10:57:43 am
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

nWorkers = 5
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

# coverFlags = StaticNets.get_coverage_conf_int_par_boot_ergm(model, N, parDgp, quantilesVals, nSampleDgp = 100, nSampleBands=500)

# @show mean(coverFlags, dims=2)

function get_coverage_conf_int_par_boot_ergm_parallel(model, N, parDgp, quantilesVals; nSampleDgp=100, nSampleBands = 100, estimates=nothing)

    # sample from dgp
    if isnothing(estimates)
        matsSample = StaticNets.sample_ergm(model, N, parDgp, nSampleDgp)
        estimates = reduce(hcat,[StaticNets.estimate(model, matsSample[:,:,n]) for n in 1:nSample])
    end
    
    coverFlags =falses(size(estimates))

    coverFlags = @sync @distributed  hcat for n in 1:nSampleDgp
        confInt = StaticNets.get_conf_int_ergm_par_boot(model, N, estimates[:, n], quantilesVals; nSample=nSampleBands)

        [StaticNets.is_between(parDgp[p], confInt[p,:]) for p in 1:model.nErgmPar]
    end

    return coverFlags
end


# coverFlags = get_coverage_conf_int_par_boot_ergm_parallel(model, N, parDgp, quantilesVals)


# @show mean(coverFlags, dims=2)


@elapsed df = collect_results( datadir("sims", "dgp&FIl_est"))
df["modelTag"] = string.(df["model"]) 

T = 300
N = 100
nSample = 50
indQuant = 1
modelTag = string(DynNets.SdErgmDirBin0Rec0_mle(scoreScalingType = "FISH_D"))

dgpSetting = DynNets.list_example_dgp_settings(DynNets.SdErgmDirBin0Rec0_mle()).dgpSetARlowlow

# # visual inspection of filter, dgp and conf bands
res = filter([:modelTag, :T, :N, :dgpSettings, :nSample, ] => (mTag,t,n, d, s) -> all((mTag==modelTag, t==T, n==N, d == dgpSetting, s == nSample)), df)[1,:]


model = res.model
using Statistics
n = 2

# coverFlags = @sync @distributed  hcat
nQuant = length(quantilesVals)
confBandsT = zeros(model.staticModel.nErgmPar, T, nSample, nQuant)
for t = 1:T
parDgp_t = res.allParDgpT[:,t,n]
parEst_t = res.allfVecT_filt[:,t,n]
confBandsT[:,t , n,:] = mapslices((x -> quantile(x, quantilesVals)), StaticNets.get_par_boot_ergm_distrib(model.staticModel, parEst_t, N; nSample = 100), dims=2)

end


mat_sample = StaticNets.sample_ergm(model, N, parEst, nSample)


p=2
plot(confBandsT[p,:,n,:])
plot(res.allParDgpT[p,:,n])
plot(res.allfVecT_filt[p,:,n])


#%% 
## Sequence for rescaled Par Bootstrap confidence intervals for SDERGM


#Estimate SD on data and filter



# Obtain Distribution of observations from the filtered path, sampling each ergm_t using the filtered parameters nSample times



# filter each 


