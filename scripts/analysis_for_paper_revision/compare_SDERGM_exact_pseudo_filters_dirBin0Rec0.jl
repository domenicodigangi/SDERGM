

"""
Script for comparison of two versions of  GasNetModelDirBin0Rec0, i.e. one parameter for total number of links and one for reciprocity, a version based on MLE and one on PMLE
"""



include("...\\..\\..\\..\\add_load_paths.jl")
using StaticNets:ergm_par_from_mean_vals
using DynNets:fooGasNetModelDirBin0Rec0_mle, sample_est_mle_pmle, dgp_missp, GasNetModelDirBin0Rec0_mle
using HelperFunDom
using PyPlot
import StaticNets: fooNetModelDirBin0Rec0, ergm_par_from_mean_vals
using Statistics
## define parameters

# To Do: 
# 1 - Verify that mean values are actually fixed by inverse formula (best verified via sample averages)
# 2 - Explore candidate regimes statically

# 3 - explore a few regimes and plot results even if they do not make sense:
#     - same scaling sparse dense within physical bounds
#     - different scaling sparse dense within physical bounds
# - Eventually add filter init to first obs
# - filter with same values


pygui(true)
model_mle = fooGasNetModelDirBin0Rec0_mle
staModel = fooNetModelDirBin0Rec0
T =200
## Compare models as misspecified filters, for one network size
Nsample = 2
percAmpl = 0.05
N = 20

denseScale = x-> 0.5 * (x^2-x)/2
sparseScale = x->x
scaleFun = denseScale


N=200
θ_0, η_0 = ergm_par_from_mean_vals(staModel, scaleFun(N), 20, N)
(scaleFun(N), scaleFun(N)/5)
resSin = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), N, Nsample)
# resAR = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "AR"), N, Nsample)
N=20
θ_0, η_0 = ergm_par_from_mean_vals(staModel, scaleFun(N), N, N)
(scaleFun(N), N/2)
resSin = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), N, Nsample)




N = 100
scaleFun = (x->((x^2 -x)/2)/2)
θ_0, η_0 = ergm_par_from_mean_vals(model, scaleFun(N), scaleFun(N)/10, N)
resSin = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), N, Nsample; plotFlag = true)



## Repeat comparison for varying network size

function scaling_comparison(model_mle::GasNetModelDirBin0Rec0_mle, Nvals, scaleFun, Nsample, T; plotFlag = false )
    
    staModel = fooNetModelDirBin0Rec0

    percAmpl = 0.1
    
    linkFun(x) = scaleFun(x)
    recFun(x) = scaleFun(x)/5
    
    rmse_mle = []
    rmse_pmle = []
    for N in Nvals 
        θ_0, η_0 = ergm_par_from_mean_vals(staModel, linkFun(N), recFun(N), N)
        res = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), N, Nsample; plotFlag = plotFlag)
        push!(rmse_pmle, mean(res.rmse_pmle, dims=2))
        push!(rmse_mle, mean(res.rmse_mle, dims=2))
    end

    return rmse_mle, rmse_pmle
end


T = 200
Nvals= round.(Int, collect(20:20:100))
Nsample = 10
denseScale = x-> 0.1 * (x^2-x)/2
sparseScale = x->x

denseRes = scaling_comparison(fooGasNetModelDirBin0Rec0_mle, Nvals, denseScale, Nsample, T, plotFlag = false )

sparseRes = scaling_comparison(fooGasNetModelDirBin0Rec0_mle, Nvals, denseScale, Nsample, T, plotFlase = false )#