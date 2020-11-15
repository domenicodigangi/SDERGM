

"""
Script for comparison of two versions of  GasNetModelDirBin0Rec0, i.e. one parameter for total number of links and one for reciprocity, a version based on MLE and one on PMLE
"""



include("...\\..\\..\\..\\add_load_paths.jl")
using StaticNets
using DynNets
using HelperFunDom
using PyPlot


## define parameters
θ_0 = -0.5
η_0 = 0.5
T =100


model_mle = fooGasNetModelDirBin0Rec0_mle




## Compare models as misspecified filters, for one network size
Nsample = 2
percAmpl = 50/100

resAR = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "AR"), 30, Nsample)
resSin = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), 60, Nsample)



## Repeat comparison for varying network size

Nvals= round.(Int, collect(range(20,stop=100,length=10)))
rmse_mle = []
rmse_pmle = []
for N in Nvals 
    res = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), N, Nsample)
    push!(rmse_pmle, mean(res.rmse_pmle, dims=2))
    push!(rmse_mle, mean(res.rmse_mle, dims=2))

end




