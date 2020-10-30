

"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""



include("...\\..\\..\\..\\add_load_paths.jl")
using StaticNets
using DynNets
using HelperFunDom
using PyPlot

sn = StaticNets
hfd = HelperFunDom


## define parameters
N=30
θ_0 = 1
η_0 = -1.5
T =100


model_mle = fooGasNetModelDirBin0Rec0_mle
model_pmle = fooGasNetModelDirBin0Rec0_pmle
indTvPar = trues(2)
B0=0.95
A0 = 0.01
vResGasPar_0 = DynNets.array2VecGasPar(model_mle, [[θ_0*(1-B0), B0, A0], [η_0*(1-B0), B0, A0] ], indTvPar)
vResGasPar = vResGasPar_0
ftot_0= zeros(Real,2)




## define tv par dgp
Nsample = 50
dgpType  ="AR"
parMatDgp_T = dgp_missp(model_mle, θ_0, η_0, percAmpl, dgpType)
resAR = sample_est_mle_pmle(model_mle, parMatDgp_T, N, Nsample)


dgpType  ="sin"
parMatDgp_T = dgp_missp(model_mle, θ_0, η_0, percAmpl, dgpType)
resSin = sample_est_mle_pmle(model_mle, parMatDgp_T, N, Nsample)







 