
__precompile__(false)
module DynNets
#export sample_bern_mat, sample_fitTV, sampleDgpNetsTv

using Distributions, StatsBase,Optim, LineSearches, StatsFuns, Roots,MLBase, GLM, LinearAlgebra, JLD,DataFrames, ForwardDiff,NLSolversBase, RCall
# using PyCall;# pygui(:qt);
# using PyPlot

abstract type GasNetModel end
abstract type GasNetModelW <: GasNetModel end
abstract type GasNetModelWcount <: GasNetModelW end
abstract type GasNetModelBin <: GasNetModel end

#constants
targetErrValDynNets = 0.01

#Relations between Static and Dynamic Models]
identify(Model::GasNetModel,UnPar::Array{<:Real,1};idType = "pinco") =
    StaticNets.identify(StaModType(Model),UnPar;idType = idType)


include("./AReg.jl")
include("./HelperFunDom.jl")
include("./StaticNets.jl")


include("./DynNets_GasNetModelBin1.jl")
include("./DynNets_GasNetModelDirBin1.jl")
include("./DynNets_DirBinGlobalPseudo.jl")
include("./DynNets_GasNetModelDirBin0Rec0.jl")

end



