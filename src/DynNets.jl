module DynNets


using Distributions
using StatsBase
using Optim
using LineSearches
using StatsFuns
using Statistics
using Roots
using MLBase
using GLM
using LinearAlgebra
using JLD
using DataFrames
using ForwardDiff
using NLSolversBase
using RCall
using PyPlot
using Logging


using ..AReg
using ..Utilities
using ..Scalings
using ..ErgmRcall
using ..StaticNets


abstract type GasNetModel end
abstract type GasNetModelW <: GasNetModel end
abstract type GasNetModelWcount <: GasNetModelW end
abstract type GasNetModelBin <: GasNetModel end
abstract type GasNetModelDirBin0Rec0 <: GasNetModel end


#constants
targetErrValDynNets = 0.01

#Relations between Static and Dynamic Models]
identify(Model::GasNetModel,UnPar::Array{<:Real,1};idType = "pinco") =
    StaticNets.identify(StaModType(Model),UnPar;idType = idType)



include("./DynNets_models/DynNets_SDModelsUtils.jl")

include("./DynNets_models/DynNets_GasNetModelBin1.jl")

include("./DynNets_models/DynNets_GasNetModelDirBin1.jl")

include("./DynNets_models/DynNets_GasNetModelDirBin0Rec0.jl")

include("./DynNets_models/DynNets_GasNetModelDirBinERGM.jl")

include("./DynNets_models/DynNets_DirBinGlobalPseudo.jl")

include("./DynNets_models/DynNets_paper_helper_funs.jl")


end



