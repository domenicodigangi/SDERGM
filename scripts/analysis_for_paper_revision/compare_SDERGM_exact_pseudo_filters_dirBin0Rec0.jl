

"""
Script for comparison of two versions of  GasNetModelDirBin0Rec0, i.e. one parameter for total number of links and one for reciprocity, a version based on MLE and one on PMLE

To Do: 
1 - Verify that mean values are actually fixed by inverse formula (best verified via sample averages) DONE
2 - Correct sin dgp to oscillate between min and max values that are obtained from average mean values

3 - explore a few regimes and plot results even if they do not make sense:
    - same scaling sparse dense within physical bounds
    - different scaling sparse dense within physical bounds
- Eventually add filter init to first obs
- filter with same values

"""



include("...\\..\\..\\..\\add_load_paths.jl")
using StaticNets:ergm_par_from_mean_vals
using DynNets:fooGasNetModelDirBin0Rec0_mle, sample_est_mle_pmle, dgp_missp, GasNetModelDirBin0Rec0_mle
using HelperFunDom
import StaticNets: fooNetModelDirBin0Rec0, ergm_par_from_mean_vals
using Statistics
using PyPlot
pygui(true)



model_mle = fooGasNetModelDirBin0Rec0_mle
staModel = fooNetModelDirBin0Rec0
T =200
## Compare models as misspecified filters, for one network size
Nsample = 2


"""
Number of links for various network size, in dense regime
"""
denseLScal(N)  = 0.1 * (N^2-N)
"""
Number of links for various network size, in sparse regime
"""
sparseLScal(N) =  1 * N 

"""
Average number of reciprocal pairs in Erdos Reny Model 
 - avgL is the expected number of links, related to the probability p of a single link being present by expL = (N^2-N) * p  

"""
erdosRenyRecScal(expL, N) = expL^2/(2*(N^2-N)) 


linkScal(N) = sparseLScal(N)
recScal(N) =  10 * erdosRenyRecScal(linkScal(N), N)


N=100
percAmpl = 0.1
θ_0, η_0 = ergm_par_from_mean_vals(staModel, linkScal(N), recScal(N), N)
(linkScal(N), recScal(N))
resSin = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "sin"), N, Nsample)
# resAR = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0, η_0, percAmpl, "AR"), N, Nsample)
N=20
percAmpl = 0.05
θ_0, η_0 = ergm_par_from_mean_vals(staModel, linkScal(N), recScal(N), N)
(linkScal(N), recScal(N))
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
Nsample
"""
Number of links for various network size, in dense regime
""" = 10
denseLScal = x-> 0.1 * (x^2-x)/2
"""
Number of links for various network size, in sparse regime
""" = 10
sparseLScal = x->x

"""
Average number of reciprocal pairs in Erdos Reny Model 
 - avgL is the expected number of links, related to the probability of a single link by expL = (N^2-N) * p  
"""

denseRes = expling_comparison(fooGasNetModelDirBin0R
avgL, "" = "
Number of links for various network size, in dense regime
"""ec0_mle, Nvals, denseLinkScal, Nsample, T, pl
"""
"""
Number of links for various network size, in sparse regime
"""ec0_mle, Nvals, denseLinkScal, Nsample, T, pl
Number of links for various network size, in dense regime
"""otFlag = false
"""
Average number of reciprocal pairs in Erdos Reny Model 
"""
 - avgL is the average number of links, expected to  

"""
Number of the probability of a single link by expL = (N^2-N) * p links for various network size, expsparse avgL, re = gime
"""otFlag = false
sparsees = scaling_comparison(fooGasN
"""
 - avgL is the average number of expected, related to  

Average the probability of a single link by expL = (N^2-N) * p number of reciprocal pairs in Erdos Reny Model 
"""etModelDirBin00le, Nvals, expseLinkScal, Nsample, T, plotFlase = faavgL, ls = e