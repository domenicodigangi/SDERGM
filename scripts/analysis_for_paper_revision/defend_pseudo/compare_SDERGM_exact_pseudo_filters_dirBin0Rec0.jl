

"""
Script for comparison of two versions of  GasNetModelDirBin0Rec0, i.e. one parameter for total number of links and one for reciprocity, a version based on MLE and one on PMLE

To Do: 

"""



using StaticNets:ergm_par_from_mean_vals
using DynNets:fooGasNetModelDirBin0Rec0_mle, sample_est_mle_pmle, dgp_missp, GasNetModelDirBin0Rec0_mle
using Utilities
import StaticNets: fooNetModelDirBin0Rec0, ergm_par_from_mean_vals,exp_val_stats, diadProbFromPars, samplSingMatCan, statsFromMat
using Statistics
using PyPlot
pygui(true)



model_mle = fooGasNetModelDirBin0Rec0_mle
staModel = fooNetModelDirBin0Rec0



"""
Number of links for various network size, in dense regime
"""
denseLScal(N)  = 0.2 * n_pox_dir_links(N)

"""
Number of links for various network size, in semiDense regime
"""
semiDenseLScal(N)  = 0.1 * n_pox_dir_links(N)/sqrt(N)

"""
Number of links for various network size, in sparse regime
"""
sparseLScal(N) =  4 * N 

"""
Average number of reciprocal pairs in Erdos Reny Model 
 - avgL is the expected number of links, related to the probability p of a single link being present by expL = (N^2-N) * p  

"""
erdosRenyRecScal(expL, N) = expL^2/(2*(N^2-N)) 

## Search an appropriate dgp, i.e. one that stays away form physical bounds



function scaled_sin_sample_est_mle_pmle(alphaMeanVal::Function, N, T, Nsample; regimeString ="", plotDgpOrigin=false, plotFlag=false)
    
    betaMeanVal(N) = alphaMeanVal(N)/5 

    amplAlpha = 0.3
    amplBeta = 0.5

    nCycles = 2

    get_alpha_beta_seq(N) = hcat(dgpSin( (1-amplAlpha/2)*alphaMeanVal(N), (1+amplAlpha/2)*alphaMeanVal(N), nCycles, T; phase = 10*rand()), dgpSin( (1-amplBeta/2)*betaMeanVal(N), (1+amplBeta/2)*betaMeanVal(N), nCycles, T; phase = 10*rand()))



    alpha_beta_to_theta_eta(α, β, N) = collect(ergm_par_from_mean_vals(staModel, α*n_pox_dir_links(N), β*n_pox_dir_links(N), N))

    theta_eta_to_alpha_beta(θ, η, N) =  collect(exp_val_stats(staModel, θ, η, N))./n_pox_dir_links(N)

    get_theta_eta_seq(N) = reduce(hcat, [alpha_beta_to_theta_eta(ab[1], ab[2], N) for ab in eachrow(get_alpha_beta_seq(N))])


    if plotDgpOrigin
        startAlphaBetaSeq = get_alpha_beta_seq(N)
        thetaEtaSeq = permutedims(reduce(hcat, [alpha_beta_to_theta_eta(ab[1], ab[2], N) for ab in eachrow(startAlphaBetaSeq)]))

        fig, ax = subplots(2,2) 
        ax[1,1].plot(startAlphaBetaSeq[:,1])
        ax[1,1].set_title("Alpha Dgp")
        ax[2,1].plot(startAlphaBetaSeq[:,2])
        ax[2,1].set_title("Beta Dgp")
        ax[1,2].plot(thetaEtaSeq[:,1])
        ax[1,2].set_title("Theta Dgp")
        ax[2,2].plot(thetaEtaSeq[:,2])
        ax[2,2].set_title("Eta Dgp")
        fig.tight_layout()
        fig.suptitle("$regimeString Dgp for N = $N ")
    end
        
    # In case we want to add a check before filtering, start from this commented lines
    # minNPairs = 10
    # lBoundBeta = minNPairs/n_pox_dir_links(N)
    # uBoundBeta = alphaMeanVal(N)/2 -  minNPairs/n_pox_dir_links(N)
    # startLRSeq = get_alpha_beta_seq(N).*n_pox_dir_links(N)
    # endAlphaBetaSeq = permutedims( reduce(hcat, [theta_eta_to_alpha_beta(θ, η, N) for (θ, η) in eachcol(get_theta_eta_seq(N))]) )
    # θ_0, η_0 = get_theta_eta_seq(N)[:,1]
    # A_vec = [statsFromMat(staModel, samplSingMatCan(staModel, diadProbFromPars(staModel, [θ_0, η_0]), N)) for i=1:100]

    return sample_est_mle_pmle(model_mle, get_theta_eta_seq(N), N, Nsample; plotFlag=plotFlag, regimeString=regimeString)
end

# Check that misspecified filters run without evident issues on few samples at N extrema

denseAlphaScal(N) = denseLScal(N)/n_pox_dir_links(N)
sparseAlphaScal(N) = sparseLScal(N)/n_pox_dir_links(N)
T = 200
scaled_sin_sample_est_mle_pmle(denseAlphaScal, 10, T, 2, regimeString = "Dense", plotDgpOrigin=true)
scaled_sin_sample_est_mle_pmle(sparseAlphaScal, 10, T, 2, regimeString = "Sparse", plotDgpOrigin=true)

scaled_sin_sample_est_mle_pmle(denseAlphaScal, 100, T, 2,regimeString = "Dense", plotDgpOrigin=true)
scaled_sin_sample_est_mle_pmle(sparseAlphaScal, 100, T, 2, regimeString = "Sparse", plotDgpOrigin=true)


"""
for a single misspecified dgp compare mle and mple rmse for varying network sizes 
"""
function scaling_comparison(model_mle::GasNetModelDirBin0Rec0_mle, Nvals, alphaScal, Nsample, T; plotFlag = false, regimeString="")
    
    staModel = fooNetModelDirBin0Rec0

    rmse_mle = []
    rmse_pmle = []
    mean_vals = []
    
    for N in Nvals 

        res = scaled_sin_sample_est_mle_pmle(alphaScal, N, T, Nsample, regimeString=regimeString; plotFlag=plotFlag)

        push!(rmse_pmle, res.rmse_pmle)
        push!(rmse_mle, res.rmse_mle)

    end

    return rmse_mle, rmse_pmle, mean_vals
end



using Serialization

T = 200
# Nvals= [20, 100]# round.(Int, collect(20:20:100))
Nvals=  round.(Int, collect(20:10:100))
Nsample = 50

savePath = (@__DIR__ )* "\\revision_data"

scaling_res = scaling_comparison(model_mle, Nvals, sparseAlphaScal, Nsample, T; plotFlag = false, regimeString = "Sparse")
serialize(savePath * "\\sparse_Nsample_$(Nsample)_precise_filter_init.jls", scaling_res)

scaling_res =scaling_comparison(model_mle, Nvals, denseAlphaScal, Nsample, T; plotFlag = false, regimeString = "Dense")
serialize(savePath * "\\dense_Nsample_$(Nsample)_precise_filter_init.jls", scaling_res)




Nsample = 50
rmse_mle_sparse, rmse_pmle_sparse, mean_vals_sparse = deserialize(savePath * "\\sparse_Nsample_$(Nsample)_precise_filter_init.jls")

rmse_mle_dense, rmse_pmle_dense, mean_vals_dense = deserialize(savePath * "\\dense_Nsample_$(Nsample)_precise_filter_init.jls")


fig, ax = subplots(2,2)
ax[1,1].plot(Nvals, reduce(hcat,rmse_mle_sparse)[1,:], "-.b")
ax[1,1].plot(Nvals, reduce(hcat,rmse_pmle_sparse)[1,:], "-.r")
ax[2,1].plot(Nvals, reduce(hcat,rmse_mle_sparse)[2,:], "-.b")
ax[2,1].plot(Nvals, reduce(hcat,rmse_pmle_sparse)[2,:], "-.r")
ax[1,1].set_title("Sparse regime")

ax[1,2].plot(Nvals, reduce(hcat,rmse_mle_dense)[1,:], "-.b")
ax[1,2].plot(Nvals, reduce(hcat,rmse_pmle_dense)[1,:], "-.r")
ax[2,2].plot(Nvals, reduce(hcat,rmse_mle_dense)[2,:], "-.b")
ax[2,2].plot(Nvals, reduce(hcat,rmse_pmle_dense)[2,:], "-.r")
ax[1,2].set_title("Dense  regime")



# N=80
# linkScal(N) = semiDenseLScal(N)
# recScal(N) = min(5 + 10 * erdosRenyRecScal(linkScal(N), N), linkScal(N)/2 -1)
# percAmpl = 0.5
# θ_0_min, η_0_min = ergm_par_from_mean_vals(staModel, linkScal(N)*(1-percAmpl), recScal(N)*(1-percAmpl), N)
# θ_0_max, η_0_max = ergm_par_from_mean_vals(staModel, linkScal(N)*(1+percAmpl), recScal(N)*(1+percAmpl), N)
# (linkScal(N), recScal(N))

# η_0_minMax = sort([η_0_min, η_0_max])
# θ_0_minMax = sort([θ_0_min, θ_0_max])
# resSin = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0_minMax, η_0_minMax, "sin"), N, Nsample)
# resAR = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0_minMax, η_0_minMax, "AR"), N, Nsample)

