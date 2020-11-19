

"""
Script for comparison of two versions of  GasNetModelDirBin0Rec0, i.e. one parameter for total number of links and one for reciprocity, a version based on MLE and one on PMLE

To Do: 

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
Number of links for various network size, in semiDense regime
"""
semiDenseLScal(N)  = 0.1 * (N^2-N)/sqrt(N)
"""
Number of links for various network size, in sparse regime
"""
sparseLScal(N) =  1 * N 

"""
Average number of reciprocal pairs in Erdos Reny Model 
 - avgL is the expected number of links, related to the probability p of a single link being present by expL = (N^2-N) * p  

"""
erdosRenyRecScal(expL, N) = expL^2/(2*(N^2-N)) 




## Repeat comparison for varying network size

function scaling_comparison(model_mle::GasNetModelDirBin0Rec0_mle, Nvals, linkScal, Nsample, T; plotFlag = false, dgpShape = "sin", regimeString="", percAmpl = 0.2 )
    
    staModel = fooNetModelDirBin0Rec0

    recScal(N) = 5 + 10 * erdosRenyRecScal(linkScal(N), N)

    rmse_mle = []
    rmse_pmle = []
    mean_vals = []
    
    for N in Nvals 

        #  fluctuations in the number of links
        minL = linkScal(N)*(1-percAmpl)
        maxL = linkScal(N)*(1+percAmpl)

        # fluctuations in the number of rec pairs within the physical bounds
        minR = recScal(N)*(1-percAmpl)
        minR = max(maxL - (N^2-N)/2, minR)
        minR = min(minL/2-1, minR)

        maxR = recScal(N)*(1+percAmpl)
        maxR = min(minL/2-1, maxR)
        maxR = max(maxL - (N^2-N)/2, maxR)


        minR>maxR ? ( x=deepcopy(minR); minR=deepcopy(maxR); maxR=x) : () 

        titleString = " minL=$(round(minL, digits=2)), maxL=$(round(maxL, digits=2)), minR=$(round(minR, digits=2)), maxR=$(round(maxR, digits=2)), regime=$regimeString)"


        θ_0_min, η_0_min = ergm_par_from_mean_vals(staModel, minL, minR, N)
        θ_0_max, η_0_max = ergm_par_from_mean_vals(staModel, maxL, maxR, N)
        
        push!(mean_vals, (linkScal(N), recScal(N)) )

        η_0_minMax = sort([η_0_min, η_0_max])
        θ_0_minMax = sort([θ_0_min, θ_0_max])
        res = sample_est_mle_pmle(model_mle, dgp_missp(model_mle, T, θ_0_minMax, η_0_minMax, dgpShape), N, Nsample; plotFlag=plotFlag, regimeString=titleString)


        push!(rmse_pmle, res.rmse_pmle)
        push!(rmse_mle, res.rmse_mle)
    end

    return rmse_mle, rmse_pmle, mean_vals
end



using Serialization

T = 200
# Nvals= [20, 100]# round.(Int, collect(20:20:100))
Nvals=  round.(Int, collect(20:20:100))
Nsample = 50

savePath = (@__DIR__ )* "\\revision_data"
percAmpl = 0.1

rmse_mle_sparse, rmse_pmle_sparse, mean_vals_sparse = scaling_comparison(model_mle, Nvals, sparseLScal, Nsample, T; plotFlag = true, regimeString = "sparse", percAmpl = percAmpl)
serialize(savePath * "\\sparse_Nsample_$(Nsample)_$(percAmpl).jls", (rmse_mle_sparse, rmse_pmle_sparse, mean_vals_sparse))

rmse_mle_semDense, rmse_pmle_semDense, mean_vals_semDense = scaling_comparison(model_mle, Nvals, semiDenseLScal, Nsample, T; plotFlag = true, regimeString = "semi - dense", percAmpl = percAmpl)
serialize(savePath * "\\semiDense_Nsample_$(Nsample)_$(percAmpl).jls", (rmse_mle_semDense, rmse_pmle_semDense, mean_vals_semDense))

rmse_mle_dense, rmse_pmle_dense, mean_vals_dense = scaling_comparison(model_mle, Nvals, denseLScal, Nsample, T; plotFlag = true, regimeString = "dense", percAmpl = percAmpl)

serialize(savePath * "\\dense_Nsample_$(Nsample)_$(percAmpl).jls", (rmse_mle_dense, rmse_pmle_dense, mean_vals_dense))


close()



Nsample = 50
percAmpl = 0.5 # 0.3, 0.5
rmse_mle_sparse, rmse_pmle_sparse, mean_vals_sparse = deserialize(savePath * "\\sparse_Nsample_$(Nsample)_$(percAmpl).jls")
rmse_mle_semDense, rmse_pmle_semDense, mean_vals_semDense = deserialize(savePath * "\\semiDense_Nsample_$(Nsample)_$(percAmpl).jls")
rmse_mle_dense, rmse_pmle_dense, mean_vals_dense = deserialize(savePath * "\\dense_Nsample_$(Nsample)_$(percAmpl).jls")


fig, ax = subplots(2,3)
ax[1,1].plot(Nvals, reduce(hcat,rmse_mle_sparse)[1,:], "-.b")
ax[1,1].plot(Nvals, reduce(hcat,rmse_pmle_sparse)[1,:], "-.r")
ax[2,1].plot(Nvals, reduce(hcat,rmse_mle_sparse)[2,:], "-.b")
ax[2,1].plot(Nvals, reduce(hcat,rmse_pmle_sparse)[2,:], "-.r")
ax[1,1].set_title("Sparse regime percAmpl = $percAmpl")

ax[1,2].plot(Nvals, reduce(hcat,rmse_mle_semDense)[1,:], "-.b")
ax[1,2].plot(Nvals, reduce(hcat,rmse_pmle_semDense)[1,:], "-.r")
ax[2,2].plot(Nvals, reduce(hcat,rmse_mle_semDense)[2,:], "-.b")
ax[2,2].plot(Nvals, reduce(hcat,rmse_pmle_semDense)[2,:], "-.r")
ax[1,2].set_title("Semi-dense  regime")

ax[1,3].plot(Nvals, reduce(hcat,rmse_mle_dense)[1,:], "-.b")
ax[1,3].plot(Nvals, reduce(hcat,rmse_pmle_dense)[1,:], "-.r")
ax[2,3].plot(Nvals, reduce(hcat,rmse_mle_dense)[2,:], "-.b")
ax[2,3].plot(Nvals, reduce(hcat,rmse_pmle_dense)[2,:], "-.r")
ax[1,3].set_title("Dense  regime")



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

