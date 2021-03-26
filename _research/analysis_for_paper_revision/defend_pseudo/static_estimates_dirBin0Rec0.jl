
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""



import StaticNets: NetModelDirBin0Rec0(), ergm_par_from_mean_vals,diadProbFromPars , samplSingMatCan, statsFromMat, get_mple, estimate,NetModelDirBin0Rec0, exp_val_stats
using ErgmRcall
using PyPlot
using RCall
using Statistics

using Utilities

ErgmRcall.clean_start_RCall()
ergmTermsString = "edges +  mutual"
R"""options(warn=-1) """
model = NetModelDirBin0Rec0()
    
##-------------------- Test and COmpare MLE and MPLE estimates
#Compare for a single value of the parameters
N=100
Ldgp = n_pox_dir_links(N)/4
Rdgp = n_pox_pairs(N)/4-100

θ_0, η_0 = Tuple(estimate(model, Ldgp, Rdgp, N))

nSample = 100
diadProb = diadProbFromPars(model, [θ_0, η_0])
A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]
statsVec = reduce(hcat,[statsFromMat(model, A) for A in A_vec])
pygui(true)

fig, ax = subplots(2,1)
ax[1].hist(statsVec[1,:], 30)
ax[1].vlines(Ldgp, 0, 1, "k")
ax[2].hist(statsVec[2,:], 30)
ax[2].vlines(Rdgp, 0, 1, "k")



parMpleSingle = reduce(hcat, [get_mple(A,ergmTermsString) for A in A_vec])
parMleSingle = reduce(hcat,[estimate(model, A) for A in A_vec])
fix, ax = subplots(2,1)
ax[1].hist(parMleSingle[2,isfinite.(parMleSingle[2,:])], density=true, color="b" )
ax[1].hist(parMleSingle[1,isfinite.(parMleSingle[1,:])], density=true, color="b" )
ax[1].vlines(θ_0, 0, 1, "k")
ax[1].vlines(η_0, 0, 1, "k")
ax[1].set_title("mle")
ax[2].hist(parMpleSingle[2,isfinite.(parMpleSingle[2,:])], density=true, color="r" )
ax[2].hist(parMpleSingle[1,isfinite.(parMpleSingle[1,:])], density=true, color="r" )
ax[2].vlines(θ_0, 0, 1, "k")
ax[2].vlines(η_0, 0, 1, "k")
ax[2].set_title("mple")

function mle_pmle_comparison_var_Net_size(model::NetModelDirBin0Rec0, netSizes, linkFun, recFun; nSample=30)
    # compare the two estimators over different  par values for the dgp
    nPars = length(netSizes)
    parMmle = zeros(nPars, 2, nSample)
    parMple = zeros(nPars, 2, nSample)
    parDgp = zeros(nPars, 2)
    model = NetModelDirBin0Rec0()
    for i=1:nPars
        N = netSizes[i]
        θ_0, η_0 = ergm_par_from_mean_vals(model, linkFun(N), recFun(N), N)
        parDgp[i,:] = [θ_0, η_0]
        diadProb = diadProbFromPars(model, [θ_0, η_0])
        A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]
        parMmle[i,:,:] = reduce(hcat,[estimate(model, A) for A in A_vec])
        parMple[i,:,:] = reduce(hcat, [get_mple(A,ergmTermsString) for A in A_vec])
    end
    return parDgp, parMmle, parMple
end

# Compare for multiple values of the parameters
netSizes = collect(15:5:40)

parDgp, parMle, parMple = mle_pmle_comparison_var_Net_size(model, netSizes, x->x, x->x/5 )

[parMle[i, :, :] ]

i=1
bias_mle = reduce(hcat, [mean(parDgp[i] .- parMle[i,:,:], dims=2 ) for i=1:length(parMle[:,1,1])] )
bias_mple = reduce(hcat, [mean(parDgp[i] .- parMple[i,:,:], dims=2 ) for i=1:length(parMle[:,1,1])] )
rmse_mle = reduce(hcat, [sqrt.(mean((parDgp[i] .- parMle[i,:,:]).^2, dims=2 )) for i=1:length(parMle[:,1,1])] )
rmse_mple = reduce(hcat, [sqrt.(mean((parDgp[i] .- parMple[i,:,:]).^2, dims=2) ) for i=1:length(parMle[:,1,1])] )

xplot = parDgp[2,:]
fix, ax = subplots(2,2)
ax[1].plot(xplot, bias_mle', "--")
ax[1].set_title("bias mle")
ax[2].plot(xplot, rmse_mle',"-")
ax[2].set_title("rmse mle")
ax[3].plot(xplot, bias_mple', "--")
ax[3].set_title("bias mple")
ax[4].plot(xplot, rmse_mple',"-")
ax[4].set_title("rmse mple")


plot(parMle[2,:], parMple[2,:], ".")
plot(parMle[1,:], parMple[1,:], ".")



##----------------------------- Explore logLikelihoods around the estimates

function pseudo_loglikelihood_ddg(Model::NetModelDirBin0Rec0, A, par)
    θ, η = par
    L, R, N = statsFromMat(Model, A)

    return L * θ + R*η - sum(log.(1 .+ exp.(2θ .+ (η).*A )) ) 
end



θ_0 = 1
η_0 = -0.5
par_dgp = [θ_0, η_0]
N = 200
fig, ax = subplots(2,1)
for i=1:3
    A = samplSingMatCan(model, diadProbFromPars(model, [θ_0, η_0]), N) 
    mle = estimate(model, A)
    mple = get_mple(A,ergmTermsString)
    println(mle.-mple)

    changeStat, response, weights = decomposeMPLEmatrix(get_change_stats(A,ergmTermsString))



    ## Are the two functions similar around the estimates?
    x_vals = LinRange(-0.5,0.5,500) 
    par_vec_1 =  [mle .+ [0, x] for x in x_vals]
    par_vec_2 =  [mle .+ [x, 0] for x in x_vals]


    cost_diff_si = (logLikelihood(model, A, mle) - pseudo_loglikelihood_strauss_ikeda(mle, changeStat, response, weights)) ./N^2
    cost_diff_ddg = (logLikelihood(model, A, mle) - pseudo_loglikelihood_ddg(model, A, mle))./N^2 

    suptitle("exact vs pseudo log'likelihoods, /constant")
    ax[1].plot(x_vals,[logLikelihood(model, A, par) for par in par_vec_1]./N^2, "-b")
    ax[1].plot(x_vals, -cost_diff_si .+[pseudo_loglikelihood_strauss_ikeda(par, changeStat, response, weights) for par in par_vec_1]./N^2, "--r")
    #ax[1].legend(["log likelihood", "pesudo log likelihood"])
    ax[1].grid()
    ax[2].plot(x_vals,[logLikelihood(model, A, par) for par in par_vec_2]./N^2, "-b")
    ax[2].plot(x_vals,cost_diff_si .+[pseudo_loglikelihood_strauss_ikeda(par, changeStat, response, weights) for par in par_vec_2]./N^2, "--r")
    ax[2].grid()
end

function grad(Model::NetModelDirBin0Rec0, A, par)
    θ, η = par
    L, R, N = statsFromMat(Model, A) 
    z = 1 + 2*exp(θ) + exp(2*θ+η)
    g_θ =  L - (N*(N-1)/2) *(2*exp(θ)+ 2* exp(2*θ+η))/z
    g_η =  R - (N*(N-1)/2) *( exp(2*θ+η))/z
    return [g_θ, g_η]
end

logl(x) = logLikelihood(model, A, x)
pseudologl(x) = pseudo_loglikelihood_strauss_ikeda(x, changeStat, response, weights)

using ForwardDiff
g_logl(x) = ForwardDiff.gradient(logl,x)
g_pseudologl(x) = ForwardDiff.gradient(pseudologl,x)
h_logl(x) = ForwardDiff.hessian(logl,x)
h_pseudologl(x) = ForwardDiff.hessian(pseudologl,x)
# compare gradients
fig, ax = subplots(2,1)
suptitle("exact vs pseudo gradients along theta")
ax[1].plot(x_vals,[g_logl(par) for par in par_vec_1])
ax[1].plot(x_vals,[g_pseudologl(par) for par in par_vec_1])
ax[1].legend(["log likelihood", "pesudo log likelihood"])
ax[1].grid()
ax[2].plot(x_vals,[g_logl(par) for par in par_vec_2])
ax[2].plot(x_vals,[g_pseudologl(par) for par in par_vec_2])
ax[2].grid()
# compare hessians
fig, ax = subplots(2,2)
suptitle("exact vs pseudo hessians along theta move")
h_logl_vec = [h_logl(par) for par in par_vec_1]
h_pseudologl_vec = [h_pseudologl(par) for par in par_vec_1]
ax[1].plot(x_vals,[h[1,1] for h in h_logl_vec])
ax[1].plot(x_vals,[h[1,1] for h in h_pseudologl_vec])
ax[2].plot(x_vals,[h[2,1] for h in h_logl_vec])
ax[2].plot(x_vals,[h[2,1] for h in h_pseudologl_vec])
ax[4].plot(x_vals,[h[2,2] for h in h_logl_vec])
ax[4].plot(x_vals,[h[2,2] for h in h_pseudologl_vec])
ax[1].legend(["log likelihood", "pesudo log likelihood"])
ax[1].grid()
ax[2].grid()
ax[3].grid()
ax[4].grid()

fig, ax = subplots(2,2)
suptitle("exact vs pseudo hessians along eta move")
h_logl_vec = [h_logl(par) for par in par_vec_2]
h_pseudologl_vec = [h_pseudologl(par) for par in par_vec_2]
ax[1].plot(x_vals,[h[1,1] for h in h_logl_vec])
ax[1].plot(x_vals,[h[1,1] for h in h_pseudologl_vec])
ax[2].plot(x_vals,[h[2,1] for h in h_logl_vec])
ax[2].plot(x_vals,[h[2,1] for h in h_pseudologl_vec])
ax[4].plot(x_vals,[h[2,2] for h in h_logl_vec])
ax[4].plot(x_vals,[h[2,2] for h in h_pseudologl_vec])
ax[1].legend(["log likelihood", "pesudo log likelihood"])
ax[1].grid()
ax[2].grid()
ax[3].grid()
ax[4].grid()



# Closed form gradient vs AD gradient
fig, ax = subplots(2,1)
ax[1].plot(x_vals,[g_logl(par)[1] for par in par_vec_1])
ax[1].plot(x_vals,[grad(model, A,par)[1] for par in par_vec_1])
ax[2].plot(x_vals,[g_logl(par)[2] for par in par_vec_1])
ax[2].plot(x_vals,[grad(model, A,par)[2] for par in par_vec_1])










#
