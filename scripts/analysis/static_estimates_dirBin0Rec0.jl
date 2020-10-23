
"""
Test script for dirBin0Rec0 model: one parameter for total number of links and one for reciprocity
"""

include("../../add_load_paths.jl")

using StaticNets
using ergmRcall
using PyPlot
using RCall


ergmRcall.clean_start_RCall()
ergmTermsString = "edges +  mutual"
R"""options(warn=-1) """
model = fooNetModelDirBin0Rec0
    
##-------------------- Test and COmpare MLE and MPLE estimates
#Compare for a single value of the parameters
nSample = 50
N=100
θ_0 = 3
η_0 = -5
par_dgp = [θ_0, η_0] 
diadProb = diadProbFromPars(model, [θ_0, η_0])
A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]
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

function mle_pmle_comparison(model::NetModelDirBin0Rec0, parDgpVec; nSample=100)
    # compare the two estimators over different  par values for the dgp
    nPars = length(parDgpVec)
    par_mle = zeros(nPars, 2, nSample)
    par_mple = zeros(nPars, 2, nSample)
    model = fooNetModelDirBin0Rec0
    for i=1:nPars
        par_dgp = parDgpVec[i] 
        diadProb = diadProbFromPars(model, [θ_0, η_0])
        A_vec = [samplSingMatCan(model, diadProb, N) for i=1:nSample]
        par_mle[i,:,:] = reduce(hcat,[estimate(model, A) for A in A_vec])
        par_mple[i,:,:] = reduce(hcat, [get_mple(A,ergmTermsString) for A in A_vec])
    end
    return par_mle, par_mple
end

# Compare for multiple values of the parameters
θ_0 = 3
η_0 = -5
parDgpVec = [[θ_0, η_0 + delta] for delta in -2:0.2:2 ]
par_mle, par_mple = mle_pmle_comparison(model, parDgpVec)

i=1
bias_mle = reduce(hcat, [mean(parDgpVec[i] .- par_mle[i,:,:], dims=2 ) for i=1:length(par_mle[:,1,1])] )
bias_mple = reduce(hcat, [mean(parDgpVec[i] .- par_mple[i,:,:], dims=2 ) for i=1:length(par_mle[:,1,1])] )
rmse_mle = reduce(hcat, [sqrt.(mean((parDgpVec[i] .- par_mle[i,:,:]).^2, dims=2 )) for i=1:length(par_mle[:,1,1])] )
rmse_mple = reduce(hcat, [sqrt.(mean((parDgpVec[i] .- par_mple[i,:,:]).^2, dims=2) ) for i=1:length(par_mle[:,1,1])] )

xplot = reduce(hcat,parDgpVec)[2,:]
fix, ax = subplots(2,2)
ax[1].plot(xplot, bias_mle', "--")
ax[1].set_title("bias mle")
ax[2].plot(xplot, rmse_mle',"-")
ax[2].set_title("rmse mle")
ax[3].plot(xplot, bias_mple', "--")
ax[3].set_title("bias mple")
ax[4].plot(xplot, rmse_mple',"-")
ax[4].set_title("rmse mple")


plot(par_mle[2,:], par_mple[2,:], ".")
plot(par_mle[1,:], par_mple[1,:], ".")



##----------------------------- Explore logLikelihoods around the estimates
θ_0 = 3
η_0 = -5
par_dgp = [θ_0, η_0]
N = 100
A = samplSingMatCan(model, diadProbFromPars(model, [θ_0, η_0]), N) 
mle = estimate(model, A)
mple = get_mple(A,ergmTermsString)
println(mle.-mple)

changeStat, response, weights = decomposeMPLEmatrix(get_change_stats(A,ergmTermsString))


# pseudolikelihood function from change statistics following strauss and ikeda
logit(x) = log(x/(1-x))
inv_logit(x) = 1/(1+exp(-x))


function pseudo_loglikelihood_strauss_ikeda(par, changeStat, response, weights)
    logit_P = sum(par.*changeStat', dims=1)      
    P = inv_logit.(logit_P)    
    logPVec = log.([response[i] == zero(response[i]) ? 1 - P[i] : P[i] for i=1:length(response) ])
    logPTot = sum(logPVec.*weights)
    return logPTot
end
function pseudo_loglikelihood_ddg(Model::NetModelDirBin0Rec0, A, par)
    θ, η = par
    N=size(A)[1]
    L, R = statsFromMat(Model, A)

    return L * θ + R*η - sum(log.(1 .+ exp.(2θ .+ (η).*A )) ) 
end
function logLikelihood(Model::NetModelDirBin0Rec0, A, par)
    θ, η = par
    N=size(A)[1]
    z = 1 + 2*exp(θ) + exp(2*θ+η)
    L, R = statsFromMat(Model, A)
    return L * θ + R*η - (N*(N-1)/2)*log(z)
end

## Are the two functions similar around the estimates?
x_vals = LinRange(-0.5,0.5,500) 
par_vec_1 =  [par_dgp .+ [0, x] for x in x_vals]
par_vec_2 =  [par_dgp .+ [x, 0] for x in x_vals]


cost_diff_si = logLikelihood(model, A, par_dgp) - pseudo_loglikelihood_strauss_ikeda(par_dgp, changeStat, response, weights) 
cost_diff_ddg = logLikelihood(model, A, par_dgp) - pseudo_loglikelihood_ddg(model, A, par_dgp) 

fig, ax = subplots(2,1)
suptitle("exact vs pseudo log'likelihoods, /constant")
ax[1].plot(x_vals,[logLikelihood(model, A, par) for par in par_vec_1])
# ax[1].plot(x_vals,cost_diff_si .+[pseudo_loglikelihood_strauss_ikeda(par, changeStat, response, weights) for par in par_vec_1])
#ax[1].legend(["log likelihood", "pesudo log likelihood"])
ax[1].grid()
ax[2].plot(x_vals,[logLikelihood(model, A, par) for par in par_vec_2])
# ax[2].plot(x_vals,cost_diff_si .+[pseudo_loglikelihood_strauss_ikeda(par, changeStat, response, weights) for par in par_vec_2])
ax[2].grid()



function grad(Model::NetModelDirBin0Rec0, A, par)
    θ, η = par
    L, R = statsFromMat(Model, A) 
    N=size(A)[1]
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
