
# script that wants to numerically test chatteris diaconis (misspelled with 99% prob)
# for the estimates of beta, fitness,ergm (many names..) in the DirBin1 case
using Utilities,AReg,StaticNets,JLD,MLBase,StatsBase#,DynNets using PyCall; pygui(:qt); using PyPlot

Nsample = 100


# load my data
halfPeriod = false
## Load dataj
fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
loadFilePartialName = "Weekly_eMid_Data_from_"
halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
@load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)

matY_T = YeMidWeekly_T


function logLikelihoodGammaVectorObs(μ_vec::Array{<:Real,1},Y_vec::Array{<:Real,1};α::T where T<:Real = 1.0)
    #check that the vector contains only non zero elements
    #prod(Y_vec .>0)?():error()
    ratioVec = Y_vec./μ_vec
    if α!=1
        error("Need to add the term counting the likelihood of zero observations")
        loglike =  α * sum(log.(ratioVec)) - sum(ratioVec) - sum(log.(Y_vec)) - log(gamma(α))
    else
        loglike = - sum(log.(μ_vec) ) - sum(ratioVec)
    end
    return loglike
end


function logLikeGammaDirMatNodePar(parR::Array{<:Real,1},parC::Array{<:Real,1},Y_mat::Array{<:Real,2};α::T where T<:Real = 1.0)

    # create vectors of nonzero observations needed to compute the likelihood of
    # gamma distributions for a matrix with parameters μ_i + μ_j
    A_mat = Y_mat.>0
    indVec =  find(A_mat)
    indR,indC = ind2sub(A_mat,indVec)
    Nnnz = sum(A_mat)
    Y_vec = zeros(Real,Nnnz)
    μ_vec= zeros(Real,Nnnz)
    for n=1:Nnnz
         i,j = indR[n],indC[n]
        #println((i,j))
        Y_vec[n] = Y_mat[i,j]
        μ_vec[n] = parR[i] + parC[j]
    end
    return logLikelihoodGammaVectorObs(μ_vec,Y_vec;α = α )
end



Y_mat = YeMidWeekly_T[:,:,20]
A_mat = Y_mat.>0
N = length(Y_mat[1,:])
logLikeGammaDirMatNodePar(ones(N),ones(N),Y_mat)

using Optim
tol = eps()*10
opt = Optim.Options(  g_tol = tol,
                 x_tol = tol,
                 f_tol = tol,
                 iterations = 5000,
                 show_trace = true,#false,#
                 show_every=1)

algo = Newton(; alphaguess = LineSearches.InitialHagerZhang(),
                 linesearch = LineSearches.BackTracking())

algo = NewtonTrustRegion(; initial_delta = 0.5,
                delta_hat = 1.0,
                eta = 0.1,
                rho_lower = 0.25,
                rho_upper = 0.75)# Newton(;alphaguess = Optim.LineSearches.InitialStatic(),


function objfun(vUnPar::Array{<:Real,1})
    reParR = exp.(vUnPar[1:N])
    reParC = exp.(vUnPar[1+N:end])
    return  - logLikeGammaDirMatNodePar(reParR,reParC,Y_mat)
end

parRC_0 = ones(2N)
ADobjfun = TwiceDifferentiable(objfun, parRC_0; autodiff = :forward);

objfun(parRC_0)

optim_out = optimize(ADobjfun,parRC_0,algo,opt)

Optim.minimizer(optim_out)


#
