""
# functions that allow to sample sequences of ergms with different parameters' values from R package ergm
# and estimate the ergm

module ErgmRcall

using DataFrames
using RCall

using RCall



function clean_start_RCall()
    # Clean R enviroment and load needed packages
    R"rm(list = ls())
    library(ergm)
    library(sna)
    library(coda)
    library(network)
    sessionInfo()"
   R"set.seed(0)"
end
export clean_start_RCall


function __init__()
    R"""options(warn=-1) """

    R"""local({r <- getOption("repos")
        r["CRAN"] <- "http://cran.r-project.org" 
        options(repos=r)})"""

    clean_start_RCall()
end


function sampleErgmRcall(parDgpT,N,Nsample,formula_ergm_str)
    """
    Function that samples from an ergm defined by formula_ergm_str (according to
    the notation of R package ergm)
    in order to work it needs RCall.jl insatalled and the packages listed at the
    beginning of the module installed in R
    """

    clean_start_RCall()
    T = size(parDgpT)[2]
    # For each t, sample the ERGM with parameters corresponding to the DGP at time t
    # and run a single snapeshot estimate
    @rput T; @rput parDgpT;@rput N;@rput Nsample

    reval("formula_ergm = net ~ "*formula_ergm_str)

     #create an empty network, the formula defining ergm, sample the ensemble and
     # store the sufficient statistics and change statistics in R
     R"
     net <- network.initialize(N)
       sampledMat_T_R =    array(0, dim=c(N,N,T,Nsample))
       changeStats_T_R = list()
       stats_T_R = list()
        for(t in 1:T){
            changeStats_t_R = list()
            stats_t_R = list()
            print(t)
            for(n in 1:Nsample){
                 net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parDgpT[,t],control = control.simulate.formula(MCMC.burnin = 100000))
                 sampledMat_T_R[,,t,n] <- as.matrix.network( net)
                 print(c(t,n,parDgpT[,t]))
                 chStat_t <- ergmMPLE(formula_ergm)
                 changeStats_t_R[[n]] <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
                 stats_t_R[[n]] <- summary(formula_ergm)
                 }
                  changeStats_T_R[[t]] <-changeStats_t_R
                  stats_T_R[[t]] <- stats_t_R
             }"



    # import sampled networks in julia
    sampledMat_T = BitArray( @rget(sampledMat_T_R))
    changeStats_T = @rget changeStats_T_R;# tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  changeStats_T[t];end;changeStats_T = tmp
    stats_T = @rget(stats_T_R); #tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp
    return sampledMat_T, changeStats_T, stats_T
end
export sampleErgmRcall

function estimateErgmRcall(sampledMat_T_R , formula_ergm_str)
    """
    Function that estimates a sequence of ergm defined by formula_ergm_str (according to
    the notation of R package ergm)
    in order to work it needs RCall.jl insatalled and the packages listed at the
    beginning of the module installed in R
    """


    clean_start_RCall()

    T = size(sampledMat_T_R )[3]
    # For each t, sample the ERGM with parameters corresponding to the DGP at time t
    # and run a single snapeshot estimate
    @rput T; @rput sampledMat_T_R ;@rput N

    reval("formula_ergm = net ~ "*formula_ergm_str)

     #create an empty network, the formula defining ergm, sample the ensemble and
     # store the sufficient statistics and change statistics in R

     R"
      estParSS_T_R =    list()
        for(t in 1:T){
            net <- as.network.matrix(   sampledMat_T_R[,,t,n])
            estParSS_t_R =    list()
            changeStats_t_R = list()
            stats_t_R = list()
            print(t)
            for(n in 1:Nsample){
                print(parDgpT[,t])
                 net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parDgpT[,t],control = control.simulate.formula(MCMC.burnin = 100000))
                 tmp <- ergm(formula_ergm)#,estimate = 'MPLE')#)#
                 estParSS_t_R[[n]] <- tmp[[1]]
                 print(c(t,n))

                 print(estParSS_t_R[[n]])
                 }
                  estParSS_T_R[[t]] <- estParSS_t_R
             }"



    # import sampled networks in julia
    estParSS_T = @rget(estParSS_T_R);#tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = estParSS_T[t]; end ; estParSS_T = tmp
    estParSS_T = Utilities.collapseArr3(estParSS_T)
    return estParSS_T
end
export estimateErgmRcall




function get_change_stats(A::Matrix{T} where T<:Integer, ergmTermsString::String)
    # given a matrix returns the change statistics wrt to a given formula
    @rput A
    reval("formula_ergm = net ~ "* ergmTermsString)
    # store the sufficient statistics and change statistics in R
    R"""
        net <- network(A)
        chStat_t <- ergmMPLE(formula_ergm)
        changeStats_t_R <- cbind(chStat_t$response, chStat_t$predictor,     chStat_t$weights)
            """

    changeStats = @rget changeStats_t_R;# tmp = Array{Array{Float64,2}}(T); for 
    return changeStats
end
export get_change_stats

function get_mple(A::Matrix{T} where T<:Integer, ergmTermsString::String)
    @rput A
    reval("formula_ergm = net ~ "* ergmTermsString)
    # store the sufficient statistics and change statistics in R
    R"""
        net <- network(A)
        mple_t_R <- ergmMPLE(formula_ergm, output="fit")$coef
        print(mple_t_R)
        """

    mple_t = @rget mple_t_R;# tmp = Array{Array{Float64,2}}(T); for 
    return mple_t
end
export get_mple

function get_static_mple(A_T::Array{T,3} where T<:Integer, ergmTermsString::String)
    clean_start_RCall()
    T = size(sampledMat_T_R )[3]
    # For each t, sample the ERGM with parameters corresponding to the DGP at time t
    # and run a single snapeshot estimate
    @rput T
    @rput A_T
    reval("formula_ergm = net ~ "* ergmTermsString)
    # store the sufficient statistics and change statistics in R
    R"""
    net <- network(A_T[,,1])
    mple_t_R <- ergmMPLE(formula_ergm)   
    dfPred = data.frame( mple_t_R$predictor )
    dfWeights = data.frame( mple_t_R$weights )
    dfResp = data.frame( mple_t_R$response )
    for (t in 2:T){
        net <- network(A_T[,,t])
        mple_t_R <- ergmMPLE(formula_ergm)
        dfPred <- rbind(dfPred, data.frame(mple_t_R$predictor))
        dfWeights <- rbind(dfWeights, data.frame(mple_t_R$weights))
        dfResp <- rbind(dfResp, data.frame(mple_t_R$response))
        df <- rbind(df, data.frame(mple_t_R))
        }
        mple_static_R <- glm(dfResp$mple_t_R.response ~ . -1, data = dfPred, weights = dfWeights$mple_t_R.weights,family="binomial")$coefficients
        """

    mple_static = @rget mple_static_R;# tmp = Array{Array{Float64,2}}(T); for 

    return mple_static
end

function get_static_mple(changeStats_T::Array{T,2} where T<:Real, ergmTermsString::String)
    clean_start_RCall()
    # For each t, sample the ERGM with parameters corresponding to the DGP at time t
    # and run a single snapeshot estimate
    @rput changeStats_T
    reval("formula_ergm = net ~ "* ergmTermsString)
    # store the sufficient statistics and change statistics in R
    R"""
        mple_static_R <- glm(changeStats_T[,1] ~ . -1, data = data.frame(changeStats_T[,c(2:(ncol(changeStats_T)-1))]), weights =changeStats_T[,ncol(changeStats_T)],family="binomial")$coefficients
        """

    mple_static = @rget mple_static_R;# tmp = Array{Array{Float64,2}}(T); for 

    return mple_static
end
function get_static_mple(changeStats_T::Array{Array{T,2},1} where T<:Real, ergmTermsString::String)
    return get_static_mple(reduce(vcat,changeStats_T), ergmTermsString)
end
export get_static_mple



function decomposeMPLEmatrix(ergmMPLEmatrix)  
    changeStat = ergmMPLEmatrix[:,2:3]
    response = ergmMPLEmatrix[:,1]
    weights = ergmMPLEmatrix[:,4]
    return changeStat, response, weights 
end
export decomposeMPLEmatrix

end
