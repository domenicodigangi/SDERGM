
# functions that allow to sample sequences of ergms with different parameters' values from R package ergm
# and estimate the ergm

module ergmRcall

using Distributions, StatsBase,Optim, LineSearches, StatsFuns, Roots,MLBase
using LinearAlgebra
using JLD,DataFrames
using StaNets,HelperFunDom
using ForwardDiff,NLSolversBase
using PyCall;# pygui(:qt);
using PyPlot
using RCall
using HelperFunDom
using DynNets

R"rm(list = ls())
  library(statnet)
   library(ergm)
   library(sna)
   library(coda)
   library(network)
   sessionInfo()"

R"set.seed(0)"

function sampleErgmRcall(parMatDgp_T,N,Nsample,formula_ergm_str)
    """
    Function that samples from an ergm defined by formula_ergm_str (according to
    the notation of R package ergm)
    in order to work it needs RCall.jl insatalled and the packages listed at the
    beginning of the module installed in R
    """

    # Clean R enviroment and load needed packages
    R"rm(list = ls())
      library(statnet)
       library(ergm)
       library(sna)
       library(coda)
       library(network)
       sessionInfo()"

    R"set.seed(0)"

    T = size(parMatDgp_T)[2]
    # For each t, sample the ERGM with parameters corresponding to the DGP at time t
    # and run a single snapeshot estimate
    @rput T; @rput parMatDgp_T;@rput N;@rput Nsample

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
                 net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parMatDgp_T[,t],control = control.simulate.formula(MCMC.burnin = 100000))
                 sampledMat_T_R[,,t,n] <- as.matrix.network( net)
                 print(c(t,n,parMatDgp_T[,t]))
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

function estimateErgmRcall(sampledMat_T_R , formula_ergm_str)
    """
    Function that estimates a sequence of ergm defined by formula_ergm_str (according to
    the notation of R package ergm)
    in order to work it needs RCall.jl insatalled and the packages listed at the
    beginning of the module installed in R
    """


        # Clean R enviroment and load needed packages
    R"rm(list = ls())
      library(statnet)
       library(ergm)
       library(sna)
       library(coda)
       library(network)
       sessionInfo()"

    R"set.seed(0)"

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
                print(parMatDgp_T[,t])
                 net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = parMatDgp_T[,t],control = control.simulate.formula(MCMC.burnin = 100000))
                 tmp <- ergm(formula_ergm)#,estimate = 'MPLE')#)#
                 estParSS_t_R[[n]] <- tmp[[1]]
                 print(c(t,n))

                 print(estParSS_t_R[[n]])
                 }
                  estParSS_T_R[[t]] <- estParSS_t_R
             }"



    # import sampled networks in julia
    estParSS_T = @rget(estParSS_T_R);#tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = estParSS_T[t]; end ; estParSS_T = tmp
    estParSS_T = HelperFunDom.collapseArr3(estParSS_T)
    return estParSS_T
end




function sdergmDGPpaper_Rcall( Model::DynNets.GasNetModelDirBinGlobalPseudo, indTvPar::BitArray{1},T::Int,N::Int;vConstPar ::Array{<:Real,1} = zeros(Real,2),
                                ftot_0::Array{<:Real,1} = zeros(Real,2),UM = [-3 , 0.05],B_Re  = 0.9, A_Re  = 0.01 )
    """ SDERGM DGP  for density and GWESP as described in the paper
    Heavily relying on R to smaple from the ERGM at each t, and it assumes that
    all needed packages have been loaded in R via Rcall
     """
    indTvPar = trues(2)
    ftot_0 = zeros(Real,2)
    indTargPar = falses(2)
    NTargPar = 0


    NergmPar = 2#
    NTvPar   = sum(indTvPar)

    function initialize_pars(vParOptim_0)
        last = 0
        for i=1:NergmPar
            if indTvPar[i]
                if indTargPar[i]
                    vParOptim_0[last+1:last+2] = [ B_Re; A_Re]
                    last+=2
                else
                    vParOptim_0[last+1:last+3] = [UM[i]*(1-B_Re) ; B_Re; A_Re]
                    last+=3
                end
            else
                vParOptim_0[last+1] = UM[i]
                last+=1
            end
        end
        return vParOptim_0
    end
    vResGasPar = initialize_pars(zeros(NergmPar + NTvPar*2 -NTargPar))
    # Organize parameters of the GAS update equation
    Wvec = vResGasPar[1:3:3*NTvPar]
    Bvec = vResGasPar[2:3:3*NTvPar]
    Avec = vResGasPar[3:3:3*NTvPar]


    # start values equal the unconditional mean,and  constant ones remain equal to the unconditional mean, hence initialize as:
    UMallPar = zeros(Real,NergmPar)
    UMallPar[indTvPar] =  Wvec ./ (1 .- Bvec)
    if !prod(indTvPar) # if not all parameters are time varying
        UMallPar[.!indTvPar] = vConstPar
    end
    #println(UMallPar)
    fVecT = ones(Real,NergmPar,T)
    obsT = fill(zeros(2,2),T)

    statsT = ones(Real,NergmPar,T)

    sum(ftot_0)==0  ?    ftot_0 = UMallPar : ()# identify(Model,UMallNodesIO)

    ftot_t = copy(ftot_0)
    if NTvPar==0
        I_tm1 = ones(1,1)
    else
        I_tm1 = Float64.(Diagonal{Real}(I,NTvPar))
    end
    formula_ergm_str = " edges +  gwesp(decay = 0.25,fixed = TRUE)"
    ftot_t
    for t=1:T

        ftot_t_R = ftot_t# cat(,ftot_t,dims = 2)
        # Sample the observations using RCall -----------------------
        @rput ftot_t_R ; @rput N;
        reval("formula_ergm = net ~ "*formula_ergm_str)
         #create an empty network, the formula defining ergm, sample the ensemble and
         # store the sufficient statistics and change statistics in R
         R"
         net <- network.initialize(N)
         net <- simulate(formula_ergm, nsim = 1, seed = sample(1:100000000,1), coef = as.numeric(ftot_t_R), control = control.simulate.formula(MCMC.burnin = 100000))
         sampledMat_t_R <- as.matrix.network( net)
         chStat_t <- ergmMPLE(formula_ergm)
         changeStats_t_R <- cbind(chStat_t$response, chStat_t$predictor,chStat_t$weights)
         stats_t_R <- summary(formula_ergm)
        "
        # import sampled networks in julia
        sampledMat_t = BitArray( @rget(sampledMat_t_R))
        changeStats_t = @rget changeStats_t_R;# tmp = Array{Array{Float64,2}}(T); for t=1:T tmp[t] =  changeStats_T[t];end;changeStats_T = tmp
        #print(size(changeStats_t))
        stats_t = @rget(stats_t_R); #tmp = zeros(Nterms,T); for t=1:T tmp[:,t] = stats_T[t];end;stats_T = tmp

        #-----------------------------------------------------------------
        obs_t = changeStats_t # vector of in and out degrees
        obsT[t] = obs_t
        statsT[:,t] = stats_t
        #print((t,I_tm1))
        ftot_t,loglike_t,I_tm1,~ = DynNets.updatedGasPar(Model,obs_t,ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
        fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
    end


    return obsT,statsT,fVecT
    end




end
