struct  GasNetModelDirBin0Rec0_pmle <: GasNetModelDirBin0Rec0
     """ A gas model based on pseudolikelihood (as objective function) for
            Directed binary networks and probability depending on a generic vector
            of global statistics each associated with a time varying parameters.
         """
     obsT:: Array{Array{Float64,2},1}# for each t we have a matrix of change statistics: L and R and the network size N
     Par::Array{Array{<:Real,1},1} #  each time varying parameter has 3 static parameters in this specification.
                         # if a parameter is constant than only 1 gas parameteris present
     indTvPar :: BitArray{1} #  what parameters are time varying   ?
     scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
     # choices see function scalingMatGas
  end

fooGasPar = [ones(2), ones(2), ones(2)]
fooChangeStat = zeros(4,2)
fooGasNetModelDirBin0Rec0_pmle = GasNetModelDirBin0Rec0_pmle( [fooChangeStat, fooChangeStat], fooGasPar, trues(2), "")

export fooGasNetModelDirBin0Rec0_pmle

struct  GasNetModelDirBin0Rec0_mle <: GasNetModelDirBin0Rec0
     """ A gas model based on pseudolikelihood (as objective function) for
            Directed binary networks and probability depending on a generic vector
            of global statistics each associated with a time varying parameters.
         """
     obsT:: Array{Array{Float64,1},1}# for each t we have 2 statistics: L and R and the network size N
     Par::Array{Array{<:Real,1},1} #  each time varying parameter has 3 static parameters in this specification.
                         # if a parameter is constant than only 1 gas parameteris present
     indTvPar :: BitArray{1} #  what parameters are time varying   ?
     scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
     # choices see function scalingMatGas
  end

fooGasPar = [ones(2), ones(2), ones(2)]
fooGasNetModelDirBin0Rec0_mle = GasNetModelDirBin0Rec0_mle( [zeros(2), zeros(2), 10 .* ones(Int, 2)], fooGasPar, trues(2), "")

export fooGasNetModelDirBin0Rec0_mle

statsFromMat(Model::T where T<: GasNetModelDirBin0Rec0, A ::Matrix{<:Real}) = StaticNets.statsFromMat(StaticNets.fooNetModelDirBin0Rec0, A ::Matrix{<:Real}) 


# Relations between Static and Dynamical models: conventions on storage for
# parameters and observations
StaModType(model::T where T<: GasNetModelDirBin0Rec0 ) = StaticNets.fooNetModelDirBin0Rec0# to be substituted with a conversion mechanism

# options and conversions of parameters for optimization
setOptionsOptim(model::T where T<: GasNetModelDirBin0Rec0) = setOptionsOptim(fooGasNetModelDirBin1)

array2VecGasPar(model::T where T<: GasNetModelDirBin0Rec0,
                ArrayGasPar, indTvPar :: BitArray{1}) =
                array2VecGasPar(fooGasNetModelDirBinGlobalPseudo,
                                ArrayGasPar, indTvPar )

vec2ArrayGasPar(model::T where T<: GasNetModelDirBin0Rec0, VecGasPar::Array{<:Real,1},
                indTvPar :: BitArray{1}) =
                vec2ArrayGasPar(fooGasNetModelDirBinGlobalPseudo, VecGasPar,
                                indTvPar)

restrictGasPar(model::T where T<: GasNetModelDirBin0Rec0, vecUnGasPar::Array{<:Real,1},
                        indTvPar :: BitArray{1}) =
                        restrictGasPar(fooGasNetModelDirBinGlobalPseudo,
                                        vecUnGasPar, indTvPar )

unRestrictGasPar(model::T where T<: GasNetModelDirBin0Rec0,
                        vecReGasPar::Array{<:Real,1},
                        indTvPar :: BitArray{1})=
                        unRestrictGasPar(fooGasNetModelDirBinGlobalPseudo ,
                                            vecReGasPar,  indTvPar )

#Gas Filter Functions
function identify(model::T where T<: GasNetModelDirBin0Rec0,parIO::Array{<:Real,1})
    # "Given a vector of parameters, return the transformed vector that verifies
    # an identification condition. Do not do anything if the model is identified."
    # #set the first of the in parameters equal to one (Restricted version)
    # N,parI,parO = splitVec(parIO)
    # idType = "equalIOsums"#"firstZero"#
    # if idType == "equalIOsums"
    #     Δ = sum(parI[isfinite.(parI)]) - sum(parO[isfinite.(parO)])
    #     shift = Δ/(2N)
    # elseif idType == "firstZero"
    #     shift = parI[1]
    # end
    # parIO = [ parI - shift ;  parO + shift ]
    return parIO
end

function scalingMatGas(model::T where T<: GasNetModelDirBin0Rec0,expMat::Array{<:Real,2},I_tm1::Array{<:Real,2})
    "Return the matrix required for the scaling of the score, given the expected
     matrix and the Scaling matrix at previous time. "
    if uppercase(model.scoreScalingType) == ""
        scalingMat = 1 #
    elseif uppercase(model.scoreScalingType) == "FISHER-EWMA"
        error()
        # λ = 0.5
        #
        # I = expMat.*(1-expMat)
        # diagI = sum(I,dims = 2)
        # [I[i,i] = diagI[i] for i=1:length(diagI) ]
        # I_t =  λ*I + (1-λ) *I_tm1
        # scalingMat = sqrt(I) ##
    elseif uppercase(model.scoreScalingType) == "FISHER-DIAG"
        # display(expMat)
        #  I = expMat.*(1-expMat)
        # scalingMat = zeros(Real,2.*size(expMat))
        # diagScalIn = sqrt.(sum(I,dims = 2))
        # N = length(diagScalIn)
        # [scalingMat[i,i] = diagScalIn[i] for i=1:N ]
        # diagScalOut = sqrt.(sum(I,dims = 1))
        # # display(diagScalIn)
        # # display(diagScalOut)
        #
        # [scalingMat[N+i,N+i] = diagScalOut[i] for i=1:N ]
        error()
    end
    return scalingMat
end

function target_function_t(model::GasNetModelDirBin0Rec0_mle, obs_t, par)
    L, R, N = obs_t
    return StaticNets.logLikelihood( StaticNets.fooNetModelDirBin0Rec0, L, R, N, par)
end


function updatedGasPar( model::T where T<: GasNetModelDirBin0Rec0, obs_t,
                         ftot_t::Array{<:Real,1}, I_tm1::Array{<:Real,2},
                         indTvPar::BitArray{1}, Wgas::Array{<:Real,1},
                         Bgas::Array{<:Real,1}, Agas::Array{<:Real,1};
                         matrixScaling=false)
    #= likelihood and gradients depend on all the parameters (ftot_t), but
    only the time vaying ones (f_t) are to be updated=#
    NergmPar = 2
    NtvPar = sum(indTvPar)
    f_t = ftot_t[indTvPar] #Time varying ergm parameters

    target_fun_t(x) = target_function_t(model, obs_t, x)

    target_fun_val_t = target_fun_t(f_t)
    grad_tot_t = ForwardDiff.gradient(target_fun_t,f_t)

    # No rescaling
    I_t = I_tm1

    s_t = grad_tot_t[indTvPar]
    f_tp1 = Wgas .+ Bgas.* f_t .+ Agas.*s_t
    ftot_tp1 = copy(ftot_t)
    ftot_tp1[indTvPar] = f_tp1 #of all parameters udate dynamic ones with GAS
    return ftot_tp1, target_fun_val_t, I_t, grad_tot_t
end


function gasFilter( model::T where T<: GasNetModelDirBin0Rec0, vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1}; vConstPar ::Array{<:Real,1} = zeros(Real,2),
                    obsT=model.obsT, ftot_0::Array{<:Real,1} = zeros(Real,2), dgpNT = (0,0))
    """GAS Filter the Dynamic Fitnesses from the Observed degrees, given the GAS parameters
     given T observations for the degrees in TxN vector degsT
     """

     #Per i modelli con pseudolikelihood this funciton allows only filtering

    sum(dgpNT) == 0   ?     dgp = false : dgp = true
    NergmPar = 2#
    NTvPar   = sum(indTvPar)
    if dgp
        T = dgpNT[2]
        N_const = dgpNT[1]
        N_T = ones(Int, T) .* N_const 
    else
        T= length(obsT);
    end
    # Organize parameters of the GAS update equation
    Wvec = vResGasPar[1:3:3*NTvPar]
    Bvec = vResGasPar[2:3:3*NTvPar]
    Avec = vResGasPar[3:3:3*NTvPar]


    # start values equal the unconditional mean,and  constant ones remain equal to the unconditional mean, hence initialize as:
    UMallPar = zeros(Real,NergmPar)
    UMallPar[indTvPar] =  Wvec ./ (1 .- Bvec)
    if !all(indTvPar) # if not all parameters are time varying
        UMallPar[.!indTvPar] = vConstPar
    end
    #println(UMallPar)
    fVecT = ones(Real,NergmPar,T)
    sVecT = ones(Real,NergmPar,T)

    sum(ftot_0)==0  ?    ftot_0 = UMallPar : ()# identify(model,UMallNodesIO)

    ftot_t = copy(ftot_0)
    if NTvPar==0
        I_tm1 = ones(1,1)
    else
        I_tm1 = Float64.(Diagonal{Real}(I,NTvPar))
    end
    loglike = 0
    if dgp
        A_T = zeros(Int8, N, N, T)
    end
    fVecT[:,1] = ftot_t
    for t=2:T
    #    println(t)
        if dgp
            diadProb = StaticNets.diadProbFromPars(StaticNets.fooNetModelDirBin0Rec0, ftot_t )
            A_t = StaticNets.samplSingMatCan(StaticNets.fooNetModelDirBin0Rec0, diadProb, N)
            A_T[:, : , t] = A_t
            obs_t = StaticNets.statsFromMat(StaticNets.fooNetModelDirBin0Rec0, A_t)

        else   
            obs_t = obsT[t]
        end

        #print((t,I_tm1))
        ftot_t,loglike_t,I_tm1,grad_t = updatedGasPar(model,obs_t, ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
        fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
        sVecT[:,t] = grad_t #store the filtered parameters from previous iteration
        loglike += loglike_t
    end
    if dgp
        return fVecT, A_T, sVecT
    else
        return fVecT, loglike, sVecT
    end
end

# Estimation
function setOptionsOptim(model::T where T<: GasNetModelDirBin0Rec0)
    "Set the options for the optimization required in the estimation of the model.
    For the optimization use the Optim package."
    tol = eps()*100
    maxIter = 50
    opt = Optim.Options(  g_tol = tol,
                     x_tol = tol,
                     f_tol = tol,
                     iterations = maxIter,
                     show_trace = true,#false,#
                     show_every=1)

    algo = NewtonTrustRegion(; initial_delta = 0.1,
                    delta_hat = 0.2,
                    eta = 0.1,
                    rho_lower = 0.25,
                    rho_upper = 0.75)
    algo = Newton(; alphaguess = LineSearches.InitialHagerZhang(),
    linesearch = LineSearches.BackTracking())
      return opt, algo
end

function static_estimate(model::GasNetModelDirBin0Rec0_mle, statsT)
    L_mean  = mean([stat[1] for stat in statsT ])
    R_mean  = mean([stat[2] for stat in statsT ])
    N_mean  = mean([stat[3] for stat in statsT ])
    
    println(N_mean)
    staticPars = StaticNets.estimate(StaticNets.fooNetModelDirBin0Rec0, L_mean, R_mean, N_mean )
    return staticPars
end

function estimate(model::T where T<: GasNetModelDirBin0Rec0; obsT = model.obsT, indTvPar::BitArray{1}=trues(2), indTargPar::BitArray{1} = indTvPar, UM:: Array{<:Real,1} = zeros(2), ftot_0 :: Array{<:Real,1} = zeros(2), hess_opt_flag = false)
    "Estimate the GAS and static parameters  "

    T = length(obsT);
    NergmPar = 2 #
    NTvPar = sum(indTvPar)
    NTargPar = sum(indTargPar)

    # UM is a vector with target values for dynamical ones. Parameters
    # if not given as input use the static estimates
    # single static estimate
    staticPars = static_estimate(model, obsT)

    if prod(UM.== 0 )&(!prod(.!indTvPar))
        UM = staticPars
    end

    # ftot_0 is a vector with initial values (to be used in the SD iteration)
    # if not given as input estimate on first 3 observations
    if prod(ftot_0.== 0 )&(!prod(.!indTvPar))
        ftot_0 =  static_estimate(model, obsT[1:10])
    end
    #UM = ftot_0

    optims_opt, algo = setOptionsOptim(model)

    # #set the starting points for the optimizations
    B0_Re  = 0.98; B0_Un = log(B0_Re ./ (1 .- B0_Re ))
    ARe_min =0.000001
    A0_Re  = 0.0005 ; A0_Un = log(A0_Re  .-  ARe_min)
    # starting values for the vector of parameters that have to be optimized

    function initialize_pars(vParOptim_0)
        last = 0
        for i=1:NergmPar
            if indTvPar[i]
                if indTargPar[i]
                    vParOptim_0[last+1:last+2] = [ B0_Un; A0_Un]
                    last+=2
                else
                    vParOptim_0[last+1:last+3] = [UM[i]*(1 .- B0_Re) ; B0_Un; A0_Un]
                    last+=3
                end
            else
                vParOptim_0[last+1] = UM[i]
                last+=1
            end
        end
        return vParOptim_0
    end
    vParOptim_0 = initialize_pars(zeros(NergmPar + NTvPar*2 - NTargPar))
    @show(vParOptim_0)

    function divideCompleteRestrictPar(vecUnPar::Array{<:Real,1})

        # vecUnPar is a vector of unrestricted parameters that need to be optimized.
        # add some elements to take into account targeting, divide into GAs and
        # costant parameters, restrict the parameters to appropriate Utilitiesains
        vecReGasParAll = zeros(Real,3NTvPar )
        vecConstPar = zeros(Real,NergmPar-NTvPar)
        # add w determined by B values to targeted parameters
        lastInputInd = 0
        lastGasInd = 0
        lastConstInd = 0
        #extract the vector of gas parameters, addimng w from targeting when needed
        for i=1:NergmPar
            if indTvPar[i]
                if indTargPar[i]
                    B =  1 ./ (1 .+ exp.( .- vecUnPar[lastInputInd+1]))
                    vecReGasParAll[lastGasInd+1] = UM[i]*(1 .- B) # w
                    vecReGasParAll[lastGasInd+2] = B #B
                    vecReGasParAll[lastGasInd+3] =  ARe_min   .+  exp(vecUnPar[lastInputInd + 2]) # A
                    lastInputInd +=2
                    lastGasInd +=3
                else
                    vecReGasParAll[lastGasInd+1] = vecUnPar[lastInputInd  + 1]
                    vecReGasParAll[lastGasInd+2] =  1 ./ (1 .+ exp.( .- vecUnPar[lastInputInd + 2]))
                    vecReGasParAll[lastGasInd+3] = ARe_min   .+  exp(vecUnPar[lastInputInd + 3])
                    lastInputInd +=3
                    lastGasInd +=3
                end
            else
                vecConstPar[lastConstInd+1] = vecUnPar[lastInputInd  + 1]
                lastInputInd +=1
                lastConstInd +=1
            end
        end
    return vecReGasParAll,vecConstPar
    end
    # objective function for the optimization
    function objfunGas(vecUnPar::Array{<:Real,1})# a function of the groups parameters
        #vecUnGasPar,vecConstPar =  divideParVec(vecUnPar)

        vecReGasParAll,vecConstPar = divideCompleteRestrictPar(vecUnPar)

        oneInADterms  = (StaticNets.maxLargeVal + vecUnPar[1])/StaticNets.maxLargeVal

        foo, target_fun_val_T, foo1 = gasFilter( model,  vecReGasParAll, indTvPar; obsT = obsT, vConstPar =  vecConstPar, ftot_0 = ftot_0 .* oneInADterms)

        #println(vecReGasPar)
         return - target_fun_val_T
    end
    #Run the optimization
    if uppercase(model.scoreScalingType) == "FISHER-EWMA"
        ADobjfunGas = objfunGas
    else
        ADobjfunGas = TwiceDifferentiable(objfunGas, vParOptim_0; autodiff = :forward);
    end

    println(objfunGas(vParOptim_0))
    #error()
    optim_out2  = optimize(ADobjfunGas,vParOptim_0 ,algo,optims_opt)
    outParAllUn = Optim.minimizer(optim_out2)
    vecAllParGasHat, vecAllParConstHat = divideCompleteRestrictPar(outParAllUn)

    @show(vecAllParGasHat)
    @show(vecAllParConstHat)
    function reshape_results(vecAllParGasHat)
        arrayAllParHat = fill(Float64[],NergmPar)
        lastGasInd = 0
        lastConstInd = 0
        for i=1:NergmPar
            if indTvPar[i]
                arrayAllParHat[i] = vecAllParGasHat[lastGasInd+1:lastGasInd+3]
                lastGasInd += 3
            else
                arrayAllParHat[i] = vecAllParConstHat[lastConstInd+1]*ones(1)
                lastConstInd+=1
            end
        end
        return arrayAllParHat
    end

    arrayAllParHat = reshape_results(vecAllParGasHat)
    conv_flag =  Optim.converged(optim_out2)

    # println(optim_out2)
    if hess_opt_flag
        #if required return the hessian for the restricted parameters computed at MLE
        #Total likelihood as a function of the restricted parameters
        function likeFun(vecReGasParAll::Array{<:Real,1})
            oneInADterms  = (maxLargeVal + vecReGasParAll[1])/maxLargeVal
            foo,loglikelValue = gasFilter(model,vecReGasParAll,indTvPar;
                                            obsT = changeStats_T,vConstPar =  vecAllParConstHat,ftot_0 = ftot_0 .* oneInADterms)
            #println(vecReGasPar)
             return - loglikelValue
        end

        likeFun(vecAllParGasHat)
        print( ForwardDiff.gradient(likeFun,vecAllParGasHat))
        hess_opt =  ForwardDiff.hessian(likeFun,vecAllParGasHat)
        return  arrayAllParHat, conv_flag,UM , ftot_0 , hess_opt
    else
        return  arrayAllParHat, conv_flag,UM , ftot_0
    end
end

#--------Pseudologlikelihood SDERGM

function change_stats(model::T where T<: GasNetModelDirBin0Rec0, A_T::Array{Matrix{<:Real},1})
    return [StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A) for A in A_T]
end

function change_stats(model::T where T<: GasNetModelDirBin0Rec0, A_T::Array{<:Real,3})
    return [StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A_T[:,:,t]) for t in 1:size(A_T)[3]]
end




function static_estimate(model::GasNetModelDirBin0Rec0_pmle, A_T)
    staticPars = ErgmRcall.get_static_mple(A_T, "edges +  mutual")
    # pmle mean estimate
    return staticPars
end

function target_function_t(model::GasNetModelDirBin0Rec0_pmle, obs_t, par)
    changeStat, response, weights = ErgmRcall.decomposeMPLEmatrix(obs_t)
    return StaticNets.pseudo_loglikelihood_strauss_ikeda( StaticNets.fooNetModelDirBin0Rec0, par, changeStat, response, weights)
end


function dgp_missp(model::GasNetModelDirBin0Rec0, T, θ_0_minMax, η_0_minMax, dgpType)
    parMatDgp_T = zeros(2,T)
    Nsteps1= 2

    if dgpType=="sin"
        phase = 10rand()
        parMatDgp_T[1,:] = dgpSin(θ_0_minMax[1], θ_0_minMax[2], Nsteps1,T; phase = phase)# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
        parMatDgp_T[2,:] = dgpSin(η_0_minMax[1], η_0_minMax[2], Nsteps1,T;phase= phase + 2/T)# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
    elseif dgpType=="steps"
        parMatDgp_T[1,:] = randSteps(θ_0_minMax[1], θ_0_minMax[2], Nsteps1,T)
        parMatDgp_T[2,:] = randSteps(η_0_minMax[1], η_0_minMax[2], Nsteps1,T)
    elseif dgpType=="AR"
        B = 0.95
        sigma = 0.1
        parMatDgp_T[1,:] = dgpAR(mean(θ_0_minMax),B,sigma,T; minMax=θ_0_minMax )
        parMatDgp_T[2,:] = dgpAR(mean(η_0_minMax),B,sigma,T; minMax = η_0_minMax )
    end
    return parMatDgp_T
end


function sample_dgp(model::GasNetModelDirBin0Rec0, parMatDgp_T::Matrix, N )
    T = size(parMatDgp_T)[2]
    A_T_dgp = zeros(Int8, N, N, T)
    for t=1:T
        diadProb = StaticNets.diadProbFromPars(StaticNets.fooNetModelDirBin0Rec0, parMatDgp_T[:,t] )
        A_T_dgp[:,:,t] = StaticNets.samplSingMatCan(StaticNets.fooNetModelDirBin0Rec0, diadProb, N)
    end
    return A_T_dgp
end




function sample_est_mle_pmle(model::GasNetModelDirBin0Rec0, parMatDgp_T, N, Nsample; plotFlag = true, regimeString="") 
    model_mle = fooGasNetModelDirBin0Rec0_mle
    model_pmle = fooGasNetModelDirBin0Rec0_pmle

    T = size(parMatDgp_T)[2]
    indTvPar = trues(2)
   
    if plotFlag
        fig1, ax_mle = subplots(2,1)
        fig2, ax_pmle = subplots(2,1)
    end

    rmse(x::Matrix) = sqrt.(mean(x.^2, dims=2))

    rmse_mle = zeros(2,Nsample)
    rmse_pmle = zeros(2,Nsample)
    vEstSd_mle = zeros(6, Nsample)
    vEstSd_pmle = zeros(6, Nsample)

    for n=1:Nsample
        ## sample dgp
        A_T_dgp = sample_dgp(model_mle, parMatDgp_T, N)
        stats_T_dgp = [statsFromMat(model_mle, A_T_dgp[:,:,t]) for t in 1:T ]
        change_stats_T_dgp = change_stats(model_pmle, A_T_dgp)


        ## estimate SD
        estPar_pmle, conv_flag,UM_mple , ftot_0_mple = estimate(model_pmle; indTvPar=indTvPar,indTargPar=indTvPar, obsT = change_stats_T_dgp)
        vResEstPar_pmle = DynNets.array2VecGasPar(model_pmle, estPar_pmle, indTvPar)
        fVecT_filt_p , target_fun_val_T_p, sVecT_filt_p = gasFilter( model_pmle,  vResEstPar_pmle, indTvPar; obsT = change_stats_T_dgp, ftot_0 = ftot_0_mple)
        vEstSd_pmle[:,n] = vResEstPar_pmle

        estPar_mle, conv_flag,UM_mle , ftot_0_mle = estimate(model_mle; indTvPar=indTvPar, indTargPar=indTvPar, obsT = stats_T_dgp)
        vResEstPar_mle = DynNets.array2VecGasPar(model_mle, estPar_mle, indTvPar)
        fVecT_filt , target_fun_val_T, sVecT_filt = gasFilter( model_mle,  vResEstPar_mle, indTvPar; obsT = stats_T_dgp, ftot_0=ftot_0_mle)
        vEstSd_mle[:,n] = vResEstPar_mle

        if plotFlag
            ax_mle[1].plot(fVecT_filt[1,:], "b", alpha =0.5)
            ax_mle[2].plot(fVecT_filt[2,:], "b", alpha =0.5)
            ax_pmle[1].plot(fVecT_filt_p[1,:], "r", alpha =0.5)
            ax_pmle[2].plot(fVecT_filt_p[2,:], "r", alpha =0.5)
        end

        rmse_mle[:,n] = rmse(fVecT_filt.- parMatDgp_T)
        rmse_pmle[:,n] =rmse(fVecT_filt_p.- parMatDgp_T)
    end
    
    avg_rmse_pmle = round.(mean(drop_nan_col(rmse_pmle), dims=2), digits=3)
    avg_rmse_mle = round.(mean(drop_nan_col(rmse_mle), dims=2), digits=3)

    if plotFlag
        ax_mle[1].plot(parMatDgp_T[1,:], "k")
        ax_mle[2].plot(parMatDgp_T[2,:], "k")
        ax_pmle[1].plot(parMatDgp_T[1,:], "k")
        ax_pmle[2].plot(parMatDgp_T[2,:], "k")

        ax_mle[1].set_title("MLE-SDERGM  N= $N , θ rmse = $(avg_rmse_mle[1]) " * regimeString)   
        ax_mle[2].set_title("MLE-SDERGM  N= $N , η rmse = $(avg_rmse_mle[2]) " * regimeString)   
        
        ax_pmle[1].set_title("PMLE-SDERGM  N= $N , θ rmse = $(avg_rmse_pmle[1])  " * regimeString)   
        ax_pmle[2].set_title("PMLE-SDERGM  N= $N , η rmse = $(avg_rmse_pmle[2])  " * regimeString)   
                
        fig1.tight_layout()
        fig2.tight_layout()
    end
    
    return (vEstSd_mle=vEstSd_mle, vEstSd_pmle=vEstSd_pmle, rmse_mle=avg_rmse_mle, rmse_pmle=avg_rmse_pmle)
end

