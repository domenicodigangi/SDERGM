
Base.@kwdef struct  GasNetModelDirBin0Rec0_mle <: GasNetModelDirBin0Rec0
     indTvPar :: BitArray{1} = trues(2) #  what parameters are time varying   ?
     scoreScalingType::String = "" # String that specifies the rescaling of the score. For a list of possible choices see function scalingMatGas
end
export GasNetModelDirBin0Rec0_mle


name(x::GasNetModelDirBin0Rec0_mle) = "GasNetModelDirBin0Rec0_mle"
export name


statsFromMat(Model::GasNetModelDirBin0Rec0_mle, A ::Matrix{<:Real}) = StaticNets.statsFromMat(StaticNets.fooNetModelDirBin0Rec0, A ::Matrix{<:Real}) 


StaModType(model::T where T<: GasNetModelDirBin0Rec0 ) = StaticNets.fooNetModelDirBin0Rec0# to be substituted with a conversion mechanism



number_ergm_par(model::T where T <:GasNetModelDirBin0Rec0) = 2


"""
Given the flag of constant parameters, a starting value for their unconditional means (their constant value, for those constant), return a starting point for the optimization
"""
function starting_point_optim(model::T where T <:GasNetModelDirBin0Rec0, indTvPar, UM; indTargPar =  falses(100))
    
    nTvPar = sum(indTvPar)
    NTargPar = sum(indTargPar)
    nErgmPar = length(indTvPar)
    
    # #set the starting points for the optimizations
    B0_Re  = 0.98; B0_Un = log(B0_Re ./ (1 .- B0_Re ))
    ARe_min =0.00001
    A0_Re  = 0.1 ; A0_Un = log(A0_Re  .-  ARe_min)
    
    # starting values for the vector of parameters that have to be optimized
    vParOptim_0 = zeros(nErgmPar + nTvPar*2 - NTargPar)
    last = 0
    for i=1:nErgmPar
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
    return vParOptim_0, ARe_min
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


function target_function_t(model::GasNetModelDirBin0Rec0_mle, obs_t, f_t)

    L, R, N = obs_t

    ll = StaticNets.logLikelihood( StaticNets.fooNetModelDirBin0Rec0, L, R, N, f_t)

    return ll
end


function target_function_t_grad(model::T where T<: GasNetModelDirBin0Rec0, obs_t, f_t)

    target_fun_t(x) = target_function_t(model, obs_t, x)
    
    grad_tot_t = ForwardDiff.gradient(target_fun_t, f_t)

    return grad_tot_t
end


function target_function_t_hess(model::T where T<: GasNetModelDirBin0Rec0, obs_t, f_t)

    target_fun_t(x) = target_function_t(model, obs_t, x)
    
    hess_tot_t = ForwardDiff.hessian(target_fun_t, f_t)

    return hess_tot_t
end


function updatedGasPar( model::T where T<: GasNetModelDirBin0Rec0, obs_t, ftot_t::Array{<:Real,1}, I_tm1::Array{<:Real,2}, indTvPar::BitArray{1}, Wgas::Array{<:Real,1}, Bgas::Array{<:Real,1}, Agas::Array{<:Real,1})
    
    
    #= likelihood and gradients depend on all the parameters (ftot_t), but
    only the time vaying ones (f_t) are to be updated=#
    
    target_fun_val_t = target_function_t(model, obs_t, ftot_t)
    
    grad_tot_t = target_function_t_grad(model, obs_t, ftot_t)
    
    hess_tot_t = target_function_t_hess(model, obs_t, ftot_t)

    I_updt = -hess_tot_t
    ewmaWeight = 0.95
    I_t = ewmaWeight.*I_updt .+ (1-ewmaWeight).*I(size(I_updt)[1])

    f_t = ftot_t[indTvPar] #Time varying ergm parameters
    s_t = grad_tot_t[indTvPar]./I_t[I(size(I_t)[1])]

    f_tp1 = Wgas .+ Bgas.* f_t .+ Agas.*s_t

    ftot_tp1 = copy(ftot_t)
    ftot_tp1[indTvPar] = f_tp1 #of all parameters udate dynamic ones with GAS
    
    #any(isnan.(ftot_tp1)) ? error() : ()

    return ftot_tp1, target_fun_val_t, I_t, grad_tot_t
end


function score_driven_filter( model::T where T<: GasNetModelDirBin0Rec0, N, vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1}; vConstPar ::Array{<:Real,1} = zeros(Real,2), obsT=zeros(2,2), ftot_0::Array{<:Real,1} = zeros(Real,2), dgpNT = (0,0))

    sum(dgpNT) == 0 ? dgp = false : dgp = true

    nErgmPar = 2#
    NTvPar   = sum(indTvPar)

    if dgp
        T = dgpNT[2]
        N_const = dgpNT[1]
        N_T = ones(Int, T) .* N_const 
        N = N_const
    else
        T= length(obsT);
    end

    # Organize parameters of the GAS update equation
    Wvec = vResGasPar[1:3:3*NTvPar]
    Bvec = vResGasPar[2:3:3*NTvPar]
    Avec = vResGasPar[3:3:3*NTvPar]

    # start values equal the unconditional mean,and  constant ones remain equal to the unconditional mean, hence initialize as:
    UMallPar = zeros(Real,nErgmPar)
    UMallPar[indTvPar] =  Wvec ./ (1 .- Bvec)
    if !all(indTvPar) # if not all parameters are time varying
        UMallPar[.!indTvPar] = vConstPar
    end

    fVecT = ones(Real,nErgmPar,T)
    sVecT = ones(Real,nErgmPar,T)
    scalMatT = ones(Real,nErgmPar, nErgmPar,T)
    logLikeVecT = ones(Real,T)

    sum(ftot_0)==0 ? ftot_0 = UMallPar : ()# identify(model,UMallNodesIO)
    
    
    if NTvPar==0
        I_tm1 = ones(1,1)
    else
        I_tm1 = Float64.(Diagonal{Real}(I,NTvPar))
    end

    loglike = 0

    if dgp
        A_T = zeros(Int8, N, N, T)
    end

    fVecT[:,1] = ftot_0

    ftot_tp1 = zeros(nErgmPar)
    for t=1:T
    #    println(t)
        if dgp
            diadProb = StaticNets.diadProbFromPars(StaticNets.fooNetModelDirBin0Rec0, fVecT[:, t] )

            A_t = StaticNets.samplSingMatCan(StaticNets.fooNetModelDirBin0Rec0, diadProb, N)
            
            A_T[:, : , t] = A_t
            
            obs_t = statsFromMat(model, A_t)

        else   
            obs_t = obsT[t]
        end

        #obj fun at time t is objFun(obs_t, f_t)

        # predictive step
        ftot_tp1, logLikeVecT[t], scalMatT[:,:,t], sVecT[:,t] = updatedGasPar(model, obs_t, fVecT[:, t], I_tm1, indTvPar, Wvec, Bvec, Avec)

        I_tm1 = scalMatT[:,:,t]

        # would the following be the update step instead ??
        #fVecT[:,t], loglike_t, I_tm1, grad_t = updatedGasPar(model,N, obs_t, fVecT[:,t-1], I_tm1, indTvPar, Wvec, Bvec, Avec)

        if t!=T
            fVecT[:,t+1] = ftot_tp1
        end
    end

    #fVecT = hcat(fVecT[:, 2:end], ftot_tp1)

    if dgp
        return fVecT, A_T, sVecT, scalMatT
    else
        return fVecT, logLikeVecT, sVecT, scalMatT
    end
end


function setOptionsOptim(model::T where T<: GasNetModelDirBin0Rec0)
    "Set the options for the optimization required in the estimation of the model.
    For the optimization use the Optim package."
    tol = eps()*10
    maxIter = 150
    opt = Optim.Options(  g_tol = 1e-8,
                     x_tol = tol,
                     x_abstol = tol,
                     x_reltol = tol,
                     f_tol = tol,
                     f_reltol = tol,
                     f_abstol = tol,
                     iterations = maxIter,
                     show_trace = false,#false,#
                     show_every=5)

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
    
    staticPars = StaticNets.estimate(StaticNets.fooNetModelDirBin0Rec0, L_mean, R_mean, N_mean )
    return staticPars
end


function estimate(model::T where T<: GasNetModelDirBin0Rec0, N, obsT; indTvPar::BitArray{1}=trues(2), indTargPar::BitArray{1} = indTvPar, UM:: Array{<:Real,1} = zeros(2), ftot_0 :: Array{<:Real,1} = zeros(2), vParOptim_0 =zeros(2) )
    "Estimate the GAS and static parameters  "

    T = length(obsT);
    nErgmPar = 2 #
    NTvPar = sum(indTvPar)
    NTargPar = sum(indTargPar)
    @show (N, T)

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
        ftot_0 =  static_estimate(model, obsT[1:5])
    end
    #UM = ftot_0

    optims_opt, algo = setOptionsOptim(model)


    vParOptim_0_tmp, ARe_min = starting_point_optim(model, indTvPar, UM; indTargPar = indTargPar)

    if sum(vParOptim_0) == 0
        vParOptim_0 = vParOptim_0_tmp
    end

  
    function divideCompleteRestrictPar(vecUnPar::Array{<:Real,1})

        # vecUnPar is a vector of unrestricted parameters that need to be optimized.
        # add some elements to take into account targeting, divide into GAs and
        # costant parameters, restrict the parameters to appropriate Utilitiesains
        vecReGasParAll = zeros(Real,3NTvPar )
        vecConstPar = zeros(Real,nErgmPar-NTvPar)
        # add w determined by B values to targeted parameters
        lastInputInd = 0
        lastGasInd = 0
        lastConstInd = 0
        #extract the vector of gas parameters, addimng w from targeting when needed
        for i=1:nErgmPar
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

        foo, logLikeVecT, foo1 = score_driven_filter( model, N, vecReGasParAll, indTvPar; obsT = obsT, vConstPar =  vecConstPar, ftot_0 = ftot_0 .* oneInADterms)

         return - sum(logLikeVecT)
    end
    #Run the optimization
    if uppercase(model.scoreScalingType) == "FISHER-EWMA"
        ADobjfunGas = objfunGas
    else
        ADobjfunGas = TwiceDifferentiable(objfunGas, vParOptim_0; autodiff = :forward);
    end

    @show objfunGas(restrict_all_par(model, indTvPar, vParOptim_0))
    @show vParOptim_0
    optim_out2  = optimize(ADobjfunGas,vParOptim_0 ,algo,optims_opt)
    outParAllUn = Optim.minimizer(optim_out2)
    vecAllParGasHat, vecAllParConstHat = divideCompleteRestrictPar(outParAllUn)

    @show(optim_out2)

    @show(vecAllParGasHat)
    @show(vecAllParConstHat)

    function reshape_results(vecAllParGasHat)
        arrayAllParHat = fill(Float64[],nErgmPar)
        lastGasInd = 0
        lastConstInd = 0
        for i=1:nErgmPar
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
   
    return  arrayAllParHat, conv_flag,UM , ftot_0
   
end

#region  Pseudologlikelihood SDERGM

Base.@kwdef struct  GasNetModelDirBin0Rec0_pmle <: GasNetModelDirBin0Rec0
     indTvPar :: BitArray{1} = trues(2) #  what parameters are time varying   ?
     scoreScalingType::String = "" # String that specifies the rescaling of the score. For a list of possible choices see function scalingMatGas
end
export GasNetModelDirBin0Rec0_pmle


name(x::GasNetModelDirBin0Rec0_pmle) = "GasNetModelDirBin0Rec0_pmle"


function change_stats(model::T where T<: GasNetModelDirBin0Rec0, A_T::Array{Matrix{<:Real},1})
    return [StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A) for A in A_T]
end


function change_stats(model::T where T<: GasNetModelDirBin0Rec0, A_T::Array{<:Real,3})
    return [StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A_T[:,:,t]) for t in 1:size(A_T)[3]]
end


statsFromMat(Model::GasNetModelDirBin0Rec0_pmle, A ::Matrix{<:Real}) = StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A)


function static_estimate(model::GasNetModelDirBin0Rec0_pmle, A_T)
    staticPars = ErgmRcall.get_static_mple(A_T, "edges +  mutual")
    # pmle mean estimate
    return staticPars
end


function target_function_t(model::GasNetModelDirBin0Rec0_pmle, obs_t, par)
 
    changeStat, response, weights = ErgmRcall.decomposeMPLEmatrix(obs_t)
 
    pll = StaticNets.pseudo_loglikelihood_strauss_ikeda( StaticNets.fooNetModelDirBin0Rec0, par, changeStat, response, weights)

    return pll
end

#endregion

#region misspecified dgps
import ..StaticNets:ergm_par_from_mean_vals

alpha_beta_to_theta_eta(α, β, N) = collect(ergm_par_from_mean_vals(fooNetModelDirBin0Rec0, α*n_pox_dir_links(N), β*n_pox_dir_links(N), N))


theta_eta_to_alpha_beta(θ, η, N) =  collect(exp_val_stats(fooNetModelDirBin0Rec0, θ, η, N))./n_pox_dir_links(N)



get_theta_eta_seq_from_alpha_beta(alpha_beta_seq, N) = reduce(hcat, [alpha_beta_to_theta_eta(ab[1], ab[2], N) for ab in eachcol(alpha_beta_seq)])


function beta_min_max_from_alpha_min(minValAlpha, N;  minValBeta = minValAlpha/5 )
    # min beta val such that the minimum exp value rec links is not close zero   
     

    # min beta val such that the maximum exp value rec links is not close to the physical upper bound 
    
    physicalUpperBound = minValAlpha/2
    maxValBeta =  physicalUpperBound - minValBeta
    
    minValBeta >= maxValBeta ? error((minValBeta, maxValBeta)) : ()

    return minValBeta, maxValBeta 
end


function dgp_misspecified(model::GasNetModelDirBin0Rec0, dgpType, N, T;  minValAlpha = 0.25, maxValAlpha = 0.3, nCycles = 2, phaseshift = 0.1, plotFlag=false, phaseAlpha = 0, sigma = 0.01, B = 0.95, A=0.001, maxAttempts = 5000, indTvPar=trues(number_ergm_par(model)))

    minValBeta, maxValBeta =  beta_min_max_from_alpha_min(minValAlpha, N)


    @show minValBeta, maxValBeta
    @show minValAlpha, maxValAlpha
    
    minValBetaSin, maxValBetaSin = minValBeta, maxValBeta# .* [1,1.5]
    betaConst = minValBeta#(maxValBeta - minValBeta)/2
    #phaseAlpha = rand()  * 2π
    phaseBeta = phaseAlpha + phaseshift * 2π

    α_β_parDgpT = zeros(2,T)
    
    for n=1:maxAttempts
        Nsteps1= 2
        if dgpType=="SIN"
            α_β_parDgpT[1,:] = dgpSin(minValAlpha, maxValAlpha, nCycles, T; phase = phaseAlpha)# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
            α_β_parDgpT[2,:] .= dgpSin(minValBetaSin, maxValBetaSin, nCycles, T;phase= phaseBeta )# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
        elseif dgpType=="steps"
            α_β_parDgpT[1,:] = randSteps(α_β_minMax[1], α_β_minMax[2], Nsteps1,T)
            α_β_parDgpT[2,:] = randSteps(η_0_minMax[1], η_0_minMax[2], Nsteps1,T)
        elseif dgpType=="AR"
            α_β_minMax = zeros(2)
            meanValAlpha = (minValAlpha+maxValAlpha)/2
            meanValBeta = (minValBeta+maxValBeta)/2
            α_β_parDgpT[1,:] = dgpAR(meanValAlpha,B,sigma,T; minMax = α_β_minMax )
            α_β_parDgpT[2,:] = dgpAR(meanValBeta,B,sigma,T; minMax = α_β_minMax )
        end

        if dgpType == "SD"
            α_β_minMax = zeros(2)
            meanValAlpha = (minValAlpha+maxValAlpha)/2
            meanValBeta = (minValBeta+maxValBeta)/2
    
            UM = [meanValAlpha, meanValBeta]

            vUnPar, ~ = DynNets.starting_point_optim(model, indTvPar, UM)
            vResParDgp = DynNets.restrict_all_par(model, indTvPar, vUnPar)

            vResParDgp[2:3:end] .= B
            vResParDgp[3:3:end].= A

            fVecT, ~, ~ = DynNets.score_driven_filter( model, N, vResParDgp, indTvPar; dgpNT = (N,T))

            global θ_η_parDgpT = fVecT

            break

        else
            try 
                global θ_η_parDgpT = get_theta_eta_seq_from_alpha_beta(α_β_parDgpT, N)
                break
            catch
                Logging.@warn("Dgp attempt failed")
            end
        end

    end

    if plotFlag
        fig, ax = subplots(2,2)
        ax[1,1].plot(α_β_parDgpT[1,:], "k")
        ax[2,1].plot(α_β_parDgpT[2,:], "k")
    end
    if plotFlag
        ax[1,2].plot(θ_η_parDgpT[1,:], "k")
        ax[2,2].plot(θ_η_parDgpT[2,:], "k")
    end
    return θ_η_parDgpT
end


function sample_dgp(model::GasNetModelDirBin0Rec0, parDgpT::Matrix, N )
    T = size(parDgpT)[2]
    A_T_dgp = zeros(Int8, N, N, T)
    for t=1:T
        diadProb = StaticNets.diadProbFromPars(StaticNets.fooNetModelDirBin0Rec0, parDgpT[:,t] )
        A_T_dgp[:,:,t] = StaticNets.samplSingMatCan(StaticNets.fooNetModelDirBin0Rec0, diadProb, N)
    end
    return A_T_dgp
end


function sample_est_mle_pmle(model::GasNetModelDirBin0Rec0, parDgpT, N, Nsample; plotFlag = true, regimeString="")

    model_mle = GasNetModelDirBin0Rec0_mle()
    model_pmle = GasNetModelDirBin0Rec0_pmle()

    T = size(parDgpT)[2]
    indTvPar = trues(2)
   
    if plotFlag
        fig1, ax_mle = subplots(2,1)
        fig2, ax_pmle = subplots(2,1)
    end

    rmse(x::Matrix) = sqrt.(mean(x.^2, dims=2))

    rmse_mle = zeros(2,Nsample)
    rmse_pmle = zeros(2,Nsample)
    vEstSd_mle = zeros(8, Nsample)
    vEstSd_pmle = zeros(8, Nsample)

    for n=1:Nsample
        ## sample dgp
        A_T_dgp = sample_dgp(model_mle, parDgpT, N)
        stats_T_dgp = [statsFromMat(model_mle, A_T_dgp[:,:,t]) for t in 1:T ]
        change_stats_T_dgp = change_stats(model_pmle, A_T_dgp)


        ## estimate SD
        estPar_pmle, conv_flag,UM_mple , ftot_0_mple = estimate(model_pmle, N,  change_stats_T_dgp; indTvPar=indTvPar,indTargPar=indTvPar)
        vResEstPar_pmle = DynNets.array2VecGasPar(model_pmle, estPar_pmle, indTvPar)
        fVecT_filt_p , logLikeVecT, sVecT_filt_p = score_driven_filter( model_pmle, N, vResEstPar_pmle, indTvPar; obsT = change_stats_T_dgp, ftot_0 = ftot_0_mple)
        vEstSd_pmle[:,n] = vcat(vResEstPar_pmle, ftot_0_mple)

        estPar_mle, conv_flag,UM_mle , ftot_0_mle = estimate(model_mle, N, stats_T_dgp; indTvPar=indTvPar, indTargPar=indTvPar)
        vResEstPar_mle = DynNets.array2VecGasPar(model_mle, estPar_mle, indTvPar)
        fVecT_filt , logLikeVecT, sVecT_filt = score_driven_filter( model_mle, N, vResEstPar_mle, indTvPar; obsT = stats_T_dgp, ftot_0=ftot_0_mle)
        vEstSd_mle[:,n] = vcat(vResEstPar_mle, ftot_0_mle)

        if plotFlag
            ax_mle[1].plot(fVecT_filt[1,:], "b", alpha =0.5)
            ax_mle[2].plot(fVecT_filt[2,:], "b", alpha =0.5)
            ax_pmle[1].plot(fVecT_filt_p[1,:], "r", alpha =0.5)
            ax_pmle[2].plot(fVecT_filt_p[2,:], "r", alpha =0.5)
        end

        rmse_mle[:,n] = rmse(fVecT_filt.- parDgpT)
        rmse_pmle[:,n] =rmse(fVecT_filt_p.- parDgpT)
    end
    
    avg_rmse_pmle = round.(mean(drop_nan_col(rmse_pmle), dims=2), digits=3)
    avg_rmse_mle = round.(mean(drop_nan_col(rmse_mle), dims=2), digits=3)

    if plotFlag
        ax_mle[1].plot(parDgpT[1,:], "k")
        ax_mle[2].plot(parDgpT[2,:], "k")
        ax_pmle[1].plot(parDgpT[1,:], "k")
        ax_pmle[2].plot(parDgpT[2,:], "k")

        ax_mle[1].set_title("MLE-SDERGM  N= $N , θ rmse = $(avg_rmse_mle[1]) " * regimeString)   
        ax_mle[2].set_title("MLE-SDERGM  N= $N , η rmse = $(avg_rmse_mle[2]) " * regimeString)   
        
        ax_pmle[1].set_title("PMLE-SDERGM  N= $N , θ rmse = $(avg_rmse_pmle[1])  " * regimeString)   
        ax_pmle[2].set_title("PMLE-SDERGM  N= $N , η rmse = $(avg_rmse_pmle[2])  " * regimeString)   
                
        fig1.tight_layout()
        fig2.tight_layout()
    end
    
    return (; vEstSd_mle, vEstSd_pmle, avg_rmse_mle, avg_rmse_pmle)
end

#endregion

#region Uncertainties filtered parameters


function A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model::GasNetModelDirBin0Rec0, N, obsT, vEstSdResPar, indTvPar, ftot_0)

    T = length(obsT)
    nPar = length(vEstSdResPar)
    gradT = zeros(nPar, T)
    hessT = zeros(nPar, nPar, T)

    vecUnParAll = unrestrict_all_par(model, indTvPar, vEstSdResPar)    


    for t = 1:T
        function obj_fun_t(xUn)

            xRe = restrict_all_par(model, indTvPar, xUn)

            vecSDParRe, vConstPar = divide_SD_par_from_const(model, indTvPar, xRe)

            oneInADterms  = (StaticNets.maxLargeVal + vecSDParRe[1])/StaticNets.maxLargeVal

            fVecT_filt, logLikeVecT, ~ = DynNets.score_driven_filter( model, N,  vecSDParRe, indTvPar; obsT = obsT[1:t], vConstPar =  vConstPar, ftot_0 = ftot_0 .* oneInADterms)
        
            return - logLikeVecT[end]
        end


        obj_fun_t(vecUnParAll)

        gradT[:,t] = deepcopy(ForwardDiff.gradient(obj_fun_t, vecUnParAll))
        hessT[:,:,t] =  deepcopy(ForwardDiff.hessian(obj_fun_t, vecUnParAll))
    end

    # function obj_fun_T(xUn)

    #     vecSDParUn, vConstPar = DynNets.divide_SD_par_from_const(model, indTvPar, xUn)

    #     vecSDParRe = DynNets.restrict_SD_static_par(model, vecSDParUn)

    #     oneInADterms  = (StaticNets.maxLargeVal + vecSDParRe[1])/StaticNets.maxLargeVal

    #     fVecT_filt, target_fun_val_T, ~ = DynNets.score_driven_filter( model, N,  vecSDParRe, indTvPar; obsT = obsT, vConstPar =  vConstPar, ftot_0 = ftot_0 .* oneInADterms)

    #     return - target_fun_val_T
    # end
    # hessCumT =  ForwardDiff.hessian(obj_fun_T, vecUnParAll)
    # HessSum = hessCumT./(T-2)

    OPGradSum = sum([gt * gt' for gt in eachcol(gradT[:,2:end])] )
    HessSum = dropdims(sum(hessT[:,:,2:end], dims=3 ), dims=3)

    return OPGradSum, HessSum
end


function white_estimate_cov_mat_static_sd_par(model::GasNetModelDirBin0Rec0,  N,obsT, indTvPar, ftot_0, vEstSdResPar; returnAllMats=false)

    T = length(obsT)
    nErgmPar = number_ergm_par(model)
    errorFlag = false

   
    OPGradSum, HessSum = A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model, N, obsT, vEstSdResPar, indTvPar, ftot_0)

    parCovHat = pinv(HessSum) * OPGradSum * pinv(HessSum)
    
    parCovHatPosDef, minEigenVal = make_pos_def(parCovHat)
   
    if minEigenVal < 0 
        minEiegOPGrad = minimum(eigen(OPGradSum).values)
        minEiegHess = minimum(eigen(HessSum).values)
        Logging.@info("Outer Prod minimum eigenvalue $(minEiegOPGrad) , hessian minimum eigenvalue $(minEiegHess)")

        # if the negative eigenvalue of the cov mat is due to a negative eigenvalue of the hessian, do not use that estimate
        if minEiegHess < 0 
            errorFlag = true
        end
    end

    mvSDUnParEstCov = Symmetric(parCovHatPosDef)
    if returnAllMats
        return mvSDUnParEstCov, errorFlag, OPGradSum, HessSum
    else
        return mvSDUnParEstCov, errorFlag 
    end
end


function divide_in_B_A_mats_as_if_all_TV(model::GasNetModelDirBin0Rec0, indTvPar, vEstSdResPar)
  
    nTvPar = sum(indTvPar)
    nErgmPar = number_ergm_par(model)
    B = zeros(nTvPar)
    A = zeros(nTvPar)
    
    lastInd = 0
    for e in 1:nErgmPar
        if indTvPar[e]
            bInd = lastInd+2
            aInd = lastInd+3
            B[e] = vEstSdResPar[bInd] 
            A[e] = vEstSdResPar[aInd] 
            lastInd += 3
        else
            lastInd += 1
        end

    end
    return B, A
end


function distrib_filtered_par_from_mv_normal(model::GasNetModelDirBin0Rec0, N, obsT, indTvPar, ftot_0, vEstSdResPar, mvSDUnParEstCov; nSample = 1000)
        
    T = length(obsT)
    nErgmPar = number_ergm_par(model)
    nTvPar = sum(indTvPar)
    any(.!indTvPar) ? error("assuming all parameters TV. otherwise need to set filter unc to zero for static ones") : ()    


    distribFilteredSD = zeros(nSample, nErgmPar,T)
    filtCovHatSample = zeros(nErgmPar, T, nSample)

    vEstSdUnPar = unrestrict_all_par(model, indTvPar, vEstSdResPar)
    # @show vEstSdUnPar
    # @show mvSDUnParEstCov

    zeroMeanSample =  rand(MvNormal(mvSDUnParEstCov), nSample ) 
    sampleUnParAll = zeroMeanSample .+ vEstSdUnPar

    for n=1:nSample
        vResPar = restrict_all_par(model, indTvPar, sampleUnParAll[:,n])

        vecSDParRe, vConstPar = divide_SD_par_from_const(model, indTvPar, vResPar)

        distribFilteredSD[n, :, :] , ~, ~, scalMatT = DynNets.score_driven_filter( model, N,  vecSDParRe, indTvPar; obsT = obsT, ftot_0=ftot_0, vConstPar=vConstPar)

        BMatSD, AMatSD = divide_in_B_A_mats_as_if_all_TV(model, indTvPar, vResPar)

        constFiltUncCoeff = (BMatSD.^(-1)).*AMatSD
        
        for t in 1:T 
            filtCovHatSample[:,t,n] = constFiltUncCoeff./ scalMatT[:,:,t][I(nTvPar)]
        end

    end

    # remove samples that resulted in inconsistent P_t
    indBadSamples = [any(filtCovHatSample[:,:,n].<0) | any(.!isfinite.(filtCovHatSample[:,:,n])) for n in 1:nSample]

    @show mean(indBadSamples)

    if mean(indBadSamples) > 0.5 
        Logging.@error("Sampling from mvNormal resulted in $(mean(indBadSamples)*100) % of bad samples. Rejected") 
        errFlag = true
    else
        errFlag = false
    end

    filtCovHatSampleGood = filtCovHatSample[:,:,.!indBadSamples]

    return distribFilteredSD, filtCovHatSampleGood, mvSDUnParEstCov, errFlag
end


"""
Use parametric bootstrap to sample from the distribution of static parameters of the SD filter
"""
function par_bootstrap_distrib_filtered_par(model::GasNetModelDirBin0Rec0, N, obsT, indTvPar, ftot_0, vEstSdResPar; nSample = 50)

    T = length(obsT)
    nErgmPar = number_ergm_par(model)
     
     
    nStaticSdPar = length(vEstSdResPar)
    vUnParEstDistrib = zeros(length(vEstSdResPar), nSample)
    distribFilteredSD = zeros(nSample, nErgmPar,T)
    filtCovHatSample = zeros(nErgmPar, nSample)
    errFlagVec = falses(nSample)

    for n=1:nSample
        
        # try
            # sample SD dgp
            fVecTDgp, A_T_dgp, ~ = DynNets.score_driven_filter( model, N, vEstSdResPar, indTvPar; dgpNT = (N,T))
            obsTNew = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

            # Estimate SD static params from SD DGP
            arrayAllParHat, conv_flag,UM, ftot_0 = estimate(model, N, obsTNew; indTvPar=indTvPar, indTargPar=falses(2), ftot_0=ftot_0)

            vResPar = DynNets.array2VecGasPar(model, arrayAllParHat, indTvPar)
    
            vUnPar = unrestrict_all_par(model, indTvPar, deepcopy(vResPar))

            vUnParEstDistrib[:,n] =  vUnPar

            vecSDParRe, vConstPar = divide_SD_par_from_const(model, indTvPar, vResPar)

            distribFilteredSD[n, :, :] , ~, ~ = DynNets.score_driven_filter( model, N,  vecSDParRe, indTvPar; obsT = obsT, ftot_0=ftot_0, vConstPar=vConstPar)

            BMatSD, AMatSD = divide_in_B_A_mats_as_if_all_TV(model, indTvPar, vResPar)

            filtCovHat = (BMatSD.^(-1)).*AMatSD
            filtCovHat[.!indTvPar] .= 0
            filtCovHatSample[:,n] = filtCovHat
            
        # catch
        #     errFlagVec[n] = true
        # end

     
    end

    mvSDUnParEstCov = cov(vUnParEstDistrib')

    return distribFilteredSD, filtCovHatSample, errFlagVec, mvSDUnParEstCov
end



function conf_bands_buccheri(model::GasNetModelDirBin0Rec0, obsT, indTvPar, fVecT_filt, distribFilteredSD, filtCovHatSample, quantilesVals::Vector{Vector{Float64}}; nSample = 500, )


    T = length(obsT)
    nErgmPar = number_ergm_par(model)

    parUncVarianceT = zeros(nErgmPar,T)
    filtUncVarianceT = zeros(nErgmPar,T)
 
    distribFilteredSD[isnan.(distribFilteredSD)] .= 0
    fVecT_filt[isnan.(fVecT_filt)] .= 0

    any(isnan.(filtCovHatSample)) ? error() : ()

    filtCovDiagHatMean = mean(filtCovHatSample, dims=3)

    #for each time compute the variance of the filtered par under the normal distrib of static par

    for k=1:nErgmPar
        if indTvPar[k]
            indAmongTv = sum(indTvPar[1:k-1]) +1 
            for t=1:T
                a_t_vec = distribFilteredSD[:,indAmongTv,t]
                aHat_t = fVecT_filt[k,t] 
        
                # add filtering and parameter unc
                diff = a_t_vec[isfinite.(a_t_vec)] .- aHat_t
                extr = quantile(diff, [0.01, 0.99])
                if extr[1] != extr[2]
                    diffNoOutl = filter(x->extr[1] < x < extr[2], diff)
                    parUncVarianceT[k, t] = var(diffNoOutl)
                else    
                    parUncVarianceT[k, t] = var(diff)
                end                
                isnan(parUncVarianceT[k, t]) ? (@show a_t_vec; @show aHat_t; @show diff; @show extr; @show t; @show k; error()) : ()

                filtUncVarianceT[k, t] = filtCovDiagHatMean[indAmongTv,t]
            end
        else         
            indAmongAll = sum(indTvPar[1:k-1])*3 +  sum(.!indTvPar[1:k-1]) + 1                        
            parUncVarianceT[k, :] .= mvSDUnParEstCov[indAmongAll,indAmongAll]
            filtUncVarianceT[k, t] = 0
        end
    end


    nBands = length(quantilesVals)
    nQuant = 2*nBands
  
    confBandsPar = repeat(fVecT_filt, outer=(1,1,nBands, 2))
    confBandsParFilt = repeat(fVecT_filt, outer=(1,1,nBands, 2))

    for p =1:number_ergm_par(model)
        for t=1:T
            for b=1:nBands
                length(quantilesVals[b]) ==2 ? () : error()
                confBandsPar[p, t, b,  :] = quantile.(Normal(fVecT_filt[p, t], sqrt(parUncVarianceT[p,t])), quantilesVals[b])
                try
                    confBandsParFilt[p, t, b,  :] = quantile.(Normal(fVecT_filt[p, t], sqrt(parUncVarianceT[p,t] + filtUncVarianceT[p,t])), quantilesVals[b])
                catch e
                    @show (parUncVarianceT[p,t] , filtUncVarianceT[p,t])
                    @show filtCovDiagHatMean

                error(e)
                end
            end
        end
    end
    return confBandsParFilt, confBandsPar
end


function conf_bands_par_uncertainty_blasques(model::GasNetModelDirBin0Rec0, obsT, fVecT_filt, distribFilteredSD, quantilesVals::Vector{Vector{Float64}}; nSample = 500)
    
    T = length(obsT)
    nErgmPar = number_ergm_par(model)

    parUncVarianceT = zeros(nErgmPar,T)

    nBands = length(quantilesVals)
    nQuant = 2*nBands   
    confBandsPar = repeat(fVecT_filt, outer=(1,1,nBands, 2))

    distribFilteredSD[isnan.(distribFilteredSD)] .= 0
    fVecT_filt[isnan.(fVecT_filt)] .= 0

    #for each time compute the variance of the filtered par under the normal distrib of static par
    
    for k=1:nErgmPar
        for t=1:T
            for b=1:nBands
                length(quantilesVals[b]) ==2 ? () : error()
             
                a_t_vec = distribFilteredSD[:,k,t]
                
                # add filtering and parameter unc
                confBandsPar[k, t, b, :] = Statistics.quantile(a_t_vec[.! isnan.(a_t_vec)], quantilesVals[b])
            end
        end
    end
   
    return confBandsPar
end


function conf_bands_coverage(parDgpT, confBands )

    nErgmPar, T, nBands, nQuant  = size(confBands)
    nQuant ==2 ? () : error()
    nParDgp, TDgp = size(parDgpT)
    nParDgp == nErgmPar ? () : error()
    T == TDgp ? () : error()

    isCovered = falses(nErgmPar, T, nBands)
    for b in 1:nBands
        for p in 1:nErgmPar 
            for t in 1:T
                ub = confBands[p, t, b, 1] 
                lb = confBands[p, t, b, 2]
                ub<lb ? error("wrong bands ordering") : ()
                isCov = lb <= parDgpT[p, t] <= ub 
                isCovered[p, t, b] = isCov
            end
        end
    end
    return isCovered
end

function plot_filtered_and_conf_bands(model::GasNetModelDirBin0Rec0, N, fVecT_filt, confBands1; confBands2 =zeros(2,2), parDgpT=zeros(2,2), nameConfBand1="1", nameConfBand2="2")

    T = size(fVecT_filt)[2]

    nBands = size(confBands1)[3]

    fig1, ax = subplots(2,1)
    for p in 1:2
        x = 1:T
        bottom = minimum(fVecT_filt[p,:])
        top = maximum(fVecT_filt[p,:])
        if parDgpT != zeros(2,2)
            ax[p].plot(x, parDgpT[p,:], "k", alpha =0.5)
            bottom = minimum(parDgpT[p,:])
            top = maximum(parDgpT[p,:])
        end
        delta = top-bottom
        margin = 0.5*delta
        ax[p].plot(x, fVecT_filt[p,:], "b", alpha =0.5)
        # ax[p].set_ylim([bottom - margin, top + margin])
        for b in 1:nBands
            if sum(confBands1) !=0
                ax[p].fill_between(x, confBands1[p, :, b,1], y2 =confBands1[p,:,b, 2],color =(0.9, 0.2 , 0.2, 0.1), alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
            end
            if confBands2 != zeros(2,2)
                ax[p].plot(x, confBands2[p, :, b, 1], "-g", alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
                ax[p].plot(x, confBands2[p,:, b, 2], "-g", alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
            end
        end
        ax[p].grid()
        
    end

    if sum(confBands1) !=0
        cov1 = round(mean(conf_bands_coverage(parDgpT, confBands1)), digits=2)
    else
        cov1 = 0
    end
    titleString = "$(name(model)), N = $N, T=$T, \n  $nameConfBand1 = $cov1"

    if confBands2 != zeros(2,2)
        cov2 = round(mean(conf_bands_coverage(parDgpT, confBands2)), digits=2)
        titleString = titleString * " $nameConfBand2 = $cov2"
    end

    ax[1].set_title(titleString)
end


function filter_and_conf_bands(model::GasNetModelDirBin0Rec0, A_T_dgp, quantilesVals::Vector{Vector{Float64}}; indTvPar = model.indTvPar, parDgpT=zeros(2,2), plotFlag=false, parUncMethod = "WHITE-MLE")
    
    N = size(A_T_dgp)[1]
    T = size(A_T_dgp)[3]
    obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

    estSdResPar, conv_flag, UM_mple, ftot_0 = estimate(model, N, obsT; indTvPar=indTvPar, indTargPar=falses(2))


    vEstSdResPar = array2VecGasPar(model, estSdResPar, indTvPar)

    fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter(model, N,  vEstSdResPar, indTvPar; obsT = obsT, ftot_0 = ftot_0)
    

    if contains(parUncMethod, "PAR-BOOTSTRAP")

        distribFilteredSD, filtCovHatSample, errFlagVec, mvSDUnParEstCov = par_bootstrap_distrib_filtered_par(model, N, obsT, indTvPar, ftot_0, vEstSdResPar)
        
        mean(errFlagVec) > 0.5 ? errFlagEstCov=true : errFlagEstCov = false

        if parUncMethod[14:end] == "SAMPLE"
            
        elseif parUncMethod[14:end] == "COV-MAT"

            distribFilteredSD, filtCovHatSample, mvSDUnParEstCov, errFlagMvNormSample = distrib_filtered_par_from_mv_normal(model, N, obsT, indTvPar, ftot_0, vEstSdResPar, mvSDUnParEstCov)

        end

    elseif parUncMethod == "WHITE-MLE"

         mvSDUnParEstCov, errFlagEstCov = white_estimate_cov_mat_static_sd_par(model, N, obsT, indTvPar, ftot_0, vEstSdResPar)
         
        distribFilteredSD, filtCovHatSample, mvSDUnParEstCov, errFlagMvNormSample = distrib_filtered_par_from_mv_normal(model, N, obsT, indTvPar, ftot_0, vEstSdResPar, mvSDUnParEstCov)
    end

    confBandsFiltPar,  confBandsFiltParPar = conf_bands_buccheri(model, obsT, indTvPar, fVecT_filt, distribFilteredSD, filtCovHatSample, quantilesVals)

    confBandsPar = conf_bands_par_uncertainty_blasques(model, obsT, fVecT_filt, distribFilteredSD, quantilesVals)

    if plotFlag
        plot_filtered_and_conf_bands(model, N, fVecT_filt, confBandsFiltPar ; parDgpT=parDgpT, nameConfBand1= "$parUncMethod - Filter+Par", nameConfBand2= "Par", confBands2=confBandsFiltParPar)

    end

    errFlag = errFlagEstCov & errFlagMvNormSample

    return obsT, vEstSdResPar, fVecT_filt, confBandsFiltPar, confBandsPar, errFlag, vEstSdResPar, mvSDUnParEstCov, distribFilteredSD
end


function conf_bands_coverage(model::GasNetModelDirBin0Rec0, dgpType, dgpOptions, T, N, nSampleCoverage, quantilesVals::Vector{Vector{Float64}}; plotFlag=false)

    nBands = length(quantilesVals)
    nErgmPar = number_ergm_par(model)


    # obs have different types for different models. storing them might require some additional steps
    #allObsT = Array{Array{Array{Float64,2},1},1}(undef, nSampleCoverage)
    allvEstSdResPar = zeros(3*nErgmPar, nSampleCoverage)
    allfVecT_filt = zeros(nErgmPar, T, nSampleCoverage)
    allConfBandsFiltPar = zeros(nErgmPar, T, nBands, 2, nSampleCoverage)
    allConfBandsPar = zeros(nErgmPar, T, nBands, 2, nSampleCoverage)
    allErrFlags = falses(nSampleCoverage)
    allCoverFiltPar = zeros(nErgmPar, T, nBands, nSampleCoverage)
    allCoverPar = zeros(nErgmPar, T, nBands, nSampleCoverage)

    for k=1:nSampleCoverage
        
        parDgpT = DynNets.dgp_misspecified(model, dgpType, N, T;  dgpOptions...)

        A_T_dgp = sample_dgp(model, parDgpT,N)
        

        ~, allvEstSdResPar[:,k], allfVecT_filt[:,:,k], allConfBandsFiltPar[:,:,:,:,k], allConfBandsPar[:,:,:,:,k], allErrFlags[k], ~, ~, ~ = filter_and_conf_bands(model, A_T_dgp, quantilesVals; plotFlag =plotFlag, parDgpT = parDgpT, parUncMethod = "WHITE-MLE")

        
        allCoverFiltPar[:, :, :, k] = conf_bands_coverage(parDgpT,  allConfBandsFiltPar[:,:,:,:,k])

        allCoverPar[:, :, :, k] = conf_bands_coverage(parDgpT,  allConfBandsPar[:,:,:,:,k])
        
    end

    Logging.@info("The fraction of estimates that resulted in errors is $(mean(allErrFlags)) ")

    return allCoverFiltPar, allCoverPar, allvEstSdResPar, allfVecT_filt, allConfBandsFiltPar,allConfBandsPar, allErrFlags
end



#endregion