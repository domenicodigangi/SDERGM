
Base.@kwdef struct  GasNetModelDirBin0Rec0_mle <: GasNetModelDirBin0Rec0
     indTvPar :: BitArray{1} = trues(2) #  what parameters are time varying   ?
     scoreScalingType::String = "" # String that specifies the rescaling of the score. For a list of possible choices see function scalingMatGas
end
export GasNetModelDirBin0Rec0_mle



name(x::GasNetModelDirBin0Rec0_mle) = "GasNetModelDirBin0Rec0_mle"
export name

statsFromMat(Model::GasNetModelDirBin0Rec0_mle, A ::Matrix{<:Real}) = StaticNets.statsFromMat(StaticNets.fooNetModelDirBin0Rec0, A ::Matrix{<:Real}) 


StaModType(model::T where T<: GasNetModelDirBin0Rec0 ) = StaticNets.fooNetModelDirBin0Rec0# to be substituted with a conversion mechanism



array2VecGasPar(model::T where T<: GasNetModelDirBin0Rec0, ArrayGasPar, indTvPar :: BitArray{1}) = array2VecGasPar(fooGasNetModelDirBinGlobalPseudo, ArrayGasPar, indTvPar )


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
    ARe_min =0.00000000001
    A0_Re  = 0.000005 ; A0_Un = log(A0_Re  .-  ARe_min)
    
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


"""
Given vecAllPar divide it into a vector of Score Driven parameters and one of costant parameters
"""
function divide_SD_par_from_const(model::T where T <:GasNetModel, indTvPar,  vecAllPar::Array{<:Real,1})

    nTvPar = sum(indTvPar)
    nErgmPar = length(indTvPar)

    vecSDParAll = zeros(Real,3nTvPar )
    vConstPar = zeros(Real,nErgmPar-nTvPar)

    lastInputInd = 0
    lastIndSD = 0
    lastConstInd = 0
    #extract the vector of gas parameters, addimng w from targeting when needed
    for i=1:nErgmPar
        if indTvPar[i] 
            vecSDParAll[lastIndSD+1] = vecAllPar[lastInputInd + 1]
            vecSDParAll[lastIndSD+2] = vecAllPar[lastInputInd + 2]
            vecSDParAll[lastIndSD+3] = vecAllPar[lastInputInd + 3]
            lastInputInd +=3
            lastIndSD +=3
        else
            vConstPar[lastConstInd+1] = vecAllPar[lastInputInd  + 1]
            lastInputInd +=1
            lastConstInd +=1
        end
    end
    return vecSDParAll, vConstPar
end


function merge_SD_par_and_const(model::T where T <:GasNetModel, indTvPar,  vecSDPar::Array{<:Real,1}, vConstPar)

    nTvPar = sum(indTvPar)
    nConstPar = sum(.!indTvPar)
    nConstPar == length(vConstPar) ? () : error()

    nErgmPar = length(indTvPar)
    nAllPar = 3*nTvPar + nConstPar

    vecAllPar = zeros(Real, nAllPar)

    lastIndAll = 0
    lastIndSD = 0
    lastIndConst = 0
    for i=1:nErgmPar
        if indTvPar[i] 
            
            vecAllPar[lastIndAll+1] = vecSDPar[lastIndSD + 1]
            vecAllPar[lastIndAll+2] = vecSDPar[lastIndSD + 2]
            vecAllPar[lastIndAll+3] = vecSDPar[lastIndSD + 3]

            lastIndAll +=3
            lastIndSD +=3
        else
            vecAllPar[lastIndAll+1] = vConstPar[lastIndConst + 1]
                        
            lastInputInd +=1
            lastConstInd +=1
        end
    end
    return vecAllPar
end


"""
Restrict the  Score Driven parameters  to appropriate link functions to ensure that they remain in the region where the SD dynamics is well specified (basically 0<=B<1  A>=0)
"""
function restrict_SD_static_par(model::T where T <:GasNetModel, vecUnSDPar::Array{<:Real,1})

    nSDPar = length(vecUnSDPar)
    nTvPar, rem = divrem(nSDPar,3)

    rem == 0 ? () : error()

    arrayOfVecsReSd = [ [vecUnSDPar[i], link_R_in_0_1(vecUnSDPar[i+1]), link_R_in_R_pos(vecUnSDPar[i+2]) ] for i in 1:3:nSDPar]

    vecReSDPar = reduce(vcat, arrayOfVecsReSd)

    return vecReSDPar
end


"""
Restrict the  Score Driven parameters  to appropriate link functions to ensure that they remain in the region where the SD dynamics is well specified (basically 0<=B<1  A>=0)
"""
function unrestrict_SD_static_par(model::T where T <:GasNetModel, vecReSDPar::Array{<:Real,1})

    nSDPar = length(vecReSDPar)
    nTvPar, rem = divrem(nSDPar,3)

    rem == 0 ? () : error()

    arrayOfVecsUnSd = [ [vecReSDPar[i], inv_link_R_in_0_1(vecReSDPar[i+1]), inv_link_R_in_R_pos(vecReSDPar[i+2]) ] for i in 1:3:nSDPar]

    vecUnSDPar = reduce(vcat, arrayOfVecsUnSd)

    return vecUnSDPar
end



function unrestrict_all_par(model, indTvPar, vAllPar)
    vSDRe, vConst = divide_SD_par_from_const(model, indTvPar, vAllPar)

    vSDUn = unrestrict_SD_static_par(model, vSDRe)

    merge_SD_par_and_const(model, indTvPar, vSDUn, vConst)
end


function restrict_all_par(model, indTvPar, vAllPar)
    vSDUn, vConst = divide_SD_par_from_const(model, indTvPar, vAllPar)

    vSDRe = restrict_SD_static_par(model, vSDUn)

    merge_SD_par_and_const(model, indTvPar, vSDRe, vConst)
end

"""
To be filled in in case we want to explore conf bands in targeted models
"""
function target_unc_mean(UM, indTargPar)

end


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


function target_function_t(model::GasNetModelDirBin0Rec0_mle, obs_t, f_t)

    L, R, N = obs_t

    return StaticNets.logLikelihood( StaticNets.fooNetModelDirBin0Rec0, L, R, N, f_t)
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


function updatedGasPar( model::T where T<: GasNetModelDirBin0Rec0, obs_t, ftot_t::Array{<:Real,1}, I_tm1::Array{<:Real,2}, indTvPar::BitArray{1}, Wgas::Array{<:Real,1}, Bgas::Array{<:Real,1}, Agas::Array{<:Real,1};matrixScaling=false)
    
    
    #= likelihood and gradients depend on all the parameters (ftot_t), but
    only the time vaying ones (f_t) are to be updated=#
    

    target_fun_val_t = target_function_t(model, obs_t, ftot_t)
    
    grad_tot_t = target_function_t_grad(model, obs_t, ftot_t)

    # No rescaling
    I_t = I_tm1

    f_t = ftot_t[indTvPar] #Time varying ergm parameters
    s_t = grad_tot_t[indTvPar]
    f_tp1 = Wgas .+ Bgas.* f_t .+ Agas.*s_t
    ftot_tp1 = copy(ftot_t)
    ftot_tp1[indTvPar] = f_tp1 #of all parameters udate dynamic ones with GAS
    return ftot_tp1, target_fun_val_t, I_t, grad_tot_t
end


function score_driven_filter( model::T where T<: GasNetModelDirBin0Rec0, vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1}; vConstPar ::Array{<:Real,1} = zeros(Real,2), obsT=model.obsT, ftot_0::Array{<:Real,1} = zeros(Real,2), dgpNT = (0,0))

    sum(dgpNT) == 0 ? dgp = false : dgp = true

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

    fVecT = ones(Real,NergmPar,T)
    sVecT = ones(Real,NergmPar,T)

    sum(ftot_0)==0 ? ftot_0 = UMallPar : ()# identify(model,UMallNodesIO)
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

    for t=1:T-1
    #    println(t)
        if dgp
            diadProb = StaticNets.diadProbFromPars(StaticNets.fooNetModelDirBin0Rec0, ftot_t )

            A_t = StaticNets.samplSingMatCan(StaticNets.fooNetModelDirBin0Rec0, diadProb, N)
            
            A_T[:, : , t] = A_t
            
            obs_t = statsFromMat(model, A_t)

        else   
            obs_t = obsT[t]
        end

        ftot_t,loglike_t,I_tm1,grad_t = updatedGasPar(model,obs_t, ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
        fVecT[:,t+1] = ftot_t #store the filtered parameters from previous iteration
        sVecT[:,t+1] = grad_t #store the filtered parameters from previous iteration
        loglike += loglike_t
    end
    if dgp
        return fVecT, A_T, sVecT
    else
        return fVecT, loglike, sVecT
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
                     show_trace = true,#false,#
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
    
    println(N_mean)
    staticPars = StaticNets.estimate(StaticNets.fooNetModelDirBin0Rec0, L_mean, R_mean, N_mean )
    return staticPars
end


function estimate(model::T where T<: GasNetModelDirBin0Rec0, obsT; indTvPar::BitArray{1}=trues(2), indTargPar::BitArray{1} = indTvPar, UM:: Array{<:Real,1} = zeros(2), ftot_0 :: Array{<:Real,1} = zeros(2))
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
        ftot_0 =  static_estimate(model, obsT[1:2])
    end
    #UM = ftot_0

    optims_opt, algo = setOptionsOptim(model)


    vParOptim_0, ARe_min = starting_point_optim(model, indTvPar, UM; indTargPar = indTargPar)
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

        foo, target_fun_val_T, foo1 = score_driven_filter( model,  vecReGasParAll, indTvPar; obsT = obsT, vConstPar =  vecConstPar, ftot_0 = ftot_0 .* oneInADterms)

        #println(vecReGasPar)
         return - target_fun_val_T
    end
    #Run the optimization
    if uppercase(model.scoreScalingType) == "FISHER-EWMA"
        ADobjfunGas = objfunGas
    else
        ADobjfunGas = TwiceDifferentiable(objfunGas, vParOptim_0; autodiff = :forward);
    end

    @show objfunGas(vParOptim_0)
    optim_out2  = optimize(ADobjfunGas,vParOptim_0 ,algo,optims_opt)
    outParAllUn = Optim.minimizer(optim_out2)
    vecAllParGasHat, vecAllParConstHat = divideCompleteRestrictPar(outParAllUn)

    @show(optim_out2)
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
 
    return StaticNets.pseudo_loglikelihood_strauss_ikeda( StaticNets.fooNetModelDirBin0Rec0, par, changeStat, response, weights)
end

#endregion

#region misspecified dgps
import ..StaticNets:ergm_par_from_mean_vals

alpha_beta_to_theta_eta(α, β, N) = collect(ergm_par_from_mean_vals(fooNetModelDirBin0Rec0, α*n_pox_dir_links(N), β*n_pox_dir_links(N), N))


theta_eta_to_alpha_beta(θ, η, N) =  collect(exp_val_stats(fooNetModelDirBin0Rec0, θ, η, N))./n_pox_dir_links(N)



get_theta_eta_seq_from_alpha_beta(alpha_beta_seq, N) = reduce(hcat, [alpha_beta_to_theta_eta(ab[1], ab[2], N) for ab in eachcol(alpha_beta_seq)])


function beta_min_max_from_alpha_min(minValAlpha, N; minPairsNumberDiffSupBound = minValAlpha*n_pox_dir_links(N)/5, minPairsNumber = minValAlpha*n_pox_dir_links(N)/5 )
    # min beta val such that the minimum exp value rec links is not close zero   
    minValBeta = minPairsNumber / n_pox_dir_links(N)

    # min beta val such that the maximum exp value rec links is not close to the physical upper bound 
    
    physicalUpperBound = minValAlpha/2
    maxValBeta =  physicalUpperBound - minPairsNumberDiffSupBound /n_pox_dir_links(N) 
    
    minValBeta >= maxValBeta ? error((minValBeta, maxValBeta)) : ()

    return minValBeta, maxValBeta 
end


function dgp_misspecified(model::GasNetModelDirBin0Rec0, dgpType, N, T;  minValAlpha = 0.25, maxValAlpha = 0.3, nCycles = 2, phaseshift = 0.1, plotFlag=false, phaseAlpha = 0, sigma = 0.01, B = 0.95)

    minValBeta, maxValBeta =  beta_min_max_from_alpha_min(minValAlpha, N)
    @show minValBeta, maxValBeta
    @show minValAlpha, maxValAlpha
    
    minValBetaSin, maxValBetaSin = minValBeta, maxValBeta# .* [1,1.5]
    betaConst = minValBeta#(maxValBeta - minValBeta)/2
    #phaseAlpha = rand()  * 2π
    phaseBeta = phaseAlpha + phaseshift * 2π


    α_β_parDgpT = zeros(2,T)
    Nsteps1= 2

    if dgpType=="sin"
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

    if plotFlag
        fig, ax = subplots(2,2)
        ax[1,1].plot(α_β_parDgpT[1,:], "k")
        ax[2,1].plot(α_β_parDgpT[2,:], "k")
    end
    θ_η_parDgpT = get_theta_eta_seq_from_alpha_beta(α_β_parDgpT, N)

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
        estPar_pmle, conv_flag,UM_mple , ftot_0_mple = estimate(model_pmle, change_stats_T_dgp; indTvPar=indTvPar,indTargPar=indTvPar)
        vResEstPar_pmle = DynNets.array2VecGasPar(model_pmle, estPar_pmle, indTvPar)
        fVecT_filt_p , target_fun_val_T_p, sVecT_filt_p = score_driven_filter( model_pmle,  vResEstPar_pmle, indTvPar; obsT = change_stats_T_dgp, ftot_0 = ftot_0_mple)
        vEstSd_pmle[:,n] = vcat(vResEstPar_pmle, ftot_0_mple)

        estPar_mle, conv_flag,UM_mle , ftot_0_mle = estimate(model_mle, stats_T_dgp; indTvPar=indTvPar, indTargPar=indTvPar)
        vResEstPar_mle = DynNets.array2VecGasPar(model_mle, estPar_mle, indTvPar)
        fVecT_filt , target_fun_val_T, sVecT_filt = score_driven_filter( model_mle,  vResEstPar_mle, indTvPar; obsT = stats_T_dgp, ftot_0=ftot_0_mle)
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


function A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model, obsT, vecUnParAll, indTvPar, ftot_0)

    T = length(obsT)
    nPar = length(vecUnParAll)
    gradT = zeros(nPar, T)
    hessT = zeros(nPar, nPar, T)
    
    for t = 2:T
        function obj_fun_t(xUn)

            vecSDParUn, vConstPar = divide_SD_par_from_const(model, indTvPar, xUn)

            vecSDParRe = restrict_SD_static_par(model, vecSDParUn)

            oneInADterms  = (StaticNets.maxLargeVal + vecSDParRe[1])/StaticNets.maxLargeVal

            fVecT_filt, target_fun_val_T, ~ = DynNets.score_driven_filter( model,  vecSDParRe, indTvPar; obsT = obsT[1:t-1], vConstPar =  vConstPar, ftot_0 = ftot_0 .* oneInADterms)
        
            return - DynNets.target_function_t(model, obsT[t-1], fVecT_filt[:,end])
        end


        obj_fun_t(vecUnParAll)

        gradT[:,t] = ForwardDiff.gradient(obj_fun_t, vecUnParAll)
        hessT[:,:,t] =  ForwardDiff.hessian(obj_fun_t, vecUnParAll)
    end

    # function obj_fun_T(xUn)

    #     vecSDParUn, vConstPar = DynNets.divide_SD_par_from_const(model, indTvPar, xUn)

    #     vecSDParRe = DynNets.restrict_SD_static_par(model, vecSDParUn)

    #     oneInADterms  = (StaticNets.maxLargeVal + vecSDParRe[1])/StaticNets.maxLargeVal

    #     fVecT_filt, target_fun_val_T, ~ = DynNets.score_driven_filter( model,  vecSDParRe, indTvPar; obsT = obsT, vConstPar =  vConstPar, ftot_0 = ftot_0 .* oneInADterms)

    #     return - target_fun_val_T
    # end
    # hessCumT =  ForwardDiff.hessian(obj_fun_T, vecUnParAll)
    # HessSum = hessCumT./(T-2)

    OPGradSum = sum([gt * gt' for gt in eachcol(gradT[:,2:end])] )
    HessSum = dropdims(sum(hessT[:,:,2:end], dims=3 ), dims=3)

    return OPGradSum, HessSum
end


function white_estimate_cov_mat_static_sd_par(model, obsT, indTvPar, ftot_0, vEstSdResPar)
    T = length(obsT)
    nErgmPar = number_ergm_par(model)
    errorFlag = false
    # sample parameters in unrestricted space
    vecUnParAll = unrestrict_all_par(model, indTvPar, vEstSdResPar)

    OPGradSum, HessSum = A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model, obsT, vecUnParAll, indTvPar, ftot_0)

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

    mvNormalCov = Symmetric(parCovHatPosDef)
    return mvNormalCov, errorFlag 
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



function var_filtered_par_from_filt_and_par_unc(model, obsT, indTvPar, ftot_0, vEstSdResPar, fVecT_filt; nSample = 1000)

    T = length(obsT)
    nErgmPar = number_ergm_par(model)
    
    
    mvNormalCov, errFlag = white_estimate_cov_mat_static_sd_par(model, obsT, indTvPar, ftot_0, vEstSdResPar)

    # sample parameters in unrestricted space
    vecUnParAll = unrestrict_all_par(model, indTvPar, vEstSdResPar)

     parUncVarianceT = zeros(nErgmPar,T)
    filtUncVarianceT = zeros(nErgmPar,T)
 
    if !errFlag
        sampleUnParAll = rand(MvNormal(zeros(6), mvNormalCov), nSample )

        sampleResParAll = reduce(hcat,[restrict_all_par(model, indTvPar,vecUnParAll.+ sampleUnParAll[:,i]) for i in 1:size(sampleUnParAll)[2]])

        nErgmPar = number_ergm_par(model)
        distribFilteredSD = zeros(nSample, nErgmPar,T)
        filtCovHatSample = zeros(nErgmPar, nSample)

        for n=1:nSample
            vResPar = sampleResParAll[:,n]

            vecSDParRe, vConstPar = divide_SD_par_from_const(model, indTvPar, vResPar)

            distribFilteredSD[n, :, :] , ~, ~ = DynNets.score_driven_filter( model,  vecSDParRe, indTvPar; obsT = obsT, ftot_0=ftot_0, vConstPar=vConstPar)

            BMatSD, AMatSD = divide_in_B_A_mats_as_if_all_TV(model, indTvPar, vResPar)

            filtCovHat = (BMatSD.^(-1)).*AMatSD
            filtCovHat[.!indTvPar] .= 0
            filtCovHatSample[:,n] = filtCovHat
            
        end

        if mean(.!isfinite.(distribFilteredSD)) > 1e-3 
            errFlag = true
        else
            errFlag =false
            distribFilteredSD[isnan.(distribFilteredSD)] .= 0
            fVecT_filt[isnan.(fVecT_filt)] .= 0

            filtCovDiagHatMean = mean(filtCovHatSample, dims=2)

            #for each time compute the variance of the filtered par under the normal distrib of static par
        
            for k=1:nErgmPar
                if indTvPar[k]
                    indAmongTv = sum(indTvPar[1:k-1]) +1 
                    for t=1:T
                        a_t_vec = distribFilteredSD[:,indAmongTv,t]
                        aHat_t = fVecT_filt[k,t] 
                
                        # add filtering and parameter unc
                        parUncVarianceT[k, t] = var(a_t_vec[isfinite.(a_t_vec)] .- aHat_t) 
                        isnan(parUncVarianceT[k, t]) ? (@show a_t_vec; @show aHat_t; error()) : ()

                        filtUncVarianceT[k, t] = filtCovDiagHatMean[indAmongTv]
                    end
                else         
                    indAmongAll = sum(indTvPar[1:k-1])*3 +  sum(.!indTvPar[1:k-1]) + 1                        
                    parUncVarianceT[k, :] .= mvNormalCov[indAmongAll,indAmongAll]
                    filtUncVarianceT[k, t] = 0
                end
            end
        end
    end
    return parUncVarianceT, filtUncVarianceT, errFlag
end


function filter_and_conf_bands(model, A_T_dgp, quantilesVals; indTvPar = model.indTvPar,  plotFlag = false, parDgpT=zeros(2,2))
    
    N = size(A_T_dgp)[1]
    T = size(A_T_dgp)[3]
    obsT = [statsFromMat(model, A_T_dgp[:,:,t]) for t in 1:T ]

    estSdResPar, conv_flag, UM_mple, ftot_0 = estimate(model, obsT; indTvPar=indTvPar, indTargPar=falses(2))


    vEstSdResPar = array2VecGasPar(model, estSdResPar, indTvPar)

    fVecT_filt , target_fun_val_T, sVecT_filt = DynNets.score_driven_filter(model,  vEstSdResPar, indTvPar; obsT = obsT, ftot_0 = ftot_0)
    
    vecUnParAll = unrestrict_all_par(model, indTvPar, vEstSdResPar)

    parUncVarianceT, filtUncVarianceT, errFlag = var_filtered_par_from_filt_and_par_unc(model, obsT, indTvPar, ftot_0, vEstSdResPar, fVecT_filt)
    
    nQuant = length(quantilesVals)
    nBands, r = divrem(nQuant,2)
    r>0 ? error() : ()

    confQuantPar = repeat(fVecT_filt, outer=(1,1,nQuant))
    confQuantParFilt = repeat(fVecT_filt, outer=(1,1,nQuant))

    for p =1:number_ergm_par(model)
        for t=1:T
            confQuantPar[p, t, :] = quantile.(Normal(fVecT_filt[p, t], sqrt(parUncVarianceT[p,t])), quantilesVals)
            confQuantParFilt[p, t, :] = quantile.(Normal(fVecT_filt[p, t], sqrt(parUncVarianceT[p,t] + filtUncVarianceT[p,t])), quantilesVals)
        end
    end

    if plotFlag 
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
            ax[p].set_ylim([bottom - margin, top + margin])
            for b in 1:nBands
                ax[p].fill_between(x, confQuantParFilt[p, :, b], y2 =confQuantParFilt[p,:, end-b+1],color =(0.9, 0.2 , 0.2, 0.1), alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
                ax[p].plot(x, confQuantPar[p, :, b], "-g", alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
                ax[p].plot(x, confQuantPar[p,:, end-b+1], "-g", alpha = 0.2*b/nBands  )#, color='b', alpha=.1)
            end
            ax[p].grid()
            
        end
        ax[1].set_title("$(name(model)), N = $N, T=$T")
        
    end
    return obsT, vEstSdResPar, fVecT_filt, confQuantPar, confQuantParFilt, errFlag
end

function conf_bands_coverage(model, parDgpT, N; nSampleCoverage=100, quantilesVals = [0.975, 0.95, 0.05, 0.025])

    T = size(parDgpT)[2]
    nQuant = length(quantilesVals)
    nBands, check = divrem(nQuant,2)
    check!=0 ? error("quantiles should be eaven to define bands") : ()
    nErgmPar = number_ergm_par(model)
    # obs have different types for different models. storing them might require some additional steps
    #allObsT = Array{Array{Array{Float64,2},1},1}(undef, nSampleCoverage)
    allvEstSdResPar = zeros(3*nErgmPar, nSampleCoverage)
    allfVecT_filt = zeros(nErgmPar, T, nSampleCoverage)
    allConfBandsParFilt = zeros(nErgmPar, T, nQuant, nSampleCoverage)
    allErrFlags = falses(nSampleCoverage)
    allCover = zeros(nErgmPar, T, nBands, nSampleCoverage)

    for k=1:nSampleCoverage
        
        A_T_dgp = sample_dgp(model, parDgpT,N)
        allObsT, allvEstSdResPar[:,k], allfVecT_filt[:,:,k], ~, allConfBandsParFilt[:,:,:,k], allErrFlags[k] = filter_and_conf_bands(model, A_T_dgp, quantilesVals)

        for b in 1:nBands
            for p in 1:nErgmPar 
                for t in 1:T
                    ub = allConfBandsParFilt[p, t, b, k] 
                    lb = allConfBandsParFilt[p, t, end-b+1, k]
                    ub<lb ? error("wrong bands ordering") : ()
                    isCovered = lb <= parDgpT[p, t] <= ub 
                    allCover[p, t, b, k] = isCovered
                end
            end
        end
        
    end

    Logging.@info("The fraction of estimates that resulted in errors is $(mean(allErrFlags)) ")

    return allCover, allvEstSdResPar, allfVecT_filt, allConfBandsParFilt, allErrFlags
end



#endregion