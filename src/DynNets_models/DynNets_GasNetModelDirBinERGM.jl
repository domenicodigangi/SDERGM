


Base.@kwdef struct  GasNetModelDirBinErgmPml <: GasNetModelBin
    staticModel::NetModelDirBinErgmPml
    indTvPar :: BitArray{1}  #  what parameters are time varying   ?
    scoreScalingType::String = "" # String that specifies the rescaling of the score. For a list of possible choices see function scalingMatGas
end
export GasNetModelDirBinErgmPml


GasNetModelDirBinErgmPml(ergmTermsString, indTvPar) = GasNetModelDirBinErgmPml(NetModelDirBinErgmPml(ergmTermsString), indTvPar)

GasNetModelDirBinErgmPml(ergmTermsString) = GasNetModelDirBinErgmPml(ergmTermsString, trues)


function Base.getproperty(x::GasNetModelDirBinErgmPml, p::Symbol)
    if p in [:staticModel, :indTvPar, :scoreScalingType]
        return getfield(x, p)
    else 
        return getproperty(getfield(x,:staticModel), p)
    end
end


name(x::GasNetModelDirBinErgmPml) = "SD" * StaticNets.name(x.staticModel)
export name


statsFromMat(model::GasNetModelDirBinErgmPml, A ::Matrix{<:Real}) = StaticNets.statsFromMat(model.staticModel, A ::Matrix{<:Real}) 



function scalingMatGas(model::T where T<: GasNetModelDirBinErgmPml,expMat::Array{<:Real,2},I_tm1::Array{<:Real,2})
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


function target_function_t(model::GasNetModelDirBinErgmPml, obs_t, f_t)

    L, R, N = obs_t

    return StaticNets.logLikelihood( StaticNets.fooNetModelDirBin0Rec0, L, R, N, f_t)
end


function target_function_t_grad(model::T where T<: GasNetModelDirBinErgmPml, obs_t, f_t)

    target_fun_t(x) = target_function_t(model, obs_t, x)
    
    grad_tot_t = ForwardDiff.gradient(target_fun_t, f_t)

    return grad_tot_t
end


function target_function_t_hess(model::T where T<: GasNetModelDirBinErgmPml, obs_t, f_t)

    target_fun_t(x) = target_function_t(model, obs_t, x)
    
    hess_tot_t = ForwardDiff.hessian(target_fun_t, f_t)

    return hess_tot_t
end


function predict_score_driven_par( model::T where T<: GasNetModelDirBinErgmPml, obs_t, ftot_t::Array{<:Real,1}, I_tm1::Array{<:Real,2}, indTvPar::BitArray{1}, Wgas::Array{<:Real,1}, Bgas::Array{<:Real,1}, Agas::Array{<:Real,1};matrixScaling=false)
    
    
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



function setOptionsOptim(model::T where T<: GasNetModelDirBinErgmPml)
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



function static_estimate(model::GasNetModelDirBinErgmPml, A_T)
    staticPars = ErgmRcall.get_static_mple(A_T, model.ergmTermsString)
    # pmle mean estimate
    return staticPars
end


function estimate(model::T where T<: GasNetModelDirBinErgmPml, obsT; indTvPar::BitArray{1}=trues(2), indTargPar::BitArray{1} = indTvPar, UM:: Array{<:Real,1} = zeros(2), ftot_0 :: Array{<:Real,1} = zeros(2))
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


function change_stats(model::T where T<: GasNetModelDirBinErgmPml, A_T::Array{Matrix{<:Real},1})
    return [StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A) for A in A_T]
end


function change_stats(model::T where T<: GasNetModelDirBinErgmPml, A_T::Array{<:Real,3})
    return [StaticNets.change_stats(StaticNets.fooNetModelDirBin0Rec0, A_T[:,:,t]) for t in 1:size(A_T)[3]]
end


function target_function_t(model::GasNetModelDirBinErgmPml, obs_t, par)
 
    changeStat, response, weights = ErgmRcall.decomposeMPLEmatrix(obs_t)
 
    return StaticNets.pseudo_loglikelihood_strauss_ikeda( StaticNets.fooNetModelDirBin0Rec0, par, changeStat, response, weights)
end



#region Uncertainties filtered parameters


function A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model::GasNetModelDirBinErgmPml, obsT, vecUnParAll, indTvPar, ftot_0)

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


function white_estimate_cov_mat_static_sd_par(model::GasNetModelDirBinErgmPml, obsT, indTvPar, ftot_0, vEstSdResPar)
    T = length(obsT)
    nErgmPar = model.nErgmPar
    
    # sample parameters in unrestricted space
    vecUnParAll = unrestrict_all_par(model, indTvPar, vEstSdResPar)

    OPGradSum, HessSum = A0_B0_est_for_white_cov_mat_obj_SD_filter_time_seq(model, obsT, vecUnParAll, indTvPar, ftot_0)

    parCovHat = pinv(HessSum) * OPGradSum * pinv(HessSum)
    
    parCovHatPosDef, minEigenVal = make_pos_def(Symmetric(parCovHat))
    mvNormalCov = Symmetric(parCovHatPosDef)
    return mvNormalCov, minEigenVal 
end


function divide_in_B_A_mats_as_if_all_TV(model::GasNetModelDirBinErgmPml, indTvPar, vEstSdResPar)
  
    nTvPar = sum(indTvPar)
    nErgmPar = model.nErgmPar
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


function conf_bands_par_uncertainty_naive(model::GasNetModelDirBinErgmPml, obsT, vecUnParAll, indTvPar, ftot_0, vEstSdResPar; nSample = 500, quantilesVals = [0.975, 0.025])
 
    T = length(obsT)
    nErgmPar = model.nErgmPar
    
    
    mvNormalCov, minEigenVal = white_estimate_cov_mat_static_sd_par(model, obsT, indTvPar, ftot_0, vEstSdResPar)
    
    confBand = zeros(length(quantilesVals), nErgmPar,T)

    if minEigenVal < -10
        errFlag =true
    else
        errFlag =false
        sampleUnParAll = rand(MvNormal(zeros(6),mvNormalCov), nSample )

        sampleResParAll = reduce(hcat,[restrict_all_par(model, indTvPar,vecUnParAll.+ sampleUnParAll[:,i]) for i in 1:size(sampleUnParAll)[2]])

        nErgmPar = model.nErgmPar
        distribFilteredSD = zeros( nSample, nErgmPar,T)

        for n=1:nSample
             vResPar = sampleResParAll[:,n]

            vecSDParRe, vConstPar = divide_SD_par_from_const(model, indTvPar, vResPar)

            distribFilteredSD[n, :, :] , ~, ~ = DynNets.score_driven_filter( model,  vecSDParRe, indTvPar; obsT = obsT, ftot_0=ftot_0, vConstPar=vConstPar)
        end

        # Compute confidence bands as sequence of quantiles for each tine
        for n=1:nSample
            for t=1:T
                for k=1:nErgmPar
                    filt_t = distribFilteredSD[:,k,t]
                    confBand[:,k,t] = Statistics.quantile(filt_t[.!isnan.(filt_t)],quantilesVals)
                end
            end
        end
    end
    return confBand, errFlag
end


function var_filtered_par_from_filt_and_par_unc(model::GasNetModelDirBinErgmPml, obsT, indTvPar, ftot_0, vEstSdResPar, fVecT_filt; nSample = 1000)

    T = length(obsT)
    nErgmPar = model.nErgmPar
    
    
    mvNormalCov, minEigenVal = white_estimate_cov_mat_static_sd_par(model, obsT, indTvPar, ftot_0, vEstSdResPar)

    # sample parameters in unrestricted space
    vecUnParAll = unrestrict_all_par(model, indTvPar, vEstSdResPar)

     parUncVarianceT = zeros(nErgmPar,T)
    filtUncVarianceT = zeros(nErgmPar,T)
    if minEigenVal < -10
        errFlag =true
    else
        sampleUnParAll = rand(MvNormal(zeros(6), mvNormalCov), nSample )

        sampleResParAll = reduce(hcat,[restrict_all_par(model, indTvPar,vecUnParAll.+ sampleUnParAll[:,i]) for i in 1:size(sampleUnParAll)[2]])

        nErgmPar = model.nErgmPar
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


function filter_and_conf_bands(model::GasNetModelDirBinErgmPml, A_T_dgp, quantilesVals; indTvPar = model.indTvPar,  plotFlag = false, parDgpT=zeros(2,2))
    
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

    for p =1:model.nErgmPar
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

function conf_bands_coverage(model::GasNetModelDirBinErgmPml, parDgpT, N; nSampleCoverage=100, quantilesVals = [0.975, 0.95, 0.05, 0.025])

    T = size(parDgpT)[2]
    nQuant = length(quantilesVals)
    nBands, check = divrem(nQuant,2)
    check!=0 ? error("quantiles should be eaven to define bands") : ()
    nErgmPar = model.nErgmPar
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

    fractErr = mean(allErrFlags)

    return allCover, allvEstSdResPar, allfVecT_filt, allConfBandsParFilt, fractErr
end



#endregion