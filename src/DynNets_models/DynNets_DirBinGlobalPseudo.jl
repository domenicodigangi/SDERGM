#-------------------------- ergm pseudo likelihood

struct  GasNetModelDirBinGlobalPseudo# <: GasNetModelBin
     """ A gas model based on pseudolikelihood (as objective function) for
            Directed binary networks and probability depending on a generic vector
            of global statistics each associated with a time varying parameters.
         """
     obsT::  Array{Array{Float64,2}}# for each t one matrix. Each matrix has the response in the first column
                                # i.e. Adjacency matrix,the predictors (change stats) in the following
                                # columns and, in the last column, the multiplicity of each row, i.e
                                # combination of response and values for the predictors
     Par::Array{Array{<:Real,1},1} #  each time varying parameter has 3 static parameters in this specification.
                         # if a parameter is constant than only 1 gas parameteris present
     indTvPar :: BitArray{1} #  what parameters are time varying   ?

     scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
     # choices see function scalingMatGas
end


#inizializza senza specificare assenza di scaling
GasNetModelDirBinGlobalPseudo( obsT , Par   ) = GasNetModelDirBinGlobalPseudo( obsT,Par ,trues(length(Par[:,1])) ,"")


#Initialize by observations only
GasNetModelDirBinGlobalPseudo(obsT::  Array{Array{Real,2}}) =   (Npar =  length(obsT[1][1,:]) - 2;  GasNetModelDirBinGlobalPseudo(obsT,zeros(Npar,3) ))


GasNetModelDirBinGlobalPseudo(obsT:: Array{Real,3},scoreScalingType::String) =   (N = round(length(obsT[:,1])/2); GasNetModelDirBinGlobalPseudo(obsT,zeros(3Npar),scoreScalingType) )

fooGasPar = [ones(3), ones(3), ones(3)]

fooGasNetModelDirBinGlobalPseudo = GasNetModelDirBinGlobalPseudo(repeat([ones(Float64,10,30)],outer = 5),fooGasPar)

# Relations between Static and Dynamical models: conventions on storage for
# parameters and observations
StaModType(Model::GasNetModelDirBinGlobalPseudo ) = StaticNets.fooNetModelDirBinGlobalPseudo# to be substituted with a conversion mechanism

# options and conversions of parameters for optimization
setOptionsOptim(Model::GasNetModelDirBinGlobalPseudo) = setOptionsOptim(fooGasNetModelDirBin1)


function restrictGasPar(Model::GasNetModelDirBinGlobalPseudo,
                        vecUnGasPar::Array{<:Real,1},
                        indTvPar :: BitArray{1})
     "From the Unrestricted values of the parameters return the restricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, e.g. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     NTvPar = sum(indTvPar)
     NergmPar = length(indTvPar)
     vecReGasPar =  zeros(Real, size(vecUnGasPar))
     last_i = 1
     for p in 1:NergmPar
        vecReGasPar[last_i] = vecUnGasPar[last_i]
        if indTvPar[p]
            vecReGasPar[last_i + 1] = 1 ./ (1 .+ exp.(.-vecUnGasPar[last_i+1]))
            vecReGasPar[last_i + 2] = exp.(vecUnGasPar[last_i+2])
            last_i = last_i + 3
        else
            last_i = last_i + 1
        end
     end

     return vecReGasPar
end


function unRestrictGasPar(Model::GasNetModelDirBinGlobalPseudo,
                        vecReGasPar::Array{<:Real,1},
                        indTvPar :: BitArray{1})
 "From the restricted values of the parameters return the unrestricted ones
 takes as inputs a vector of parameters that can take any value in R and returns
 the appropriately contrained values, i.e. the coefficient of autoregressive
 component in the GAS update rule has to be positive and <1."
 # nei modelli pseudolikelihood ho usato un diverso ordinamento rispetto al resto dei ccasi
     NTvPar = sum(indTvPar)
     NergmPar = length(indTvPar)
     vecUnGasPar =  zeros(Real, size(vecReGasPar))
     last_i = 1
     for p in 1:NergmPar
        vecUnGasPar[last_i] = vecReGasPar[last_i]
        if indTvPar[p]
            vecUnGasPar[last_i + 1] = log.(vecReGasPar[last_i+1] ./ (1 .- vecReGasPar[last_i+1] ))
            vecUnGasPar[last_i + 2] = log.(vecReGasPar[last_i+2])
            last_i = last_i + 3
        else
            last_i = last_i + 1
        end
     end

     return vecUnGasPar
end


#Gas Filter Functions
function identify(Model::GasNetModelDirBinGlobalPseudo,parIO::Array{<:Real,1})
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


function scalingMatGas(Model::GasNetModelDirBinGlobalPseudo,expMat::Array{<:Real,2},I_tm1::Array{<:Real,2})
    "Return the matrix required for the scaling of the score, given the expected
     matrix and the Scaling matrix at previous time. "
    if uppercase(Model.scoreScalingType) == ""
        scalingMat = 1 #
    elseif uppercase(Model.scoreScalingType) == "FISHER-EWMA"

        # λ = 0.8
        #
        # I = expMat.*(1-expMat)
        # diagI = sum(I,dims = 2)
        # [I[i,i] = diagI[i] for i=1:length(diagI) ]
        # I_t =  λ*I + (1-λ) *I_tm1
        # scalingMat = sqrt(I) ##
    elseif uppercase(Model.scoreScalingType) == "FISHER-DIAG"
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
    end
    return scalingMat
end


function predict_score_driven_par( Model::GasNetModelDirBinGlobalPseudo, obs_t::Array{<:Real,2},
                         ftot_t::Array{<:Real,1}, I_tm1::Array{<:Real,2},
                         indTvPar::BitArray{1}, Wgas::Array{<:Real,1},
                         Bgas::Array{<:Real,1}, Agas::Array{<:Real,1};
                         matrixScaling=false)
     #= likelihood and gradients depend on all the parameters (ftot_t), but
     only the time vaying ones (f_t) are to be updated=#
     NergmPar = length(obs_t[1,:])-2
     NtvPar = sum(indTvPar)
     f_t = ftot_t[indTvPar] #Time varying ergm parameters

     indPres = (obs_t[:,1]).>0 # change stats of matrix elements that are present
     grad_t = zeros(Real,NtvPar)
     Idiag_t = ones(Real,NtvPar)
     I_t = ones(Real,NtvPar,NtvPar)

     #for each ergm parameter
     tmpMatPar = ftot_t' .* obs_t[:,2:end-1]
     p_ij =  exp.( .- sum(tmpMatPar,dims = 2))
     mult =obs_t[:,end]
     logpseudolike_t =  .-  sum(mult .* log.(1 .+ 1 ./ p_ij ))
     π = (1 ./ (1 .+ p_ij  ))
     #for the static Parameters
     for i =1:sum(.!indTvPar)
         p = findall(.!indTvPar)[i]
         δ_p = obs_t[:,p+1]
         tmp1 = sum( mult[indPres] .* δ_p[indPres] )
         #@show(tmp1)
         logpseudolike_t += ftot_t[p]*tmp1
     end

     # for all the time varying parameters
     for i =1: NtvPar
         p = findall(indTvPar)[i]
         δ_p = obs_t[:,p+1]
         tmp1 = sum( mult[indPres] .* δ_p[indPres] )
         #@show(tmp1)
         logpseudolike_t += ftot_t[p]*tmp1

         #println( (1 ./ (1 + exp.(sum(tmpMatPar,2))))  )
            # display(tmp1)
         grad_t[i] = tmp1  .-  sum( mult.*( δ_p.* π))
         tmp2 = ((δ_p.^2) .*  π .*  (1 .- π) )
         #println(size(tmp2))
         Idiag_t[i] =  (  sum( mult.*tmp2))
         I_t[i,i] = Idiag_t[i]
         for j=i+1: NtvPar
             pj = findall(indTvPar)[j]
             δ_pj = obs_t[:,pj+1]
             tmp3 =((δ_p.*δ_pj) .*  π .*  (1 .- π) )
             #println(size(tmp2))
             I_t[i,j] =  (  sum( mult.*tmp3))
            # @show(I_t)
             I_t[j,i] = I_t[i,j]
         end
     end
#     I_t = Idiag_t
     if sum(indTvPar)>=1
         λ = 0.2
         #I_t = λ* I_t +( 1-λ).*eye(Real,sum(indTvPar)) #
         I_t = λ* I_tm1  .+  (1-λ) * I_t #Moving average
        # println(size(I_t))
        ##@show(grad_t)
        #@show(I_t)
        if matrixScaling
            s_t =  sqrt(I_t)\ grad_t  # identify(Model,grad_t )
        else
            s_t = ( grad_t) ./ sqrt.(Idiag_t) # identify(Model,grad_t )
        end
     else
        s_t = grad_t
     end
     #s_t = ( grad_t ./ sqrt.(Idiag_t) )
     # updating direction without rescaling (for the moment)
     #I_t = scalingMatGas(Model,ones(2,2),I_tm1)
     #display(diag(I_t))
     # display(size(gradIO_t))
     # display(size(exp_mat_t ))
     #@show(grad_t)
     #println((grad_t,scal_t))
     #println(scal_t)
     # GAS update for next iteration
     #@show(s_t)
     f_tp1 = Wgas .+ Bgas.* f_t .+ Agas.*s_t
     ftot_tp1 = copy(ftot_t)
     ftot_tp1[indTvPar] = f_tp1 #of all parameters udate dynamic ones with GAS
     return ftot_tp1, logpseudolike_t, I_t, grad_t
  
end


function score_driven_filter( Model::GasNetModelDirBinGlobalPseudo,
                                vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1};vConstPar ::Array{<:Real,1} = zeros(Real,2),
                                obsT:: Array{Array{Float64,2},1}=Model.obsT ,ftot_0::Array{<:Real,1} = zeros(Real,2))
    """GAS Filter the Dynamic Fitnesses from the Observed degrees, given the GAS parameters
     given T observations for the degrees in TxN vector degsT
     """

     #Per i modelli con pseudolikelihood this funciton allows only filtering
    NergmPar = length(obsT[1][1,:])-2#
    NTvPar   = sum(indTvPar)
    T= length(obsT);
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

    sum(ftot_0)==0  ?    ftot_0 = UMallPar : ()# identify(Model,UMallNodesIO)

    ftot_t = copy(ftot_0)
    if NTvPar==0
        I_tm1 = ones(1,1)
    else
        I_tm1 = Float64.(Diagonal{Real}(I,NTvPar))
    end
    loglike = 0

    for t=1:T
    #    println(t)
        obs_t = obsT[t] # vector of in and out degrees
        #print((t,I_tm1))
        ftot_t,loglike_t,I_tm1,~ = predict_score_driven_par(Model,obs_t,ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
        fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
        loglike += loglike_t
    end
    return fVecT, loglike
end


function gasScoreSeries(Model::GasNetModelDirBinGlobalPseudo,
                                 vTvPar::Array{<:Real,2};obsT:: Array{Array{Float64,2},1}=Model.obsT)
     # se come input viene data una serie temporale di parametr, allora usa quella
     NergmPar = length(obsT[1][1,:])-2#
     T= length(obsT);

     I_tm1 = ones(NergmPar,NergmPar)
     loglike = 0
     scoreVec_T = ones(Real,NergmPar,T)
     rescScoreVec_T = ones(Real,NergmPar,T)
     I_Store_T = ones(Real,NergmPar,T)
     indTvPar = trues(NergmPar)
     Wvec = zeros(Real,NergmPar)
     Bvec = zeros(Real,NergmPar)
     Avec = zeros(Real,NergmPar)
     matScal = true #ci sono dei problemi con il matrix scaling. quando lo uso
     # per calcolare
     for t=1:T
         #    println(t)
         obs_t = obsT[t] #
         ~,~,I_tm1,grad_t = predict_score_driven_par(Model,obs_t,vTvPar[:,t],I_tm1,indTvPar,
                                             Wvec,Bvec,Avec;matrixScaling = matScal)
         scoreVec_T[:,t] = grad_t #store the filtered parameters from previous iteration
         #I_Store_T[:,t] = I_tm1
         if matScal
             rescScoreVec_T[:,t] = I_tm1\grad_t
         else
             rescScoreVec_T[:,t] = grad_t ./ diag(I_tm1)
         end

     end

     return scoreVec_T,rescScoreVec_T,I_Store_T
end


gasScoreSeries(Model::GasNetModelDirBinGlobalPseudo, vStaticPar::Array{<:Real,1};obsT:: Array{Array{Float64,2},1}=Model.obsT) = gasScoreSeries(Model,repeat(vStaticPar,outer = (1,length(obsT)) );obsT= obsT)


# Estimation
function setOptionsOptim(Model::GasNetModelDirBinGlobalPseudo)
    "Set the options for the optimization required in the estimation of the model.
    For the optimization use the Optim package."
    tol = eps()*100
    maxIter = 5000
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


function estimate(Model::GasNetModelDirBinGlobalPseudo;indTvPar::BitArray{1}=trues(length(Model.obsT[1][1,:])-2), indTargPar::BitArray{1} = indTvPar, UM :: Array{<:Real,1} = zeros(length(Model.obsT[1][1,:])-2),ftot_0 :: Array{<:Real,1} = zeros(length(Model.obsT[1][1,:])-2),changeStats_T = Model.obsT, hess_opt_flag = false)


    # se targeting è vero allora usa le start Um com target

    #ftot_0 = zeros(length(Model.obsT[1][1,:])-2)
    #changeStats_T = Model.obsT
    #check that targeting is required only for time varying Parameters
    sum(indTargPar[.!indTvPar])!= 0 ? error() : ()

    T = length(Model.obsT);
    NergmPar = length(Model.obsT[1][1,:])-2 #
    NTvPar = sum(indTvPar)
    NTargPar = sum(indTargPar)

    # UM is a vector with target values for dynamical ones. Parameters
    # if not given as input use the static estimates
    if prod(UM.== 0 )&(!prod(.!indTvPar))
        tmpParStat,~ = estimate(Model;indTvPar = BitArray(undef,NergmPar))
        parStat = zeros(NergmPar); for i=1:NergmPar parStat[i] = tmpParStat[i][1]; end
        UM = parStat
        @show(UM)
    end

    # ftot_0 is a vector with initial values (to be used in the SD iteration)
    # if not given as input estimate on first 3 observations
    if prod(ftot_0.== 0 )&(!prod(.!indTvPar))
        tmpParStat,~ = estimate(Model;indTvPar = BitArray(undef,NergmPar),changeStats_T = changeStats_T[1:5 ])
        parStat = zeros(NergmPar); for i=1:NergmPar parStat[i] = tmpParStat[i][1]; end
        ftot_0 = parStat
        @show(ftot_0)
    end
    #UM = ftot_0

    optims_opt, algo = setOptionsOptim(Model)

    # #set the starting points for the optimizations
    B0_Re  = 0.9; B0_Un = log(B0_Re ./ (1 .- B0_Re ))
    ARe_min =0.001
    A0_Re  = 0.005 ; A0_Un = log(A0_Re  .-  ARe_min)
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

        oneInADterms  = (maxLargeVal + vecUnPar[1])/maxLargeVal
        foo,loglikelValue = score_driven_filter(Model,vecReGasParAll,indTvPar;
                                        obsT = changeStats_T,vConstPar =  vecConstPar,ftot_0 = ftot_0 .* oneInADterms)
        #println(vecReGasPar)
         return - loglikelValue
    end
    #Run the optimization
    if uppercase(Model.scoreScalingType) == "FISHER-EWMA"
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
            foo,loglikelValue = score_driven_filter(Model,vecReGasParAll,indTvPar;
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


function pValStatic_SDERGM(model::GasNetModelDirBinGlobalPseudo;
                            obsT:: Array{Array{Float64,2},1}=model.obsT,
                            estPar = zeros(2,2))
    Nterms =  length(obsT[1][1,:])-2#
    T=length(obsT)
    scoreAutoc_Info = zeros(Nterms,4)
     sum(estPar) == 0   ?     twoSteps=false : twoSteps = true
    if !twoSteps
        estParStatic,~,~, ftot_0= estimate(model;UM = zeros(Nterms),indTvPar = falses(Nterms))
        constParVec2 = zeros(Nterms); for i=1:Nterms constParVec2[i] = estParStatic[i][1]; end
        #constParVec2[2] = 0.55
        grad_T,s_T,I_T = gasScoreSeries(model,constParVec2)
        estPar = constParVec2
    else
         grad_T,s_T,I_T = gasScoreSeries(model,estPar)
    end
    pValVec = zeros(Nterms)
    for i=1:Nterms
        g_T_t = Float64.(grad_T[i,2:end])
        s_T_tm1 = Float64.(s_T[i,1:end-1])
        g_T_tm1 = Float64.(grad_T[i,1:end-1])
        g_delta_t = Float64.(grad_T[(1:Nterms).!=i,2:end])
        I_T_t = Float64.(I_T[i,2:end])
        y = ones(T-1)
        if twoSteps
            ols = lm([   g_T_t  g_T_tm1.*g_T_t],y)
        else
            ols = lm([  g_delta_t' g_T_t  s_T_tm1.*g_T_t],y)
        end
        yHat = predict(ols)
        yMean = 1
        RSS = sum((yHat .- yMean).^2)
        ESS = T - RSS
        pval =  ccdf(Chisq(1), ESS)
        pValVec[i] = pval
        scoreAutoc_Info[i,1] = autocor(Float64.(grad_T[i,:]),[1])[1]
        scoreAutoc_Info[i,2] = autocor(Float64.(s_T[i,:]),[1])[1]
        scoreAutoc_Info[i,3] = cor(s_T_tm1, g_T_t)
        scoreAutoc_Info[i,4] = pval
    end
    return pValVec,scoreAutoc_Info,estPar
end


function dgp_paper_SDERGM(model::GasNetModelDirBinGlobalPseudo;dgpType = "sin", T = 50, N = 50,
                            Nterms = 2, Nsteps1 = 1 ,Nsteps2 = 0)
     #Generate the matrix of parameters values for each t, i.e. the DGP path for
     #the TV parametrs
     parDgpT = zeros(Nterms,T)
     minpar1 = -3#-5
     maxpar1 = -2.5# -1.5
     minpar2 = 0.05#0.02
     maxpar2 = 0.5#0.7
     if dgpType =="sin"
         parDgpT[1,:] = dgpSin(minpar1,maxpar1,Nsteps1,T)# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
         parDgpT[2,:] = dgpSin(minpar2,maxpar2,Nsteps2,T) #

     elseif dgpType=="steps"
         parDgpT[1,:] = randSteps(minpar1,maxpar1,Nsteps1,T)# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
         parDgpT[2,:] = randSteps(minpar2,maxpar2,Nsteps2,T)#
         #parDgpT[2,:] = 0.25
     elseif dgpType=="AR"
         B = 0.95
         sigma = 0.1
         parDgpT[1,:] = dgpAR(minpar1,B,sigma,T,minMax=[minpar1,maxpar1])
         Nsteps1==0 ? parDgpT[1,:] =minpar1 : ()# -3# randSteps(0.05,0.5,2,T) #1.5#.00000000000000001
         parDgpT[2,:] = dgpAR((maxpar2 + minpar2)/2,B,sigma,T;minMax = [minpar2,maxpar2])#
         nonoInds = parDgpT[2,:].<0
         if sum(nonoInds)>0
             parDgpT[2,:] = parDgpT[2,:] - minimum(parDgpT[2,nonoInds])
         end
         Nsteps2==0 ? parDgpT[2,:] =maxpar2 : ()
          load_fold = "./data/estimatesTest/sdergmTest/R_MCMC_estimates/"
         @load(load_fold * "ARpath_edges_and_GWESP.jld",parDgpT_store)
         parDgpT = parDgpT_store'
         parInd = 1
         # close()
         # subplot(1,2,1);plot(1:T,parDgpT[parInd,:],"k",linewidth=5)
         # parInd = 2
         # subplot(1,2,2);plot(1:T,parDgpT[parInd,:],"k",linewidth=5)
     end

     return parDgpT
end


function logLike_t(Model::GasNetModelDirBinGlobalPseudo, obsT,
                  vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1};
                  ftot_0::Array{<:Real,1} = zeros(Real,2))
      """
            loglikelihood of the last observartion, depends on the static
            parameters beacause of all previous iterations
      """
      NergmPar = length(obsT[1][1,:])-2#
      NTvPar   = sum(indTvPar)
      T= length(obsT);
      # Organize parameters of the GAS update equation
      NTvPar = sum(indTvPar)
      Wvec = zeros(Real,NTvPar)
      Bvec = zeros(Real,NTvPar)
      Avec = zeros(Real,NTvPar)
      vConstPar = zeros(Real,NTvPar)
      lastTv_i = 1
      lastConst_i = 1
      lastTot_i = 1
      for p in 1:NergmPar
         if indTvPar[p]
             Wvec[lastTv_i] = vResGasPar[lastTot_i]
             Bvec[lastTv_i] = vResGasPar[lastTot_i+1]
             Avec[lastTv_i] = vResGasPar[lastTot_i+2]
             lastTv_i = lastTv_i + 1
             lastTot_i = lastTot_i + 3
         else
             vConstPar[lastConst_i] = vResGasPar[lastTot_i]
             lastConst_i = lastConst_i + 1
             lastTot_i = lastTot_i + 1
         end
      end
      # start values equal the unconditional mean,and  constant ones remain equal to the unconditional mean, hence initialize as:
      UMTVPar =  Wvec ./ (1 .- Bvec)
      fVecT = ones(Real,NergmPar,T)
      if sum(ftot_0)==0
          ftot_0 =  zeros(Real,NergmPar)
          ftot_0[indTvPar] = UMTVPar
          ftot_0[indTvPar] = vConstPar
      end

      ftot_t = copy(ftot_0)
      if NTvPar==0
            I_tm1 = ones(1,1)
      else
            I_tm1 = Float64.(Diagonal{Real}(I,NTvPar))
      end
      loglike_t = 0
      for t=1:T
            obs_t = obsT[t] # vector of in and out degrees
            ftot_t, loglike_t, I_tm1, ~ = predict_score_driven_par(Model, obs_t, ftot_t,
                                                        I_tm1, indTvPar, Wvec,
                                                        Bvec, Avec)
            fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
      end
      return  loglike_t::T where T <:Real
end


function logLike_T(Model::GasNetModelDirBinGlobalPseudo, obsT,
                  vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1};
                  ftot_0::Array{<:Real,1} = zeros(Real,2))
      NergmPar = length(obsT[1][1,:])-2#
      T= length(obsT);
      NTvPar = sum(indTvPar)
      Wvec = zeros(Real,NTvPar)
      Bvec = zeros(Real,NTvPar)
      Avec = zeros(Real,NTvPar)
      vConstPar = zeros(Real,NTvPar)
      lastTv_i = 1
      lastConst_i = 1
      lastTot_i = 1
      for p in 1:NergmPar
         if indTvPar[p]
             Wvec[lastTv_i] = vResGasPar[lastTot_i]
             Bvec[lastTv_i] = vResGasPar[lastTot_i+1]
             Avec[lastTv_i] = vResGasPar[lastTot_i+2]
             lastTv_i = lastTv_i + 1
             lastTot_i = lastTot_i + 3
         else
             vConstPar[lastConst_i] = vResGasPar[lastTot_i]
             lastConst_i = lastConst_i + 1
             lastTot_i = lastTot_i + 1
         end
      end

      # start values equal the unconditional mean,and  constant ones remain equal to the unconditional mean, hence initialize as:
      UMallPar = zeros(Real,NergmPar)
      UMallPar[indTvPar] =  Wvec ./ (1 .- Bvec)
      if !prod(indTvPar) # if not all parameters are time varying
            UMallPar[.!indTvPar] = vConstPar
      end
      fVecT = ones(Real,NergmPar,T)
      sum(ftot_0)==0  ?    ftot_0 = UMallPar : ()# identify(Model,UMallNodesIO)
      ftot_t = copy(ftot_0)
      if NTvPar==0
            I_tm1 = ones(1,1)
      else
            I_tm1 = Float64.(Diagonal{Real}(I,NTvPar))
      end
      loglike_T = zero(Real)
      for t=1:T
            obs_t = obsT[t] # vector of in and out degrees
            ftot_t,loglike_t,I_tm1,~ = predict_score_driven_par(Model,obs_t,ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
            fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
            loglike_T += loglike_t
      end

      return  loglike_T::T where T <:Real
end
