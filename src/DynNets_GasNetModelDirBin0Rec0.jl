

struct  GasNetModelDirBin0Rec0 <: GasNetModelBin
     """ A gas model based on pseudolikelihood (as objective function) for
            Directed binary networks and probability depending on a generic vector
            of global statistics each associated with a time varying parameters.
         """
     obsT:: Array{Array{Float64,1},1}# for each t we have 2 statistics: L and R
     N::Int
     Par::Array{Array{<:Real,1},1} #  each time varying parameter has 3 static parameters in this specification.
                         # if a parameter is constant than only 1 gas parameteris present
     indTvPar :: BitArray{1} #  what parameters are time varying   ?
     scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
     # choices see function scalingMatGas
  end

#inizializza senza specificare assenza di scaling
GasNetModelDirBin0Rec0( obsT, N,  Par   ) =
        GasNetModelDirBin0Rec0( obsT, N, Par ,trues(length(Par[:,1])) ,"")
#Initialize by observations only
GasNetModelDirBin0Rec0(obsT::  Array{Array{Real,2}}) =   (Npar =  length(obsT[1][1,:]) - 2;
                                                GasNetModelDirBin0Rec0(obsT,zeros(Npar,3) ))

GasNetModelDirBin0Rec0(obsT:: Array{Real,3},scoreScalingType::String) =   (N = round(length(obsT[:,1])/2);
                                                GasNetModelDirBin0Rec0(obsT,zeros(3Npar),scoreScalingType) )
fooGasPar = [ones(2), ones(2), ones(2)]
fooGasNetModelDirBin0Rec0 = GasNetModelDirBin0Rec0( [zeros(2), zeros(2)], 10, fooGasPar)

# Relations between Static and Dynamical models: conventions on storage for
# parameters and observations
StaModType(Model::GasNetModelDirBin0Rec0 ) = StaNets.fooNetModelDirBin0Rec0# to be substituted with a conversion mechanism

# options and conversions of parameters for optimization
setOptionsOptim(Model::GasNetModelDirBin0Rec0) = setOptionsOptim(fooGasNetModelDirBin1)

array2VecGasPar(Model::GasNetModelDirBin0Rec0,
                ArrayGasPar, indTvPar :: BitArray{1}) =
                array2VecGasPar(fooGasNetModelDirBinGlobalPseudo,
                                ArrayGasPar, indTvPar )

vec2ArrayGasPar(Model::GasNetModelDirBin0Rec0, VecGasPar::Array{<:Real,1},
                indTvPar :: BitArray{1}) =
                vec2ArrayGasPar(fooGasNetModelDirBinGlobalPseudo, VecGasPar,
                                indTvPar)

restrictGasPar(Model::GasNetModelDirBin0Rec0, vecUnGasPar::Array{<:Real,1},
                        indTvPar :: BitArray{1}) =
                        restrictGasPar(fooGasNetModelDirBinGlobalPseudo,
                                        vecUnGasPar, indTvPar )

unRestrictGasPar(Model::GasNetModelDirBin0Rec0,
                        vecReGasPar::Array{<:Real,1},
                        indTvPar :: BitArray{1})=
                        unRestrictGasPar(fooGasNetModelDirBinGlobalPseudo ,
                                            vecReGasPar,  indTvPar )

#Gas Filter Functions
function identify(Model::GasNetModelDirBin0Rec0,parIO::Array{<:Real,1})
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

function scalingMatGas(Model::GasNetModelDirBin0Rec0,expMat::Array{<:Real,2},I_tm1::Array{<:Real,2})
    "Return the matrix required for the scaling of the score, given the expected
     matrix and the Scaling matrix at previous time. "
    if uppercase(Model.scoreScalingType) == ""
        scalingMat = 1 #
    elseif uppercase(Model.scoreScalingType) == "FISHER-EWMA"
        error()
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
        error()
    end
    return scalingMat
 end



function updatedGasPar( Model::GasNetModelDirBin0Rec0, obs_t, N,
                         ftot_t::Array{<:Real,1}, I_tm1::Array{<:Real,2},
                         indTvPar::BitArray{1}, Wgas::Array{<:Real,1},
                         Bgas::Array{<:Real,1}, Agas::Array{<:Real,1};
                         matrixScaling=false)
     #= likelihood and gradients depend on all the parameters (ftot_t), but
     only the time vaying ones (f_t) are to be updated=#
     NergmPar = 2
     NtvPar = sum(indTvPar)
     f_t = ftot_t[indTvPar] #Time varying ergm parameters

     L_t, R_t = obs_t
     logLike_t = StaNets.logLikelihood(StaNets.fooNetModelDirBin0Rec0, L_t, R_t, ftot_t)

     θ_t = f_t[1]
     η_t = f_t[2]

     den = 1 + 2*exp(θ_t) + exp(2*(θ_t + η_t))
     grad_tot_t = [ L_t - N*(N-1) *(exp(θ_t) + exp(2*(θ_t + η_t)))/den,
                    R_t - N*(N-1) *(exp(2*(θ_t + η_t)))/den]

     # No rescaling
     I_t = I_tm1

     s_t = grad_tot_t[indTvPar]

     f_tp1 = Wgas .+ Bgas.* f_t .+ Agas.*s_t
     ftot_tp1 = copy(ftot_t)
     ftot_tp1[indTvPar] = f_tp1 #of all parameters udate dynamic ones with GAS
     return ftot_tp1, logLike_t, I_t, grad_tot_t
  end


function gasFilter( Model::GasNetModelDirBin0Rec0,
                    vResGasPar::Array{<:Real,1}, indTvPar::BitArray{1};vConstPar ::Array{<:Real,1} = zeros(Real,2),
                    obsT=Model.obsT, N = Model.N, ftot_0::Array{<:Real,1} = zeros(Real,2), dgpNT = (0,0))
    """GAS Filter the Dynamic Fitnesses from the Observed degrees, given the GAS parameters
     given T observations for the degrees in TxN vector degsT
     """

     #Per i modelli con pseudolikelihood this funciton allows only filtering

    sum(dgpNT) == 0   ?     dgp = false : dgp = true
    NergmPar = 2#
    NTvPar   = sum(indTvPar)
    if dgp
        T = dgpNT[2]
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
    if dgp
        A_T = zeros(Int8, N, N, T)
    end
    for t=1:T
    #    println(t)
        if dgp
            diadProb = StaNets.diadProbFromPars(StaNets.fooNetModelDirBin0Rec0, ftot_t )
            A_t = StaNets.samplSingMatCan(StaNets.fooNetModelDirBin0Rec0, diadProb, N)
            A_T[:, : , t] = A_t
            obs_t = StaNets.statsFromMat(StaNets.fooNetModelDirBin0Rec0, A_t)
        else
            obs_t = obsT[t]
        end

        #print((t,I_tm1))
        ftot_t,loglike_t,I_tm1,~ = updatedGasPar(Model,obs_t, N, ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
        fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
        loglike += loglike_t
    end
    if dgp
        return fVecT, A_T
    else
        return fVecT, loglike
    end
    end

# Estimation
function setOptionsOptim(Model::GasNetModelDirBin0Rec0)
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

function estimate(Model::GasNetModelDirBin0Rec0; indTvPar::BitArray{1}=trues(length(Model.obsT[1][1,:])-2),
                    indTargPar::BitArray{1} = indTvPar, UM :: Array{<:Real,1} = zeros(length(Model.obsT[1][1,:])-2),
                    ftot_0 :: Array{<:Real,1} = zeros(length(Model.obsT[1][1,:])-2), obsT = Model.obsT)
    "Estimate the GAS and static parameters  "

    T = length(Model.obsT);
    NergmPar = 2 #
    NTvPar = sum(indTvPar)
    NTargPar = sum(indTargPar)

    # UM is a vector with target values for dynamical ones. Parameters
    # if not given as input use the static estimates
    # single static estimate
    L_mean  = mean([stat[1] for stat in obsT ])
    R_mean  = mean([stat[2] for stat in obsT ])
    parStat = estimate(fooNetModelDirBin0Rec0, L_mean, R_mean, N )
    if prod(UM.== 0 )&(!prod(.!indTvPar))
        UM = parStat
    end

    # ftot_0 is a vector with initial values (to be used in the SD iteration)
    # if not given as input estimate on first 3 observations
    if prod(ftot_0.== 0 )&(!prod(.!indTvPar))
        ftot_0 = parStat
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
        # costant parameters, restrict the parameters to appropriate HelperFunDomains
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
        foo,loglikelValue = gasFilter(Model,vecReGasParAll,indTvPar;
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
            foo,loglikelValue = gasFilter(Model,vecReGasParAll,indTvPar;
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
