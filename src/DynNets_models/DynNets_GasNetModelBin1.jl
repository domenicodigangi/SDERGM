#------------UNDIRECTED NETWORKS ----------------------------------------------------------

struct  GasNetModelBin1 <: GasNetModelBin
    """
    A Logistic Gas model for binary networks and probability depending only
            on time varying parameters.
            p^{(t)}_{ij}  = logistic( theta^{(t)}_i + theta^{(t)}_j)
            # Bin1 Model with one unconditional mean per node and one gas parameter per group of nodes
            # Model Definition and parameters conversion functions
            """
    # node - specific W and group Specific (A B)
    obsT::Array{<:Real,2}
    Par::Array{Array{<:Real,1},1} #Both nodes and Group specific parameters (Gas and non Gas)
    # a group of nodes has the same GAS parameters, but different fitnesses
    groupsInds::Array{Array{<:Real,1},1} #I assume that nodes with the same Wgroup have the same ABgroup!!!
    scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
    # choices see function scalingMatGas
 end

GasNetModelBin1(obsT::Array{<:Real,2},Par::Array{Array{<:Real,1},1},groupsInds::Array{Array{<:Real,1},1}) =
            GasNetModelBin1(obsT,Par,groupsInds,"")
GasNetModelBin1(obsT::Array{<:Real,2}) = GasNetModelBin1(obsT,
                                                            [ zeros(length(obsT[1,:])),
                                                             0.9  * ones(1),
                                                             0.01 * ones(1)],
                                                             [ones(Int64,length(obsT[1,:])) ,ones(1)],"" )

GasNetModelBin1(obsT::Array{<:Real,2},scoreScalingType::String) = GasNetModelBin1(obsT,
                                                            [ zeros(length(obsT[1,:])),
                                                             0.9  * ones(1),
                                                             0.01 * ones(1)],
                                                             [ones(Int64,length(obsT[1,:])) ,ones(1)],scoreScalingType)

fooGasNetModelBin1 = GasNetModelBin1(ones(3,3))
# Functions that are obtained directly from the static version of the model

StaModType(Model::GasNetModelBin1) = StaticNets.fooNetModelBin1# to be substituted with a conversion mechanism
# StaPar2DynPar(Model::GasNetModelBin1,RePar::Array{<:Real,1}) = log.(RePar)
# DynPar2StaPar(Model::GasNetModelBin1,UnPar::Array{<:Real,1}) = exp.(UnPar)
linSpacedPar(Model::GasNetModelBin1,Nnodes::Int;NgroupsW = Nnodes, deltaN::Int=3,graphConstr = true) =   StaticNets.linSpacedPar(StaModType(Model),Nnodes;Ngroups = NgroupsW,deltaN=deltaN,graphConstr =graphConstr);


# options and conversions of parameters for optimization
function setOptionsOptim(Model::GasNetModelBin1)
    "Set the options for the optimization required in the estimation of the model.
    For the optimization use the Optim package."
    tol = eps()*1000
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
      return opt, algo
 end
function array2VecGasPar(Model::GasNetModelBin1,ArrayGasPar::Array{Array{Real,1},1})
         VecGasPar = [ArrayGasPar[1];ArrayGasPar[2];ArrayGasPar[3];ArrayGasPar[4]]
         return VecGasPar
     end
function vec2ArrayGasPar(Model::GasNetModelBin1,VecGasPar::Array{<:Real,1})
     T,N2 = size(Model.obsT)
     groupsInds = Model.groupsInds
     NGW = length(unique(groupsInds[1]))
     GBA = length( unique(groupsInds[2][.!(groupsInds[2].==0)]) )
     ArrayGasPar =  Array{Array{Real,1},1}(3)
     ArrayGasPar[1] = VecGasPar[1:NGW]
     ArrayGasPar[2] = ones(GBA).*VecGasPar[NGW+1:NGW+GBA]
     ArrayGasPar[3] = ones(GBA).*VecGasPar[NGW+GBA+1:NGW+2GBA]
     return ArrayGasPar
     end
function restrictGasPar(Model::GasNetModelBin1,vecUnGasPar::Array{<:Real,1})
     "From the Unrestricted values of the parameters return the restricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, e.g. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     groupsInds = Model.groupsInds
     NGW = length(unique(groupsInds[1]))
     GBA = length( unique(groupsInds[2][.!(groupsInds[2].==0)]) )
     W_un = vecUnGasPar[1:NGW]
     diag_B_Un = vecUnGasPar[NGW+1:NGW+GBA]
     diag_A_Un = vecUnGasPar[NGW+ GBA+1:NGW+2GBA]
     W_Re = W_un
     diag_A_Re = exp.(diag_A_Un)
     diag_B_Re = 1 ./ (1+exp.(.-diag_B_Un))
     vecRePar =  [W_Re; diag_B_Re; diag_A_Re]
     return vecRePar
     end
function unRestrictGasPar( Model::GasNetModelBin1,vecReGasPar::Array{<:Real,1})
     "From the restricted values of the parameters return the unrestricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, i.e. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     groupsInds = Model.groupsInds
     NGW = length(unique(groupsInds[1]))
     GBA = length( unique(groupsInds[2][.!(groupsInds[2].==0)]) )

     W_Re = vecReGasPar[1:NGW]
     diag_B_Re = vecReGasPar[NGW+1:NGW+GBA]
     diag_A_Re = vecReGasPar[NGW+GBA+1:end]

     W_Un = W_Re
     diag_A_Un = log.(diag_A_Re)
     diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
     vecUnPar =  [W_Re;diag_B_Un; diag_A_Un]
     return vecUnPar
     end

#Gas Filter Functions

function scalingMatGas(Model::GasNetModelBin1,expMat::Array{<:Real,2},I_tm1::Array{<:Real,2})
    "Return the matrix required for the scaling of the score, given the expected
     matrix and the Scaling matrix at previous time. "
    if uppercase(Model.scoreScalingType) == ""
        scalingMat = 1 #
    elseif uppercase(Model.scoreScalingType) == "FISHER-EWMA"
        λ = 0.8

        I = expMat.*(1 .- expMat)
        diagI = sum(I,dims = 2)
        I = zeros(expMat)
        [I[i,i] = diagI[i] for i=1:length(diagI) ]
        I_t =  λ*I + (1-λ) *I_tm1
        scalingMat = sqrt(I) ##
    end
    return scalingMat
 end
function predict_score_driven_par( Model::GasNetModelBin1,N::Int,degs_t::Array{<:Real,1},
                        ftot_t::Array{<:Real,1},I_tm1::Array{<:Real,2},indTvNodes::BitArray{1},
                        Wgas::Array{<:Real,1},Bgas::Array{<:Real,1},Agas::Array{<:Real,1})
    "Return the GAS updated parameters, with the scaling matrix and the
    loglikelihood of the current observation."
    #= likelihood and gradients depend on all the parameters (ftot_t), but
    only the time vaying ones (f_t) are to be updated=#
    f_t = ftot_t[indTvNodes] #Time varying fitnesses

    thetas_mat_t_exp, exp_mat_t = StaticNets.expMatrix2(StaticNets.fooNetModelBin1,ftot_t)
    # thetas_mat_t_exp = exp.(Symmetric(ftot_t  .+  ftot_t'))
    # StaticNets.putZeroDiag!(thetas_mat_t_exp)
    # exp_mat_t = thetas_mat_t_exp ./ (1+thetas_mat_t_exp)

    exp_deg_t = sum(exp_mat_t,dims = 2)
    #The following part is a faster version of logLikelihood_t()

    loglike_t = sum(ftot_t.*degs_t) -  sum(UpperTriangular(log.(1 + thetas_mat_t_exp))) #  sum(log.(1 + thetas_mat_t_exp))
    # compute the score of observations at time t wrt parameters at tm1
    grad_t = degs_t[indTvNodes]  .-  exp_deg_t[indTvNodes] # only for TV paramaters
    I_t = scalingMatGas(Model,exp_mat_t,I_tm1)
    s_t = I_t\grad_t
    # GAS update for next iteration
    f_tp1 = Wgas + Bgas.* f_t + Agas.*s_t
    f_t = f_tp1

    ftot_t[indTvNodes] = f_t #of all parameters udate dynamic ones with GAS
    ftot_tp1 = ftot_t
    return ftot_tp1,I_t,loglike_t,s_t
 end
function score_driven_filter( Model::GasNetModelBin1,
                                vResGasPar::Array{<:Real,1};
                                obsT::Array{<:Real,2}=Model.obsT,
                                groupsInds:: Array{Array{<:Real,1},1}=Model.groupsInds,
                                dgpTN = (0,0))
    "Run the GAS Filter  for the input model given the GAS Parameters and
     observations, and return the TV parameters with the loglikelihood.
     Alternatively, if dgpTN != (0,0) then use it as a Data
     Generating Process and return the TV parameters and the Observations."

     #se parametri per dgp sono dati allora campiona il dgp, altrimenti filtra
    sum(dgpTN) == 0   ?     dgp = false : dgp = true
    if dgp
        T,N = dgpTN
     Y_T = zeros(Int8,T,N,N)
    else
     T,N = size(obsT)
    end

    NGW = length(unique(groupsInds[1]))# Number of W groups
    GBA = length( unique(groupsInds[2][.!(groupsInds[2].==0)]) )# Number of BA groups
    ABgroupsIndNodes = groupsInds[2][groupsInds[1]]#
    indTvNodes =   .!(ABgroupsIndNodes[1:N] .==0)
    # Organize parameters of the GAS update equation
    WGroups = vResGasPar[1:NGW] #W parameters for each group
    W_all = WGroups[groupsInds[1]]# W par for each node
    BgasGroups  = vResGasPar[NGW+1:NGW+GBA]# B par for each group

    AgasGroups  = vResGasPar[NGW+GBA+1:NGW+2GBA]# A par for each group
    #distribute nodes among groups with equal parameters
    Bgas = BgasGroups[ABgroupsIndNodes[indTvNodes]]#assign 1 B par to each TV node
    Agas = AgasGroups[ABgroupsIndNodes[indTvNodes]]#assign 1 A par to each TV node
    Wgas   = W_all[indTvNodes] #  The W par of TV nodes
     # start values equal the unconditional mean, but  constant ones remain equal to the unconditional mean, hence initialize as:
    UMallNodes = W_all
    UMallNodes[indTvNodes] =  Wgas ./ (1 .- Bgas)

    fVecT = ones(Real,T,N)
    ftot_t = UMallNodes
    I_tm1 = UniformScaling(N) # Intialization of the scaling Mat
    loglike = 0
    for t=1:T
        if dgp
            expMat_t = StaticNets.expMatrix(StaModType(Model),ftot_t )# exp.(ftot_t)
            Y_T[t,:,:] =StaticNets.samplSingMatCan(StaModType(Model),expMat_t)
            degs_t = dropdims(sum(Y_T[t,:,:],2),dims = 2)
        else
            degs_t = obsT[t,:] # vector of in and out degrees
        end
        ftot_t,I_tm1,loglike_t = predict_score_driven_par(Model,N,degs_t,ftot_t,I_tm1,indTvNodes,Wgas,Bgas,Agas)
        fVecT[t,:] = ftot_t #store the filtered parameters from previous iteration
        loglike += loglike_t
    end
    if dgp
        return Y_T , fVecT
    else
        return fVecT, loglike
    end
    end
sampl(Mod::GasNetModelBin1,T::Int) = (  N = length(Mod.groupsInds[1]) ;
                                        tmpTuple = score_driven_filter(Mod,[Mod.Par[1];Mod.Par[2];Mod.Par[3]];  dgpTN = (T,N) );
                                        degs_T = dropdims(sum(tmpTuple[1],dims = 2),2);
                                        (GasNetModelBin1(degs_T,Mod.Par,Mod.groupsInds,Mod.scoreScalingType),tmpTuple[1],tmpTuple[2]) )

# Estimation
function estSingSnap(Model::GasNetModelBin1,degs_t::Array{<:Real,1}; groupsInds = Model.groupsInds, targetErr::Real=targetErrValDynNets)
    hatUnPar,~,~ = StaticNets.estimate(StaModType(Model); deg = degs_t ,groupsInds = groupsInds[1], targetErr =  targetErr)
    return hatUnPar
 end
function estimateSnapSeq(Model::GasNetModelBin1; degsT::Array{<:Real,2}=Model.obsT,
                     groupsInds = Model.groupsInds,targetErr::Real=targetErrValDynNets)
    T,N = size(degsT)
    NGW = length(unique(groupsInds[1]))
    UnParT = zeros(T,NGW)

    for t = 1:T
        #println(t)
        degs_t = degsT[t,:]
        UnParT[t,:] = estSingSnap(Model,degs_t; groupsInds = groupsInds ,targetErr =  targetErr)
    end
    #Remove infinites
    UnParT[1,findall(.!isfinite.(UnParT[1,:]))] = 0
    for t=2:T
        infInd = findall(.!isfinite.(UnParT[t,:]))
        UnParT[t,infInd] = UnParT[t-1,infInd]
    end
    return UnParT
 end
function estimate(Model::GasNetModelBin1;start_values = [zeros(Real,10),zeros(Real,3),zeros(Real,3) ])
    T,N = size(Model.obsT)
    groupsInds = Model.groupsInds
    NGW = length(unique(groupsInds[1]))
    GBA = length( unique(groupsInds[2][.!(groupsInds[2].==0)]) )
    ABgroupsIndNodes = groupsInds[2][groupsInds[1]]
    indTvNodes =   .!(ABgroupsIndNodes[1:N] .==0)
    optims_opt, algo = setOptionsOptim(Model)
    #set the starting points for the optimizations if not given
    if sum([sum(start_values[i]) for i =1:3 ])==0
         B0_ReGroups = 0.9 .* ones(GBA)
         A0_ReGroups = 0.01 .* ones(GBA)
         BgasNodes_0   = B0_ReGroups[ ABgroupsIndNodes[indTvNodes]]
         AgasNodes_0   = A0_ReGroups[ ABgroupsIndNodes[indTvNodes]]
         # choose the starting points for W in order to fix the Unconditional Means
         # equal to the static estimates for the average degrees(over time)
         MeanDegsT =  dropdims(mean(Model.obsT,1),dims= 1)
         UMstaticGroups = estSingSnap(Model ,MeanDegsT; groupsInds = groupsInds )
         conv_flag1 = true

         GBA>0  ?    ( W0_Groups = UMstaticGroups.*(1 .- B0_ReGroups[groupsInds[2]])) : W0_Groups = UMstaticGroups
    else
         W0_Groups = start_values[1]
         B0_ReGroups = start_values[2]
         A0_ReGroups = start_values[3]
         conv_flag1 = true
    end
    if GBA == 0  return [UMstaticGroups, zeros(1), zeros(1)] ,  conv_flag1 end;
    # Estimate  jointly all  Gas parameters
    vGasParAll0_Un = unRestrictGasPar(Model,[W0_Groups ;B0_ReGroups;A0_ReGroups ] ) # vGasParGroups0_Un#
    # objective function for the optimization
    function objfunGas(vparUn::Array{<:Real,1})# a function of the groups parameters
              foo,loglikelValue = score_driven_filter( Model,restrictGasPar(Model,vparUn))
           return - loglikelValue
    end

    #Run the optimization without AD if a matrix sqrt is required in the filter
    if uppercase(Model.scoreScalingType) == "FISHER-EWMA"
        ADobjfunGas = objfunGas
    else
        ADobjfunGas = TwiceDifferentiable(objfunGas, vGasParAll0_Un; autodiff = :forward);
    end

    optim_out2  = optimize(ADobjfunGas,vGasParAll0_Un,algo,optims_opt)
    outParAll = Optim.minimizer(optim_out2)

    vGasParGroupsHat_Re = restrictGasPar(Model,outParAll)
    arrayGasParHat_Re = vec2ArrayGasPar(Model,vGasParGroupsHat_Re )
    conv_flag = conv_flag1*Optim.converged(optim_out2)
    #println(optim_out2)
    return  arrayGasParHat_Re, conv_flag
    end
function estimateTargeting(Model::GasNetModelBin1)
    "Estimate the GAS parameters of the GAS dynamic Undirected fitness model,
    using a targeting on the unconditional means. Then return the GASparameters "
    T,N = size(Model.obsT)
    groupsInds = Model.groupsInds
    NGW = length(unique(groupsInds[1]))
    GBA = length( unique(groupsInds[2][.!(groupsInds[2].==0)]) )
    ABgroupsIndNodes = groupsInds[2][groupsInds[1]]
    indTvNodes =   .!(ABgroupsIndNodes[1:N] .==0)

    optims_opt, algo = setOptionsOptim(Model)
    #target the unconditional mean
    ParTvSnap = estimateSnapSeq(Model) # sequence of snapshots estimate
    uncMeansWGroups = dropdims(mean(ParTvSnap,1),dims = 1)
    #functions needed for the specific targeting case
    function WGroupsFromB(B::Array{<:Real,1})
         BGroupsW = B[groupsInds[2]]
         WGroups = uncMeansWGroups.*(1 .- BGroupsW)
         return WGroups
    end
    function restrAB(vecUnPar::Array{<:Real,1}) #restrict a vector of only A and B
         diag_B_Un = vecUnPar[1:GBA]
         diag_A_Un = vecUnPar[GBA+1:end]

         diag_A_Re = exp.(diag_A_Un)
         diag_B_Re = 1 ./ (1+exp.(.-diag_B_Un))
         vecRePar =  [diag_B_Re; diag_A_Re]
         return vecRePar
         end
    function UnrestrAB(vecRePar::Array{<:Real,1}) #restrict a vector of only A and B
         diag_B_Re = vecRePar[1:GBA]
         diag_A_Re = vecRePar[GBA+1:end]
         diag_A_Un = log.(diag_A_Re)
         diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
         vecUnPar =  [diag_B_Un; diag_A_Un]
         return vecUnPar
    end
    if GBA == 0  return [uncMeans, zeros(1), zeros(1)] , Optim.converged(optim_out) end;
    # Estimate As and Bs


    # #set the starting points for the optimizations
    B0_ReGroups = 0.99 .* ones(GBA)
    A0_ReGroups = 0.07 .* ones(GBA)
    vGasParAll0_Un = UnrestrAB([B0_ReGroups;A0_ReGroups ] ) # vGasParGroups0_Un#
    # objective function for the optimization
    function objfunGas(vecUnAB::Array{<:Real,1})# a function of the groups parameters
         vecReAB = restrAB(vecUnAB)
         ReB = vecReAB[1:GBA]
         WGroups = WGroupsFromB(ReB)
         vecReGasPar = [WGroups;vecReAB]
         foo,loglikelValue = score_driven_filter( Model,vecReGasPar )
         return - loglikelValue
    end
    #Run the optimization
    if uppercase(Model.scoreScalingType) == "FISHER-EWMA"
        ADobjfunGas = objfunGas
    else
        ADobjfunGas = TwiceDifferentiable(objfunGas, vGasParAll0_Un; autodiff = :forward);
    end

    optim_out2  = optimize(ADobjfunGas,vGasParAll0_Un,algo,optims_opt)
    outParAB = Optim.minimizer(optim_out2)
    outReAB = restrAB(outParAB)
    vGasParGroupsHat_Re = [WGroupsFromB(outReAB[1:GBA]);outReAB]
    arrayGasParHat_Re = vec2ArrayGasPar(Model,vGasParGroupsHat_Re )
    conv_flag =  Optim.converged(optim_out2)
    # println(optim_out2)
    return  arrayGasParHat_Re, conv_flag
    end
