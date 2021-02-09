#DIRECTED NETWORKS------------
# --------------------------Sender Receiver Effect
struct  GasNetModelDirBin1 <: GasNetModelBin
     """ A Logistic Gas model for Directed binary networks and probability depending only
             on time varying parameters.
             ''p^{(t)}_{ij}  = logistic( theta^{(t)}_i + theta^{(t)}_j)''
             Bin1 Model with one unconditional mean per node and one gas parameter per group of nodes
             # Model Definition and parameters conversion functions
         """
     # node - specific W and group Specific (A B)
     obsT:: Array{<:Real,2} # In degrees and out Degrees for each t . The first dimension is always time. Tx2N
     Par::Array{Array{<:Real,1},1} #[[InW;outW], Agroups , Bgroups]  Group specific parameters
     #Nodes in a group have the same GAS parameters, but different fitnesses

     # each node belongs to a W group, each W group is assigned 2 AB groups, in order
     #to allow different parameters for in and our. The AB groups In and Out can be the same

     groupsInds::Array{Array{Int,1},1} #[[Par_2_groupsW, Wgroups_2_groupsAB ]  (_2_ stands for "assignement to" )
     # In and out W parameters are assumed to be  different, while In and Out AB
     # are be equal. Moreover, I assume that nodes with the same Wgroup have the same ABgroups
     scoreScalingType::String # String that specifies the rescaling of the score. For a list of possible
     # choices see function scalingMatGas
  end

#inizializza senza specificare assenza di scaling
GasNetModelDirBin1( obsT , Par  , groupsInds) =
        GasNetModelDirBin1( obsT,Par, groupsInds,"")
#Initialize by observations only
GasNetModelDirBin1(obsT:: Array{<:Real,2}) =   (N = round(length(obsT[:,1])/2);
                                                GasNetModelDirBin1(obsT,
                                                    [ zeros(Int(2N)) ,0.9  * ones(1), 0.01 * ones(1)],
                                                    [Int.(1:2N),Int.(ones(Int(2N)))],"") )

GasNetModelDirBin1(obsT:: Array{<:Real,2}, scoreScalingType::String, group_string::String) =
    (N = round(length(obsT[:,1])/2);
    if group_string == "ONE_PAR_EACH"
        return GasNetModelDirBin1(obsT,
            [ zeros(Int(2N)) ,0.9  * ones(1), 0.01 * ones(1)],
            [Int.(1:2N), Int.(1:2N)], scoreScalingType)
    elseif group_string == "ONE_PAR_ALL"
        return GasNetModelDirBin1(obsT,
            [ zeros(Int(2N)) ,0.9  * ones(1), 0.01 * ones(1)],
            [Int.(1:2N),Int.(ones(Int(2N)))],scoreScalingType)
    end)
GasNetModelDirBin1(obsT:: Array{<:Real,2},scoreScalingType::String) =   (N = round(length(obsT[:,1])/2);
                                                GasNetModelDirBin1(obsT,
                                                    [ zeros(Int(2N)) ,0.9  * ones(1), 0.01 * ones(1)],
                                                    [Int.(1:2N),Int.(ones(Int(2N)))],scoreScalingType) )
fooGasNetModelDirBin1 = GasNetModelDirBin1(ones(100,20))

# Relations between Static and Dynamical models: conventions on storage for
# parameters and observations
StaModType(Model::GasNetModelDirBin1) = StaticNets.fooNetModelDirBin1# to be substituted with a conversion mechanism
linSpacedPar(Model::GasNetModelDirBin1,Nnodes::Int;NgroupsW = Nnodes, deltaN::Int=3,graphConstr = true) =
        StaticNets.linSpacedPar(StaModType(Model),Nnodes;Ngroups = NgroupsW,deltaN=deltaN,graphConstr =graphConstr);

# options and conversions of parameters for optimization
function setOptionsOptim(Model::GasNetModelDirBin1)
    "Set the options for the optimization required in the estimation of the model.
    For the optimization use the Optim package."
    tol = eps()*1000
    opt = Optim.Options(  g_tol = tol,
                     x_tol = tol,
                     f_tol = tol,
                     iterations = 500,
                     show_trace = true,#false,#
                     show_every=1)#,
    #                 extended_trace = true)

    algo = Newton(; alphaguess = LineSearches.InitialHagerZhang(),
                     linesearch = LineSearches.HagerZhang())
    algo = NewtonTrustRegion(; initial_delta = 0.1,
                    delta_hat = 0.8,
                    eta = 0.1,
                    rho_lower = 0.25,
                    rho_upper = 0.75)# Newton(;alphaguess = Optim.LineSearches.InitialStatic(),
       return opt, algo
 end

function array2VecGasPar(Model::GasNetModelDirBin1,ArrayGasPar::Array{<:Array{Float64,1},1})
         VecGasPar = [ArrayGasPar[1];ArrayGasPar[2];ArrayGasPar[3]]
         return VecGasPar
     end

function NumberOfGroups(Model::GasNetModelDirBin1,groupsInds::Array{<:Array{<:Real,1},1})
    NGW = length(unique(groupsInds[1]))
    GBA = length( unique( groupsInds[2][.!(groupsInds[2].==0)] ) )
     return NGW, GBA
 end
function NumberOfGroupsAndABindNodes(Model::GasNetModelDirBin1,groupsInds::Array{<:Array{<:Real,1},1})
    NGW = length(unique(groupsInds[1]))
    GBA = length( unique( [groupsInds[2][.!(groupsInds[2].==0)]  ]) )
    ABgroupsIndNodesIO = groupsInds[2][groupsInds[1]]
    indTvNodesIO =   (.!(ABgroupsIndNodesIO .==0))#if at least one group ind is zero then both in and out tv par are constant
    return NGW, GBA, ABgroupsIndNodesIO, indTvNodesIO
 end
 function NumberOfGroupsAndABindNodes(Model::GasNetModelDirBin1,groupsInds::Array{Array{Int,1},1})
     NGW = length(unique(groupsInds[1]))
     GBA = length( unique( groupsInds[2][.!(groupsInds[2].==0)] ) )
     ABgroupsIndNodesIO = groupsInds[2][groupsInds[1]]
     indTvNodesIO =   (.!(ABgroupsIndNodesIO .==0))#if at least one group ind is zero then both in and out tv par are constant
     return NGW, GBA, ABgroupsIndNodesIO, indTvNodesIO
  end
function vec2ArrayGasPar(Model::GasNetModelDirBin1,VecGasPar::Array{<:Real,1};groupsInds::Array{<:Array{<:Real,1},1} = Model.groupsInds)
     println(size(VecGasPar))
     NGW, GBA = NumberOfGroups(Model,groupsInds)
     ArrayGasPar =  [VecGasPar[1:NGW],ones(GBA).*VecGasPar[NGW+1:NGW+GBA],ones(GBA).*VecGasPar[NGW+GBA+1:NGW + 2GBA]]
     return ArrayGasPar
     end
function restrictGasPar(Model::GasNetModelDirBin1,vecUnGasPar::Array{<:Real,1})
     "From the Unrestricted values of the parameters return the restricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, e.g. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     groupsInds = Model.groupsInds
     NGW, GBA= NumberOfGroups(Model,groupsInds)
     W_un = vecUnGasPar[1:2NGW]
     diag_B_Un = vecUnGasPar[2NGW+1:2NGW+GBA]
     diag_A_Un = vecUnGasPar[2NGW+ GBA+1:2NGW+2GBA]
     W_Re = W_un
     diag_A_Re = exp.(diag_A_Un)
     diag_B_Re = 1 ./ (1+exp.(.-diag_B_Un))
     vecRePar =  [W_Re; diag_B_Re; diag_A_Re]
     return vecRePar
     end
function unRestrictGasPar( Model::GasNetModelDirBin1,vecReGasPar::Array{<:Real,1})
     "From the restricted values of the parameters return the unrestricted ones
     takes as inputs a vector of parameters that can take any value in R and returns
     the appropriately contrained values, i.e. the coefficient of autoregressive
     component in the GAS update rule has to be positive and <1."
     groupsInds = Model.groupsInds
     NGW, GBA= NumberOfGroups(Model,groupsInds)
     W_Re = vecReGasPar[1:2NGW]
     diag_B_Re = vecReGasPar[2NGW+1:2NGW+GBA]
     diag_A_Re = vecReGasPar[2NGW+GBA+1:end]

     W_Un = W_Re
     diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
     diag_A_Un = log.(diag_A_Re)

     vecUnPar =  [W_Re;diag_B_Un; diag_A_Un]
     return vecUnPar
     end

#Gas Filter Functions
identify(Model::GasNetModelDirBin1,parIO::Array{Tp,1} where Tp<:Real; idType = "equalIOsums")=
 StaticNets.identify(StaticNets.fooNetModelDirBin1,parIO;  idType =idType )

function scalingMatGas(Model::GasNetModelDirBin1,expMat::Array{<:Real,2};
                        I_tm1:: Union{UniformScaling, Matrix}=UniformScaling(2))
    "Return the matrix required for the scaling of the score, given the expected
     matrix and the Scaling matrix at previous time. "
    if Model.scoreScalingType == ""
        scalingMat = 1#
    elseif Model.scoreScalingType == "FISHER-EWMA"

        # λ = 0.8
        #
        # I = expMat.*(1-expMat)
        # diagI = sum(I,2)
        # [I[i,i] = diagI[i] for i=1:length(diagI) ]
        # I_t =  λ*I + (1-λ) *I_tm1
        # scalingMat = sqrt(I) ##
    elseif  Model.scoreScalingType == "FISHER-DIAG"
        # display(expMat)
         I = expMat.*(1 .- expMat)
        # scalingMat = zeros(Real,2 .* size(expMat))
        diagScalIn = sqrt.(sumSq(I, 2))
        N = length(diagScalIn)
        # [scalingMat[i,i] = diagScalIn[i] for i=1:N ]
        diagScalOut = sqrt.(sumSq(I, 1))
        # display(diagScalIn)
        # display(diagScalOut)

        # [scalingMat[N+i,N+i] = diagScalOut[i] for i=1:N ]
        scalingMat = Diagonal(vcat(diagScalIn, diagScalOut))
        #print(size(scalingMat))
    end
    return scalingMat
 end
function predict_score_driven_par( Model::GasNetModelDirBin1,N::Int,degsIO_t::Array{<:Real,1},
                         ftotIO_t::Array{<:Real,1},I_tm1:: Union{UniformScaling, Matrix},
                         indTvNodesIO::Union{BitArray{1},Array{Bool}},
                         WgasIO::Array{<:Real,1},BgasIO::Array{<:Real,1},AgasIO::Array{<:Real,1})
     #= likelihood and gradients depend on all the parameters (ftot_t), but
     only the time vaying ones (f_t) are to be updated=#

     fIO_t = ftotIO_t[indTvNodesIO] #Time varying fitnesses
     thetas_mat_t_exp, exp_mat_t = StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,ftotIO_t)
     exp_degIO_t = [sumSq(exp_mat_t,2);sumSq(exp_mat_t,1)]
         #The following part is a faster version of logLikelihood_t()
     loglike_t = sum(ftotIO_t .* degsIO_t) -  sum(log.(1 .+ thetas_mat_t_exp))
     # compute the score of observations at time t wrt parameters at tm1
     gradIO_t = degsIO_t  .-  exp_degIO_t # only for TV paramaters
     # updating direction without rescaling (for the moment)
     I_t = scalingMatGas(Model,exp_mat_t;I_tm1=I_tm1)
     # display(diag(I_t))
     # display(size(gradIO_t))
     # display(size(exp_mat_t ))
      sIO_t =   I_t\gradIO_t
      sIO_t[gradIO_t.==0] .= 0

     # GAS update for next iteration
     fIO_tp1 = WgasIO .+ BgasIO .* fIO_t .+ AgasIO .* sIO_t[indTvNodesIO]
     fIO_t = fIO_tp1
     ftotIO_tp1 = copy(ftotIO_t)
     ftotIO_tp1 = fIO_t #of all parameters udate dynamic ones with GAS
     ftotIO_tp1 =identify(Model,ftotIO_tp1 )
     return ftotIO_tp1,loglike_t,gradIO_t
  end
function score_driven_filter( Model::GasNetModelDirBin1,
                                vResGasPar::Array{<:Real,1};vConstPar::Array{<:Real,1} = zeros(Real,2),
                                obsT::Array{<:Real,2}=Model.obsT, ftotIO_0::Array{<:Real,1} = zeros(2),
                                groupsInds::Array{<:Array{<:Real,1},1}=Model.groupsInds,
                                dgpNT = (0,0))
    """GAS Filter the Dynamic Fitnesses from the Observed degrees, given the GAS parameters
     given T observations for the degrees in TxN vector degsT
     """
     #se parametri per dgp sono dati allora campiona il dgp, altrimenti filtra
    sum(dgpNT) == 0   ?     dgp = false : dgp = true
    if dgp
        N,T = dgpNT
        Y_T = zeros(Int8,N,N,T)
    else
        N2,T = size(obsT);N = round(Int,N2/2)
    end

    NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model,groupsInds)

    # ABgroups Inds for each par (In and Out). IO means that we have a 2N vector
    # of indices: one index for each In par and below one index for each Out par
    #println(NumberOfGroupsAndABindNodes(Model,groupsInds))

    # Organize parameters of the GAS update equation
    WGroupsIO = vResGasPar[1:NGW]
    #    StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,WGroupsIO )
    W_allIO = WGroupsIO[groupsInds[1]]

    BgasGroups  = vResGasPar[NGW+1:NGW+GBA]
    AgasGroups  = vResGasPar[NGW+GBA+1:NGW+2GBA]
    #distribute nodes among groups with equal parameters
    AgasIO = AgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
    BgasIO = BgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]

    WgasIO   = W_allIO[indTvNodesIO]


    if all(ftotIO_0.==0)
        # start values equal the unconditional mean,
        # but  constant ones remain equal to the unconditional mean, hence initialize as
        ftotIO_0  =  WgasIO ./ (1 .- BgasIO)
    end

    ftotIO_t =  identify(Model, ftotIO_0)
    I_tm1 = UniformScaling(N)
    loglike = 0

    fIOVecT = ftotIO_t
    for t=1:T-1

        if dgp
            expMat_t = StaticNets.expMatrix(StaModType(Model),ftotIO_t )# exp.(ftot_t)
            Y_T[:,:,t] = StaticNets.samplSingMatCan(StaModType(Model),expMat_t)
            degsIO_t = [dropdims(sum(Y_T[:,:,t],dims = 2),dims = 2); dropdims(sum(Y_T[:,:,t],dims = 1),dims = 1)]
        else
            degsIO_t = obsT[:,t] # vector of in and out degrees
        end
        #t==150  ?     println(WgasIO) : ()
        ftotIO_tp1,loglike_t = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t,I_tm1,
                                            indTvNodesIO,WgasIO,BgasIO,AgasIO)
        #fIOVecT[:,t+1] = ftotIO_tp1 #store the filtered parameters from previous iteration
        fIOVecT = hcat(fIOVecT, ftotIO_tp1 )

        ftotIO_t = ftotIO_tp1
        loglike += loglike_t
    end
    if dgp
        return Y_T , fIOVecT
    else
        return fIOVecT, loglike
    end
    end

score_driven_filter(Model::GasNetModelDirBin1) = score_driven_filter(Model,array2VecGasPar(Model,Model.Par))
function gasScoreSeries( Model::GasNetModelDirBin1,
                                dynPar_T::Array{<:Real,2};
                                obsT::Array{<:Real,2}=Model.obsT)
    N2,T = size(obsT);N = round(Int,N2/2)
    # Organize parameters of the GAS update equation
    sIO_T = ones(Real,2N,T)
    gradIO_T = ones(Real,2N,T)
    # identify(Model,UMallNodesIO)
    indTvNodesIO  = trues(2N)
    for t=1:T
        ftotIO_t = dynPar_T[:,t]
        degsIO_t = obsT[:,t] # vector of in and out degrees

             fIO_t = ftotIO_t[indTvNodesIO] #Time varying fitnesses
             thetas_mat_t_exp, exp_mat_t = StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,ftotIO_t)
             exp_degIO_t = [sumSq(exp_mat_t,2);sumSq(exp_mat_t,1)]

             #The following part is a faster version of logLikelihood_t()
             loglike_t = sum(ftotIO_t.*degsIO_t)  .-   sum(log.(1 .+ thetas_mat_t_exp))
             # compute the score of observations at time t wrt parameters at tm1
             # THIS IS A NEW VERSION!!!!!!!!!!!!!!!!!!!!!!1------------------------------------------------------
             gradIO_t = sIO_t = identify(Model, (degsIO_t .- exp_degIO_t ) ) # only for TV paramaters
             # updating direction without rescaling (for the moment)
             I_t = scalingMatGas(Model,exp_mat_t)
             # display(diag(I_t))
             # display(size(gradIO_t))
             # display(size(exp_mat_t ))
              sIO_t = identify(Model, I_t\gradIO_t )
              sIO_t[gradIO_t.==0] .= 0
              sIO_T[:,t] =  sIO_t
              gradIO_T[:,t] =  gradIO_t
    end
    return sIO_T,gradIO_T
    end

sampl(Mod::GasNetModelDirBin1,T::Int)=( N = length(Mod.groupsInds[1]) ;
                                        tmpTuple = score_driven_filter(Mod,[Mod.Par[1];Mod.Par[2];Mod.Par[3]];  dgpNT = (N,T) );
                                        degsIO_T = [sumSq(tmpTuple[1],3) sumSq(tmpTuple[1],2)];
                                        (GasNetModelDirBin1(degsIO_T,Mod.Par,Mod.groupsInds,Mod.scoreScalingType),tmpTuple[1],tmpTuple[2]) )

# Estimation

function estSingSnap(Model::GasNetModelDirBin1, degs_t::Array{<:Real,1}; groupsInds = Model.groupsInds, targetErr::Real=targetErrValDynNets)
    hatUnPar,~,~ = StaticNets.estimate(StaModType(Model); degIO = degs_t ,groupsInds = groupsInds[1], targetErr =  targetErr)
    return hatUnPar
 end
function estimateSnapSeq(Model::GasNetModelDirBin1; degsIO_T::Array{<:Real,2}=Model.obsT,
                            targetErr::Real=1e-5,identPost=false,identIter=false)
    #this funciton does  not allow groups estimate a sequence of single node's
    # fitnesses
    return  StaticNets.estimate(StaticNets.SnapSeqNetDirBin1(degsIO_T),targetErr = targetErr ,identPost=identPost,identIter= identIter)
 end

function estimateTarg(Model::GasNetModelDirBin1; SSest::Array{<:Real,2} =zeros(2,2),
             groupsInds::Array{<:Array{<:Real,1},1} = Model.groupsInds)
    "Estimate the GAS parameters of the GAS  model,
    using a targeting on the unconditional means. Then return the GASparameters "
    # Not ready for groups indices
    # only one B and one A
    N2,T = size(Model.obsT);N = round(Int,N2/2)
    println(N)
    NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model,groupsInds)

    optims_opt, algo = setOptionsOptim(Model)

    #target the unconditional means from the sequence of single snapshots estimate
    if sum(SSest) == 0
        uncMeansIO =  meanSq(estimateSnapSeq(Model),2)
    else
        uncMeansIO =  meanSq(SSest,2)
    end

    uncMeansIO = identify(Model,uncMeansIO)
    ftotIO_0  = meanSq(SSest[:, 1:5],2)


    #functions needed if I want to enable grouping of nodes again
    # function WGroupsFromB(B::Array{<:Real,1})
    #      BGroupsW_IO = B[[groupsInds[2];groupsInds[2] ] ]
    #      WGroupsIO = uncMeansWGroupsIO.*(1-BGroupsW_IO)
    #      return WGroupsIO
    # end
    function restrAB(vecUnPar::Array{<:Real,1}) #restrict a vector of only A and B
         diag_B_Un = vecUnPar[1:GBA]
         diag_A_Un = vecUnPar[GBA+1:end]

         diag_B_Re = 1 ./ (1 .+ exp.( .- diag_B_Un))
         diag_A_Re = exp.(diag_A_Un)
         vecRePar =  [diag_B_Re; diag_A_Re]
         return vecRePar
         end
    function UnrestrAB(vecRePar::Array{<:Real,1}) #restrict a vector of only A and B
         diag_B_Re = vecRePar[1:GBA]
         diag_A_Re = vecRePar[GBA+1:end]
         diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
         diag_A_Un = log.(diag_A_Re)
         vecUnPar =  [diag_B_Un; diag_A_Un]
         return vecUnPar
    end
    # #set the starting points for the optimizations
    B0_Re  = 0.99
    A0_Re_H  = 0.1
    A0_Re_L  = 0.000001
    if true && GBA>1
        #calcola autocorrelazioni dello score e usale per definire A_0
        meanDegsIO = meanSq(Model.obsT,2)
        constParsIO,~,~ = StaticNets.estimate(StaticNets.fooNetModelDirBin1,degIO = meanDegsIO)
         fooPar = [zeros(2N), zeros(2N) ,zeros(2N)]
          indsGroups = [Int.(1:2N),Int.(1:2N)]
        modGasDirBin1 = DynNets.GasNetModelDirBin1(Model.obsT, fooPar, indsGroups,"FISHER-DIAG")
        sIO_T,gIO_T =  gasScoreSeries(modGasDirBin1,repeat(constParsIO,1,T);obsT = Model.obsT)
        indsLikelyTV = dropdims(autocor(Float64.(gIO_T'),[1]).>0.2, dims = 1)
        Avec = ones(Real,GBA).*A0_Re_L
        println(size(indsLikelyTV))
        Avec[indsLikelyTV] .= A0_Re_H
    else
        Avec = ones(Real,GBA).*A0_Re_L
    end
    vGasParAll0_Un = UnrestrAB( [ones(Real,GBA).*B0_Re;Avec ]) # vGasParGroups0_Un#

    # objective function for the optimization
    function objfunGas(vecUnAB::Array{<:Real,1})# a function of the groups parameters
         ReB = 1 ./ (1 .+ exp.(.-vecUnAB[1:GBA]))
         ReA = exp.(vecUnAB[GBA+1:end])
         WNodes = uncMeansIO.*((maxLargeVal .+ ReB[1])/maxLargeVal)
        # println(typeof(WNodes))
         WNodes[indTvNodesIO] = uncMeansIO[indTvNodesIO].*(1 .- ReB[ABgroupsIndNodesIO[indTvNodesIO]])

         #StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,WNodes)
         vecReGasPar = [WNodes;ReB;ReA ]
         #StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,vecReGasPar[1:NGW])
         #NON POSSONO STARE NELLO STESSO VETTORE PARAMETRI DA OTTIMIZZARE DI TIPO FORWARDDIFF E TYPES NORMALI  ?      ?      ?      ?
         foo,loglikelValue = score_driven_filter( Model,vecReGasPar;
                                        groupsInds = groupsInds, ftotIO_0=ftotIO_0)
         return - loglikelValue
    end
    diff_mode = "forward"
    if diff_mode == "forward"
        if uppercase(Model.scoreScalingType) == "FISHER-EWMA"
            ADobjfunGas = objfunGas
        else
            ADobjfunGas = TwiceDifferentiable(objfunGas, vGasParAll0_Un; autodiff = :finite);
        end
    elseif diff_mode in ["reverse",  "backward"]

    end
    println(objfunGas(vGasParAll0_Un))
    optim_out2  = optimize(ADobjfunGas,vGasParAll0_Un,algo,optims_opt)
    outParAB = Optim.minimizer(optim_out2)
    outReAB = restrAB(outParAB)
    vGasParGroupsHat_Re = zeros(NGW + 2GBA)
    vGasParGroupsHat_Re[(1:NGW)[.!indTvNodesIO]] = uncMeansIO[.!indTvNodesIO]
    #println( uncMeansIO[.!indTvNodesIO])
    vGasParGroupsHat_Re[(1:NGW)[indTvNodesIO]] = uncMeansIO[indTvNodesIO].*(1 .- outReAB[1:GBA][ABgroupsIndNodesIO[indTvNodesIO]])
    vGasParGroupsHat_Re[1+NGW:end] = outReAB
    #println(vGasParGroupsHat_Re)
    arrayGasParHat_Re = vec2ArrayGasPar(Model,vGasParGroupsHat_Re;groupsInds = groupsInds )
    conv_flag =  Optim.converged(optim_out2)
    # println(optim_out2)
    return  arrayGasParHat_Re, conv_flag
end

function estimateOld(Model::GasNetModelDirBin1;start_values = [zeros(Real,10),zeros(Real,3),zeros(Real,3) ])
    T,N2 = size(Model.obsT);N = round(Int,N2/2)
    groupsInds = Model.groupsInds
    NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model,groupsInds)

    optims_opt, algo = setOptionsOptim(Model)


    #set the starting points for the optimizations if not given
    if sum([sum(start_values[i]) for i =1:3 ])==0
         B0_ReGroups = 0.9 * ones(GBA)
         A0_ReGroups = 0.01 * ones(GBA)

         BgasNodes_0   = B0_ReGroups[ABgroupsIndNodesIO[indTvNodesIO]]
         AgasNodes_0   = A0_ReGroups[ABgroupsIndNodesIO[indTvNodesIO]]
         # choose the starting points for W in order to fix the Unconditional Means
         # equal to the static estimates for the average degrees(over time)
         MeanDegsT =  meanSq(Model.obsT,1)

         UMstaticGroups = estSingSnap(Model ,MeanDegsT; groupsInds = groupsInds )
         conv_flag1 = true

         GBA>0  ?    ( W0_Groups = UMstaticGroups.*(1 .- B0_ReGroups[[groupsInds[2];groupsInds[3]]])) : W0_Groups = UMstaticGroups
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

function estimateTargOld(Model::GasNetModelDirBin1)
    "Estimate the GAS parameters of the GAS  model,
    using a targeting on the unconditional means. Then return the GASparameters "
    T,N2 = size(Model.obsT);N = round(Int,N2/2)
    groupsInds = Model.groupsInds
    NGW,GBA = NumberOfGroups(Model,groupsInds)
    optims_opt, algo = setOptionsOptim(Model)

    #target the unconditional mean
    ParTvSnap = estimateSnapSeq(Model) # sequence of snapshots estimate
    uncMeansWGroupsIO = dropdims(mean(ParTvSnap,1), dims =1)
    #functions needed for the specific targeting case
    function WGroupsFromB(B::Array{<:Real,1})
         BGroupsW_IO = B[[groupsInds[2];groupsInds[2] ] ]
         WGroupsIO = uncMeansWGroupsIO.*(1 .- BGroupsW_IO)
         return WGroupsIO
    end
    function restrAB(vecUnPar::Array{<:Real,1}) #restrict a vector of only A and B
         diag_B_Un = vecUnPar[1:GBA]
         diag_A_Un = vecUnPar[GBA+1:end]

         diag_B_Re = 1 ./ (1 .+ exp.(.-diag_B_Un))
         diag_A_Re = exp.(diag_A_Un)
         vecRePar =  [diag_B_Re; diag_A_Re]
         return vecRePar
         end
    function UnrestrAB(vecRePar::Array{<:Real,1}) #restrict a vector of only A and B
         diag_B_Re = vecRePar[1:GBA]
         diag_A_Re = vecRePar[GBA+1:end]
         diag_B_Un = log.(diag_B_Re ./ (1 .- diag_B_Re ))
         diag_A_Un = log.(diag_A_Re)
         vecUnPar =  [diag_B_Un; diag_A_Un]
         return vecUnPar
    end
    if GBA == 0  return [uncMeansWGroupsIO, zeros(1), zeros(1)] , Optim.converged(optim_out) end;
    # Estimate As and Bs



    # #set the starting points for the optimizations
    B0_ReGroups = 0.9 * ones(GBA)
    A0_ReGroups = 0.01 * ones(GBA)
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

function forecastEvalGasNetDirBin1(obsNet_T::BitArray{3}, gasParEstOnTrain::Array{Array{Float64,1},1},Ttrain::Int;thVarConst = 0.0005,mat2rem=falses(obsNet_T[:,:,1]) )
    # this function evaluates the forecasting performances of the GAS model over the train sample.
    # to do so uses gas parameters previously estimated and observations at time t of the test sample
    # to forecast tv par at t+1, then using these estimates and t+1 observations goes to t+2 and so on

    @show N2 = length(gasParEstOnTrain[1]);N = round(Int,N2/2)
    T = length(obsNet_T[1,1,:])
    @show Ttest = T-Ttrain
    testObsNet_T = obsNet_T[:,:,Ttrain+1:end]
    degsIO_T = [sumSq(obsNet_T,2); sumSq(obsNet_T,1)]



    #disregard predictions for degrees that are constant in the training sample
    # gas shoud manage but AR fitness no and I dont want to advantage gas
    # to do so i need to count the number of links that are non constant (nnc)
    Nc,isConIn,isConOut = StaticNets.defineConstDegs(degsIO_T[:,1:Ttrain];thVarConst = thVarConst )
    noDiagIndnnc = putZeroDiag(((.!isConIn).*(.!isConOut')).* (.!mat2rem)    )

    @show Nlinksnnc =sum(noDiagIndnnc)
    # forecast fitnesses using Gas parameters and observations
    foreFit,~ = score_driven_filter( DynNets.GasNetModelDirBin1(degsIO_T),[gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])

    TRoc = Ttest-1
    #storage variables
    foreVals = 100 * ones(Nlinksnnc*(TRoc))
    realVals = trues( Nlinksnnc*(TRoc))

    lastInd=1

    for t=2:TRoc
        adjMat = testObsNet_T[:,:,t]   #[testObsNet_T[50:end,:,t];testObsNet_T[1:49,:,t]] #a shuffle to test that the code does not not work
        adjMat_tm1 = testObsNet_T[:,:,t-1]
        indZeroR = sum(adjMat_tm1,dims = 2).==0
        indZeroC = sum(adjMat_tm1,dims = 1).==0
        indZeroMat = indZeroR.*indZeroC
    #    println(sum(indZeroMat)/(length(indZeroC)^2))
        expMat = StaticNets.expMatrix(StaticNets.fooNetModelDirBin1,foreFit[:,Ttrain+1:end][:,t])

        #expMat[indZeroMat] = 0
        foreVals[lastInd:lastInd+Nlinksnnc-1] = expMat[noDiagIndnnc]
        realVals[lastInd:lastInd+Nlinksnnc-1] = adjMat[noDiagIndnnc]
        lastInd += Nlinksnnc
    end
    return realVals,foreVals,foreFit
 end

function forecastEvalGasNetDirBin1(obsMat_allT::BitArray{3}, Ttrain::Int;thVarConst = 0.0005,mat2rem=falses(obsNet_T[:,:,1]) )
    # if all observations and the length of training sample are given as input, then
    # estimate the gas par and run the one step forecast
    degsIO_T = [sumSq(obsMat_allT,2);sumSq(obsMat_allT,1)]
    N2,T = size(degsIO_T);N = round(Int,N2/2)
    Ttest = T - Ttrain
    #Define the train model and estimate
    modGasDirBin1_eMidTrain = DynNets.GasNetModelDirBin1(degsIO_T[:,1:Ttrain])
    estTargDirBin1_eMidTrain,~ = DynNets.estimateTarg(modGasDirBin1_eMidTrain)
    gasParEstOnTrain = estTargDirBin1_eMidTrain
    #produce a vector of real values and probabilities from 1 step ahead forecasting
    foreEval = forecastEvalGasNetDirBin1(obsMat_allT,gasParEstOnTrain,Ttrain,thVarConst = thVarConst,mat2rem =mat2rem)
    tp,fp,AUC = Utilities.rocCurve(foreEval[1],foreEval[2])
    return tp,fp,AUC,foreEval
 end

function multiSteps_predict_score_driven_par( Model::GasNetModelDirBin1,N::Int,degsIO_t_in::Array{<:Real,1},
                         ftotIO_t0_input::Array{<:Real,1},
                         WgasIO_in::Array{<:Real,1},BgasIO_in::Array{<:Real,1},AgasIO_in::Array{<:Real,1},Nsample::Int,Nsteps::Int )
     #Simula Nsteps nel futuro per Nsample volte
     degsIO_t0 = copy(degsIO_t_in)
     ftotIO_t0 = copy(ftotIO_t0_input)
     WgasIO = copy(WgasIO_in)
     BgasIO = copy(BgasIO_in)
     AgasIO = copy(AgasIO_in)
     ftotIO_T = zeros(length(ftotIO_t0),Nsteps)

     indTvNodesIO = Model.groupsInds[2].!=0
     expMat_T_means = zeros(N,N,Nsteps)
     expMat_T_Y = zeros(N,N,Nsteps)
     ftotIO_t1  ,~ = predict_score_driven_par(Model,N,degsIO_t0,ftotIO_t0,eye(N,N),
            indTvNodesIO,WgasIO,BgasIO,AgasIO)
     ftotIO_T[:,1] = copy(ftotIO_t1)
     for i=1:Nsample
         ftotIO_t1_tmp = copy(ftotIO_t1)
         expMat_t = StaticNets.expMatrix(StaModType(Model),ftotIO_t1_tmp )
         expMat_T_means[:,:,1] = expMat_t
         Y_t = StaticNets.samplSingMatCan(StaModType(Model),expMat_t)
         degsIO_t = [dropdims(sum(Y_t,dims = 2),dims =2); dropdims(sum(Y_t,dims = 1), dims =1)]
         #matrice media dell matrici campionate
         expMat_T_Y[:,:,1] = expMat_T_Y[:,:,1]  .+  Y_t ./ Nsample
         #il primo step viene fatto usando l'osservzione reale
         ftotIO_tp1  ,~ = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t1_tmp,eye(N,N),
                            indTvNodesIO,WgasIO,BgasIO,AgasIO)
         ftotIO_t = copy(ftotIO_tp1)
         ftotIO_T[:,1] = copy(ftotIO_t)
         expMat_t = StaticNets.expMatrix(StaModType(Model),ftotIO_t )
        # println(sum(ftotIO_t0))

         for t=2:Nsteps
            # fai update dei parametri
             ftotIO_tp1  ,~ = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t,eye(expMat_t),
                                                 indTvNodesIO,WgasIO,BgasIO,AgasIO)
             ftotIO_t = copy(ftotIO_tp1)
             # da ora stai usando i parametri  t
             ftotIO_T[:,t] = ftotIO_T[:,t] .+ ftotIO_t ./ Nsample


             expMat_t = StaticNets.expMatrix(StaModType(Model),ftotIO_t )
            # t==2  ?     println(mean(expMat_t)) : ()
             #matrice media delle matrici medie
             expMat_T_means[:,:,t] = expMat_T_means[:,:,t]  .+  expMat_t ./ Nsample

             Y_t = StaticNets.samplSingMatCan(StaModType(Model),expMat_t)
             #matrice media dell matrici campionate
             expMat_T_Y[:,:,t] = expMat_T_Y[:,:,t]  .+  Y_t ./ Nsample
             degsIO_t = [dropdims(sum(Y_t,dims = 2), dims =2); dropdims(sum(Y_t,dims = 1), dims =1)]

         end

     end
     return expMat_T_means, ftotIO_T,expMat_T_Y
  end

function multiStepsForecastExpMat( Model::GasNetModelDirBin1,obsNet_T::BitArray{3},
      gasParEstOnTrain::Array{Array{Float64,1},1},Nsample::Int,Nsteps::Int,Ttrain::Int;meanOfMeans=false)

      N2 = length(gasParEstOnTrain[1]);N = round(Int,N2/2)
      T = length(obsNet_T[1,1,:])
      degsIO_T = [sumSq(obsNet_T,2);sumSq(obsNet_T,1)]
      foreFitGas1,~ = DynNets.score_driven_filter( Model,
      [gasParEstOnTrain[1];gasParEstOnTrain[2];gasParEstOnTrain[3]])
      Ttest = T-Ttrain

      expMat_T_means = zeros(N,N,T)
      expMat_T_Y = zeros(N,N,T)
      foreFit_T = zeros(2N,T)
      expMat_T_meanPar = zeros(N,N,T)
      for t=Ttrain +1:T
          tFit = t-Nsteps
          tmpAllMat_means,tmpFit,tmpAllMat_Y = multiSteps_predict_score_driven_par(Model,N,degsIO_T[:,tFit],foreFitGas1[:,tFit],
                                gasParEstOnTrain[1],gasParEstOnTrain[2],gasParEstOnTrain[3],Nsample,Nsteps)
          expMat_T_means[:,:,tFit+Nsteps] = tmpAllMat_means[:,:,Nsteps]
          expMat_T_Y[:,:,tFit+Nsteps] = tmpAllMat_Y[:,:,Nsteps]
          foreFit_T[:,tFit+Nsteps] = tmpFit[:,end]
           expMat_T_meanPar[:,:,tFit+Nsteps] =   StaticNets.expMatrix(StaModType(Model),tmpFit[:,end])
      end
      return expMat_T_means,foreFit_T,expMat_T_Y,expMat_T_meanPar
  end

function multiStepsForecastExpMat_roll( Model::GasNetModelDirBin1,obsNet_T::BitArray{3},
      gasEst_rolling_input::Array{Array{Array{Float64,1},1},1},gasFiltAndForeFitFromRollEst_input::Array{Float64,2},Nsample::Int,Nsteps::Int,Ttrain::Int;meanOfMeans=false)


      gasEst_rolling = copy(gasEst_rolling_input)
      gasFiltAndForeFitFromRollEst = copy(gasFiltAndForeFitFromRollEst_input)
      N2 = length(gasEst_rolling[1][1]);N = round(Int,N2/2)
      T = length(obsNet_T[1,1,:])
      degsIO_T = [sumSq(obsNet_T,2);sumSq(obsNet_T,1)]
      Ttest = T-Ttrain

      expMat_T_means = zeros(N,N,T)
      expMat_T_Y = zeros(N,N,T)
      foreFit_T = zeros(2N,T)
      expMat_T_meanPar = zeros(N,N,T)
      # forecast the parameters at each t between Ttrain+1 and T using the observation at t-N_Nsteps
      # and the static gas estimates obtained from  t-Nsteps-Ttrain:T-Nsteps
      for t=Ttrain +1:T
          tEst = max(t-Nsteps - Ttrain,1)
          tObs = t-Nsteps-1


          tmpAllMat_means,tmpFit,tmpAllMat_Y = multiSteps_predict_score_driven_par(Model,N,degsIO_T[:,tObs],gasFiltAndForeFitFromRollEst[:,tObs+1],
                                gasEst_rolling[tEst][1],gasEst_rolling[tEst][2],gasEst_rolling[tEst][3],Nsample,Nsteps)
          expMat_T_means[:,:,tObs+Nsteps] = tmpAllMat_means[:,:,Nsteps]
          expMat_T_Y[:,:,tObs+Nsteps] = tmpAllMat_Y[:,:,Nsteps]
          foreFit_T[:,tObs+Nsteps] = tmpFit[:,end]
           expMat_T_meanPar[:,:,tObs+Nsteps] =   StaticNets.expMatrix(StaModType(Model),tmpFit[:,end])
      end
      return expMat_T_means,foreFit_T,expMat_T_Y,expMat_T_meanPar
  end

function logLike_t(Model::GasNetModelDirBin1, obsT, vReGasPar)
 groupsInds = Model.groupsInds
      N2,T = size(obsT);N = round(Int,N2/2)
      NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model, groupsInds)
      # Organize parameters of the GAS update equation
      WGroupsIO = vReGasPar[1:NGW]
      #    StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,WGroupsIO )
      W_allIO = WGroupsIO[groupsInds[1]]
      BgasGroups  = vReGasPar[NGW+1:NGW+GBA]
      AgasGroups  = vReGasPar[NGW+GBA+1:NGW+2GBA]
      #distribute nodes among groups with equal parameters
      AgasIO = AgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      BgasIO = BgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      WgasIO   = W_allIO[indTvNodesIO]
      UMallNodesIO = W_allIO
      UMallNodesIO  =  WgasIO ./ (1 .- BgasIO)
      ftotIO_t =  identify(Model,UMallNodesIO)
      I_tm1 = UniformScaling(N)
      loglike_t = zero(Real)
      for t=1:T-1
            degsIO_t = obsT[:,t] # vector of in and out degrees
            ftotIO_tp1,loglike_t = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t,I_tm1,
                                                  indTvNodesIO,WgasIO,BgasIO,AgasIO)
            ftotIO_t = ftotIO_tp1
      end
      return  loglike_t::T where T <:Real
end

function logLike_T(Model::GasNetModelDirBin1, obsT, vReGasPar)
      groupsInds = Model.groupsInds
      N2,T = size(obsT);N = round(Int,N2/2)
      NGW,GBA,ABgroupsIndNodesIO,indTvNodesIO = NumberOfGroupsAndABindNodes(Model, groupsInds)
      # Organize parameters of the GAS update equation
      WGroupsIO = vReGasPar[1:NGW]
      #    StaticNets.expMatrix2(StaticNets.fooNetModelDirBin1,WGroupsIO )
      W_allIO = WGroupsIO[groupsInds[1]]
      BgasGroups  = vReGasPar[NGW+1:NGW+GBA]
      AgasGroups  = vReGasPar[NGW+GBA+1:NGW+2GBA]
      #distribute nodes among groups with equal parameters
      AgasIO = AgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      BgasIO = BgasGroups[ABgroupsIndNodesIO[indTvNodesIO]]
      WgasIO   = W_allIO[indTvNodesIO]
      # start values equal the unconditional mean, but  constant ones remain equal to the unconditional mean, hence initialize as:
      UMallNodesIO = W_allIO
      UMallNodesIO  =  WgasIO ./ (1 .- BgasIO)
      ftotIO_t =  identify(Model,UMallNodesIO)
      I_tm1 = UniformScaling(N)
      loglike_T = zero(Real)
      for t=1:T-1
            degsIO_t = obsT[:,t] # vector of in and out degrees
            ftotIO_tp1,loglike_t = predict_score_driven_par(Model,N,degsIO_t,ftotIO_t,I_tm1,
                                                  indTvNodesIO,WgasIO,BgasIO,AgasIO)
            ftotIO_t = ftotIO_tp1
            loglike_T += loglike_t
      end
      return  loglike_T::T where T <:Real
end
