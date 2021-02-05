

ErgmRcall.clean_start_RCall()
R"""options(warn=-1) """

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
        ftot_t,loglike_t,I_tm1,~ = DynNets.predict_score_driven_par(Model,obs_t,ftot_t,I_tm1,indTvPar,Wvec,Bvec,Avec)
        fVecT[:,t] = ftot_t #store the filtered parameters from previous iteration
    end


    return obsT,statsT,fVecT
end
export sdergmDGPpaper_Rcall
