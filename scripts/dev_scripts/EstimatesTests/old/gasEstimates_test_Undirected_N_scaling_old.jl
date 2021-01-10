using StatsBase, JLD,  Hwloc,  DynNets, StatsFuns
initWorkers()
@everywhere using DynNets, StatsFuns
N_est = 4
NGWpoints = 3
Tvals = [50 250 1250]# [ 100 200  300  ]
N_Tind = length(Tvals)[1]
scoreRescType = ""#"FISHER-EWMA"#
useStartVal = false# true# false
Nmin = 10
Nspacing =  5
Nmax = 20
Nvals = Vector(Nmin:Nspacing:Nmax)
N_Nind = length(Nvals)

est_times = SharedArray{Float64}(zeros(N_Tind,N_Nind,NGWpoints,N_est))
est_conv_flag = SharedArray{Bool,4}(N_Tind,N_Nind,NGWpoints,N_est)
est_flag = SharedArray{Bool,4}(N_Tind,N_Nind,NGWpoints,N_est)
GWvalStore  = Array{Array{Int,1},1}(N_Nind)

@progress "N" for indN = N_Nind:-1:1
#indN = 1
N = Nvals[indN] ;
GWvals = unique(Vector{Int}(floor.(linspace(1,N,NGWpoints))))
GWvalStore[indN] = GWvals
N_GWind = length(GWvals)

@progress "GW"  for indGW =N_GWind:-1:1
#indGW = 2
GW = GWvals[indGW]
GBAdgp =1; GWdgp = GW ;
GBAest =1; GWest = GW;
#Distribuisci i nodi tra i gruppi in parti uguali

groupsindTsDgp = [[DynNets.distributeAinVecN(Array{Int}(1:GWdgp),N)];
                    [DynNets.distributeAinVecN(Array{Int}(1:GBAdgp),GWdgp)]]
groupsindTsEst = [[DynNets.distributeAinVecN(Array{Int}(1:GWest),N)];
                   [DynNets.distributeAinVecN(Array{Int}(1:GBAest),GWest)]]


#Initialize the variable to store
realUm = DynNets.linSpacedPar(DynNets.fooGasNetModelBin1, N;NgroupsW = GW )[1]
sim_B = Array(linspace(0.4,0.9,GBAdgp+1)[2:end]) #
sim_A = Array(linspace(0.02,0.2,GBAdgp+1)[2:end]) #
sim_W = realUm.*(1-sim_B[groupsindTsDgp[2]])
[sim_W,sim_B,sim_A]


@progress "T" for indT = N_Tind:-1:1
#  indT=1
T = Tvals[indT]
@sync @parallel for n = 1:N_est
#n=1

        Model2Sample = DynNets.GasNetModelBin1(ones(T,N),[sim_W,sim_B,sim_A],
                                                    groupsindTsDgp,scoreRescType)
        Model2Sample,Y_T,Fitness_T  =  DynNets.sampl(Model2Sample,T)
        sampleDegsT = squeeze(sum(Y_T,2),2)

        Model2Est =  DynNets.GasNetModelBin1(sampleDegsT,[zeros(GWest),zeros(GBAest),zeros(GBAest)],
                                                    groupsindTsEst,scoreRescType)
        #estimate
        useStartVal ?   start_val = [sim_W, sim_B, sim_A] : start_val = zeros(3,3)
        est_times[indT,indN,indGW,n] = @elapsed  hat_pars , conv_flag =
                            DynNets.estimate(Model2Est,start_values = start_val )
        est_conv_flag[indT,indN,indGW,n] = conv_flag
        est_flag[indT,indN,indGW,n] = true


    end
end
end
    # save
    save_fold = "./data/estimatesTest/"
    file_name = "EstTest_varN_Nest_$(N_est)_$(Nmin)_$(Nspacing)_$(Nmax)_varGW_$(NGWpoints)_T_$(Tvals)_corectSpec"* scoreRescType *".jld"
    save_path = save_fold*file_name#
    @save(save_path,est_conv_flag,est_times,est_flag,Nvals,GWvalStore ,useStartVal)
    println([N,sum(est_times[:,indN,:,:]),sum(est_conv_flag[:,indN,:,:])/sum(est_flag[:,indN,:,:])])
end

##














 ##
