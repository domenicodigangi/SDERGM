using StatsBase, JLD,  Hwloc

counts = Hwloc.histmap(Hwloc.topology_load())
ncores = counts[:Core]
if nworkers()<ncores
    addprocs(ncores   - nprocs())
end
using DynNets
using StatsFuns
@everywhere using DynNets
@everywhere using StatsFuns
N_est = 300
Tvals =  [1000]
N_ind = length(Tvals)[1]


N = 10 ;
GBAdgp =1; GWdgp = N ;
GBAest =1; GWest = 10; useStartVal = false

#Distribuisci i nodi tra i gruppi in parti uguali

groupsIndsDgp = [[DynNets.distributeAinVecN(Array{Int}(1:GWdgp),N)];
                    [DynNets.distributeAinVecN(Array{Int}(1:GBAdgp),GWdgp)]]
groupsIndsEst = [[DynNets.distributeAinVecN(Array{Int}(1:GWest),N)];
                   [DynNets.distributeAinVecN(Array{Int}(1:GBAest),GWest)]]

#the values for the unconidtional means to be used in the DGP are obtained estimating
#an equally spaced degree sequence
est_times = Array{Float64}(zeros(N_est))
est_conv_flag = falses(N_est)
est_flag = falses(N_est)
# Set the value of the parameters common to all replicate of DGP+Estimation
realUm = DynNets.linSpacedFitnesses(N;Ngroups = GWdgp,mode = "UNDIRECTED")
sim_B = Array(linspace(0.4,0.9,GBAdgp+1)[2:end]) #
sim_A = Array(linspace(0.02,0.2,GBAdgp+1)[2:end]) #
sim_W = realUm.*(1-sim_B[groupsIndsDgp[2]])
#Initialize the variable to store
sizesParSim = [GWdgp,GBAdgp,GBAdgp]
sizesParEst = [GWest,GBAdgp,GBAest]
for i =1:3
    ex1 = "estPar_$(i)  = SharedArray{Float64}(zeros(N_ind,N_est, $(sizesParEst[i]) ))"
    ex2 = "simPar_$(i) = SharedArray{Float64}(zeros(N_ind,N_est, $(sizesParSim[i])))"
    eval(parse(ex1))
    eval(parse(ex2))
 end
allObs = SharedArray{Int16,4}(N,Tvals[end],N_ind,N_est)
est_times = SharedArray{Float64}(zeros(N_ind,N_est))
est_conv_flag = SharedArray{Bool,2}(N_ind,N_est)
est_flag = SharedArray{Bool,2}(N_ind,N_est)
@progress "T"  for ind = 1:N_ind
  # ind=1
    T = Tvals[ind]
     @sync @parallel   for n = 1:N_est
#  n=1
        simPar_3[ind,n,1:GBAdgp] = sim_A
        simPar_2[ind,n,1:GBAdgp] = sim_B
        simPar_1[ind,n,1:GWdgp] = sim_W
        Model2Sample = DynNets.GasNetModelBin1(ones(T,N),[sim_W,sim_B,sim_A],
                                                    groupsIndsDgp)
        Y_T,fit_T  =  DynNets.sampl(Model2Sample,T)
        sampleDegsT = squeeze(sum(Y_T,2),2)
        allObs[:,1:T,ind,n] = sampleDegsT'
        Model2Est =  DynNets.GasNetModelBin1(sampleDegsT,[zeros(GWest),zeros(GBAest),zeros(GBAest)],
                                                    groupsIndsEst)
        #estimate
        useStartVal ?   start_val = [sim_W, sim_B, sim_A] : start_val = zeros(3,3)
        est_times[ind,n] = @elapsed  hat_pars , conv_flag =
                            DynNets.estimate(Model2Est,start_values = start_val )
        estPar_1[ind,n,:] = hat_pars[1]
        estPar_2[ind,n,:] = hat_pars[2]
        estPar_3[ind,n,:] = hat_pars[3]
        est_conv_flag[ind,n] = conv_flag
        est_flag[ind,n] = true


    end

    # save
    save_fold = "./data/estimatesTest/"
    file_name = "EstTest_$(N)_$(GWdgp)_$(GBAdgp)_$(Tvals)_$(N_est)_$(useStartVal)_Miss_$(GWest)_$(GBAest).jld"
    save_path = save_fold*file_name#
    @save(save_path,groupsIndsEst,groupsIndsDgp,allObs,simPar_1,simPar_2,simPar_3,estPar_1,estPar_2,estPar_3,est_conv_flag,est_times)
    println([N,GWdgp,GWest,T,sum(est_times[ind,:]),sum(est_conv_flag[ind,:])/sum(est_flag[ind,:])])
end

##















 ##
