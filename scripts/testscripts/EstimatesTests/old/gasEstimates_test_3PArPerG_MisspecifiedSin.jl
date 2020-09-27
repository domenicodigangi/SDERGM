using StatsBase, StatsFuns, JLD,  Hwloc

counts = Hwloc.histmap(Hwloc.topology_load())
ncores = counts[:Core]
if nworkers()<ncores
    addprocs(ncores   - nprocs())
end

using DynNets
@everywhere using DynNets
@everywhere using StatsFuns
N_est = 500
Tvals = [50 300  ]#500#
N_ind = length(Tvals)[1]


N = 10 ;
GBAest =1; GWest = 10; useStartVal = false

groupsIndsEst = [[DynNets.distributeAinVecN(Array{Int}(1:GWest),N)];
                   [DynNets.distributeAinVecN(Array{Int}(1:GBAest),GWest)]]




est_times = Array{Float64}(zeros(N_est))
est_conv_flag = falses(N_est)
est_flag = falses(N_est)

sizesPar = [N,GBAest,GBAest]
for i =1:3
    ex1 = "estPar_$(i)  = SharedArray{Float64}(zeros(N_ind,N_est, $(sizesPar[i]) ))"
    eval(parse(ex1))
end

ObsAllT = Array{Array{Int16,3},1}(N_ind)
ParTvDgpAllT = Array{Array{Float64,3},1}(N_ind)

#define save path and list variables to save
save_fold = "./data/estimatesTest/"
    file_name = "EstTest_$(N)_$(Tvals)_$(N_est)_$(useStartVal)_Sin_Miss_$(GWest)_$(GBAest).jld"
    save_path = save_fold*file_name#



est_times = SharedArray{Float64}(zeros(N_ind,N_est))
est_conv_flag = SharedArray{Bool,2}(N_ind,N_est)
est_flag = SharedArray{Bool,2}(N_ind,N_est)
#
@progress "T"  for ind = 1:N_ind
#   ind=1

    T = Tvals[ind]
    ObsTmp = SharedArray{Int16,3}(Tvals[ind],N,N_est)
    ParTvDgpTmp = SharedArray{Float64,3}(Tvals[ind],N,N_est)

    @sync @parallel     for n = 1:N_est
#     n=1


        Y_T,fit_T  = DynNets.sampleDynNetsBin1(T,"FIT-MULTI-SIN";N=N)
        sampleDegsT = squeeze(sum(Y_T,2),2)

        ParTvDgpTmp[:,:,n] = fit_T
        ObsTmp[:,:,n] = sampleDegsT

        Model2Est =  DynNets.GasNetModelBin1_UmBAgr(sampleDegsT,[zeros(GWest),zeros(GBAest),zeros(GBAest)],
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

   ObsAllT[ind] = ObsTmp
   ParTvDgpAllT[ind] = ParTvDgpTmp

    # save

    @save(save_path,groupsIndsEst,ObsAllT,ParTvDgpAllT ,estPar_1,estPar_2,estPar_3,est_conv_flag,est_times)
    println([N,GWest,T,sum(est_times[ind,:]),sum(est_conv_flag[ind,:])/sum(est_flag[ind,:])])
end

##















 ##
