using StatsBase, StatsFuns, JLD,  Hwloc

counts = Hwloc.histmap(Hwloc.topology_load())
ncores = counts[:Core]
if nworkers()<ncores
    addprocs(ncores   - nprocs())
end

using DynNets
@everywhere using DynNets
@everywhere using StatsFuns
N_est = 900
Tvals = [50 500 5000]#500#
N_ind = length(Tvals)[1]


N = 10 ;
GBAest =1; GWest = N; useStartVal = false

groupsIndsEst = DynNets.distributeAinVecN(ones(Int64,1),N)




est_times = Array{Float64}(zeros(N_est))
est_conv_flag = falses(N_est)
est_flag = falses(N_est)

tmpdegsT = ones(10,N)

#ridefinizioni che servono per allocare
tmpModelEst =  DynNets.GasNetModelBin1(tmpdegsT,[zeros(GWest),zeros(GBAest),zeros(GBAest)],
                                    groupsIndsEst)
# tmpModelDgp =  DynNets.GasNetModelBin1_UmBAgr(tmpdegsT,[zeros(GWdgp),zeros(GBAdgp),zeros(GBAdgp)],
#                                     groupsIndsDgp)
# tmpModelDgp.Par[1] = ones(GWdgp);tmpModelDgp.Par[2]  = ones(GBAdgp);tmpModelDgp.Par[3]  = ones(GBAdgp)
size(tmpModelEst.Par[1])[1]
# Set the value of the parameters common to all replicate of DGP+Estimation
sampleDegsT = 0; n=0;fit_T = 0;
#Initialize the variable to store
for i =1: size(tmpModelEst.Par)[1]
    ex1 = "estPar_$(i)  = SharedArray{Float64}(zeros(N_ind,N_est, size(tmpModelEst.Par[$(i)])[1] ))"
    eval(parse(ex1))
 end

ObsAllT = Array{Array{Int16,3},1}(N_ind)
ParTvDgpAllT = Array{Array{Float64,3},1}(N_ind)

#define save path and list variables to save
save_fold = "./data/estimatesTest/"
    file_name = "EstTest_$(N)_$(Tvals)_$(N_est)_$(useStartVal)_Sin_Miss_Targeting_GBA$(GBAest).jld"
    save_path = save_fold*file_name#



est_times = SharedArray{Float64}(zeros(N_ind,N_est))
est_conv_flag = SharedArray{Bool,2}(N_ind,N_est)
est_flag = SharedArray{Bool,2}(N_ind,N_est)

@progress "T"  for ind = 1:N_ind
#   ind=1

    T = Tvals[ind]
    ObsTmp = SharedArray{Int16,3}(Tvals[ind],N,N_est)
    ParTvDgpTmp = SharedArray{Float64,3}(Tvals[ind],N,N_est)

    @sync @parallel   for n = 1:N_est
#     n=1


        Y_T,fit_T  = DynNets.sampleDynNetsBin1(T,"FIT-MULTI-SIN";N=N)
        sampleDegsT = squeeze(sum(Y_T,2),2)

        ParTvDgpTmp[:,:,n] = fit_T
        ObsTmp[:,:,n] = sampleDegsT

        Model2Est =  DynNets.GasNetModelBin1(sampleDegsT,[N,zeros(GBAest),zeros(GBAest)],
                                                    groupsIndsEst)
        #estimate
        useStartVal ?   start_val = [sim_W, sim_B, sim_A] : start_val = zeros(3,3)
        est_times[ind,n] = @elapsed  hat_pars , conv_flag =
                            DynNets.estimateTargeting(Model2Est )

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
