
# Sample a gas fitness DGP and estimate the parameters for different lenghts of
# time series

using StatsBase
using JLD

initWorkers()

using DynNets
using StatsFuns
@everywhere using DynNets
@everywhere using StatsFuns
N_est = 300

N = 10 ; Ngroups = 1; GTV =Ngroups; useStartVal = true
NnodesPerGroup = floor(N/Ngroups);groupsInds=zeros(Int64,N)
#Distribuisci i nodi tra i gruppi in parti uguali
for i =0:Ngroups-2 groupsInds[1 + Int64(NnodesPerGroup*i):Int64(NnodesPerGroup*(i+1))] = i+1 end;
groupsInds[find(groupsInds .==0)] = Ngroups
indTvGroups = falses(Ngroups);indTvGroups[1:GTV] = true
Tvals = [1000]
N_ind = length(Tvals)[1]

#the values for the unconidtional means to be used in the DGP are obtained estimating
#an equally spaced degree sequence
tmpDegs = Array{Real,1}(linspace(2,N-2,N))'; tmpdegsT = repmat(tmpDegs,2,1)#ones(N)*N/2
tmpModel =  DynNets.GasNetModelBin1(tmpdegsT,[zeros(N),zeros(GTV),zeros(GTV)],
                                    zeros(Int,N),
                                    falses(1))

#

#Initialize the variable to store
for i =1: size(tmpModel.Par)[1]
    ex1 = "estPar_$(i)  = SharedArray{Float64}(zeros(N_ind,N_est, size(tmpModel.Par[$(i)])[1] ))"
    ex2 = "simPar_$(i) = SharedArray{Float64}(zeros(N_ind,N_est, size(tmpModel.Par[$(i)])[1] ))"
    eval(parse(ex1))
    eval(parse(ex2))
 end
 est_times = SharedArray{Float64}(zeros(N_ind,N_est))
 est_conv_flag = SharedArray{Bool,2}(N_ind,N_est)
 est_flag = SharedArray{Bool,2}(N_ind,N_est)

# Set the value of the parameters common to all replicate of DGP+Estimation
tmp,~ = DynNets.estimate(tmpModel )
realUm = tmp[1]
sim_B = Array(linspace(0.4,0.9,GTV+1)[2:end]) #
sim_A = Array(linspace(0.02,0.2,GTV+1)[2:end]) #
sim_W = realUm.*(1-sim_B)
@progress "T"  for ind = 1:N_ind
 # ind=1
   T = Tvals[ind]
    @sync @parallel   for n = 1:N_est
  #n=1

        simPar_3[ind,n,1:GTV] = sim_A
        simPar_2[ind,n,1:GTV] = sim_B
        simPar_1[ind,n,1:N] = sim_W
        Y_T,fit_T  =  DynNets.sampleDynNetsBin1(T,"Fit-Multi-gas";N=N,
                                                    indTvGroups=indTvGroups,
                                                    groupsInds=groupsInds,
                                                    AgasGroups =sim_A,
                                                    BgasGroups = sim_B,
                                                    WgasNodes = sim_W )
        degsT = squeeze(sum(Y_T,2),2)
        model_data = DynNets.GasNetModelBin1(   degsT,
                                                [zeros(N),zeros(Ngroups),zeros(Ngroups)],
                                                groupsInds,
                                                indTvGroups)
        #estimate
        useStartVal ?   start_val = [sim_W, sim_B, sim_A] : start_val = zeros(3,3)
        est_times[ind,n] = @elapsed  hat_pars , conv_flag =
                            DynNets.estimate(model_data,start_values = start_val )

        indTvNodes =   .!(groupsInds .==0);
        estPar_1[ind,n,:] = hat_pars[1]
        estPar_2[ind,n,:] = hat_pars[2]
        estPar_3[ind,n,:] = hat_pars[3]
        est_conv_flag[ind,n] = conv_flag
        est_flag[ind,n] = true

    end

    # save
    save_fold = "./data/estimatesTest/"
    file_name = "EstTest_$(N)_$(GTV)_$(Tvals)_$(N_est)_$(useStartVal).jld"
    save_path = save_fold*file_name#
        @save(save_path,groupsInds,indTvGroups,simPar_1,simPar_2,simPar_3,estPar_1,estPar_2,estPar_3,est_conv_flag,est_times)
    println([T,sum(est_times[ind,:]),sum(est_conv_flag[ind,:])/sum(est_flag[ind,:])])
end

##


 ##
