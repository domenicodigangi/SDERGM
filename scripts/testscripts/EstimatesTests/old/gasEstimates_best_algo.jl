
# Sample a gas fitness DGP and estimate the parameters for different lenghts of
# time series
using StatsFuns
using StatsBase
using JLD
N_est = 1

N = 50 ; Ngroups = 1; GTV =Ngroups; useStartVal = false
NnodesPerGroup = floor(N/Ngroups);groupsInds=zeros(Int64,N)
#Distribuisci i nodi tra i gruppi in parti uguali
for i =0:Ngroups-2 groupsInds[1 + Int64(NnodesPerGroup*i):Int64(NnodesPerGroup*(i+1))] = i+1 end;
groupsInds[find(groupsInds .==0)] = Ngroups
indTvGroups = falses(Ngroups);indTvGroups[1:GTV] = true
T = 200

#the values for the unconidtional means to be used in the DGP are obtained estimating
#an equally spaced degree sequence
tmpDegs = Array{Real,1}(linspace(2,N-2,N) )'; tmpdegsT = repmat(tmpDegs,2,1)  #ones(N)*N/2
tmpModel =  DynNets.GasNetModelBin1(tmpdegsT,[zeros(N),zeros(GTV),zeros(GTV)],
                                    ones(Int,N))
est_times =  Array{Float64}(zeros( N_est))
est_conv_flag =  falses(N_est)
est_flag =  falses(N_est)

# Set the value of the parameters common to all replicate of DGP+Estimation
tmp,~ = DynNets.estimate(tmpModel )
realUm = tmp[1]
sim_B = Array(linspace(0.4,0.9,GTV+1)[2:end]) #
sim_A = Array(linspace(0.02,0.2,GTV+1)[2:end]) #
sim_W = realUm.*(1-sim_B[groupsInds])
 #  for n = 1:N_est
n = 1
        Y_T,fit_T  =  DynNets.sampleDynNetsBin1(T,"Fit-Multi-gas";N=N,
                                                    groupsInds=groupsInds,
                                                    AgasGroups =sim_A,
                                                    BgasGroups = sim_B,
                                                    WgasNodes = sim_W )
        degsT = squeeze(sum(Y_T,2),2)
        model_data = DynNets.GasNetModelBin1(   degsT,
                                                [zeros(N),zeros(Ngroups),zeros(Ngroups)],
                                                groupsInds)
        #estimate
        useStartVal ?   start_val = [sim_W, sim_B, sim_A] : start_val = zeros(3,3)
        est_times[n] = @elapsed  hat_pars , conv_flag =
                            DynNets.estimate(model_data,start_values = start_val )
        indTvNodes =   .!(groupsInds .==0);
        est_conv_flag[n] = conv_flag
        est_flag[n] = true
#    end
println([N,T,sum(est_times),sum(est_conv_flag)/sum(est_flag)])




 ##
