using StatsBase, JLD,  Hwloc,  DynNets, StatsFuns,Utilities

#initWorkers()
#@everywhere  using DynNets, StatsFuns, Utilities
N_est = 2
NGWpoints = 1
Tvals =  [  600 400 200 ]#[50 250 1250]#
N_Tind = length(Tvals)[1]
scoreRescType = ""#"FISHER-EWMA"#
useStartVal = false# true# false
N_Nmin = 500
Nspacing =  500
Nmax = 10000
Nvals =  Vector(Nmin:Nspacing:Nmax)
Nind = length(Nvals)

est_times = SharedArray{Float64}(zeros(N_Tind,N_Nind,NGWpoints,N_est)) #zeros(N_Tind,N_Nind,NGWpoints,N_est)#
est_conv_flag =  SharedArray{Bool,4}(N_Tind,N_Nind,NGWpoints,N_est) #falses(N_Tind,N_Nind,NGWpoints,N_est) #
est_flag = SharedArray{Bool,4}(N_Tind,N_Nind,NGWpoints,N_est) #falses(N_Tind,N_Nind,NGWpoints,N_est) #
GWvalStore  = Array{Array{Int,1},1}(N_Nind)


@everywhere function loopEval(T,N,simParVec, useStartVal,groupsindTsDgp,scoreRescType,GBAest,GWest,groupsindTsEst)
    Model2Sample = DynNets.GasNetModelBin1(ones(T,N),simParVec,
                                                groupsindTsDgp,scoreRescType)
    Model2Sample,Y_T,Fitness_T  =  DynNets.sampl(Model2Sample,T)
    sampleDegsT = squeeze(sum(Y_T,2),2)

    Model2Est =  DynNets.GasNetModelBin1(sampleDegsT,[zeros(GWest),zeros(GBAest),zeros(GBAest)],
                                                groupsindTsEst,scoreRescType)
    #estimate
    useStartVal ?   start_val = simParVec : start_val = zeros(3,3)
    estTime = @elapsed  hat_pars , conv_flag =
                        DynNets.estimate(Model2Est,start_values = start_val )
    return estTime,conv_flag
end



save_fold = "./data/estimatesTest/compComplexity/"
file_name = "EstTest_varN_Nest_$(N_est)_$(Nmin)_$(Nspacing)_$(Nmax)_varGW_$(NGWpoints)_T_$(Tvals)_corectSpec"* scoreRescType *".jld"
save_path = save_fold*file_name#

#

GW = 1
@progress "N" for indN = 1:N_Nind
#indN = 1
N = Nvals[indN]
#GWvals = unique(Vector{Int}(floor.(linspace(1,N,NGWpoints))))
#GWvalStore[indN] = GWvals

#@progress "GW"  for indGW = NGWpoints:-1:1
indGW = 1
#GW =  GWvals[indGW]
GBAdgp =1; GWdgp = GW ;
GBAest =1; GWest = GW;
#Distribuisci i nodi tra i gruppi in parti uguali
groupsindTsDgp = [[distributeAinVecN(Array{Int}(1:GWdgp),N)];
                    [distributeAinVecN(Array{Int}(1:GBAdgp),GWdgp)]]
groupsindTsEst = [[distributeAinVecN(Array{Int}(1:GWest),N)];
                   [distributeAinVecN(Array{Int}(1:GBAest),GWest)]]


#Initialize the variable to store
realUm = DynNets.linSpacedPar(DynNets.fooGasNetModelBin1, N;NgroupsW = GW )[1]
sim_B = Array(linspace(0.4,0.9,GBAdgp+1)[2:end]) #
sim_A = Array(linspace(0.02,0.2,GBAdgp+1)[2:end]) #
sim_W = realUm.*(1-sim_B[groupsindTsDgp[2]])
simParVec = [sim_W,sim_B,sim_A]



@progress "T" for indT = 1:N_Tind
#  indT=1
T = Tvals[indT]
#@sync @parallel
for n = 1:N_est
#n=1
est_times[indT,indN,indGW,n] , est_conv_flag[indT,indN,indGW,n]  = loopEval(T,N,simParVec, useStartVal,groupsindTsDgp,scoreRescType,GBAest,GWest,groupsindTsEst)
est_flag[indT,indN,indGW,n] = true
println(11)
#end
# save
end
@save(save_path,est_conv_flag,est_times,est_flag,Nvals,GWvalStore ,useStartVal)
println([N,GW,T,sum(est_times[indT,indN,indGW,:]),sum(est_conv_flag[indT,indN,indGW,:])/sum(est_flag[indT,indN,indGW,:])])

end

end

##














 ##
