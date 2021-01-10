
using Plots, JLD, Utilities
plotly()

N_est = 5
NGWpoints = 3
Tvals =  [  600 400 200 ]#[50 250 1250]#
N_Tind = length(Tvals)[1]
scoreRescType = ""#"FISHER-EWMA"#
useStartVal = false# true# false
save_fold = "./data/estimatesTest/compComplexity/"
##




Nmin = 10; Nspacing =  10; Nmax = 70
file_name = "EstTest_varN_Nest_$(N_est)_$(Nmin)_$(Nspacing)_$(Nmax)_varGW_$(NGWpoints)_T_$(Tvals)_corectSpec"* scoreRescType *".jld"; save_path = save_fold*file_name#
@load(save_path,est_conv_flag,est_times,est_flag,Nvals,GWvalStore ,useStartVal)

NvalsAll = Nvals
est_timesAll = est_times
##


Nmin = 80; Nspacing =  10; Nmax = 100
file_name = "EstTest_varN_Nest_$(N_est)_$(Nmin)_$(Nspacing)_$(Nmax)_varGW_$(NGWpoints)_T_$(Tvals)_corectSpec"* scoreRescType *".jld"; save_path = save_fold*file_name#
@load(save_path,est_conv_flag,est_times,est_flag,Nvals,GWvalStore ,useStartVal)

NindSelect = 1
NvalsAll = [Nvals[1];NvalsAll ]
est_timesAll = cat(2, est_times[1:end,NindSelect:NindSelect,1:end,1:end],est_timesAll)

##

GWvalsLeg = ["1 par" ; "N/2 par"; " N par"]

est_times[1,1,1:end,1:end]
#est_times[indT,indN,indGW,n]
##
h = plot(); indGW = 2
 for t = 1:3
    est_timesAll[t,1:end,indGW,1:end]
    plot!(h,NvalsAll,meanSq(est_timesAll[t,1:end,indGW,1:end],2),label = "$(Tvals[t])")
 end
 plot(h, title = GWvalsLeg[indGW])

## add more single par estimates
Nmin = 500; Nspacing =  500; Nmax = 10000; NGWpoints = 1;N_est = 1
file_name = "EstTest_varN_Nest_$(N_est)_$(Nmin)_$(Nspacing)_$(Nmax)_varGW_$(NGWpoints)_T_$(Tvals)_corectSpec"* scoreRescType *".jld"; save_path = save_fold*file_name#
save_path
@load(save_path,est_conv_flag,est_times,est_flag,Nvals,GWvalStore ,useStartVal)
tmp = est_times[:,:,1,1]


##
NvalsAll_GW1 = [NvalsAll;Nvals]
est_timesAll_1GW = [meanSq(est_timesAll[:,:,1,:],3) tmp]



h = plot(); indGW = 1
 nnz = est_timesAll_1GW[1,:].!=0
 for t = 1:3
    est_timesAll_1GW[t,nnz]
    plot!(h,NvalsAll_GW1[nnz],est_timesAll_1GW[t,nnz],label = "$(Tvals[t])")
 end
 plot(h, title = GWvalsLeg[indGW])


#




















##
