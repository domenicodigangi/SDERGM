using  JLD,Dates
## Load Data
startDate = Date("2012-03-12")
endDate = Date("2015-02-27")
save_fold =   "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
file_nameStart = "Weekly_eMid_Estimates"
file_nameEnd = "_from_$(Dates.format(startDate,"yyyy_mm_dd"))"*
        "_to_" *"$(Dates.format(endDate,"yyyy_mm_dd"))"*".jld"

save_path = save_fold*file_nameStart* "_W1_" *file_nameEnd#
@load(save_path,estSSW1,estGasW1,estGasTargW1,estSSW1_1GW,estGasW1_1GW,AeMidWeekly_T, YeMidWeekly_T ,weekInd,degsIO_T)

# test netmodel direceted binary with covariates on eMid Data
using  StaticNets ,Utilities


AeMidWeekly_T = BitArray(YeMidWeekly_T )
AmatsEMid = AeMidWeekly_T[:,:,2:end]
#use the presence or absence of a link as regressor
regs0eMid = AeMidWeekly_T[:,:,1:end-1]
#regs0eMid = YeMidWeekly_T[:,:,1:end-1]#use the weight of previous time link





























##
