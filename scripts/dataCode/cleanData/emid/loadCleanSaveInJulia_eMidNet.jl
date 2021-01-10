using StatsBase , JLD,  Dates, PyPlot, DelimitedFiles
        using  Utilities

        # script that loads the csv file with emid data and cleans it:
        #
        # --- aggregate the ON loans in weekly adjacency matrices
        # --- selecti only banks that are active at least a given fraction of the weeks
        # --- Store also a half period version of the data, starting after LTRO began


        halfPeriod = false # Consider only data after LTRO began ?
        ####--- Load Data from txt files
        ##load Dates
        loadPath = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/csvFiles/banksNames.txt"
        banksIDs = Array{Int}(readdlm(loadPath,','))
        NallBanks = length(unique(banksIDs))
        N = NallBanks
        loadPath = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/csvFiles/dates.txt"
        tmp = Array{Int}(readdlm(loadPath,','))
        alldatesONeMid = Array{Date,1}(undef,length(tmp) ) # vector of the dates in the data
        for i=1:length(tmp) alldatesONeMid[i] = Date("$(tmp[i])","yyyymmdd") end # convert it in dates type
        TdaysAll = length(unique(alldatesONeMid))
        halfPeriod ? startDate = Date("2012-03-12") : startDate=alldatesONeMid[1] #
        endDate = Date("2015-02-27")
        Dates.dayname(startDate) # The first day in the data is a
        timesToKeep = (alldatesONeMid .>=startDate) .&(alldatesONeMid .<=endDate)
        Tdays = sum(timesToKeep)
        datesONeMid = alldatesONeMid[timesToKeep]
        ## Load Matrices
        loadPath = "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/csvFiles/eMid_all_T.txt"

tmp = readdlm(loadPath,',')
YeMidDaylyAll_T = reshape(tmp,NallBanks,NallBanks,TdaysAll)
YeMidDaylyAll_T = YeMidDaylyAll_T[:,:,timesToKeep]
## Aggregate Matrices weekly
fracWeek =  [tmp.value for tmp in (datesONeMid .- datesONeMid[1])]/7# to wich week each day belongs to ?
[Dates.dayname.(datesONeMid) fracWeek floor.(fracWeek)]
weekInd = 1 .+ floor.(Int,fracWeek)
weekInd[end] == length(unique(weekInd)) ? () : error()
Tweeks = weekInd[end]
YeMidWeeklyAll_T = zeros(NallBanks,NallBanks,Tweeks)
for t = 1:Tweeks YeMidWeeklyAll_T[:,:,t] = sum(YeMidDaylyAll_T[:,:,weekInd .== t],dims=3) end
date_first_day_week = [datesONeMid[weekInd .== t][1] for t =1:Tweeks  ]# t = 1:Tweeks  = sum(YeMidDaylyAll_T[:,:,weekInd .== t],dims=3) end
days = day.(date_first_day_week)
months = month.(date_first_day_week)
years = year.(date_first_day_week)

#Some banks' statistics:
AeMidWeeklyAll_T = YeMidWeeklyAll_T[:,:,:] .> 0
degsI_T = sumSq(AeMidWeeklyAll_T,2)
degsO_T = sumSq(AeMidWeeklyAll_T,1)
strI_T = sumSq(YeMidWeeklyAll_T,2)
strO_T = sumSq(YeMidWeeklyAll_T,1)
# Select only banks active for more than a fixed proportion of the times
activityThreshold = 0.05# 0.9999
isBankiActiveAtTimet =  (degsI_T.>0) .| (degsO_T.>0)
activityBanks = sumSq(isBankiActiveAtTimet,2)./(length(isBankiActiveAtTimet[1,:]))
inactiveBanks = activityBanks .<= activityThreshold; NinactiveBanks = sum(inactiveBanks)
N = NallBanks - NinactiveBanks

##

#compute the strengths and degrees for the sub network of active banks!!

YeMidWeekly_T = YeMidWeeklyAll_T[.!inactiveBanks,.!inactiveBanks,:]
AeMidWeekly_T = YeMidWeekly_T[:,:,:] .> 0
degsI_T = sumSq(AeMidWeekly_T,2)
degsO_T = sumSq(AeMidWeekly_T,1)
strI_T = sumSq(YeMidWeekly_T,2)
strO_T = sumSq(YeMidWeekly_T,1)
degsIO_T = [degsI_T; degsO_T]
strIO_T = [strI_T; strO_T]
~,meandegsI,meandegsO = splitVec(meanSq(degsIO_T,2))

sum(strI_T)-sum(strO_T) >1e-12 ? error() : ()
save_fold =   "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
file_nameStart = "Weekly_eMid_Data"
file_nameEnd = "_from_$(Dates.format(startDate,"yyyy_mm_dd"))"*
        "_to_" *"$(Dates.format(endDate,"yyyy_mm_dd"))"*".jld"

# save Data
save_path = save_fold*file_nameStart  *file_nameEnd#
@save(save_path, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,
datesONeMid ,degsIO_T,strIO_T, years, months, days)
