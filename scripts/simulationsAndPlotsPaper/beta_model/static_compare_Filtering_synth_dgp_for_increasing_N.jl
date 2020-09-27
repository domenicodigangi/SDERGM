
# The purpose of this script is to compare the filtering capability of the single snapshots
# estimates and the gas filter on a given

using StatsBase, JLD,  Hwloc,  DynNets, StaNets, StatsFuns,HelperFunDom

# STATIC ESTIMATES INCREASING N  ---------------------------------------------
#need to choose how to distribute the unconditional means. I want to stay
# in the case of chaterjee Diaconis Sly, hence have fixed maximum value for the
# parameters that does not increase with N
fractionParTV = 0.05;

Nrange = 50:50:100#100:150#800
Nsample = 5
degMin = 30
degMax = 30
function simulAndEstVarN(model::NetModelDirBin1,Nrange::StepRange{Int,Int},Nsample::Int,degMin::Int,degMax::Int)
    density=zeros(length(Nrange),3)
    dgpFits = Array{Array{Float64,1},2}(length(Nrange),3 )
    #generate multiple versions of the degree sequences for the dgp
     for n = 1:length(Nrange)
        println(n)
        Nn = Nrange[n]
        fits,degs =  dgpFitVarN(Nn,degMin,degMax;exponent = 1)
        density[n,1] = sum(degs[1:Nn])/(Nn^2-Nn)
        dgpFits[n,1] = fits
        fits,degs =  dgpFitVarN(Nn,degMin,Nn-degMin;exponent = 0.4 )
        density[n,2] = sum(degs[1:Nn])/(Nn^2-Nn)
        dgpFits[n,2] = fits
        fits,degs =  dgpFitVarN(Nn,degMin,Nn-degMin;exponent = 1 )
        density[n,3] = sum(degs[1:Nn])/(Nn^2-Nn)
        dgpFits[n,3] = fits
    end
    estFits = Array{Array{Float64,2},2}(length(Nrange),3 )
    @time for n = 1:length(Nrange)

        Nn = Nrange[n]
        for i=1:size(dgpFits)[2]
            # sample the dgpFits
            matSamp =  StaNets.sampl(NetModelDirBin1(dgpFits[n,i]),Nsample;  parGroupsIO = dgpFits[n,i])
            degsSamp = [sumSq(matSamp,2);sumSq(matSamp,1)]
            #println(degsSamp   )
            tmpEst = zeros(2Nn,Nsample)
            for s=1:Nsample
                tmpEst[:,s],~ = StaNets.estimate(fooNetModelDirBin1,degIO = degsSamp[:,s])

            end
            estFits[n,i] = tmpEst
        end
            println(n)
    end
 return dgpFits,estFits
end

function simulAndEstVarN(model::NetModelDirBin1,dgpDegs::Array{<:Real,1},Nsample::Int)
    # if the degree sequence is given as input use that
    dgpFits,~,~ = meanFit ,~ = StaNets.estimate(StaNets.NetModelDirBin1(dgpDegs),degIO = dgpDegs  )
    N2 = length(dgpDegs)
    estFits = Array{Float64,2}(N2, Nsample )
    estFits = zeros(N2,Nsample)
    matSamp = StaNets.sampl(NetModelDirBin1(dgpFits),Nsample;  parGroupsIO = dgpFits)
    degsSamp = [sumSq(matSamp,2);sumSq(matSamp,1)]
        for s=1:Nsample
            estFits[:,s],~ = StaNets.estimate(fooNetModelDirBin1,degIO = degsSamp[:,s])
            println(s)
        end
 return dgpFits,estFits
end

# sample and estimate networks of increasing size AND STORE DATA
#@time dgpFits,estFits = simulAndEstVarN(fooNetModelDirBin1,Nrange,Nsample,degMin,degMax)

# Load data and plot results
save_fold = "./data/estimatesTest/asympTest/"

data = load(save_fold*"asympTestStaticLarge.jld")
estFits = data["estFits"]
dgpFits = data["dgpFits"]
N,Scal = size(dgpFits)
Sam = size(estFits[1,1])[2]
expDegsDgp = copy(dgpFits)
for n=1:N,scal=1:Scal expDegsDgp[n,scal] = expValStats(fooNetModelDirBin1,dgpFits[n,scal]) end
expDegsEst = copy(estFits)
for n=1:N,scal=1:Scal,sam =1:Sam  expDegsEst[n,scal][:,sam] = expValStats(fooNetModelDirBin1,estFits[n,scal][:,sam]) end
# plot the Degrees in the dgp
for n=1:N
    scal = 3
    plot(expDegsDgp[n,scal][1:midInd(expDegsDgp[n,scal])])
 end
 grid()
# plot densities
for scal=[1,3]
     dens = [sum(expDegsDgp[n,scal])/(length(expDegsDgp[n,scal])*(length(expDegsDgp[n,scal])-1)) for n=1:N]
     size = [length(expDegsDgp[n,scal]) for n=1:N]
     plot(size,dens ,".-")
  end
  grid()
  xlabel("number of nodes")
  ylabel("density")
#Plot mean relative bias over ensembles
for n=1:N
    scal = 3

 diffDegs = (expDegsDgp[n,scal] .- expDegsEst[n,scal])./expDegsDgp[n,scal]
  sampMenDiff = meanSq(diffDegs,2)
 plot(expDegsDgp[n,scal][1:midInd(expDegsDgp[n,scal])],n*0.2 + sampMenDiff[1:midInd(sampMenDiff)],".")
 end
 yticks(0.2:0.2:N*0.2)
 grid()
for n=1:N
     scal = 3
  diffDegs = (expDegsDgp[n,scal] .- expDegsEst[n,scal])./expDegsDgp[n,scal]
  sampStdDiff = squeeze(std(diffDegs,2),2)
  plot(expDegsDgp[n,scal][1:midInd(sampStdDiff)],n*0.2 + sampStdDiff[1:midInd(sampStdDiff)],".")
 end
  yticks(0.2:0.2:N*0.2)
  grid()
# load eMid data and run the estimates tests on realistically spaced fitnesses
LoadData =true
 halfPeriod=false
 fold_Path =  "/home/Domenico/Dropbox/Dynamic_Networks/data/emid_data/juliaFiles/"
 loadFilePartialName = "Weekly_eMid_Data_from_"
 halfPeriod? periodEndStr =  "2012_03_12_to_2015_02_27.jld": periodEndStr =  "2009_06_22_to_2015_02_27.jld"
 @load(fold_Path*loadFilePartialName*periodEndStr, AeMidWeekly_T,banksIDs,inactiveBanks, YeMidWeekly_T,weekInd,datesONeMid ,degsIO_T,strIO_T)
Sam = 500
degsIO_T = [sumSq(AeMidWeekly_T,2);sumSq(AeMidWeekly_T,1)]
#plot(sumSq(degsIO_T.>0,1)
meanDegsEmid = meanSq(degsIO_T,2)
dgpDegs =  round.(Int,meanDegsEmid)
dgpFits,estFits = simulAndEstVarN(fooNetModelDirBin1,dgpDegs,Sam)
#plot empirical cumulative distribution
dgpDegDistr = ecdf(dgpDegs)
 plot(0:0.1:maximum(dgpDegs),dgpDegDistr(0:0.1:maximum(dgpDegs)),".-")
 grid()
 xticks(0:maximum(dgpDegs))
 xlabel("Degree")
 ylabel("emprical distribution")
expDegsDgp = expValStats(fooNetModelDirBin1,dgpFits)
expDegsEst = copy(estFits)
for sam =1:Sam  expDegsEst[:,sam] = expValStats(fooNetModelDirBin1,estFits[:,sam]) end
diff =  (expDegsDgp .- expDegsEst)
diffRel = diff./expDegsDgp
 diffRel[expDegsDgp.==0,:] = diff[expDegsDgp.==0,:]
# plot mean and std of relative error degrees estimates over the ensemble
fig, ax = subplots(2,1)
 ax[1,1][:plot](expDegsDgp,mean(diffRel,2),".")
 xticks(1:ceil(maximum(expDegsDgp)))
 ax[2,1][:plot](expDegsDgp,std(diffRel,2),".")
 ax[1,1][:grid]("on")
 ax[2,1][:grid]("on")
 ax[1,1][:set_xticks](0:ceil(maximum(expDegsDgp)))
 ax[1,1][:set_ylabel]("mean of relative error degs")
 ax[2,1][:set_xticks](0:ceil(maximum(expDegsDgp)))
 ax[2,1][:set_ylabel]("std of relative error degs")

indPlot = 7
plt[:hist](diffRel[indPlot,:],bins = 10,normed =true);
 grid()
 title("deg = $(round(Int,expDegsDgp[indPlot]))")
plt[:hist](diffRel[indmax(expDegsDgp),:],bins = 18,normed =true);
 grid()
 title("deg = $(round(Int,maximum(expDegsDgp)))")

 #
