
# The purpose of this script is to compare the filtering capability of the single snapshots
# estimates and the gas filter on a given

using StatsBase, JLD,  Hwloc,  DynNets, StaticNets, StatsFuns,HelperFunDom

# static Estimates for increasing N ---------------------------------------------

#need to choose how to distribute the unconditional means. I want to stay
# in the case of chaterjee Diaconis Sly, hence have fixed maximum value for the
# parameters that does not increase with N
fractionParTV = 0.05;

Nrange = 50:100:800
Nsample = 50
degMin = 5
degMax = 30

midInd(vec::Array{<:Real,1}) = round(Int,length(vec)/2)

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


##
