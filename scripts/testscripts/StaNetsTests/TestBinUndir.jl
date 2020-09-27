N = 10
Ngroups = N
groupsInds = distributeAinVecN(Vector(1:Ngroups),N)
tmpDic = countmap(groupsInds)
#number of nodes per group
NgroupsMemb = sortrows([[i for i in keys(tmpDic)] [  tmpDic[i] for i in keys(tmpDic) ] ])

deltaN = 3
groupsPar =  ones(Ngroups)
degGroups =  DegSeq2graphDegSeq(fooNetModelBin1, Array{Real,1}(linspace(1+deltaN,N - round(0.1*N)-deltaN,N)))
#tmpdegs = [8, 4, 7, 8, 5, 5, 6, 5, 6, 0]
#sort([7, 5, 8, 9, 8, 5, 4, 2, 8, 4],rev = true)# DegSeq2graphDegSeq(fooNetModelBin1, sort([9, 8, 4, 7, 3, 6, 4, 7, 4, 2],rev = true))
Mod = NetModelBin1(degGroups,groupsPar,groupsInds)
groupsPar, nIterm,estMod = estimate(Mod)
firstOrderCond(Mod;deg = tmpdegs,par = groupsPar, groupsInds = groupsInds )
expMatrix(Mod,groupsPar[groupsInds])

## test sampling
groupsPar[10] -= 10
Nsample = 1000
sam = sampl(estMod,Nsample)
meanMatSam = squeeze(mean(sam,3),3)
degsSam = squeeze(sum(meanMatSam,2),2)
degGroupsSam = zeros(Ngroups) ; for i =1:N degGroupsSam[groupsInds[i]] += degsSam[i] end
relErrSam = maximum(abs.(( degGroupsSam - degGroups )./degGroups ))
## Test estimate umbiasedness
estParSam = zeros(N,Nsample)
@time for n=1:Nsample
    degs_n = expValStatsFromMat(fooNetModelBin1,sam[:,:,n] )
    estParSam[:,n], ~,~ = estimate(Mod; deg= degs_n)
end

##
using Plots, AverageShiftedHistograms, Distributions, StatPlots
plotly()
listPlots=""
    for n=1:N
        err = ((estParSam .- groupsPar))[n,:]
        rng_space = -3:.1:3
        n==10 ? (rng_space = -300:.1:300) : ()
        tmpAsh = ash(err;rng = rng_space)
        pA1 = plot(tmpAsh)
        expr = "PW$(n) = plot!(Normal(0,1),linewidth = 5) "
        eval(parse(expr))
        listPlots *= "PW$(n) , "
    end

exprLast = "plot(" * listPlots  *" titlefont = font(12) ,titleloc=:center,size=(1350,600),legend = :none )"
 eval(parse(exprLast))


##








#
