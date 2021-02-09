N = 10
NG = N#5

fooNetModelDirBin1
deltaN = 3

tmpdegsIO = round.([ Vector(linspace(deltaN,N-1-deltaN,N));Vector(linspace(deltaN,N-1-deltaN,N))])
tmpdegsIO  =   DegSeq2graphDegSeq(fooNetModelDirBin1,tmpdegsIO)
linSpacedPar(fooNetModelDirBin1,N; Ngroups = NG, deltaN=deltaN)
## test graphicability functions


##
groupsInds = distributeAinVecN(Vector(1:NG),N)
Mod = NetModelDirBin1(deg,zeros(N),groupsInds)
groupsNames = unique(groupsInds)
NG = length(groupsNames)
tmpDic = countmap(groupsInds)
#number of nodes per group
NGMemb = sortrows([[i for i in keys(tmpDic)] [  tmpDic[i] for i in keys(tmpDic) ] ])
# agrregate degree of each group
degGroups = zeros(2NG) ; for i =1:N degGroups[groupsInds[i]] += deg[i]; degGroups[NG + groupsInds[i]] += deg[N + i];  end
L = sum(degGroups,1)[1]
L/(N^2-N)
groupsPar = ones(2NG)
firstOrderCond(Mod;degIO = deg,parGroupsIO = groupsPar, groupsInds = groupsInds )
expMatrix(Mod,groupsPar)

targetErrValStaticNets = 0.005
estPar,estIt,estMod = estimate(Mod)

##
sam = sampl(estMod,1000)

meanMatSam = squeeze(mean(sam,3),3)
degsSam = expValStatsFromMat(estMod,meanMatSam)
degGroupsSam = zeros(2NG) ; for i =1:N degGroupsSam[groupsInds[i]] += degsSam[i]; degGroupsSam[NG + groupsInds[i]] += degsSam[ N + i] end
relErrSam = maximum(abs.(( degGroupsSam - degGroups )./degGroups )[degGroups.!=0])

##
















#
