
using Utilities

@time for i=1:1000
    mat = rand(10,10).>0.5
     mat = putZeroDiag(mat)
     f = transitivityDirBin_funErgm(mat)
     i,j = (2,1)
     delta = transitivityDirBin_funErgm_delta(mat,i,j)
     mat[i,j] = false ; f1 =  transitivityDirBin_funErgm(mat);
     mat[i,j] = true ; f2 =  transitivityDirBin_funErgm(mat);
     (f2 - f1) != delta ? error():()
end




function pseudoLogLikelihood(θ::Real,A::BitArray{2},statF::statErgmFun)
     N1,N2 = size(A)
     statVal = statF.f(A) #the value of the statistics on the real networks
     # compute a matrix of changes in the statistics when each element is toggled from 0 to 1
     # come lo faccio in modo efficiente??
     deltaFun(l,m) = statF.delta(A,l,m)
     logPsLike = 0
     for i=1:N1,j=1:N2
         if i!=j
            π_ij = 1/(1+ exp(- θ * deltaFun(i,j) ))
            if A[i,j] #is non zero
                logPsLike = logPsLike + log(π_ij)
            else
                logPsLike = logPsLike + log( 1 - π_ij )
            end
        end
     end
     return logPsLike
 end

testStat = statErgmFun(transitivityDirBin_funErgm,transitivityDirBin_funErgm_delta)

mat = rand(10,10).>0.5
pseudoLogLikelihood(0.0,mat,testStat)
