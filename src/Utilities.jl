module Utilities

using Logging
using MLBase
using LinearAlgebra
using PyCall; pygui(:qt); using PyPlot

using ..AReg


logit(x) = log(x/(1-x))
export logit


inv_logit(x) = 1/(1+exp(-x))
export inv_logit

link_R_in_0_1(x) = inv_logit(x)
export link_R_in_0_1

inv_link_R_in_0_1(x) = logit(x)
export inv_link_R_in_0_1

link_R_in_R_pos(x) = exp(x) 
export link_R_in_R_pos

inv_link_R_in_R_pos(x) = log(x) 
export inv_link_R_in_R_pos

function make_pos_def( M; negOrSmallTh = 100*eps(), smallVal = 100*eps(), warnTh = 0.01) 
    eig = eigen(M)
    minEigenVal = minimum(eig.values)

    if minEigenVal < negOrSmallTh
        Logging.@info("negative eigenvalue $minEigenVal")
        absMinEigVal = abs(minEigenVal)

        deltaEigenVal = absMinEigVal + smallVal

        absMinEigVal > warnTh ? Logging.@warn("large negative eigenvalue $minEigenVal") : ()

        eig.values[eig.values.< smallVal] .+= deltaEigenVal

        M = Matrix(eig)
    end
    
    return M, minEigenVal
end
export make_pos_def



drop_nan_col(x::Matrix) = x[:, .!dropdims(any(isnan.(x), dims=1), dims=1) ]
export drop_nan_col


drop_nan_row(x::Matrix) = x[.!dropdims(any(isnan.(x), dims=2), dims=2), : ]
export drop_nan_row


#Funzioni che uso spesso
sumSq(array :: AbstractArray,Dim::Int) = dropdims(sum(array,dims = Dim),dims =Dim)
sumSq(array :: AbstractArray) = sumSq(array,1)
export sumSq


meanSq(array::AbstractArray,Dim::Int) = dropdims(mean(array,dims = Dim),dims =Dim)
export meanSq


splitVec(Vec::AbstractArray{<:Any,1}) = (Nhalf = round(Int,length(Vec)/2); (Nhalf,Vec[1:Nhalf],Vec[Nhalf+1:2Nhalf]) )
export splitVec


splitMat(Mat::AbstractArray{<:Any,2}) = (Nhalf = round(Int,length(Mat[1,:])/2); (Nhalf,Mat[:,1:Nhalf],Mat[:,Nhalf+1:2Nhalf]) )
export splitMat


"""put to zero diagonal of matrices"""
putZeroDiag!(Matr::Array{T,2} where T<: Real) =   for i=1:length(Matr[:,1]) Matr[i,i] = 0 end
export putZeroDiag!


#putZeroDiag!(Matr::Symmetric{T,Array{T,2}} where T<: Real) =   for i=1:length(Matr[:,1]) Matr[i,i] = 0 end
putZeroDiag(Matr::Array{T,2} where T<: Real) = (tmp = Matr; putZeroDiag!(tmp); tmp)
#putZeroDiag(Matr::Symmetric{T,Array{T,2}} where T<: Real) = (tmp = Matr; putZeroDiag!(tmp); tmp)
putZeroDiag!(Matr::BitArray{2}) =   for i=1:length(Matr[:,1]) Matr[i,i] = false end
putZeroDiag(Matr::BitArray{2}) = (tmp = Matr; putZeroDiag!(tmp); tmp)
export putZeroDiag


putZeroDiag_no_mut(Matr::Array{T,2} where T<: Real) = (Matr - Diagonal(Diagonal(Matr)))
export putZeroDiag_no_mut


function sortrowsUtilities(A::Matrix{<:Number},rowInd::Int;rev = true)
    tmp  = sortperm(A[:,rowInd],rev=rev)
    return A[tmp,:]
end
export sortrowsUtilities


#define an adjacency matrix that has falses on diagonals: central upper and lower
tmpAjacMatTridiagFalses(N::Int) = (tmp = trues(N,N);for i=1:N tmp[i,i] = false; i>1 ? tmp[i,i-1]= 0 : (); i<N ? tmp[i,i+1]=0 : ()  end;tmp)
export tmpAjacMatTridiagFalses


function tmpArrayAdjmatrEqual(N::Int,T::Int;matGenFun::Function = tmpAjacMatTridiagFalses)
    out = trues(N,N,T)
    mat = matGenFun(N)
    for t=1:T
        out[:,:,t] = mat
    end
    return out
end
export tmpArrayAdjmatrEqual


function distributeAinB!(A::Array{<:Real,1},B::Array{<:Real,1})
    NA = length(A)
    NB = length(B)
    GroupSize  = floor(NB/NA)
    for i =0:NA-1 B[1 + Int64( GroupSize *i):Int64(GroupSize*(i+1))] = A[i+1] end;
    B[Int64(GroupSize*(NA-1))+1:end] = A[NA];
end
export distributeAinB!


function distributeAinVecN(A::Array{<:Real,1},N::I where I<: Int64)
    if length(A) == 0
        return zeros(Int64,N);
    else
        vecN = zeros(typeof(A[1]),N);
        distributeAinB!(A,vecN);
        return vecN
    end
end
distributeAinVecN(A::UnitRange,N::I where I<: Int64) = distributeAinVecN(Array{Int,1}(A),N)
export distributeAinVecN


function rocCurve(realVals::BitArray{1},foreVals::Vector{Float64}; Th::Vector{Float64}=Vector(0:0.00001:1),plotFlag = true)
    #Given binary observations and estimated probabilityes plot the roc and return true positives and false positives
    rocOut = roc(realVals,foreVals,Th)
    Nth = length(Th)
    tpr = [rocOut[i].tp/rocOut[i].p for i=1:Nth]
    fpr = [rocOut[i].fp/rocOut[i].n for i=1:Nth]
    auc = sum(-Base.diff(fpr).*tpr[1:end-1])
    if plotFlag
        plot(fpr,tpr)
        xlim([0,1])
        ylim([0,1])
        title("AUC =  $(round(auc,3))")
        grid()
    end
    return tpr, fpr,auc
end
export rocCurve


function collapseArr3(arr)
     N1,N2,N3 = length(arr),length(arr[1]),length(arr[1][1])
     tmp = zeros(N1,N2,N3)
     for i=1:N1,j=1:N2,k=1:N3
         tmp[i,j,k] = arr[i][j][k]
     end
     tmp = permutedims(tmp,[2,1,3])
     return tmp
end
export collapseArr3


function randSteps(startVal,endVal,Nsteps::Int,T::Int;rand=true)
     out = zeros(T)
     if Nsteps>1
         heights = linspace(startVal,endVal,Nsteps)
         if rand
             heights = shuffle(heights)
         end
         Tstep = round(Int,T/Nsteps)
         last = 0
         for i=1:Nsteps-1
             out[1+last:last+Tstep] = heights[i]
             last +=Tstep
         end
         out[last:end] = heights[Nsteps]
     elseif Nsteps == 1
         out =startVal
     end
     return out
end
export randSteps


function dgpSin(minVal::Ty  where Ty <:Real  ,maxVal::Ty  where Ty <:Real ,Ncycles,T::Int; phase =0)
     out = zeros(T)
     maxVal<=minVal ? error() : ()
     medVal = (maxVal + minVal)./2
     amplitude = abs.(maxVal-minVal)
     out = (medVal .+ (amplitude/2) .* sin.(2π*Ncycles/T .*  Vector(1:T) .+ phase))
     return out
end
export dgpSin


function dgpAR(mu,B,sigma,T::Int;minMax=zeros(2),scaling = "uniform")
  ar1 = AReg.ARp([mu*(1-B) , B],sigma,[mu])
  path = AReg.simulate(ar1,T)

  if minMax !=zeros(2)
      min = minimum(minMax)
      max = maximum(minMax)
      if scaling == "uniform"
          minPath = minimum(path)
          maxPath = maximum(path)
          Δ = (max-min)/( maxPath - minPath)
          rescPath = min .+ (path .- minPath).*Δ
      elseif scaling == "nonlinear"
          rescPath = min + (max - min)* 1 ./(1 + exp.(path))
      end

  else
      rescPath = path
  end
  return rescPath
end
export dgpAR


end
