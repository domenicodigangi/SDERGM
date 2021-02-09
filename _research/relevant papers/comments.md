# DDG's comments on relevant papers

##  "A framework for the comparison of maximum pseudo-likelihood and maximum likelihood estimation of exponential family random graph models"


- For quantitative comparisons they only look at RMSE. No attempt is made to quantitatively disentangle bias from variance of the estimators. MLEs look biased in natural parametrization.

- They use wrong estimators for the estimator's variances (personally I would stop here..) : For each estimator, a standard error estimate is derived from the estimated curvature of the corresponding log-likelihood or log-pseudo-likelihood ((5), (10) and (11)). We refer to these as “perceived” standard errors because these are the values formally derived as the standard approximations to the true standard errors from asymptotic arguments that have not been justified for these models.

- I do not see a strong case against MPLE in the natural parametrization :
        - From their boxplots, MPLE has higher variance  but lower bias than MLE 

- Inference in the mean value parametrization is completely unfair on the MPLE. Seems that they proceed as follows:
   1) Get MPLE in natural parametrization  
   2) Estimate the variance of this estimator with the curvature of the pseudolikelihood. This is known to be incorrect!
   3) Transform natural par estimate in a mean value one by sampling and taking the mean
   4) Obtain the variance of the mean val MPLE as the inverse of the (incorrect) variance of the natural par MPLE

    - Why should anyone even consider MPLE in the mean val parametrization ?
    - We dont's consider SDERGM in the mean val parametrization ? 
    - Need to check if we rescaled in the numerical results for P-SDERGM. In any case there is a parameter in front that helps "reducing the error". 




- In the "computational details" part they clearly state that MPLE is always used as starting point of the MLE. That is understandable but does not make the case for a fair comparison. Most importantly:
    - The ergm default is to use the MPLE fit to begin
the MLE fit. This approachwas used for the samples from the original
model. The increased transitivity samples, however,were more
difficult to fit, so the MPLE estimate proved to be too far from the
MLE to lead consistently to converged estimates over large numbers
of sampled networks. Two modifications were introduced to
address this problem. First, the initial parameter estimates were
computed by adding 90% of the MPLE estimate to 10% of the model
parameters from the original model. This had the effect of providing
some correction in cases where the MPLE estimates were very
far from the MLE. Second, the model fit was accomplished using
two ergm calls. The first was less precise and with smaller sample
size, and aimed at producing a rough initial set of parameter estimates
closer to the MLE than its original values. The second ergm
call was started at the estimate produced by the first, and involved
a larger sample size in the interest of producing a more precise final
set of parameter estimates. Note that these sophistications robustify
the estimation of the MLE for the purposes of automation over
the 1000 networks, but do not effect the ultimateMLE itself.   
    

## "Overview of Composite Likelihoods Methods" 
- Notation: m number of elements of vector random variable (as a consequence in many cases it is equalto the  number of composite likelhood's terms), n number of observations 
- For maximum composite likelihoods estimators (MPLE is a particular case) the "information matrix identity" does not hold. 
- The variance of MPLE estimators is asymptotically G^-1, where G = E[h] * (var[s])^-1 * E[h] is  the Godambe or sandwitch information matrix 
- Not clear how to estimate var[s]