#' MCMC for Bayesian Spatial Model using Probit link
#'
#' @param path      Path to csv file.
#' @param nsamples      Number of samples to run.
#' @param burn      Number of samples to burn before saving parameter samples.
#' @param thin      Rate to save parameter samples.
compareidealK1D <- function(votes, missing, group, demleader, repleader, nsamples = 30000L, burn = 1000L, thin = 1L, printevery = 250L, samplezeta = TRUE, startconfig = 0, varalpprop = 0.06) {
    
    missing    = as.matrix(missing)
    votes      = as.matrix(votes)
    group      = as.matrix(group)
    samplezeta = as.numeric(samplezeta)
    ngroups    = length(unique(group))
    
    
    print(is.integer(nsamples))
    
    if(ngroups!=(max(group)+1))   stop("groups must be labeled continuously starting at 0")
    if(any(dim(missing)!=dim(y))) stop("The dimensions of votes and missing must match")
    if(dim(y)[2]!=length(group))  stop("The length of groups does not match the number of votes")
    if(nsamples <= burn)          stop("nsamples must be greater than burn")
    if(nsamples <= printevery)    stop("nsamples must be greater than printevery")
    if(!samplezeta%in%c(0,1))      stop("samplezeta must be a logical variable")
    if(!startconfig%in%c(0,1))     stop("startconfig can only be 0 (single cluster for all) or 1 (K clusters for all)")    

     .Call(`_compareKgroups`, ngroups, votes, missing, demleader, repleader, group, nsamples, burn, thin, printevery, samplezeta, startconfig, varalpprop)
   
}

