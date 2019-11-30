\name{sampling}
\alias{sampling}
\title{MCMC for Bayesian Spatial Model}
\usage{
sampling(path, nsamples, burn, thin)
}
\arguments{
\item{path}{Path to csv file.}

\item{nsamples}{Number of samples to run.}

\item{burn}{Number of samples to burn before saving parameter samples.}

\item{thin}{Rate to save parameter samples.}

\item{link}{Choose the link function. "probit" or "logit".}
}
\value{
Returns a list of matrices consisting of the sampled perameters.
}
