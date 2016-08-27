\name{MBSGSSS}
\alias{MBSGSSS}
\title{
Multivariate Bayesian Sparse Group Selection with Spike and Slab priors
}
\description{
Run a gibbs sampler for a Multivariate Bayesian sparse group selection model with spike and slab prior. This function is designed for a regression model with multivariate response, where the design matrix has a group structure.
}
\usage{
MBSGSSS(Y, X, group_size, pi0 = 0.5, pi1 = 0.5,
a1 = 1, a2 = 1, c1 = 1,c2 = 1, pi_prior = TRUE,
niter = 10000, burnin = 5000, d = 3,
num_update = 100, niter.update = 100)
}

\arguments{
\item{Y}{A numerical vector representing the univariate response variable.}
  \item{X}{A matrix respresenting the design matrix of the linear regression model.}
 \item{group_size}{Integer vector representing the size of the groups of the  design matrix \code{X}}
  \item{pi0}{Initial value for pi0 which will be updated if \code{pi_prior}="TRUE"}
  \item{pi1}{Initial value for pi1 which will be updated if \code{pi_prior}="TRUE"}
 \item{a1}{First shape parameter of the conjugate beta hyper-prior for \code{pi_0}. Default is 1.}
  \item{a2}{Second shape parameter of the conjugate beta prior for \code{pi_0}. Default is 1.}
  \item{c1}{First shape parameter of the conjugate beta hyper-prior for \code{pi_1}. Default is 1.}
  \item{c2}{Second shape parameter of the conjugate beta prior for \code{pi_1}. Default is 1.}
 \item{pi_prior}{Logical. If "TRUE" beta priors are used for pi0 and pi1}
  \item{niter}{Number of iteration for the Gibbs sampler.}
\item{burnin}{Number of burnin iteration}
 \item{d}{Degree of freedom of the inverse Wishart prior of the covariance matrix of the response variable. By default \code{d} is set to 3.}
  \item{num_update}{Number of update regarding the scaling of the shrinkage parameter lambda which is calibrated by a Monte Carlo EM algorithm}
  \item{niter.update}{Number of itertion regarding the scaling of the shrinkage parameter lambda which is calibrated by a Monte Carlo EM algorithm}
}

\references{
B. Liquet, K. Mengersen, A. Pettitt and M. Sutton. (2016). Bayesian Variable Selection Regression Of Multivariate Responses For Group Data. \emph{Submitted in Bayesian Analysis}.
}
\author{
Benoit Liquet and Matthew Sutton.
}



\seealso{
\code{\link{MBGLSS}}
}
\examples{
\dontrun{
## Simulation of datasets X and Y with group variables
data1 = gen_data_Multi(nsample = 120, ntrain = 80)
data1 = Mnormalize(data1)

true_model <- data1$true_model
X <- data1$X
Y<- data1$Y
train_idx <- data1$train_idx
gsize <- data1$gsize
niter <- 2000
burnin <- 1000

model <- MBSGSSS(Y,X,niter=niter,burnin=burnin,group_size=gsize,
num_update = 50,niter.update = 50)
model$pos_median[,1]!=0
}}
\keyword{regression}
\keyword{multivariate}
