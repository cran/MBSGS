\name{BSGSSS}
\alias{BSGSSS}
\title{
Bayesian Sparse Group Selection with Spike and Slab priors
}
\description{
Run a gibbs sampler for a Bayesian sparse group selection model with spike and slab priors. This function is designed for an univariate response model and when the design matrix has a group structure.
}
\usage{
BSGSSS(Y, X, group_size, niter = 10000, burnin = 5000,
pi0 = 0.5, pi1 = 0.5, num_update = 100, niter.update = 100,
alpha = 0.1, gamma = 0.1, a1 = 1, a2 = 1, c1 = 1, c2 = 1,
pi_prior = TRUE)
}

\arguments{
\item{Y}{A numerical vector representing the univariate response variable.}
  \item{X}{A matrix respresenting the design matrix of the linear regression model.}
 \item{group_size}{Integer vector representing the size of the groups of the  design matrix \code{X}}
 \item{niter}{Number of iteration for the Gibbs sampler.}
  \item{burnin}{Number of burnin iteration}
  \item{pi0}{Initial value for pi0 which will be updated if \code{pi_prior}="TRUE"}
  \item{pi1}{Initial value for pi1 which will be updated if \code{pi_prior}="TRUE"}
 \item{num_update}{Number of update regarding the scaling of the shrinkage parameter lambda which is calibrated by a Monte Carlo EM algorithm}
  \item{niter.update}{Number of itertion regarding the scaling of the shrinkage parameter lambda which is calibrated by a Monte Carlo EM algorithm}
 \item{alpha}{Shape parameter of the Inverse Gamma prior on the variance of the noise for the linear regression model.}
  \item{gamma}{Scale parameter of the Inverse Gamma prior on the variance of the noise for the linear regression model.}
  \item{a1}{First shape parameter of the conjugate beta hyper-prior for \code{pi_0}. Default is 1.}
  \item{a2}{Second shape parameter of the conjugate beta prior for \code{pi_0}. Default is 1.}
  \item{c1}{First shape parameter of the conjugate beta hyper-prior for \code{pi_1}. Default is 1.}
  \item{c2}{Second shape parameter of the conjugate beta prior for \code{pi_1}. Default is 1.}
  \item{pi_prior}{Logical. If "TRUE" beta priors are used for pi0 and pi1}
  }

\value{
\code{BSGSSS} returns a list that contains the following components:

\item{pos_mean}{The posterior mean estimate of the regression coefficients}
\item{pos_median}{The posterior mean estimate of the regression coefficients}
\item{coef}{A matrix with the regression coefficients sampled at each iteration}
}

\references{
B. Liquet, K. Mengersen, A. Pettitt and M. Sutton. (2016). Bayesian Variable Selection Regression Of Multivariate Responses For Group Data. \emph{Submitted in Bayesian Analysis}.

Xu, X. and Ghosh, M. (2015). Bayesian Variable Selection and Estimation for Group Lasso.  \emph{Bayesian Analysis}, 10(4): 909-936.

}
\author{
Benoit Liquet, Matthew Sutton and Xiaofan Xu.
}



\seealso{
\code{\link{BGLSS}}
}
\examples{

## Simulation of datasets X and Y with group variables
set.seed(1)
data1 = gen_data_uni(nsample = 120,cor.var=0.5, ntrain = 80)
data1 = normalize(data1)

true_model <- data1$true_model
X <- data1$X
Y<- data1$Y
train_idx <- data1$train_idx
gsize <- data1$gsize
## We recommend to set niter=50000, burnin=10000
## num_update = 100 and niter.update = 100 
## to reach convergence

model <- BSGSSS(Y[,1],X,niter=500,burnin=100,group_size=gsize,
num_update = 20,niter.update = 20)
model$pos_median!=0
true_model
}
\keyword{regression}
