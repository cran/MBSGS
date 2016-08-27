\name{MBGLSS}
\alias{MBGLSS}
\title{
Multivariate Bayesian Group Lasso with Spike and Slab prior
}
\description{
Run a gibbs sampler for a Multivariate Bayesian group lasso model with spike and slab prior. This function is designed for a regression model with multivariate response, where the design matrix has a group structure.
}
\usage{
MBGLSS(Y, X, niter = 10000, burnin = 5000, group_size,
a = 1, b = 1, num_update = 100, niter.update = 100,
verbose = FALSE, pi_prior = TRUE, pi = 0.5,
d = 3, update_tau = TRUE, option.update = "global")
}
\arguments{
\item{Y}{A numerical vector representing the univariate response variable.}
  \item{X}{A matrix respresenting the design matrix of the linear regression model.}
  \item{niter}{Number of iteration for the Gibbs sampler.}
  \item{burnin}{Number of burnin iteration}
  \item{group_size}{Integer vector representing the size of the groups of the  design matrix \code{X}}
  \item{a}{First shape parameter of the conjugate beta prior for \code{pi_0}. Default is 1.}
  \item{b}{Second shape parameter of the conjugate beta prior for \code{pi_0}. Default is 1.}
  \item{num_update}{Number of update regarding the scaling of the shrinkage parameter lambda which is calibrated by a Monte Carlo EM algorithm}
  \item{niter.update}{Number of itertion regarding the scaling of the shrinkage parameter lambda which is calibrated by a Monte Carlo EM algorithm}
  \item{verbose}{Logical. If "TRUE" iterations are displayed.}
  \item{pi_prior}{Logical. If "TRUE" a beta prior is used for pi}
  \item{pi}{Initial value for pi_0 which will be updated if \code{pi_prior}="TRUE""}
  \item{d}{Degree of freedom of the inverse Wishart prior of the covariance matrix of the response variable. By default \code{d} is set to 3.}
\item{update_tau}{Logical. If "TRUE" then a Monte Carlo EM algorithm is used to update lambda}
\item{option.update}{Two options are proposed for updating lambda. A "Local" update or a "Global" update}
}

\value{
\code{BSGSSS} returns a list that contains the following components:

\item{pos_mean}{The posterior mean estimate of the regression coefficients}
\item{pos_median}{The posterior mean estimate of the regression coefficients}
\item{coef}{A matrix with the regression coefficients sampled at each iteration}
}


\references{
B. Liquet, K. Mengersen, A. Pettitt and M. Sutton. (2016). Bayesian Variable Selection Regression Of Multivariate Responses For Group Data. \emph{Submitted in Bayesian Analysis}.
}
\author{
Benoit Liquet and Matthew Sutton.
}



\seealso{
\code{\link{MBSGSSS}}
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

model <- MBGLSS(Y,X,niter,burnin,gsize,num_update = 100,
niter.update = 100)
model$pos_median[,1]!=0
}}
\keyword{regression}
\keyword{multivariate}
