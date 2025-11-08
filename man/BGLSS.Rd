\name{BGLSS}
\alias{BGLSS}
\title{
Bayesian Group Lasso with Spike and Slab prior
}
\description{
Run a gibbs sampler for a Bayesian group lasso model with spike and slab prior. This function is designed for an univariate response model and when the design matrix has a group structure.
}
\usage{
BGLSS(Y, X, niter = 10000, burnin = 5000, group_size, a = 1,
b = 1, num_update = 100, niter.update = 100, verbose = FALSE,
alpha = 0.1, gamma = 0.1, pi_prior = TRUE, pi = 0.5,
update_tau = TRUE, option.weight.group = FALSE,
option.update = "global", lambda2_update = NULL)
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
  \item{alpha}{Shape parameter of the Inverse Gamma prior on the variance of the noise for the linear regression model.}
  \item{gamma}{Scale parameter of the Inverse Gamma prior on the variance of the noise for the linear regression model.}
  \item{pi_prior}{Logical. If "TRUE" a beta prior is used for pi}
  \item{pi}{Initial value for pi_0 which will be updated if \code{pi_prior}="TRUE""}
  \item{update_tau}{Logical. If "TRUE" then a Monte Carlo EM algorithm is used to update lambda}
  \item{option.weight.group}{If "TRUE" then the group size is used for shrinkage penalty purpose.}
  \item{option.update}{Two options are proposed for updating lambda. A "Local" update or a "Global" update}
  \item{lambda2_update}{Value of the square of lambda when \code{update_tau}="FALSE"}
}

\value{
\code{BGLSS} returns a list that contains the following components:

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
\code{\link{BSGSSS}}
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

model <- BGLSS(Y[,1],X,niter=500,burnin=100,group_size=gsize,
num_update = 20,niter.update = 20)
model$pos_median!=0
true_model

}
\keyword{regression}
