BSGSSS = function(Y, X, group_size,niter = 10000, burnin = 5000,pi0 = 0.5, pi1 = 0.5,num_update = 100, niter.update =100,
                    alpha = 1e-1, gamma = 1e-1, a1 = 1, a2 = 1,
                    c1 = 1, c2 = 1, pi_prior = TRUE)
{
  ####################################
  # Create and Initialize parameters #
  ####################################

  # book keeping
  n = length(Y)
  m = dim(X)[2]
  ngroup = length(group_size)

  # initializing parameters
  b = rep(1, m)
  b_prob = rep(0.5, ngroup)
  tau = rep(1, m)
  tau_prob = rep(0.5, m)
  beta = rep(1, m)
  sigma2 = 1
  s2 = 1

  ###############################
  # EM for t #
  ############
  fit_for_t = BSGSSS_EM_t(Y, X, group_size,num_update = num_update, niter = niter.update)
  t = tail(fit_for_t$t_path, 1)

  ###############################
  # avoid replicate computation #
  ###############################

  YtY = crossprod(Y, Y)
  XtY = crossprod(X, Y)
  XtX = crossprod(X, X)

  # group
  XktY = vector(mode='list', length=ngroup)
  XktXk = vector(mode='list', length=ngroup)
  XktXmk = vector(mode='list', length=ngroup)

  begin_idx = 1
  for(i in 1:ngroup)
  {
    end_idx = begin_idx + group_size[i] - 1
    Xk = X[,begin_idx:end_idx]
    XktY[[i]] = crossprod(Xk, Y)
    XktXk[[i]] = crossprod(Xk, Xk)
    XktXmk[[i]] = crossprod(Xk, X[,-(begin_idx:end_idx)])
    begin_idx = end_idx + 1
  }

  # individual
  YtXi = rep(0, m)
  XmitXi = array(0, dim=c(m,m-1))
  XitXi = rep(0, m)
  for(j in 1:m)
  {
    YtXi[j] = crossprod(Y, X[,j])
    XmitXi[j,] = crossprod(X[,-j], X[,j])
    XitXi[j] = crossprod(X[,j], X[,j])
  }

  #####################
  # The Gibbs Sampler #
  #####################


  # initialize coefficients vector
  coef = array(0, dim=c(m, niter-burnin))
  bprob = array(-1, dim=c(ngroup, niter-burnin))
  tauprob = array(-1, dim=c(m, niter-burnin))

  for (iter in 1:niter)
  {
    # number of non zero groups
    n_nzg = 0
    # number of non zero taus
    n_nzt = 0

    ### Generate tau ###
    unif_samples = runif(m)
    idx = 1
    for(g in 1:ngroup)
    {
      for (j in 1:group_size[g])
      {
        M = (YtXi[idx] * b[idx] - crossprod(beta[-idx], XmitXi[idx,]) * b[idx])/sigma2
        N = 1/s2 + XitXi[idx] * b[idx]^2 / sigma2
        # we need to aviod 0 * Inf
        tau_prob[idx] = pi1/( pi1 + 2 * (1-pi1) * (s2*N)^(-0.5)*
                                exp(pnorm(M/sqrt(N),log.p=TRUE) + M^2/(2*N)) )
        if(unif_samples[idx] < tau_prob[idx]) {
          tau[idx] = 0
        } else {
          tau[idx] = rtruncnorm(1, mean=M/N, sd=sqrt(1/N), a=0)
          n_nzt = n_nzt + 1
        }
        idx = idx + 1
      }
    }

    # generate b and compute beta
    begin_idx = 1
    unif_samples_group = runif(ngroup)
    for(g in 1:ngroup)
    {
      end_idx = begin_idx + group_size[g] - 1
      # Covariance Matrix for Group g
      Vg = diag(tau[begin_idx:end_idx])
      #
      Sig = solve(crossprod(Vg, XktXk[[g]]) %*% Vg / sigma2 + diag(group_size[g]))
      dev = (Vg %*% (XktY[[g]] - XktXmk[[g]] %*% beta[-c(begin_idx:end_idx)]))/sigma2
      # probability for bg to be 0 vector
      b_prob[g] = pi0 / ( pi0 + (1-pi0) * det(Sig)^(0.5) * exp(crossprod(dev,Sig)%*%dev/2) )
      if(unif_samples_group[g] < b_prob[g])
        b[begin_idx:end_idx] = 0
      else
      {
        b[begin_idx:end_idx] = t(rmnorm(1, mean=Sig%*%dev, varcov=Sig))
        n_nzg = n_nzg + 1
      }
      # Compute beta
      beta[begin_idx:end_idx] = Vg %*% b[begin_idx:end_idx]
      # Update index
      begin_idx = end_idx + 1
    }


    # store coefficients vector beta
    if(iter - burnin > 0)
    {
      coef[,iter-burnin] = beta
      bprob[,iter-burnin] = b_prob
      tauprob[,iter-burnin] = tau_prob
    }

    # Generate sigma2
    shape = n/2 + alpha
    scale = (YtY + crossprod(beta,XtX)%*%beta - 2*crossprod(beta, XtY))/2 + gamma
    sigma2 = rinvgamma(1, shape=shape, scale=scale)

    if(pi_prior == TRUE) {
      # Generate pi0
      pi0 = rbeta(1, ngroup - n_nzg + a1, n_nzg + a2)

      # Generate pi1
      pi1 = rbeta(1, m - n_nzt + c1, n_nzt + c2)
    }

    # Update s2
    s2 = rinvgamma(1, shape = 1 + n_nzt/2, scale = t + sum(tau)/2)

  }

  # output is the posterior mean at first (iter-burnin) iteration, iter = burnin+1, ..., niter

  pos_mean = apply(coef, 1, mean)
  pos_median = apply(coef, 1, median)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef)
}


