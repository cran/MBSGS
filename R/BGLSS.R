BGLSS = function(Y, X, niter = 10000, burnin = 5000, group_size, a=1, b=1,num_update = 100, niter.update =100,
                 verbose = FALSE, alpha=1e-1,gamma=1e-1, pi_prior=TRUE, pi=0.5, update_tau=TRUE,option.weight.group=FALSE,option.update="global",lambda2_update=NULL)
{
  ####################################
  # Create and Initialize parameters #
  ####################################
  n = length(Y)
  p = dim(X)[2]
  ngroup = length(group_size)
  # Initialize parameters
  tau2 = rep(1, ngroup)
  sigma2 = 1
  l = rep(0, ngroup)
  beta = vector(mode='list', length=ngroup)
  for(i in 1:ngroup) beta[[i]]=rep(0, group_size[i])
  Z = rep(0, ngroup)

  ##################################
  # Compute lambda2 via EM         #
  ##################################
  if (update_tau==TRUE){
  fit_for_lambda2 = BGLSS_EM_lambda(Y, X, num_update = num_update, niter = niter.update,group_size = group_size,option.update=option.update,option.weight.group=option.weight.group)
  lambda2 = apply(fit_for_lambda2$lambda2_path,2,tail,1)}else
  {
    lambda2 <- lambda2_update
  }

  #print(c("lambda2",option.update))

  ###############################
  # avoid duplicate computation #
  ###############################

  YtY = t(Y) %*% Y
  XtY = t(X) %*% Y
  XtX = t(X) %*% X

  XktY = vector(mode='list', length=ngroup)
  XktXk = vector(mode='list', length=ngroup)
  XktXmk = vector(mode='list', length=ngroup)

  begin_idx = 1
  for(i in 1:ngroup)
  {
    end_idx = begin_idx + group_size[i] - 1
    Xk = X[,begin_idx:end_idx]
    XktY[[i]] = t(Xk) %*% Y
    XktXk[[i]] = t(Xk) %*% Xk
    XktXmk[[i]] = t(Xk) %*% X[,-(begin_idx:end_idx)]
    begin_idx = end_idx + 1
  }

  #####################
  # The Gibbs Sampler #
  #####################

  # burnin
  coef = array(0, dim=c(p, niter-burnin))
  coef_tau = array(0, dim=c(ngroup, niter))
  for (iter in 1:niter)
  {
    # print the current number of iteration
    if (verbose == TRUE) {print(iter)}

    # Update beta's
    for(i in 1:ngroup)
    {
      bmk = c()
      for(j in 1:ngroup)
      {
        if(j!=i) bmk = c(bmk, beta[[j]])
      }
      f1 = XktY[[i]] - XktXmk[[i]] %*% bmk
      f2 = XktXk[[i]]+1/tau2[i]*diag(nrow=group_size[i])
      f2_inverse = solve(f2)
      mu = f2_inverse %*% f1
       ### Main part
      l[i] = pi/(pi+(1-pi)*(tau2[i])^(-group_size[i]/2)*det(f2)^(-1/2)*exp(t(f1)%*%mu/(2*sigma2)))
      maxf <- max(f2)
      trythis <- (-group_size[i]/2)*log(tau2[i]) + (-1/2)*log(det(f2/maxf)) + (-dim(f2)[1]/2)*log(maxf) + t(f1)%*%mu/(2*sigma2)
      l[i] = pi/(pi+(1-pi)*exp(trythis))

      if(runif(1)<l[i])
      {
        beta[[i]] = rep(0, group_size[i])
        Z[i] = 0
      }
      else
      {
        beta[[i]] = rmnorm(1, mean=mu, varcov=sigma2*f2_inverse)
        Z[i] = 1
      }
    }

  # Update tau2's



    if(update_tau) {

      if (option.weight.group== FALSE){
        for(i in 1:ngroup)
        {
          if(Z[i]==0){tau2[i] = rgamma(1, shape=(group_size[i]+1)/2, rate=lambda2[i]/2)}
          else{tau2[i] = 1/rig(1, mean=sqrt(lambda2[i]*sigma2/sum(beta[[i]]^2)), scale = 1/(lambda2[i]))}
        }
      }else{
      for(i in 1:ngroup)
      {
        if(Z[i]==0){tau2[i] = rgamma(1, shape=(group_size[i]+1)/2, rate=lambda2[i]*(group_size[i])/2)}
        else{tau2[i] = 1/rig(1, mean=sqrt(lambda2[i]*group_size[i]*sigma2/sum(beta[[i]]^2)), scale = 1/(group_size[i]*lambda2[i]))}
      coef_tau[i,iter]  = tau2[i]

      }

      }
    }


    # Update sigma2
    s=0
    for(i in 1:ngroup)
    {
      s = s + sum(beta[[i]]^2)/tau2[i]
    }
    beta_vec = c()
    for(j in 1:ngroup) beta_vec = c(beta_vec, beta[[j]])
    if(iter > burnin)
      coef[,iter-burnin] = beta_vec
    sigma2 = rinvgamma(1, shape=(n-1)/2 + sum(Z*group_size)/2 + alpha,
                       scale=(YtY-2*t(beta_vec)%*%XtY+t(beta_vec)%*%XtX%*%beta_vec+s)/2 + gamma)

    # Update pi
    if(pi_prior==TRUE)
    pi = rbeta(1, shape1=a+ngroup-sum(Z), shape2=b+sum(Z))
  }

  # output the posterior mean and median as our estimator
  pos_mean = apply(coef, 1, mean)
  pos_median = apply(coef, 1, median)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef)
}




