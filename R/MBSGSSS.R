# Gibbs Sampler for Multivariate Bayesian Sparse Group Lasso with ZIMP
# we use a truncated normal prior for standard deviation \tau
MBSGSSS = function(Y, X, group_size, pi0 = 0.5, pi1 = 0.5, a1 = 1, a2 = 1,
                    c1 = 1, c2 = 1, pi_prior = TRUE, niter = 10000,burnin = 5000,d=3,num_update = 100, niter.update =100)
{

  ####################################
  # Create and Initialize parameters #
  ####################################
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(X)[2]
  ###Multivariate linear model
  
  if(p<n){
	model <- lm(Y~X)
	k <- mean(apply(model$residuals,2,var))}else{
	var.res <- NULL
	for(l in 1:q){
  	mydata <- data.frame(y=Y[,l],X)
  	min.model <- lm(y~-1,data=mydata)
  	formula.model <- formula(lm(y~0+.,data=mydata))
  	model <- stepAIC(min.model,direction='forward',scope=formula.model,steps=n-1,trace=FALSE)
  	var.res <- c(var.res,var(residuals(model)))	
}
	k <- mean(var.res)
	}
  
  
  Q = k*diag(q)
  ngroup = length(group_size)
  # Initialize parameters
  #tau2 = rep(1, ngroup)
  Sigma <- diag(q)
  #sigma2 = 1
  l = rep(0, ngroup)
  b = vector(mode='list', length=ngroup)
  beta = vector(mode='list', length=ngroup)
  for(i in 1:ngroup){ beta[[i]]=matrix(1,ncol=q,nrow=group_size[i])
  b[[i]]=matrix(1,ncol=q,nrow=group_size[i])}

  Z = rep(0, ngroup)
  m = dim(X)[2]


  # initializing parameters
  b_prob = rep(0.5, ngroup)
  tau = rep(1, m)
  tau_prob = rep(0.5, m)


  #sigma2 = 1
  s2 = 1

  ###############################
  # EM for t #
  ############
  fit_for_t =MBSGSSS_EM_t(k,Y, X, group_size=group_size,num_update = num_update, niter=niter.update)
  t = tail(fit_for_t$t_path, 1)
  #t <- 0.7940479
  #print(t)
  #stop("after EM")
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



  YtXi <- vector(mode='list', length=m)
  XmitXi = array(0, dim=c(m,m-1))
  XitXi = rep(0, m)
  for(j in 1:m)
  {
    YtXi[[j]] = crossprod(Y, X[,j])
    XmitXi[j,] = crossprod(X[,-j], X[,j])
    XitXi[j] = crossprod(X[,j], X[,j])
  }

  #####################
  # The Gibbs Sampler #
  #####################

  # number of iterations
  #niter = 10000
  # burnin
  #burnin = 5000
  # initialize coefficients vector

  coef = vector("list", length=q)
  for(jj in 1:q){
    coef[[jj]] <- array(0, dim=c(p, niter-burnin))
  }



  #coef = array(0, dim=c(m, niter-burnin))
  bprob = array(-1, dim=c(ngroup, niter-burnin))
  tauprob = array(-1, dim=c(m, niter-burnin))

  for (iter in 1:niter)
  {
    #print(iter)
    Z = rep(0, ngroup)
    # number of non zero groups
    n_nzg = 0
    # number of non zero taus
    n_nzt = 0

    ### Generate tau ###
    unif_samples = runif(m)
    idx = 1
    Sigma_inv <- solve(Sigma)
    Beta <- do.call(rbind, beta)
    matb <- do.call(rbind, b)
    for(g in 1:ngroup)
    {
      ### verifier b et beta tilde avec le group jjjj
      ####
      for (j in 1:group_size[g])
      {
        M = Sigma_inv%*%(YtXi[[idx]]-t(Beta[-idx,])%*%matrix(XmitXi[idx,],ncol=1,nrow=m-1))%*%matrix(matb[idx,],nrow=1,ncol=q)
        bsquare <- t(matrix(matb[idx,],ncol=q,nrow=1))%*%matrix(matb[idx,],ncol=q,nrow=1)
        vgj_square = 1/(sum(diag(XitXi[idx]*Sigma_inv%*%bsquare))+1/s2)
        ugj = sum(diag(M))*vgj_square

        ##M = (YtXi[idx] * b[idx] - crossprod(beta[-idx], XmitXi[idx,]) * b[idx])/sigma2
        ##N = 1/s2 + XitXi[idx] * b[idx]^2 / sigma2
        # we need to aviod 0 * Inf
        tau_prob[idx] = pi1/( pi1 + 2 * (1-pi1) * ((s2)^(-0.5))*(vgj_square)^(0.5)*
                                exp(pnorm(ugj/sqrt(vgj_square),log.p=TRUE) + ugj^2/(2*vgj_square)))
        if(unif_samples[idx] < tau_prob[idx]) {
          tau[idx] = 0
        } else {
          tau[idx] = rtruncnorm(1, mean=ugj, sd=sqrt(vgj_square), a=0)
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
      #Sig = solve(crossprod(Vg, XktXk[[g]]) %*% Vg / sigma2 + diag(group_size[g]))
      Sig = solve(crossprod(Vg, XktXk[[g]]) %*% Vg  + diag(group_size[g]))
      dev = (Vg %*% (XktY[[g]] - XktXmk[[g]] %*% Beta[-c(begin_idx:end_idx),]))

      # probability for bg to be 0 vector
      num <- 0.5*sum(diag(solve(Sigma)%*%crossprod(dev,Sig)%*%dev))
      exp.num <- exp(num)
      if(exp.num == Inf) exp.num <- 1.797693e+307
      b_prob[g] = pi0 / ( pi0 + (1-pi0) * det(Sig)^(q/2) * exp.num)
      if(unif_samples_group[g] < b_prob[g]){
        matb[begin_idx:end_idx,] = 0
       b[[g]] = matrix(0,ncol=q,nrow=group_size[g])
       Z[g] = 0}
      else
      {
#        b[begin_idx:end_idx] = mvrnorm(1, mu = Sig%*%dev, Sigma = Sig)
        #print(dim(dev))
        #print(dim(Sig))
        #print(dim((Sig%*%dev)))
        #print(dim(kronecker(Sigma,Sig)))
        #print(c(begin_idx,end_idx))
        matb[begin_idx:end_idx,] = matrix(rmnorm(1, mean=c(Sig%*%dev), varcov=kronecker(Sigma,Sig)),ncol=q,nrow=group_size[g],byrow = FALSE)
        b[[g]] = matb[begin_idx:end_idx,]
        Z[g] = 1
        n_nzg = n_nzg + 1
      }
      # Compute beta
      Beta[begin_idx:end_idx,] = Vg %*% matb[begin_idx:end_idx,]
      beta[[g]] <- Beta[begin_idx:end_idx,]
      # Update index
      begin_idx = end_idx + 1
    }

      ### save results about beta
      if(iter > burnin){
        for(jj in 1:q){
          coef[[jj]][,iter-burnin] <- Beta[,jj]
        }
        bprob[,iter-burnin] = b_prob
        tauprob[,iter-burnin] = tau_prob
      }


      # Update Sigma  (to check value of the hyper-parameter d=3)
      ddf <- d + n +sum(group_size*Z)
      #Beta <- do.call(rbind, beta)
      #D_tau <- diag(rep(1/tau2,times=group_size))
      res <- Y-X%*%Beta
      #S <- t(res)%*%res+t(matb)%*%D_tau%*%matb + Q
      S <- t(res)%*%res+t(matb)%*%matb + Q
      Sigma <- riwish(ddf, S)


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
  # mean = t(apply(coef, 1, function(x) cumsum(x)/c(1:(niter-burnin))))
    pos_mean = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, mean))),ncol=q,nrow=p,byrow=FALSE)
    pos_median = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, median))),ncol=q,nrow=p,byrow=FALSE)

  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef)
}

