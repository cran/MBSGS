MBGLSS = function(Y, X, niter = 10000, burnin = 5000, group_size, a=1, b=1,num_update = 100, niter.update = 100,
                 verbose = FALSE, pi_prior=TRUE, pi=0.5, d=3,update_tau=TRUE,option.update="global")
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
  tau2 = rep(1, ngroup)
  Sigma <- diag(q)
  #sigma2 = 1
  l = rep(0, ngroup)
  beta = vector(mode='list', length=ngroup)
  for(i in 1:ngroup) beta[[i]]=matrix(0,ncol=q,nrow=group_size[i])
  Z = rep(0, ngroup)

  ##################################
  # Compute lambda2 via EM         #
  ##################################
  ### TO DO UPDATE of LAMBDA
  fit_for_lambda2 = MBGLSS_EM_lambda(k,Y, X,num_update = num_update, niter = niter.update, group_size = group_size,option.update=option.update)

  lambda2 = apply(fit_for_lambda2$lambda2_path,2,tail,1)

  ###############################
  # avoid duplicate computation #
  ###############################

  YtY = t(Y) %*% Y
  XtY = t(X) %*% Y
  XtX = t(X)%*% X

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

  coef = vector("list", length=q)
  for(jj in 1:q){
    coef[[jj]] <- array(0, dim=c(p, niter-burnin))
  }


  for (iter in 1:niter)
  {
    # print the current number of iteration
    if (verbose == TRUE) {print(iter)}

    # Update beta's
    for(i in 1:ngroup)
    {
      bmk = NULL
      for(j in 1:ngroup)
      {
        if(j!=i) bmk = rbind(bmk, beta[[j]])
      }
      f1 = XktY[[i]] - XktXmk[[i]] %*% bmk
      f2 = XktXk[[i]]+1/tau2[i]*diag(nrow=group_size[i]) ### Sigma_g inverse
      f2_inverse = solve(f2) ### Sigma_g
      mu = f2_inverse %*% f1
      maxf <- max(f2)

      trythis <- (-q*group_size[i]/2)*log(tau2[i]) + (-q/2)*log(det(f2/maxf)) + (-q*dim(f2)[1]/2)*log(maxf) + 0.5*sum(diag(solve(Sigma)%*%t(mu)%*%f2%*%mu))

      l[i] = pi/(pi+(1-pi)*exp(trythis))

      if(runif(1)<l[i])
      {
        beta[[i]] = matrix(0,ncol=q,nrow=group_size[i])
        Z[i] = 0
      }
      else
      {

        beta[[i]] = matrix(rmnorm(1, mean=c(mu), varcov=kronecker(Sigma,f2_inverse)),ncol=q,nrow=group_size[i],byrow = FALSE)
        Z[i] = 1
      }
    }


    # Update tau2's
    if(update_tau) {
      for(i in 1:ngroup)
        {
          if(Z[i]==0){tau2[i] = rgamma(1, shape=(q*group_size[i]+1)/2, rate=group_size[i]*lambda2[i]/2)}
          else{tau2[i] = 1/rig(1, mean=sqrt((lambda2[i]*group_size[i])/sum(diag(beta[[i]]%*%solve(Sigma)%*%t(beta[[i]])))), scale = 1/(lambda2[i]*group_size[i]))}
        }
     }





    # Update Sigma  (to check value of the hyper-parameter d=3)
    ddf <- d + n +sum(group_size*Z)
    Beta <- do.call(rbind, beta)
    D_tau <- diag(rep(1/tau2,times=group_size))
    res <- Y-X%*%Beta
    S <- t(res)%*%res+t(Beta)%*%D_tau%*%Beta + Q
    Sigma <- riwish(ddf, S)

    ### save results about beta
    if(iter > burnin){
      for(jj in 1:q){
      coef[[jj]][,iter-burnin] <- Beta[,jj]
      }
    }

    # Update pi
    if(pi_prior==TRUE)
      pi = rbeta(1, shape1=a+ngroup-sum(Z), shape2=b+sum(Z)) ### email for confirmation
  }

  # output the posterior mean and median as our estimator

  pos_mean = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, mean))),ncol=q,nrow=p,byrow=FALSE)
  pos_median = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, median))),ncol=q,nrow=p,byrow=FALSE)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef)
}




