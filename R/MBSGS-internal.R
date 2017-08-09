gen_data_Multi = function( nsample = 120, ntrain = 80) {
  # Generate design matrix
  mu = rep(0, 20)
  Sigma = diag(rep(1,20)) + array(0.5, dim=c(20,20))
  X = mvrnorm(nsample, mu=mu, Sigma=Sigma)

  b1 = c(0.3, -1, 0, 0.5, 0.01, rep(0,5), rep(0.8, 5), rep(0,5))
  b2 = c(0.2, -1.1, 0, 0.6, 0.02, rep(0,5), rep(0.7, 5), rep(0,5))
  b3 = c(0.1, -1.2, 0, 0.7, 0.03, rep(0,5), rep(0.6, 5), rep(0,5))
  beta <- matrix(c(b1, b2, b3), ncol = 3, nrow = 20, byrow = FALSE)

  Sigma <- matrix(c(1, 0.95, 0.3, 0.95, 1, 0.5, 0.3, 0.5, 1), ncol = 3, nrow = 3, byrow = FALSE) # equation (18)
  c <- 4
  # Generate response
  Y = X %*% beta + mvrnorm(nsample, mu = rep(0, 3),Sigma = c*Sigma)
  # Group size
  gsize = c(5, 5, 5, 5)

  # Seperate training set and test set
  train_idx = sample(nsample, size = ntrain)

  # Return value
  list(X=X, Y=Y, gsize=gsize, true_model=(beta[,1]!=0),
       train_idx = train_idx)
}




gen_data_uni=  function(nsample = 120, ntrain = 80,cor.var=0,
                     inter = 0, sd = 1) {
  # Generate design matrix
  mu = rep(0, 20)
  Sigma = diag(rep(1,20)) + array(cor.var, dim=c(20,20))
  X = mvrnorm(nsample, mu=mu, Sigma=Sigma)

  beta = c(0.3, -1, 0, 0.5, 0.01, rep(0,5), rep(0.8, 5), rep(0,5))

  # Generate response
  Y = cbind(1, X) %*% c(inter, beta) + rnorm(n = nsample, sd = sd)

  # Group size
  gsize = c(5, 5, 5, 5)

  # Seperate training set and test set
  train_idx = sample(nsample, size = ntrain)

  # Return value
  list(X=X, Y=Y, gsize=gsize, true_model=(beta!=0),
       train_idx = train_idx)
}






# Centralize and scale data
normalize = function(data) {
  # data must be a data frame with 'X' and 'Y' elements

  # Centralize Y
  data$Y = data$Y - mean(data$Y[data$train_idx])

  # Normalize X
  data$X = scale(data$X, center = apply(data$X[data$train_idx,], 2, mean),
                 scale = apply(data$X[data$train_idx,], 2, sd))

  # Return value
  data
}


Mnormalize = function(data) {
  # data must be a data frame with 'X' and 'Y' elements

  # Centralize Y
  data$Y = data$Y - matrix(apply(data$Y[data$train_idx,], 2, mean),ncol=3,nrow=dim(data$Y)[1],byrow=TRUE)

  # Normalize X
  data$X = scale(data$X, center = apply(data$X[data$train_idx,], 2, mean),
                 scale = apply(data$X[data$train_idx,], 2, sd))

  # Return value
  data
}


BSGSSS_EM_t = function(Y, X, group_size, niter = 100, num_update = 100, pi0 = 0.5, pi1 = 0.5,
                       alpha = 1e-1, gamma = 1e-1, a1 = 1, a2 = 1,
                       c1 = 1, c2 = 1, pi_prior = TRUE, t = 1)
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

  t_path = rep(-1, num_update)

  for (update in 1:num_update) {
    # initialize coefficients vector
    coef = array(0, dim=c(m, niter))
    bprob = array(-1, dim=c(ngroup, niter))
    tauprob = array(-1, dim=c(m, niter))
    s2_vec = rep(-1, niter)
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
      coef[,iter] = beta
      bprob[,iter] = b_prob
      tauprob[,iter] = tau_prob

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
      s2_vec[iter] = s2

    }
    # Update t
    t = 1 / (mean(1/s2_vec))
    t_path[update] = t
  }

  # output
  pos_mean = apply(coef, 1, mean)
  pos_median = apply(coef, 1, median)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef, t_path = t_path)
}



BGLSS_EM_lambda = function(Y, X, num_update = 100, niter = 100, group_size, a=1, b=1,
                           verbose = FALSE, delta=0.001, alpha=1e-1,
                           gamma=1e-1, pi_prior=TRUE, pi=0.5,option.update="global",option.weight.group=FALSE)
{
  ####################################
  # Create and Initialize parameters #
  ####################################
  n = length(Y)
  p = dim(X)[2]
  ngroup = length(group_size)
  # initialize parameters
  tau2 = rep(1, ngroup)
  sigma2 = 1
  lambda2 = 1
  matlambda2 = rep(1,ngroup)
  lambda2_path = rep(-1, num_update)
  matlambda2_path = matrix(-1,ncol=ngroup,nrow=num_update)
  l = rep(0, ngroup)
  beta = vector(mode='list', length=ngroup)
  for(i in 1:ngroup) beta[[i]]=rep(0, group_size[i])
  Z = rep(0, ngroup)

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

  for (update in 1:num_update) {
    # print(c("updadte=",update))
    coef = array(0, dim=c(p, niter))
    tau2_each_update = array(0, dim=c(ngroup, niter))
    for (iter in 1:niter)    {
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
      if (option.weight.group== FALSE){
        for(i in 1:ngroup)
        {
          if(Z[i]==0){tau2[i] = rgamma(1, shape=(group_size[i]+1)/2, rate=matlambda2[i]/2)}
          else{tau2[i] = 1/rig(1, mean=sqrt(matlambda2[i]*sigma2/sum(beta[[i]]^2)), scale = 1/(matlambda2[i]))}
        }
      }else{
        for(i in 1:ngroup)
        {
          if(Z[i]==0){tau2[i] = rgamma(1, shape=(group_size[i]+1)/2, rate=matlambda2[i]*(group_size[i])/2)}
          else{tau2[i] = 1/rig(1, mean=sqrt(matlambda2[i]*group_size[i]*sigma2/sum(beta[[i]]^2)), scale = 1/(group_size[i]*matlambda2[i]))}

        }
      }
      tau2_each_update[,iter] = tau2

      # Update sigma2
      s=0
      for(i in 1:ngroup)
      {
        s = s + sum(beta[[i]]^2)/tau2[i]
      }
      beta_vec = c()
      for(j in 1:ngroup) beta_vec = c(beta_vec, beta[[j]])
      coef[,iter] = beta_vec
      sigma2 = rinvgamma(1, shape=(n-1)/2 + sum(Z*group_size)/2 + alpha,
                         scale=(YtY-2*t(beta_vec)%*%XtY+t(beta_vec)%*%XtX%*%beta_vec+s)/2 + gamma)

      # Update pi
      if(pi_prior==TRUE)

        pi = rbeta(1, shape1=a+ngroup-sum(Z), shape2=b+sum(Z))

    }
    # Update lambda
    tau2_mean = apply(tau2_each_update, 1, mean)

    matlambda2 = (group_size+1)/(tau2_mean*group_size)

    lambda2 = (p + ngroup) / sum(tau2_mean*group_size)


    if(option.update=="global") matlambda2 <- rep(lambda2,ngroup)
    if(option.weight.group==FALSE) matlambda2 <- rep((p + ngroup) / sum(tau2_mean),ngroup)
    matlambda2_path[update,] = matlambda2



  }

  # output the posterior mean and median as our estimator
  pos_mean = apply(coef, 1, mean)
  pos_median = apply(coef, 1, median)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef, lambda2_path = matlambda2_path)
}


MBGLSS_EM_lambda = function(k,Y, X, num_update = 100, niter = 100, group_size, a=1, b=1,
                            verbose = FALSE, output=1, pi_prior=TRUE, pi=0.5 , d=3,option.update="global")
{

  ####################################
  # Create and Initialize parameters #
  ####################################
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(X)[2]
  
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

  sigma2 = 1
  lambda2 = 1
  matlambda2 = rep(1,ngroup)
  lambda2_path = rep(-1, num_update)
  matlambda2_path = matrix(-1,ncol=ngroup,nrow=num_update)


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

  for (update in 1:num_update) {
    coef = vector("list", length=q)
    for(jj in 1:q){
      coef[[jj]] <- array(0, dim=c(p, niter))
    }
    tau2_each_update = array(0, dim=c(ngroup, niter))

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


      for(i in 1:ngroup)
      {
        if(Z[i]==0){tau2[i] = rgamma(1, shape=(q*group_size[i]+1)/2, rate=group_size[i]*matlambda2[i]/2)}
        else{tau2[i] = 1/rig(1, mean=sqrt((matlambda2[i]*group_size[i])/sum(diag(beta[[i]]%*%solve(Sigma)%*%t(beta[[i]])))), scale = 1/(matlambda2[i]*group_size[i]))
        }
      }

      tau2_each_update[,iter] = tau2

      # Update Sigma  (to check value of the hyper-parameter d=3)
      ddf <- d + n +sum(group_size*Z)
      Beta <- do.call(rbind, beta)
      D_tau <- diag(rep(1/tau2,times=group_size))
      res <- Y-X%*%Beta
      S <- t(res)%*%res+t(Beta)%*%D_tau%*%Beta + Q
      Sigma <- riwish(ddf, S)

      ### save results about beta

      for(jj in 1:q){
        coef[[jj]][,iter] <- Beta[,jj]
      }

      # Update pi
      if(pi_prior==TRUE)
        pi = rbeta(1, shape1=a+ngroup-sum(Z), shape2=b+sum(Z))
    }
    # Update lambda
    tau2_mean = apply(tau2_each_update, 1, mean)

    matlambda2 = (q*group_size+1)/(tau2_mean*group_size)

    lambda2 = (p*q + ngroup) / sum(tau2_mean*group_size)
    if(option.update=="global") matlambda2 <- rep(lambda2,ngroup)
    matlambda2_path[update,] = matlambda2
  }

  # output the posterior mean and median as our estimator
  pos_mean = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, mean))),ncol=q,nrow=p,byrow=FALSE)
  pos_median = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, median))),ncol=q,nrow=p,byrow=FALSE)

  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef, lambda2_path = matlambda2_path)
}


MBSGSSS_EM_t = function(k,Y, X, group_size, niter = 100, num_update = 100, pi0 = 0.5, pi1 = 0.5,
                        alpha = 1e-1, gamma = 1e-1, a1 = 1, a2 = 1,
                        c1 = 1, c2 = 1, pi_prior = TRUE, t = 1,d=3)
{
  ####################################
  # Create and Initialize parameters #
  ####################################
  n = dim(Y)[1]
  q = dim(Y)[2]
  p = dim(X)[2]
  
  Q = k*diag(q)
  ngroup = length(group_size)
  # Initialize parameters

  Sigma <- diag(q)

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
  #YtXi = array(0,dim=c(q,m))
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

  t_path = rep(-1, num_update)

  for (update in 1:num_update) {
    # initialize coefficients vector
    coef = array(0, dim=c(m, niter))
    bprob = array(-1, dim=c(ngroup, niter))
    tauprob = array(-1, dim=c(m, niter))
    s2_vec = rep(-1, niter)
    coef = vector("list", length=q)
    for(jj in 1:q){
      coef[[jj]] <- array(0, dim=c(p, niter))
    }
    for (iter in 1:niter)
    {
      # print(c("EM",iter))
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
        ###
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
          Z[g] = 0}else
          {

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

      # store coefficients vector beta
      ### save results about beta

      for(jj in 1:q){
        coef[[jj]][,iter] <- Beta[,jj]
      }
      bprob[,iter] = b_prob
      tauprob[,iter] = tau_prob



      # Update Sigma  (to check value of the hyper-parameter d=3)
      ddf <- d + n +sum(group_size*Z)

      res <- Y-X%*%Beta

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
      s2_vec[iter] = s2

    }

    # Update t
    t = 1 / (mean(1/s2_vec))
    t_path[update] = t
  }

  # output
  # output is the posterior mean at first (iter-burnin) iteration, iter = burnin+1, ..., niter

  pos_mean = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, mean))),ncol=q,nrow=p,byrow=FALSE)
  pos_median = matrix(unlist(lapply(coef, FUN=function(x) apply(x, 1, median))),ncol=q,nrow=p,byrow=FALSE)
  list(pos_mean = pos_mean, pos_median = pos_median, coef = coef, t_path = t_path)
}





