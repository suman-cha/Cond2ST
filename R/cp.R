# =============================================================================
# R/cp.R -- Helpers for the classifier-based two-sample test (CP / DCP).
#
# Internal building blocks consumed by debiased_test, CP_test, and DCP_test
# in R/hypothesis_tests.R. Not intended to be called directly by experiment
# scripts.
# =============================================================================

## functions

getUstat <- function(g1, v1, v2){
  ## to calculate \sum_j \hat{U}_j/n21

  n1 <- length(g1); n2 <- length(v2)
  zeta <- runif(n2)
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1
  inner_sum <- sapply(v1, function(x){mean(zeta*(x==v2))})
  u1 <- mean(g1*((1-ecdf1)+inner_sum))
  return(u1/mean(g1))
}

## calculate the empirical cdf under 1st sample
getEcdf1 <- function(v1, v2){
  fun <- ecdf(v2)
  ecdf1 <- fun(v1)
  ecdf2 <- sapply(v1, function(x){return(mean(x==v2))})
  return(list(ecdf1=ecdf1, ecdf2=ecdf2))
}

## calculate the empirical cdf under 2nd sample
getEcdf2 <- function(v){
  fun <- ecdf(v)
  ecdf1 <- fun(v)
  ecdf2 <- sapply(v, function(x){return(mean(x==v))})
  return(list(ecdf1=ecdf1, ecdf2=ecdf2))
}

## calculate the asymptotic variance under sample 2 (cancellation one g)
## ecdf1: P(V_2<=v), ecdf2: P(V_2=v)
getAsyVar2 <- function(g2, v2, m1){
  m2 <- length(g2)
  ecdfVal <- getEcdf2(v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  s1 <- mean(g2*(1-ecdf1)^2)-1/4
  s2 <- mean(g2*ecdf2^2)/3
  s3 <- mean(g2*(1-ecdf1)*ecdf2)
  ss <- s1+s2+s3
  ss <- ifelse(ss>0, ss, 0)
  sigma1 <- ss/m1 + 1/12/m2
  sigma2 <- (mean(g2)-1)/m1
  sigma2 <- ifelse(sigma2>0, sigma2, 0)
  sigma12 <- (mean(g2*(1-ecdf1))+mean(g2*ecdf2)/2 - 1/2)/m1
  sigma <- sigma1 + sigma2/4 - sigma12
  return(list(sigma=sigma, sigma1=sigma1, sigma2=sigma2, sigma12=sigma12))
}


## calculate the asymptotic variance under 1st sample
getAsyVar1 <- function(g1, v1, v2){
  m1 <- length(g1); m2 <- length(v2)
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  ss <- mean(g1^2*(1-ecdf1)^2)+mean(g1^2*ecdf2^2)/3+mean(g1^2*(1-ecdf1)*ecdf2)-
    c(1/4, mean(g1*(1-ecdf1))+1/2*mean(g1*ecdf2)-1/4, (mean(g1*(1-ecdf1))+mean(g1*ecdf2)/2)^2)
  ss <- ifelse(ss>0, ss, 0)
  sigma1 <- ss/m1 + 1/12/m2
  sigma2 <- (mean(g1^2)-c(1, 2*mean(g1)-1, mean(g1)^2))/m1
  sigma2 <- ifelse(sigma2>0, sigma2, 0)
  sigma12 <- (mean(g1^2*(1-ecdf1))+mean(g1^2*ecdf2)/2 -
                c(1/2, mean(g1)/2+mean(g1*(1-ecdf1)+g1*ecdf2/2)-1/2, mean(g1)*mean(g1*(1-ecdf1)+g1*ecdf2/2)))/m1
  sigma <- sigma1 + sigma2/4 - sigma12
  return(list(sigma=sigma, sigma1=sigma1, sigma2=sigma2, sigma12=sigma12))
}

## calculate the final statistic
getFinalStat <- function(g1, g2, v1, v2){
  m1 <- length(g1); m2 <- length(g2)
  U <- getUstat(g1, v1, v2)
  sigma.1 <- getAsyVar1(g1, v1, v2)$sigma # length: 3
  sigma.2 <- getAsyVar2(g2, v2, m1)$sigma
  sigma.gm <- sqrt(sigma.1*sigma.2)
  sigma.hm <- 2/(1/sigma.1+1/sigma.2)
  z1 <- (U-0.5)/sqrt(sigma.1)
  z2 <- (U-0.5)/sqrt(sigma.2)
  z.gm <- (U-0.5)/sqrt(sigma.gm)
  z.hm <- (U-0.5)/sqrt(sigma.hm)
  return(list(U=U, z1=z1, z2=z2, z.gm=z.gm, z.hm=z.hm, sigma.1=sigma.1, sigma.2=sigma.2,
              sigma.gm=sigma.gm, sigma.hm=sigma.hm))
}

NNfun <- function(x, z, xpre1, xpre2, nnrep=10, hidden.layers=NA, acfun = 'sigmoid', optim.type='sgd',
                  n.epochs=500, learn.rates=0.001, L1=0){

  n1 <- dim(xpre1)[1]; n2 <- dim(xpre2)[1]
  prob1 <- matrix(0, nrow = n1, ncol = nnrep)
  prob2 <- matrix(0, nrow = n2, ncol = nnrep)

  for(i in 1:nnrep){
    # print(paste0('Neural Net:', i))
    fit_nn <- neuralnetwork(x, z, hidden.layers = hidden.layers, optim.type = optim.type,
                            val.prop=0, learn.rates = learn.rates, L1=L1,
                            n.epochs = n.epochs, activ.functions = acfun, verbose = 0)

    prob1[,i] <- predict(fit_nn, xpre1)$probabilities[,2]
    prob2[,i] <- predict(fit_nn, xpre2)$probabilities[,2]
  }

  prob1.fit <- rowMeans(prob1)
  prob2.fit <- rowMeans(prob2)
  prob1.fit[prob1.fit<0.01] <- 0.01; prob1.fit[prob1.fit>0.99] <- 0.99
  prob2.fit[prob2.fit<0.01] <- 0.01; prob2.fit[prob2.fit>0.99] <- 0.99

  return(list(prob1.fit=prob1.fit, prob2.fit=prob2.fit))
}


KLR.CV <- function(x, y, xpre, ypre, lambdaseq, sigmaseq){
  lambda_len <- length(lambdaseq); sigma_len <- length(sigmaseq)
  centropy <- matrix(0, lambda_len, sigma_len)
  ypre <- as.numeric(as.character(ypre))
  minerror <- Inf
  for(i in 1:lambda_len){
    lambda <- lambdaseq[i]
    for(j in 1:sigma_len){
      print(paste0('KLR CV: lambda', i, 'sigma', j))
      sigma <- sigmaseq[j]
      cv_data <- constructData(x, y)
      cv_data <- shuffleData(cv_data)
      klr_learner <- constructKlogRegLearner()
      # rbf: exp(-sigma\|x-y\|^2)
      params <- list(kernel='rbfdot', sigma=sigma, lambda=lambda/getN(cv_data), tol=10e-6, maxiter=500)
      fit_klr <- klr_learner$learn(cv_data, params)
      K = kernelMult(fit_klr$kernel, xpre, fit_klr$data, fit_klr$alpha) # predict
      prob_pre = 1 / (1 + exp(-as.vector(K))) # predicted probabilities
      # cross entropy error
      prob_pre[prob_pre<0.01] <- 0.01; prob_pre[prob_pre>0.99] <- 0.99
      centropy[i, j] <- mean(-ypre*log(prob_pre)-(1-ypre)*log(1-prob_pre))
      if(centropy[i,j]<minerror){
        minerror <- centropy[i,j]
        iind <- i; jind <- j
        prob <- prob_pre
      }
    }
  }
  prob[prob<0.01] <- 0.01; prob[prob>0.99] <- 0.99
  lambda <- lambdaseq[iind]; sigma <- sigmaseq[jind]
  return(list(prob=prob, lambda=lambda, sigma=sigma, centropy=centropy))

}

## calculate the cross entropy error
## y:label, 0 or 1; p: estimated probability for label 1
centropy <- function(p, y){
  return(mean(-y*log(p)-(1-y)*log(1-p)))
}

gerror <- function(g1, g2){
  sum(abs(g1/sum(g1)-g2/sum(g2)))
}

## calculate relevant quantities in assumptions
getCor <- function(g1,v1,v2,gorac1,vorac1,vorac2){
  mean_hatg_minus_g <- mean(g1-gorac1)
  mean_abs_hatg_minus_g <- mean(abs(g1-gorac1))
  ecdfVal <- getEcdf1(v1, v2)
  ecdf1 <- ecdfVal$ecdf1; ecdf2 <- ecdfVal$ecdf2
  mean_hatd <- mean(1-ecdf1 + ecdf2/2)
  var_hatg_minus_g <- var(g1-gorac1)
  second_hatd <- mean(1-ecdf1+ecdf2/3)
  var_hatd <- second_hatd - mean_hatd^2
  second_hatg_minus_g_hatd <- mean((g1-gorac1)*((1-ecdf1)+ecdf2/2))

  e1 <- second_hatg_minus_g_hatd - mean_hatg_minus_g*mean_hatd # covariance between \hat G- G and \hat D
  e2 <- mean_hatg_minus_g
  e3 <- mean(gorac1*((1-ecdf1)+ecdf2/2)) - mean_hatd # covariance between G and \hat D

  # correlation between \hat G- G and \hat D
  rho1 <- e1 / sqrt(var_hatg_minus_g*var_hatd)
  rho2 <- mean_hatg_minus_g/mean_abs_hatg_minus_g
  # correlation between G and \hat D
  var_g <- var(gorac1)
  rho3 <- e3/sqrt(var_hatd*var_g)

  ecdfValorac <- getEcdf1(vorac1, vorac2)
  ecdforac1 <- ecdfValorac$ecdf1; ecdforac2 <- ecdfValorac$ecdf2
  e4 <- mean(gorac1*((1-ecdf1)+ecdf2/2)) - mean(gorac1*((1-ecdforac1)+ecdforac2/2)) # mean of G(\hat D - D)

  return(c('rho1'=rho1, 'rho2'=rho2, 'rho3'=rho3, 'e1'=e1, 'e2'=e2, 'e3'=e3, 'e4'=e4))
}
