#' MLE Parameters Trio
#'
#' This function computes the MLE of prior parameters
#' in the Gamma-Poisson model (GPM)
#'
#' @param x  The $n$ by $m$ matrix of microbial abundance.
#'
#' @return A trio of two estimates of gamma-poisson parameters and the coefficient of zero-inflated model.
#'
#' @example
#'
#'
#' @export
#' 

is_dispersion <- function(v){
  # evaluate dispersion
  # v is an m-length vector, denoting the abundance of a specific microbe
  return(var(v)>mean(v))
}

is_overdispersion <- function(v, alpha=0.05){
  # testing overdispersion
  # v is an m-length vector, denoting the abundance of a specific microbe
  test_dat <- data.frame(count = v)
  fit_test <- glm(count~., data=test_dat, family=poisson)
  disp_test <- AER::dispersiontest(fit_test, trafo=2)
  return(disp_test$p.value <= alpha)
}

is_zero_infla <- function(v, n_sample=NULL, replicates=NULL, max_prop_zero=2/3){
  # evaluate zero inflation
  # v is an m-length vector, denoting the abundance of a specific microbe
  replicates <- ifelse(is.null(replicates), 3, replicates)
  n_sample <- ifelse(is.null(n_sample), length(v)/replicates, n_sample)
  record <- FALSE
  max_zero <- replicates * max_prop_zero
  sudo_m <- diag(n_sample)
  for(k in 1:ncol(sudo_m)){
    sudo_v <- rep(sudo_m[,k], each=replicates)
    num_0  <- sum((v * sudo_v) == 0) - sum(sudo_v==0)
    record <- ifelse(num_0>max_zero,TRUE,record)
  }
  return(record)
}


calc_mle_gpm <- function(x, n_sample=NULL, replicates=NULL, max_prop_zero=2/3){
  param_trio <- data.frame()
  for(i in 1:nrow(x)){
    if (is_zero_infla(as.numeric(x[i,]), n_sample, replicates, max_prop_zero=2/3)){
      # fit zero-inflated model
      test_dat <- data.frame(count = as.numeric(x[i,]))
      tryCatch({      
        fit_test <- pscl::zeroinfl(count~1|1, link='logit', dist='negbin',data=test_dat)
        mu  <- exp(fit_test$coefficients$count[['(Intercept)']])
        k   <- fit_test$theta
        logit_pi <- fit_test$coefficients$zero[['(Intercept)']]
        pi0 <- exp(logit_pi)/(1+exp(logit_pi))
        }, error=function(e){
          cat("ERROR :",conditionMessage(e), "when analyzing microbe ", rownames(x)[i], "\n")
          })
    } else if (is_overdispersion(as.numeric(x[i,]), alpha=0.05)) {
      # fit negative binomial 
      test_dat <- data.frame(count = as.numeric(x[i,]))
      tryCatch({
        fit_test <- VGAM::vglm(count~1, negbinomial, data=test_dat)
        mu  <- exp(fit_test@coefficients[['(Intercept):1']])
        k   <- exp(fit_test@coefficients[['(Intercept):2']])
        pi0 <- 0
        
        # an alternative function to fit negative binomial
        #fit_test <- MASS::glm.nb(count~1, data=test_dat, 
        #                         control=glm.control(maxit=100,trace=FALSE))
        #mu  <- exp(fit_test$coefficients[['(Intercept)']])
        #k   <- fit_test$theta
        #pi0 <- 0
        
        }, error=function(e){
          cat("ERROR :",conditionMessage(e), "when analyzing microbe ", rownames(x)[i], "\n")
          })
    } else {
      # fit poisson
      test_dat <- data.frame(count = as.numeric(x[i,]))
      fit_test <- glm(count~1, data=test_dat, family=poisson(), trace=FALSE)
      mu  <- exp(fit_test$coefficients[['(Intercept)']])
      k   <- Inf
      pi0 <- 0
    }
    param_trio <- rbind(param_trio,
                        data.frame(mu  = mu,
                                   k   = k,
                                   pi0 = pi0))
    print(paste0('process ',round(i/nrow(x)*100, digits=2),'%' ))
  }
  rownames(param_trio) <- rownames(x)
  print('completed!!!')
  return(param_trio)
}

calc_ll_poisson <- function(v, theta){
  return(sum(- theta + v * log(theta) - lgamma(v + 1)))
}

id_overdispersion <- function(x){
  is_od <- apply(x, 1, function(v) is_overdispersion(as.numeric(v), alpha=0.05))
  id_od <- rownames(x)[is_od]
  id_nod <- rownames(x)[!is_od]
  lis_out <- list(id_od, id_nod)
  names(lis_out) <- c('od','nod')
  return(lis_out)
}

ll_mu_k <- function(x){
  # x is n by m abundance matrix
  p_value <- c()
  for(i in 1:nrow(x)){
    test_dat <- data.frame(count = as.numeric(x[i,]))
    # is zero-inflated or not
    if (is_zero_infla(test_dat$count)){
      # fit zero-inflated poisson
      fit_zp <- pscl::zeroinfl(count~1|1, link='logit', dist='poisson',data=test_dat)
      ll_zp <- fit_zp$loglik
      # fit zero-inflated gamma poisson
      fit_zgp <- pscl::zeroinfl(count~1|1, link='logit', dist='negbin',data=test_dat)
      ll_zgp <- fit_zgp$loglik
      # h0 prefer zero-inflated poisson
      # h0 prefer k-wing
      # mu-wing if p_value less than 0.05
      p_lrt <- pchisq(2*(ll_zgp-ll_zp), df = 3-2, lower.tail = FALSE)
      p_value <- c(p_value, p_lrt)
      
    } else {
      # test over-dispersion
      # h0 prefer non over-dispersion
      # h0 prefer k-wing
      # mu-wing if p_value less than 0.05
      fit_test <- glm(count~., data=test_dat, family=poisson)
      disp_test <- AER::dispersiontest(fit_test, trafo=2)
      p_value <- c(p_value, disp_test$p.value)
    }
  }
  return(p_value)
}

classify_mu_k <- function(x){
  # x is n by m abundance matrix
  id_mu <- c()
  id_k  <- c()
  for(i in 1:nrow(x)){
    
    test_dat <- data.frame(count = as.numeric(x[i,]))
    # over-dispersion test
    fit_p <- glm(count~., data=test_dat, family=poisson)
    disp_test <- AER::dispersiontest(fit_p, trafo=2)
    p_od  <- disp_test$p.value
    # k coefficient test
    fit_gp <- VGAM::vglm(count~1, negbinomial, data=test_dat)
    tbl_gp <- summary(fit_gp)
    p_gpk <- tbl_gp@coef3['(Intercept):2','Pr(>|z|)']
    
    if((p_od<=0.05)&(p_gpk<=0.05)){
      id_mu <- c(id_mu, rownames(x)[i])
    } else if ((p_od>0.05)&(p_gpk>0.05)) {
      id_k  <- c(id_k,  rownames(x)[i])
    } else {
      # is zero-inflated or not
      if (is_zero_infla(test_dat$count)){
        # fit zero-inflated poisson
        fit_zp <- pscl::zeroinfl(count~1|1, link='logit', dist='poisson',data=test_dat)
        ll_zp <- fit_zp$loglik
        # fit zero-inflated gamma poisson
        fit_zgp <- pscl::zeroinfl(count~1|1, link='logit', dist='negbin',data=test_dat)
        ll_zgp <- fit_zgp$loglik
        # h0 prefer zero-inflated poisson
        # h0 prefer k-wing
        # mu-wing if p_value less than 0.05
        p_lrt <- pchisq(2*(ll_zgp-ll_zp), df = 3-2, lower.tail = FALSE)
        if(p_lrt <= 0.05) {
          id_mu <- c(id_mu, rownames(x)[i])
        } else {
          id_k  <- c(id_k,  rownames(x)[i])
        }
      } else {
        # fit gamma poisson
        fit_gp <- VGAM::vglm(count~1, negbinomial, data=test_dat)
        ll_gp <- fit_gp@criterion$loglikelihood
        # fit poisson
        fit_p <- glm(count~1, data=test_dat, family=poisson(), trace=FALSE)
        mu  <- exp(fit_p$coefficients[['(Intercept)']])
        ll_p <- calc_ll_poisson(test_dat$count, theta=mu)
        # h0 prefer poisson
        # h0 prefer k-wing
        # mu-wing if p_value less than 0.05
        p_lrt <- pchisq(2*(ll_gp-ll_p), df = 2-1, lower.tail = FALSE)
        if(p_lrt <= 0.05) {
          id_mu <- c(id_mu, rownames(x)[i])
        } else {
          id_k  <- c(id_k,  rownames(x)[i])
        }
      }
    }
  }
  lis_out <- list(id_mu, id_k)
  names(lis_out) <- c('mu','k')
  return(lis_out)
}



