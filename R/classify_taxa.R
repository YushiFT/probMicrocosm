#' Categories of microbial taxa
#'
#' This function classify microbial taxa into six categories
#' RT: rare taxa 
#' MT: moderate taxa
#' AT: abundant taxa
#' CRT: conditionally rare taxa
#' CAT: conditionally abundant taxa
#' CRAT: conditionally rare or abundant taxa
#'
#' @param x  The $n$ by $m$ matrix of microbial abundance.
#'
#' @return A vector of microbial taxa belonging to each category.
#'
#' @example
#'
#'
#' @export
#' 

is_rt <- function(x, t_r=0.001){
  # extract rare taxa
  # rare taxa has relative abundance no greater than 0.01% in all samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_rt <- c()
  f_sub <- function(v){
    record_rt <- sum(v<=t_r)
    return(record_rt==length(v))
  }
  is_rt_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_rt_vec==TRUE])
}

is_mt <- function(x, t_r=0.001, t_a=0.01){
  # extract moderate taxa
  # moderate taxa has relative abundance between 0.01% and 1% in all samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_mt <- c()
  f_sub <- function(v){
    record_mt <- sum((v>t_r)&(v<=t_a))
    return(record_mt==length(v))
  }
  is_mt_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_mt_vec==TRUE])
}

is_at <- function(x, t_a=0.01){
  # extract abundant taxa
  # abundant taxa has relative abundance greater than 1% in all samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_at <- c()
  f_sub <- function(v){
    record_at <- sum(v>t_a)
    return(record_at==length(v))
  }
  is_at_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_at_vec==TRUE])
}

is_crt <- function(x, t_r=0.001, t_a=0.01){
  # extract conditionally rare taxa
  # conditionally rare taxa has relative abundance no greater than 1% in all samples 
  # and no greater than 0.01% in some samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_crt <- c()
  f_sub <- function(v){
    record_1 <- sum(v<=t_a)
    record_2 <- sum((v>t_r)&(v<=t_a))
    return( (record_1==length(v))&(record_2>0)&(record_2<length(v)) )
  }
  is_crt_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_crt_vec==TRUE])
}

is_cat <- function(x, t_r=0.001, t_a=0.01){
  # extract conditionally abundant taxa
  # conditionally abundant taxa has relative abundance greater than 0.01% in all samples 
  # and greater than 1% in some samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_cat <- c()
  f_sub <- function(v){
    record_1 <- sum(v>t_r)
    record_2 <- sum(v>t_a)
    return( (record_1==length(v))&(record_2>0)&(record_2<length(v)) )
  }
  is_cat_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_cat_vec==TRUE])
}

is_crat <- function(x, t_r=0.001, t_a=0.01){
  # extract conditionally rare or abundant taxa
  # conditionally rare or abundant taxa has relative abundance varying from no greater than 0.01% 
  # to greater than 1%
  x_prop <- x/colSums(x) # transform to relative abundance
  id_crat <- c()
  f_sub <- function(v){
    record_1 <- sum(v<=t_r)
    record_2 <- sum(v>t_a)
    return( (record_1>0)&(record_2>0) )
  }
  is_crat_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_crat_vec==TRUE])
}

is_generalist <- function(x){
  # generalists with B > b_g
  # moderate generalists with B between b_s and b_g
  x_prop <- x/colSums(x) # transform to relative abundance
  mean_p <- rowMeans(x_prop)
  x_sub <- x[mean_p >= 2e-5,]
  p <- x_sub/rowSums(x_sub)
  b_j <- 1/rowSums(p^2)
  q1 <- as.numeric(quantile(b_j)[2])
  q2 <- as.numeric(quantile(b_j)[3])
  q3 <- as.numeric(quantile(b_j)[4])
  b_g <- q3 + 1.5*(q3-q1)
  return(rownames(x_sub)[b_j>b_g])
}

is_specialists <- function(x, b_s=1.2){
  # specialists with B < b_s
  # generalists with B > b_g
  # moderate generalists with B between b_s and b_g
  x_prop <- x/colSums(x) # transform to relative abundance
  mean_p <- rowMeans(x_prop)
  x_sub <- x[mean_p >= 2e-5,]
  p <- x_sub/rowSums(x_sub)
  b_j <- 1/rowSums(p^2)
  return(rownames(x_sub)[b_j<b_s])
}


cal_nicheb <- function(x){
  # specialists with B < b_s
  # generalists with B > b_g
  # moderate generalists with B between b_s and b_g
  p <- x/rowSums(x)
  b_j <- 1/rowSums(p^2)
  return(b_j)
}

###########################################
# so many categorizing rules
# fluctuate the thresholds
###########################################
# method a ##################################################
is_rt_a <- function(x){
  # extract rare taxa according to method a
  # rare taxa has relative abundance less than 1% in all samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_rt <- c()
  f_sub <- function(v){
    record_rt <- sum(v<0.01)
    return(record_rt==length(v))
  }
  is_rt_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_rt_vec==TRUE])
}

is_at_a <- function(x){
  # extract rare taxa according to method a
  # rare taxa has relative abundance no less than 1% in any sample
  x_prop <- x/colSums(x) # transform to relative abundance
  id_at <- c()
  f_sub <- function(v){
    record_at <- sum(v>=0.01)
    return(record_at>0)
  }
  is_at_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_at_vec==TRUE])
}

# method b ####################################################
is_rt_b <- function(x){
  # extract regional rare taxa according to method b
  # regional rare taxa has mean relative abundance less than 0.001%
  x_prop <- x/colSums(x) # transform to relative abundance
  v_prop <- rowMeans(x_prop)
  return(rownames(x)[v_prop<0.00001])
}

is_at_b <- function(x){
  # extract regional abundant taxa according to method b
  # regional abundant taxa has mean relative abundance greater than 0.1%
  x_prop <- x/colSums(x) # transform to relative abundance
  v_prop <- rowMeans(x_prop)
  return(rownames(x)[v_prop>0.001])
}

is_it_b <- function(x){
  # extract regional intermediate taxa according to method b
  # regional intermediate taxa has mean relative abundance between 0.001% and 0.1%
  x_prop <- x/colSums(x) # transform to relative abundance
  v_prop <- rowMeans(x_prop)
  return(rownames(x)[(v_prop>=0.00001)&(v_prop<=0.001)])
}

# method c ##################################################
is_rt_c <- function(x){
  # extract rare taxa according to method a
  # rare taxa has relative abundance less than 0.1% in all samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_rt <- c()
  f_sub <- function(v){
    record_rt <- sum(v<0.001)
    return(record_rt==length(v))
  }
  is_rt_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_rt_vec==TRUE])
}

is_at_c <- function(x){
  # extract rare taxa according to method a
  # abundant taxa has relative abundance no less than 0.1% in half or more samples
  # and occurred in more than 80% samples
  x_prop <- x/colSums(x) # transform to relative abundance
  id_at <- c()
  f_sub <- function(v){
    record_at_1 <- sum(v>=0.001)
    record_at_2 <- sum(v>0)
    return((record_at_1>0.5*length(v)) & (record_at_2>0.8*length(v)))
  }
  is_at_vec <- apply(x_prop, 1, f_sub)
  return(rownames(x)[is_at_vec==TRUE])
}





