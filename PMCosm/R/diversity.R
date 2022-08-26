#' Diversity Evaluation
#'
#' A set of functions that assist to calculate diversity indeces including
#' alpha diversity index
#' gamma diversity index
#' beta diversity index
#' 
#'
#' @param x  The $n$ by $m$ matrix of microbial abundance.
#'
#' @return 
#'
#' @example
#'
#'
#' @export
#' 

count_to_preabs <- function(x){
  # transform counting data to presence-absence data
  f_cell <- function(x_ij){
    return(ifelse(x_ij>0, 1, 0))
  }
  return(apply(x, c(1,2), f_cell))
}

calc_commu_beta <- function(x){
  # describe beta diversity at the community level
  # x is presence-absence matrix
  b_out <- data.frame()
  sum_stat <- function(v){
    u=round(mean(v)+1.96*sd(v)/sqrt(length(v)), digits=3)
    l=round(mean(v)-1.96*sd(v)/sqrt(length(v)), digits=3)
    dat_out <- data.frame(m=round(mean(v),digits=3),
                          ci=paste0('(',l,',',u,')'))
    return(dat_out)
  }
  b_lis <- c('w','-1','c','wb','r','I',
             'e','t','me','m','-2','co',
             'cc','g','-3','l','19','hk',
             'sim','z')
  for(b in b_lis){
    b_out <- rbind(b_out, sum_stat(as.vector(vegan::betadiver(x, b))))
  }
  rownames(b_out) <- paste0('beta_',b_lis)
  return(b_out)
}


