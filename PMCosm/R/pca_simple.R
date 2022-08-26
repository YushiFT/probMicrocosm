#' Principle Components Analysis
#'
#' This function conducts PCA for microbial data set
#'
#' @param x  The $n$ by $m$ matrix of microbial abundance with columns denoting sample sites.
#'
#' @return A list of principle components, loadings, and proportion of variance explained by each pc.
#'
#' @example
#'
#'
#' @export
#' 

pca <- function(x, center=TRUE, scale=FALSE) {
  x <- t(scale(t(x), center=center, scale=scale))
  x <- x/sqrt(nrow(x)-1)
  s <- svd(x)
  loading <- s$u
  colnames(loading) <- paste0("Loading", 1:ncol(loading))
  rownames(loading) <- rownames(x)
  pc <- diag(s$d) %*% t(s$v)
  rownames(pc) <- paste0("PC", 1:nrow(pc))
  colnames(pc) <- colnames(x)
  pve <- s$d^2 / sum(s$d^2)
  return(list(pc=pc, loading=loading, pve=pve))
  }

