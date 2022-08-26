#' Annotate ribosome RNA operons
#'
#' This function annotates taxa with rrn count
#'
#' @param tax  The $n$ by $3$ matrix of taxa id, family, and genus.
#'
#' @return A vector of mean value of rrn copy number.
#'
#' @example
#'
#'
#' @export
#' 

rrn_anno <- function(tax, rrn){
  rrn$anno <- rrn$name
  rrn_genus  <- rrn[rrn$rank=='genus',]
  rrn_family <- rrn[rrn$rank=='family',]
  rrn_order  <- rrn[rrn$rank=='order',]
  rrn_class  <- rrn[rrn$rank=='class',]
  rrn_phylum <- rrn[rrn$rank=='phylum',]
  
  tax_genus  <- tax[tax$G!="", ]
  tax_family <- tax[(tax$G=="") & (tax$F!=""), ]
  tax_order  <- tax[(tax$G=="") & (tax$F=="") & (tax$O!=""), ]
  tax_class  <- tax[(tax$G=="") & (tax$F=="") & (tax$O=="") & (tax$C!=""), ]
  tax_phylum <- tax[(tax$G=="") & (tax$F=="") & (tax$O=="") & (tax$C=="") & (tax$P!=""), ]
  
  tax_genus$anno  <- tax_genus$G
  tax_family$anno <- tax_family$F
  tax_order$anno  <- tax_order$O
  tax_class$anno  <- tax_class$C
  tax_phylum$anno  <- tax_phylum$P
  
  tax_genus  <- tax_genus[,c('ID','anno')]
  tax_family <- tax_family[,c('ID','anno')]
  tax_order  <- tax_order[,c('ID','anno')]
  tax_class  <- tax_class[,c('ID','anno')]
  tax_phylum <- tax_phylum[,c('ID','anno')]
  
  out_genus  <- merge(tax_genus,  rrn_genus,  by='anno')
  out_family <- merge(tax_family, rrn_family, by='anno')
  out_order  <- merge(tax_order,  rrn_order,  by='anno')
  out_class  <- merge(tax_class,  rrn_class,  by='anno')
  out_phylum <- merge(tax_phylum, rrn_phylum, by='anno')
  
  out_rrn <- rbind(out_genus[,c('ID','mean')], out_family[,c('ID','mean')],
                   out_order[,c('ID','mean')], out_class[,c('ID','mean')],
                   out_phylum[,c('ID','mean')])
  return(out_rrn)
}

calc_comm_rrn <- function(rrn_anno){
  # calculate community-level rrn copy number
  # the first two columns of rrn_anno are
  # ID, Mean
  # return a vector with length of the community size
  top <- colSums(rrn_anno[,3:ncol(rrn_anno)])
  bot <- colSums(rrn_anno[,3:ncol(rrn_anno)]/rrn_anno$mean)
  return(top/bot)
}

