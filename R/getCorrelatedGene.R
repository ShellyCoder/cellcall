#' get relation among specific TF
#' @param data a dataframe with row of gene and column of sample
#' @param cell_type specific cell type
#' @param tf choose one specific TF
#' @param target_list choose target genes corresponding TF
#' @param pValue set the pValue of spearman
#' @param corValue set the corValue of spearman
#' @param topGene use topTargetCor of candidate genes which has firlter by above parameters, default is 1, means 100%
#' @return result target genes in this regulon
#' @importFrom psych corr.test

getCorrelatedGene <- function(data, cell_type="", tf, target_list,
  pValue=0.05, corValue=0, topGene=0.1
){
  my_Expr <- data
  expr_tmp <- my_Expr[,which(colnames(my_Expr)==cell_type)]
  if(tf %in% rownames(expr_tmp)){
    # library(psych)
    p_tmp<-psych::corr.test(t(expr_tmp[tf,]),t(expr_tmp[target_list,]),adjust = "none",use = "pairwise", method = "spearman")
    tmp_pValue = p_tmp$p
    tmp_corr = p_tmp$r
    tmp_pValue[is.na(tmp_pValue)]=1  ## NA是两gene有标准差为零的情况，因此设P为1，即不显著
    tmp_corr[is.na(tmp_corr)]=0  ## NA是两gene有标准差为零的情况，因此设P为0，即不显著
    corr_gene_tmp <- as.data.frame(t(tmp_corr[1,tmp_pValue<pValue,drop=F]))

    if(nrow(corr_gene_tmp)>1){
      corr_gene_tmp <- corr_gene_tmp[corr_gene_tmp>corValue,1,drop=F]
      corr_target_tmp <- rownames(corr_gene_tmp[order(corr_gene_tmp[,1],decreasing = T),1,drop=F])
      topTarget = floor(topGene * length(corr_target_tmp))

      if(length(corr_target_tmp)>(topTarget+1)){
        # print(topTarget)
        corr_target_tmp <- corr_target_tmp[2:(topTarget+1)]
      }
    }else{
      corr_target_tmp <- vector()
    }
  }else{
    corr_target_tmp <- vector()
  }
  return(corr_target_tmp)
}





