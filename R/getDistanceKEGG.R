#' abstract path from kegg pathway file
#' @param data a dataframe store L-R-TF relation in KEGG and the distance between LR and TF
#' @param method the method to calculate multiple path corresponding specific L-R-TF, default is "mean"
#' @return the distance between LR and TF
#' @importFrom stringr str_split
#' @importFrom stats median quantile
#' @importFrom magrittr %>%

getDistanceKEGG <- function(data, method){
  # library(stringr)
  triple_relation <- data
  l_r_inter <- unique(triple_relation[,5:6])
  dist_lr_tf <- matrix(data = 0, nrow = nrow(l_r_inter), ncol = length(unique(triple_relation$TF_Symbol)))
  dist_lr_tf <- as.data.frame(dist_lr_tf)
  rownames(dist_lr_tf) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames(dist_lr_tf) <- unique(triple_relation$TF_Symbol)

  for (i in 1:nrow(triple_relation)) {
    index_row_tmp <- paste(triple_relation[i,5],triple_relation[i,6],sep = "-")
    index_col_tmp <- triple_relation[i,7]
    tmp_mat <- str_split(triple_relation[i,4],",")[[1]] %>% str_split('_',simplify = T)
    if(method=="mean"){
      mean_dist <- mean(as.numeric(tmp_mat[,2]))
    }else if(method=="median"){
      mean_dist <- median(as.numeric(tmp_mat[,2]))
    }
    dist_lr_tf[index_row_tmp, index_col_tmp] <- mean_dist
  }
  return(dist_lr_tf)
}





