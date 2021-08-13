#' get score of Triple relation corresponding specific cellA-cellB
#' @param object a Cellwave objects
#' @param sender_cell the cell type of sender cell
#' @param recevier_cell the cell type of recevier cell
#' @param slot plot the graph with the data of specific slot
#' @param org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @importFrom utils read.table head
#' @importFrom stringr str_split
#' @importFrom dplyr filter
#' @export

LR2TF <- function(object, sender_cell, recevier_cell, slot="expr_l_r_log2_scale", org="Homo sapiens", IS_core=TRUE){
  options(stringsAsFactors = F) 

  myData <- object@data[[slot]]

  detect_gene <- rownames(object@data$expr_mean)
  regulons_matrix <- object@data$regulons_matrix
  DistanceKEGG <- object@data$DistanceKEGG
  expr_mean <- object@data$expr_mean
  sender_cell <- sender_cell
  recevier_cell <- recevier_cell
  # top_n <- top_n

  if(org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs.txt", package="cellcall")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)

    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_extended.txt", package="cellcall")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }

    f.tmp <- system.file("extdata", "tf_target.txt", package="cellcall")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology.txt", package="cellcall")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)

    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology_extended.txt", package="cellcall")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }

    f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellcall")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }

  interest_group <- myData[,paste(sender_cell, recevier_cell, sep = "-"), drop=F]
  interest_df <- interest_group[order(interest_group,decreasing = T),,drop=F]
  interest_df <- interest_df[which(interest_df[,1]!=0),,drop=F]

  # library(stringr)
  my_ligand_receptor_df <- str_split(rownames(interest_df),"-",simplify = T)
  my_tfSet <- c()
  for (n in 1:nrow(interest_df)) {
    sender_tmp <- my_ligand_receptor_df[n,1]
    receiver_tmp <- my_ligand_receptor_df[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")

    info_tmp <- dplyr::filter(triple_relation, Ligand_Symbol==sender_tmp & Receptor_Symbol==receiver_tmp)[,c("Ligand_Symbol", "Receptor_Symbol", "TF_Symbol")]
    tfs_tmp <- info_tmp$TF_Symbol[info_tmp$TF_Symbol %in% detect_gene]
    tfs_tmp <- tfs_tmp[regulons_matrix[tfs_tmp, recevier_cell]>0]
    my_tfSet <- c(my_tfSet, tfs_tmp)
  }
  my_tfSet <- unique(my_tfSet)

  sankey_matrix <- matrix(0, ncol = 5)
  l_r_tf <- matrix(0, ncol = length(my_tfSet), nrow = nrow(interest_df))
  l_r_tf <- data.frame(l_r_tf)
  rownames(l_r_tf) <- rownames(interest_df)
  colnames(l_r_tf) <- my_tfSet
  for (n in 1:nrow(interest_df)) {
    print(n)
    sender_tmp <- my_ligand_receptor_df[n,1]
    receiver_tmp <- my_ligand_receptor_df[n,2]
    sender_val <- expr_mean[sender_tmp,sender_cell]
    receiver_val <- expr_mean[receiver_tmp,recevier_cell]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")

    info_tmp <- dplyr::filter(triple_relation, Ligand_Symbol==sender_tmp & Receptor_Symbol==receiver_tmp)[,c("Ligand_Symbol", "Receptor_Symbol", "TF_Symbol")]
    my_tfSet_tmp <- info_tmp$TF_Symbol[info_tmp$TF_Symbol %in% detect_gene]
    my_tfSet_tmp <- my_tfSet_tmp[regulons_matrix[my_tfSet_tmp, recevier_cell]>0]

    tfVal_tmp <- regulons_matrix[my_tfSet_tmp,recevier_cell]
    distance2w_tmp <- (1/DistanceKEGG[row_index,my_tfSet_tmp])
    w_tmp<- distance2w_tmp/sum(distance2w_tmp)
    l_r_val_tmp <- myData[row_index,paste(sender_cell, recevier_cell, sep = "-")]

    # l_r_tf[row_index,my_tfSet_tmp] <- as.numeric(tfVal_tmp * w_tmp * l_r_val_tmp)
    sankey_tmp <- str_split( paste(sender_tmp, receiver_tmp, my_tfSet_tmp, l_r_val_tmp, tfVal_tmp * w_tmp )," ",simplify = TRUE)
    sankey_matrix <- rbind(sankey_matrix, sankey_tmp)

  }
  # l_r_tf <- l_r_tf[,apply(l_r_tf, 2, function(x){sum(x!=0)})>0] 
  sankey_matrix <- data.frame(sankey_matrix)
  colnames(sankey_matrix) <- c("Ligand", "Receptor",	"TF",	"weight1",	"weight2")
  sankey_matrix <- sankey_matrix[-1,]
  sankey_matrix$weight1 <- as.numeric(sankey_matrix$weight1)
  sankey_matrix$weight2 <- as.numeric(sankey_matrix$weight2)

  # object@reductions$LR_TF <- l_r_tf
  object@reductions$sankey <- sankey_matrix
  return(object)
}




