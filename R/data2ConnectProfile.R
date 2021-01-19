#' create a Cellwave objects
#' @param object a Cellwave objects
#' @param probs Percentile of gene expression in one cell type to represents this cell type
#' @param use.type the type of compute, default is "median"
#' @param pValueCor firlter target gene of TF with spearson, p > pValueCor, default is 0.05
#' @param CorValue firlter target gene of TF with spearson, value > CorValue, default is 0.1
#' @param topTargetCor use topTargetCor of candidate genes which has firlter by above parameters, default is 1, means 100%
#' @param p.adjust gsea pValue of regulons with BH adjusted threshold, default is 0.05
#' @param method "weighted", "max", "mean", of which "weighted" is default. choose the proper method to score downstream activation of ligand-receptor all regulons of given ligand-receptor relation
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @return the result dataframe of \code{cell communication}
#' @importFrom stringr str_split
#' @importFrom stats quantile median
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom utils read.table head
#' @export

ConnectProfile <- function(object, pValueCor=0.05, CorValue=0.1, topTargetCor=1, method="weighted", p.adjust=0.05, use.type="median", probs = 0.75, Org = 'Homo sapiens'){
  Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv, 修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
  options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor

  if(Org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "ligand_receptor_TFs.txt", package="cellwave")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    f.tmp <- system.file("extdata", "tf_target.txt", package="cellwave")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(Org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "ligand_receptor_TFs_homology.txt", package="cellwave")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellwave")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }

  triple_relation$pathway_ID <- NULL
  print(triple_relation[1:4,])
  complex_tmp <- triple_relation$Receptor_Symbol[grep(",",triple_relation$Receptor_Symbol)] %>% unique()
  tmp_complex_symbol <- triple_relation$Receptor_Symbol[grep(",",triple_relation$Receptor_Symbol)] %>% unique() %>% str_split(",") %>% unlist %>% unique()
  all.gene.needed <- unique(as.character(c(triple_relation$Ligand_Symbol, triple_relation$Receptor_Symbol, triple_relation$TF_Symbol, target_relation$TF_Symbol, target_relation$Target_Symbol,tmp_complex_symbol)))
  # triple_relation[1:4,1:4]
  # target_relation[1:4,1:4]
  my_Expr <- object@data$withoutlog
  colnames(my_Expr) <- as.character(object@meta.data$celltype)
  my_Expr[1:4,1:4]
  detect_gene <- rownames(my_Expr)

  expr_set <- my_Expr[intersect(detect_gene, all.gene.needed),]  ## 表达矩阵输入
  detect_gene <- rownames(expr_set)
  cell_type = unique(colnames(expr_set))
  expr.fc <- object@data$withoutlog[detect_gene,]
  colnames(expr.fc) <- colnames(expr_set)

  ## 读入complex_list ,在res 最后几行添加
  complex_matrix <- matrix(ncol = length(colnames(expr_set)))
  complex_matrix <- as.data.frame(complex_matrix)
  colnames(complex_matrix) <- colnames(expr_set)
  myrownames <- c()

  complex <- complex_tmp
  for(i in 1:length(complex)){
    i_tmp = strsplit(complex[i], ',')
    # print(i_tmp)
    if( sum(i_tmp[[1]] %in% detect_gene) == length(i_tmp[[1]]) ){
      tmp_df <- expr_set[i_tmp[[1]],]
      tmp_mean <- colMeans(tmp_df)
      tmp_index <- unique(unlist(apply(tmp_df, 1,function(x) {which(x==0)}))) # complex任意一个表达为0, 该complex不存在
      tmp_mean[tmp_index] <- 0

      # print(res_tmp)
      complex_matrix <- rbind(complex_matrix, tmp_mean)
      myrownames <- c(myrownames, complex[i])
    }
  }

  complex_matrix <- complex_matrix[-1,]

  ## 把complex的联合表达值加上
  if(nrow(complex_matrix) > 0){
    rownames(complex_matrix) <- myrownames
    expr_set <- rbind(expr_set, complex_matrix)
  }

  expr_set <- expr_set[apply(expr_set, 1, function(x){sum(x!=0)})>0,] ##删除没检测到的gene
  detect_gene <- rownames(expr_set)
  # expr_set[1:4,1:4]

  ### 计算每个细胞类型的 gene 中位数
  print("step1: compute means of gene")
  expr_mean <- matrix(nrow = nrow(expr_set), ncol = length(cell_type))
  myColnames <- c()
  for (i in 1:length(cell_type)) {
    myCell <- cell_type[i]
    myMatrix <- expr_set[,colnames(expr_set)==myCell,drop=F]
    if(use.type=="mean"){
      myMatrix_mean <- as.numeric(apply(myMatrix, 1, mean))
    }else if(use.type=="median"){
      quantil.tmp <- as.numeric(apply(myMatrix, 1, function(x){
          quantile(x, probs = probs,names=FALSE)
      }))
      mean.tmp <- rowMeans(myMatrix)
      mean.tmp[which(quantil.tmp==0)]<-0 # 3/4 以上全为0的gene均值也为 0
      myMatrix_mean <- mean.tmp## 每个细胞类型的gene求均值
    }
    expr_mean[,i] <- myMatrix_mean
    myColnames <- c(myColnames, myCell)
    # print(myCell)
  }
  expr_mean <- data.frame(expr_mean)
  colnames(expr_mean) <- myColnames
  rownames(expr_mean) <- rownames(expr_set)

  expr_mean <- expr_mean[apply(expr_mean, 1, function(x){sum(x!=0)})>0,]
  detect_gene <- rownames(expr_mean)

  # gene rank for each cell type for following GSEA
  # source("./myProject/code/project/GSEA.R")
  if(use.type=="median"){
    # fc.list <- mylog2foldChange(inData = expr.fc, cell.type = cell_type, method="mean", probs = probs)
    fc.list <- mylog2foldChange.diy(inData = expr.fc, cell.type = cell_type, method="median", probs = probs)
  }

  #  tf -> spearson cor ( first filter) -> gsea enrich testing targets trend of expression
  print("step2: filrter tf-gene with correlation, then score regulons")
  # source("./myProject/code/project/getCorrelatedGene.R")
  tfs_set <- unique(triple_relation$TF_Symbol)
  regulons_matrix <- matrix(data = 0, nrow = length(tfs_set), ncol = length(cell_type))

  regulons_matrix <- as.data.frame(regulons_matrix)
  rownames(regulons_matrix) <- tfs_set
  colnames(regulons_matrix) <- cell_type
  my_minGSSize <- 5

  for (i in cell_type) {
    print(i)
    tf_val <- lapply(tfs_set, function(x) {
      if(x %in% detect_gene){
        # Chip-seq实验数据 -> TF-Target
        # print(x)
        targets <- target_relation[which(target_relation$TF_Symbol==x),2]
        targets <- targets[targets %in% detect_gene]
        if(length(targets)<=0){
          return(0)
        }
        # 表达值相关(spearman) -> TF-Target [ source("./project/getCorrelatedGene.R") ]
        corGene_tmp <- getCorrelatedGene(data = expr_set,cell_type = i,tf=x, target_list=targets, pValue=pValueCor, corValue=CorValue,topGene=topTargetCor)
        common_targets_tmp <- intersect(corGene_tmp, targets)
        # print(length(common_targets_tmp))
        if(length(common_targets_tmp)==0){
          return(0)
        }
        # print(length(common_targets_tmp))
        gene.name.tmp <- common_targets_tmp
        term_gene_list.tmp <- data.frame(term.name=rep(1, length(gene.name.tmp)), gene=gene.name.tmp)

        if(length(gene.name.tmp)<my_minGSSize){
          return(0)
        }
        # print(gene.name.tmp)
        tryCatch({
          nes.tmp <- getGSEA(term_gene_list = term_gene_list.tmp,
                             FC_OF_CELL = fc.list[[i]], minGSSize=my_minGSSize, maxGSSize=500)
          if(length(nes.tmp@result$NES)>0 & length(nes.tmp@result$p.adjust)>0 & expr_mean[x,i]>0){
            if(nes.tmp@result$p.adjust<p.adjust & nes.tmp@result$NES>0){
              tf.val.enriched <- nes.tmp@result$NES
              return(tf.val.enriched)
            }else{
              return(0)
            }
          }else{
            return(0)
          }
        },error=function(e){
          return(0)
        })

      }else{
        return(0)
      }
    })
    tf_val <- unlist(tf_val)
    regulons_matrix[,i] <- tf_val
    print(sum(tf_val>0))
  }

  gsea.list <- list()
  gsea.genes.list <- list()
  for (i in cell_type) {
    print(i)
    print(length(which(regulons_matrix[,i]!=0)))
    if(length(which(regulons_matrix[,i]!=0))!=0){
      tfs_set.tmp <- tfs_set[which(regulons_matrix[,i]!=0)]
      tf_val.df <- do.call(rbind, lapply(tfs_set.tmp, function(x) {
        # Chip-seq实验数据 -> TF-Target
        targets <- target_relation[which(target_relation$TF_Symbol==x),2]
        targets <- targets[targets %in% detect_gene]
        if(length(targets)<=0){
          return(NULL)
        }
        # 表达值相关(spearman) -> TF-Target [ source("./project/getCorrelatedGene.R") ]
        corGene_tmp <- getCorrelatedGene(data = expr_set,cell_type = i,tf=x, target_list=targets, pValue=pValueCor, corValue=CorValue,topGene=topTargetCor)
        common_targets_tmp <- intersect(corGene_tmp, targets)
        # print(length(common_targets_tmp))
        if(length(common_targets_tmp)==0){
          return(NULL)
        }
        # print(length(common_targets_tmp))
        gene.name.tmp <- common_targets_tmp
        if(length(gene.name.tmp)<my_minGSSize){
          return(NULL)
        }

        term_gene_list.tmp <- data.frame(term.name=rep(x, length(gene.name.tmp)), gene=gene.name.tmp)
        return(term_gene_list.tmp)
      }))
      nes.tmp <- getGSEA(term_gene_list = tf_val.df, FC_OF_CELL = fc.list[[i]], minGSSize=my_minGSSize, maxGSSize=500)

      gsea.list <- c(gsea.list, list(nes.tmp))
    }else{
      gsea.list <- c(gsea.list, list(vector()))
    }
  }
  names(gsea.list) <- cell_type


  print("step3: get distance between receptor and tf in pathway")
  # source("./myProject/code/project/getDistanceKEGG.R")
  DistanceKEGG <- getDistanceKEGG(data = triple_relation,method = "mean")

  print("step4: score downstream activation of ligand-receptor all regulons of given ligand-receptor relation (weighted, max, or mean) ####")
  l_r_inter <- unique(triple_relation[,5:6])
  expr_r_regulons <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type)) ## 细胞通讯的大类信息（细胞类型） ：A->A,A->B,A->C,,,,C->C
  expr_r_regulons <- as.data.frame(expr_r_regulons)
  rownames(expr_r_regulons) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames(expr_r_regulons) <- cell_type

  for (n in 1:nrow(l_r_inter)) {
    sender_tmp <- l_r_inter[n,1]
    receiver_tmp <- l_r_inter[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    # print(n)

    val_tmp = 0
    if( sum(l_r_inter[n,] %in% detect_gene)==2 ){

      info_tmp <- dplyr::filter(triple_relation, Ligand_Symbol==sender_tmp & Receptor_Symbol==receiver_tmp)[,5:7]
      tfs_tmp <- info_tmp$TF_Symbol[info_tmp$TF_Symbol %in% detect_gene]
      if(length(tfs_tmp) > 0){
        regulon_tmp_df <- regulons_matrix[tfs_tmp,]
        if(method=='max'){
          expr_r_regulons[row_index,] = as.numeric(apply(regulon_tmp_df, 2, function(x){max(x)}))
        }else if(method=="weighted"){
          distance2w_tmp <- (1/DistanceKEGG[row_index,tfs_tmp])
          w_tmp<- distance2w_tmp/sum(distance2w_tmp)
          expr_r_regulons[row_index,] = as.numeric(apply(regulon_tmp_df, 2, function(x){
            sum(w_tmp*x)
          }))
        }else if(method=="mean"){
          expr_r_regulons[row_index,] = as.numeric(apply(regulon_tmp_df, 2, function(x){mean(x)}))
        }
      }
    }
  }

  print("step5: softmax for ligand")
  # softmax for ligand
  ligand_symbol <- unique(triple_relation$Ligand_Symbol)
  softmax_ligand <- expr_mean[intersect(ligand_symbol, detect_gene),]
  colnames(softmax_ligand) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_ligand)

  softmax_ligand <- do.call(rbind,lapply(1:nrow(softmax_ligand), function(i){
    softmax_ligand[i,]/rowCounts[i]
  }))

  # softmax for receptor
  receptor_symbol <- unique(triple_relation$Receptor_Symbol)
  softmax_receptor <- expr_mean[intersect(receptor_symbol, detect_gene),]
  colnames(softmax_receptor) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_receptor)

  softmax_receptor <- do.call(rbind,lapply(1:nrow(softmax_receptor), function(i){
    softmax_receptor[i,]/rowCounts[i]
  }))

  #  l-r in cell type level
  print("step6: score ligand-receptor relation (weighted, max, or mean) ####")

  l_r_inter <- unique(triple_relation[,5:6])
  expr_l_r <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type)^2) ## 细胞通讯的大类信息（细胞类型） ：A->A,A->B,A->C,,,,C->C
  expr_l_r <- as.data.frame(expr_l_r)
  rownames(expr_l_r) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  myColnames <- character()
  for (i in cell_type) {
    for (j in cell_type) {
      myColnames <- c(myColnames, paste(i,j,sep = "-"))
    }
  }
  colnames(expr_l_r) <- myColnames

  for (n in 1:nrow(l_r_inter)) {
    sender_tmp <- l_r_inter[n,1]
    receiver_tmp <- l_r_inter[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    # print(n)
    for (i in cell_type) {
      for (j in cell_type) {
        myColnames <- c(myColnames, paste(i,j,sep = "-"))
        val_tmp = 0
        if( sum(l_r_inter[n,] %in% detect_gene)==2 ){
          sender_val <- expr_mean[sender_tmp,i]
          receiver_val <- expr_mean[receiver_tmp,j]
          tf_val <- expr_r_regulons[row_index,j]

          if(tf_val > 0 & sender_val>0 & receiver_val >0){
            # val_tmp <- sender_val^2 + receiver_val^2 + tf_val  #sender_val^2 + receiver_val^2 + tf_val^2
            sender_val_weighted <- softmax_ligand[sender_tmp, i]
            receiver_val_weighted <- softmax_receptor[receiver_tmp, j]
            # val_tmp <- sender_val_weighted * (receiver_val + tf_val*(sender_val+receiver_val)/(sender_val+receiver_val+tf_val))
            val_tmp <- 100*(sender_val_weighted^2 + receiver_val_weighted^2) * tf_val
            # print(val_tmp)
          }else{
            val_tmp = 0
          }
        }else{
          val_tmp = 0
        }
        col_index_tmp <- paste(i,j,sep = "-")
        expr_l_r[n,col_index_tmp] <- val_tmp
      }
    }
  }

  expr_l_r <- expr_l_r[apply(expr_l_r, 1, function(x){sum(x!=0)})>0,] ##删除没检测到的L-R
  # expr_l_r <- expr_l_r[,which(colSums(expr_l_r)!=0)] ##删除没检测到的cell-M_cell-N
  expr_l_r <- as.data.frame(expr_l_r)

  # test = as.data.frame(exp(tmp2))
  expr_l_r_log2 <- log2(expr_l_r+1)
  expr_l_r_log2_scale <- (expr_l_r_log2-min(expr_l_r_log2))/(max(expr_l_r_log2)-min(expr_l_r_log2))


  Result <- list(expr_mean = expr_mean,
                 regulons_matrix = regulons_matrix,
                 gsea.list = gsea.list,
                 fc.list = fc.list,
                 expr_r_regulons = expr_r_regulons,
                 softmax_ligand = softmax_ligand,
                 softmax_receptor = softmax_receptor,
                 expr_l_r =  expr_l_r,
                 expr_l_r_log2 = expr_l_r_log2,
                 expr_l_r_log2_scale = expr_l_r_log2_scale,
                 DistanceKEGG= DistanceKEGG)

  return(Result)
}










