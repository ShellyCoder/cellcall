#' enrich communication relation on the pathway
#' @param data a dataframe of communication score with row LR and column cellA-cellB
#' @param cella_cellb explore the LR between sender cellA and receiver cellB, eg: "A-B"
#' @param Org choose the species source of gene, only "Homo sapiens" in this version.
#' @return the dataframe with column: Pvalue, Jaccard, NES and pathway
#' @importFrom utils read.table
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom stats phyper sd
#' @export

getHyperPathway <- function(data, cella_cellb, Org="Homo sapiens"){
  n<- data
  Org <- Org
  # library(stringr)
  ## 统计每条通路中涉及到的三元关系
  if(Org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "ligand_receptor_TFs.txt", package="cellwave")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(Org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "ligand_receptor_TFs_homology.txt", package="cellwave")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }

  tmp <- triple_relation
  tmp$triple <- paste(triple_relation$Ligand_Symbol, triple_relation$Receptor_Symbol, triple_relation$TF_Symbol, sep = "-")
  list.tmp <- apply(tmp, 1, function(x){
    str_split(x[5],",",simplify = F) %>% unlist()  -> path.tmp
    tmp.df.part <- data.frame(path=path.tmp, triple=rep(as.character(x[9]), length(path.tmp)))
  })
  path.triple.tmp.df <- do.call(rbind, list.tmp)


  ## 统计每个细胞间存在的三元关系
  path.list.tmp <- lapply(rownames(n), function(x){
    names.tmp <- colnames(n)[n[x,]>0]

    if(length(names.tmp)>0){
      # print(x)
      l.tmp <- str_split(x,"-")[[1]][1]
      r.tmp <- str_split(x,"-")[[1]][2]
      a.tmp <- triple_relation %>% filter(Ligand_Symbol==l.tmp & Receptor_Symbol==r.tmp) %>% .[,c('pathway_ID','TF_Symbol')]
      recevier.tmp <- str_split(names.tmp, '-', simplify = T)[,2]
      df.tmp <- mt@data$regulons_matrix[a.tmp$TF_Symbol, recevier.tmp, drop=F]
      colnames(df.tmp) <- names.tmp
      df.tmp$path <- a.tmp$pathway_ID
      df.tmp$triple <- paste(x, a.tmp$TF_Symbol, sep = '-')
      df.longer.tmp <- df.tmp %>% tidyr::pivot_longer(cols = names.tmp, names_to = "cc", values_to = "count") %>% as.data.frame()
      df.longer.tmp <- df.longer.tmp[df.longer.tmp$count>0,]
      list.tmp <- apply(df.longer.tmp, 1, function(x){
        str_split(x[1],",",simplify = F) %>% unlist()  -> path.tmp
        tmp.df.part <- data.frame(triple=as.character(x[2]), path=path.tmp, cc=rep(as.character(x[3]), length(path.tmp)), count=rep(as.character(x[4]), length(path.tmp)))
      })
      list.tmp.df <- do.call(rbind, list.tmp)
      return(list.tmp.df)
    }
  })
  path.list.tmp.df <- do.call(rbind, path.list.tmp)

  apply(path.list.tmp.df, 1, function(x){
    row_index_tmp <- paste(str_split(x[1], '-')[[1]][1], str_split(x[1], '-')[[1]][2], sep = '-')
    col_index_tmp <- as.character(x[3])
    return(n[row_index_tmp, col_index_tmp])
  }) -> score
  path.list.tmp.df$lrscore <- score
  colnames(path.list.tmp.df)[4] <- "enrichment_score"

  cc.tmp <- cella_cellb
  cc.tmp.triple <- unique(filter(path.list.tmp.df, cc==cc.tmp)[,1:2,drop=T])
  if(nrow(cc.tmp.triple)==0){
    res.df <- data.frame(Pvalue=1, Jaccard=0, NES=0, pathway="hsa04330")
    rownames(res.df) <- "hsa04330"
    return(res.df)
  }

  ## 计算cc.tmp涉及到每个通路的超几何P值
  ## 预先定义好函数
  myHypergeometric <- function(intersect_mirna=intersect_mirna, mrna_has_mirna=mrna_has_mirna, all_mirna=all_mirna, lncrna_has_mirna=lncrna_has_mirna){
    q = intersect_mirna-1 # 抽到的白球个数-1，这里就是lncRNA和mrna共有miRNA个数-1，
    # 减1的原因这里算的是累计概率值，然后设置lower.tail = FALSE，最后得到的（1-累计概率值）
    m = mrna_has_mirna    # 白球的个数，这里就是具体某一个mrna拥有的miRNA个数
    n = all_mirna-mrna_has_mirna # 黑球的个数，这里就是所有的miRNA-具体某一个mrna拥有的miRNA个数,这里的所有miRNA其实就是背景miRNA的个数
    k = lncrna_has_mirna  # 从球袋里面无放回的拿球的个数，这里就是具体某一个lncRNA拥有的miRNA个数
    stats::phyper(q=q, m=m, n=n, k=k, log = FALSE, lower.tail = FALSE) # lower.tail = FALSE使得最后结果是（1-累计概率值）,也就是P[X > q],也就是P[X >= intersect_mirna]
  }

  cc.tmp.pathway <- unique(cc.tmp.triple$path)
  do.call(rbind, lapply(cc.tmp.pathway, function(x){
    unique(filter(path.triple.tmp.df, path==x)[,2]) -> pathway.triple.tmp
    unique(filter(cc.tmp.triple, path==x)[,1]) -> cc.triple.tmp
    intersect_mirna <- length(intersect(pathway.triple.tmp, cc.triple.tmp))
    mrna_has_mirna <- length(pathway.triple.tmp)
    all_mirna <- length(unique(path.triple.tmp.df$triple))
    lncrna_has_mirna <- length(unique(cc.tmp.triple$triple))

    ## p值和jaccard系数
    Hyper.pValue.tmp <- myHypergeometric(intersect_mirna = intersect_mirna,
                                         mrna_has_mirna = mrna_has_mirna,
                                         all_mirna = all_mirna,
                                         lncrna_has_mirna = lncrna_has_mirna)
    Jaccard.tmp <- length(intersect(pathway.triple.tmp, cc.triple.tmp))/length(unique(c(pathway.triple.tmp, cc.triple.tmp)))
    return(c(Hyper.pValue.tmp, Jaccard.tmp))
  })) -> res.df
  res.df <- as.data.frame(res.df)
  colnames(res.df) <- c("Pvalue", "Jaccard")
  rownames(res.df) <- cc.tmp.pathway

  ## all.path.Jaccard.tmp是特定cella_cellb设计到的三元关系在所有pathway的jaccard系数，作为后续矫正score的背景

  all.path.tmp <- unique(path.triple.tmp.df$path)
  all.path.Jaccard.tmp <- do.call(rbind, lapply(all.path.tmp, function(x){
    unique(filter(path.triple.tmp.df, path==x)[,2]) -> pathway.triple.tmp
    unique(filter(cc.tmp.triple, path==x)[,1]) -> cc.triple.tmp

    Jaccard.tmp <- length(intersect(pathway.triple.tmp, cc.triple.tmp))/length(unique(c(pathway.triple.tmp, cc.triple.tmp)))
    return(Jaccard.tmp)
  })) %>% unlist()

  ## 标准化jaccard系数
  mean.tmp <- base::mean(all.path.Jaccard.tmp) ## 均值
  sd.tmp <- stats::sd(all.path.Jaccard.tmp) ## 标准差
  res.df$NES <- (res.df$Jaccard-mean.tmp)/sd.tmp

  res.df$pathway <- rownames(res.df)

  return(res.df)
}
