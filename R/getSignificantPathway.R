#' enrich communication relation on the pathway
#' @param data a dataframe of communication score with row LR and column cellA-cellB
#' @param cella_cellb explore the LR between sender cellA and receiver cellB, eg: "A-B"
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @param Org choose the species source of gene, only "Homo sapiens" in this version.
#' @return the dataframe with column: Pvalue, Jaccard, NES and pathway
#' @importFrom utils read.table
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom stats phyper sd
#' @export

getHyperPathway <- function(data, object, cella_cellb, IS_core=TRUE, Org="Homo sapiens"){
  n<- data
  Org <- Org
  mt <- object

  if(Org == 'Homo sapiens'){
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
  }else if(Org == 'Mus musculus'){
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

  tmp <- triple_relation
  tmp$triple <- paste(triple_relation$Ligand_Symbol, triple_relation$Receptor_Symbol, triple_relation$TF_Symbol, sep = "-")
  list.tmp <- apply(tmp, 1, function(x){
    str_split(x[5],",",simplify = F) %>% unlist()  -> path.tmp
    tmp.df.part <- data.frame(path=path.tmp, triple=rep(as.character(x[9]), length(path.tmp)))
  })
  path.triple.tmp.df <- do.call(rbind, list.tmp)

  path.list.tmp <- lapply(rownames(n), function(x){
    names.tmp <- colnames(n)[n[x,]>0]

    if(length(names.tmp)>0){
      # print(x)
      l.tmp <- str_split(x,"-")[[1]][1]
      r.tmp <- str_split(x,"-")[[1]][2]
      a.tmp <- triple_relation %>% dplyr::filter(Ligand_Symbol==l.tmp & Receptor_Symbol==r.tmp) %>% .[,c('pathway_ID','TF_Symbol')]
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
  cc.tmp.triple <- unique(dplyr::filter(path.list.tmp.df, cc==cc.tmp)[,1:2,drop=T] %>% as.data.frame())
  if(nrow(cc.tmp.triple)==0){
    res.df <- data.frame(Pvalue=1, Jaccard=0, NES=0, pathway="hsa04330")
    rownames(res.df) <- "hsa04330"
    return(res.df)
  }

  myHypergeometric <- function(intersect_mirna=intersect_mirna, mrna_has_mirna=mrna_has_mirna, all_mirna=all_mirna, lncrna_has_mirna=lncrna_has_mirna){
    q = intersect_mirna-1 
    m = mrna_has_mirna  
    n = all_mirna-mrna_has_mirna 
    k = lncrna_has_mirna 
    stats::phyper(q=q, m=m, n=n, k=k, log = FALSE, lower.tail = FALSE)
  }

  cc.tmp.triple$triple <- as.character(cc.tmp.triple$triple)
  cc.tmp.triple$path <- as.character(cc.tmp.triple$path)

  cc.tmp.pathway <- unique(cc.tmp.triple$path)
  do.call(rbind, lapply(cc.tmp.pathway, function(x){
    unique(dplyr::filter(path.triple.tmp.df, path==x)[,2]) -> pathway.triple.tmp
    unique(dplyr::filter(cc.tmp.triple, path==x)[,1]) -> cc.triple.tmp
    intersect_mirna <- length(intersect(pathway.triple.tmp, cc.triple.tmp))
    mrna_has_mirna <- length(pathway.triple.tmp)
    all_mirna <- length(unique(path.triple.tmp.df$triple))
    lncrna_has_mirna <- length(unique(cc.tmp.triple$triple))

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

  all.path.tmp <- unique(path.triple.tmp.df$path)
  all.path.Jaccard.tmp <- do.call(rbind, lapply(all.path.tmp, function(x){
    unique(dplyr::filter(path.triple.tmp.df, path==x)[,2]) -> pathway.triple.tmp
    unique(dplyr::filter(cc.tmp.triple, path==x)[,1]) -> cc.triple.tmp

    Jaccard.tmp <- length(intersect(pathway.triple.tmp, cc.triple.tmp))/length(unique(c(pathway.triple.tmp, cc.triple.tmp)))
    return(Jaccard.tmp)
  })) %>% unlist()

  mean.tmp <- base::mean(all.path.Jaccard.tmp) 
  sd.tmp <- stats::sd(all.path.Jaccard.tmp)
  res.df$NES <- (res.df$Jaccard-mean.tmp)/sd.tmp

  res.df$pathway <- rownames(res.df)

  return(res.df)
}
