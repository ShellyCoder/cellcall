#' transform count to RPKM or TPM within Hsapiens and Mmusculus
#' @param data a dataframe with row of gene and column of sample
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @param toType whether RPKM or TPM of result you wanted, eg "RPKM", "TPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @return the RPKM or TPM of \code{data}
counts2normalized_smartseq2 <- function(data, Org, toType, scale.factor=10^6){
  if(Org=="Homo sapiens"){
    f.tmp <- system.file("extdata", "homo/transcriptmaxLength.Rdata", package="cellcall")
    load(f.tmp)
  }else if(Org=="Mus musculus"){
    f.tmp <- system.file("extdata", "mmus/transcriptmaxLength.Rdata", package="cellcall")
    load(f.tmp)
  }
  # head(g_l)
  a <- data
  a[1:4,1:4]
  ng=intersect(rownames(a),g_l$symbol) 

  exprSet=a[ng,]
  lengths=g_l[match(ng,g_l$symbol),2]
  head(lengths)
  head(rownames(exprSet))

  exprSet[1:4,1:4]
  if(toType=="RPKM"){
    total_count<- colSums(exprSet)
    head(total_count)
    head(lengths)

    rpkm <- t(do.call( rbind,
                       lapply(1:length(total_count),
                              function(i){
                                10^9*exprSet[,i]/lengths/total_count[i]
                              }) ))
    return(rpkm)
  }else if(toType=="TPM"){
    rpk <- exprSet/lengths
    total_count<- colSums(rpk)
    tpm <- t(do.call( rbind,
                      lapply(1:length(total_count),
                             function(i){
                               scale.factor*rpk[,i]/total_count[i]
                             }) ))
    tpm <- data.frame(tpm, stringsAsFactors = FALSE)
    colnames(tpm) <- colnames(exprSet)
    rownames(tpm) <- rownames(exprSet)

    return(tpm)
  }

}

#' transform count to CPM
#' @param data a dataframe with row of gene and column of sample
#' @param toType whether RPKM or TPM of result you wanted, eg "CPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @return the CPM of \code{data}
counts2normalized_10X <- function(data, toType, scale.factor=10^6){
  exprSet <- data
  if(toType=="CPM"){
    total_count<- colSums(exprSet)

    cpm <- t(do.call( rbind,
                       lapply(1:length(total_count),
                              function(i){
                                scale.factor*exprSet[,i]/total_count[i]
                              }) ))
    cpm <- data.frame(cpm, stringsAsFactors = FALSE)
    colnames(cpm) <- colnames(exprSet)
    rownames(cpm) <- rownames(exprSet)
    return(cpm)
  }
}
