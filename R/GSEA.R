#' calculate foldchange for each celltype
#' @param inData a dataframe of gene expression
#' @param cell.type the cell type which you want to calculate the foldchange
#' @param method default is "median",eg "median" or "mean"
#' @param probs the quantile for median calculate
#' @return the foldchange value of each celltype in \code{cellwave object}
#' @importFrom stats median

mylog2foldChange.diy<-function(inData, cell.type, method="median", probs = 0.75) # median or mean
{
  re.list <- list()
  cell.type <- cell.type
  print(cell.type)
  for (c in cell.type) {
    print(c)
    cell_fc_labels<-colnames(inData)
    cell_fc_labels[cell_fc_labels==c]<-0
    cell_fc_labels[cell_fc_labels!='0']<-1

    classLabel <- cell_fc_labels
    sampleIdsCase<-which(classLabel==0);#0 tumer
    sampleIdsControl<-which(classLabel==1);#1 normal
    probeFC <- as.numeric(apply(inData, 1, function(x){
      # print(x)
      if(method=="mean"){
        (mean(as.numeric(x[sampleIdsCase]))+1)/(mean(as.numeric(x[sampleIdsControl]))+1);
      }else if(method=="median"){
        quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)

        a.tmp <- quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)+1
        b.tmp <- quantile(as.numeric(x[sampleIdsControl]), probs = probs, names=FALSE)+1
        a.tmp/b.tmp
      }
    }))

    probeFC<-log(probeFC,base=2);
    fc_res<-probeFC;

    fc_res[is.infinite(fc_res)]<-0
    fc_res[is.na(fc_res)]<-0
    res = data.frame(gene_id = as.vector(rownames(inData)), log2fc = fc_res)
    # print(dim(res))
    re.list <- c(re.list, list(res))
  }
  names(re.list) <- cell.type
  return(re.list)
}

#' calculate enrich for each celltype with all regulons
#' @param term_gene_list all regulons of genes denotes geneSet
#' @param FC_OF_CELL the cell type which you want to calculate the foldchange
#' @param minGSSize set the min size of each geneSet
#' @param maxGSSize set the max size of each geneSet
#' @return the enrich result of each celltype with all regulons
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr arrange desc

getGSEA<-function(term_gene_list, FC_OF_CELL, minGSSize=5, maxGSSize=500) {
  term_gene <- term_gene_list
  # library(dplyr)
  FC_OF_CELL <- FC_OF_CELL
  geneList.sort <- dplyr::arrange(FC_OF_CELL, desc(FC_OF_CELL$log2fc))

  ### geneList 是全表达谱排序后的foldchange  data(geneList)
  ##  TERM2GENE 一条感兴趣的term中所有基因
  ##
  ##  TERM2GENE is a data.frame with first column of term ID and second column of corresponding mapped gene
  ##  and TERM2NAME is a data.frame with first column of term ID and second column of corresponding term name.
  ##  TERM2NAME is optional

  rownames(geneList.sort) <- geneList.sort[,1]
  gene.name <- geneList.sort[,1]
  geneList.sort$log2fc <- as.numeric(geneList.sort$log2fc)
  geneList<-geneList.sort[,-1]
  names(geneList) <- gene.name

  # if most fc are not high enough will to error like: GSEA Error in if (abs(max.ES) > abs(min.ES))
  # you can set exponent = 0
  y<-GSEA(geneList = geneList, TERM2GENE = term_gene, minGSSize = minGSSize,exponent = 1,
          maxGSSize = maxGSSize, pvalueCutoff = 1, pAdjustMethod = "BH",
          verbose = FALSE)
  return(y)
}







