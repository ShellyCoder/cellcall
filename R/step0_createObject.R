#' The CellInter Class
#' The CellInter object is the center of each cell-cell communication analysis with the scRNA-seq.
#' @slot data  The normalized expression matrix, mean matrix, gsea result, lr score matrix, and lr score matrix(log-scale)
#' @slot meta.data cell type, barcode, nFeature and count of each cell
#' @slot reductions data.frame of L-R-TF with specific cellA-cellB, later for sankey plot
#' @slot project the name of this analysis
#' @slot Org the species of this scRNA-seq
#' @slot version the version of this R package
#' @import methods
#' @name CellInter-class
#' @aliases CellInter-class
#' @exportClass CellInter
#' @export
#'
CellInter <- setClass("CellInter", slots = list(data = "list",
                                                meta.data = "data.frame",
                                                reductions = "list",
                                                project = "character",
                                                Org = "character",
                                                version = "character"),

                      prototype = list(data = list(),
                                       meta.data = data.frame(),
                                       reductions = list(),
                                       project = "Microenvironment",
                                       Org = "Hsapiens",
                                       version = "1.0")
)


#' create a Cellwave objects
#' @param data a dataframe with row of gene and column of sample
#' @param min.feature Include cells where at least this many features are detected
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE_CELLTYPE, set this to "_" to separate the cell name into its component parts for picking the relevant field.
#' @param names.field BARCODE_CELLTYPE in the input matrix, set names.field to 2 to set the initial identities to CELLTYPE.
#' @param project Project name for the this object
#' @param source the type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @return the value of \code{cellwave object}
#' @import graphics
#' @import methods
#' @importFrom stringr str_split
#' @export

# 定义CreateNichConObject的现实，并指定参数类型为nichcon对象
CreateNichConObject <- function(data,
                   min.feature = 3,
                   names.field = 1,
                   names.delim = "_",
                   project = "Microenvironment",
                   source = "UMI", # "UMI" or "fullLength"
                   scale.factor = 10^6,
                   Org = "Homo sapiens" # "Homo sapiens" or "Mus musculus"
                   )
{

  nichcon <- new("CellInter",project = "Microenvironment",Org = Org)

  if(!is.data.frame(data)){
    stop("data must be a dataframe.")
  }

  sampleID <- colnames(data)
  if(sum(duplicated(sampleID))>0){
    stop("The cell ID should be unique.")
  }

  cell_type <- stringr::str_split(colnames(data), names.delim, simplify = T)[,names.field]
  if(length(grep("-",unique(cell_type)))>0){
    stop("The cell ID can't contain '-'.")
  }
  if(length(grep("_",unique(cell_type)))>0){
    stop("The cell ID can't contain '_'.")
  }

  nFeature <- apply(data, 2, function(x) {sum(x>0)})
  nCounts <- colSums(data)

  init.meta.data <- data.frame(sampleID = sampleID, celltype = cell_type, nFeature = nFeature, nCounts = nCounts)

  # source("./myProject/code/project/counts2rpm.R")
  if(source=="UMI"){
    data_cpm <- counts2normalized_10X(data, toType = "CPM", scale.factor=scale.factor)
    # data_cpm_log <- log2(data_cpm+1)
  }else if(source=="fullLength"){
    data_cpm <- counts2normalized_smartseq2(data, Org, "TPM", scale.factor=scale.factor)
    # data_cpm_log <- log2(data_cpm+1)
  }else if(source=="TPM"){
    data_cpm <- data
    # data_cpm_log <- log2(data+1)
  }else if(source=="CPM"){
    data_cpm <- data
    # data_cpm_log <- log2(data+1)
  }

  data.list = list()
  # data.list = c(data.list, list(count = data, data = data_cpm_log, withoutlog = data_cpm))
  data.list = c(data.list, list(count = data, withoutlog = data_cpm))

  nichcon@meta.data <- init.meta.data
  nichcon@data <- data.list

  return(nichcon)
}

