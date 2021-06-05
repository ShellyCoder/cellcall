#' Create Nichobject from Seurat object
#' @param Seurat.object  The Seurat object which stores the expression matrix and cell type information
#' @param slot The name of slot which contains expression matrix, default "counts".
#' @param cell_type The name of specific column which contains cell type information, default "orig.ident".
#' @param source the type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @param Org the species of this scRNA-seq
#' @importFrom Seurat GetAssayData
#' @export
#'
CreateObject_fromSeurat <- function(Seurat.object, slot="counts", 
                                    cell_type="orig.ident", data_source="UMI", scale.factor=10^6,
                                    Org="Homo sapiens"){
  Seurat.object <- Seurat.object
  slot <- slot # counts, data
  cell_type <- cell_type
  your_source <- data_source # "UMI", "fullLength", "TPM", or "CPM".
  
  is_Seurat <- is(Seurat.object, "Seurat")
  
  if(is_Seurat){
    if(require(Seurat)){
      
      tryCatch({
        data <-  as.matrix(Seurat::GetAssayData(object = Seurat.object, slot = slot, assay = "RNA"))
      },error=function(e){
        stop("there is no slot: ", slot," in your seurat object in RNA assay")
      })
      
      tryCatch({
        myCelltype <- as.character(Seurat.object@meta.data[,cell_type])
      },error=function(e){
        stop("there is no columns: ", cell_type," in meta.data of Seurat.")
      })
      
      data <- as.data.frame(data)
      colnames(data) <- paste(1:length(myCelltype), myCelltype, sep = "_")
      
      mt <- CreateNichConObject(data=data, min.feature = 0,
                                names.field = 2,
                                names.delim = "_",
                                source = your_source, # fullLength, UMI, TPM
                                scale.factor = scale.factor,
                                Org = Org,
                                project = "Microenvironment")
      
      return(mt)
      
    }else{
      stop("No package Seurat! Fail to extract data from Seurat object.")
    }
  }else{
    stop("It's not Seurat object!")
  }
  
}








