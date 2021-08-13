#' get CommuProfile from a Cellwave objects
#' @param object a Cellwave objects
#' @param probs Percentile of gene expression in one cell type to represents this cell type
#' @param use.type the type of compute, default is "median"
#' @param pValueCor firlter target gene of TF with spearson, p > pValueCor, default is 0.05
#' @param CorValue firlter target gene of TF with spearson, value > CorValue, default is 0.1
#' @param topTargetCor use topTargetCor of candidate genes which has firlter by above parameters, default is 1, means 100%
#' @param p.adjust gsea pValue of regulons with BH adjusted threshold, default is 0.05
#' @param method "weighted", "max", "mean", of which "weighted" is default. choose the proper method to score downstream activation of ligand-receptor all regulons of given ligand-receptor relation
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @return the value of \code{cellwave object}
#' @import graphics
#' @import methods
#' @export

TransCommuProfile <- function(object,
                              pValueCor = 0.05,
                              CorValue = 0.1,
                              topTargetCor=1,
                              p.adjust=0.05,
                              use.type="median",
                              probs = 0.75,
                              method="weighted",  
                              Org = 'Homo sapiens', 
                              IS_core = TRUE
                              )
{
    is_myObject <- is(object, "CellInter")
    if(is_myObject){
      profile <- ConnectProfile(object,
                     pValueCor = pValueCor,
                     CorValue = CorValue,
                     topTargetCor = topTargetCor,
                     p.adjust=p.adjust,
                     use.type=use.type,
                     probs = probs,
                     method = method,
                     IS_core = IS_core,
                     Org = Org)

      object@data = c(object@data, profile)
      return(object)
    }else{
      stop("object should be CellInter type.")
    }
}
