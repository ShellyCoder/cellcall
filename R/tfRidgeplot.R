##' plot ridge graph
##' @param x gseaResult object
##' @param showCategory number of categories for plotting
##' @param fill one of "pvalue", "p.adjust", "qvalue"
##' @param core_enrichment whether only using core_enriched genes
##' @param orderBy The order of the Y-axis, default is NES, other colnames is ok, eg: "ID", "Description", "setSize", "enrichmentScore", "p.adjust".
##' @param decreasing logical. Should the orderBy order be increasing or decreasing?
##' @importFrom ggplot2 scale_fill_gradientn aes_string scale_fill_continuous xlab ylab guide_colorbar
##' @importFrom ggridges geom_density_ridges
##' @importFrom DOSE theme_dose geneInCategory
##' @export
ridgeplot.DIY <- function(x, showCategory = targetTFindex, fill="p.adjust",
                                 core_enrichment = TRUE,
                                 orderBy = "NES", decreasing = FALSE) {
  if (!is(x, "gseaResult"))
    stop("currently only support gseaResult")

  if (!fill %in% colnames(x@result)) {
    stop("'fill' variable not available ...")
  }

  if (orderBy !=  'NES' && !orderBy %in% colnames(x@result)) {
    message('wrong orderBy parameter; set to default `orderBy = "NES"`')
    orderBy <- "NES"
  }
  n <- showCategory
  if (core_enrichment) {
    gs2id <- geneInCategory(x)[n]
  } else {
    gs2id <- x@geneSets[x$ID[n]]
  }

  gs2val <- lapply(gs2id, function(id) {
    res <- x@geneList[id]
    res <- res[!is.na(res)]
  })

  nn <- names(gs2val)
  i <- match(nn, x$ID)
  nn <- x$Description[i]

  j <- order(x@result[[orderBy]][i], decreasing = decreasing)
  len <- sapply(gs2val, length)
  gs2val.df <- data.frame(category = rep(nn, times=len),
                          color = rep(x[i, fill], times=len),
                          value = unlist(gs2val))

  colnames(gs2val.df)[2] <- fill
  gs2val.df$category <- factor(gs2val.df$category, levels=nn[j])

  ggplot(gs2val.df, aes_string(x="value", y="category", fill=fill)) +
    ggridges::geom_density_ridges() +
    scale_fill_continuous(low="red", high="blue", name = fill,
                          guide=guide_colorbar(reverse=TRUE)) +
    xlab(NULL) + ylab(NULL) +  DOSE::theme_dose()
}
