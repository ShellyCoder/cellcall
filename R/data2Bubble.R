#' melt enrich dataframe for bubble graph
#' @param pathway.hyper.list a list of pathway enrichment pvalue and NES
#' @return the result for pbubble graph
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join inner_join
#' @importFrom utils read.table
#' @export

getForBubble <- function(pathway.hyper.list = pathway.hyper.list){
  ## 合并不同cc间的p值，jaccard系数和NES值
  # library(reshape2)
  pathway.hyper.df <- Reduce(function(x,y) full_join(x, y, by=c("pathway")), pathway.hyper.list, accumulate =FALSE)
  rownames(pathway.hyper.df) <- pathway.hyper.df$pathway
  pathway.hyper.df$pathway <- NULL

  ## 拆分后，获得不同cc间的p值
  pathway.hyper.df.pvalue <- pathway.hyper.df[, seq(from=1, to=ncol(pathway.hyper.df), by=3)]
  colnames(pathway.hyper.df.pvalue) <- colnames(n)
  rownames(pathway.hyper.df.pvalue) <- rownames(pathway.hyper.df)
  pathway.hyper.df.pvalue[is.na(pathway.hyper.df.pvalue)] <- 1

  ## 拆分后，获得不同cc间的NES值
  pathway.hyper.df.NES <- pathway.hyper.df[, seq(from=3, to=ncol(pathway.hyper.df), by=3)]
  colnames(pathway.hyper.df.NES) <- colnames(n)
  rownames(pathway.hyper.df.NES) <- rownames(pathway.hyper.df)
  pathway.hyper.df.NES[is.na(pathway.hyper.df.NES)] <- 0

  ## 删除不涉及的pathway
  pathway.hyper.df.NES <- pathway.hyper.df.NES[rowSums(pathway.hyper.df.NES)!=0,]
  pathway.hyper.df.pvalue <- pathway.hyper.df.pvalue[rowSums(pathway.hyper.df.NES)!=0,]

  f.tmp <- system.file("extdata", "KEGG_SYMBOL_ID.txt", package="cellwave")
  pathway.info <- read.table(f.tmp, sep = '\t', quote = "", header = FALSE, stringsAsFactors = F)
  colnames(pathway.info) <- c('id', 'name', "main.object.name", "sub.object.name")
  rownames(pathway.info) <- pathway.info$id

  ## melt 不同cc间的p值 适应 ggplot 绘图格式
  pathway.hyper.df.pvalue$pathway <- rownames(pathway.hyper.df.pvalue)
  pathway.hyper.df.pvalue$pathway <- pathway.info[pathway.hyper.df.pvalue$pathway, 'name']
  pvalue.melted <- melt(pathway.hyper.df.pvalue, id=c('pathway'),variable.name="cc",value.name="pvalue")

  ## melt 不同cc间的NES值 适应 ggplot 绘图格式
  pathway.hyper.df.NES$pathway <- rownames(pathway.hyper.df.NES)
  pathway.hyper.df.NES$pathway <- pathway.info[pathway.hyper.df.NES$pathway, 'name']
  NES.melted <- melt(pathway.hyper.df.NES, id=c('pathway'),variable.name="cc",value.name="NES")
  NES.melted$NES[NES.melted$NES > 2]=2
  NES.melted$NES[NES.melted$NES < 0]=0

  ## 合并上两步结果
  myPub.df <- inner_join(pvalue.melted, NES.melted, by=c('pathway', 'cc'))

  return(myPub.df)
}
