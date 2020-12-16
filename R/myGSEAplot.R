#' plot enrichment graph
#' @param gsea.list a list of enrichment result from cellwave
#' @param myCelltype the cell type of receiver cell
#' @param fc.list foldchange list in the cellwave object
#' @param geneSetID the character of TF symbol, only significant activated can be inspected
#' @param selectedGeneID default is NULL, label the position of specific gene in FC flow.
#' @param mycol the color of each TF. the length is consistent with geneSetID
#' @importFrom enrichplot plot_grid
#' @importFrom ggplot2 ggplot geom_line scale_color_manual aes_ xlab ylab geom_hline theme_bw theme element_blank element_rect element_text margin geom_linerange scale_y_continuous geom_rect scale_fill_gradientn geom_segment geom_bar scale_fill_manual unit element_line
#' @importFrom ggrepel geom_text_repel
#' @export

getGSEAplot <- function(gsea.list, myCelltype, fc.list, geneSetID, selectedGeneID = NULL,
                        mycol = NULL)
  {
  Sys.setenv(LANGUAGE = "en") #显示英文报错信息
  options(stringsAsFactors = FALSE) #禁止chr转成factor
  fc.df <- fc.list[[myCelltype]]
  #按照foldchange排序
  fc.df.sorted <- fc.df[order(fc.df$log2fc, decreasing = T),]

  gsea.object <- gsea.list[[myCelltype]]

  # 突出显示感兴趣的基因 (在图中想标记出来的)
  # selectedGeneID <- c("ABCB1", "BAGE", "PRSS58")

  # 自定义足够多的颜色
  if(is.null(mycol)){
    mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9",
               "#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D",
               "#7CC767")
  }
  x <- gsea.object
  geneList <- position <- NULL ## to satisfy codetool

  #合并多条通路的数据
  gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
  gsdata$gsym <- rep(fc.df$gene_id, length(geneSetID))

  #########################################################################
  ##  running score - First Part
  #####################################################################################
  p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
    scale_color_manual(values = mycol) +

    #scale_x_continuous(expand=c(0,0)) + #两侧不留空
    geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + #在0的位置画虚线
    ylab("Enrichment\n Score") +

    theme_bw() + #背景为白色
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "top", legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent")) + # legend 位置

    theme(axis.text.y=element_text(size = 12, face = "bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b= 0, l=.2, unit="cm"))
  #p.res

  #########################################################################
  ##  rank lines - Second Part
  #####################################################################################
  # i的目的是给不同通路的线设定不同的y轴位置, 比如第一个通路的rank线0-1,第二个通路的rank线1-2,,,
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  #head(gsdata)

  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) +
    scale_color_manual(values = mycol) + #用自定义的颜色

    theme_bw() +
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    #scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  #p2
  #########################################################################
  ##  gradient rect - Third Part
  #####################################################################################
  test <- gsdata[1:floor(nrow(gsdata)/length(geneSetID)),]
  test$xend = test$x+1

  p.grad <- ggplot(test)+ geom_rect(aes(xmin = x,xmax = xend , ymin = 0 , ymax = 1, fill=geneList))+
    scale_fill_gradientn(colours = c("#035BFD", "#397EFC", "#5B94FB","white", "#F77A7C", "#F45557", "#FB0407"), limits = c(-max(abs((test$geneList))),max(abs((test$geneList)))))+
    theme_bw() + # 背景色
    theme(panel.grid = element_blank()) + #不画网格

    theme(legend.position = "none",
          plot.margin = margin(t=0, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    #scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))

  #########################################################################
  ##  run foldchange - Fourth Part
  #####################################################################################
  df2 <- p.res$data
  df2$y <- p.res$data$geneList[df2$x]
  df2$gsym <- p.res$data$gsym[df2$x]
  #head(df2)

  # 提取感兴趣的基因的变化倍数
  if(!is.null(selectedGeneID)){
    selectgenes <- data.frame(gsym = selectedGeneID)
    selectgenes <- merge(selectgenes, df2, by = "gsym")
    selectgenes <- unique(selectgenes[,c("x","y","gsym")])
    head(selectgenes)

    p.pos <- ggplot(selectgenes, aes(x, y, fill = "black", color = "black", label = gsym)) +
      geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                   color = "#80b1d3") +
      geom_bar(position = "dodge", stat = "identity",  width=0.5) +
      scale_fill_manual(values = "black", guide=FALSE) + #用自定义的颜色
      scale_color_manual(values = "black", guide=FALSE) + #用自定义的颜色

      #scale_x_continuous(expand=c(0,0)) +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +

      theme_bw() +
      theme(axis.text.y=element_text(size = 12, face = "bold"),
            panel.grid = element_blank()) +

      # 显示感兴趣的基因的基因名
      geom_text_repel(data = selectgenes,
                      show.legend = FALSE, #不显示图例
                      direction = "x", #基因名横向排列在x轴方向
                      ylim = c(2, NA), #基因名画在-2下方
                      angle = 90, #基因名竖着写
                      size = 2.5, box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines")) +
      theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
    #p.pos
  }else{
    p.pos <- ggplot() + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                   color = "#80b1d3") +

      #scale_x_continuous(expand=c(0,0)) +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + #在0的位置画虚线
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +

      theme_bw() +
      theme(axis.text.y=element_text(size = 12, face = "bold"),
            panel.grid = element_blank())+  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  }

  #########################################################################
  ##  arrange all Parts together
  #####################################################################################
  rel_heights <- c(1.5, .5, .2, 1.5) # 上中下三个部分的比例
  plotlist <- list(p.res, p2, p.grad, p.pos)
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(size = 12, face = "bold"))

  plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

  # 保存到PDF文件
  # ggsave("GSEA_multi_pathways.pdf", width=6, height=5)
}




