#' get score of Triple relation corresponding specific cellA-cellB
#' @param object a Cellwave objects or a dataframe about score of Triple relation corresponding specific cellA-cellB
#' @param fontSize the font size of text.
#' @param nodeWidth the node Width of sankey graph.
#' @param nodePadding the padding of node.
#' @param height the height of graph, default is NULL.
#' @param width the width of graph, default is 1200.
#' @param sinksRight boolean. If TRUE, the last nodes are moved to the right border of the plot.
#' @param DIY.color boolean. If TRUE, set the parameter color.DIY with your color-setting, default is FALSE.
#' @param color.DIY a color dataframe, rownames is cell type, value is color, default is NULL.
#' @importFrom networkD3 sankeyNetwork JS saveNetwork
#' @importFrom dplyr filter
#' @importFrom stats aggregate
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom jsonlite toJSON
#' @export

LRT.Dimplot <- function(object, fontSize = 12, nodeWidth = 30,nodePadding = 5,
                        height = NULL, width = 1200, sinksRight=FALSE,
                        DIY.color = FALSE, color.DIY = NULL){
  is_myObject <- is(object, "CellInter")
  # if(!is_myObject){
  #   stop("object should be CellInter type.")
  # }
  if(!is_myObject){
    df <- object
  }else{
    df <- object@reductions$sankey
  }

  if(is.numeric(df$weight1) & is.numeric(df$weight2)){
  }else{
    stop("The weight1 and weight2 should be numeric")
  }

  df <- dplyr::filter(df, weight1 !=0 & weight2 !=0)
  df$Ligand <- paste("sender:", df$Ligand, sep = "")
  df$Receptor <- paste("receiver:", df$Receptor, sep = "")
  df$TF <- paste("TF:", df$TF, sep = "")
  a <- df
  sum_tf <- aggregate(a[,5],by=list(a$Ligand,a$Receptor),FUN=sum) # 统计每对ligand-receptor对应所有TF的和
  colnames(sum_tf) <- c("Ligand", "Receptor", 'Score')
  # aa <- apply(a, 1, function(x){
  #   tmp = dplyr::filter(sum_tf, Ligand == x[1] & Receptor == x[2])
  #   res_tmp = as.numeric(x[1,5])/as.numeric(tmp[,3][1])
  #
  #   print(class(x[1,5]))
  #   return(res_tmp)
  # })
  aa <- c()
  for (i in 1:nrow(a)) {
    tmp = dplyr::filter(sum_tf, Ligand == a[i,1] & Receptor == a[i,2])
    res_tmp = as.numeric(a[i,5])/as.numeric(tmp[,3][1])
    aa <- c(aa, res_tmp)
  }

  a$weight2 <- as.numeric(aa)
  a$value <- a$weight1*a$weight2
  a <- a[,c(1,2,3,6)]

  Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限  options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor

  # library(networkD3)
  # 加载数据
  mydata1 <- data.frame(ligand=as.character(a$Ligand),
                        receptor=as.character(a$Receptor),
                        tf=as.character(a$TF),
                        value=as.numeric(a$value),stringsAsFactors=FALSE)
  mydata.raw <- mydata1

  if(DIY.color){
    color.df <- color.DIY# 输入一个如下的color对应矩阵
    #  node   color
    #  CCL1   #223D6C
    #  NFKBIA #D20A13
    #  FGF20  #FFD121
    #  CSF2   #088247

  }else{
    # library(RColorBrewer)
    cellInter.color <- colorRampPalette(RColorBrewer::brewer.pal(n = 12, name = "Paired"))(24)
    ligand <- unique(mydata.raw$ligand)
    receptor <- unique(mydata.raw$receptor)
    tf <- unique(mydata.raw$tf)
    color.df <- data.frame(node = c(ligand, receptor, tf),
                           color = c(sample(cellInter.color,length(ligand),replace = TRUE),
                                     rep('#C3AFD3',length(receptor),replace = TRUE),
                                     sample(cellInter.color,length(tf),replace = TRUE) )
    )
  }

  mydata2<-aggregate(mydata1[,4],by=list(mydata1$ligand,mydata1$receptor),FUN=sum)
  names(mydata2) <- c("source","target","value")
  mydata2$group <- mydata2$source

  mydata1 <- mydata1[,-1]
  names(mydata1) <- c("source","target","value")
  mydata1$group <- mydata1$target

  fist_index <- nrow(mydata2)
  second_index <- nrow(mydata1)
  mydata <- rbind(mydata2,mydata1)

  Sankeylinks <- mydata
  Sankeynodes <- data.frame(name=unique(c(Sankeylinks$source,Sankeylinks$target)))
  Sankeynodes$index <- 0:(nrow(Sankeynodes) - 1)
  Sankeylinks <- merge(Sankeylinks,Sankeynodes,by.x="source",by.y="name")  #替换source列
  Sankeylinks <- merge(Sankeylinks,Sankeynodes,by.x="target",by.y="name")  #替换target列
  Sankeylinks <- merge(Sankeylinks,Sankeynodes,by.x="group",by.y="name")  #替换group列
  Sankeydata <- Sankeylinks[,c(5,6,4,7)]  #只保留索引值列，用于绘图
  names(Sankeydata) <- c("Source","Target","Value","group")
  Sankeydata$group <- as.character(Sankeydata$group)
  Sankeyname <- Sankeynodes[,1,drop=FALSE]

  Sankeyname$group <-as.character(0:(nrow(Sankeyname)-1))

  myColor <- merge(Sankeyname, color.df, by.x = "name", by.y = "node")
  color_scale <- data.frame(
    range = as.character(myColor$color),
    domain = as.character(myColor$group),
    # nodes = energy$nodes,
    stringsAsFactors = FALSE
  )

  sank <- sankeyNetwork(Links=Sankeydata,
                        Nodes=Sankeyname,
                        NodeGroup = "group",
                        LinkGroup = "group",
                        Source ="Source",
                        Target = "Target",
                        Value = "Value",
                        NodeID = "name",
                        colourScale = networkD3::JS( sprintf(
                          'd3.scaleOrdinal()
                        .domain(%s)
                        .range(%s)',
                          jsonlite::toJSON(color_scale$domain),
                          jsonlite::toJSON(color_scale$range)
                        )),
                        fontSize = fontSize,
                        nodeWidth = nodeWidth,
                        nodePadding = nodePadding,
                        height = height,
                        width = width,
                        sinksRight=sinksRight
  )

  # saveNetwork(sank, "test.pdf")
  return(sank)
}


