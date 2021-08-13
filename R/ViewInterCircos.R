#' plot circle graph with communication profile
#' @param object a Cellwave objects
#' @param font the size of font
#' @param cellColor a color dataframe, rownames is cell type, value is color
#' @param lrColor a color vector denotes the color of ligand and receptor, containing two elements, default is c('#D92E27', "#35C6F4")
#' @param order.vector default is null, a celltype vector with the order you want in the circle graph
#' @param trackhight1 Height of the outer track
#' @param trackhight2 Height of the inner track
#' @param linkcolor.from.sender logical value, whether the color of line correspond with color of sender cell
#' @param linkcolor one color you want link to be, only if parameter linkcolor.from.sender=FALSE
#' @param arr.type Type of the arrows, default value is big.arrow There is an additional option triangle
#' @param arr.length Length of the arrows, measured in 'cm'. If arr.type is set to big.arrow, the value is percent to the radius of the unit circle.
#' @param DIY logical value, if TRUE, the parameter object should be a dataframe, and set slot="expr_l_r_log2_scale". otherwise object should be a Cellwave objects.
#' @param slot plot the graph with the data of specific slot
#' @param gap.degree between two neighbour sectors. It can be a single value or a vector. If it is a vector, the first value corresponds to the gap after the first sector.
#' @param track.margin2 affect current track
#' @importFrom grid pushViewport unit upViewport viewport gpar
#' @importFrom graphics plot.new par
#' @importFrom gridBase gridOMI
#' @importFrom circlize circos.clear circos.par circos.initialize circos.trackPlotRegion circos.link get.cell.meta.data highlight.sector colorRamp2 rand_color
#' @importFrom stringr str_split
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr filter
#' @importFrom ComplexHeatmap draw Legend packLegend
#' @export

ViewInterCircos <- function(object, font = 2, cellColor ,lrColor = NULL, order.vector = NULL,
                            trackhight1 = 0.05, linkcolor.from.sender = TRUE, linkcolor = NULL,
                            arr.type = "big.arrow",arr.length = 0.04, DIY = TRUE, gap.degree = NULL,
                            trackhight2 = 0.032, track.margin2 = c(0.01,0.12),slot="expr_l_r_log2_scale")
{

  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                        just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)


  circos.clear()  

  # cell.padding: the padding between sector and other sector

  if(is.null(gap.degree)){
    circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.01,0,0.01,0))
  }else{
    circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.01,0,0.01,0), gap.degree=gap.degree)
  }

  # library(stringr)
  if(DIY){
    a <- colSums(object)
  }else{
    a <- colSums(object@data[[slot]])
  }


  b <- stringr::str_split(names(a), "-", simplify = T)
  c <- data.frame(b, stringsAsFactors = FALSE)
  c$x3 <- as.numeric(a)
  test <- c
  colnames(test) <- c( "cell_from", "cell_to", "n")

  test$id1 <- paste("sender", test$cell_from, sep = "_", 1:nrow(test))
  test$id2 <- paste("recevier", test$cell_to, sep = "_", 1:nrow(test))

  a_tmp <- test[,c(1,4)]
  colnames(a_tmp) <- c('celltype', 'id')
  a_tmp$arrowType <-  "sender"

  b_tmp <- test[,c(2,5)]
  colnames(b_tmp) <- c('celltype', 'id')
  b_tmp$arrowType <-  "recevier"

  ab_tmp <- rbind(a_tmp,b_tmp)

  ab_tmp <- data.frame(ab_tmp, stringsAsFactors = FALSE)
  ab_tmp <- ab_tmp[order(ab_tmp$celltype,ab_tmp$id,decreasing = TRUE),]

  sector_id <- ab_tmp$id
  fa = ab_tmp$id
  if(is.null(order.vector)){
    fa = factor(fa,levels = fa)
  }else{
    fa.df <- str_split(fa, "_", simplify = T) %>% as.data.frame() %>% magrittr::set_colnames(c('sender_or_receiver', 'clltype', "index_number"))
    my.levels <- do.call(rbind,lapply(order.vector, function(x){
      fa.df[which(fa.df$clltype==x),]
    })) %>% apply(1, function(x){
      paste(x[1],x[2],x[3],sep='_')
    }) %>% unlist %>% as.character()

    fa = factor(fa,levels = my.levels)
  }
  circos.initialize(factors = fa, xlim = c(0,1)) # 

  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = trackhight1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
      # print(sector.index)
      # print(xlim)
      # print(ylim)
    }
  )

  if(is.null(cellColor)){
    cell_color <- data.frame(color = rand_color(length(unique(ab_tmp$celltype)), luminosity = "light"), stringsAsFactors = FALSE)
    rownames(cell_color) <- unique(ab_tmp$celltype)
  }else{
    cell_color <- cellColor
  }
  cell_type <- unique(ab_tmp$celltype)
  for(i in 1:length(cell_type)){
    myCell_Type <-  cell_type[i]
    mySector_id <- as.character(ab_tmp[ab_tmp$celltype==myCell_Type,'id'])
    myColor <- as.character(cell_color[myCell_Type,1])
    highlight.sector(mySector_id, track.index = 1,
                     text = myCell_Type, text.vjust = -1,niceFacing = T, font = font, col = myColor)

  }

  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = trackhight2,
    bg.border = NA,
    track.margin = track.margin2,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )

  if(is.null(lrColor)){
    ligand_receptor <- c('#D92E27', "#35C6F4")  #color of sender and receiver
  }else{
    ligand_receptor <- lrColor
  }

  for(i in 1:length(cell_type)){
    myCell_Type <-  cell_type[i]
    mySector_id_tmp <- ab_tmp[ab_tmp$celltype==myCell_Type,c('arrowType', 'id')]
    my_ligand <- dplyr::filter(mySector_id_tmp, arrowType == 'sender')
    my_receptor <- dplyr::filter(mySector_id_tmp, arrowType == 'recevier')

    if(length(as.character(my_ligand[,'id']))>0){
      highlight.sector(as.character(my_ligand[,'id']), track.index = 2,
                       text = '', niceFacing = F,col = ligand_receptor[1],text.col = 'white')

    }

    if(length(as.character(my_receptor[,'id']))>0){
      highlight.sector(as.character(my_receptor[,'id']), track.index = 2,
                       text = '', niceFacing = F,col = ligand_receptor[2],text.col = 'white')

    }

  }


  test$weighted_n <- (test$n-min(test$n))/(max(test$n)-min(test$n))
  test <- test[order(test$weighted_n, decreasing = F),]
  # print(test$weighted_n)
  # test$weighted_n <- as.numeric(test$n/sum(test$n))

  min_n <- min(test$weighted_n)
  max_n <- max(test$weighted_n)
  mean_n <- mean(min_n, max_n)

  if(linkcolor.from.sender){

    color.function <- apply(cell_color, 1, function(x){
      col_fun = colorRamp2(c(0, 1), c("#FFFFFF", x))
    })

    for(i in 1:nrow(test)){
      my_Line_color = as.character(cell_color[test[i,1],1])
      print(test$weighted_n[i])
      circos.link(sector.index1 = test[i,'id1'],
                  point1 = c(0,1),
                  sector.index2 = test[i,'id2'],
                  point2 = c(0,1),
                  directional = 1,
                  arr.type = arr.type,
                  # arr.width = 0,
                  arr.length = arr.length,
                  col=color.function[[test[i,'cell_from']]](test$weighted_n[i])
      )

    }

    upViewport()
    list.obj <- list()
    # discrete
    lgd_points = Legend(at = c('ligand', "receptor"), type = "points", gap = unit(2, "mm"),
                        legend_gp = gpar(col = lrColor), title_position = "topleft",
                        title = "")
    # discrete
    lgd_lines = Legend(at = rownames(cell_color), type = "lines", gap = unit(2, "mm"),
                       legend_gp = gpar(col = cell_color$color, lwd = 2), title_position = "topleft",
                       title = "Cell type")

    list.obj <- c(list.obj, list(lgd_lines))
    list.obj <- c(list.obj, list(lgd_points))
    # # continuous
    for (f in color.function) {
      lgd_links = Legend(at = c(0, 0.5, 1), col_fun = f, gap = unit(1, "mm"),direction="horizontal",
                         title_position = "topleft", title = "")
      list.obj <- c(list.obj, lgd_links)
    }

    lgd_list_vertical = packLegend(list = list.obj,
                                   gap = unit(2, "mm") # ligend distance
    )

    draw(lgd_list_vertical, x = circle_size, just = "left")
  }else{

    if(is.null(linkcolor)){
      col_fun = colorRamp2(c(0, 1), linkcolor)
    }else{
      col_fun = colorRamp2(c(0, 1), c("#FFFFFF", "#f349eb"))
    }

    for(i in 1:nrow(test)){
      my_Line_color = as.character(cell_color[test[i,1],1])
      print(test$weighted_n[i])
      circos.link(sector.index1 = test[i,'id1'],
                  point1 = c(0,1),
                  sector.index2 = test[i,'id2'],
                  point2 = c(0,1),
                  directional = 1,
                  arr.type = arr.type,
                  # arr.width = 0,
                  arr.length = arr.length,
                  col=col_fun(test$weighted_n[i])
      )

    }
    upViewport()

    # discrete
    lgd_points = Legend(at = c('ligand', "receptor"), type = "points", gap = unit(2, "mm"),
                        legend_gp = gpar(col = lrColor), title_position = "topleft",
                        title = "")
    # discrete
    lgd_lines = Legend(at = rownames(cell_color), type = "lines", gap = unit(2, "mm"),
                       legend_gp = gpar(col = cell_color$color, lwd = 2), title_position = "topleft",
                       background="#FFFFFF", title = "Cell type")
    # # continuous
    lgd_links = Legend(at = c(0, 0.5, 1), col_fun = col_fun, gap = unit(2, "mm"), direction="vertical",
                       background="#FFFFFF", title_position = "topleft", title = "Adjusted Score")

    lgd_list_vertical = packLegend(lgd_lines, lgd_points, lgd_links,
                                   gap = unit(4, "mm") # ligend distance
                                   )

    draw(lgd_list_vertical, x = circle_size, just = "left")
  }


}




