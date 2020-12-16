#' plot pheatmap graph with communication profile
#' @param object a Cellwave objects
#' @param slot plot the graph with the data of specific slot
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param treeheight_row the height of a tree for rows, if these are clustered. Default value 0 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered. Default value 50 points.
#' @param cluster_rows boolean values determining if rows should be clustered or hclust object,
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object.
#' @param fontsize base fontsize for the plot
#' @param angle_col angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315)
#' @param color vector of colors used in heatmap.
#' @param main the title of the plot, default is "score".
#' @importFrom pheatmap pheatmap
#' @export

viewPheatmap <- function(object, slot="expr_l_r_log2_scale", show_rownames = T, show_colnames = T,
                         treeheight_row=0, treeheight_col=50, cluster_rows = T, cluster_cols = F,
                         fontsize = 12, angle_col = "45", color=NULL, main="score"){
  # library(pheatmap)
  dat <- object@data[[slot]]
  pheatmap::pheatmap(dat, show_rownames = show_rownames, show_colnames = show_colnames,
           treeheight_row=treeheight_row, treeheight_col=treeheight_col,
           cluster_rows = cluster_rows, cluster_cols = cluster_cols, fontsize = fontsize, angle_col = angle_col,
           # color = colorRampPalette(colors = c('#FFFFFF','#ffffb2','#fd8d3c','#e31a1c'))(1000),
           main=main)
}



