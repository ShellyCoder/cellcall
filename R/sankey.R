#' plot circle graph with communication profile
#' @param df a dataframe with five or four columns depending on the parameter isGrandSon
#' @param axes if triple realtion of sankey, set axes=1:3, otherwise bi is 1:2. default 1:3
#' @param mycol a vector of character, denotes the color of each node
#' @param nudge_x a vector of numeric, denotes the horizontal position of each node label
#' @param font.size the font size of node label
#' @param boder.col the color of node border
#' @param isGrandSon If FALSE, the third axe inherits the first two axes and only consider about relation instead of score. otherwise every axe only inherit the first one axes ggtitle and consider about score.
#' @param set_alpha set the alpha of color in the node, a numeric bwtween 0-1.
#' @importFrom ggplot2 ggplot aes geom_text scale_x_continuous theme element_blank element_text ggtitle after_stat theme_bw xlab ylab scale_fill_manual guides scale_x_discrete
#' @importFrom ggalluvial geom_flow geom_stratum to_lodes_form StatStratum
#' @export

sankey_graph <- function(df, axes, mycol=NULL, nudge_x=NULL, font.size=5,
                         boder.col = "black", isGrandSon = FALSE, set_alpha = 1){

  options(stringsAsFactors = FALSE)
  diy_stratum <- ggalluvial::StatStratum

if(isGrandSon){
  
  subdf <- df
  colnames(subdf) <- c("Ligand_symbol", "Receptor_symbol", "TF", "value")
  subdf <- as.data.frame(subdf)

  if(is.null(mycol)){
    p <- ggplot(as.data.frame(subdf),
                aes(y = value,
                    axis1 = Ligand_symbol,
                    axis2 = Receptor_symbol,
                    axis3 = TF)) +
      ggalluvial::geom_flow(stat = "alluvium",width = 1/8,aes(fill = Ligand_symbol), alpha=set_alpha) +
      ggalluvial::geom_stratum(width = 1/8, reverse = T,alpha = .9, size =0.001) +
      geom_text(stat = diy_stratum, size = font.size, aes(label = after_stat(stratum)),
                reverse = T) +
      scale_x_continuous(breaks = 1:3, labels = c("Ligand", "Receptor", "TF")) +

      theme(legend.position = "bottom",
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(size = font.size, face = "bold", color = "black")) +

      xlab("") + ylab("") +
      theme_bw() +
      theme(panel.grid =element_blank()) +
      theme(panel.border = element_blank()) + 
      theme(axis.line = element_blank(),axis.ticks = element_blank(),
            axis.text.y = element_blank()) + 
      ggtitle("")
  }else{
    p <- ggplot(as.data.frame(subdf),
                aes(y = value,
                    axis1 = Ligand_symbol,
                    axis2 = Receptor_symbol,
                    axis3 = TF)) + 
      scale_fill_manual(values = mycol) +
      ggalluvial::geom_flow(stat = "alluvium",width = 1/8,aes(fill = Ligand_symbol), alpha=set_alpha) +
      ggalluvial::geom_stratum(width = 1/8, reverse = T,alpha = .9, size =0.001) +
      geom_text(stat = diy_stratum, size = font.size, aes(label = after_stat(stratum)),
                reverse = T) +
      scale_x_continuous(breaks = 1:3, labels = c("Ligand", "Receptor", "TF")) +

      theme(legend.position = "bottom", 
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(size = font.size, face = "bold", color = "black")) +

      xlab("") + ylab("") +
      theme_bw() + 
      theme(panel.grid =element_blank()) +
      theme(panel.border = element_blank()) + 
      theme(axis.line = element_blank(),axis.ticks = element_blank(), 
            axis.text.y = element_blank()) + 
      ggtitle("")
  }

}else{
  UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                             axes = axes,
                             id = "Cohort")

  if(is.null(mycol)){
    p<-ggplot(UCB_lodes,
              aes(x = x, stratum = stratum, alluvium = Cohort,
                  fill = stratum, label = stratum)) +
      scale_x_discrete(expand = c(0, 0)) +
      ggalluvial::geom_flow(width = 1/8, alpha=set_alpha) + 
      ggalluvial::geom_stratum(alpha = .9,width = 1/6, size =0.001, color = boder.col) + 
      xlab("") + ylab("") +
      theme_bw() + 
      theme(panel.grid =element_blank()) + 
      theme(panel.border = element_blank()) + 
      theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) +
      ggtitle("")+
      guides(fill = FALSE)
  }else{
    p<-ggplot(UCB_lodes,
              aes(x = x, stratum = stratum, alluvium = Cohort,
                  fill = stratum, label = stratum)) +
      scale_x_discrete(expand = c(0, 0)) +
      ggalluvial::geom_flow(width = 1/8, alpha=set_alpha) +
      ggalluvial::geom_stratum(alpha = .9,width = 1/6, size =0.001, color = boder.col) + 

      scale_fill_manual(values = mycol) +

      xlab("") + ylab("") +
      theme_bw() + 
      theme(panel.grid =element_blank()) +
      theme(panel.border = element_blank()) +
      theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + 
      ggtitle("")+
      guides(fill = FALSE)
  }

  if(is.null(nudge_x)){
    p <- p + geom_text(stat = diy_stratum, size = font.size, color="black")
  }else{
    p <- p + geom_text(stat = diy_stratum, size = font.size, color="black", vjust="right", nudge_x=nudge_x)
  }
}
  return(p)
}

