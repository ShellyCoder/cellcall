#' plot bubble graph
#' @param dat a dataframe of melted communication score with columns: pathway, cc, pvalue, NES, logPvalue
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradientn scale_size labs theme element_text element_rect
#' @export

plotBubble <- function(dat){
    myPub.df <- dat
    myPub.df$logPvalue <- -log10(myPub.df$pvalue)
    myPub.df$logPvalue[ myPub.df$logPvalue>=2 ] = 2
    myPub.df$logPvalue[ myPub.df$logPvalue>=1.3 & myPub.df$logPvalue<2 ] = 1.3
    myPub.df$logPvalue[ myPub.df$logPvalue<1.3 ] = 1

    p = ggplot(myPub.df, aes(x = cc, y = pathway))

    pbubble = p + geom_point(aes(size=logPvalue,colour=NES))


    pr = pbubble+scale_colour_gradientn(colours=c('#1400FE', '#009FFE', '#D00B5D', '#FB6203','#F60404'),
                                        values = c(0,1/4,1/2,3/4,1)) +
      scale_size("Pvalue", breaks=c(1,1.3,2), range=c(1,5),
                 labels = c("0.1", "0.05", "<0.01")) 

    pr = pr + labs(color="NES",
                   size="-log10 p-value",
                   x="CC",y="Pathway",title="")

    pr = pr + theme(axis.text.x = element_text(size = 8, color = "black", face = "plain", vjust = 1.0, hjust = 1, angle = 45),
                    panel.background = element_rect(fill = "white", colour = "black", size = 1))
    return(pr)
}


