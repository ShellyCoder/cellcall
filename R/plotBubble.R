#' plot bubble graph
#' @param dat a dataframe of melted communication score with columns: pathway, cc, pvalue, NES, logPvalue
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradientn scale_size labs theme element_text element_rect
#' @export

plotBubble <- function(dat){
    myPub.df <- dat
    ## 提前对数据进行处理，来印射size属性
    myPub.df$logPvalue <- -log10(myPub.df$pvalue)
    myPub.df$logPvalue[ myPub.df$logPvalue>=2 ] = 2
    myPub.df$logPvalue[ myPub.df$logPvalue>=1.3 & myPub.df$logPvalue<2 ] = 1.3
    myPub.df$logPvalue[ myPub.df$logPvalue<1.3 ] = 1

    p = ggplot(myPub.df, aes(x = cc, y = pathway))
    # 修改点的大小
    # 展示三维数据
    pbubble = p + geom_point(aes(size=logPvalue,colour=NES))
    # 设置渐变色
    # ls("package:ggplot2", pattern="^scale_fill.+")

    pr = pbubble+scale_colour_gradientn(colours=c('#1400FE', '#009FFE', '#D00B5D', '#FB6203','#F60404'),
                                        values = c(0,1/4,1/2,3/4,1)) +
      scale_size("Pvalue", breaks=c(1,1.3,2), range=c(1,5),
                 labels = c("0.1", "0.05", "<0.01")) ## 设置几个刻度, 1, 1.3, 2, 3

    # 绘制p气泡图
    pr = pr + labs(color="NES",
                   size="-log10 p-value",
                   x="CC",y="Pathway",title="")

    ## 修改x轴的文字  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
    pr = pr + theme(axis.text.x = element_text(size = 12, color = "black", face = "plain", vjust = 1.0, hjust = 1, angle = 45),
                    panel.background = element_rect(fill = "white", colour = "black", size = 1))
    return(pr)
}


