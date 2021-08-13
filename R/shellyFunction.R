#' get score of Triple relation corresponding specific cellA-cellB
#' @param df a dataframe consist of five columns: "Ligand", "Receptor", "TF", "weight1", "weight2"
#' @return a dataframe consist of four columns: "Ligand", "Receptor", "TF", "weight"
#' @importFrom dplyr filter
#' @export

trans2tripleScore <- function(df){
  colnames(df) <- c("Ligand", "Receptor", "TF", "weight1", "weight2")
  df <- dplyr::filter(df, weight1 !=0 & weight2 !=0)
  df$Ligand <- paste("sender:", df$Ligand, sep = "")
  df$Receptor <- paste("receiver:", df$Receptor, sep = "")
  df$TF <- paste("TF:", df$TF, sep = "")
  a <- df
  sum_tf <- aggregate(a[,5],by=list(a$Ligand,a$Receptor),FUN=sum) 
  colnames(sum_tf) <- c("Ligand", "Receptor", 'Score')
  
  aa <- c()
  for (i in 1:nrow(a)) {
    tmp = dplyr::filter(sum_tf, Ligand == a[i,1] & Receptor == a[i,2])
    res_tmp = as.numeric(a[i,5])/as.numeric(tmp[,3][1])
    aa <- c(aa, res_tmp)
  }

  a$weight2 <- as.numeric(aa)
  a$value <- a$weight1*a$weight2
  a <- a[,c(1,2,3,6)]

  return(a)
}
