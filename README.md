<a name="PNqrO"></a>
# CellCall: Integrating paired ligand-receptor and transcription factor activities for cell-cell communication


<a name="adf749d9"></a>
## Updated information of CellCall
<a name="58715c1d"></a>
#### 2021/08/02 -- The research of CellCall is online. Please cite us with https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkab638/6332819.
<a name="58715c1b"></a>
#### 2021/06/25 -- Update the manual more comprehensively.
<a name="58715c1a"></a>
#### 2021/04/15 -- Update the function getHyperPathway.
<a name="9b526a59"></a>
#### 2021/02/02 -- Increase the number of LR datasets reference to 1141.
<a name="79d2b28b"></a>
#### 2021/02/01 -- The reference LR datasets change to, core and extended, two parts. And rename the tool to CellCall.
<a name="6378b244"></a>
#### 2020/12/11 -- The tool CellWave is online.


<a name="c8f560dd"></a>
## 1. Introduction to CellCall


<a name="JqFcu"></a>
### 1.1 Overview of CellCall
CellCall is a toolkit to infer intercellular communication networks and internal regulatory signals by integrating intracellular and intercellular signaling. (1) CellCall collects ligand-receptor-transcript factor (L-R-TF) axis datasets based on KEGG pathways. (2) According to prior knowledge of L-R-TF interactions, CellCall infers intercellular communication by combining the expression of ligands/receptors and downstream TF activities for certain L-R pairs. (3) CellCall embeds a pathway activity analysis method to identify the crucial pathways involved in communications between certain cell types. (4) CellCall offers a rich suite of visualization options (Circos plot, Sankey plot, bubble plot, ridge plot, etc.) to intuitively present the analysis results. The overview figure of CellCall is shown as follows.

![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168160448-8ee48425-8406-45a1-8680-3c220176e73b.png#align=left&display=inline&height=634&margin=%5Bobject%20Object%5D&name=image.png&originHeight=634&originWidth=714&size=211419&status=done&style=none&width=714#id=XTMe1&originHeight=634&originWidth=714&originalType=binary&ratio=1&status=done&style=none)<br />

<a name="27dc18d1"></a>
### 1.2 Installing R package
To install this package, start R (version 3.6 or higher) and enter:
```
library(devtools)
devtools::install_github("ShellyCoder/cellcall")
```
If you encounter the following error -- ERROR: dependency is not available for package 'cellcall', install corresponding R package. And appropriate version is in the section 4.<br />

<a name="689d6784"></a>
## 2. Main functions of CellCall
CellCall provides a variety of functions including intercellular communication analysis, pathway activity analysis and a rich suite of visualization tools to intuitively present the results of the analysis (including Heatmap, Circos plot, Bubble plot, Sankey plot, TF enrichment plot and Ridge plot).<br />

<a name="LP04w"></a>
### 2.1 Intercellular communication analysis


<a name="945a6b52"></a>
#### 2.1.1 Load data
The format of the input file is as follow table:<br />1. The row names: gene symbols.<br />2. The column names: cell IDs. The colnames can't contain punctuation such as commas, periods, dashes, etc. Using underline to connect barcoder and cell type is recommended. Take the input format below as an example, the column name is made up of index and cell type. Users should set names.field=2 and names.delim="_" in the function CreateNichConObject(). After that, cell type information is obtained and stored in the S4 object for later analysis. Because method in this paper depends on the cell type information, obtaining celltype information correctly is important.<br />3. Other place: the expression values (counts or TPM) for a gene in a cell.

|  | 1_ST | 2_ST | 3_ST | 4_SSC | 5_SSC | 6_SPGed | 7_SPGed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TSPAN6 | 2.278 | 2.031 | 0.000 | 12.385 | 0.000 | 0.553 | 24.846 |
| TNMD | 9.112 | 6.031 | 0.000 | 0.000 | 11.615 | 10.518 | 0.000 |
| DPM1 | 0.000 | 0.000 | 21.498 | 4.246 | 7.382 | 0.000 | 2.385 |
| SCYL3 | 5.983 | 1.215 | 0.000 | 0.518 | 2.386 | 4.002 | 14.792 |

This instruction may take the in-house dataset included in the package as an example. User can load the dataset with command following in the code box. There are 366 single cells and 35,135 genes that were performed with the scRNA sequencing.
```r
f.tmp <- system.file("extdata", "example_Data.Rdata", package="cellcall")
load(f.tmp)

## gene expression stored in the variable in.content
dim(in.content)
in.content[1:4, 1:4]
table(str_split(colnames(in.content), "_", simplify = T)[,2])
```
We next use the expression dataframe  to create a CreateNichCon object with the function CreateNichConObject as the code in part 2.1.2. The object serves as a container that contains both data (like the expression dataframe) and analysis (like score, or enrichment results) for a single-cell dataset.<br />

<a name="e815f697"></a>
#### 2.1.2 Create object
CellCall use the expression dataframe to create an S4 object by the function CreateNichConObject, The line of code is shown in the code box. The S4 object serves as a container that contains both data (such as the expression dataframe) and analysis results (such score and enrichment results) for a project (see section 3 for details).
```r
 mt <- CreateNichConObject(data=in.content, min.feature = 3,
                            names.field = 2,
                            names.delim = "_",
                            source = "TPM",
                            scale.factor = 10^6,
                            Org = "Homo sapiens",
                            project = "Microenvironment")
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **data** | A dataframe with row of gene and column of sample and the value must be numeric. Meanwhile what matters is that the colnames of dataframe should be in line with the paramter 'names.delim' and 'names.field', the former for pattern to splite every colnames, the latter for setting which index in splited colnames is cell type information.<br />The function can get the 'CELLTYPE' information from the colnames 'BARCODE_CLUSTER_CELLTYPE' with names.delim="_" and names.field='3', and then stored in slot meta.data of CreateNichCon.<br />Cell type annotation from every cell is essential for scoring cell communication. If the colnames of data don't coincide with the paramter 'names.delim' and 'names.field', CreateNichCon object may fail to create. |
| **min.feature** | Include cells where enough features equalling min.feature are detected. It's a preprocess which is the same as Seurat and set min.feature=0, if you don't want to filter cell. This parameter depends on the sequencing technology of the input data. |
| **names.delim** | Set the pattern to splite column name into vector. If the column name of the input matrix is BARCODE_CLUSTER_CELLTYPE, set names.delim="_" to get CELLTYPE of BARCODE_CLUSTER_CELLTYPE with names.field=3. |
| **names.field** | Set the index of column name vector which is splited by parameter names.delim to get cell type information. If the column name of the input matrix is BARCODE_CLUSTER_CELLTYPE, set names.field=3 to get CELLTYPE of BARCODE_CLUSTER_CELLTYPE with names.delim="_". |
| **source** | The type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM". When the source of input data is  "TPM" or "CPM", no transformation on the data. Otherwise, we transform the data to TPM with the parameter source="fullLength" and to CPM with source="UMI". |
| **scale.factor** | Sets the scale factor for cell-level normalization, default "10^6", if the parameter is "UMI" or "fullLength". Otherwise this parameter doesn't work. |
| **Org** | Set the species source of gene, eg "Homo sapiens", "Mus musculus". This parameter decides the paired ligand-receptor dataset and the transcript length which is needed in "TPM" transformation. |
| **project** | Sets the project name for the CreateNichCon object. |


<br />

<a name="oEar5"></a>
#### 2.1.3 Infer the cell-cell communication score
The communication score of an L-R interaction between cell types is evaluated by integrating the L2- norm of the L-R interaction and the activity score of the downstream TFs. The code is shown in the code box.
```r
mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **ob​ject** | A Cellcall S4 object, the result of function CreateNichConObject(). |
| **pValueCor** | Set the threshold of spearson Correlation significance between target gene and TF, ( significance  < pValueCor, default is 0.05 ). |
| **CorValue** | Set the threshold of spearson Correlation Coefficient between target gene and TF, ( Coefficient > CorValue, default is 0.1 ). |
| **topTargetCor** | Set the rank of candidate genes which has firlter by spearson Correlation, default is 1, that means 100% filtered candidate genes will be used. |
| **p.adjust** | Set the threshold of regulons's GSEA pValue which adjusted by Benjamini & Hochberg, default is 0.05. |
| **use.type** | With parameter "median", CellCall set the mean value of gene as zero, when the percentile of gene expression in one celltype below the parameter "probs". The other choice is "mean" and means that we not concern about the percentile of gene expression in one celltype but directly use the mean value. |
| **probs** | Set the percentile of gene expression in one celltype to represent mean value, when use.type="median". |
| **method** | Choose the proper method to score downstream activation of all regulons of given ligand-receptor relation. Candidate values are "weighted", "max", "mean", of which "weighted" is default. |
| **Org** | Choose the dataset source of this project, eg "Homo sapiens", "Mus musculus". |
| **IS_core** | Logical variable, whether use core reference LR data with high confidence or include extended datasets on the basis of core reference. |

<a name="Hi4AT"></a>
### 2.2 Pathway activity analysis
CellCall embeds a pathway activity analysis method to help explore the main pathways involved in communication between certain cells. The code is shown in the code box.<br />

```r
n <- mt@data$expr_l_r_log2_scale

pathway.hyper.list <- lapply(colnames(n), function(i){
    print(i)
    tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
    return(tmp)
})
```
**getHyperPathway():**

| **Arguments** | **Detail** |
| --- | --- |
| **data** | A dataframe of communication score where row name is ligand-receptor and column names is cellA-cellB, stored in the data$expr_l_r_log2_scale slot of S4 object. |
| **ob​ject** | A Cellcall S4 object, the result of function CreateNichConObject() and TransCommuProfile(). |
| **cella_cellb** | If explore the pathway enriched by paired ligand-receptor dataset between sender cellA and receiver cellB, user can set cella_cellb="A-B". |
| **Org** | Choose the dataset source of this project, eg "Homo sapiens", "Mus musculus". |
| **IS_core** | Logical variable, whether use core reference LR data with high confidence or include extended datasets on the basis of core reference. |

For pathway activity analysis, Bubble plot is adopted to present the analysis results. Function of getForBubble is used to merge the data and plotBubble is used to draw the bubble plot.
```r
myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
p <- plotBubble(myPub.df)
```
**getForBubble():**

| **Arguments** | **Detail** |
| --- | --- |
| **pathway.hyper.list** | A list of enrichment result of function getHyperPathway(). |
| **cella_cellb** | If explore the pathway enriched by paired ligand-receptor dataset between sender cellA and receiver cellB, user can set cella_cellb="A-B". |

![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608044718538-54a82161-c714-4389-80f9-e715429c0754.png#align=left&display=inline&height=685&margin=%5Bobject%20Object%5D&name=image.png&originHeight=685&originWidth=1047&size=161481&status=done&style=none&width=1047#from=url&id=CmUMN&margin=%5Bobject%20Object%5D&originHeight=685&originWidth=1047&originalType=binary&ratio=2&status=done&style=none)
<a name="ztJhz"></a>
### 2.3 Visualization
CellCall offers a rich suite of visualization tools to intuitively present the results of the analysis, including heatmap, Circos plot, Bubble plot, Sankey plot, TF enrichment plot and Ridge plot.
<a name="TsPXv"></a>
#### 2.3.1 Circle plot
Circle plot is adopted to present the global cell-cell communications between cell types.<br />Setting the color and name of each cell type:
```r
  cell_color <- data.frame(color=c('#e31a1c','#1f78b4',
                                   '#e78ac3','#ff7f00'), stringsAsFactors = FALSE)
  rownames(cell_color) <- c("SSC", "SPGing", "SPGed", "ST")
```

<br />Plotting circle with CellCall object:
```r
ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 2,
                order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168608941-00510457-8fa4-425e-bb73-57b1c8718f01.png#align=left&display=inline&height=475&margin=%5Bobject%20Object%5D&name=image.png&originHeight=475&originWidth=559&size=74933&status=done&style=none&width=559#id=Z16r5&originHeight=475&originWidth=559&originalType=binary&ratio=1&status=done&style=none)

Plotting circle with DIY dataframe of mt@data$expr_l_r_log2_scale:
```r
  ViewInterCircos(object = mt@data$expr_l_r_log2_scale, font = 2, 
                  cellColor = cell_color,
                  lrColor = c("#F16B6F", "#84B1ED"),
                  arr.type = "big.arrow",arr.length = 0.04,
                  trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                  linkcolor.from.sender = TRUE,
                  linkcolor = NULL, gap.degree = 2,
                  order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                  trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = T)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168681351-261eaa95-e761-4f08-9e2f-76462fa223b4.png#align=left&display=inline&height=452&margin=%5Bobject%20Object%5D&name=image.png&originHeight=452&originWidth=534&size=72472&status=done&style=none&width=534#id=tqJ8P&originHeight=452&originWidth=534&originalType=binary&ratio=1&status=done&style=none)<br />**ViewInterCircos():**

| **Arguments** | **Detail** |
| --- | --- |
| **object** | A Cellcall S4 object, the result of function CreateNichConObject() and TransCommuProfile(). |
| **font** | The size of font. |
| **cellColor** | A color dataframe, rownames is cell type, value is color. |
| **lrColor** | A color vector denotes the color of ligand and receptor, containing two elements, default is c('#D92E27', "#35C6F4"). |
| **order.vector** | Default is null, a celltype vector with the order you want in the circle graph. |
| **trackhight1** | Height of the outer track. |
| **linkcolor.from.sender** | Logical value, whether the color of line correspond with color of sender cell. |
| **linkcolor** | One color you want link to be, only if parameter linkcolor.from.sender=FALSE. |
| **arr.type** | Type of the arrows, default value is big.arrow There is an additional option triangle. |
| **arr.length** | Length of the arrows, measured in 'cm'. If arr.type is set to big.arrow, the value is percent to the radius of the unit circle. |
| **DIY** | Logical value, if TRUE, the parameter object should be a dataframe, and set slot="expr_l_r_log2_scale". otherwise object should be a Cellwave objects. |
| **gap.degree** | Between two neighbour sectors. It can be a single value or a vector. If it is a vector, the first value corresponds to the gap after the first sector. |
| **trackhight2** | Height of the inner track. |
| **track.margin2** | Set the margin of current track, a numeric vector. |
| **slot** | Plot the graph with the data of specific slot |



<a name="c2fcc6f5"></a>
#### 2.3.2 Pheatmap plot
Pheatmap plot is adopted to present the detailed communication scores for the L-R interactions between different cell types.
```r
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 12,angle_col = "45", 	
             main="score")
```
**viewPheatmap():**

| **Arguments** | **Detail** |
| --- | --- |
| **object** | A Cellcall S4 object, the result of function CreateNichConObject() and TransCommuProfile(). |
| **slot** | Set the slot of data which is used to plot the graph. |
| **show_rownames** | Boolean parameter specifying if row names are be shown. |
| **show_colnames** | Boolean parameter specifying if column names are be shown. |
| **treeheight_row** | Set the height of a tree for rows, if these are clustered. Default value 0 points. |
| **treeheight_col** | Set the height of a tree for columns, if these are clustered. Default value 50 points. |
| **cluster_rows** | Boolean values determining if rows should be clustered. |
| **cluster_cols** | Boolean values determining if columns should be clustered. |
| **fontsize** | Base fontsize for the plot. |
| **angle_col** | Set the angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315). |
| **color** | Vector of colors used in heatmap. |
| **main** | Set the title of the plot, default is "score". |


<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168729707-830a042f-afbd-4a67-a0fa-f3804abc6aac.png#align=left&display=inline&height=477&margin=%5Bobject%20Object%5D&name=image.png&originHeight=477&originWidth=1185&size=112516&status=done&style=none&width=1185#id=EJbnc&originHeight=477&originWidth=1185&originalType=binary&ratio=1&status=done&style=none)<br />

<a name="b101b966"></a>
#### 2.3.3 Sankey plot
Sankey plot is adopted to present the detailed L-R-TF axis for the communications between different cell types.
```r
  mt <- LR2TF(object = mt, sender_cell="ST", recevier_cell="SSC",
              slot="expr_l_r_log2_scale", org="Homo sapiens")
  head(mt@reductions$sankey)
```
There are three ways for users to draw personalized Sankey plot:<br />(1) The color depends on ligand and TF (by function LRT.Dimplot)
```r
  if(!require(networkD3)){
  		BiocManager::install("networkD3")
  }
  
  sank <- LRT.Dimplot(mt, fontSize = 8, nodeWidth = 30, height = NULL, width = 1200, 		 
                      sinksRight=FALSE, DIY.color = FALSE)
  networkD3::saveNetwork(sank, "~/ST-SSC_full.html")
```
**LRT.Dimplot():**

| **Arguments** | **Detail** |
| --- | --- |
| **object** | A Cellcall S4 object, the result of function CreateNichConObject() and TransCommuProfile(). |
| **fontSize** | Set the font size of text in the graph. |
| **nodeWidth** | Set the node width of sankey graph. |
| **nodePadding** | Set the padding of node. |
| **height** | Set the height of graph, default is NULL. |
| **width** | Set the width of graph, default is 1200. |
| **sinksRight** | Boolean parameter. If TRUE, the last nodes are moved to the right border of the plot. |
| **DIY.color** | Boolean parameter. If TRUE, set the parameter color. DIY with your color-setting, default is FALSE. |
| **color.DIY** | A color dataframe, rownames is cell type, value is color, default is NULL. |

The first pillar is ligand，the second pillar is receptor and the last pillar is TF.<br />And the color of left and right flow are consistent with ligand and TF respectively.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608037011420-8a08e2c5-73f0-4d9a-bc49-87119980fcda.png#align=left&display=inline&height=287&margin=%5Bobject%20Object%5D&name=image.png&originHeight=714&originWidth=1855&size=282995&status=done&style=none&width=746#id=WOHOu&originHeight=714&originWidth=1855&originalType=binary&ratio=1&status=done&style=none)

---

(2) The color depends on ligand and receptor (by function sankey_graph with isGrandSon = FALSE)<br />

```r
library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0) ## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  ## transform weight1 and weight2 to one value (weight)
head(tmp.df)

## set the color of node in sankey graph
mycol.vector = c('#5d62b5','#29c3be','#f2726f','#62b58f','#bc95df', '#67cdf2', '#ffc533', '#5d62b5', '#29c3be')  
elments.num <-  tmp.df %>% unlist %>% unique %>% length()
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))
```


```r
sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], nudge_x = NULL, font.size = 4, boder.col="white", isGrandSon = F)
```

<br />The first pillar is ligand，the second pillar is receptor and the last pillar is TF.<br />And the color of left and right flow are consistent with ligand and receptor respectively.<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168925812-ab05a489-d7f4-497f-b517-939486808a11.png#align=left&display=inline&height=488&margin=%5Bobject%20Object%5D&name=image.png&originHeight=488&originWidth=547&size=183644&status=done&style=none&width=547#id=sLu3g&originHeight=488&originWidth=547&originalType=binary&ratio=1&status=done&style=none)

---

(3) The color depends on ligand (by function sankey_graph with isGrandSon = TRUE)<br />

```r
library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0)  ## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  ## transform weight1 and weight2 to one value (weight)

## set the color of node in sankey graph
mycol.vector = c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
elments.num <-  length(unique(tmp.df$Ligand))
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], 
             isGrandSon = TRUE, nudge_x = nudge_x, font.size = 2, boder.col="white", 			
             set_alpha = 0.8)
```

<br />The first pillar is ligand，the second pillar is receptor and the last pillar is TF.<br />And the color of left and right flow are consistent with ligand.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608038695456-2ce50409-aaf1-4484-bfda-bde39ec7738f.png#align=left&display=inline&height=562&margin=%5Bobject%20Object%5D&name=image.png&originHeight=562&originWidth=1025&size=322055&status=done&style=none&width=1025#id=I8hDo&originHeight=562&originWidth=1025&originalType=binary&ratio=1&status=done&style=none)<br />**sankey_graph():**

| **Arguments** | **Detail** |
| --- | --- |
| **df** | A dataframe with five or four columns depending on the parameter isGrandSon. |
| **axes** | If plot triple realtion of sankey, set axes=1:3, otherwise bipartite realtion is 1:2. Default 1:3. |
| **mycol** | A vector of character, denotes the color of each node. |
| **nudge_x** | A vector of numeric, denotes the horizontal position of each node label. |
| **font.size** | Set the font size of node label. |
| **boder.col** | Set the color of node border. |
| **isGrandSon** | If FALSE, the flow inherits it's source axe and only consider about relation instead of score. Otherwise every axe only inherit the first one axes ggtitle and consider about score. |
| **set_alpha** | Set the alpha of color in the node, a numeric bwtween 0-1. |

<a name="e4de6250"></a>
#### 2.3.4 TF enrichment plot
TF enrichment plot is adopted to present the TF activities in receiver cells.<br />**Obtain the gene sets (TGs of TF):**<br />For some biologists, they pay more attention to the TGs of the TF. This tool provide options to present or export the details of the TG list, stored in the NichConObject@data$gsea.list$cell_type@geneSets.
```r
mt@data$gsea.list$SSC@geneSets
```
The figure presents a part of result stored in the NichConObject@data$gsea.list$cell_type@geneSets.<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1622793918329-07ac2149-03ff-4243-ba2f-06345a112b69.png#clientId=uc7e6a175-7c8e-4&from=paste&height=235&id=bJOWi&margin=%5Bobject%20Object%5D&name=image.png&originHeight=469&originWidth=980&originalType=binary&ratio=1&size=114458&status=done&style=none&taskId=u62af5ea2-5c00-4ce1-8da3-2dac8a7047a&width=490)<br />**Show all TFs in the SSC:**
```r
  ssc.tf <- names(mt@data$gsea.list$SSC@geneSets)
  ssc.tf
```
<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608100011863-2c349ce5-f197-4528-a167-2eec6004524c.png#align=left&display=inline&height=122&margin=%5Bobject%20Object%5D&name=image.png&originHeight=122&originWidth=1317&size=19698&status=done&style=none&width=1317#id=ofHIX&margin=%5Bobject%20Object%5D&originHeight=122&originWidth=1317&originalType=binary&ratio=1&status=done&style=none)<br />**Draw the TF enrichment plot:**
```r
getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID=c("CREBBP", "ESR1", "FOXO3"), 
            myCelltype="SSC", fc.list=mt@data$fc.list,  
            selectedGeneID = mt@data$gsea.list$SSC@geneSets$CREBBP[1:10],
            mycol = NULL)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608100142219-a62cfef9-fe81-4a2e-8222-73d30770a122.png#align=left&display=inline&height=682&margin=%5Bobject%20Object%5D&name=image.png&originHeight=682&originWidth=692&size=75484&status=done&style=none&width=692#id=wOTtk&originHeight=682&originWidth=692&originalType=binary&ratio=1&status=done&style=none)<br />**getGSEAplot():**

| **Arguments** | **Detail** |
| --- | --- |
| gsea.list | A list of enrichment result from cellcall. |
| myCelltype | The cell type of receiver cell. |
| fc.list | The foldchange list in the cellcall object. |
| geneSetID | A character of TF symbol to draw enrichment plot, only significant activated can be inspected. |
| selectedGeneID | Default is NULL, label the position of specific gene in FC flow. |
| mycol | Set the color of each TF. The length is consistent with geneSetID. |

<a name="l1Mx3"></a>
#### 2.3.5 Ridge plot
Ridge plot is adopted to present FC distribution of TGs of activated TFs.
```
## gsea object
egmt <- mt@data$gsea.list$SSC

## filter TF
egmt.df <- data.frame(egmt)
head(egmt.df[,1:6])
flag.index <- which(egmt.df$p.adjust < 0.05)

ridgeplot.DIY(x=egmt, fill="p.adjust", showCategory=flag.index, core_enrichment = T,
                orderBy = "NES", decreasing = FALSE)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608099926349-53777dc6-352b-484e-a3a2-9a4c519906a1.png#align=left&display=inline&height=686&margin=%5Bobject%20Object%5D&name=image.png&originHeight=686&originWidth=688&size=117664&status=done&style=none&width=688#id=h5lZy&originHeight=686&originWidth=688&originalType=binary&ratio=1&status=done&style=none)**ridgeplot.DIY():**

| **Arguments** | **Detail** |
| --- | --- |
| **x** | A gseaResult object. |
| **showCategory** | Set the number of categories for plotting. |
| **fill** | Choose one of "pvalue", "p.adjust", "qvalue" attribute to be the color of ridge. |
| **core_enrichment** | Whether only using core_enriched genes. |
| **orderBy** | The order of the Y-axis, default is NES, or other colnames, eg: "ID", "Description", "setSize", "enrichmentScore", "p.adjust". |
| **decreasing** | Logical variable. Should the orderBy order be increasing or decreasing? |

<a name="kPWcN"></a>
### 2.4 Interface with Seurat
```r
test <- CreateObject_fromSeurat(Seurat.object=Seurat.object, 
                                slot="counts", 
                                cell_type="orig.ident",
                                data_source="UMI",
                                scale.factor = 10^6, 
                                Org = "Homo sapiens")
```
**Arguments:**

| **Arguments** | **Detail** |
| --- | --- |
| **Seurat.object** | The Seurat object which stores the expression matrix and cell type information. |
| **slot** | The name of slot which contains expression matrix, default "counts". |
| **cell_type** | The name of specific column which contains cell type information, default "orig.ident". |
| **data_source** | The type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM". |
| **scale.factor** | Sets the scale factor for cell-level normalization, default "10^6". The parameter is only for "UMI" or "fullLength", otherwise it doesn't work. |
| **Org** | Set the species source of gene, eg "Homo sapiens", "Mus musculus". This decide which ligand-receptor reference dataset be used. |

<a name="e56ue"></a>
## 3. Structure of S4 object
<a name="ad42O"></a>
### 3.1 all data slots stored in the S4 object of R package
More detailed results in the intermediate process are stored in the S4 object of R package for some customized analyses. And explaination of every slot  is in the table below.<br />


| **Slot(@data)** | **Detail** |
| --- | --- |
| **count** | Raw input matrix input in the function CreateNichConObject |
| **withoutlog** | The matrix transformed from raw input, TPM or CPM |
| **expr_mean** | The mean value of gene in different cell |
| **regulons_matrix** | The normalized enrichment value of tanscriptional factor in different cell |
| **gsea.list** | The enrichment result and target genes of tanscriptional factor in different cell |
| **fc.list** | The fold change value between specific cell type and others |
| **expr_r_regulons** | The sum value of normalized enrichment value of tanscriptional factor downstreaming specific receptor |
| **softmax_ligand** | Softmax value of the ligand expression across all cell types |
| **softmax_receptor** | Softmax value of the receptor expression across all cell types |
| **expr_l_r** | The score of ligand-receptor in cellA-cellB |
| **expr_l_r_log2** | The score of ligand-receptor in cellA-cellB with log transform |
| **expr_l_r_log2_scale** | The score of ligand-receptor in cellA-cellB with log transform and scale to [0,1] |
| **DistanceKEGG** | The distance between receptor and tf in one pathway in KEGG |

| **Slot(@meta.data)** | **Detail** |
| --- | --- |
| **sampleID** | The cell id of all cell. |
| **celltype** | The metadata of cell type information. |
| **nFeature** | The number of gene which value is greater than 0 in every cell. |
| **nCounts** | The sum of expression value in every cell. |

<a name="X7o5u"></a>
### 3.2 present or export the details of the TG list
For some biologists, they pay more attention to the TGs of the TF. This tool provide options to present or export the details of the TG list, stored in the NichConObject@data$gsea.list$cell_type@geneSets.
```r
mt@data$gsea.list$SSC@geneSets
```
The figure presents a part of result stored in the NichConObject@data$gsea.list$cell_type@geneSets.<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1622793918329-07ac2149-03ff-4243-ba2f-06345a112b69.png#clientId=uc7e6a175-7c8e-4&from=paste&height=235&id=ud5973156&margin=%5Bobject%20Object%5D&name=image.png&originHeight=469&originWidth=980&originalType=binary&ratio=1&size=114458&status=done&style=none&taskId=u62af5ea2-5c00-4ce1-8da3-2dac8a7047a&width=490)<br />​<br />
<a name="b4d8b003"></a>
## 4. My session info


```
• Session info ------------------------------------------------------------------------------------------------
setting  value
version  R version 3.6.0 (2019-04-26)
os       Windows 10 x64
system   x86_64, mingw32
ui       RStudio
language en
• Packages ----------------------------------------------------------------------------------------------------
package         * version    date       lib source
AnnotationDbi     1.48.0     2019-10-29 [1] Bioconductor
assertthat        0.2.1      2019-03-21 [1] CRAN (R 3.6.3)
backports         1.1.7      2020-05-13 [1] CRAN (R 3.6.3)
Biobase           2.46.0     2019-10-29 [1] Bioconductor
BiocGenerics      0.32.0     2019-10-29 [1] Bioconductor
BiocManager       1.30.10    2019-11-16 [1] CRAN (R 3.6.3)
BiocParallel      1.20.1     2019-12-21 [1] Bioconductor
bit               1.1-15.2   2020-02-10 [1] CRAN (R 3.6.2)
bit64             0.9-7      2017-05-08 [1] CRAN (R 3.6.2)
blob              1.2.1      2020-01-20 [1] CRAN (R 3.6.3)
callr             3.4.3      2020-03-28 [1] CRAN (R 3.6.3)
cellcall        * 0.0.0.9000 2021-02-01 [1] local
circlize          0.4.10     2020-06-15 [1] CRAN (R 3.6.3)
cli               2.0.2      2020-02-28 [1] CRAN (R 3.6.3)
clue              0.3-57     2019-02-25 [1] CRAN (R 3.6.3)
cluster           2.0.8      2019-04-05 [1] CRAN (R 3.6.0)
clusterProfiler   3.14.3     2020-01-08 [1] Bioconductor
colorspace        1.4-1      2019-03-18 [1] CRAN (R 3.6.3)
ComplexHeatmap    2.2.0      2019-10-29 [1] Bioconductor
cowplot           1.0.0      2019-07-11 [1] CRAN (R 3.6.3)
crayon            1.3.4      2017-09-16 [1] CRAN (R 3.6.3)
data.table        1.12.8     2019-12-09 [1] CRAN (R 3.6.3)
DBI               1.1.0      2019-12-15 [1] CRAN (R 3.6.3)
desc              1.2.0      2018-05-01 [1] CRAN (R 3.6.3)
devtools        * 2.3.1      2020-07-21 [1] CRAN (R 3.6.3)
digest            0.6.25     2020-02-23 [1] CRAN (R 3.6.3)
DO.db             2.9        2020-08-23 [1] Bioconductor
DOSE              3.12.0     2019-10-29 [1] Bioconductor
dplyr           * 1.0.0      2020-05-29 [1] CRAN (R 3.6.3)
ellipsis          0.3.1      2020-05-15 [1] CRAN (R 3.6.3)
enrichplot        1.6.1      2019-12-16 [1] Bioconductor
europepmc         0.4        2020-05-31 [1] CRAN (R 3.6.3)
fansi             0.4.1      2020-01-08 [1] CRAN (R 3.6.3)
farver            2.0.3      2020-01-16 [1] CRAN (R 3.6.3)
fastmatch         1.1-0      2017-01-28 [1] CRAN (R 3.6.0)
fgsea             1.12.0     2019-10-29 [1] Bioconductor
fs                1.4.2      2020-06-30 [1] CRAN (R 3.6.3)
generics          0.1.0      2020-10-31 [1] CRAN (R 3.6.3)
GetoptLong        1.0.2      2020-07-06 [1] CRAN (R 3.6.3)
ggalluvial        0.12.1     2020-08-10 [1] CRAN (R 3.6.0)
ggforce           0.3.2      2020-06-23 [1] CRAN (R 3.6.3)
ggplot2           3.3.2      2020-06-19 [1] CRAN (R 3.6.3)
ggplotify         0.0.5      2020-03-12 [1] CRAN (R 3.6.3)
ggraph            2.0.3      2020-05-20 [1] CRAN (R 3.6.3)
ggrepel           0.8.2      2020-03-08 [1] CRAN (R 3.6.3)
ggridges          0.5.2      2020-01-12 [1] CRAN (R 3.6.3)
GlobalOptions     0.1.2      2020-06-10 [1] CRAN (R 3.6.3)
glue              1.4.1      2020-05-13 [1] CRAN (R 3.6.3)
GO.db             3.10.0     2020-08-23 [1] Bioconductor
GOSemSim          2.12.1     2020-03-19 [1] Bioconductor
graphlayouts      0.7.0      2020-04-25 [1] CRAN (R 3.6.3)
gridBase          0.4-7      2014-02-24 [1] CRAN (R 3.6.3)
gridExtra         2.3        2017-09-09 [1] CRAN (R 3.6.3)
gridGraphics      0.5-0      2020-02-25 [1] CRAN (R 3.6.3)
gtable            0.3.0      2019-03-25 [1] CRAN (R 3.6.3)
hms               0.5.3      2020-01-08 [1] CRAN (R 3.6.3)
htmltools         0.5.0      2020-06-16 [1] CRAN (R 3.6.3)
htmlwidgets       1.5.1      2019-10-08 [1] CRAN (R 3.6.3)
httr              1.4.1      2019-08-05 [1] CRAN (R 3.6.3)
igraph            1.2.5      2020-03-19 [1] CRAN (R 3.6.3)
IRanges           2.20.2     2020-01-13 [1] Bioconductor
jsonlite          1.7.0      2020-06-25 [1] CRAN (R 3.6.3)
knitr             1.29       2020-06-23 [1] CRAN (R 3.6.3)
labeling          0.3        2014-08-23 [1] CRAN (R 3.6.0)
lattice           0.20-38    2018-11-04 [1] CRAN (R 3.6.0)
lifecycle         0.2.0      2020-03-06 [1] CRAN (R 3.6.3)
magrittr        * 2.0.1      2020-11-17 [1] CRAN (R 3.6.3)
MASS              7.3-51.4   2019-03-31 [1] CRAN (R 3.6.0)
Matrix            1.2-17     2019-03-22 [1] CRAN (R 3.6.0)
memoise           1.1.0      2017-04-21 [1] CRAN (R 3.6.3)
mnormt            1.5-7      2020-04-30 [1] CRAN (R 3.6.3)
munsell           0.5.0      2018-06-12 [1] CRAN (R 3.6.3)
networkD3       * 0.4        2017-03-18 [1] CRAN (R 3.6.3)
nlme              3.1-139    2019-04-09 [1] CRAN (R 3.6.0)
pheatmap          1.0.12     2019-01-04 [1] CRAN (R 3.6.3)
pillar            1.4.4      2020-05-05 [1] CRAN (R 3.6.3)
pkgbuild          1.0.8      2020-05-07 [1] CRAN (R 3.6.3)
pkgconfig         2.0.3      2019-09-22 [1] CRAN (R 3.6.3)
pkgload           1.1.0      2020-05-29 [1] CRAN (R 3.6.3)
plyr              1.8.6      2020-03-03 [1] CRAN (R 3.6.3)
png               0.1-7      2013-12-03 [1] CRAN (R 3.6.0)
polyclip          1.10-0     2019-03-14 [1] CRAN (R 3.6.0)
prettyunits       1.1.1      2020-01-24 [1] CRAN (R 3.6.3)
processx          3.4.3      2020-07-05 [1] CRAN (R 3.6.3)
progress          1.2.2      2019-05-16 [1] CRAN (R 3.6.3)
ps                1.3.3      2020-05-08 [1] CRAN (R 3.6.3)
psych             1.9.12.31  2020-01-08 [1] CRAN (R 3.6.3)
purrr             0.3.4      2020-04-17 [1] CRAN (R 3.6.3)
qvalue            2.18.0     2019-10-29 [1] Bioconductor
R6                2.4.1      2019-11-12 [1] CRAN (R 3.6.3)
RColorBrewer      1.1-2      2014-12-07 [1] CRAN (R 3.6.0)
Rcpp              1.0.4.6    2020-04-09 [1] CRAN (R 3.6.3)
remotes           2.2.0      2020-07-21 [1] CRAN (R 3.6.3)
reshape2          1.4.4      2020-04-09 [1] CRAN (R 3.6.3)
rjson             0.2.20     2018-06-08 [1] CRAN (R 3.6.0)
rlang             0.4.6      2020-05-02 [1] CRAN (R 3.6.3)
roxygen2        * 7.1.1      2020-06-27 [1] CRAN (R 3.6.3)
rprojroot         1.3-2      2018-01-03 [1] CRAN (R 3.6.3)
RSQLite           2.2.0      2020-01-07 [1] CRAN (R 3.6.3)
rstudioapi        0.11       2020-02-07 [1] CRAN (R 3.6.3)
rvcheck           0.1.8      2020-03-01 [1] CRAN (R 3.6.3)
S4Vectors         0.24.4     2020-04-09 [1] Bioconductor
scales            1.1.1      2020-05-11 [1] CRAN (R 3.6.3)
sessioninfo       1.1.1      2018-11-05 [1] CRAN (R 3.6.3)
shape             1.4.4      2018-02-07 [1] CRAN (R 3.6.0)
stringi           1.4.6      2020-02-17 [1] CRAN (R 3.6.2)
stringr         * 1.4.0      2019-02-10 [1] CRAN (R 3.6.3)
testthat          2.3.2      2020-03-02 [1] CRAN (R 3.6.3)
tibble            3.0.1      2020-04-20 [1] CRAN (R 3.6.3)
tidygraph         1.2.0      2020-05-12 [1] CRAN (R 3.6.3)
tidyr             1.1.0      2020-05-20 [1] CRAN (R 3.6.3)
tidyselect        1.1.0      2020-05-11 [1] CRAN (R 3.6.3)
triebeard         0.3.0      2016-08-04 [1] CRAN (R 3.6.3)
tweenr            1.0.1      2018-12-14 [1] CRAN (R 3.6.3)
urltools          1.7.3      2019-04-14 [1] CRAN (R 3.6.3)
usethis         * 1.6.1      2020-04-29 [1] CRAN (R 3.6.3)
vctrs             0.3.1      2020-06-05 [1] CRAN (R 3.6.3)
viridis           0.5.1      2018-03-29 [1] CRAN (R 3.6.3)
viridisLite       0.3.0      2018-02-01 [1] CRAN (R 3.6.3)
withr             2.4.1      2021-01-26 [1] CRAN (R 3.6.0)
xfun              0.19       2020-10-30 [1] CRAN (R 3.6.3)
xml2              1.3.2      2020-04-23 [1] CRAN (R 3.6.3)
[1] I:/R/R-3.6.0/library
```
