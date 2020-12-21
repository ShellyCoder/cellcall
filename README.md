# cellwave: inference of intercellular networks from single-cell transcriptomics

<a name="wMi2P"></a>
## 1. Introduction to CellWave
<a name="I3XMe"></a>
### 1.1 workflow
The figure below shows a graphical representation of the CellWave workflow.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608018839628-0d742450-a88b-41dd-b712-84195dd8ea8a.png#align=left&display=inline&height=281&margin=%5Bobject%20Object%5D&name=image.png&originHeight=281&originWidth=846&size=81307&status=done&style=none&width=846)<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608018860688-c3991704-21a5-47a2-a6a0-82c09a4f9d85.png#align=left&display=inline&height=259&margin=%5Bobject%20Object%5D&name=image.png&originHeight=259&originWidth=846&size=129208&status=done&style=none&width=846)
<a name="zuiad"></a>
### 1.2 how to install R package
```
library(devtools)
devtools::install_github("shellylab/cellwave")
```
If you encounter the following error -- ERROR: dependency * are not available for package 'cellwave', installing * package manually to install dependency is a good choice.
<a name="v9J0k"></a>
## 2. Main functionalities of CellWave
Specific functionalities of this package include (use data included in the package):
<a name="GFKvX"></a>
### 2.1 assessing how well ligands expressed by a sender cell interact with the receptor of receiver cell.
<a name="fhhux"></a>
#### 2.1.1 load data
The colnames can't contain punctuation such as commas, periods, dashes, etc. Using underline to connect barcoder_celltype is recommended.
```
  f.tmp <- system.file("extdata", "example_Data.Rdata", package="cellwave")
  load(f.tmp)
  
  ## gene expression stored in the variable in.content
  dim(in.content)
  in.content[1:4, 1:4]
  table(str_split(colnames(in.content), "_", simplify = T)[,2])
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608020082823-4b9c7385-56cf-4eb0-89e8-ab76167005fc.png#align=left&display=inline&height=261&margin=%5Bobject%20Object%5D&name=image.png&originHeight=261&originWidth=741&size=22973&status=done&style=none&width=741)
<a name="d8DIn"></a>
#### 2.1.2 createobject
What's important is **the parameter** as followed:<br />**names.delim**  For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE_CELLTYPE, set this to "_" to separate the cell name into its component parts for picking the relevant field.<br />**source**  the type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM". If you don't want any transformation in the cellwave, using CPM data as input and set source = "CPM".<br />**Org**  choose the species source of gene, eg "Homo sapiens", "Mus musculus". This parameter matters following ligand-receptor-tf resource.
```
  mt <- CreateNichConObject(data=in.content, min.feature = 3,
                            names.field = 2,
                            names.delim = "_",
                            source = "TPM", # fullLength, UMI, TPM
                            scale.factor = 10^6,
                            Org = "Homo sapiens",
                            project = "Microenvironment")
```
What's in the **NichConObject ？**
```
mt@data$count: raw data
mt@data$withoutlog: data proceeded by cellwave
mt@meta.data: metadata of the data, sampleID, celltype, etc
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608031777700-0b7382f9-865b-43c0-9965-8b9c75906bf7.png#align=left&display=inline&height=318&margin=%5Bobject%20Object%5D&name=image.png&originHeight=318&originWidth=790&size=31888&status=done&style=none&width=790)
<a name="bulKQ"></a>
#### 2.1.3 compute the score
What's important is **the parameter** as followed:<br />**names.delim**  For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE_CELLTYPE, set this to "_" to separate the cell name into its component parts for picking the relevant field.
```
mt <- TransCommuProfile(mt,
                          pValueCor = 0.05,
                          CorValue = 0.1,
                          topTargetCor=1,
                          p.adjust = 0.05,
                          use.type="median",
                          probs = 0.75,
                          method="weighted",
                          Org = 'Homo sapiens')
```
What's new in the **NichConObject ？**<br />
```
mt@data$expr_l_r: raw score
mt@data$expr_l_r_log2: log2(raw score+1)
mt@data$expr_l_r_log2_scale: do max,min transform to expr_l_r_log2
mt@data$gsea.list: result of TF-activation
```
<a name="jDh7s"></a>
### 2.2 visualization of the ligand-receptor-TF model
<a name="ZpHeZ"></a>
#### 2.2.1 overveiw in circle plot 
set the color of each cell type
```
  cell_color <- data.frame(color=c('#e31a1c','#1f78b4',
                                   '#e78ac3','#ff7f00'), stringsAsFactors = FALSE)
  rownames(cell_color) <- c("SSC", "SPGing", "SPGed", "ST")
```
plot circle with cellwave object 
```
ViewInterCircos(object = mt, font = 2, cellColor = cell_color, lrColor = c("#F16B6F", "#84B1ED"),
                  arr.type = "big.arrow",arr.length = 0.04,
                  trackhight1 = 0.05, slot="expr_l_r_log2_scale",linkcolor.from.sender = TRUE,
                  linkcolor = NULL, gap.degree = 2,order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                  trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608035884476-ec0ec869-3246-4260-9b8b-d922cbb2683f.png#align=left&display=inline&height=402&margin=%5Bobject%20Object%5D&name=image.png&originHeight=402&originWidth=589&size=75656&status=done&style=none&width=589)<br />plot circle with DIY dataframe of mt@data$expr_l_r_log2_scale 
```
  ViewInterCircos(object = mt@data$expr_l_r_log2_scale, font = 2, cellColor = cell_color,
                  lrColor = c("#F16B6F", "#84B1ED"),
                  arr.type = "big.arrow",arr.length = 0.04,
                  trackhight1 = 0.05, slot="expr_l_r_log2_scale",linkcolor.from.sender = TRUE,
                  linkcolor = NULL, gap.degree = 2,order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                  trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = T)
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608036102646-62d254cb-af6d-4154-bdc6-e7c09d4c72d1.png#align=left&display=inline&height=401&margin=%5Bobject%20Object%5D&name=image.png&originHeight=401&originWidth=594&size=75665&status=done&style=none&width=594)
<a name="uCkIS"></a>
#### 2.2.2 overveiw in pheatmap plot 
```
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,show_colnames = T,
             treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 12,angle_col = "45",
             main="score")
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608036336849-fbfe4e8e-f542-4f07-8572-8206dc0d9b4e.png#align=left&display=inline&height=623&margin=%5Bobject%20Object%5D&name=image.png&originHeight=623&originWidth=1820&size=195839&status=done&style=none&width=1820)
<a name="uyYcQ"></a>
#### 2.2.3 inspect Ligand-Receptor-TF in specific cellA-cellB
There are three types to show this triple relation.<br />Funtion LR2TF to measure the triple relation between specific cells.
```
  mt <- LR2TF(object = mt, sender_cell="ST", recevier_cell="SSC",
              slot="expr_l_r_log2_scale", org="Homo sapiens")
  head(mt@reductions$sankey)
```
First type, function LRT.Dimplot.
```
  sank <- LRT.Dimplot(mt, fontSize = 8, nodeWidth = 30, height = NULL, width = 1200, sinksRight=FALSE, DIY.color = FALSE)
  saveNetwork(sank, "~/ST-SSC_full.html")
```
The first pillar is ligand，the second pillar is receptor，the last pillar is tf.<br />And the color of left and right flow is consistent with ligand and receptor respectively.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608037011420-8a08e2c5-73f0-4d9a-bc49-87119980fcda.png#align=left&display=inline&height=287&margin=%5Bobject%20Object%5D&name=image.png&originHeight=714&originWidth=1855&size=282995&status=done&style=none&width=746)

---

Second type, function sankey_graph with isGrandSon = FALSE.
```
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
```
sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], nudge_x = NULL,
font.size = 4, boder.col="white", isGrandSon = F)
```
The first pillar is ligand，the second pillar is receptor，the last pillar is tf.<br />And the color of left and right flow is consistent with ligand and receptor respectively.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608038131311-57d827cd-6a0f-4cb1-aa68-87ec48b80f47.png#align=left&display=inline&height=604&margin=%5Bobject%20Object%5D&name=image.png&originHeight=604&originWidth=1192&size=277773&status=done&style=none&width=1192)

---

Third type, function sankey_graph with isGrandSon = TRUE.
```
library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0)  ## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  ## transform weight1 and weight2 to one value (weight)

## set the color of node in sankey graph
mycol.vector = c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')
elments.num <-  length(unique(tmp.df$Ligand))
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))

sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], isGrandSon = TRUE,
              nudge_x = nudge_x, font.size = 2, boder.col="white", set_alpha = 0.8)
```
The first pillar is ligand，the second pillar is receptor，the last pillar is tf.<br />And the color of left and right flow is consistent with one node (ligand or receptor).<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608038695456-2ce50409-aaf1-4484-bfda-bde39ec7738f.png#align=left&display=inline&height=562&margin=%5Bobject%20Object%5D&name=image.png&originHeight=562&originWidth=1025&size=322055&status=done&style=none&width=1025)
<a name="4oHTG"></a>
### 2.3 enrichment in pathway
getHyperPathway to perform enrichment, getForBubble to merge data for graph and plotBubble produce the bubble plot.
```
n <- mt@data$expr_l_r_log2_scale

pathway.hyper.list <- lapply(colnames(n), function(i){
    print(i)
    tmp <- getHyperPathway(data = n, cella_cellb = i, Org="Homo sapiens")
    return(tmp)
})

myPub.df <- getForBubble(pathway.hyper.list)
p <- plotBubble(myPub.df)
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608044718538-54a82161-c714-4389-80f9-e715429c0754.png#align=left&display=inline&height=685&margin=%5Bobject%20Object%5D&name=image.png&originHeight=685&originWidth=1047&size=161481&status=done&style=none&width=1047)
<a name="toPeb"></a>
### 2.4 activation of TF in receiver cell
<a name="MY3Ne"></a>
#### 2.4.1 ridge plot
plot enrichment result of TF (filter or not). 
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
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608099926349-53777dc6-352b-484e-a3a2-9a4c519906a1.png#align=left&display=inline&height=686&margin=%5Bobject%20Object%5D&name=image.png&originHeight=686&originWidth=688&size=117664&status=done&style=none&width=688)
<a name="tumnl"></a>
#### 2.4.2 gsea plot
```
  ssc.tf <- names(mt@data$gsea.list$SSC@geneSets)
  ssc.tf
```
Show all TF have result in the SSC.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608100011863-2c349ce5-f197-4528-a167-2eec6004524c.png#align=left&display=inline&height=122&margin=%5Bobject%20Object%5D&name=image.png&originHeight=122&originWidth=1317&size=19698&status=done&style=none&width=1317)
```
  getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID=c("CREBBP", "ESR1", "FOXO3"), myCelltype="SSC",
              fc.list=mt@data$fc.list,  selectedGeneID = mt@data$gsea.list$SSC@geneSets$CREBBP[1:10],
              mycol = NULL)
```
![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608100142219-a62cfef9-fe81-4a2e-8222-73d30770a122.png#align=left&display=inline&height=682&margin=%5Bobject%20Object%5D&name=image.png&originHeight=682&originWidth=692&size=75484&status=done&style=none&width=692)<br />
<br />

