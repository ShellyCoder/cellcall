<a name="PNqrO"></a>
# CellCall: Integrating paired ligand-receptor and transcription factor activities for cell-cell communication


<a name="adf749d9"></a>
## Updated information of CellCall
<a name="58715c1a"></a>
#### 2021/04/15 -- Change the function getHyperPathway for the debug of factor variable.
<a name="9b526a59"></a>
#### 2021/02/02 -- Increase the number of LR datasets reference to 1141.
<a name="79d2b28b"></a>
#### 2021/02/01 -- The reference LR datasets change to, core and extended, two parts. And rename the tool CellCall not CellWave.
<a name="6378b244"></a>
#### 2020/12/11 -- The tool CellWave is online.


<a name="c8f560dd"></a>
## 1. Introduction to CellCall


<a name="a22c6860"></a>
### 1.1 workflow

<br />The figure below shows a graphical representation of the CellCall workflow.<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168160448-8ee48425-8406-45a1-8680-3c220176e73b.png#align=left&display=inline&height=634&margin=%5Bobject%20Object%5D&name=image.png&originHeight=634&originWidth=714&size=211419&status=done&style=none&width=714#id=XTMe1&originHeight=634&originWidth=714&originalType=binary&status=done&style=none)<br />

<a name="27dc18d1"></a>
### 1.2 how to install R package


```
library(devtools)
devtools::install_github("ShellyCoder/cellcall")
```

<br />If you encounter the following error -- ERROR: dependency _ are not available for package 'cellcall', installing _ package manually to install dependency, in refer to the namespace file,  is a good choice.<br />

<a name="689d6784"></a>
## 2. Main functionalities of CellCall

<br />Specific functionalities of this package include (use data included in the package):<br />

<a name="e8321f2b"></a>
### 2.1 assessing how well ligands expressed by a sender cell interact with the receptor of receiver cell.


<a name="945a6b52"></a>
#### 2.1.1 load data
The format of the counts file is as follow table:<br />1. The row names: gene symbols.<br />2. The column names: cell IDs. The colnames can't contain punctuation such as commas, periods, dashes, etc. Using underline to connect barcoder_celltype is recommended.<br />3. Other place: the expression values (counts or TPM) for a gene in a cell

|  | 1_ST | 2_ST | 3_ST | 4_SSC | 5_SSC | 6_SPGed | 7_SPGed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TSPAN6 | 2.278 | 2.031 | 0.000 | 12.385 | 0.000 | 0.553 | 24.846 |
| TNMD | 9.112 | 6.031 | 0.000 | 0.000 | 11.615 | 10.518 | 0.000 |
| DPM1 | 0.000 | 0.000 | 21.498 | 4.246 | 7.382 | 0.000 | 2.385 |
| SCYL3 | 5.983 | 1.215 | 0.000 | 0.518 | 2.386 | 4.002 | 14.792 |

This instruction may take the in-house dataset included in the package as an example. User can load the dataset with command following in the code box. There are 366 single cells and 35,135 genes that were performed with the scRNA sequencing.<br />​<br />
```
  f.tmp <- system.file("extdata", "example_Data.Rdata", package="cellcall")
  load(f.tmp)
  
  ## gene expression stored in the variable in.content
  dim(in.content)
  in.content[1:4, 1:4]
  table(str_split(colnames(in.content), "_", simplify = T)[,2])
```
​

We next use the expression dataframe  to create a CreateNichCon object with the function CreateNichConObject as the code in part 2.1.2. The object serves as a container that contains both data (like the expression dataframe) and analysis (like score, or enrichment results) for a single-cell dataset.<br />

<a name="e815f697"></a>
#### 2.1.2 createobject
```r
 mt <- CreateNichConObject(data=in.content, min.feature = 3,
                            names.field = 2,
                            names.delim = "_",
                            source = "TPM", # fullLength, UMI, TPM
                            scale.factor = 10^6,
                            Org = "Homo sapiens",
                            project = "Microenvironment")
```
What's important is **the parameter setting** as followed:

1. **data**

A dataframe with row of gene and column of sample and the value must be numeric. Meanwhile what matters is that the colnames of dataframe should be in line with the paramter 'names.delim' and 'names.field', the former for pattern to splite every colnames, the latter for setting which index in splited colnames is cell type information. A names.delim, "_", and a names.field, '3', get the information 'CELLTYPE' from the colnames, 'BARCODE_CLUSTER_CELLTYPE',  stored in CreateNichCon@meta.data. If the colnames of data don't coincide with the paramter 'names.delim' and 'names.field', CreateNichCon object may fail to create.

2. **min.feature**

Include cells where at least this many features are detected. It's a preprocess which is the same as Seurat and set min.feature=0, if you don't want to filter cell. This parameter depends on the sequencing technology of the input data, e.g. Smart-seq2 and 10x Genomics.

3. **names.field**

For the initial identity class for each cell, choose this field from the cell's name. E.g. If your cells are named as BARCODE_CLUSTER_CELLTYPE in the input matrix, set names.field to 3 to set the initial identities to CELLTYPE.

4. **names.delim**

For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE_CELLTYPE, set this to "_" to separate the cell name into its component parts for picking the relevant field.

5. **project**

Sets the project name for the CreateNichCon object.

6. **source**

The type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM". When the source of input data is  "TPM" or "CPM", no transformation on the data. Otherwise, we transfer the data to "TPM" with "fullLength" data source and to "CPM" with "UMI" data source.

7. **scale.factor**

Sets the scale factor for cell-level normalization, default "10^6", if the parameter is "UMI" or "fullLength". Otherwise this parameter doesn't work.

8. **Org**

Set the species source of gene, eg "Homo sapiens", "Mus musculus". This decide which ligand-receptor reference dataset we use and considering the transcript length for "TPM" transformation.<br />
<br />What's in the **NichConObject **after create object**？Please refer to the part 3.1 for more detailed information.**

| **slot** | **detail** |
| --- | --- |
| count | raw input matrix input in the function CreateNichConObject |
| withoutlog | the matrix transformed from raw input, TPM or CPM |
| meta.data | metadata of the data, sampleID, celltype, etc |



<a name="74f90243"></a>
#### 2.1.3 compute the score


```
mt <- TransCommuProfile(mt,
                          pValueCor = 0.05,
                          CorValue = 0.1,
                          topTargetCor=1,
                          p.adjust = 0.05,
                          use.type="median",
                          probs = 0.75,
                          method="weighted",
                          IS_core = TRUE,
                          Org = 'Homo sapiens')
```

<br />What's new in the **NichConObject after **TransCommuProfile function**？Please refer to the part 3.1 for more detailed information.**

| **slot** | **detail** |
| --- | --- |
| expr_mean | the mean value of gene in different cell |
| regulons_matrix | the normalized enrichment value of tanscriptional factor in different cell |
| gsea.list | the enrichment result and target genes of tanscriptional factor in different cell |
| fc.list  | the fold change value between specific cell type and others |
| expr_r_regulons | the sum value of normalized enrichment value of tanscriptional factor downstreaming specific receptor |
| softmax_ligand | softmax value of the ligand expression across all cell types |
| softmax_receptor | softmax value of the receptor expression across all cell types |
| expr_l_r | the score of ligand-receptor in cellA-cellB |
| expr_l_r_log2 | the score of ligand-receptor in cellA-cellB with log transform |
| expr_l_r_log2_scale | the score of ligand-receptor in cellA-cellB with log transform and scale to [0,1] |
| DistanceKEGG | the distance between receptor and tf in one pathway in KEGG |



<a name="4131104a"></a>
### 2.2 visualization of the ligand-receptor-TF model


<a name="4d1f704f"></a>
#### 2.2.1 overveiw in circle plot

<br />set the color of each cell type<br />

```
  cell_color <- data.frame(color=c('#e31a1c','#1f78b4',
                                   '#e78ac3','#ff7f00'), stringsAsFactors = FALSE)
  rownames(cell_color) <- c("SSC", "SPGing", "SPGed", "ST")
```

<br />plot circle with cellcall object<br />

```
ViewInterCircos(object = mt, font = 2, cellColor = cell_color, lrColor = c("#F16B6F", "#84B1ED"),
                  arr.type = "big.arrow",arr.length = 0.04,
                  trackhight1 = 0.05, slot="expr_l_r_log2_scale",linkcolor.from.sender = TRUE,
                  linkcolor = NULL, gap.degree = 2,order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                  trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = FALSE)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168608941-00510457-8fa4-425e-bb73-57b1c8718f01.png#align=left&display=inline&height=475&margin=%5Bobject%20Object%5D&name=image.png&originHeight=475&originWidth=559&size=74933&status=done&style=none&width=559#id=Z16r5&originHeight=475&originWidth=559&originalType=binary&status=done&style=none)<br />plot circle with DIY dataframe of mt@data$expr_l_r_log2_scale<br />

```
  ViewInterCircos(object = mt@data$expr_l_r_log2_scale, font = 2, cellColor = cell_color,
                  lrColor = c("#F16B6F", "#84B1ED"),
                  arr.type = "big.arrow",arr.length = 0.04,
                  trackhight1 = 0.05, slot="expr_l_r_log2_scale",linkcolor.from.sender = TRUE,
                  linkcolor = NULL, gap.degree = 2,order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                  trackhight2 = 0.032, track.margin2 = c(0.01,0.12), DIY = T)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168681351-261eaa95-e761-4f08-9e2f-76462fa223b4.png#align=left&display=inline&height=452&margin=%5Bobject%20Object%5D&name=image.png&originHeight=452&originWidth=534&size=72472&status=done&style=none&width=534#id=tqJ8P&originHeight=452&originWidth=534&originalType=binary&status=done&style=none)<br />

<a name="c2fcc6f5"></a>
#### 2.2.2 overveiw in pheatmap plot


```
viewPheatmap(object = mt, slot="expr_l_r_log2_scale", show_rownames = T,show_colnames = T,
             treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 12,angle_col = "45",
             main="score")
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168729707-830a042f-afbd-4a67-a0fa-f3804abc6aac.png#align=left&display=inline&height=477&margin=%5Bobject%20Object%5D&name=image.png&originHeight=477&originWidth=1185&size=112516&status=done&style=none&width=1185#id=EJbnc&originHeight=477&originWidth=1185&originalType=binary&status=done&style=none)<br />

<a name="b101b966"></a>
#### 2.2.3 inspect Ligand-Receptor-TF in specific cellA-cellB

<br />There are three types to show this triple relation.<br />Funtion LR2TF to measure the triple relation between specific cells.<br />

```
  mt <- LR2TF(object = mt, sender_cell="ST", recevier_cell="SSC",
              slot="expr_l_r_log2_scale", org="Homo sapiens")
  head(mt@reductions$sankey)
```

<br />First type, function LRT.Dimplot.<br />

```
  if(!require(networkD3)){
  		BiocManager::install("networkD3")
  }
  
  sank <- LRT.Dimplot(mt, fontSize = 8, nodeWidth = 30, height = NULL, width = 1200, sinksRight=FALSE, DIY.color = FALSE)
  networkD3::saveNetwork(sank, "~/ST-SSC_full.html")
```

<br />The first pillar is ligand，the second pillar is receptor，the last pillar is tf.<br />And the color of left and right flow is consistent with ligand and receptor respectively.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608037011420-8a08e2c5-73f0-4d9a-bc49-87119980fcda.png#align=left&display=inline&height=287&margin=%5Bobject%20Object%5D&name=image.png&originHeight=714&originWidth=1855&size=282995&status=done&style=none&width=746#id=WOHOu&originHeight=714&originWidth=1855&originalType=binary&status=done&style=none)

---

Second type, function sankey_graph with isGrandSon = FALSE.<br />

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

<br />The first pillar is ligand，the second pillar is receptor，the last pillar is tf.<br />And the color of left and right flow is consistent with ligand and receptor respectively.<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1612168925812-ab05a489-d7f4-497f-b517-939486808a11.png#align=left&display=inline&height=488&margin=%5Bobject%20Object%5D&name=image.png&originHeight=488&originWidth=547&size=183644&status=done&style=none&width=547#id=sLu3g&originHeight=488&originWidth=547&originalType=binary&status=done&style=none)

---

Third type, function sankey_graph with isGrandSon = TRUE.<br />

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

<br />The first pillar is ligand，the second pillar is receptor，the last pillar is tf.<br />And the color of left and right flow is consistent with one node (ligand or receptor).<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608038695456-2ce50409-aaf1-4484-bfda-bde39ec7738f.png#align=left&display=inline&height=562&margin=%5Bobject%20Object%5D&name=image.png&originHeight=562&originWidth=1025&size=322055&status=done&style=none&width=1025#id=I8hDo&originHeight=562&originWidth=1025&originalType=binary&status=done&style=none)<br />

<a name="5e449248"></a>
### 2.3 enrichment in pathway

<br />getHyperPathway to perform enrichment, getForBubble to merge data for graph and plotBubble produce the bubble plot.<br />

```
n <- mt@data$expr_l_r_log2_scale

pathway.hyper.list <- lapply(colnames(n), function(i){
    print(i)
    tmp <- getHyperPathway(data = n, object = mt, cella_cellb = i, Org="Homo sapiens")
    return(tmp)
})

myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
p <- plotBubble(myPub.df)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608044718538-54a82161-c714-4389-80f9-e715429c0754.png#align=left&display=inline&height=685&margin=%5Bobject%20Object%5D&name=image.png&originHeight=685&originWidth=1047&size=161481&status=done&style=none&width=1047#id=YbZcY&originHeight=685&originWidth=1047&originalType=binary&status=done&style=none)<br />

<a name="ac88ba82"></a>
### 2.4 activation of TF in receiver cell


<a name="8431efc8"></a>
#### 2.4.1 ridge plot

<br />plot enrichment result of TF (filter or not).<br />

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

<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608099926349-53777dc6-352b-484e-a3a2-9a4c519906a1.png#align=left&display=inline&height=686&margin=%5Bobject%20Object%5D&name=image.png&originHeight=686&originWidth=688&size=117664&status=done&style=none&width=688#id=unTPP&originHeight=686&originWidth=688&originalType=binary&status=done&style=none)<br />

<a name="e4de6250"></a>
#### 2.4.2 gsea plot


```
  ssc.tf <- names(mt@data$gsea.list$SSC@geneSets)
  ssc.tf
```

<br />Show all TF have result in the SSC.<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608100011863-2c349ce5-f197-4528-a167-2eec6004524c.png#align=left&display=inline&height=122&margin=%5Bobject%20Object%5D&name=image.png&originHeight=122&originWidth=1317&size=19698&status=done&style=none&width=1317#id=Zcli0&originHeight=122&originWidth=1317&originalType=binary&status=done&style=none)<br />

```
  getGSEAplot(gsea.list=mt@data$gsea.list, geneSetID=c("CREBBP", "ESR1", "FOXO3"), myCelltype="SSC",
              fc.list=mt@data$fc.list,  selectedGeneID = mt@data$gsea.list$SSC@geneSets$CREBBP[1:10],
              mycol = NULL)
```

<br />![image.png](https://cdn.nlark.com/yuque/0/2020/png/1705105/1608100142219-a62cfef9-fe81-4a2e-8222-73d30770a122.png#align=left&display=inline&height=682&margin=%5Bobject%20Object%5D&name=image.png&originHeight=682&originWidth=692&size=75484&status=done&style=none&width=692#id=wOTtk&originHeight=682&originWidth=692&originalType=binary&status=done&style=none)<br />

<a name="e56ue"></a>
## 3. Structure of S4 object
<a name="ad42O"></a>
### 3.1 all data slots stored in the S4 object of R package
More detailed results in the intermediate process are stored in the S4 object of R package for some customized analyses. And explaination of every slot  is in the table below.

| **slot** | **detail** |
| --- | --- |
| count | raw input matrix input in the function CreateNichConObject |
| withoutlog | the matrix transformed from raw input, TPM or CPM |
| expr_mean | the mean value of gene in different cell |
| regulons_matrix | the normalized enrichment value of tanscriptional factor in different cell |
| gsea.list | the enrichment result and target genes of tanscriptional factor in different cell |
| fc.list  | the fold change value between specific cell type and others |
| expr_r_regulons | the sum value of normalized enrichment value of tanscriptional factor downstreaming specific receptor |
| softmax_ligand | softmax value of the ligand expression across all cell types |
| softmax_receptor | softmax value of the receptor expression across all cell types |
| expr_l_r | the score of ligand-receptor in cellA-cellB |
| expr_l_r_log2 | the score of ligand-receptor in cellA-cellB with log transform |
| expr_l_r_log2_scale | the score of ligand-receptor in cellA-cellB with log transform and scale to [0,1] |
| DistanceKEGG | the distance between receptor and tf in one pathway in KEGG |

<a name="X7o5u"></a>
### 3.2 present or export the details of the TG list
For some biologists, they pay more attention to the TGs of the TF. This tool provide options to present or export the details of the TG list, stored in the NichConObject@data$gsea.list$cell_type@geneSets.
```r
mt@data$gsea.list$SSC@geneSets
```
The figure presents a part of result stored in the NichConObject@data$gsea.list$cell_type@geneSets.<br />![image.png](https://cdn.nlark.com/yuque/0/2021/png/1705105/1622793918329-07ac2149-03ff-4243-ba2f-06345a112b69.png#clientId=uc7e6a175-7c8e-4&from=paste&height=235&id=ud5973156&margin=%5Bobject%20Object%5D&name=image.png&originHeight=469&originWidth=980&originalType=binary&size=114458&status=done&style=none&taskId=u62af5ea2-5c00-4ce1-8da3-2dac8a7047a&width=490)<br />​<br />
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
