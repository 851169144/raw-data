setwd('e:/work/hd/scRNA/')

### 补充做个拟时序，可以和前面的内容放在一起
scRNAsub=readRDS('scRNA_monocyte.RDS')


## 拟时序分析 先关掉再开    suggest 4.1.3
library(Seurat)
#没有monocle要先安装 BiocManager::install
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNAsub@assays$RNA@counts)
data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修改
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因，不修改
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)


#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
plot1

##Pseudotime轨迹图
plot2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
plot2


plot4 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
plot4

seurat_obj=readRDS('hdWGCNA_object.rds')
library(hdWGCNA)
library(tidyverse)
# 可以根据自己的需要扩大个数
hub_df <- GetHubGenes(seurat_obj)
### !!!!!
gene=hub_df[hub_df$module %in% c('red','green'),]
gene=gene$gene_name

my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[gene,],
                                                 # num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

my_pseudotime_cluster 



##### 准备富集和PPI 需要的gene list
library(hdWGCNA)
library(tidyverse)
seurat_obj=readRDS('hdWGCNA_object.rds')
# 可以根据自己的需要扩大个数
hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)
gene=hub_df[hub_df$module %in% c('red','blue','yellow'),]
gene=gene$gene_name

### 人不用运行！！！
library(Hmisc)
gene=toupper(gene)
gene=capitalize(gene)


write.csv(gene,file ='hdWGCNA_gene_for_PPI.csv',quote = F)


##############
FeaturePlot(scRNAsub,features = 'Arg1',split.by = 'tissue_type')
#https://metascape.org/gp/index.html
#http://genemania.org/


##### 药物预测，分子对接准备，提上来讲###########

#https://ctdbase.org/
# 药物-蛋白-疾病
dev.off()

##在卒中中，ARG1的抑制剂

setwd('e:/work/hd')
chemi=read.csv(file ='CTD_383_diseases_20230101092921.csv',check.names = F)
chemi=chemi[chemi$`Disease Name`=='Stroke',]
chemi=chemi$`Inference Network`
chemi=stringr::str_split(chemi,pattern = '\\|')
chemi=chemi[[1]]

intera=read.csv(file ='CTD_383_ixns_20230101093249.csv',check.names = F)

intera=intera[intera$`Chemical Name` %in% chemi,]


### https://www.rcsb.org/
### https://pubchem.ncbi.nlm.nih.gov/
## https://www.dockeasy.cn/DockCompound


# parameter_file AD4_parameters.dat
