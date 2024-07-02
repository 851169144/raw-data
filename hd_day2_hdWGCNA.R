# 推荐本地安装
#devtools::install_local('hdWGNCA.zip')
#BiocManager::install('harmony',update=F,ask=F)
library(hdWGCNA)
#加载单细胞分析包
library(Seurat)
#加载作图包
library(tidyverse)
library(cowplot)
library(patchwork)
#加载共表达网络分析包
#BiocManager::install('WGCNA',update=F,ask=F)
library(WGCNA)
gc()

setwd('./scRNA/')
#设置随机种子
set.seed(12345)

rm(list = ls())
## 读取上节课的数据
scRNA=readRDS('./scRNA_T.RDS')

#过滤出至少在5%的细胞中表达的基因
scRNA <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Bio_com" # the name of the hdWGCNA experiment
)


#构建metacells!!这一步非常重要，WGCNA对数据的稀疏性非常敏感，与普通转录组的WGCNA分析相比
# 单细胞的稀疏矩阵的解决方法是WGCNA流程的问题核心

# construct metacells  in each group
scRNA<- MetacellsByGroups(
  seurat_obj = scRNA,k=20,
  max_shared = 10,
  # group.by一般关注的是组织类型和细胞类型!
  group.by = c("tissue_type",'T_celltype'), # 也可以选择别的groupby
  ident.group = 'tissue_type' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
scRNA <- NormalizeMetacells(scRNA)
metacell_obj <- GetMetacellObject(scRNA)


#转置表达矩阵
# 安全起见，另起一个对
seurat_obj  <- SetDatExpr(
  scRNA,
  group_name = "CD4_EM", # 选择感兴趣恶的细胞类群！
  group.by='T_celltype' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)


#选择softpower
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # 这边不用转置了，前面已转置
)


# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
pdf(file="FIG 2 hdwgcna.pdf", width=10, height=8)
wrap_plots(plot_list, ncol=2)
dev.off()


#查看powerTable
power_table <- GetPowerTable(seurat_obj)
head(power_table)


#构建共表达网络
softpower=8  # 根据自己的改!!
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=softpower,
  setDatExpr = F,overwrite_tom = T)


dev.off()

#可视化WGCNA网络
pdf(file="FIG 2 hdwcgna_2.pdf", width=8, height=4)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

#(可选)获取TOM矩阵，可用来进行其他高级分析
TOM <- GetTOM(seurat_obj)


#计算模块协调特征
#记得scale一下 or else harmony throws an error:
seurat_obj <- Seurat::ScaleData(
  seurat_obj,
  features = GetWGCNAGenes(seurat_obj),
  
)
# 计算ME，根据组织类型分组
# harmony必须biocManager安装，不可以用github安装！！！
library(harmony)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="tissue_type" #harmony对象
)



seurat_obj <- ModuleConnectivity(seurat_obj)

# plot genes ranked by kME for each module
#可视化每个模块中，按照kME打分的基因PlotKMEs(seurat_obj, ncol=5)

# 获取hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)

head(hub_df)
pdf(file="FIG 2 hdwgcna3.pdf", width=8, height=4)
PlotKMEs(seurat_obj)
dev.off()
#记得保存上面hdWGNCA关键分析过程！
saveRDS(seurat_obj, file='hdWGCNA_object.rds')

dev.off()


####------------一些可视化-----------------------
## 模块间的相关性
library(igraph)
library(qgraph)
# 载入保存的

seurat_obj=readRDS('hdWGCNA_object.rds')

# 画模块间相关性图

pdf(file="FIG 2 hdwgcna4cor.pdf", width=6, height=6)
ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2)
dev.off()



## 0，2，4monocytes到底与哪些模块相关
# 由于识别到了每个模块的hub基因，可以去计算hub基因打分
# compute gene scoring for the top 25 hub genes by kME for each module
# (方法一)with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# (方法二)with UCell method #推荐这种方法
# 由于Ucell刚刚更新，所以4.1.x的同学请用本地安装,依赖包自行安装
#devtools::install_local("UCell-1.3.zip")
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# featureplot
# 瞄一眼
DimPlot(scRNA, label=TRUE,split.by = 'tissue_type',group.by = "T_celltype") 

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE ,# order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf(file="FIG 3 hdwgcna5.pdf", width=8, height=4)
wrap_plots(plot_list)
dev.off()

### dotplot
# get hMEs from seurat object

MEs <- GetMEs(seurat_obj, harmonized=TRUE)

mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods,group.by = "T_celltype")

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
pdf(file="FIG 2 hdwgcna6.pdf", width=6, height=8)
p
dev.off()

hub_df <- GetHubGenes(seurat_obj)
### !!!!!
gene=hub_df[hub_df$module %in% c('red','green'),]
gene=gene$gene_name
save(gene,file ='key_hdwgcna_gene.Rdata')
