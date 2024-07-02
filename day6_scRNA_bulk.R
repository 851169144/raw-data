### 下游分析（轨迹、细胞通讯、代谢、bulk）

## 回到单细胞看表达-----------------
scRNA=readRDS('./scRNA_anno.RDS')
load('table1.Rdata')

## 关键基因
gene=unique(table1$Gene)



library(Seurat)
library(viridis)
DotPlot(scRNA,features = gene,cols = c('#dadada','#bc3c29'))


scRNA_T=readRDS('./scRNA_T.RDS')

library(viridis)
FeaturePlot(scRNA_T,features = 'APOBEC3G',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))
FeaturePlot(scRNA_T,features = 'YWHAQ',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'))


##  轨迹的基因表达分析----------------------------------

library(mixtools)
#devtools::install_github("SGDDNB/GeneSwitches")
library(GeneSwitches)
#install.packages('fastglm')
library(SingleCellExperiment)



## time and expression in CD8_CM
scRNA_cm=subset(scRNA_T,T_celltype=='CD8_CM')

cellinfo <- scRNA_cm@meta.data

## 构建SingleCellExperiment对象
sce <- as.SingleCellExperiment(scRNA_cm)

## run
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP', 
                           start.clus = c(3,5), shrink = 0.2)





dev.off()
## 可视化
cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)

## 下面这行关键，否则容易报错！！（特别是Matrx>1.5-0的同学）
igraph::igraph.options(sparsematrices = FALSE)

## 曲线折线仍选一个即可
#lines(SlingshotDataSet(sce_slingshot), lwd=2,col = 'black')#,type = 'lineages'
lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = 'black')


legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)

### 开关基因（驱动基因分析）
allexpdata <- as.matrix(scRNA_cm@assays$RNA@data);dim(allexpdata)
allcells<-colData(sce_slingshot);dim(allcells)

allcells$slingPseudotime_1

cells <- allcells[!is.na(allcells$slingPseudotime_1),];dim(cells)
expdata <- allexpdata[,rownames(cells)];dim(expdata)

#filter genes expressed in less than 5 cells
#过滤少于五个细胞表达的基因
expdata <- expdata[apply(expdata > 0,1,sum) >= 5,];dim(expdata)

rd_UMAP <- Embeddings(object = scRNA_cm, reduction = "umap");dim(rd_UMAP)#原 object = seu3obj.integrated
rd_UMAP <- rd_UMAP[rownames(cells), ];dim(rd_UMAP)
all(rownames(rd_UMAP) == colnames(expdata))

library(mixtools)
library(GeneSwitches)
library(SingleCellExperiment)


## create SingleCellExperiment object with log-normalized single cell data
## 使用对数规范化的单细胞数据创建SingleCellExperiment对象
sce <- SingleCellExperiment(assays = List(expdata = expdata))
## 添加伪时间信息
colData(sce)$Pseudotime <- cells$slingPseudotime_1
## 添加降维，例如 PCA、UMAP、tSNE
reducedDims(sce) <- SimpleList(UMAP = rd_UMAP)
sce_p1 <- sce
### 检查二值化阈值
h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.2, col="blue")}

## 二值化分析
sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = 0.2)
# sce_p1 <- binarize_exp(sce_p1, ncores = 3)
sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE, zero_ratio = 0.65, ds_cutoff = 0.65)

table(rowData(sce_p1)$prd_quality)

# 过滤出开关基因
sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92);dim(sg_allgenes)


sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92,
                                genelists = gs_genelists);dim(sg_gtypes)#, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),]);dim(sg_vis)

## 自己关注的基因
gl <- unique(table1$Gene)
intersect(sg_vis$geneID, gl)
sg_my <- rowData(sce_p1)[gl,];head(sg_my)
sg_my$feature_type <- "Mendelian genes"
sg_vis <- rbind(sg_vis, sg_my)
plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3.5)

## R2大于0，上调型开关基因，R2小于0 ，下调型开关基因
a=sce_p1@assays@data$expdata['APOBEC3G',] ## !!!
b=sce_p1$Pseudotime

df=data.frame(gene=a,time=b)

ggstatsplot::ggscatterstats(data=df,x='time',y='gene')



####细胞通讯--------阴阳性群
gc()

scRNA=readRDS('./scRNA_anno.RDS')
scRNA_other=subset(scRNA,celltype != 'T_cells')
rm(scRNA)
gc()
scRNA_T=readRDS('./scRNA_T.RDS')


gc()

##!!!
scRNA_CM=subset(scRNA_T,T_celltype=='CD8_CM')
scRNA_CM$gene_group=ifelse(scRNA_CM@assays$RNA@counts['APOBEC3G',]>0,'APOBEC3G+CD8_CM','APOBEC3G-CD8_CM')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD8_CM')

# 加列
scRNA_other$gene_group =scRNA_other$celltype
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_chat=merge(scRNA_CM,c(scRNA_other,scRNA_otherT))

rm(scRNA_CM,scRNA_otherT)
rm(scRNA_T)
rm(scRNA_other)
gc()

##
scRNA_chat_SLE=subset(scRNA_chat,tissue_type=='SLE')

set.seed(123)
a=sample(1:ncol(scRNA_chat_SLE),2000)
scRNA_chat_SLE=scRNA_chat_SLE[,a]

meta =scRNA_chat_SLE@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_SLE@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "gene_group")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, sources.use = c('APOBEC3G+CD8_CM','APOBEC3G-CD8_CM'),
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('APOBEC3G+CD8_CM','APOBEC3G-CD8_CM'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble


### 代谢相关------------------

scRNA_T=readRDS('./scRNA_T.RDS')

gc()
scRNA_CM=subset(scRNA_T,T_celltype=='CD8_CM')
scRNA_CM$gene_group=ifelse(scRNA_CM@assays$RNA@counts['APOBEC3G',]>0,'APOBEC3G+CD8_CM','APOBEC3G-CD8_CM')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD8_CM')

# 加列
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_metab=merge(scRNA_CM,c(scRNA_otherT))

gc()

rm(scRNA_chat)
gc()
scRNA_metab_SLE=subset(scRNA_metab,tissue_type=='SLE')

set.seed(123)
a=sample(1:ncol(scRNA_metab_SLE),2000)
scRNA_metab_SLE=scRNA_metab_SLE[,a]




#### scMetabolism评估巨噬细胞代谢活性
#devtools::install_github("YosefLab/VISION")
#devtools::install_github("wu-yc/scMetabolism")

library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_metab_SLE<-sc.metabolism.Seurat(obj = scRNA_metab_SLE, method = 'AUCell', imputation = F, ncores = 2, metabolism.type = "KEGG")

input.pathway <- rownames(scRNA_metab_SLE@assays[["METABOLISM"]][["score"]])[61:90]
DotPlot.metabolism(obj =scRNA_metab_SLE,
                   pathway = input.pathway, phenotype = "gene_group", norm = "y")


gc()
##差异基因------------------------------------
library(Seurat)
## Idents()
Idents(scRNA_CM)=scRNA_CM$gene_group
df=FindAllMarkers(scRNA_CM,only.pos = T,logfc.threshold =0.25)
write.csv(df,'APOBEC_marker.csv',quote = F)

## 富集分析怎么做,负值csv中的基因列到网站，参见下面教程
#https://mp.weixin.qq.com/s/ClHOFvw3GSM9wvmIPip4VA


## bulk---------------------------------------------
gc()

library(data.table)
rt=fread('./GSE112087_counts-matrix-EnsembIDs-GRCh37.p10.txt',data.table = F)
rownames(rt)=rt$V1
rt$V1=NULL


#IOBR有一系列依赖包
#https://mp.weixin.qq.com/s/nVziQeInS-4QxVNPQCilVQ
#   devtools::install_github("IOBR/IOBR")
library(IOBR)
#
#gc()
#rm(scRNA_chat_SLE)
#rm(scRNA_CM)
#rm(scRNA_T)
gc()

rt=count2tpm(rt,idType = 'Ensembl',org = 'hsa')

rt=as.data.frame(rt)

max(rt)
rt=log2(rt+1)

a1=grep('SLE',colnames(rt))

exp1=rt[,a1]
exp2=rt[,-a1]

rt=cbind(exp2,exp1)


load('table1.Rdata')

data=rt[unique(table1$Gene),]

### 注释文件
anno=data.frame(row.names =colnames(rt),group=c(rep('Healthy',58),
                                               rep('SLE',62)))
pheatmap::pheatmap(data,cluster_cols = F,
                   scale = 'row',show_colnames = F,annotation_col = anno)

df=data.frame(gene=as.numeric(rt['APOBEC3G',]),group=anno$group)
ggpubr::ggboxplot(data = df, x = 'group',y='gene',color = 'group',palette = 'jco',notch = T,size = 1)+
  stat_compare_means()

