# 关键细胞亚群的深入分析

#BiocManager::install("slingshot")
## 读取数据
scRNA_T=readRDS('./scRNA_T.RDS')
library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

## 细胞轨迹分析----------------------------

## 载入示例数据
table(scRNA_T$T_celltype)

cellinfo <- scRNA_T@meta.data


## 构建SingleCellExperiment对象
sce <- as.SingleCellExperiment(scRNA_T)

## run
sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP', 
                           start.clus = c(3,5), shrink = 0.2)



#lin1 <- getLineages(sce_slingshot, 
#                    clusterLabels = "seurat_clusters", 
#                   start.clus ="0",#可指定起始细胞簇，用处不大
#                    end.clus="5",#可指定终点细胞簇,用处不大
#                    reducedDim = "UMAP")



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
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(0.5,2,4), pch = 16)





# 细胞通讯----------------------
# 在SLE和衰老中分别做，看看功能是不是一样的

## cell-cell chat----
scRNA=readRDS('./scRNA_anno.RDS')
scRNA_other=subset(scRNA,celltype != 'T_cells')
rm(scRNA)

# 加一列！！！
scRNA_other$T_celltype =scRNA_other$celltype
scRNA_chat=merge(scRNA_other,scRNA_T)
rm(scRNA_T)
rm(scRNA_other)
gc()

scRNA_chat_old=subset(scRNA_chat,tissue_type=='OLD')


## 抽2000细胞做
set.seed(123)
a=sample(1:ncol(scRNA_chat_old),2000)
scRNA_chat_old=scRNA_chat_old[,a]

meta =scRNA_chat_old@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_old@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "T_celltype")

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
                 weight.scale = T, label.edge= T, sources.use = 'CD4_EM',
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                            sources.use = 'CD4_EM',
                           remove.isolate = FALSE)+coord_flip()
p_bubble

library(CellChat)
cellchat <- updateCellChat(cellchat)

## SLE
scRNA_chat_SLE=subset(scRNA_chat,tissue_type=='ASC')

set.seed(123)
#a=sample(1:ncol(scRNA_chat_SLE),2000,replace = T)
scRNA_chat_SLE=scRNA_chat_SLE[,]

meta =scRNA_chat_SLE@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_SLE@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "T_celltype")

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
                 weight.scale = T, label.edge= T, sources.use = 'CD4_EM',
                 title.name = "Number of interactions")
dev.off()

pp_bubble= netVisual_bubble(cellchat,
                           sources.use = 'CD4_EM',
                           remove.isolate = FALSE)+coord_flip()
pp_bubble

