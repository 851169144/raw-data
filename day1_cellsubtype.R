# 发现关键细胞亚群

library(Seurat)
gc()
# 读取数据
SLE1=Read10X('./ASC/')

CT1=Read10X('./CT1/')
CT3=Read10X('./CT3/')
CT3=CT3[[1]]
OLD1=Read10X('./OLD1/')
OLD1=OLD1[[1]]
OLD2=Read10X('./OLD2/')
OLD2=OLD2[[1]]

rm(scRNA)


## 建立seurat对象
SLE1=CreateSeuratObject(SLE1,project = 'ASC')
#SLE2=CreateSeuratObject(SLE2,project = 'SLE2')
CT1=CreateSeuratObject(CT1,project = 'CT1')
CT3=CreateSeuratObject(CT3,project = 'CT3')
OLD1=CreateSeuratObject(OLD1,project = 'OLD1')
OLD2=CreateSeuratObject(OLD2,project = 'OLD2')

### 合并scRNA1,scRNA2,scRNA3
## 没有正常样本则不加scRNA1
### 更多样本可以自行加入
#scRNA <- merge(scRNA1, y = c(scRNA2,scRNA3))
scRNA <- merge(SLE1,y=c(CT1,CT3,OLD1,OLD2))

# metadata为样本信息，我们需要定义分组
scRNA@meta.data$tissue_type=scRNA@meta.data$orig.ident
#install.packages('stringr')
# 去除数字,获得干净的分组
scRNA@meta.data$tissue_type=stringr::str_remove(scRNA@meta.data$tissue_type,'[0-9]')

# 质控
##计算质控指标
#计算细胞中线粒体核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
#head(scRNA@meta.data)
library(ggplot2)
col.num <- length(levels(scRNA@active.ident))
# 过滤前
VlnPlot(scRNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 


##设置质控标准，比较随意
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200
maxGene=4000
pctMT=10

##数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
col.num <- length(levels(scRNA@active.ident))
VlnPlot(scRNA,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
        cols =rainbow(col.num), 
        pt.size = 0.1, 
        ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#降维聚类#########################
library(Seurat)
library(tidyverse)
library(patchwork)


scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")

### harmony去批次
#BiocManager::
library(harmony)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNA_harmony,reduction = 'harmony')

# 一定要指定“harmony”！
scRNA <- FindNeighbors(scRNA_harmony, dims = 1:10, reduction = "harmony")
scRNA <- FindClusters(scRNA)
scRNA <- RunUMAP(scRNA, dims = 1:10,reduction = 'harmony')

# 去批次成功
DimPlot(scRNA,split.by = 'tissue_type')

rm(scRNA_harmony)
#BiocManager
library(SingleR)
# 人用下面
#refdata <- SingleR::HumanPrimaryCellAtlasData()
refdata <-get(load("refdata_singleR.Rdata"))
# 鼠用下面
#refdata <- SingleR::MouseRNAseqData()

library(Seurat)
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)


scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}


Idents(scRNA)=scRNA$celltype

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

DimPlot(scRNA,split.by = 'tissue_type',cols = col_vector[20:30])

## 保存数据
saveRDS(scRNA,file = 'scRNA_anno.RDS')


### 根据原文 我们直接提取T亚群
scRNA_T=subset(scRNA,celltype=='T_cells')
library(Seurat)

# 提T细胞亚群,重新降维聚类
scRNAsub<- FindVariableFeatures(scRNA_T, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub<- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNAsub)


## 重新harmony
library(harmony)
set.seed(1000)
scRNAsub <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")
DimPlot(scRNAsub, reduction = "harmony", group.by = "orig.ident")

ElbowPlot(scRNAsub,reduction = 'harmony')


scRNAsub <- FindNeighbors(scRNAsub, reduction = 'harmony',dims = 1:10)
scRNAsub <- FindClusters(scRNAsub)
scRNAsub <- RunUMAP(scRNAsub, reduction = 'harmony',dims = 1:10)

DimPlot(scRNAsub, label=TRUE,split.by = 'tissue_type') 

dev.off()

## 细胞比例图
library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$seurat_clusters,scRNAsub@meta.data$tissue_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:10] 

# 作图
prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop

# T细胞亚群注释参考
# https://blog.csdn.net/weixin_52505487/article/details/126687526
Idents(scRNAsub)=scRNAsub$seurat_clusters
T_marker <- c("CCR7","LEF1", "TCF7",'SELL','KLF2', #CD4_Naive
                "ANXA1", "CXCR4", "IL2", #CD4_EM
                "BCL6", "CXCR5","ICA1", #CD4_FH
                'IL23R',"CCR6",'CAPG','RORC','IL17A', #TH17
                'FOXP3','IL2RA','IL1R2',#CD4_REG
              'CD8A','GZMK',  # CD8_EM
              'GZMA','CCL5',  #CD8_CM
               'HAVCR2','PDCD1','LAG3', # CD8_exhau
                'EPCAM','CD19','CD3E') 
DotPlot(scRNAsub, 
        features = T_marker,
        group.by = "seurat_clusters") + coord_flip()


## 不要求非常准确
#'CD4_Naive','CD8_CM','CD8_EM',
#'CD4_Naive','Th17','CD4_Naive',
#'CD4_CM','CD4_Naive','CD8_EM',
#'CD4_REG','CD8_Naive','CD8_Exhau',
#'CD8_Exhau'

T_celltype=c('CD4_Naive',
             'CD4_Naive',
             'CD4_EM',
             'CD4_CM',
             'CD8_CM',
             'CD8_EM',
             'CD8_Naive',
             'CD4_Naive',
             'CD8_EM',
             'CD4_EM',
             'CD4_EM',
             'CD4_Naive',
             'CD4_Tfr')



Idents(scRNAsub) <- scRNAsub@meta.data$seurat_clusters
names(T_celltype) <- levels(scRNAsub)
scRNAsub<- RenameIdents(scRNAsub, T_celltype)

scRNAsub@meta.data$T_celltype <- Idents(scRNAsub)


#设置idents主要识别标签
Idents(scRNAsub)=scRNAsub@meta.data$T_celltype

colors=c('#313c63','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#bc3c29')
DimPlot(scRNAsub, group.by="T_celltype", label=T, label.size=5,cols = col_vector[30:50],
        pt.size = 1,split.by = 'tissue_type')

## 细胞比例图再次

library(reshape2)
library(ggplot2)
prop_df <- table(scRNAsub@meta.data$T_celltype,scRNAsub@meta.data$tissue_type) %>% melt()
colnames(prop_df) <- c("Cluster","Sample","Number")
prop_df$Cluster <- factor(prop_df$Cluster)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#处理后有73种差异还比较明显的颜色，基本够用
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

sample_color <- col_vector[1:10] 

prop <- ggplot(data = prop_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 
prop



saveRDS(scRNAsub,file = 'scRNA_T.RDS')
