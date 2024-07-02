load('merge.normalzie.Rdata')
rt=outTab[,group_list=='ASC']

cluster=read.table('./cluster.txt',row.names = 1)
colnames(cluster)='Cluster'

cluster$Cluster=paste0('Cluster',cluster$Cluster)

rt=rt[,rownames(cluster)]

cluster$Cluster1=ifelse(cluster$Cluster=='Cluster1',1,0)
cluster$Cluster2=ifelse(cluster$Cluster=='Cluster2',1,0)

cluster$Cluster=NULL

gc()
#BiocManager::install("impute")
library(WGCNA)
# 可换5000
WGCNA_matrix = t(rt[order(apply(rt,1,mad), decreasing = T)[1:2000],])
datExpr0 <- WGCNA_matrix  ## top mad genes
datExpr0 <- as.data.frame(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # 返回TRUE则继续
# （可选）如果存在太多的缺失值
if (!gsg$allOK){
  # 把含有缺失值的基因或样本打印出来
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # 去掉那些缺失值
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

## 样本过滤前
sampleTree = hclust(dist(datExpr0), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

## 根据图片挑选cutheight,我们不丢样本了!!!!!
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust) # 0代表切除的，1代表保留的
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]

## 更新anno
cluster=cluster[rownames(datExpr),,drop=F]

## 样本过滤后
sampleTree = hclust(dist(datExpr), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers (after)", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


datTraits =cluster
datExpr=datExpr[rownames(datTraits),]
sampleNames = rownames(datExpr)
# 能全部对上
traitRows = match(sampleNames, rownames(datTraits))  

###power值散点图
#allowWGCNAThreads()
enableWGCNAThreads()   #多线程工作
powers = c(1:30)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,blockSize = 100000)

dev.off()
par(mfrow = c(1,2))
cex1 = 0.85
###拟合指数与power值散点图

P1<-plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red") #可以修改
###平均连通性与power值散点图
P2<-plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.1,col="red") #可以修改
pdf(file="FIG  WCGNA 1.pdf", width=8, height=4)
P1+P2
dev.off()
###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate #最佳power值
# 发现不合适就自定义，此处我自定义!!!!!
#softPower=6
adjacency = adjacency(datExpr, power = softPower)

net = blockwiseModules(datExpr, power = softPower,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PRTOM",
                       verbose = 3)
# 显示模块数量以及各自包含的基因数目
# 0表示未分入任何模块的基因
# 1是最大的模块，往后依次降序排列，分别对应各自模块的基因
table(net$colors)


mergedColors = labels2colors(net$colors)
mergedColors
##手动保存
pdf(file="FIG WGCNA 2.pdf", width=8, height=4)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 用color labels重新计算MEs（Module Eigengenes:模块的第一主成分）
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") #（这是重点）计算ME和表型相关性
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


# 设置热图上的文字（两行数字：第一行是模块与各种表型的相关系数；
# 第二行是p值）
# signif 取有效数字
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# 然后对moduleTraitCor画热图
pdf(file="FIG WCGNA 3.pdf", width=8, height=12)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()
##把turquoise挑出来
# 选择导出模块
module = c('yellow')
# 选择模块中基因/探针
probes = names(datExpr)
inModule = (moduleColors %in% module)
modProbes_POS = probes[inModule]
#modprobes可以后续分析
modProbes_POS
write.table(modProbes_POS,file ='MEyellow.txt',row.names = F,col.names = F,quote=F)

#GO和KEGG----
library("org.Hs.eg.db")
rt=read.table("MEyellow.txt",sep="\t",check.names=F,header=F)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,1,entrezID=entrezIDs)
colnames(out)=c("symbol","logFC","entrezID")
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)

#GO----
dir.create("WCGNNA_GOKEGG")
setwd("./WCGNNA_GOKEGG")
library("clusterProfiler")
library("org.Hs.eg.db")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

cor=rt$logFC
gene=rt$entrezID
names(cor)=gene

#GO????????
kk <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)

#??״ͼ
tiff(file="GO_barplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory =20)
dev.off()

#??ͼ
tiff(file="GO_dotplot.tiff",width = 35,height = 20,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk,showCategory = 20)
dev.off()
#KEGG----
library("clusterProfiler")
library("KEGG.db")
rt=read.table("id.txt",sep="\t",header=T,check.names=F)
rt=rt[is.na(rt[,"entrezID"])==F,]

geneFC=rt$logFC
gene=rt$entrezID
names(geneFC)=gene

#kegg????????
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05,use_internal_data = T)
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)

#??״ͼ
tiff(file="KEGG_barplot.tiff",width = 30,height = 20,units ="cm",compression="lzw",bg="white",res=600)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()

#??ͼ
tiff(file="KEGG_dotplot.tiff",width = 30,height = 20,units ="cm",compression="lzw",bg="white",res=600)
dotplot(kk, showCategory = 20)
dev.off()
#BiocManager::install("pathview")
#ͨ·ͼ
library("pathview")
keggxls=read.table("KEGG.txt",sep="\t",header=T)
for(i in keggxls$ID){
  pv.out <- pathview(gene.data = geneFC, pathway.id = i, species = "hsa", out.suffix = "pathview")
}




########FGSEA----
## 导入基因集
#BiocManager::install('fgsea')
library(fgsea)
library(msigdbr)
msigdbr_species()
a=msigdbr_collections()
# 假设做鼠，人就Homo sapiens
m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = 'BP')
# 变list
BP <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

head(BP)

# 构造预制函数，内置流程是先排序差异基因，针对GO_BP的fgsea
preranked_BP <- function(x) {
  ranks <- x %>% 
    na.omit()%>%
    mutate(ranking=logFC)
  ranks <- ranks$ranking
  names(ranks) <- rownames(x)
  head(ranks, 10)
  set.seed(123456)
  BP_x <- fgsea(pathways = BP, 
                stats = ranks,
                minSize=10,
                maxSize=500,
                nperm=1000)
  
  BP_x$pathway<-gsub("GOBP_","",BP_x$pathway)
  BP_x$pathway<-gsub("_"," ",BP_x$pathway)
  return(BP_x)
}

# 获得各亚型特征的通路
## cluster1的marker基因和fgsea
# 再次确认
identical(colnames(rt),rownames(cluster))

library(limma)
# 先分组信息使用cluster1列

group_list=cluster$Cluster1
group_list=ifelse(group_list==0,'other','Cluster1')
group_list=factor(group_list,levels=c('other','Cluster1'))
design=model.matrix(~ group_list)

fit=lmFit(rt,design)
fit=eBayes(fit) 
c1_allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1) 

# 使用预置函数直接进行fgsea
library(dplyr)
BP_c1 <- preranked_BP(c1_allDiff)
sig_BP_c1 <- BP_c1 %>% filter(abs(NES)>1 & pval<0.01)
sig_BP_c1 <- sig_BP_c1[order(sig_BP_c1$NES,decreasing = T),]

## cluster2类的marker基因和fgsea
# 先分组信息使用cluster2列

group_list=cluster$Cluster2
group_list=ifelse(group_list==0,'other','Cluster2')
group_list=factor(group_list,levels=c('other','Cluster2'))
design=model.matrix(~ group_list)


fit=lmFit(rt,design)
fit=eBayes(fit) 
c2_allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1)



# 使用预置函数直接进行fgsea
library(dplyr)
BP_c2 <- preranked_BP(c2_allDiff)
sig_BP_c2 <- BP_c2 %>% filter(abs(NES)>1 & pval<0.01)
sig_BP_c2 <- sig_BP_c2[order(sig_BP_c2$NES,decreasing = T),]



## cluster3类的marker基因和fgsea------
# 先分组信息使用cluster2列

group_list=cluster$Cluster3
group_list=ifelse(group_list==0,'other','Cluster3')
group_list=factor(group_list,levels=c('other','Cluster3'))
design=model.matrix(~ group_list)


fit=lmFit(rt,design)
fit=eBayes(fit) 
c3_allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=1)



# 使用预置函数直接进行fgsea
library(dplyr)
BP_c3 <- preranked_BP(c3_allDiff)
sig_BP_c3 <- BP_c3 %>% filter(abs(NES)>1 & pval<0.01)
sig_BP_c3 <- sig_BP_c3[order(sig_BP_c3$NES,decreasing = T),]




# 合并作图-----
##A和B合并 suffixes后缀，一个是A，一个是B

merged <- merge(sig_BP_c1[sample(1:nrow(sig_BP_c1),20),c(1,5)], sig_BP_c2[sample(1:nrow(sig_BP_c2),20),c(1,5)], by = "pathway" , all = T,
                suffixes = c(".c1",".c2"))

## 重复合并，合并C
#merged2 <- merge(merged,sig_BP_c3[sample(1:nrow(sig_BP_c3),20),c(1,5)], by = "pathway" , all = T,
              #   suffixes = c(".c1",".mixed"))

merged3=merged
colnames(merged3) <- c("pathway","Cluster1", "Cluster2")
merged3[is.na(merged3)] <- 0
merged3=as.data.frame(merged3)
rownames(merged3) <- merged3$pathway
rownames(merged3)<- tolower(rownames(merged3))
merged3[,1] <- NULL

library(pheatmap)
pdf(file="clus_ph.pdf", width=12, height=8)
pheatmap(merged3, cluster_cols = F, cluster_rows = T, border_color=NA, 
         cellwidth =30,  color = colorRampPalette(c(rep("Darkblue",1), "white", rep("red",1)))(1000) ,
         breaks = seq(-2,2,length.out = 1000), main = "Biological process enrichment", 
         angle_col = 45, fontsize=14, fontsize_row = 10)
dev.off()
save(c1_allDiff,c2_allDiff,file ='alldiff_cluster.Rdata')



####### 输出protemap
load('alldiff_cluster.Rdata')
c2_sig=c2_allDiff[c2_allDiff$logFC>0.5 & c2_allDiff$adj.P.Val<0.05,]
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
gene = bitr(rownames(c2_sig), fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(SYMBOL = rownames(c2_sig),logFC=c2_sig$logFC)

gene_df <- merge(gene,gene_df,by="SYMBOL")
gene_df <- dplyr::distinct(gene_df,UNIPROT,.keep_all=TRUE)

prote=gene_df[,-1]
write.table(prote,file ='prote_c2.tsv',sep = '\t',row.names = F,col.names = F,quote = F)


c1_sig=c1_allDiff[c1_allDiff$logFC>0.5 & c1_allDiff$adj.P.Val<0.05,]
library(clusterProfiler)
library(org.Hs.eg.db)
library(clusterProfiler)
gene = bitr(rownames(c1_sig), fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(SYMBOL = rownames(c1_sig),logFC=c1_sig$logFC)

gene_df <- merge(gene,gene_df,by="SYMBOL")
gene_df <- dplyr::distinct(gene_df,UNIPROT,.keep_all=TRUE)

prote=gene_df[,-1]
write.table(prote,file ='prote_c1.tsv',sep = '\t',row.names = F,col.names = F,quote = F)
#https://www.proteomaps.net/