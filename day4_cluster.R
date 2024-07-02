load('merge.normalzie.Rdata')
rt=outTab[,group_list=='ASC']
load('gene_rfe_rf.Rdata')

rt=rt[rownames(rt) %in% gene,]

#BiocManager::install("ConsensusClusterPlus")
#聚类
library(ConsensusClusterPlus)
rt=as.matrix(rt)
maxK=10
results=ConsensusClusterPlus(rt,
                             maxK=maxK,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title='final_cluster',
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123456,
                             plot="pdf",)

#一致性打分
dir.create('calscore')
setwd('calscore/')
calcICL(results, title="consensusScore", plot="pdf")
#输出结果
clusterNum=2      #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
#setwd('./stroke_mito/')
setwd('..')
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)



### 作热图
load('merge.normalzie.Rdata')
rt=outTab[,group_list=='ASC']
load('gene_rfe_rf.Rdata')

rt=rt[rownames(rt) %in% gene,]

cluster=read.table('cluster.txt',row.names = 1)
cluster=cluster[order(cluster$V2),,drop=F]
colnames(cluster)='Cluster'
rt=rt[,rownames(cluster)]
cluster$Cluster=factor(cluster$Cluster)
pheatmap::pheatmap(rt,scale = 'row',cluster_cols = F,border_color = NA,
                   breaks = seq(-2,2,length.out = 100),annotation_col = cluster,
                   show_colnames = F)


####炎症因子
load('merge.normalzie.Rdata')
rt=outTab[,group_list=='ASC']
gene=read.table('BIOCARTA_INFLAM_PATHWAY.v7.5.1.gmt')[3:31]
gene=as.data.frame(t(gene))
gene=gene$V1
rt1=rt[rownames(rt)%in% gene,]
#identical(colnames(rt1),rownames(anno))


cluster=read.table('cluster.txt')
colnames(cluster)=c('Sample','Cluster')


rt1=as.data.frame(t(rt1))
rt1=rt1[cluster$Sample,]

rt1$group=cluster$Cluster


library(tidyverse)
rt2=tidyr::pivot_longer(rt1,cols = -c('group'),names_to = "cell",values_to = 'expression')

rt2$group=factor(rt2$group)
library(ggpubr)
pdf(file="FIG 3 cluster_immf.pdf", width=8, height=4)
ggboxplot(
  rt2,
  x = "cell",
  y = "expression",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "expression", palette=c("#20854e","#E18727",'#bc3c29')
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
