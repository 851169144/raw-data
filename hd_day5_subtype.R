## 确保自己的路径
library(NMF)

seurat_obj=readRDS('./scRNA/hdWGCNA_object.rds')
library(hdWGCNA)
library(tidyverse)
# 可以根据自己的需要扩大个数
hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)
gene=hub_df[hub_df$module %in% c('red','blue','yellow'),]
gene=gene$gene_name
gene=read.table('alldiff_sig.csv', header=T, sep=",", check.names=T, row.names=1)
gene=rownames(gene)

rm(list=ls())
gc()
## 人并不用运行！！！！
library(Hmisc)
gene=toupper(gene)
gene=capitalize(gene)

load('merge.normalzie.Rdata')
pdata2=outTab
group_list2=group_list
exprSet=outTab
### gsva 评估Bulk 层面 是否有差异
load('lasso_genes.Rdata')
gene=ss
gene_set=list(gene)
names(gene_set)='MachineLearning_GENE'


library(genefilter)
library(GSVA)
library(Biobase)

# gsva方法
gsva_matrix<- gsva(as.matrix(exprSet), gene_set,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(t(gsva_matrix))
gsva_matrix$group=group_list2

library(ggpubr)
library(ggsci)
ggboxplot(gsva_matrix,x='group',y='MachineLearning_GENE',color = 'group',
          add='jitter')+ggsci::scale_color_jco()+stat_compare_means()



### lasso基因与免疫细胞相关性####################
# https://mp.weixin.qq.com/s/qpsRzkHXJlf-Hr7vXhoZPA
## devtools::install_github("IOBR/IOBR")
library(IOBR)

#！！！！ group_list
exprSet2=exprSet[,group_list2=='ASC']
mcp=IOBR::deconvo_mcpcounter(eset = exprSet2)
rownames(mcp)=mcp$ID
mcp$ID=NULL

####lasso gene读取
load('lasso_genes.Rdata')

rt=exprSet2[rownames(exprSet2) %in% genes,]
rt=as.data.frame(t(rt))
# 确保一样
identical(rownames(rt),rownames(mcp))


library(psych)
pp <- corr.test(rt,mcp,method="spearman",adjust = "fdr")

cor <- pp$r
pvalue <- pp$p

myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}

library(reshape2)
heatmap <- melt(cor)
colnames(heatmap)=c('sample','gene','cor')  

heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
  mutate(signif = sapply(pvalue, function(x) myfun(x)))


ggplot(heatmap,aes(sample,gene,col=cor))+
  geom_tile(color="grey70",fill="white",size=1)+
  geom_point(aes(size = abs(cor)),shape=15) + 
  geom_text(aes(label=signif),size=6,color="white",
            hjust=0.5,vjust=0.7)+
  labs(x = NULL,y = NULL,color=NULL) + 
  scale_color_gradient2(low = 'blue',high = 'red')+
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(),
        axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
        axis.ticks.y=element_blank(),
        panel.border = element_rect(fill=NA,color="grey70",
                                    size=2, linetype="solid")) +
  scale_size(range=c(1,10),guide=NULL)+
  guides(color = guide_colorbar(direction = "vertical",
                                reverse = F,barwidth = unit(.5, "cm"),
                                barheight = unit(15, "cm")))



################   NMF分型（疾病样本） ###################
library(NMF)

## !!!
exprSet2=exprSet[,group_list2=='ASC']

data=exprSet2[rownames(exprSet2) %in% genes,]

ranks=1:10
estim.coad <- nmf(data,rank= ranks, nrun=100)
plot(estim.coad)
seed = 2022
nmf.rank2 <- nmf(data, 
                 rank =3, 
                 nrun=100,
                 seed = seed, 
                 method = "brunet")
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
index <- extractFeatures(nmf.rank2,"max") 
sig.order <- unlist(index)
NMF.Exp.rank2 <- data[sig.order,]
NMF.Exp.rank2 <- na.omit(NMF.Exp.rank2) 
group <- predict(nmf.rank2) 
table(group)
dev.off()
consensusmap(nmf.rank2,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank2)]),
            annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
#consensusmap(nmf.rank2)
group_nmf= as.data.frame(group)


write.table(group, file="cluster_nmf.txt", sep="\t", quote=F, col.names = F)



### 比较临床信息（个性化）
####!注意修改
pdata=pdata2[group_list2=='IS',]

clinical=data.frame(time=pdata$`time after stroke (h):ch1`,group=group_nmf$group)
clinical$group=as.character(clinical$group)
library(ggstatsplot)
ggbarstats(clinical,x='time',y='group')

#### 如果变为连续型变量
clinical$time=as.numeric(clinical$time)
library(ggpubr)
ggboxplot(clinical,x='group',y='time',color = 'group')+stat_compare_means()

### 热图
data1=data[,group_nmf$group==1]
data2=data[,group_nmf$group==2]
### 三型
data3=data[,group_nmf$group==3]
data=cbind(data1,data2,data3)

data=na.omit(data)
anno=group_nmf
anno$group=paste0('Subtype',anno$group)
pheatmap::pheatmap(data,cluster_cols = F,scale = 'row',show_colnames = F,
                   border_color = NA,breaks = seq(-3,3,length.out = 100),
                   annotation_col = anno)

### 亚型间的比较
gsva_matrix=gsva_matrix[gsva_matrix$group=='IS',]
gsva_matrix$nmf_group=group_nmf$group
gsva_matrix$nmf_group=paste0('Subtype',gsva_matrix$nmf_group)

library(ggpubr)
library(ggsci)
ggboxplot(gsva_matrix,x='nmf_group',y='IS_Monocytes',color = 'nmf_group',
          add='jitter')+ggsci::scale_color_jco()+stat_compare_means()


############# 分型和免疫浸润######################################
library(IOBR)
mcp=IOBR::deconvo_mcpcounter(eset = exprSet2)
rownames(mcp)=mcp$ID

mcp$ID=NULL

mcp1=mcp[group_nmf$group==1,]
mcp2=mcp[group_nmf$group==2,]

## 三型
mcp3=mcp[group_nmf$group==3,]

rt=rbind(mcp1,mcp2,mcp3)
rt=as.data.frame(t(rt))
anno=group_nmf
anno$group=paste0('Subtype',anno$group)
pheatmap::pheatmap(rt,cluster_cols = F,scale = 'row',show_colnames = F,
                   border_color = NA,breaks = seq(-3,3,length.out = 100),
                   annotation_col = anno)


rt=as.data.frame(t(rt))
rt$nmf_group=c(rep('Subtype1',nrow(mcp1)),rep('Subtype2',nrow(mcp2)),
               rep('Subtype3',nrow(mcp3)))

rt=gather(rt,key='Celltype',value='Abundance',-c(nmf_group))

#a=grep(zhuhong',rt1$ID)


library(ggsci)
library(ggplot2)
library(ggpubr)
ggboxplot(
  rt,
  x = "Celltype",
  y = "Abundance",
  color = "black",
  fill = "nmf_group",
  xlab = "",
  ylab = "Cell composition",title = 'All patients'
) +
  stat_compare_means(
    aes(group = nmf_group),
    label = "p.signif", 
    method = "t.test",
    hide.ns = T,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))+scale_fill_jco()


