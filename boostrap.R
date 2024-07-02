rm(list=ls())
gc()

#----训练集boostrap----


load('GSE59867_46_111.Rdata')
pdata1=pdata
group_list1=c(rep('CT',46),rep('ASC',111))
exp1=rt
#gene=read.table("importanceGene.RF.txt", header=T, sep="\t", check.names=F)
#gene=as.vector(gene[,1])
load('gene_rfe_rf.Rdata')
#gene=ss
rt=rt[rownames(rt) %in% gene,]
gene=rownames(rt)
rt=as.data.frame(t(rt))
group_list1=as.data.frame(group_list1)
rt$group=group_list1$group_list1

rt1=rt[,colnames(rt)!="group"]
rt1 <- apply(rt1, 2, scale)
rt1 <- as.data.frame(rt1)
rt1$group=rt$group
rownames(rt1)=rownames(rt)

rt=rt1


rt$group=ifelse(rt$group=='CT',0,1)

glm=glm(group ~.,data = rt,family= binomial(link='logit'))

OR=exp(glm$coefficients)[2:(length(gene)+1)]
OR_CI=exp(confint(glm,level = 0.95))
OR.95L=OR_CI[2:(length(gene)+1),1]
OR.95H=OR_CI[2:(length(gene)+1),2]


df=data.frame(gene=gene,OR=OR,OR.95L=OR.95L,OR.95H=OR.95H)


#install.packages("forestploter")
library(forestploter)

df$`OR (95% CI)` <- ifelse(is.na(df$OR), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   df$OR,df$OR.95L,df$OR.95H))


df$"" = paste(rep(" ",57), collapse = " ")
plot <- forest(df[, c(1, 6, 5)],
               est = df$OR,
               lower =df$OR.95L,
               upper =df$OR.95H,ref_line = 1,
               ci_column = 2)
pdf(file="FIG 3 train_log.pdf", width=8, height=4)
plot
dev.off()
pred=predict(glm,rt,type = 'response')
pred_df=as.data.frame(pred)
save(pred_df,file ='pred_logistic.Rdata')



library(pROC)
pdf(file="FIG 3 train_ROC.pdf", width=6, height=6)
plot.roc(rt$group,pred,print.auc=T,print.thres=T)
dev.off()
###bootstrap验证
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data) 
  # 得到预测的ROC
  glmROC=roc(rt$group,as.numeric(glm.pred))
  plot.roc(rt$group,glm.pred,add=T,col='#0072b5')
  # 返回ROC的AUC
  return(glmROC$auc)
} 

library(boot) 
set.seed(123456) 


# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)



#
print(results)
pdf(file="FIG 3 train_boot2.pdf", width=12, height=6)
plot(results)
dev.off()
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))

rm(outtab)

###boot 灵敏度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(sen)
} 

library(boot) 
set.seed(123456) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
pdf(file="FIG 3 train_boot_灵敏度.pdf", width=12, height=6)
plot(results)
dev.off()
## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))



###boot 特异度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(spe)
} 

library(boot) 
set.seed(123456) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
pdf(file="FIG 3 train_boot_特异度.pdf", width=12, height=6)
plot(results)
dev.off()
## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))



#install.packages('regplot')
library(regplot)
regplot(glm,observation = rt[2,],interval = 'confidence',clickable = T,title = 'Nomogram')


library(rms)
glm2=lrm(group~., data = rt,x=T,y=T)

cal<-calibrate(glm2, method = 'boot', B=1000, data = rt)
par(oma=c(3,3,3,3)) 
par(mar=c(6,5,4,3) + 0.1) 
plot(cal,
     xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",main = "Calibration Curve")

#install.packages('set')
#devtools::install_github("yikeshu0611/do")
#if (!require(ggDCA)) {
#  devtools::install_github("yikeshu0611/ggDCA")
#}
#install.packages("rmda")
library(rmda)
library(ggDCA)
library(ggplot2)
library(rms)
library(caret)


set.seed(123)
colnames(rt)

model1 <- decision_curve(group ~ CD247, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
model2 <- decision_curve(group ~ PDCD1LG2, data = rt,
                         thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)
colnames(rt)
#!!!!!!!!!!!!!!!!!!
model0 <- decision_curve(group ~ CD247 + `HLA-DRA` + LCK + CD3E + `HLA-DPA1` + `HLA-DQA1` + CD3D + PDCD1LG2 + CD3G,
                         data = rt,thresholds = seq(0, 0.4, by = 0.005), bootstraps = 10)

par(oma=c(1,1,1,1)) 
par(mar=c(6,5,4,3) + 0.1) 
plot_decision_curve( list(model1,model2,model0),
                     # curve.names = c("Baseline model", "Full model"),
                     # col = c("blue", "red"),
                     confidence.intervals = F,  #remove confidence intervals
                     cost.benefit.axis = FALSE, #remove cost benefit axis
                     legend.position = "bottomleft") #add the legend





#----测试集boostrap----
load('GSE62646_14_28.Rdata')
#pdata2=pdata
group_list2=c(rep('CT',14),rep('ASC',28))
#exp2=rt
#rt=exp2
#gene=read.table("importanceGene.RF.txt", header=T, sep="\t", check.names=F)
#gene=as.vector(gene[,1])
#load('key_marker_gene.Rdata')
#gene=ss
load('gene_rfe_rf.Rdata')
rt=rt[rownames(rt) %in% gene,]
gene=rownames(rt)
rt=as.data.frame(t(rt))
group_list2=as.data.frame(group_list2)
rt$group=group_list2$group_list2





rt1=rt[,colnames(rt)!="group"]
rt1 <- apply(rt1, 2, scale)
rt1 <- as.data.frame(rt1)
rt1$group=rt$group
rownames(rt1)=rownames(rt)

rt=rt1


rt$group=ifelse(rt$group=='CT',0,1)
#binomial    #probit  #identity
glm=glm(group ~.,data = rt,family= binomial(link='logit'),control=list(maxit=1000))

OR=exp(glm$coefficients)[2:(length(gene)+1)]
OR_CI=exp(confint(glm,level = 0.95))
OR.95L=OR_CI[2:(length(gene)+1),1]
OR.95H=OR_CI[2:(length(gene)+1),2]


df=data.frame(gene=gene,OR=OR,OR.95L=OR.95L,OR.95H=OR.95H)


#install.packages("forestploter")
library(forestploter)

df$`OR (95% CI)` <- ifelse(is.na(df$OR), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   df$OR,df$OR.95L,df$OR.95H))


df$"" = paste(rep(" ",57), collapse = " ")
plot <- forest(df[, c(1, 6, 5)],
               est = df$OR,
               lower =df$OR.95L,
               upper =df$OR.95H,ref_line = 1,
               ci_column = 2)
pdf(file="FIG 3 txt_log.pdf", width=8, height=4)
plot
dev.off()
pred=predict(glm,rt,type = 'response')
pred_df=as.data.frame(pred)
save(pred_df,file ='txt_pred_logistic_txt.Rdata')



library(pROC)
pdf(file="FIG 3 txt_ROC.pdf", width=6, height=6)
plot.roc(rt$group,pred,print.auc=T,print.thres=T)
dev.off()
###bootstrap验证
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data) 
  # 得到预测的ROC
  glmROC=roc(rt$group,as.numeric(glm.pred))
  plot.roc(rt$group,glm.pred,add=T,col='#0072b5')
  # 返回ROC的AUC
  return(glmROC$auc)
} 

library(boot) 
set.seed(123456) 


# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)



#
print(results)
pdf(file="FIG 3 txt_boot2.pdf", width=12, height=6)
plot(results)
dev.off()
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))

rm(outtab)

###boot 灵敏度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(sen)
} 

library(boot) 
set.seed(123456) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
pdf(file="FIG 3 txt_boot_灵敏度.pdf", width=12, height=6)
plot(results)
dev.off()
## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))



###boot 特异度
rsq <- function(formula, data, indices) { 
  # 通过data得到重抽样样本
  d <- data[indices,] 
  # 通过formula得到模型
  fit.glm <- glm(formula, data=d)
  # 用boot的模型在原数据集上预测
  glm.pred <- predict(fit.glm, data,type='response') 
  print(glm.pred)
  b=ifelse(glm.pred>0.523,1,0)
  outtab=table(b,rt$group)
  sen=outtab[2,2]/(outtab[1,2]+outtab[2,2])
  spe=outtab[1,1]/(outtab[1,1]+outtab[2,1])
  return(spe)
} 

library(boot) 
set.seed(123456) 

# 运行bootstrap
results <- boot(data=rt, statistic=rsq, 
                R=1000, formula=group~.)

#
print(results)
pdf(file="FIG 3 txt_boot_特异度.pdf", width=12, height=6)
plot(results)
dev.off()
## 和我们一样0.96
mean(results$t)

# 获取置信区间
boot.ci(results, type=c("perc", "bca"))






