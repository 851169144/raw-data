

#???ð?
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(9999999)      #????????
load('GSE59867_46_111.Rdata')
pdata1=pdata
group_list1=c(rep('CT',46),rep('ASC',111))
exp1=rt
gene=read.table("./alldiff_sig.csv" ,sep= ",",header = T,row.names = 1)
gene=row.names(gene)



rt=rt[rownames(rt) %in% gene,]
gene=rownames(rt)
rt=as.data.frame(t(rt))
group_list1=as.data.frame(group_list1)
rt$Type=group_list1$group_list1

data=rt

#?????ݽ??з???
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#RF????ɭ????ģ??
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

#SVM????ѧϰģ??
mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

#XGBģ??
mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)

#GLMģ??
mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)


#????Ԥ?⺯??
p_fun=function(object, newdata){
	predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="CT", 0, 1)

#RF????ɭ????ģ??Ԥ??????
explainer_rf=explain(mod_rf, label = "RF",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_rf=model_performance(explainer_rf)
#SVM????ѧϰģ??Ԥ??????
explainer_svm=explain(mod_svm, label = "SVM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_svm=model_performance(explainer_svm)
#XGBģ??Ԥ??????
explainer_xgb=explain(mod_xgb, label = "XGB",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)
#GLMģ??Ԥ??????
explainer_glm=explain(mod_glm, label = "GLM",
                         data = test, y = yTest,
                         predict_function = p_fun,
                         verbose = FALSE)
mp_glm=model_performance(explainer_glm)

#???????ַ????Ĳв???ۼƷֲ?ͼ
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm)
print(p1)
dev.off()

#???????ַ????Ĳв?????ͼ
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, geom = "boxplot")
print(p2)
dev.off()


#????ROC????
pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
legend('bottomright',
	   c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
	     paste0('SVM: ',sprintf("%.03f",roc2$auc)),
	     paste0('XGB: ',sprintf("%.03f",roc3$auc)),
	     paste0('GLM: ',sprintf("%.03f",roc4$auc))),
	   col=c("red","blue","green","yellow"), lwd=2, bty = 'n')
dev.off()

#?????ַ??????л???????Ҫ?Է???,?õ????ַ?????????Ҫ???��?
importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)
#???ƻ?????Ҫ??ͼ??
pdf(file="importance.pdf", width=7, height=10)
plot(importance_rf[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_svm[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_xgb[c(1,(ncol(data)-8):(ncol(data)+1)),],
	 importance_glm[c(1,(ncol(data)-8):(ncol(data)+1)),])
dev.off()
#??????Ҫ???��????ߵĻ???
geneNum=10    #???û???????Ŀ
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)


#------临床-----
library(rms)
library(rmda)

load('GSE59867_46_111.Rdata')
pdata1=pdata
group_list1=c(rep('CT',46),rep('ASC',111))
exp1=rt
gene=read.table("importanceGene.RF.txt", header=T, sep="\t", check.names=F)
gene=as.vector(gene[,1])


rt=rt[rownames(rt) %in% gene,]
gene=rownames(rt)
rt=as.data.frame(t(rt))
group_list1=as.data.frame(group_list1)
rt$Type=group_list1$group_list1

data=rt[,colnames(rt)!='Type']
paste(colnames(data), collapse="+")

ddist=datadist(rt)
options(datadist="ddist")

#????ģ?ͣ?????????ͼ
lrmModel=lrm(Type~ PTP4A2+CD3D+ABI3+GZMM+CALR+DIAPH1+ANXA6+BTN3A2+NCR3+ACTB, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.0001,0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of Disease")
#????????ͼ
pdf("Nomo.pdf", width=8, height=6)
plot(nomo)
dev.off()

#????У׼????
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()

#???ƾ???????
rt$Type=ifelse(rt$Type=="CT", 0, 1)
dc=decision_curve(Type ~ PTP4A2+CD3D+ABI3+GZMM+CALR+DIAPH1+ANXA6+BTN3A2+NCR3+ACTB, 
                  data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
#????DCAͼ??
pdf(file="DCA.pdf", width=5.5, height=5.5)

plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="red",
                    confidence.intervals=FALSE,
                    standardize=FALSE)

dev.off()

