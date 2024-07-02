library(IOBR)
#load('lasso_genes.Rdata')
load('gene_rfe_rf.Rdata')
#rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)

## 测试集
load('GSE59867_46_111.Rdata')
pdata1=pdata
group_list1=c(rep('CT',46),rep('ASC',111))
exp1=rt
load('GSE62646_14_28.Rdata')
pdata2=pdata
group_list2=c(rep('CT',14),rep('ASC',28))
exp2=rt
### 确定好训练和验证集，需要尝试
train=exp1
test=exp2
rownames(train) <- trimws(rownames(train), "both")
rownames(test) <- trimws(rownames(test), "both")
colnames(train) <- trimws(colnames(train), "both")
colnames(test) <- trimws(colnames(test), "both")
############# 训练集######################################
mcp_train=IOBR::deconvo_mcpcounter(eset = train)
rownames(mcp_train)=mcp_train$ID

genes=gene
mcp_train$ID=NULL
exp_train=train[rownames(train) %in% genes,]
mcp_train=as.data.frame(t(mcp_train))


### 训练集制作第一个#######
mcp_a=mcp_train[,1]
exp_a=exp_train[,1]

table_a=matrix(rep(1, length(genes)*10),ncol =10,nrow =length(genes))


rownames(table_a)=rownames(exp_train)
colnames(table_a)=rownames(mcp_train)


for (i in 1:ncol(table_a)) {
  for (j in 1:nrow(table_a)) {
    table_a[j,i]=mcp_a[i] / exp_a[j]
  }
}
pdf(file="FIG 6 train_exm.pdf", width=8, height=8)
pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F)
dev.off()
### 训练集制作第二到最后一个#########
for (k in 2:ncol(mcp_train)) {
  mcp_k=mcp_train[,k]
  exp_k=exp_train[,k]
  table_k=matrix(rep(1,length(genes)*10),ncol =10,nrow = nrow(exp_train))
  rownames(table_k)=rownames(exp_train)
  colnames(table_k)=rownames(mcp_train)
  for (i in 1:ncol(table_k)) {
    for (j in 1:nrow(table_k)) {
      table_k[j,i]=mcp_k[i] / exp_k[j]
    }
  }
  table_a=rbind(table_a,table_k)
}


library(keras)


x_train <- array_reshape(table_a, dim = c(ncol(mcp_train), length(genes), 10))

pheatmap::pheatmap(x_train[60,,],cluster_rows = F,cluster_cols = F)




############# 测试集######################################
mcp_test=IOBR::deconvo_mcpcounter(eset = test)
rownames(mcp_test)=mcp_test$ID


mcp_test$ID=NULL
exp_test=test[rownames(test) %in% genes,]
mcp_test=as.data.frame(t(mcp_test))


### 训练集制作第一个#######
mcp_a=mcp_test[,3]
exp_a=exp_test[,3]

table_a=matrix(rep(1, length(genes)*10),ncol =10,nrow =length(genes))


rownames(table_a)=rownames(exp_test)
colnames(table_a)=rownames(mcp_test)


for (i in 1:ncol(table_a)) {
  for (j in 1:nrow(table_a)) {
    table_a[j,i]=mcp_a[i] / exp_a[j]
  }
}
pdf(file="FIG 6 txt_exm.pdf", width=8, height=8)
pheatmap::pheatmap(table_a,cluster_rows = F,cluster_cols = F)
dev.off()
### 训练集制作第二到最后一个#########
for (k in 2:ncol(mcp_test)) {
  mcp_k=mcp_test[,k]
  exp_k=exp_test[,k]
  table_k=matrix(rep(1,length(genes)*10),ncol =10,nrow = nrow(exp_test))
  rownames(table_k)=rownames(exp_test)
  colnames(table_k)=rownames(mcp_test)
  for (i in 1:ncol(table_k)) {
    for (j in 1:nrow(table_k)) {
      table_k[j,i]=mcp_k[i] / exp_k[j]
    }
  }
  table_a=rbind(table_a,table_k)
}


library(keras)
#dataset_mnist()

x_test <- array_reshape(table_a, dim = c(ncol(mcp_test), length(genes), 10))

pheatmap::pheatmap(x_test[40,,],cluster_rows = F,cluster_cols = F)




###############卷积，注意修改数目 ######################！！！！！！！！！！！

### 基因个数你不一定是11个
x_train <- array_reshape(x_train, dim = c(157,length(genes), 10, 1))
x_test <- array_reshape(x_test, dim = c(42, length(genes), 10, 1))

y_train <- ifelse(group_list1=='ASC',1,0)
y_test <- ifelse(group_list2=='ASC',1,0)


model <- keras_model_sequential() 
model %>%   
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), padding = 'same',  input_shape = c(10, 10, 1)) %>%
  layer_activation('relu') %>%
  #layer_activation(activation = 'relu') %>%
  #layer_max_pooling_2d(pool_size=c(2, 2), strides=c(2, 2)) %>%
  layer_conv_2d(filters = 16, kernel_size = c(2, 2), 
                dilation_rate = c(1,1), activation = 'softplus', padding = 'same') %>%
  layer_max_pooling_2d(pool_size=c(2, 2)) %>%
  layer_flatten() %>%
  layer_dense(64, activation = 'relu') %>%
  layer_dropout(0.5) %>%
  layer_dense(1, activation = 'sigmoid')

summary(model)

model %>% compile(
  loss = 'binary_crossentropy',
  optimizer ='adam',
  metrics = c('accuracy')
)

history <- model %>% fit(
  x_train, y_train, 
  epochs = 200
)
pdf(file="FIG 6 train_leve.pdf", width=6, height=4)
plot(history)
dev.off()
prob=model %>% predict(x_train) 
pdf(file="FIG 6 train_deep_ROC.pdf", width=6, height=6)
pROC::plot.roc(factor(y_train),as.numeric(prob),print.auc=T,print.thres=T)
dev.off()
prob_test=model %>% predict(x_test) 
pdf(file="FIG 6 txt_deep_ROC.pdf", width=6, height=6)
pROC::plot.roc(factor(y_test),as.numeric(prob_test),
               print.auc=T,print.thres=T)
dev.off()

###在最佳cutoff下 ，看测试集表现
prob_bi=ifelse(prob_test>0.859,1,0)

table(prob_bi,y_test)




####################### 神经网络，注意修改数目##########################
x_train <- array_reshape(x_train, dim = c(157, length(genes), 10))
x_test <- array_reshape(x_test, dim = c(42, length(genes), 10))

#loading keras library
library(keras)
#loading the keras inbuilt mnist dataset
# converting a 2D array into a 1D array for feeding into the MLP and normalising the matrix
train_x <- array(x_train, dim = c(dim(x_train)[1], prod(dim(x_train)[-1]))) 
test_x <- array(x_test, dim = c(dim(x_test)[1], prod(dim(x_test)[-1])))
#converting the target variable to once hot encoded vectors using keras inbuilt function
train_y=y_train
test_y=y_test
#defining a keras sequential model
model <- keras_model_sequential()
#defining the model with 1 input layer[784 neurons], 1 hidden layer[784 neurons] with dropout rate 0.4 and 1 output layer[10 neurons]
#i.e number of digits from 0 to 9
model %>% 
  layer_dense(units = length(genes)*10, input_shape = length(genes)*10) %>% 
  layer_dropout(rate=0.4)%>%
  layer_activation(activation = 'relu') %>% 
  layer_dense(units = 1) %>% 
  layer_activation(activation = 'sigmoid')
#compiling the defined model with metric = accuracy and optimiser as adam.
model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = 'adam',
  metrics = c('accuracy')
)
#fitting the model on the training dataset


history <- model %>% fit(
  train_x, train_y, 
  epochs = 200,batch_size = 10
)


history <- model %>% fit(
  x_train, y_train, 
  epochs = 200,batch_size = 10
)

plot(history)

prob=model %>% predict(train_x) 
pdf(file="FIG 6 train_NE_ROC.pdf", width=6, height=6)
pROC::plot.roc(factor(train_y),as.numeric(prob),print.auc=T,print.thres=T)
dev.off()
prob_bi=ifelse(prob>0.882,1,0)

table(prob_bi,y_train)


prob_test=model %>% predict(test_x) 
pdf(file="FIG 6 txt_NE_ROC.pdf", width=6, height=6)
pROC::plot.roc(factor(test_y),as.numeric(prob_test),
               print.auc=T,print.thres=T,smooth=T)
dev.off()
prob_bi=ifelse(prob_test>0.5,1,0)

table(prob_bi,y_test)
