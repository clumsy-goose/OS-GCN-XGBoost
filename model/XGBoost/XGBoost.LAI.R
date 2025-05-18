# 0. environment ----------------------------------------------------------------------------------------

library(xgboost)
library(data.table)
library(mlr)
library(ggpubr)
library(ggplot2)

rm(list = setdiff(ls(), lsf.str()))

# 1. Load data -----------------------------------------------------------------------------------

data<-read.csv(file = "../.././dataProcessedResult/phenotype_rpkm/LAI.csv",header = T,row.names='Entry_Treatment_Replicate') 

## Remove unwanted features
data=data[,-c(0:1)] 

# 2. XGBoost -------------------------------------------------------------------------------------

## Setting parameters

### About the data set
n=220 # number of genotypes
c=6  # number of samples per genotype

### About the run
k=20 # number of iteration
jj=k-1

### Hyperparameters
r=20 # number of rounds
colsample=0.33
eta=0.075 #0.075 #0.1
num_parallel_tree=1
subsample=0.25

### About the output structure
## Need to reset the following otherwise you'll get 'out of subscript' error
y=0
obs = matrix(nrow = k*n, ncol=c)
pred = matrix(nrow = k*n, ncol=c)
rmse = vector()
train.rmse = vector()
impt.out=vector('list',k*n)

for (i in 1:n){
  for (j in 0:jj) {
    
    test.index <-c(i,i+220,i+440,i+660,i+880,i+1100)
    testing<-data[test.index,]
    training<-data[-test.index,]
    
    #convert data frame to data table
    setDT(training) 
    setDT(testing)
    
    #using one hard encoding 
    train.trait <- training$LAI
    test.trait <- testing$LAI
    new_training <- model.matrix(~.+0,data = training[,-c("LAI"),with=F]) 
    new_testing <- model.matrix(~.+0,data = testing[,-c("LAI"),with=F])
    
    #preparing matrix 
    dtrain <- xgb.DMatrix(data = new_training,label = train.trait) 
    dtest <- xgb.DMatrix(data = new_testing,label=test.trait)
    watchlist <- list(train=dtrain, test=dtest)
    
    #user defined evaluation metric
    cor_metric <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      cor <- cor(preds,labels)
      list(metric = "cor", value = cor)}
    
    params <- list(booster = "gbtree", 
                   objective = "reg:squarederror", 
                   eta= eta,
                   gamma= 50,  
                   max_depth=6,
                   min_child_weight=1, 
                   eval_metric=cor_metric,
                   subsample=subsample, 
                   colsample_bytree=colsample,
                   num_parallel_tree=num_parallel_tree) 
    
    set.seed(j)
    bst.val<-xgb.train( params = params, 
                        data = dtrain, 
                        nrounds = r,
                        nfold = 5, 
                        showsd = T, 
                        stratified = T, 
                        print_every_n = 20, 
                        early_stop_round = 5, 
                        watchlist = watchlist,
                        maximize = F,
                        verbose = 0)
    
    y=y+1
    
    pred[y,1:c]<- predict(bst.val, dtest)
    obs[y,1:c]<-test.trait
    
    rmse[y]<- as.numeric(bst.val$evaluation_log[r,3])
    train.rmse[y]<-as.numeric(bst.val$evaluation_log[r,2])
    
    # extract important features
    importance_matrix <- xgb.importance(model = bst.val)
    impt.out[[y]]<-paste(importance_matrix$Feature,importance_matrix$Gain, sep = ",")
  }}

save(data,obs, pred, impt.out, rmse, train.rmse,k,n,jj,r,c,colsample,eta,num_parallel_tree,subsample,
     file="./LAI.RData")

## Organize cor

load(file = "./LAI.RData")

n=220
print(ncol(pred.mat))  # 实际列数
print(ncol(obs.mat))
print(ncol(pred))
pred.mat = matrix(pred, nrow=k)
obs.mat = matrix(obs, nrow=k)

# For each accession, calculate COR for each iteration 

cor.mat=matrix(nrow = jj+1,ncol = n)

for (i in 1:n){
  for (j in 1:(jj+1)){
    O=c(obs.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    P=c(pred.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    cor.mat[j,i]=cor(P,O)
  }}
sum(is.na(cor.mat))  # 统计 NA 的总数
mean(is.na(cor.mat)) # 计算 NA 占比
COR=vector()
#COR=colMeans(cor.mat)
COR <- colMeans(cor.mat, na.rm = TRUE)
mean(COR, na.rm = TRUE)  # 忽略 NA 计算全局均值
#mean(COR)   # 0.6811985

## Plot the output
apply(cor.mat,2,sd)
sd=apply(cor.mat,2,sd)
upper<-COR+sd
lower<-COR-sd

genotype <- factor(sapply(strsplit(rownames(data)[1:220], "_"), `[`, 1))

df <- data.frame(genotype=genotype,
                 Correlation=COR)

par(mar=c(5,6,4,2)+0.1)
par(mgp=c(5,6,4))

ggplot(data=df, aes(x=genotype, y=Correlation)) +
  geom_bar(stat="identity",fill="#af8dc3",
           width =0.4,position=position_dodge())+
  geom_errorbar(aes(ymin=Correlation, ymax=upper), width=.2,colour="#af8dc3",
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90,vjust=1))+
  theme_minimal()

# 3. Extract importantgenes  ----------------------------------------------------------------------

library(tidyverse)
library(data.table)

s=220 # number of genotypes, n
t=100 # number of iteration, jj+1
impt=vector("list",s)

for (i in 1:s){
  m=(i-1)*t+1
  n=t*i
  
  d2 <- as.data.frame(do.call(rbind,flatten(impt.out[m:n]))) %>% 
    separate(.,V1, into =c("Gene","Importance"), sep=",")
  
  ## convert importance from character to numeric
  d2$Importance=as.numeric(d2$Importance)
  
  ## convert data frame to data table for easy calculation
  setDT(d2)
  d2[,sum(Importance),by=Gene]
  output=d2[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
  ## calculate impt
  impt[[i]]=paste(output$Gene,output$SUM, sep = ",")
}

# Combine the impt
d3 <- as.data.frame(do.call(rbind,flatten(impt))) %>% 
  separate(.,V1, into =c("Gene","Importance"), sep=",")

d3$Importance=as.numeric(d3$Importance)

# Calculate composite score of importance
setDT(d3)
d3[,sum(Importance),by=Gene]
output=d3[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
nrow(output) #610

top_n <- 25
plot_data <- head(output, top_n)

# 绘制柱状图（纵向）
ggplot(plot_data, aes(x = reorder(Gene, SUM), y = SUM)) +
  geom_bar(stat = "identity", fill = "#55A868", width = 0.7) +
  labs(
    title = "Top Genes by Composite Importance Score (Leaf area index (%) (LAI))",
    x = "Gene",
    y = "Importance Score"
  ) +
  coord_flip() +  # 翻转坐标轴
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )
