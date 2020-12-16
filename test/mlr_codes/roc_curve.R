setwd('G:/lab/DCcomboNet/mlr_models/')
library(mlr)
library(pROC)
library(ggplot2)
library(ggthemes)
library(grid)
load('results2/ranger_result_LOO_0.2.Rdata')
# plot(density(r_ranger$pred$data$prob.P))
# lines(density(r_ranger$pred$data$prob.N), col='red')
# length(r_ranger$pred$data$prob.P)
# length(r_ranger$pred$data$prob.N)
rf_pred = getRRPredictions(r_ranger)
rf_df = data.frame(prob = rf_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
rf_df = rf_df[order(rf_df$prob, decreasing = T),]
rf_df$rank = 1:nrow(rf_df)

rf_modelroc <- roc(rf_df$label, rf_df$prob,
                   ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                   plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                   print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/ranger_result_LOO_0.2_roc.pdf')
plot(rf_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

load('results2/ksvm_result_LOO_0.2.Rdata')
ksvm_pred = getRRPredictions(r_ksvm)
ksvm_df = data.frame(prob = ksvm_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
ksvm_df = ksvm_df[order(ksvm_df$prob, decreasing = T),]
ksvm_df$rank = 1:nrow(ksvm_df)
ksvm_modelroc <- roc(ksvm_df$label, ksvm_df$prob,
                     ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/ksvm_result_LOO_0.2_roc.pdf')
plot(ksvm_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

load('results2/elasticnet_result_LOO_0.2.Rdata')
elasticnet_pred = getRRPredictions(r_elasticnet)
elasticnet_df = data.frame(prob = elasticnet_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
elasticnet_df = elasticnet_df[order(elasticnet_df$prob, decreasing = T),]
elasticnet_df$rank = 1:nrow(elasticnet_df)
elasticnet_modelroc <- roc(elasticnet_df$label, elasticnet_df$prob,
                           ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                           plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                           print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/elasticnet_result_LOO_0.2_roc.pdf')
plot(elasticnet_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()


load('results2/glm_result_LOO_0.2.Rdata')
glm_pred = getRRPredictions(r_glm)
glm_df = data.frame(prob = glm_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
glm_df = glm_df[order(glm_df$prob, decreasing = T),]
glm_df$rank = 1:nrow(glm_df)
glm_modelroc <- roc(glm_df$label, glm_df$prob,
                    ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                    plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                    print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/glm_result_LOO_0.2_roc.pdf')
plot(glm_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()



load('results2/naiveBayes_result_LOO_0.2.Rdata')
naiveBayes_pred = getRRPredictions(r_naiveBayes)
naiveBayes_df = data.frame(prob = naiveBayes_pred$data$prob.P,
                           label = c(rep('P',185),rep('N',6639-185)),rank=0)
naiveBayes_df = naiveBayes_df[order(naiveBayes_df$prob, decreasing = T),]
naiveBayes_df$rank = 1:nrow(naiveBayes_df)
naiveBayes_modelroc <- roc(naiveBayes_df$label, naiveBayes_df$prob,
                           ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                           plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                           print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/naiveBayes_result_LOO_0.2_roc.pdf')
plot(naiveBayes_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

load('results2/xgboost_result_LOO_0.2.Rdata')
xgboost_pred = getRRPredictions(r_xgboost)
xgboost_df = data.frame(prob = xgboost_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
xgboost_df = xgboost_df[order(xgboost_df$prob, decreasing = T),]
xgboost_df$rank = 1:nrow(xgboost_df)
xgboost_modelroc <- roc(xgboost_df$label, xgboost_df$prob,
                        ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/xgboost_result_LOO_0.2_roc.pdf')
plot(xgboost_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()


load('results2/nnet_result_LOO_0.2.Rdata')
nnet_pred = getRRPredictions(r_nnet)
nnet_df = data.frame(prob = nnet_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
nnet_df = nnet_df[order(nnet_df$prob, decreasing = T),]
nnet_df$rank = 1:nrow(nnet_df)
nnet_modelroc <- roc(nnet_df$label, nnet_df$prob,
                     ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/nnet_result_LOO_0.2_roc.pdf')
plot(nnet_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()


load('results2/nnTrain_result_LOO_0.2.Rdata')
nnTrain_pred = getRRPredictions(r_nnTrain)
nnTrain_df = data.frame(prob = nnTrain_pred$data$prob.P, label = c(rep('P',185),rep('N',6639-185)),rank=0)
nnTrain_df = nnTrain_df[order(nnTrain_df$prob, decreasing = T),]
nnTrain_df$rank = 1:nrow(nnTrain_df)
nnTrain_modelroc <- roc(nnTrain_df$label, nnTrain_df$prob,
                        ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/nnTraint_result_LOO_0.2_roc.pdf')
plot(nnTrain_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()


######################
#     DComboNet      #
######################

outputdir = 'G:/lab/DCcomboNet/Rpackage/'
# setwd(outputdir)
require(data.table)

model = 'L1'
cellline = NULL
eta=1
r=0.7

drugnet_feature2 = read.csv(paste0(outputdir,'OCL_LY3/drug_net/features2.csv'),stringsAsFactors = F)
drugfeature = drugnet_feature2[drugnet_feature2$integrated_score2>=0.2,]
drugpair_all = unique(drugfeature[c(1,2)])

x=0.5;y=0.5;z=0;A=0.5;B=0.5;C=0
resultdir = paste0(outputdir,'/wwin_LOOCV8/L1_model/x_',x,'_y_',y,'_z_',z,'_a_',A,'_b_',B,'_c_',C,'_r_0.7/')
a <- list.files(paste0(resultdir,model,'_result/drugrank/'))
dir <-paste0(resultdir,model,'_result/drugrank/',a)
n <- length(dir);n
b = gsub('_rank.csv','',a)


A_B = data.table(b)
drugpair_positive = data.frame(A_B[,c('A','B') := tstrsplit(b,'_',fix = T)][])[c('A','B')]

names(drugpair_positive) = c('drugseed','DrugID')
drugpair_positive$Pred2 = rep(1,nrow(drugpair_positive))
drugpair = strsplit(b[1],'_')[[1]]
drugrank <- read.csv(file = dir[1],header=T, stringsAsFactors = F)
drugseed = drugpair[1]; drug.pred = drugpair[2]
drugrank$drugseed = drugseed
drugrank$Pred1 = rep(0,nrow(drugrank))
drugrank$Pred2 = rep(0,nrow(drugrank))
drugrank[which(drugrank$DrugID == drug.pred),]$Pred1 = 1
m = drugpair_positive[which(drugpair_positive$drugseed == drugseed),2]
if(length(m)>1){for(i in 1:length(m)){
  drugrank[which(drugrank$DrugID == m[i]),]$Pred2 = 1
}}else if(m == 1){
  drugrank[which(drugrank$DrugID == m),]$Pred2 = 1
}

for(i in 2:n){
  #i = 2
  drugpair = strsplit(b[i],'_')[[1]]
  drugrank2 <- read.csv(file = dir[i],header=T, stringsAsFactors = F)
  drugseed = drugpair[1]; drug.pred = drugpair[2]
  drugrank2$drugseed = drugseed
  drugrank2$Pred1 = rep(0,nrow(drugrank2))
  drugrank2$Pred2 = rep(0,nrow(drugrank2))
  drugrank2[which(drugrank2$DrugID == drug.pred),]$Pred1 = 1
  tmp = drugpair_positive[which(drugpair_positive$drugseed == drugseed),2]
  for(j in 1:length(tmp)){
    if(length(intersect(drugrank2$DrugID,tmp[j]))!=0){
      drugrank2[which(drugrank2$DrugID == tmp[j]),]$Pred2 = 1
    }
  }
  
  drugrank = rbind(drugrank, drugrank2)
}

drugrank = drugrank[c(4,1,2,3,5,6)]

drugrank = unique(drugrank[(drugrank$Pred1 == 0 & drugrank$Pred2 != 1)|drugrank$Pred1 == 1,])
save(drugrank, file='results2/DComboNet_result.Rdata')

a=0;
b=1;
y = drugrank$Score
Ymax=max(drugrank$Score);
Ymin=min(drugrank$Score);
k=(b-a)/(Ymax-Ymin);
norY=a+k*(y-Ymin)
drugrank$norm_Score = norY

pred <- prediction( drugrank$Score, drugrank$Pred1)
perf <- performance( pred, 'tpr', 'fpr' )
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AucValue<-round(AUCValue,3)
AucValue


DComboNet_modelroc <- roc(drugrank$Pred1, drugrank$norm_Score,
                        ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                        plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                        print.auc=TRUE, show.thres=TRUE)

pdf('results2_plot/DComboNet_result_LOO_0.2_roc.pdf')
plot(DComboNet_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

pdf('results2_plot/model_roc.pdf', width=10,height=10)
par(mfrow= c(3,3))
plot(rf_modelroc,
     print.auc=TRUE,
     auc.polygon=TRUE, 
     grid=c(0.1, 0.2),
     grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", 
     print.thres=TRUE)
plot(ksvm_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(elasticnet_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(glm_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(naiveBayes_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(xgboost_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(nnet_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(nnTrain_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(DComboNet_modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()


 rf_roc = data.frame(model = 'rf', x=1-rf_modelroc$specificities, y=rf_modelroc$sensitivities)
 ksvm_roc = data.frame(model = 'ksvm', x=1-ksvm_modelroc$specificities, y=ksvm_modelroc$sensitivities)
 elasticnet_roc = data.frame(model = 'elasticnet', x=1-elasticnet_modelroc$specificities, y=elasticnet_modelroc$sensitivities)
 glm_roc = data.frame(model = 'glm', x=1-glm_modelroc$specificities, y=glm_modelroc$sensitivities)
 naiveBayes_roc = data.frame(model = 'naiveBayes', x=1-naiveBayes_modelroc$specificities, y=naiveBayes_modelroc$sensitivities)
 xgboost_roc = data.frame(model = 'xgboost', x=1-xgboost_modelroc$specificities, y=xgboost_modelroc$sensitivities)
 nnet_roc = data.frame(model = 'nnet', x=1-nnet_modelroc$specificities, y=nnet_modelroc$sensitivities)
 DComboNet_roc = data.frame(model = 'DComboNet', x=1-DComboNet_modelroc$specificities, y=DComboNet_modelroc$sensitivities)

dat = rbind(DComboNet_roc, rf_roc, ksvm_roc, elasticnet_roc, glm_roc, naiveBayes_roc, xgboost_roc, nnet_roc)
dat$model = factor(dat$model, levels = unique(dat$model))

 nnTrain_modelroc = data.frame(model = 'nnTrain', x=1-nnTrain_modelroc$specificities, y=nnTrain_modelroc$sensitivities)

p = ggplot(dat) + 
  geom_line(aes(x,y,colour = model),size = 0.3)  +  
  labs(x="1 - Specificity", y="Sensitivity") + 
  scale_colour_manual(breaks=unique(dat$model),
                      values = unique(dat$model),
                      name="Different Methods (AUC)",
                      labels=c(paste0('DComboNet (',round(DComboNet_modelroc$auc[1],3),')'),
                               paste0('rf (',round(rf_modelroc$auc[1],3),')'),
                               paste0('ksvm (',round( ksvm_modelroc$auc[1],3),')'),
                               paste0('elasticnet (',round( elasticnet_modelroc$auc[1],3),')'),
                               paste0('glm (',round(glm_modelroc$auc[1],3),')'),
                               paste0('naiveBayes (',round(naiveBayes_modelroc$auc[1],3),')'),
                               paste0('xgboost (',round(xgboost_modelroc$auc[1],3),')'),
                               paste0('nnet (',round(nnet_modelroc$auc[1],3),')'))) + 
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.grid.major=element_line(colour=NA)) + 
  theme_set(theme_bw()) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed",colour='grey',alpha=0.8) +
  theme(legend.position=c(0.8,0.2))

p

ggsave(file='results2_plot/model_roc2.pdf',p, width = 10, height = 10)

######################################################
# rf_modelroc
# ksvm_modelroc
# elasticnet_modelroc
# glm_modelroc
# xgboost_modelroc
# naiveBayes_modelroc
# nnet_modelroc
# nnTrain_modelroc
# DComboNet_modelroc

rf_measures = coords(rf_modelroc, "best", ret= 'all',transpose = FALSE);rownames(rf_measures) = 'rf';rf_measures$auc =  rf_modelroc$auc[1]
ksvm_measures = coords(ksvm_modelroc, "best", ret= 'all',transpose = FALSE);rownames(ksvm_measures) = 'ksvm';ksvm_measures$auc =  ksvm_modelroc$auc[1]
elasticnet_measures = coords(elasticnet_modelroc, "best", ret= 'all',transpose = FALSE);rownames(elasticnet_measures) = 'elasticnet';elasticnet_measures$auc =  elasticnet_modelroc$auc[1]
glm_measures = coords(glm_modelroc, "best", ret= 'all',transpose = FALSE);rownames(glm_measures) = 'glm';glm_measures$auc =  glm_modelroc$auc[1]
xgboost_measures = coords(xgboost_modelroc, "best", ret= 'all',transpose = FALSE);rownames(xgboost_measures) = 'xgboost';xgboost_measures$auc =  xgboost_modelroc$auc[1]
naiveBayes_measures = coords(naiveBayes_modelroc, "best", ret= 'all',transpose = FALSE);rownames(naiveBayes_measures) = 'naiveBayes';naiveBayes_measures$auc =  naiveBayes_modelroc$auc[1]
nnet_measures = coords(nnet_modelroc, "best", ret= 'all',transpose = FALSE);rownames(nnet_measures) = 'nnet';nnet_measures$auc =  nnet_modelroc$auc[1]
nnTrain_measures = coords(nnTrain_modelroc, "best", ret= 'all',transpose = FALSE);rownames(nnTrain_measures) = 'nnTrain';nnTrain_measures$auc =  nnTrain_modelroc$auc[1]
DComboNet_measures = coords(DComboNet_modelroc, "best", ret= 'all',transpose = FALSE);rownames(DComboNet_measures) = 'DComboNet';DComboNet_measures$auc =  DComboNet_modelroc$auc[1]

model_compairson = rbind(rf_measures, ksvm_measures, elasticnet_measures,
      glm_measures, xgboost_measures, naiveBayes_measures,
      nnet_measures, nnTrain_measures, DComboNet_measures)

model_compairson[c("auc","specificity","sensitivity","accuracy","tpr","precision","recall","F1")]

write.csv(model_compairson,'model_compairson.csv',quote=F)

model_compairson$model = rownames(model_compairson)

p1 = ggplot(model_compairson, aes(x = model, y = auc))+ 
     geom_bar(stat = "identity", fill='#88CFF0') +
     geom_text(aes(label=round(auc,3)), vjust=1.6, color="black", size=3.5)
# p1

p2 = ggplot(model_compairson, aes(x = model, y = tpr))+ 
    geom_bar(stat = "identity", fill='#88CFF0') +
     geom_text(aes(label=round(tpr,3)), vjust=1.6, color="black", size=3.5)
# p2

pdf('results2_plot/auc_tpr_res.pdf', width=8, height=8)

grid.newpage()  ##新建页面
pushViewport(viewport(layout = grid.layout(2,1))) ####将页面分成2*2矩阵
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p1, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图c
print(p2, vp = vplayout(2,1))   ###将(2,1)的位置画图b

dev.off()


save(rf_modelroc,
     ksvm_modelroc,
     elasticnet_modelroc,
     glm_modelroc,
     xgboost_modelroc,
     naiveBayes_modelroc,
     nnet_modelroc,
     nnTrain_modelroc,
     DComboNet_modelroc,
     file = 'model_measures.Rdata')

setwd('/picb/bigdata/project/FengFYM/mlr_models')
setwd('G:/lab/DCcomboNet/mlr_models/')
library(pROC)

load('model_measures.Rdata')
model_compr = data.frame(model1 = c('rf','ksvm','elasticnet','glm','xgboost','naiveBayes','nnet','nnTrain'),
                         model2 = 'DComboNet',
                         venkatraman_pvalue = 0,
                         delong_pvalue = 0)


rf_DComboNet_delong = roc.test(rf_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='rf',]$delong_pvalue = rf_DComboNet_delong$p.value

ksvm_DComboNet_delong = roc.test(ksvm_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='ksvm',]$delong_pvalue = ksvm_DComboNet_delong$p.value

elasticnet_DComboNet_delong = roc.test(elasticnet_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='elasticnet',]$delong_pvalue = elasticnet_DComboNet_delong$p.value

glm_DComboNet_delong = roc.test(glm_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='glm',]$delong_pvalue = glm_DComboNet_delong$p.value

xgboost_DComboNet_delong = roc.test(xgboost_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='xgboost',]$delong_pvalue = xgboost_DComboNet_delong$p.value

naiveBayes_DComboNet_delong = roc.test(naiveBayes_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='naiveBayes',]$delong_pvalue = naiveBayes_DComboNet_delong$p.value

nnet_DComboNet_delong = roc.test(nnet_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='nnet',]$delong_pvalue = nnet_DComboNet_delong$p.value

nnTrain_DComboNet_delong = roc.test(nnTrain_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="delong")
model_compr[model_compr$model1=='nnTrain',]$delong_pvalue = nnTrain_DComboNet_delong$p.value

#########################

rf_DComboNet_venkatraman = roc.test(rf_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='rf',]$venkatraman_pvalue = rf_DComboNet_venkatraman$p.value

ksvm_DComboNet_venkatraman = roc.test(ksvm_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='ksvm',]$venkatraman_pvalue = ksvm_DComboNet_venkatraman$p.value

elasticnet_DComboNet_venkatraman = roc.test(elasticnet_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='elasticnet',]$venkatraman_pvalue = elasticnet_DComboNet_venkatraman$p.value

glm_DComboNet_venkatraman = roc.test(glm_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='glm',]$venkatraman_pvalue = glm_DComboNet_venkatraman$p.value

xgboost_DComboNet_venkatraman = roc.test(xgboost_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='xgboost',]$venkatraman_pvalue = xgboost_DComboNet_venkatraman$p.value

naiveBayes_DComboNet_venkatraman = roc.test(naiveBayes_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='naiveBayes',]$venkatraman_pvalue = naiveBayes_DComboNet_venkatraman$p.value

nnet_DComboNet_venkatraman = roc.test(nnet_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='nnet',]$venkatraman_pvalue = nnet_DComboNet_venkatraman$p.value

nnTrain_DComboNet_venkatraman = roc.test(nnTrain_modelroc, DComboNet_modelroc, reuse.auc=FALSE,method="venkatraman")
model_compr[model_compr$model1=='nnTrain',]$venkatraman_pvalue = nnTrain_DComboNet_venkatraman$p.value

write.csv(model_compr, 'model_compr_result.csv',row.names = F, quote = F)
