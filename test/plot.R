setwd('G:/lab/Projects/p1_DComboNet/Rpackage/OCL_LY3/')

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
#library(ggplot2)
set.seed(123)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plyr)
library(pROC)
library(ggplot2)

label = read.csv('drug_net/drug_features.csv', stringsAsFactors = F)
label = unique(label[c(1,2,8)])

all_rank = read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/result_all.csv')
all_rank = all_rank[-c(1,2),]

a=0;
b=1;
y = all_rank$Score_ABplusBA
Ymax=max(all_rank$Score_ABplusBA);
Ymin=min(all_rank$Score_ABplusBA);
k=(b-a)/(Ymax-Ymin);
norY=a+k*(y-Ymin)
all_rank$norm_Score = norY


#===============================================================================
#      Figure 3 L1 model ROC curve
#===============================================================================


drug_class = read.csv('data/druginfo/drug_class.csv', header = T, stringsAsFactors = F)
# drug_class = read.csv('/picb/bigdata/project/FengFYM/DComboNetV2/OCL_LY3/druginfo/drug_class.csv', header = T, stringsAsFactors = F)
classlist = sort(unique(drug_class$class),decreasing = T)[c(2,1,3)]
drugRWRfile <- all_rank
names(drugRWRfile) <- c('drugseed','DrugID', "Score_AB","Rank_AB","Score_BA","Rank_BA","Label","Score_ABplusBA","index","Rank_ABplusBA","norm_Score")
names(drug_class) = c('drugseed','drugseed.class')
drugRWRfile1 = merge(drugRWRfile,drug_class,by = 'drugseed', all.x = T)
names(drug_class) = c('DrugID','DrugID.class')
drugRWRfile2 = unique(merge(drugRWRfile1,drug_class,by = 'DrugID', all.x = T))
#drugRWRfile2 = drugRWRfile2[c('drugseed','drugseed.class','DrugID','DrugID.class','Score',"norm_Score",'rank','Pred1','Pred2')]
drugRWRfile.class1 = drugRWRfile2[drugRWRfile2$drugseed.class == classlist[1] & drugRWRfile2$DrugID.class == classlist[1], ]
drugRWRfile.class2 = drugRWRfile2[drugRWRfile2$drugseed.class == classlist[2] & drugRWRfile2$DrugID.class == classlist[2], ]
drugRWRfile.class3 = drugRWRfile2[drugRWRfile2$drugseed.class == classlist[3] & drugRWRfile2$DrugID.class == classlist[3], ]
drugRWRfile.class4 = drugRWRfile2[drugRWRfile2$drugseed.class == classlist[1] & drugRWRfile2$DrugID.class == classlist[2] |drugRWRfile2$drugseed.class == classlist[2] & drugRWRfile2$DrugID.class == classlist[1], ]
drugRWRfile.class5 = drugRWRfile2[drugRWRfile2$drugseed.class == classlist[1] & drugRWRfile2$DrugID.class == classlist[3] |drugRWRfile2$drugseed.class == classlist[3] & drugRWRfile2$DrugID.class == classlist[1], ]
drugRWRfile.class6 = drugRWRfile2[drugRWRfile2$drugseed.class == classlist[2] & drugRWRfile2$DrugID.class == classlist[3] |drugRWRfile2$drugseed.class == classlist[3] & drugRWRfile2$DrugID.class == classlist[2], ]
# write.csv(drugRWRfile.class1,paste0(resultdir,model,'_result/ROC/chemo_chemo_RWRfile.csv',sep = ''), row.names = F)
# write.csv(drugRWRfile.class2,paste0(resultdir,model,'_result/ROC/target_target_RWRfile.csv',sep = ''), row.names = F)
# write.csv(drugRWRfile.class3,paste0(resultdir,model,'_result/ROC/other_other_RWRfile.csv',sep = ''), row.names = F)
# write.csv(drugRWRfile.class4,paste0(resultdir,model,'_result/ROC/chemo_target_RWRfile.csv',sep = ''), row.names = F)
# write.csv(drugRWRfile.class5,paste0(resultdir,model,'_result/ROC/chemo_other_RWRfile.csv',sep = ''), row.names = F)
# write.csv(drugRWRfile.class6,paste0(resultdir,model,'_result/ROC/target_other_RWRfile.csv',sep = ''), row.names = F)

#drugrank2 = drugrank[(drugrank$Label == 0 & drugrank$Pred2 != 1)|drugrank$Label == 1,]
ROC.list = list()
if(length(unique(drugRWRfile$Label )) > 1){ ROC.list[[1]] <- roc(drugRWRfile$Label,drugRWRfile$norm_Score) }else{ROC.list[[1]] = NA}
if(length(unique(drugRWRfile.class1$Label )) > 1){ROC.list[[2]] <- roc(drugRWRfile.class1$Label,drugRWRfile.class1$norm_Score)  }else{ROC.list[[2]] = NA}
if(length(unique(drugRWRfile.class2$Label )) > 1){ROC.list[[3]] <- roc(drugRWRfile.class2$Label,drugRWRfile.class2$norm_Score)  }else{ROC.list[[3]] = NA}
if(length(unique(drugRWRfile.class3$Label )) > 1){ROC.list[[4]] <- roc(drugRWRfile.class3$Label,drugRWRfile.class3$norm_Score)  }else{ROC.list[[4]] = NA}
if(length(unique(drugRWRfile.class4$Label )) > 1){ROC.list[[5]] <- roc(drugRWRfile.class4$Label,drugRWRfile.class4$norm_Score)  }else{ROC.list[[5]] = NA}
if(length(unique(drugRWRfile.class5$Label )) > 1){ROC.list[[6]] <- roc(drugRWRfile.class5$Label,drugRWRfile.class5$norm_Score)  }else{ROC.list[[6]] = NA}
if(length(unique(drugRWRfile.class6$Label )) > 1){ROC.list[[7]] <- roc(drugRWRfile.class6$Label,drugRWRfile.class6$norm_Score)  }else{ROC.list[[7]] = NA}
names(ROC.list) = c('All','chemo_chemo','target_target','other_other','chemo_target','chemo_other','target_other')


pdf(file = paste0('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/L1_multiROC.pdf',sep = ''), width = 10, height = 10)

if(length(which(is.na(ROC.list))) != 0 ){ROC.list = ROC.list[-which(is.na(ROC.list))]}else{ROC.list = ROC.list}
try(plot(ROC.list[[1]], col='black',main = 'multiROC') ,silent = T)
try(plot.roc(ROC.list[[2]], add=TRUE, col='#4FFFFF')  ,silent = T)
try(plot.roc(ROC.list[[3]], add=TRUE, col='#FF00FF')  ,silent = T)
try(plot.roc(ROC.list[[4]], add=TRUE, col='#FF3333')  ,silent = T)
try(plot.roc(ROC.list[[5]], add=TRUE, col='#3381FF')  ,silent = T)
try(plot.roc(ROC.list[[6]], add=TRUE, col='#C6A300')  ,silent = T)
try(plot.roc(ROC.list[[7]], add=TRUE, col='#00BB00')  ,silent = T)

AUC.list = list()
for(Z in 1:length(ROC.list)){AUC.list[Z] = round(ROC.list[[Z]]$auc[1],3)}
names(AUC.list) = names(ROC.list)

if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','chemo_target')))==1 ){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[4]],sep='')),
         col = c('black','#4FFFFF','#FF00FF','#3381FF'))
}else if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','chemo_target','target_other')))==1){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[4]],sep=''),
                    paste('target_other_AUC=',AUC.list[[5]],sep='')
         ), col = c('black','#4FFFFF','#FF00FF', '#3381FF','#00BB00' ))
}else if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','chemo_target','chemo_other')))==1){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[4]],sep=''),
                    paste('chemo_other_AUC=',AUC.list[[5]],sep='')
         ), col = c('black','#4FFFFF','#FF00FF', '#3381FF','#C6A300'))
}else if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','other_other','chemo_target')))==1){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('other_other_AUC=',AUC.list[[4]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[5]],sep='')),
         col = c('black','#4FFFFF','#FF00FF','#FF3333','#3381FF'))
}else if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','chemo_target','chemo_other','target_other')))==1){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[4]],sep=''),
                    paste('chemo_other_AUC=',AUC.list[[5]],sep=''),
                    paste('target_other_AUC=',AUC.list[[6]],sep='')),
         col = c('black','#4FFFFF','#FF00FF', '#3381FF','#C6A300','#00BB00'))
}else if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','other_other','chemo_target','target_other')))==1){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('other_other_AUC=',AUC.list[[4]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[5]],sep=''),
                    paste('target_other_AUC=',AUC.list[[6]],sep='')),
         col = c('black','#4FFFFF','#FF00FF','#FF3333','#3381FF','#00BB00'))
}else if(length(unique(names(AUC.list) == c('All','chemo_chemo','target_target','other_other','chemo_target','chemo_other')))==1){
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('other_other_AUC=',AUC.list[[4]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[5]],sep=''),
                    paste('chemo_other_AUC=',AUC.list[[6]],sep='')),
         col = c('black','#4FFFFF','#FF00FF','#FF3333','#3381FF','#C6A300' ))
}else {
  legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('All_AUC=',AUC.list[[1]],sep=''),
                    paste('chemo_chemo_AUC=',AUC.list[[2]],sep=''),
                    paste('target_target_AUC=',AUC.list[[3]],sep=''),
                    paste('other_other_AUC=',AUC.list[[4]],sep=''),
                    paste('chemo_target_AUC=',AUC.list[[5]],sep=''),
                    paste('chemo_other_AUC=',AUC.list[[6]],sep=''),
                    paste('target_other_AUC=',AUC.list[[7]],sep='')),
         col = c('black','#4FFFFF','#FF00FF','#FF3333','#3381FF','#C6A300','#00BB00'))
}

dev.off()

#===============================================================================
#      Figure 3 L1 model TPR/TNR
#===============================================================================
DComboNet_modelroc <- roc(drugRWRfile$Label,drugRWRfile$norm_Score,
                          ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                          plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                          print.auc=TRUE, show.thres=TRUE, print.thres.best.method='closest.topleft')

DComboNet_measures = coords(DComboNet_modelroc, "best", best.method = 'closest.topleft', 
                            ret= 'all',transpose = FALSE)

rownames(DComboNet_measures) = 'all'
DComboNet_measures$auc =  DComboNet_modelroc$auc[1]

class1_modelroc <- roc(drugRWRfile.class1$Label,drugRWRfile.class1$norm_Score,
                          ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                          plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                          print.auc=TRUE, show.thres=TRUE, 
                          print.thres.best.method='closest.topleft')
class1_measures = coords(class1_modelroc, "best", best.method = 'closest.topleft', 
                            ret= 'all',transpose = FALSE)
rownames(class1_measures) = 'chemo_chemo'
class1_measures$auc =  class1_modelroc$auc[1]

class2_modelroc <- roc(drugRWRfile.class2$Label,drugRWRfile.class2$norm_Score,
                       ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE, 
                       print.thres.best.method='closest.topleft')
class2_measures = coords(class2_modelroc, "best", best.method = 'closest.topleft', 
                         ret= 'all',transpose = FALSE)
rownames(class2_measures) = 'target_target'
class2_measures$auc =  class2_modelroc$auc[1]

class3_modelroc <- roc(drugRWRfile.class3$Label,drugRWRfile.class3$norm_Score,
                       ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE, 
                       print.thres.best.method='closest.topleft')
class3_measures = coords(class3_modelroc, "best", best.method = 'closest.topleft', 
                         ret= 'all',transpose = FALSE)
rownames(class3_measures) = 'other_other'
class3_measures$auc =  class3_modelroc$auc[1]

class4_modelroc <- roc(drugRWRfile.class4$Label,drugRWRfile.class4$norm_Score,
                       ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE, 
                       print.thres.best.method='closest.topleft')
class4_measures = coords(class4_modelroc, "best", best.method = 'closest.topleft', 
                         ret= 'all',transpose = FALSE)
rownames(class4_measures) = 'chemo_target'
class4_measures$auc =  class4_modelroc$auc[1]

class5_modelroc <- roc(drugRWRfile.class5$Label,drugRWRfile.class5$norm_Score,
                       ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE, 
                       print.thres.best.method='closest.topleft')
class5_measures = coords(class5_modelroc, "best", best.method = 'closest.topleft', 
                         ret= 'all',transpose = FALSE)
rownames(class5_measures) = 'chemo_other'
class5_measures$auc =  class5_modelroc$auc[1]

class6_modelroc <- roc(drugRWRfile.class6$Label,drugRWRfile.class6$norm_Score,
                       ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                       plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                       print.auc=TRUE, show.thres=TRUE, 
                       print.thres.best.method='closest.topleft')
class6_measures = coords(class6_modelroc, "best", best.method = 'closest.topleft', 
                         ret= 'all',transpose = FALSE)
rownames(class6_measures) = 'target_other'
class6_measures$auc =  class6_modelroc$auc[1]


DComboNet_roc = data.frame(model = 'All', x=1-DComboNet_modelroc$specificities, y=DComboNet_modelroc$sensitivities)
class1_roc = data.frame(model = 'C_C', x=1-class1_modelroc$specificities, y=class1_modelroc$sensitivities)
class2_roc = data.frame(model = 'T_T', x=1-class2_modelroc$specificities, y=class2_modelroc$sensitivities)
class3_roc = data.frame(model = 'O_O', x=1-class3_modelroc$specificities, y=class3_modelroc$sensitivities)
class4_roc = data.frame(model = 'C_T', x=1-class4_modelroc$specificities, y=class4_modelroc$sensitivities)
class5_roc = data.frame(model = 'C_O', x=1-class5_modelroc$specificities, y=class5_modelroc$sensitivities)
class6_roc = data.frame(model = 'T_O', x=1-class6_modelroc$specificities, y=class6_modelroc$sensitivities)

dat = rbind(DComboNet_roc,
            class1_roc, 
            class2_roc, 
            class3_roc, 
            class4_roc, 
            class5_roc, 
            class6_roc)
dat$model = factor(dat$model, levels = unique(dat$model))


p = ggplot(dat) +
  geom_line(aes(x,y,colour = model),size = 0.3)  +
  labs(x="1 - Specificity", y="Sensitivity") +
  scale_colour_manual(breaks=unique(dat$model),
                      values = c('black','#4FFFFF','#FF00FF','#FF3333','#3381FF','#C6A300','#00BB00'),#unique(dat$model),
                      name="Different drug classes (AUC)",
                      labels=c(paste0('All (',round(DComboNet_modelroc$auc[1],3),')'),
                               paste0('C_C (',round(class1_modelroc$auc[1],3),')'),
                               paste0('T_T (',round(class2_modelroc$auc[1],3),')'),
                               paste0('O_O (',round(class3_modelroc$auc[1],3),')'),
                               paste0('C_T (',round(class4_modelroc$auc[1],3),')'),
                               paste0('C_O (',round(class5_modelroc$auc[1],3),')'),
                               paste0('T_O (',round(class6_modelroc$auc[1],3),')'))) +
  theme(legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        panel.grid.major=element_line(colour=NA)) +
  theme_set(theme_bw()) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed",colour='grey',alpha=0.8) +
  theme(legend.position=c(0.8,0.2))

p

ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/L1_multiROC2.pdf',p, width = 10, height = 10)

subclass_performance = rbind(DComboNet_measures, 
                             class1_measures,
                             class2_measures,
                             class3_measures,
                             class4_measures,
                             class5_measures,
                             class6_measures)
subclass_performance$subclasses = c("All", "C_C","T_T","O_O","C_T","C_O","T_O")
subclass_p_performance = subclass_performance[c('subclasses','tpr','fnr')]
df_p = reshape2::melt(subclass_p_performance)
df_p$subclasses = factor(df_p$subclasses,levels =c("All", "C_C","T_T","O_O","C_T","C_O","T_O"))
df_p$variable = factor(df_p$variable, levels = c('fnr','tpr'))
p = ggplot(df_p, aes(x = subclasses , y = value, fill = variable))+
  geom_bar(stat="identity",position="fill",width=0.8, color = 'black') +
  scale_fill_manual(values=c('#91bfdb','#fc8d59'))+ 
  geom_text(stat='identity',aes(label=paste0(round((1-value),4)*100,'%')), 
            color="black", size=7,
            position=position_dodge(0),vjust=-2)+
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+ 
    ylim(0, 1)+ xlab('')+ ylab('')
p
ggsave('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/positive_performance.pdf',plot=p, width=14, height=6)

df_p2 = df_p[df_p$variable == 'tpr',]
p1 = ggplot(df_p2, aes(x = subclasses, y = value))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') + #'#88CFF0'
  geom_text(aes(label=round(value,3)), vjust=-0.8, color="black", size=3.4)+
  ylim(0,1)+
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
p1
ggsave('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/positive_performance2.pdf',plot=p1, width=10, height=3)


subclass_n_performance = subclass_performance[c('subclasses','tnr','fpr')]
df_n = reshape2::melt(subclass_n_performance)
df_n$subclasses = factor(df_n$subclasses,levels = c("All", "C_C","T_T","O_O","C_T","C_O","T_O"))
df_n$variable = factor(df_n$variable, levels = c('fpr','tnr'))
n = ggplot(df_n, aes(x = subclasses , y = value, fill = variable))+
  geom_bar(stat="identity",position="fill",width=0.8, color = 'black') + 
  scale_fill_manual(values=c('#91bfdb', '#fc8d59'))+ 
  geom_text(stat='identity',aes(label=paste0(round((1-value),4)*100,'%')), 
            color="black", size=7,
            position=position_dodge(0),vjust=-2.5)+
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))+ 
  ylim(0, 1)+ xlab('')+ ylab('')
n
ggsave('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/negative_performance.pdf',plot=n, width=14, height=6)
df_n2 = df_n[df_n$variable == 'tnr',]
n1 = ggplot(df_n2, aes(x = subclasses, y = value))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') + #'#88CFF0'
  geom_text(aes(label=round(value,3)), vjust=-0.8, color="black", size=3.4)+
  ylim(0,1)+
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25))
n1
ggsave('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/negative_performance2.pdf',plot=n1, width=10, height=3)


pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/p_n_performance.pdf', width=12, height=11)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p, vp = vplayout(1,1))  
print(n, vp = vplayout(2,1))  
dev.off()

df2 = subclass_performance[c('subclasses','tpr','tnr')]
names(df2) = c("Classes", "TPR", "TNR")
df2$Classes = factor(df2$Classes,levels = c("All", "C_C","T_T","O_O","C_T","C_O","T_O"))
df2$color = c('black','#4FFFFF','#FF00FF','#FF3333','#3381FF','#C6A300','#00BB00')
df2p = ggplot(df2, aes(x = TPR , y = TNR, color = Classes))+
  geom_point(size=6) + #aes(color=color) + 
  scale_colour_manual(values=df2$color)+
  theme_set(theme_bw())+ 
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+ 
  xlim(0.5, 1) + ylim(0.5, 1)
df2p
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/p_n_performance2.pdf',df2p, width = 9.5, height = 6)
df2p = ggplot(df2, aes(x = TPR , y = TNR, color = Classes))+
  geom_point(size=6) + #aes(color=color) + 
  scale_colour_manual(values=df2$color)+
  theme_set(theme_bw())+ 
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
df2p
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/p_n_performance3.pdf',df2p, width = 9.5, height = 6)

#===============================================================================
#      Figure 3 L1 model comparison
#===============================================================================


# ROC curves see file roc_curve.R

setwd('G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/')
options(stringsAsFactors = F)
library(mlr)
library(pROC)
library(ggplot2)
library(ggthemes)
library(grid)
library(reshape2)

model_compairson <- read.csv('model_compairson_new.csv')

model_compairson$model = factor(rownames(model_compairson), levels = rownames(model_compairson))

p1 = ggplot(model_compairson, aes(x = X, y = auc))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') + #'#88CFF0'
  geom_text(aes(label=round(auc,3)), vjust=-0.8, color="black", size=3.4)+
  ylim(0,1) + theme(axis.text.x = element_text(angle = 0, hjust = 0.8, vjust = 0.5))
# p1

p2 = ggplot(model_compairson, aes(x = X, y = tpr))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') +
  geom_text(aes(label=round(tpr,3)), vjust=-0.8, color="black", size=3.4)+
  ylim(0,1) + theme(axis.text.x = element_text(angle = 0, hjust = 0.8, vjust = 0.5))
# p2
p3 = ggplot(model_compairson, aes(x = X, y = tnr))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') +
  geom_text(aes(label=round(tnr,3)), vjust=-0.8, color="black", size=3.4 )+
  ylim(0,1) + theme(axis.text.x = element_text(angle = 0, hjust = 0.8, vjust = 0.5))
# p3
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/results2_plot/auc.pdf',p1, width = 10, height = 3)
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/results2_plot/tpr.pdf',p2, width = 10, height = 3)
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/results2_plot/tnr.pdf',p3, width = 10, height = 3)



pdf('results2_plot/auc_tpr_res_new1_1.pdf', width=14, height=3.5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p1, vp = vplayout(1,1))  
print(p2, vp = vplayout(1,2))  
print(p3, vp = vplayout(1,3))  
dev.off()

pdf('results2_plot/auc_tpr_res_new1_2.pdf', width=6, height=12)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,1))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p1, vp = vplayout(1,1))  
print(p2, vp = vplayout(2,1))  
print(p3, vp = vplayout(3,1))  
dev.off()

model_compairson$model = factor(rownames(model_compairson), levels = rownames(model_compairson))
df <- model_compairson[c('model','auc', 'tpr','tnr')]
df = reshape2::melt(df)
df$color = '#fc8d59'
df[df$variable == 'tpr',]$color = '#ffffbf'
df[df$variable == 'tnr',]$color = '#91bfdb'

p = ggplot(df, aes(x = model, y = value, fill = variable))+
  geom_bar(stat = "identity", position = 'dodge',width=0.8, color = 'black') + 
  ylim(0,1)+
  scale_fill_manual(values=c('#fc8d59','#ffffbf','#91bfdb'))+ 
  geom_text(stat='identity',aes(label=round(value,3)), color="black", size=4,
            position=position_dodge(0.8),vjust=-0.8)
p
ggsave('results2_plot/auc_tpr_res_new2.pdf',plot=p, width=18, height=5)

df2 <- model_compairson[c('X', 'tpr','tnr')]
names(df2) = c("Models", "TPR", "TNR")
df2$color = c('#34AC37','#2F469A','#7BCADE','#BC4894','#C7C5C3','#F7EB00','#00B9E3','#FF0000','black')
df2$Models = factor(df2$Models,levels = df2$Models)
df2$color = factor(df2$color, levels = df2$color)
df2p = ggplot(df2, aes(x = TPR , y = TNR, color = Models))+
  geom_point(size=9) + #aes(color=color) + 
  scale_colour_manual(values=c('#34AC37','#2F469A','#7BCADE','#BC4894','#C7C5C3','#F7EB00','#00B9E3','#FF0000','black'))+
  theme_set(theme_bw())+ 
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+ 
  xlim(0, 1) + ylim(0, 1)
df2p
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/results2_plot/tpr_tnr.pdf',df2p, width = 9.5, height = 6)
df2p = ggplot(df2, aes(x = TPR , y = TNR, color = Models))+
  geom_point(size=9) + #aes(color=color) + 
  scale_colour_manual(values=c('#34AC37','#2F469A','#7BCADE','#BC4894','#C7C5C3','#F7EB00','#00B9E3','#FF0000','black'))+
  theme_set(theme_bw())+ 
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/results2_plot/tpr_tnr2.pdf',df2p, width = 9.5, height = 6)

#===============================================================================
#      Figure 3 result heatmap
#===============================================================================
setwd('C:/Users/fengf/Documents/Veronica_files/DCcomboNet/plot_update/for_heatmap/')
options(stringsAsFactors = F)
#library(ggplot2)
set.seed(123)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plyr)

# color = unique(c(brewer.pal(12,"Set3"),
#           brewer.pal(8,'Dark2'),
#           brewer.pal(8,'Accent'),
#           brewer.pal(8,'Pastel1'),
#           brewer.pal(9,"Set1")))
# color= c(brewer.pal(8,'Dark2'),brewer.pal(9,"Set3"))
# barplot(1:length(color),col=color)
# subclass_color = sample(unique(color)[-17] ,size = 14, replace = FALSE)
# barplot(1:length(subclass_color),col=subclass_color)
# 
# 
# 
# 
subclass_color = read.csv('heatmap_subclass_color.csv',stringsAsFactors = F)[,1]
#subclass_color = c(color[ seq(from=1, to=36, by=2) ],color[ seq(from=2, to=36, by=2) ])
drug_info = read.csv('drug_info2.csv',stringsAsFactors = F)
drugclass = drug_info[c('Drug','Class')]
drugsubclass = drug_info[c('Drug','Subclass')]
names(drugclass) = c('Drug','class')
names(drugsubclass) = c('Drug','subclass')

subclass_col = data.frame(subclass = unique(drugsubclass$subclass), 
                          col = subclass_color)
drugsubclass = merge(drugsubclass,subclass_col, by ='subclass')

all_rank = read.csv("result_all_score2.csv")
drugpair = read.csv("result_all.csv")
drugpair = drugpair[drugpair$Label == 1,][c('A','B')]
drugpair_index = c(paste0(drugpair$A, "_", drugpair$B), paste0(drugpair$B, "_", drugpair$A))
all_rank1 = all_rank[c("A","B","Rank_AplusB")]
all_rank1$index = paste0(all_rank1$A, "_", all_rank1$B)
all_rank2 = all_rank1[-which(is.element(all_rank1$index, drugpair_index)) , ]
all_rank2 = all_rank2[c("A","B","Rank_AplusB")]

names(all_rank2) = c("drugseeds","drugcandidates","Rank_ABplusBA")
tmp=reshape(all_rank2, idvar="drugseeds", timevar="drugcandidates", direction="wide")
tmp=tmp[order(tmp$drugseeds, decreasing = F),]

#tmp[is.na(tmp)]=0
rownames(tmp)=gsub('Rank_ABplusBA','',tmp$drugseeds)
colnames(tmp)=gsub('Rank_ABplusBA.','',colnames(tmp))
mat = as.matrix(tmp[,-1])
rownames(mat) = tmp[,1]
colnames(mat) = colnames(tmp)[-1]

palf = colorRampPalette(c('white','#436eee'))
colors = structure(c(rep("#ff0000",5), rep("#ffc4c4",10), rep('#ffffc4',15), palf(127)), names =1:157)

fh = function(x) fastcluster::hclust(dist(x))

col_mat = data.frame(Drug = colnames(mat))
row_mat = data.frame(Drug = rownames(mat))
col_mat = merge(drugsubclass,col_mat,by='Drug')
row_mat = merge(drugsubclass,row_mat,by='Drug')
col_mat = col_mat[order(col_mat$subclass,decreasing = T),]
row_mat = row_mat[order(row_mat$subclass,decreasing = T),]
mat = mat[row_mat$Drug,col_mat$Drug]
p1 = Heatmap(mat, 
             name = "mat", 
             na_col = "gray70", 
             col = colors,
             column_title = "Candidate Drugs", 
             row_title = "Drug Seeds",
             cluster_rows = F,
             cluster_columns  = F,
             show_heatmap_legend = F,
             row_names_gp = gpar(fontsize = 9),
             column_names_gp = gpar(fontsize = 9)) 
p1
#row_order(p1)
#column_order(p1)
mat1 = mat[row_order(p1),column_order(p1)]
#mat1 = mat[drugsubclass$Drug,drugsubclass$Drug]

anno_col = data.frame(Drug = colnames(mat1))
anno_col = join(anno_col,drugclass, type='left')
anno_col = join(anno_col,drugsubclass, type='left')
#anno_col[is.na(anno_col$class),]
#ann1 <- anno_col[c(2,3)]
col_class1 = data.frame(class = unique(drugclass$class))
col_class1$col = c("#5c5c5c",
                   "#a10303",
                   "#d1cfcf" )
col_class = col_class1$col
names(col_class) = unique(col_class1$class)
anno_col1 = unique(anno_col[c(3,4)])
col_subclass = anno_col1$col
names(col_subclass) = unique(anno_col1$subclass)
col_colours <- list(class=col_class,
                    subclass = col_subclass)
rownames(anno_col) <- anno_col$Drug
colAnn <- HeatmapAnnotation(df=anno_col[-c(1,4)],
                            which="col", 
                            col=col_colours, 
                            annotation_width=unit(c(1, 1), "cm"),
                            annotation_height = unit(10,'cm'),
                            gap=unit(1, "mm"), show_annotation_name = FALSE)


anno_row1 = data.frame(Drug = rownames(mat1))
rownames(anno_row1) = rownames(mat1)
anno_row1 = join(anno_row1,drugclass)
anno_row1 = join(anno_row1,drugsubclass)
#ann2 <- anno_row1[c(2,3)]
# row_class = c("Standard Chemotherapy"="#5c5c5c",
#               "Targeted Therapy"="#a10303",
#               "Other drug" = "#d1cfcf" )

row_class1 = data.frame(class = unique(drugclass$class))
row_class1$col = c("#5c5c5c",
                   "#a10303",
                   "#d1cfcf" )
row_class = row_class1$col
names(row_class) = row_class1$class

anno_row2= unique(anno_row1[c(3,4)])
row_subclass = anno_row2$col
names(row_subclass) = anno_row2$subclass
row_colours <- list(class=row_class,
                    subclass = row_subclass)
#ann2 <- anno_row1[c(1,2)]
rownames(anno_row1) <- anno_row1$Drug
row_ha = HeatmapAnnotation(df=anno_row1[-c(1,4)], 
                           which="row", 
                           col=row_colours, 
                           annotation_width=unit(c(0.5, 1), "cm"), 
                           annotation_height = unit(1,'cm'),
                           gap=unit(1, "mm"), 
                           show_annotation_name = FALSE)

p2 = Heatmap(mat1, 
             name = "mat", 
             na_col = "gray60", 
             col = colors,
             column_title = "Candidate Drugs", 
             row_title = "Drug Seeds",
             cluster_rows = F,
             cluster_columns  = F,
             bottom_annotation = colAnn,
             right_annotation = row_ha,
             show_heatmap_legend = F,
             row_names_gp = gpar(fontsize = 4),
             column_names_gp = gpar(fontsize = 5)
) 
# p2
pdf(paste0("heatmap_all_score2.pdf"), width = 18, height = 13)
draw(p2, 
     show_heatmap_legend = F, 
     show_annotation_legend = T)
dev.off()

library(pROC)
result_cv <- read.csv('result_all.csv')
a=0;
b=1;
y = result_cv$Score_ABplusBA
Ymax=max(result_cv$Score_ABplusBA);
Ymin=min(result_cv$Score_ABplusBA);
k=(b-a)/(Ymax-Ymin);
norY=a+k*(y-Ymin)
result_cv$norm_Score = norY
DComboNet_modelroc <- roc(result_cv$Label, result_cv$norm_Score,
                          ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                          plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                          print.auc=TRUE, show.thres=TRUE, print.thres.best.method='closest.topleft')
DComboNet_measures = coords(DComboNet_modelroc, "best", best.method = 'closest.topleft', ret= 'all',transpose = FALSE);rownames(DComboNet_measures) = 'DComboNet';DComboNet_measures$auc =  DComboNet_modelroc$auc[1]
result_cv_p = result_cv[result_cv$Label == 1,]
nrow(result_cv_p[result_cv_p$norm_Score > DComboNet_measures$threshold,])
drugseeds = as.list(union(result_cv_p$A, result_cv_p$B))
result_cv_Tp <- lapply(drugseeds, function(d){
  tmp = result_cv_p[result_cv_p$A == d | result_cv_p$B == d, ]
  result = data.frame(Drug = d,
                      num_p = nrow(tmp),
                      num_Tp = nrow(tmp[tmp$norm_Score > DComboNet_measures$threshold,]))
  return(result)
})
result_cv_Tp = do.call(rbind, result_cv_Tp)
result_cv_Tp$TPR_per_drugseed <- result_cv_Tp$num_Tp/result_cv_Tp$num_p
rownames(result_cv_Tp) = result_cv_Tp$drugseed

anno_row2 = merge(anno_row1, result_cv_Tp)
rownames(anno_row2) = anno_row2$Drug
row_ha2 = rowAnnotation(TPR = anno_barplot(anno_row2$TPR_per_drugseed, gp = gpar(fill = '#fcb09d'), 
                                         bar_width = 1, height = unit(6, "cm")))

p2 = Heatmap(mat1, 
             name = "mat", 
             na_col = "gray60", 
             col = colors,
             column_title = "Candidate Drugs", 
             row_title = "Drug Seeds",
             cluster_rows = F,
             cluster_columns  = F,
             bottom_annotation = colAnn,
             right_annotation = row_ha2,
             show_heatmap_legend = F,
             row_names_gp = gpar(fontsize = 4),
             column_names_gp = gpar(fontsize = 5)
)

p2
pdf(paste0("heatmap_all_score2_barplot.pdf"), width = 18, height = 13)
draw(p2, 
     show_heatmap_legend = F, 
     show_annotation_legend = T)
dev.off()

#===============================================================================
#      Figure 4 L2 model ROC curve
#===============================================================================
# HEPG2 L1-L2 ROC curve
library(pROC)
drugrank_L2 <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L2/HEPG2_result_all.csv')



drugrank_L1 <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/result_all.csv')
drugrank_L1 <- drugrank_L1[order(drugrank_L1$Score_ABplusBA, decreasing = T), ]
drugrank_L1$Rank_ABplusBA <- 1:nrow(drugrank_L1)
drugrank_L1 <- drugrank_L1[-9]
result_p2 = drugrank_L2[drugrank_L2$Label == 1, ]
result_p2 = result_p2[-9]

names(result_p2) = c('A','B',"Score_AB_L2","Rank_AB_L2","Score_BA_L2","Rank_BA_L2","Label","Score_ABplusBA_L2","Rank_ABplusBA_L2")
names(drugrank_L1) = c('A','B',"Score_AB_L1","Rank_AB_L1","Score_BA_L1","Rank_BA_L1","Label","Score_ABplusBA_L1","Rank_ABplusBA_L1")
tmp <- merge(result_p2, drugrank_L1, by = c("A", "B", "Label"))
tmp2 <- tmp[c("A", "B", "Score_AB_L1", "Rank_AB_L1", "Score_BA_L1", "Rank_BA_L1", "Label", "Score_ABplusBA_L1", "Rank_ABplusBA_L1")]
names(tmp2) <- names(drugrank_L1)
drugrank_L1 <- rbind(tmp2, drugrank_L1[drugrank_L1$Label == 0, ])



ROC.list = list()
if(length(unique(drugrank_L1$Label )) > 1){ ROC.list[[1]] <- roc(drugrank_L1$Label ,drugrank_L1$Score_ABplusBA_L1) }else{ROC.list[[1]] = NA}
if(length(unique(drugrank_L2$Label )) > 1){ ROC.list[[2]] <- roc(drugrank_L2$Label ,drugrank_L2$Score_ABplusBA) }else{ROC.list[[2]] = NA}

names(ROC.list) = c('drugrank_L1','drugrank_L2')


# roc.test( ROC.list[[1]],  ROC.list[[2]], method="venkatraman")
# roc.test( ROC.list[[1]],  ROC.list[[2]], method="delong")

cellline = 'HEPG2'

pdf(file = 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L2/multiROC_a_b_com_HEPG2_new.pdf')

if(length(which(is.na(ROC.list))) != 0 ){
  ROC.list = ROC.list[-which(is.na(ROC.list))]
}else{ROC.list = ROC.list}
try(plot(ROC.list[[1]], col='#4FFFFF',main = paste0(cellline,' Cell Line'),lty = 1,lwd =2) ,silent = T)
try(plot.roc(ROC.list[[2]], add=TRUE, col='#FD9089'  ),silent = T)
AUC.list = list()
for(Z in 1:length(ROC.list)){AUC.list[Z] = round(ROC.list[[Z]]$auc[1],3)}
names(AUC.list) = names(ROC.list)

legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
       legend = c(paste('drugrank_L1 =',AUC.list[[1]],sep=''),
                  paste('drugrank_L2 =',AUC.list[[2]],sep='')),
       col = c('#4FFFFF','#FD9089'))

dev.off()

p_L1 <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L2/drugpair_positive_L1.csv')
p_L2 <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L2/drugpair_positive.csv')
p_L1 <- p_L1[-4]
names(p_L1) <- c("A", "B", "rank_position_L1")
names(p_L2) <- c("A", "B", "rank_position_L2", "max_rank")

res = merge(p_L1, p_L2, by = c('A','B'))
res$direction = res$rank_position_L1 - res$rank_position_L2
res$color = 'orange'
res[res$direction < 0,]$color = 'grey'
res = res[order(res$color),]
res[res$A == 'Capecitabine' & res$B == "Docetaxel", ]$color = 'red'

drugrank_L2_modelroc <- roc(drugrank_L2$Label, drugrank_L2$Score_ABplusBA,
                            ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                            print.auc=TRUE, show.thres=TRUE, print.thres.best.method='closest.topleft')
drugrank_L2_measures = coords(drugrank_L2_modelroc, "best", best.method = 'closest.topleft', ret= 'all',transpose = FALSE);rownames(DComboNet_measures) = 'DComboNet';DComboNet_measures$auc =  DComboNet_modelroc$auc[1]
drugrank_L1_modelroc <- roc(drugrank_L1$Label, drugrank_L1$Score_ABplusBA,
                            ci=TRUE, boot.n=1000, ci.alpha=0.9, stratified=FALSE,
                            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                            print.auc=TRUE, show.thres=TRUE, print.thres.best.method='closest.topleft')
drugrank_L1_measures = coords(drugrank_L1_modelroc, "best", best.method = 'closest.topleft', ret= 'all',transpose = FALSE);rownames(DComboNet_measures) = 'DComboNet';DComboNet_measures$auc =  DComboNet_modelroc$auc[1]

pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L2/HEPG2_L1_L2_p_compr_new.pdf', width = 3, height = 5)
plot(x = c(1,2),
     y = res[1,c(3,4)],
     xlim = c(0.5,2.5), ylim = c(0,10+max(c(res$rank_position_L1, res$rank_position_L2))),
     col = res[1, which(names(res)=='color')], type = 'o', pch=20,
     xlab = "Models", ylab = "Ranks")
for(i in 2:nrow(res)){
  lines(x = c(1,2),
        y = res[i,c(3,4)],
        type = 'o', pch=20, col = res[i,which(names(res)=='color')])}
# abline(h=length(which(drugrank_L2$Score_ABplusBA >= drugrank_L2_measures$threshold))/nrow(drugrank_L2)*150, col = "green", lty = 2)
# abline(h=length(which(drugrank_L1$Score_ABplusBA >= drugrank_L1_measures$threshold))/nrow(drugrank_L1)*150, col = "gray50", lty = 2)
# for plot 'HEPG2_L1_L2_p_compr3.pdf'
abline(h=length(which(drugrank_L2$Score_ABplusBA >= drugrank_L2_measures$threshold))/nrow(drugrank_L2)*150, col = "gray50", lty = 2)
dev.off()


#===============================================================================
#      Figure 4 L2 model comparison
#===============================================================================
# AUC, TPR, TNR barplots
library(ggplot2)
library(grid)
method_compr_summary <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/method_compr_overlap_summary.csv')
method_compr_summary <- read.csv('C:/Users/fengf/Documents/Veronica_files/DCcomboNet/plot_update/newplot/fig4c_d//method_compr_overlap_summary.csv')


tpr_df <- method_compr_summary[c('methods', 'TPR')]
tpr_df$FNR <- 1-tpr_df$TPR
tpr_df = reshape2::melt(tpr_df)
tpr_df$methods = factor(tpr_df$methods, levels = method_compr_summary$methods)
tpr_df$variable = factor(tpr_df$variable, levels = c('FNR','TPR'))
tpr_df = tpr_df[order(tpr_df$variable, decreasing = T),]
p1 = ggplot(tpr_df, aes(x = methods , y = value, fill=variable ))+
  geom_bar(stat="identity",position="fill",width=0.8, color = 'black') +
  scale_fill_manual(values=c('#91bfdb','#fc8d59'))+ 
    geom_text(stat='identity',aes(label=paste0((1-value)*100,'%')), 
            color="black", size=7,
            position=position_dodge(0),vjust=0)+
  theme(axis.text.x = element_text(size = 18,angle = 90), 
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  theme_set(theme_bw())+ 
  ylim(0, 1)+ xlab('')+ ylab('')
p1

tnr_df <- method_compr_summary[c('methods', 'TNR')]
tnr_df$FPR <- 1-tnr_df$TNR
tnr_df = reshape2::melt(tnr_df)
tnr_df$methods = factor(tnr_df$methods, levels = method_compr_summary$methods)
tnr_df$variable = factor(tnr_df$variable, levels = c('FPR','TNR'))
tnr_df = tnr_df[order(tnr_df$variable, decreasing = T),]
p2 = ggplot(tnr_df, aes(x = methods , y = value, fill=variable ))+
  geom_bar(stat="identity",position="fill",width=0.8, color = 'black') +
  scale_fill_manual(values=c('#91bfdb','#fc8d59'))+ 
  geom_text(stat='identity',aes(label=paste0((1-value)*100,'%')), 
            color="black", size=7,
            position=position_dodge(0),vjust=0) +
  theme(axis.text.x = element_text(size = 18,angle = 90), 
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  theme_set(theme_bw())+ 
  ylim(0, 1)+ xlab('')+ ylab('')
p2

AUC_df <- method_compr_summary[c('methods', 'AUC')]
AUC_df$methods = factor(AUC_df$methods, levels = method_compr_summary$methods)
p3 = ggplot(AUC_df, aes(x = methods , y = AUC, fill = '#91bfdb'))+
  geom_bar(stat="identity",width=0.8, color = 'black') +
  geom_text(stat='identity',aes(label=AUC), 
            color="black", size=7,
            position=position_dodge(0),vjust=0)+
  theme(axis.text.x = element_text(size = 18,angle = -90), 
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  theme_set(theme_bw())+ 
  ylim(0, 1)+ xlab('')+ ylab('')
p3
pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/TPR_TNR_barplot.pdf', width=8, height=21)
pdf('C:/Users/fengf/Documents/Veronica_files/DCcomboNet/plot_update/newplot/fig4c_d/TPR_TNR_barplot.pdf', width=8, height=15)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3,1))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p1, vp = vplayout(1,1))  
print(p2, vp = vplayout(2,1))  
print(p3, vp = vplayout(3,1))  
dev.off()

tpr_df = tpr_df[tpr_df$variable == 'TPR',]

tpr_p = ggplot(tpr_df, aes(x = methods , y = value))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') +
  geom_text(aes(label=round(value,3)), vjust=-0.8, color="black", size=3.4)+
  ylim(0,1) + theme(axis.text.x = element_text(angle = 0, hjust = 0.8, vjust = 0.5))
# p2
tnr_df = tnr_df[tnr_df$variable == 'TNR',]
tnr_p = ggplot(tnr_df, aes(x = methods, y = value))+
  geom_bar(stat = "identity", fill='#fee0d2',width=0.8, color='black') +
  geom_text(aes(label=round(value,3)), vjust=-0.8, color="black", size=3.4 )+
  ylim(0,1) + theme(axis.text.x = element_text(angle = 0, hjust = 0.8, vjust = 0.5))
# p3
pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/TPR_TNR_barplot2.pdf', width=8, height=14)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1))) 
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(tpr_p, vp = vplayout(1,1))  
print(tnr_p, vp = vplayout(2,1))  
dev.off()

tpr_tnr = method_compr_summary[c('methods', 'TPR','TNR')]
tpr_tnr_p = ggplot(tpr_tnr, aes(x = TPR , y = TNR, color = methods))+
  geom_point(size=9) + #aes(color=color) + 
  scale_colour_manual(values=c('#34AC37','#2F469A','#7BCADE','#BC4894','#C7C5C3','#F7EB00','#00B9E3','#FF0000','black'))+
  theme_set(theme_bw())+ 
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25),
        panel.grid.major=element_line(colour=NA),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))+ 
  xlim(0.3, 0.8) + ylim(0.7, 0.9)
tpr_tnr_p
ggsave(file='G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/TPR_TNR_dotplot.pdf',tpr_tnr_p, width = 9, height = 3)


# ROC curves
method_compr <- read.csv('C:/Users/fengf/Documents/Veronica_files/DCcomboNet/plot_update/newplot/fig4c_d/method_compr_overlap_confirm.csv')
ROC.list = list()
ROC.list[[1]] <- roc(method_compr$label,1/method_compr$SynGen) 
ROC.list[[2]] <- roc(method_compr$label,1/method_compr$DIGRE)  
ROC.list[[3]] <- roc(method_compr$label,1/method_compr$DrugComboRanker)
ROC.list[[4]] <- roc(method_compr$label,1/method_compr$Zhao.s.method)
ROC.list[[5]] <- roc(method_compr$label,1/method_compr$RACS)
ROC.list[[6]] <- roc(method_compr$label,1/method_compr$DComboNet)
names(ROC.list) = c('SynGen','DIGRE','DrugComboRanker','Zhao.s.method','RACS','DComboNet')


pdf(file = paste0('C:/Users/fengf/Documents/Veronica_files/DCcomboNet/plot_update/newplot/fig4c_d/method_compr_ROC.pdf',sep = ''), width = 10, height = 10)

if(length(which(is.na(ROC.list))) != 0 ){ROC.list = ROC.list[-which(is.na(ROC.list))]}else{ROC.list = ROC.list}
try(plot(ROC.list[[1]], col='#C6A300',main = 'multiROC') ,silent = T)
try(plot.roc(ROC.list[[2]], add=TRUE, col='#4FFFFF')  ,silent = T)
try(plot.roc(ROC.list[[3]], add=TRUE, col='#FF00FF')  ,silent = T)
try(plot.roc(ROC.list[[4]], add=TRUE, col='#FF3333')  ,silent = T)
try(plot.roc(ROC.list[[5]], add=TRUE, col='#3381FF')  ,silent = T)
try(plot.roc(ROC.list[[6]], add=TRUE, col='black')  ,silent = T)

AUC.list = list()
for(Z in 1:length(ROC.list)){AUC.list[Z] = round(ROC.list[[Z]]$auc[1],3)}
names(AUC.list) = names(ROC.list)

legend('bottomright', inset = .01, lty = c(1,1), cex = 0.9,
         legend = c(paste('SynGen =',AUC.list[[1]],sep=''),
                    paste('DIGRE =',AUC.list[[2]],sep=''),
                    paste('DrugComboRanker =',AUC.list[[3]],sep=''),
                    paste('Zhao.s.method =',AUC.list[[4]],sep=''),
                    paste('RACS =',AUC.list[[5]],sep=''),
                    paste('DComboNet =',AUC.list[[6]],sep='')),
         col = c('#C6A300','#4FFFFF','#FF00FF','#FF3333','#3381FF','black' ))

dev.off()


#===============================================================================
#      Figure 5 Case study - Sorafenib Heatmap
#===============================================================================
setwd('G:/lab/Projects/p1_DComboNet/DComboNet_improve/Sorafenib/HEPG2_L2_casestudy_result/')
options(stringsAsFactors = F)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(plyr)

sorafenib_rank <- read.csv('sorafenib_rank.csv')
names(sorafenib_rank) <- c('Drugseed', 'Drug',names(sorafenib_rank)[-c(1,2)])
subclass_color = read.csv('heatmap_subclass_color.csv',stringsAsFactors = F)[,1]
#subclass_color = c(color[ seq(from=1, to=36, by=2) ],color[ seq(from=2, to=36, by=2) ])
drug_info = read.csv('drug_info2.csv',stringsAsFactors = F)
drugclass = drug_info[c('Drug','Class')]
drugsubclass = drug_info[c('Drug','Subclass')]
names(drugclass) = c('Drug','class')
drugclass$class_color = "#5c5c5c"
drugclass[drugclass$class == "Targeted Therapy",]$class_color <- "#a10303"
drugclass[drugclass$class == "Other drug",]$class_color <- "#d1cfcf"

names(drugsubclass) = c('Drug','subclass')
subclass_col = data.frame(subclass = unique(drugsubclass$subclass), 
                          subclass_color = subclass_color)
drugsubclass = merge(drugsubclass,subclass_col, by ='subclass')
sorafenib_rank = merge(sorafenib_rank, drugclass, by = 'Drug')
sorafenib_rank = merge(sorafenib_rank, drugsubclass, by = 'Drug')

sorafenib_rank = sorafenib_rank[order(sorafenib_rank$Rank_ABplusBA),]
sorafenib_pred = sorafenib_rank[sorafenib_rank$Rank_ABplusBA <= 0.2*nrow(sorafenib_rank),]


old_result <- read.csv('old_result.csv')
old_result <- old_result[c(1,7:11)]
names(old_result) <- c("Drug",names(old_result)[-1])
sorafenib_pred2 <- merge(sorafenib_pred, old_result, by='Drug', all.x = T)
sorafenib_pred2 <- sorafenib_pred2[order(sorafenib_pred2$Rank_ABplusBA),]
write.csv(sorafenib_pred2,'sorafenib_pred.csv', quote=F, row.names = F)

sorafenib_pred <- read.csv('sorafenib_pred.csv')
sorafenib_pred1 <- sorafenib_pred[sorafenib_pred$Known_combinations == 1, ]
sorafenib_pred2 <- sorafenib_pred[sorafenib_pred$Known_combinations == 0, ]
sorafenib_pred1 <- sorafenib_pred1[order(sorafenib_pred1$Rank_ABplusBA),]
sorafenib_pred2 <- sorafenib_pred2[order(sorafenib_pred2$Rank_ABplusBA),]
sorafenib_pred2_1 <- sorafenib_pred2[which(sorafenib_pred2$ref_HCC == 1), ]
sorafenib_pred2_2 <- sorafenib_pred2[which(is.na(sorafenib_pred2$ref_HCC) & sorafenib_pred2$ref_othercancer == 1), ]
sorafenib_pred2_3 <- sorafenib_pred2[which(is.na(sorafenib_pred2$ref_HCC)& is.na(sorafenib_pred2$ref_othercancer)), ]
sorafenib_pred2_1 <- sorafenib_pred2_1[order(sorafenib_pred2_1$Rank_ABplusBA),]
sorafenib_pred2_2 <- sorafenib_pred2_2[order(sorafenib_pred2_2$Rank_ABplusBA),]
sorafenib_pred2_3 <- sorafenib_pred2_3[order(sorafenib_pred2_3$Rank_ABplusBA),]
sorafenib_pred <- rbind(sorafenib_pred1,
                        sorafenib_pred2_1,
                        sorafenib_pred2_2,
                        sorafenib_pred2_3)

mat <- as.matrix(sorafenib_pred[c(12,13)])
colnames(mat) <- c("HCC", "Other cancer") 
rownames(mat) <- sorafenib_pred$Drug
mat <- t(mat)

annot_df <- sorafenib_pred[c(7,8,10)]
names(annot_df) <- c("If Known", "Class", "Subclass")
rownames(annot_df) <- sorafenib_pred$Drug
ann_colors <- list(
  `If Known` = c("#8BD1A0", "#EFA29A" ),
  Class = c(
    `Standard Chemotherapy` = "#474647",
    `Targeted Therapy` = "#a10303",
    `Other drug` = "#d1cfcf"
  ),
  Subclass =c(
    `Alkylating agents` = "#80B1D3",
    `Antibiotics/antineoplastics` = "#B3DE69",
    `Antimetabolites` = "#1B9E77",
    `TKIs & RTKIs` = "#666666",
    `Others` = "#FDB462",
    `Histone Deacetylase Inhibitor` = "#7570B3",
    `MTOR inhibitors` = "#FCCDE5",
    `Mitotic inhibitors` = "#BEBADA",
    `Immune` = "#80B1D3",
    `OtherSymptoms` = "#8DD3C7"

  )
)

pdf("plot_with_legends.pdf", width = 15, height = 5)
pheatmap(mat,
         cluster_rows = F, 
         cluster_cols = F,
         breaks = c(0,1),
         color = "#828FC7", 
         na_col = "white", 
         border_color = "black",
         annotation_col = annot_df,
         annotation_colors = ann_colors,
         legend = F) 
dev.off()

pdf("plot_without_legends.pdf", width = 15, height = 2.5)
pheatmap(mat,
         cluster_rows = F, 
         cluster_cols = F,
         breaks = c(0,1),
         color = "#828FC7", 
         na_col = "white", 
         border_color = "black",
         annotation_col = annot_df,
         annotation_colors = ann_colors,
         annotation_legend = F,
         legend = F) 
dev.off()



sorafenib_rank <- read.csv('sorafenib_rank.csv')
names(sorafenib_rank) <- c('Drugseed', 'Drug',names(sorafenib_rank)[-c(1,2)])
subclass_color = read.csv('heatmap_subclass_color.csv',stringsAsFactors = F)[,1]
#subclass_color = c(color[ seq(from=1, to=36, by=2) ],color[ seq(from=2, to=36, by=2) ])
drug_info = read.csv('drug_info2.csv',stringsAsFactors = F)
drugclass = drug_info[c('Drug','Class')]
drugsubclass = drug_info[c('Drug','Subclass2')]
names(drugclass) = c('Drug','class')
drugclass$class_color = "#5c5c5c"
drugclass[drugclass$class == "Targeted Therapy",]$class_color <- "#a10303"
drugclass[drugclass$class == "Other drug",]$class_color <- "#d1cfcf"

names(drugsubclass) = c('Drug','subclass')
sorafenib_rank = merge(sorafenib_rank, drugclass, by = 'Drug')
sorafenib_rank = merge(sorafenib_rank, drugsubclass, by = 'Drug')

sorafenib_rank = sorafenib_rank[order(sorafenib_rank$Rank_ABplusBA),]


old_result <- read.csv('sorafenib_pred.csv')
old_result <- old_result[c(1,12:16)]
names(old_result) <- c("Drug",names(old_result)[-1])
sorafenib_pred <- merge(sorafenib_rank, old_result, by='Drug', all.x = T)
sorafenib_pred <- sorafenib_pred[order(sorafenib_pred$Rank_ABplusBA),]
# write.csv(sorafenib_pred2,'sorafenib_pred.csv', quote=F, row.names = F)
sorafenib_pred1 = sorafenib_pred[sorafenib_pred$Known_combinations == 1, ]
sorafenib_pred2 = sorafenib_pred[sorafenib_pred$Rank_ABplusBA <= 0.2*nrow(sorafenib_pred),]
sorafenib_pred <- unique(rbind(sorafenib_pred1, sorafenib_pred2))

sorafenib_pred1 <- sorafenib_pred[sorafenib_pred$Known_combinations == 1, ]
sorafenib_pred2 <- sorafenib_pred[sorafenib_pred$Known_combinations == 0, ]
sorafenib_pred1 <- sorafenib_pred1[order(sorafenib_pred1$Rank_ABplusBA),]
sorafenib_pred2 <- sorafenib_pred2[order(sorafenib_pred2$Rank_ABplusBA),]
sorafenib_pred2_1 <- sorafenib_pred2[which(sorafenib_pred2$ref_HCC == 1), ]
sorafenib_pred2_2 <- sorafenib_pred2[which(is.na(sorafenib_pred2$ref_HCC) & sorafenib_pred2$ref_othercancer == 1), ]
sorafenib_pred2_3 <- sorafenib_pred2[which(is.na(sorafenib_pred2$ref_HCC)& is.na(sorafenib_pred2$ref_othercancer)), ]
sorafenib_pred2_1 <- sorafenib_pred2_1[order(sorafenib_pred2_1$Rank_ABplusBA),]
sorafenib_pred2_2 <- sorafenib_pred2_2[order(sorafenib_pred2_2$Rank_ABplusBA),]
sorafenib_pred2_3 <- sorafenib_pred2_3[order(sorafenib_pred2_3$Rank_ABplusBA),]
sorafenib_pred <- rbind(sorafenib_pred1,
                        sorafenib_pred2_1,
                        sorafenib_pred2_2,
                        sorafenib_pred2_3)
sorafenib_pred[is.element(sorafenib_pred$Drug, c("Afatinib","Cetuximab","Erlotinib","Lapatinib")),]$subclass <- "EGFR inhibitor"

mat <- as.matrix(sorafenib_pred[c(11,12)])
colnames(mat) <- c("HCC", "Other cancer") 
rownames(mat) <- sorafenib_pred$Drug
mat <- t(mat)

annot_df <- sorafenib_pred[c(7,8,10)]
names(annot_df) <- c("If Known", "Class", "Subclass")
rownames(annot_df) <- sorafenib_pred$Drug
ann_colors <- list(
  `If Known` = c("#8BD1A0", "#EFA29A" ),
  Class = c(
    `Standard Chemotherapy` = "#474647",
    `Targeted Therapy` = "#a10303",
    `Other drug` = "#d1cfcf"
  ),
  Subclass =c(
    `DNAorRNA` = "#C45B14",
    `Microtubule &  Mitotic` = "#E1E1DE",
    `Histone Deacetylase Inhibitor` = "#F7ED16",
    `protein kinases inhibitor` = "#9DC5B4",
    `Proteosome inhibitor` = "#009A70",
    `Tyrosine kinase & receptor` = "#FAE1EC",
    `EGFR inhibitor` = "#EC6600",
    `mTOR kinase inhibitor` = "#E6B200",
    `PARP inhibitor` = "#2DB0CB",
    `BCL inhibitor` = "#978FC4",
    `CDK Inhibitor` = "#48B540",
    `MABs` = "#80B1D3",
    `Retinoic acid receptor` = "#793EA3",
    `OtherSymptoms` = "#F2A0E1",
    `Unknown` = "#878787"
  )
)

pdf("plot_with_legends2.pdf", width = 15, height = 5)
pheatmap(mat,
         cluster_rows = F, 
         cluster_cols = F,
         breaks = c(0,1),
         color = "#828FC7", 
         na_col = "white", 
         border_color = "black",
         annotation_col = annot_df,
         annotation_colors = ann_colors,
         legend = F) 
dev.off()

pdf("plot_without_legends2.pdf", width = 15, height = 2.5)
pheatmap(mat,
         cluster_rows = F, 
         cluster_cols = F,
         breaks = c(0,1),
         color = "#828FC7", 
         na_col = "white", 
         border_color = "black",
         annotation_col = annot_df,
         annotation_colors = ann_colors,
         annotation_legend = F,
         legend = F) 
dev.off()

#===============================================================================
#      Supplementary Figure 
#===============================================================================


setwd('G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/parameter_tuning/')
options(stringsAsFactors = F)
library(ggplot2)
net_tuning <- read.csv("net_tuning.csv", stringsAsFactors = F)
net_tuning$net <- factor(net_tuning$net, levels = net_tuning$net)
net_tuning_p = ggplot(net_tuning, aes(x=net, y=AUC)) + geom_line()+ geom_point(size=8, shape=20)+
  geom_line()+
  geom_point(size=8, shape=20) +
  geom_text(aes(label=AUC), size=6,nudge_y = -0.1,angle=90)+
  # scale_x_continuous(breaks = c(0,seq(0.1,0.9,0.1),1))+
  scale_y_continuous(limits = seq(0,1),breaks=seq(0,1,0.1))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1,size=18),
        axis.text.y=element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+
  geom_point(aes(x = "DComboNet", y = net_tuning[net_tuning$net == "DComboNet", ]$AUC),color='red',size=8)+ 
  labs(x = "Network", y = "AUC")
net_tuning_p
ggsave("net_tuning.pdf", net_tuning_p, width = 8, height = 8)


tuning_r <- read.csv("tuning_r.csv")
tuning.r = ggplot(tuning_r, aes(x=r, y=AUC)) + geom_line()+ geom_point(size=8, shape=20)+
  geom_line()+
  geom_point(size=8, shape=20) +
  geom_text(aes(label=AUC), size=6,nudge_y = -0.1,angle=90)+
  geom_vline(xintercept=c(0.7), linetype="dotted")+
  scale_x_continuous(breaks = c(0,seq(0.1,0.9,0.1),1))+
  scale_y_continuous(limits = seq(0,1),breaks=seq(0,1,0.1))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1,size=18),
        axis.text.y=element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+
  labs(x = "sigma", y = "AUC")+
  geom_point(aes(x=0.7, y=tuning_r[tuning_r$r ==0.7,]$AUC),color='red',size=8)
# tuning.r
ggsave("tuning_r.pdf", tuning.r, width = 8, height = 8)

tuning_A <- read.csv("tuning_A.csv")

tuning.A = ggplot(tuning_A, aes(x=A, y=AUC)) + 
  geom_line()+ 
  geom_point(size=8, shape=20)+
  geom_text(aes(label=AUC), size=6,nudge_y = -0.1,angle=90)+
  geom_vline(xintercept=c(0.5), linetype="dotted")+
  scale_x_continuous(breaks = c(0,seq(0.1,0.9,0.1),1))+
  scale_y_continuous(limits = seq(0,1),breaks=seq(0,1,0.1))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1,size=18),
        axis.text.y=element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+
  geom_point(aes(x=0.5, y=max(tuning_A$AUC)),color='red',size=8)+
  labs(x = "lambda_DD", y = "AUC")
# tuning.A
ggsave("tuning_A.pdf",tuning.A, width=8, height=8)

L1_random <- read.csv("G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/randomization/rand_whole_net_result_L1.csv")
L1_random_p <- ggplot(L1_random, aes(x = auc))+
  geom_density(alpha=.2, fill="grey50")+
  geom_vline(xintercept=0.825, color = 'red', size = 2)+
  scale_x_continuous(limits = seq(0,1), breaks = seq(0,1,0.1))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1,size=18),
        axis.text.y=element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))+
  labs(x = "AUC", y = "density")+
  geom_text(aes(x = 0.9, y=1, label = paste("AUC_T = 0.825\nP=", round((length(which(L1_random$auc>0.825)))/(nrow(L1_random)),3))))
L1_random_p
ggsave("G:/lab/Projects/p1_DComboNet/DComboNet_improve/LOOCV_L1/randomization/L1_random.pdf",L1_random_p, width=8, height=7)

