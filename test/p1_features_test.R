library(ROCR)
features <- read.csv('G:/lab/Projects/p1_DComboNet/plot_V3/features.csv')
features$label = 0
features[features$Label =='P',]$label = 1
pred <- prediction( features$Pathway.Similarity, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

pred <- prediction( features$Target.Distance, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

pred <- prediction( features$sim_ATC, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

pred <- prediction( features$sim_structure, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

pred <- prediction( features$sim_sideeffects, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

wilcox.test(features[features$Label == 'P',]$sim_ATC, features[features$Label == 'N',]$sim_ATC)
wilcox.test(features[features$Label == 'P',]$sim_structure, features[features$Label == 'N',]$sim_structure)
wilcox.test(features[features$Label == 'P',]$sim_sideeffects, features[features$Label == 'N',]$sim_sideeffects)
wilcox.test(features[features$Label == 'P',]$Integrated.Score, features[features$Label == 'N',]$Integrated.Score)
wilcox.test(features[features$Label == 'P',]$Target.Distance, features[features$Label == 'N',]$Target.Distance)
wilcox.test(features[features$Label == 'P',]$Pathway.Similarity, features[features$Label == 'N',]$Pathway.Similarity)



features <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/mlr_models/data/features2_HEPG2_2.csv')
features[features$TAG == 0,]$TAG = 'P'
features$label = 0
features[features$TAG =='P',]$label = 1
pred <- prediction( features$STsim, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

pred <- prediction( features$ATCsim, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue
pred <- prediction( features$SEsim, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue
pred <- prediction( features$integrated_score2, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

features$integrated_score4 = 1-(1-features$ATCsim)*(features$STsim)*(1-features$SEsim)
pred <- prediction( features$integrated_score4, features$label)
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AUCValue

features2 <- read.csv('G:/lab/Projects/p1_DComboNet/plot_V3/features.csv')

features1 = merge(features, features2, by=c('A','B'))
names(features) = c('B','A', names(features)[-c(1,2)])
features1.2 = merge(features, features2, by=c('A','B'))
features = rbind(features1, features1.2)

wilcox.test(features[features$TAG == 'P',]$ATCsim, features[features$TAG == 'N',]$ATCsim, alternative = 'greater')
wilcox.test(features[features$TAG == 'P',]$STsim, features[features$TAG == 'N',]$STsim, alternative = 'greater')
wilcox.test(features[features$TAG == 'P',]$SEsim , features[features$TAG == 'N',]$SEsim , alternative = 'greater')
wilcox.test(features[features$TAG == 'P',]$integrated_score2, features[features$TAG == 'N',]$integrated_score2, alternative = 'greater')
wilcox.test(features[features$TAG == 'P',]$Target.Distance, features[features$TAG == 'N',]$Target.Distance, alternative = 'less')
wilcox.test(features[features$TAG == 'P',]$Pathway.Similarity, features[features$TAG == 'N',]$Pathway.Similarity, alternative = 'greater')



features[is.na(features$ATCsim),]$ATCsim = mean(features[!is.na(features$ATCsim),]$ATCsim)
features[is.na(features$TS),]$TS = mean(features[!is.na(features$TS),]$TS)
features[is.na(features$SE_TS),]$SE_TS = mean(features[!is.na(features$SE_TS),]$SE_TS)
features[is.na(features$WWI),]$WWI = mean(features[!is.na(features$WWI),]$WWI)
features[is.na(features$distance),]$distance = mean(features[!is.na(features$distance),]$distance)



set.seed(1234)
res = data.frame(times = 1:1000, ATCsim_p = 0, STsim_p = 0, SEsim_p = 0, Integrated_p = 0, distance_p = 0, WWI_p = 0)
for(i in 1:1000){
  n_index = sample(which(features$TAG == 'N'), size=nrow(features[features$TAG == 'P',]))
  res[res$times == i,]$ATCsim_p = wilcox.test(features[features$TAG == 'P',]$ATCsim, features[n_index, ]$ATCsim, alternative = 'greater')$p.value
  res[res$times == i,]$STsim_p = wilcox.test(features[features$TAG == 'P',]$STsim, features[n_index, ]$STsim, alternative = 'greater')$p.value
  res[res$times == i,]$SEsim_p = wilcox.test(features[features$TAG == 'P',]$SEsim , features[n_index, ]$SEsim , alternative = 'greater')$p.value
  res[res$times == i,]$Integrated_p = wilcox.test(features[features$TAG == 'P',]$Integrated.Score, features[n_index, ]$Integrated.Score)$p.value
  res[res$times == i,]$distance_p = wilcox.test(features[features$TAG == 'P',]$Target.Distance, features[n_index, ]$Target.Distance, alternative = 'less')$p.value
  res[res$times == i,]$WWI_p = wilcox.test(features[features$TAG == 'P',]$Pathway.Similarity, features[n_index, ]$Pathway.Similarity, alternative = 'greater')$p.value
}


features$Label = factor(features$Label, levels = c('P','N'))
features[features$Target.Distance==1e-50,]$Target.Distance=NA
features[features$Pathway.Similarity==1e-50,]$Pathway.Similarity=NA
features[features$integrated_score2==1e-50,]$integrated_score2=NA
features[features$sim_ATC==1e-50,]$sim_ATC=NA
features[features$sim_structure==1e-50,]$sim_structure=NA
features[features$sim_sideeffects==1e-50,]$sim_sideeffects=NA

data1 = melt(features[c(1:5,7,8,16,17)])
p1 <- ggplot(data=data1, aes(x=TAG,y=value))+
  geom_boxplot(aes(fill=TAG))+
  geom_signif(comparisons = list(c('P','N')), test = wilcox.test, textsize = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(size=0.58),
        axis.text.y = element_text(size=0.55))

p1 = p1 + facet_wrap(~ variable, scales="free") +
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(size = 16))
p1

features <- read.csv('/picb/bigdata/project/FengFYM/DComboNet_improve/mlr_models/data/features4.csv')
features[features$Target.Distance==1e-50,]$Target.Distance=NA
features[features$Pathway.Similarity==1e-50,]$Pathway.Similarity=NA
features[features$Integrated.Score==1e-50,]$Integrated.Score=NA
features[features$sim_ATC==1e-50,]$sim_ATC=NA
features[features$sim_structure==1e-50,]$sim_structure=NA
features[features$sim_sideeffects==1e-50,]$sim_sideeffects=NA
wilcox.test(features[features$Label == 'P',]$sim_ATC, features[features$Label == 'N',]$sim_ATC)
wilcox.test(features[features$Label == 'P',]$sim_structure, features[features$Label == 'N',]$sim_structure)
wilcox.test(features[features$Label == 'P',]$sim_sideeffects, features[features$Label == 'N',]$sim_sideeffects)
wilcox.test(features[features$Label == 'P',]$Integrated.Score, features[features$Label == 'N',]$Integrated.Score)
wilcox.test(features[features$Label == 'P',]$Target.Distance, features[features$Label == 'N',]$Target.Distance)
wilcox.test(features[features$Label == 'P',]$Pathway.Similarity, features[features$Label == 'N',]$Pathway.Similarity)
set.seed(1234)
res = data.frame(times = 1:1000, ATCsim_p = 0, STsim_p = 0, SEsim_p = 0, Integrated_p = 0, distance_p = 0, WWI_p = 0)
for(i in 1:1000){
  n_index = sample(which(features$Label == 'N'), size=nrow(features[features$Label == 'P',]))
  res[res$times == i,]$ATCsim_p = wilcox.test(features[features$Label == 'P',]$sim_ATC, features[n_index, ]$sim_ATC)$p.value
  res[res$times == i,]$STsim_p =wilcox.test(features[features$Label == 'P',]$sim_structure, features[n_index, ]$sim_structure)$p.value
  res[res$times == i,]$SEsim_p = wilcox.test(features[features$Label == 'P',]$sim_sideeffects, features[n_index, ]$sim_sideeffects)$p.value
  res[res$times == i,]$Integrated_p = wilcox.test(features[features$Label == 'P',]$Integrated.Score, features[n_index, ]$Integrated.Score)$p.value
  res[res$times == i,]$distance_p = wilcox.test(features[features$Label == 'P',]$Target.Distance, features[n_index, ]$Target.Distance)$p.value
  res[res$times == i,]$WWI_p = wilcox.test(features[features$Label == 'P',]$Pathway.Similarity, features[n_index, ]$Pathway.Similarity)$p.value
}

nrow(res[res$ATCsim_p < 0.05,])
nrow(res[res$STsim_p < 0.05,])
nrow(res[res$SEsim_p < 0.05,])
nrow(res[res$Integrated_p < 0.05,])
nrow(res[res$distance_p < 0.05,])
nrow(res[res$WWI_p < 0.05,])


#===============================================================================
#===============================================================================



drugnet_feature = read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/features.csv')

drugnet_feature$Label = factor(drugnet_feature$Label, levels = c('P','N'))
drugnet_feature[drugnet_feature$Target.Distance==1e-50,]$Target.Distance=NA
drugnet_feature[drugnet_feature$Pathway.Similarity==1e-50,]$Pathway.Similarity=NA
drugnet_feature[drugnet_feature$Integrated.Score==1e-50,]$Integrated.Score=NA
drugnet_feature[drugnet_feature$sim_ATC==1e-50,]$sim_ATC=NA
drugnet_feature[drugnet_feature$sim_structure==1e-50,]$sim_structure=NA
drugnet_feature[drugnet_feature$sim_sideeffects==1e-50,]$sim_sideeffects=NA

positive_atcsim <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/positive_drug_atc_sim2.txt',sep='\t', header = T, stringsAsFactors = F)
negative_atcsim <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/negative_drug_atc_sim2.txt',sep='\t', header = T, stringsAsFactors = F)

positive_TS <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/positive_TS.txt',sep='\t', header = TRUE, stringsAsFactors = F)
negative_TS <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/negative_TS.txt',sep='\t', header = TRUE, stringsAsFactors = F)

positive_SE <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/positive_SETS.txt',sep='\t', header = TRUE, stringsAsFactors = F)
negative_SE <- read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/negative_SETS.txt',sep='\t', header = TRUE, stringsAsFactors = F)

positive_atcsim$Label = 'P'
negative_atcsim$Label = 'N'
positive_TS$Label = 'P'
negative_TS$Label = 'N'
positive_SE$Label = 'P'
negative_SE$Label = 'N'
atc = rbind(positive_atcsim, negative_atcsim)
ts = rbind(positive_TS, negative_TS)
se = rbind(positive_SE, negative_SE)

names(atc) = c('A','B','value','Label')
names(ts) = c('A','B','value','Label')
names(se) = c('A','B','value','Label')

set.seed(1234)
res = data.frame(times = 1:1000, ATCsim_p = 0, STsim_p = 0, SEsim_p = 0, Integrated_p = 0, distance_p = 0, WWI_p = 0)
for(i in 1:1000){
  n_index = sample(which(drugnet_feature$Label == 'N'), size=nrow(drugnet_feature[drugnet_feature$Label == 'P',]))
  n_index1 = sample(which(atc$Label == 'N'), size=nrow(atc[atc$Label == 'P',]))
  n_index2 = sample(which(ts$Label == 'N'), size=nrow(ts[ts$Label == 'P',]))
  n_index3 = sample(which(se$Label == 'N'), size=nrow(se[se$Label == 'P',]))
  res[res$times == i,]$ATCsim_p = wilcox.test(atc[atc$Label == 'P',]$value, atc[n_index1, ]$value, alternative = 'greater')$p.value
  res[res$times == i,]$STsim_p = wilcox.test(ts[ts$Label == 'P',]$value, ts[n_index2, ]$value, alternative = 'greater')$p.value
  res[res$times == i,]$SEsim_p = wilcox.test(se[se$Label == 'P',]$value , se[n_index3, ]$value, alternative = 'greater')$p.value
  res[res$times == i,]$Integrated_p = wilcox.test(drugnet_feature[drugnet_feature$Label == 'P',]$Integrated.Score, drugnet_feature[n_index, ]$Integrated.Score, alternative = 'greater')$p.value
  res[res$times == i,]$distance_p = wilcox.test(drugnet_feature[drugnet_feature$Label == 'P',]$Target.Distance, drugnet_feature[n_index, ]$Target.Distance, alternative = 'less')$p.value
  res[res$times == i,]$WWI_p = wilcox.test(drugnet_feature[drugnet_feature$Label == 'P',]$Pathway.Similarity, drugnet_feature[n_index, ]$Pathway.Similarity, alternative = 'greater')$p.value
}
1-nrow(res[res$ATCsim_p < 0.05,])/1000
1-nrow(res[res$STsim_p < 0.05,])/1000
1-nrow(res[res$SEsim_p < 0.05,])/1000
1-nrow(res[res$Integrated_p < 0.05,])/1000
1-nrow(res[res$distance_p < 0.05,])/1000
1-nrow(res[res$WWI_p < 0.05,])/1000


atc_real_p = wilcox.test(atc[atc$Label=='P',]$value, atc[atc$Label=='N',]$value, alternative = 'greater')$p.value
ts_real_p = wilcox.test(ts[ts$Label=='P',]$value, ts[ts$Label=='N',]$value, alternative = 'greater')$p.value
se_real_p = wilcox.test(se[se$Label=='P',]$value, se[se$Label=='N',]$value, alternative = 'greater')$p.value
integrated_real_p = wilcox.test(drugnet_feature[drugnet_feature$Label=='P',]$Integrated.Score, drugnet_feature[drugnet_feature$Label=='N',]$Integrated.Score, alternative = 'greater')$p.value
dis_real_p = wilcox.test(drugnet_feature[drugnet_feature$Label=='P',]$Target.Distance, drugnet_feature[drugnet_feature$Label=='N',]$Target.Distance, alternative = 'less')$p.value
wwi_real_p = wilcox.test(drugnet_feature[drugnet_feature$Label=='P',]$Pathway.Similarity, drugnet_feature[drugnet_feature$Label=='N',]$Pathway.Similarity, alternative = 'greater')$p.value


vals <- boxplot(data = atc, main = "NX", value ~ Label, las = 2)
vals$out
library(car)
outliers<-Boxplot(value ~ Label, data=atc)
atc2 <- atc[-as.numeric(outliers),]
Boxplot(value ~ Label, data=atc2)
boxplot(atc[atc$Label == 'P',]$value)
out = boxplot.stats(atc[atc$Label == 'P',]$value)$out

atc_p1 = atc[atc$Label == 'P',]
out = boxplot.stats(atc_p1$value)$out
atc_p2 = atc_p1[-which(is.element(out, atc_p1$value)), ]
boxplot(atc_p2$value)

atc_n1 = atc[atc$Label == 'N',]
out = boxplot.stats(atc_n1$value)$out
atc_n2 = atc_n1[-which(is.element(out, atc_n1$value)), ]
boxplot(atc_n2$value)

wilcox.test(atc_p2$value, atc_n2$value)
wilcox.test(atc_p2$value, atc_n2$value)

# pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/ATCsim_permP.pdf')
# plot(density(res$ATCsim_p))
# abline(v = atc_real_p,col='red')
# abline(v= 0.05 ,col='green')
# dev.off()
# pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/STsim_permP.pdf')
# plot(density(res$STsim_p))
# abline(v = ts_real_p,col='red')
# dev.off()
# pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/SEsim_permP.pdf')
# plot(density(res$SEsim_p))
# abline(v = se_real_p,col='red')
# dev.off()
# pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/Integrated_permP.pdf')
# plot(density(res$Integrated_p))
# abline(v = integrated_real_p,col='red')
# dev.off()
# pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/distance_permP.pdf')
# plot(density(res$distance_p))
# abline(v = dis_real_p,col='red')
# dev.off()
# pdf('G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/WWI_permP.pdf')
# plot(density(res$WWI_p))
# abline(v = wwi_real_p,col='red')
# dev.off()

library(patchwork)


p1 = ggplot(res, aes(x=ATCsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = atc_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5,show_guide=TRUE)+
  ggtitle(paste0('sim_ATC','\n Permutation_Pvalue = ',1-nrow(res[res$ATCsim_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black')) +
  xlab('P value')

p2 = ggplot(res, aes(x=STsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = ts_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_structure','\n Permutation_Pvalue = ', 1-nrow(res[res$STsim_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p3 = ggplot(res, aes(x=SEsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = se_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_sideeffect','\n Permutation_Pvalue = ', 1-nrow(res[res$SEsim_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p4 = ggplot(res, aes(x=Integrated_p)) +
  geom_density()+
  geom_vline(aes(xintercept = integrated_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_pharmacology','\n Permutation_Pvalue = ', 1-nrow(res[res$Integrated_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')


p5 = ggplot(res, aes(x=distance_p)) +
  geom_density()+
  geom_vline(aes(xintercept = dis_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Target_Distance','\n Permutation_Pvalue = ', 1-nrow(res[res$distance_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p6 = ggplot(res, aes(x=WWI_p)) +
  geom_density()+
  geom_vline(aes(xintercept = wwi_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Pathway Similarity','\n Permutation_Pvalue = ', 1-nrow(res[res$WWI_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p = (p1|p2|p3)/(p4+p5+p6)
p
ggsave( 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/features_permP.pdf',p, width = 10, height = 5)

p1 = ggplot(res, aes(x=ATCsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = atc_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5,show_guide=TRUE)+
  ggtitle(paste0('sim_ATC','\n Permutation_Pvalue = ', (nrow(res[res$ATCsim_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black')) +
  xlab('P value')

p2 = ggplot(res, aes(x=STsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = ts_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_structure','\n Permutation_Pvalue = ', (nrow(res[res$STsim_p < 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p3 = ggplot(res, aes(x=SEsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = se_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_sideeffect','\n Permutation_Pvalue = ', (nrow(res[res$SEsim_p < 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p4 = ggplot(res, aes(x=Integrated_p)) +
  geom_density()+
  geom_vline(aes(xintercept = integrated_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_pharmacology','\n Permutation_Pvalue = ', (nrow(res[res$Integrated_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')


p5 = ggplot(res, aes(x=distance_p)) +
  geom_density()+
  geom_vline(aes(xintercept = dis_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Target_Distance','\n Permutation_Pvalue = ', (nrow(res[res$distance_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p6 = ggplot(res, aes(x=WWI_p)) +
  geom_density()+
  geom_vline(aes(xintercept = wwi_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Pathway Similarity','\n Permutation_Pvalue = ', (nrow(res[res$WWI_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p = (p1|p2|p3)/(p4+p5+p6)
p
ggsave( 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/features_permP_2.pdf',p, width = 10, height = 5)


set.seed(1234)
res = data.frame(times = 1:1000, ATCsim_p = 0, STsim_p = 0, SEsim_p = 0, Integrated_p = 0, distance_p = 0, WWI_p = 0)
for(i in 1:1000){
  n_index = sample(which(drugnet_feature$Label == 'N'), size=nrow(drugnet_feature[drugnet_feature$Label == 'P',]))
  n_index1 = sample(which(atc$Label == 'N'), size=nrow(atc[atc$Label == 'P',]))
  n_index2 = sample(which(ts$Label == 'N'), size=nrow(ts[ts$Label == 'P',]))
  n_index3 = sample(which(se$Label == 'N'), size=nrow(se[se$Label == 'P',]))
  res[res$times == i,]$ATCsim_p = wilcox.test(atc[atc$Label == 'P',]$value, atc[n_index1, ]$value)$p.value
  res[res$times == i,]$STsim_p = wilcox.test(ts[ts$Label == 'P',]$value, ts[n_index2, ]$value)$p.value
  res[res$times == i,]$SEsim_p = wilcox.test(se[se$Label == 'P',]$value , se[n_index3, ]$value)$p.value
  res[res$times == i,]$Integrated_p = wilcox.test(drugnet_feature[drugnet_feature$Label == 'P',]$Integrated.Score, drugnet_feature[n_index, ]$Integrated.Score)$p.value
  res[res$times == i,]$distance_p = wilcox.test(drugnet_feature[drugnet_feature$Label == 'P',]$Target.Distance, drugnet_feature[n_index, ]$Target.Distance)$p.value
  res[res$times == i,]$WWI_p = wilcox.test(drugnet_feature[drugnet_feature$Label == 'P',]$Pathway.Similarity, drugnet_feature[n_index, ]$Pathway.Similarity)$p.value
}

(nrow(res[res$ATCsim_p >= 0.05,])+1)/1000
(nrow(res[res$STsim_p >= 0.05,])+1)/1000
(nrow(res[res$SEsim_p >= 0.05,])+1)/1000
(nrow(res[res$Integrated_p >= 0.05,])+1)/1000
(nrow(res[res$distance_p >= 0.05,])+1)/1000
(nrow(res[res$WWI_p >= 0.05,])+1)/1000


atc_real_p = wilcox.test(atc[atc$Label=='P',]$value, atc[atc$Label=='N',]$value)$p.value
ts_real_p = wilcox.test(ts[ts$Label=='P',]$value, ts[ts$Label=='N',]$value)$p.value
se_real_p = wilcox.test(se[se$Label=='P',]$value, se[se$Label=='N',]$value)$p.value
integrated_real_p = wilcox.test(drugnet_feature[drugnet_feature$Label=='P',]$Integrated.Score, drugnet_feature[drugnet_feature$Label=='N',]$Integrated.Score)$p.value
dis_real_p = wilcox.test(drugnet_feature[drugnet_feature$Label=='P',]$Target.Distance, drugnet_feature[drugnet_feature$Label=='N',]$Target.Distance)$p.value
wwi_real_p = wilcox.test(drugnet_feature[drugnet_feature$Label=='P',]$Pathway.Similarity, drugnet_feature[drugnet_feature$Label=='N',]$Pathway.Similarity)$p.value


library(patchwork)


p1 = ggplot(res, aes(x=ATCsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = atc_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5,show_guide=TRUE)+
  ggtitle(paste0('sim_ATC','\n Permutation_Pvalue = ',1-nrow(res[res$ATCsim_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black')) +
  xlab('P value')

p2 = ggplot(res, aes(x=STsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = ts_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_structure','\n Permutation_Pvalue = ', 1-nrow(res[res$STsim_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p3 = ggplot(res, aes(x=SEsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = se_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_sideeffect','\n Permutation_Pvalue = ', 1-nrow(res[res$SEsim_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p4 = ggplot(res, aes(x=Integrated_p)) +
  geom_density()+
  geom_vline(aes(xintercept = integrated_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_pharmacology','\n Permutation_Pvalue = ', 1-nrow(res[res$Integrated_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')


p5 = ggplot(res, aes(x=distance_p)) +
  geom_density()+
  geom_vline(aes(xintercept = dis_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Target_Distance','\n Permutation_Pvalue = ', 1-nrow(res[res$distance_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p6 = ggplot(res, aes(x=WWI_p)) +
  geom_density()+
  geom_vline(aes(xintercept = wwi_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Pathway Similarity','\n Permutation_Pvalue = ', 1-nrow(res[res$WWI_p < 0.05,])/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p = (p1|p2|p3)/(p4+p5+p6)
p
ggsave( 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/features_permP2.pdf',p, width = 10, height = 5)


p1 = ggplot(res, aes(x=ATCsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = atc_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5,show_guide=TRUE)+
  ggtitle(paste0('sim_ATC','\n Permutation_Pvalue = ', (nrow(res[res$ATCsim_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black')) +
  xlab('P value')

p2 = ggplot(res, aes(x=STsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = ts_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_structure','\n Permutation_Pvalue = ', (nrow(res[res$STsim_p < 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p3 = ggplot(res, aes(x=SEsim_p)) +
  geom_density()+
  geom_vline(aes(xintercept = se_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_sideeffect','\n Permutation_Pvalue = ', (nrow(res[res$SEsim_p < 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p4 = ggplot(res, aes(x=Integrated_p)) +
  geom_density()+
  geom_vline(aes(xintercept = integrated_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('sim_pharmacology','\n Permutation_Pvalue = ', (nrow(res[res$Integrated_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')


p5 = ggplot(res, aes(x=distance_p)) +
  geom_density()+
  geom_vline(aes(xintercept = dis_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Target_Distance','\n Permutation_Pvalue = ', (nrow(res[res$distance_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p6 = ggplot(res, aes(x=WWI_p)) +
  geom_density()+
  geom_vline(aes(xintercept = wwi_real_p),
             color="red", size=0.3)+
  geom_vline(aes(xintercept = 0.05),
             color="green", linetype="dashed", size=0.5)+
  ggtitle(paste0('Pathway Similarity','\n Permutation_Pvalue = ', (nrow(res[res$WWI_p >= 0.05,])+1)/1000)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "transparent",colour = 'black'))+
  xlab('P value')

p = (p1|p2|p3)/(p4+p5+p6)
p
ggsave( 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/features/plot/features_permP2_2.pdf',p, width = 10, height = 5)
