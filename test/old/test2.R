
library(DComboNet)
library(Hmisc)
library(reshape2)
options(stringsAsFactors = F)
outputdir = 'G:/lab/Projects/p1_DComboNet/药物组合/DrugCombinations_Old_files/DComboNet/OCL_LY3/'
# setwd(outputdir)
#
# fun_file=list.files(paste0(outputdir,'/DComboNet/R/'))
# fun_dir = paste0(outputdir,'/DComboNet/R/',fun_file)
# for(file in 1:length(fun_file)){
#   source(fun_dir[file])
# }

drugcandidate <-  read.csv(paste0(outputdir,"druglist_oci_ly3.csv"), sep = ',', header = T, stringsAsFactors = F)
drugpair = read.csv(paste0(outputdir,'drugpair_oci_ly3.csv'),stringsAsFactors = F)
candidateFP1 <- read.csv(paste0(outputdir,'drug_net/structure_sim/fingerprint.csv'), sep = ',',header = T, stringsAsFactors = F)
candidateFP2 <- read.csv(paste0(outputdir,'drug_net/structure_sim/PubChem.csv'), sep = ',',header = T, stringsAsFactors = F)
candidateFP3 <- read.csv(paste0(outputdir,'drug_net/structure_sim/MACCF.csv'), sep = ',',header = T, stringsAsFactors = F)
drugnetWeight = TRUE
featuretype = 'integrated_score'
drugtarget = read.csv(paste0(outputdir,'drug_gene/drug_target_oci_ly3.csv'), sep = ',',header = T, stringsAsFactors = F)
druggene = NULL
essentialgene = NULL
cellline = 'OCILY3'

# drugDEG = read.table(paste0(outputdir,'drug_gene/DEG/drug_DEG_12hrs.txt'),sep='\t',header=T,stringsAsFactors = F)
# drugDEG = drugDEG[c(1,2,4)]
cancergene = read.table(paste0(outputdir,'gene_network/OCILY3_genelist.txt'),sep='\t',header=T,stringsAsFactors = F)

dtweight = 2
dgweight = 1
dDEGweight = 1
dgpWeight = 1
dDEGpWeight = 1

model = 'L2'
x = 0.5
y = 0.5
z = 0
A = 0.5
B = 0.5
C = 0

drugDEG_overlap = read.table(paste0('G:/lab/Projects/p1_DComboNet/药物组合//DrugCombinations_Old_files/DComboNet/OCL_LY3/drug_gene/DEG/DEG_12h_0_0.05_doseoverlap.txt'),sep='\t',header=T,stringsAsFactors = F)
resultdir_overlap = paste0('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/result_overlap7/')

dir.create(resultdir_overlap)
DComboNet(load_dir = outputdir,
          resultdir = resultdir_overlap,
          model = 'L2',
          drugcandidate=drugcandidate,
          CDK_FP = NULL,
          pubchemFP = NULL,
          MACCS_FP = NULL,
          drugnetWeight = TRUE,
          featuretype = featuretype,
          drugtarget = drugtarget,
          druggene = NULL,
          dataset = NULL,
          cellline = 'OCILY3',
          treatment_time = NULL,
          drugDEG = drugDEG_overlap,
          drugDEP = NULL,
          cancergene = cancergene,
          dtweight = 2,
          dgweight = 1,
          dDEGweight = 1,
          dgpweight = 1,
          dDEpweight = 1,
          eta=1,
          r=0.7)


resultdir = resultdir_overlap

# Level2
library(Hmisc)
require(MESS)
require(pROC)
require(data.table)
require(ggplot2)
require(ROCR)
require(RCurl)


model = 'L2'



goldstandard = read.table('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_results/goldstandard.csv',header = T, sep = ',', stringsAsFactors = F)
names(goldstandard) = c('A', 'B', 'EOB', 'EOBerror')
goldstandard$snr = abs(goldstandard$EOB / goldstandard$EOBerror)
goldstandard$Label = 'additive'
goldstandard[goldstandard$EOB > 0 & goldstandard$snr >2,]$Label = 'synergistic'
goldstandard[goldstandard$EOB < 0 & goldstandard$snr >2,]$Label = 'antagonistic'
#write.csv(goldstandard,'goldstandard.csv',row.names = F, quote = F)
goldstandard$A = toupper(goldstandard$A)
goldstandard$B = toupper(goldstandard$B)
#DComboNet_results
#setwd('C:/Users/fengf/Documents/Veronica_files/DCcomboNet/Rpackage/OCL_LY3_3/')
goldstandard[goldstandard$A == "MITOMYCIN C",]$A = 'MITOMYCIN'
goldstandard[goldstandard$B == "MITOMYCIN C",]$B = 'MITOMYCIN'
goldstandard[goldstandard$A == 'RAPAMYCIN',]$A = 'SIROLIMUS'
goldstandard[goldstandard$B == 'RAPAMYCIN',]$B = 'SIROLIMUS'
model='L2'

a = list.files(paste0(resultdir,model,'_result/drugrank/'))
dir = paste0(resultdir,model,'_result/drugrank/',a)

n = length(a) ;n
b = gsub('_rank.csv','',a)
drugrank = read.csv(dir[1],stringsAsFactors = F)
drugrank$drugseed = b[1]
# drugrank$Score = (drugrank$Score - mean(drugrank$Score))/sd(drugrank$Score)
for(i in 2:n){
  drugrank1 = read.csv(dir[i],stringsAsFactors = F)
  drugrank1$drugseed = b[i]
  # drugrank1$Score = (drugrank1$Score - mean(drugrank1$Score))/sd(drugrank1$Score)

  drugrank = rbind(drugrank,drugrank1)
}
drugrank = unique(drugrank[c(4,1,2,3)])
names(drugrank) = c('drugseed','DrugID','Score','Rank')
drugrank$drugseed = toupper(drugrank$drugseed)
drugrank$DrugID = toupper(drugrank$DrugID)


#drugrank[drugrank$drugseed == 'MITOMYCIN C',]$drugseed = 'MITOMYCIN'
#drugrank[drugrank$DrugID == 'MITOMYCIN C',]$DrugID = 'MITOMYCIN'
#drugrank[drugrank$drugseed == 'SIROLIMUS',]$drugseed = 'RAPAMYCIN'
#drugrank[drugrank$DrugID == 'SIROLIMUS',]$DrugID = 'RAPAMYCIN'
#drugrank[drugrank$drugseed == 'ACLARUBICIN',]$drugseed = "ACLACINOMYCIN A"
#drugrank[drugrank$DrugID == 'ACLARUBICIN',]$DrugID = "ACLACINOMYCIN A"
names(drugrank) = c('A','B','Score','Rank')

drugrank1 = unique(goldstandard[c('A','B','Label')])
drugrank1[which(drugrank1$A == 'MITOMYCIN C'),]$A = 'MITOMYCIN'
drugrank1[which(drugrank1$B == 'MITOMYCIN C'),]$B = 'MITOMYCIN'
drugrank1 = merge(drugrank[c(1,2,3)], drugrank1, by =c('A','B'))
drugrank1$A_B = as.numeric(0)
drugrank1$B_A = as.numeric(0)
drugrank1$A_B_score = as.numeric(0)
drugrank1$B_A_score = as.numeric(0)
drugrank1$pred_tag_L2plus = 0
drugrank1$Score = 0
drugrank1$rank = 0
num_drug = max(drugrank$Rank)
pres_syn=0.1
pres_ant = 0.5


for(i in 1:nrow(drugrank1)){

  seedA = drugrank[drugrank$A == drugrank1[i,1] & drugrank$B == drugrank1[i,2],]
  seedB = drugrank[drugrank$B == drugrank1[i,1] & drugrank$A == drugrank1[i,2],]

  if(nrow(seedA)>0 & nrow(seedB)>0){

    drugrank1[i,]$A_B = round(mean(seedA$Rank),0)
    drugrank1[i,]$B_A = round(mean(seedB$Rank),0)

    drugrank1[i,]$A_B_score = mean(seedA$Score)
    drugrank1[i,]$B_A_score = mean(seedB$Score)

    if(drugrank1[i,]$A_B <= num_drug*pres_syn & drugrank1[i,]$B_A <= num_drug*pres_syn){

      drugrank1[i,]$pred_tag_L2plus = 'Syn'
      drugrank1[i,]$Score = max(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = min(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B <= num_drug*pres_syn & (drugrank1[i,]$B_A >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant))|
             ((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$A_B <= num_drug*pres_ant) & drugrank1[i,]$B_A  <= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Syn'
      drugrank1[i,]$Score = max(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = min(seedA$Rank, seedB$Rank)
    }else if(drugrank1[i,]$A_B >= num_drug*pres_ant & drugrank1[i,]$B_A >= num_drug*pres_ant){

      drugrank1[i,]$pred_tag_L2plus = 'Ant'
      drugrank1[i,]$Score = min(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = max(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B >= num_drug*pres_ant & (drugrank1[i,]$B_A >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant))|
             ((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$A_B <= num_drug*pres_ant) & drugrank1[i,]$B_A  >= num_drug*pres_ant)){

      drugrank1[i,]$pred_tag_L2plus = 'Ant'
      drugrank1[i,]$Score = min(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = max(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$B_A <= num_drug*pres_ant)|
             (drugrank1[i,]$A_B <= num_drug*pres_ant & drugrank1[i,]$B_A >= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Add'
      drugrank1[i,]$Score = mean(mean(seedA$Score), mean(seedB$Score))
      drugrank1[i,]$rank = mean(round(mean(seedA$Rank),0), round(mean(seedB$Rank),0))
    }else if((drugrank1[i,]$A_B <= num_drug*pres_syn & drugrank1[i,]$B_A >= num_drug*pres_ant)|
             (drugrank1[i,]$A_B >= num_drug*pres_ant & drugrank1[i,]$B_A <= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Add'
      drugrank1[i,]$Score = mean(mean(seedA$Score), mean(seedB$Score))
      drugrank1[i,]$rank = mean(round(mean(seedA$Rank),0), round(mean(seedB$Rank),0))

    }
  }else if(nrow(seedA)>0 & nrow(seedB)==0){

    drugrank1[i,]$A_B = round(mean(seedA$Rank),0)
    drugrank1[i,]$B_A = round(mean(seedA$Rank),0)
    drugrank1[i,]$A_B_score = mean(seedA$Score)
    drugrank1[i,]$B_A_score = mean(seedB$Score)

    if(drugrank1[i,]$A_B <= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_syn){

      drugrank1[i,]$pred_tag_L2plus = 'Syn'
      drugrank1[i,]$Score = max(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = min(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B <= num_drug*pres_syn & (drugrank1[i,]$B_A  >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant))|
             ((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$A_B <= num_drug*pres_ant) & drugrank1[i,]$B_A  <= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Syn'
      drugrank1[i,]$Score = max(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = min(seedA$Rank, seedB$Rank)
    }else if(drugrank1[i,]$A_B >= num_drug*pres_ant & drugrank1[i,]$B_A >= num_drug*pres_ant){

      drugrank1[i,]$pred_tag_L2plus = 'Ant'
      drugrank1[i,]$Score = min(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = max(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B >= num_drug*pres_ant & (drugrank1[i,]$B_A  >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant))|
             ((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$A_B <= num_drug*pres_ant) & drugrank1[i,]$B_A  >= num_drug*pres_ant)){

      drugrank1[i,]$pred_tag_L2plus = 'Ant'
      drugrank1[i,]$Score = min(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = max(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant)|
             (drugrank1[i,]$A_B <= num_drug*pres_ant & drugrank1[i,]$B_A  >= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Add'
      drugrank1[i,]$Score = mean(seedA$Score)
      drugrank1[i,]$rank = round(mean(seedA$Rank),0)
    }else if((drugrank1[i,]$A_B <= num_drug*pres_syn & drugrank1[i,]$B_A  >= num_drug*pres_ant)|
             (drugrank1[i,]$A_B >= num_drug*pres_ant & drugrank1[i,]$B_A  <= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Add'
      drugrank1[i,]$Score = mean(seedA$Score)
      drugrank1[i,]$rank = round(mean(seedA$Rank),0)
    }
  }else if(nrow(seedA)==0 & nrow(seedB)>0){

    drugrank1[i,]$A_B = round(mean(seedB$Rank),0)
    drugrank1[i,]$B_A = round(mean(seedB$Rank),0)
    drugrank1[i,]$A_B_score = mean(seedA$Score)
    drugrank1[i,]$B_A_score = mean(seedB$Score)

    if(drugrank1[i,]$A_B <= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_syn){

      drugrank1[i,]$pred_tag_L2plus = 'Syn'
      drugrank1[i,]$Score = max(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = min(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B <= num_drug*pres_syn & (drugrank1[i,]$B_A  >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant))|
             ((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$A_B <= num_drug*pres_ant) & drugrank1[i,]$B_A  <= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Syn'
      drugrank1[i,]$Score = max(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = min(seedA$Rank, seedB$Rank)
    }else if(drugrank1[i,]$A_B >= num_drug*pres_ant & drugrank1[i,]$B_A >= num_drug*pres_ant){

      drugrank1[i,]$pred_tag_L2plus = 'Ant'
      drugrank1[i,]$Score = min(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = max(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B >= num_drug*pres_ant & (drugrank1[i,]$B_A  >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant))|
             ((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$A_B <= num_drug*pres_ant) & drugrank1[i,]$B_A  >= num_drug*pres_ant)){

      drugrank1[i,]$pred_tag_L2plus = 'Ant'
      drugrank1[i,]$Score = min(seedA$Score, seedB$Score)
      drugrank1[i,]$rank = max(seedA$Rank, seedB$Rank)
    }else if((drugrank1[i,]$A_B >= num_drug*pres_syn & drugrank1[i,]$B_A  <= num_drug*pres_ant)|
             (drugrank1[i,]$A_B <= num_drug*pres_ant & drugrank1[i,]$B_A  >= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Add'
      drugrank1[i,]$Score = mean(seedB$Score)
      drugrank1[i,]$rank = round(mean(seedB$Rank),0)
    }else if((drugrank1[i,]$A_B <= num_drug*pres_syn & drugrank1[i,]$B_A  >= num_drug*pres_ant)|
             (drugrank1[i,]$A_B >= num_drug*pres_ant & drugrank1[i,]$B_A  <= num_drug*pres_syn)){

      drugrank1[i,]$pred_tag_L2plus = 'Add'
      drugrank1[i,]$Score = mean(seedB$Score)
      drugrank1[i,]$rank = round(mean(seedB$Rank),0)
    }
  }else{
    print('not predicted')
  }

}

library(ROCR)
drugrank1$label = 0
drugrank1[drugrank1$Label == 'synergistic',]$label = 1
drugrank1$pred.label = 0
drugrank1[drugrank1$pred_tag_L2plus == 'Syn',]$pred.label = 1

drugrank1$AplusB_Score = drugrank1$A_B_score + drugrank1$B_A_score
nrow(drugrank1[drugrank1$pred_tag_L2plus != 'Syn' & drugrank1$Label != 'synergistic',])
nrow(drugrank1[drugrank1$pred_tag_L2plus != 'Syn' & drugrank1$Label != 'synergistic',])/nrow(drugrank1[drugrank1$Label != 'synergistic',])
nrow(drugrank1[drugrank1$pred_tag_L2plus == 'Syn' & drugrank1$Label == 'synergistic',])
nrow(drugrank1[drugrank1$pred_tag_L2plus == 'Syn' & drugrank1$Label == 'synergistic',])/nrow(drugrank1[drugrank1$Label == 'synergistic',])

pred <- prediction( drugrank1$AplusB_Score, drugrank1$label)
perf <- performance( pred, 'tpr', 'fpr' )
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AucValue<-round(AUCValue,3)
AucValue


drugrank1 = drugrank1[order(drugrank1$AplusB_Score,decreasing = T),]
drugrank1$AplusB_Rank = 1: nrow(drugrank1)

nrow(drugrank1[drugrank1$AplusB_Rank > 20 & drugrank1$Label != 'synergistic',])
nrow(drugrank1[drugrank1$Label != 'synergistic',])
nrow(drugrank1[drugrank1$AplusB_Rank > 20 & drugrank1$Label != 'synergistic',])/nrow(drugrank1[drugrank1$Label != 'synergistic',])
nrow(drugrank1[drugrank1$AplusB_Rank <= 20 & drugrank1$Label == 'synergistic',])
nrow(drugrank1[drugrank1$Label == 'synergistic',])
nrow(drugrank1[drugrank1$AplusB_Rank <= 20 & drugrank1$Label == 'synergistic',])/nrow(drugrank1[drugrank1$Label == 'synergistic',])

tmp_s = drugrank1[drugrank1$A != 'SIROLIMUS' & drugrank1$B != 'SIROLIMUS' ,]
tmp_s = tmp_s[order(tmp_s$AplusB_Score,decreasing = T),]
tmp_s$AplusB_Rank = 1:nrow(tmp_s)
nrow(tmp_s[tmp_s$pred_tag_L2plus != 'Syn' & tmp_s$Label != 'synergistic',])
nrow(tmp_s[tmp_s$Label != 'synergistic',])
nrow(tmp_s[tmp_s$pred_tag_L2plus != 'Syn' & tmp_s$Label != 'synergistic',])/nrow(tmp_s[tmp_s$Label != 'synergistic',])
nrow(tmp_s[tmp_s$pred_tag_L2plus == 'Syn' & tmp_s$Label == 'synergistic',])
nrow(tmp_s[tmp_s$Label == 'synergistic',])
nrow(tmp_s[tmp_s$pred_tag_L2plus == 'Syn' & tmp_s$Label == 'synergistic',])/nrow(tmp_s[tmp_s$Label == 'synergistic',])

nrow(tmp_s[tmp_s$AplusB_Rank > 20 & tmp_s$Label != 'synergistic',])
nrow(tmp_s[drugrank1$Label != 'synergistic',])
nrow(tmp_s[tmp_s$AplusB_Rank > 20 & tmp_s$Label != 'synergistic',])/nrow(tmp_s[drugrank1$Label != 'synergistic',])
nrow(tmp_s[tmp_s$AplusB_Rank <= 20 & tmp_s$Label == 'synergistic',])
nrow(tmp_s[drugrank1$Label == 'synergistic',])
nrow(tmp_s[tmp_s$AplusB_Rank <= 20 & tmp_s$Label == 'synergistic',])/nrow(tmp_s[drugrank1$Label == 'synergistic',])

pred <- prediction( tmp_s$AplusB_Score, tmp_s$label)
perf <- performance( pred, 'tpr', 'fpr' )
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AucValue<-round(AUCValue,3)
AucValue


tmp_m = drugrank1[drugrank1$A != 'MITOMYCIN' & drugrank1$B != 'MITOMYCIN' ,]
tmp_m = tmp_m[order(tmp_m$AplusB_Score,decreasing = T),]
tmp_m$AplusB_Rank = 1:nrow(tmp_m)
nrow(tmp_m[tmp_m$pred_tag_L2plus != 'Syn' & tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$pred_tag_L2plus != 'Syn' & tmp_m$Label != 'synergistic',])/nrow(tmp_m[tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$pred_tag_L2plus == 'Syn' & tmp_m$Label == 'synergistic',])/nrow(tmp_m[tmp_m$Label == 'synergistic',])

nrow(tmp_m[tmp_m$AplusB_Rank > 20 & tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$AplusB_Rank > 20 & tmp_m$Label != 'synergistic',])/nrow(tmp_m[tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$AplusB_Rank <= 20 & tmp_m$Label == 'synergistic',])
nrow(tmp_m[tmp_m$Label == 'synergistic',])
nrow(tmp_m[tmp_m$AplusB_Rank <= 20 & tmp_m$Label == 'synergistic',])/nrow(tmp_m[tmp_m$Label == 'synergistic',])


pred <- prediction( tmp_m$AplusB_Score, tmp_m$label)
perf <- performance( pred, 'tpr', 'fpr' )
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AucValue<-round(AUCValue,3)
AucValue


ranking_result <- Ranking(
  drugcandidate= drugcandidate$Drug,
  model='L2',
  resultdir=resultdir_overlap,
  filename = 'ranking_test.csv',
  outputdir = resultdir_overlap,
  methods = "sum"
)
names(ranking_result) = c('B','A','Score_sum','Rank_sum')
tmp = merge(goldstandard, ranking_result, by=c('A','B'))
tmp$label=0
tmp[tmp$Label == 'synergistic',]$label = 1
pred <- prediction( tmp$Score_sum, tmp$label)
perf <- performance( pred, 'tpr', 'fpr' )
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AucValue<-round(AUCValue,3)
AucValue


tmp_m = tmp[tmp$A != 'MITOMYCIN' & tmp$B != 'MITOMYCIN' ,]
tmp_m = tmp_m[order(tmp_m$AplusB_Score, decreasing = T), ]
tmp_m$rank = 1:nrow(tmp_m)
pred <- prediction( tmp_m$Score_sum, tmp_m$label)
perf <- performance( pred, 'tpr', 'fpr' )
aucPerf <- performance( pred, 'auc' )
AUCValue<-aucPerf@y.values[[1]]
AucValue<-round(AUCValue,3)
AucValue

nrow(tmp_m[tmp_m$rank > 20 & tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$rank > 20 & tmp_m$Label != 'synergistic',])/nrow(tmp_m[tmp_m$Label != 'synergistic',])
nrow(tmp_m[tmp_m$rank <= 20 & tmp_m$Label == 'synergistic',])
nrow(tmp_m[tmp_m$Label == 'synergistic',])
nrow(tmp_m[tmp_m$rank <= 20 & tmp_m$Label == 'synergistic',])/nrow(tmp_m[tmp_m$Label == 'synergistic',])
