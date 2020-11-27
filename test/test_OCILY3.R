# library(DComboNet)
library(Hmisc)
library(reshape2)
options(stringsAsFactors = F)
outputdir = 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/Rpackage/DComboNet-master/OCI_LY3/'

drugcandidate <-  read.csv(paste0(outputdir,"druglist_oci_ly3.csv"), sep = ',', header = T, stringsAsFactors = F)
drugpair = read.csv(paste0(outputdir,'drugpair_oci_ly3.csv'),stringsAsFactors = F)
candidateFP1 <- read.csv(paste0(outputdir,'drug_net/structure_sim/fingerprint.csv'), sep = ',',header = T, stringsAsFactors = F)
candidateFP2 <- read.csv(paste0(outputdir,'drug_net/structure_sim/PubChem.csv'), sep = ',',header = T, stringsAsFactors = F)
candidateFP3 <- read.csv(paste0(outputdir,'drug_net/structure_sim/MACCF.csv'), sep = ',',header = T, stringsAsFactors = F)
drugnetWeight = TRUE
featuretype = 'integrated_score'
drugtarget = read.csv(paste0(outputdir,'drug_target_oci_ly3.csv'), sep = ',',header = T, stringsAsFactors = F)
druggene = NULL
essentialgene = NULL
cellline = 'OCILY3'

drugDEG = read.table(paste0(outputdir,'drug_DEG_12hrs.txt'),sep='\t',header=T,stringsAsFactors = F)
cancergene = read.table(paste0(outputdir,'OCILY3_genelist.txt'),sep='\t',header=T,stringsAsFactors = F)

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

resultdir_overlap = paste0('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/result_overlap8/')
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
          drugDEG = drugDEG,
          drugDEP = NULL,
          cancergene = cancergene,
          dtweight = 2,
          dgweight = 1,
          dDEGweight = 1,
          dgpweight = 1,
          dDEpweight = 1,
          eta=1,
          r=0.7)

ranking_result <- Ranking(
  drugcandidate= drugcandidate$Drug,
  model='L2',
  resultdir=resultdir_overlap,
  filename = 'ranking_test',
  outputdir = resultdir_overlap,
  methods = "sum"
)


resultdir = resultdir_overlap

# Level2

require(ROCR)
model = 'L2'


goldstandard = read.table('G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_results/goldstandard.csv',header = T, sep = ',', stringsAsFactors = F)
names(goldstandard) = c('A', 'B', 'EOB', 'EOBerror')
goldstandard$snr = abs(goldstandard$EOB / goldstandard$EOBerror)
goldstandard$Label = 'additive'
goldstandard[goldstandard$EOB > 0 & goldstandard$snr > 2,]$Label = 'synergistic'
goldstandard[goldstandard$EOB < 0 & goldstandard$snr > 2,]$Label = 'antagonistic'
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
tmp_m = tmp_m[order(tmp_m$Score_sum, decreasing = T), ]
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
