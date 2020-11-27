library(igraph)
library(data.table)
library(compiler)
library(Matrix)
library(reshape2)
library(Hmisc)
library(DComboNet)

outputdir = "/picb/bigdata/project/FengFYM/DComboNetV2/OCL_LY3/"
setwd(outputdir)
options(stringsAsFactors = F)


#outputdir = "C:/Users/fengf/Documents/Veronica_files/DCcomboNet/Rpackage/OCL_LY3/"
suppressPackageStartupMessages(library(data.table))

drugpair = read.csv(paste0(outputdir,"data/drugpair.csv"),stringsAsFactors = F)[-3]
drugcandidate <-  NULL #read.csv(paste0(outputdir,"data/druglist.csv"), sep = ",", header = T, stringsAsFactors = F)
candidateFP1 <- read.csv(paste0(outputdir,"data/fingerprints/fingerprints.csv"), sep = ",",header = T, stringsAsFactors = F)
candidateFP2 <- read.csv(paste0(outputdir,"data/fingerprints/pubchem_fingerprints.csv"), sep = ",",header = T, stringsAsFactors = F)
candidateFP3 <- read.csv(paste0(outputdir,"data/fingerprints/MACCS_fingerprints.csv"), sep = ",",header = T, stringsAsFactors = F)
drugtarget = read.csv(paste0(outputdir,"data/drugtarget.csv"), sep = ",",header = T, stringsAsFactors = F)
druggene = NULL
dt_lib = read.csv(paste0(outputdir,"data/drugtarget.csv"))
druggene = read.csv(paste0(outputdir,"data/drug_gene/drug_ipa.inPPI.csv"), sep = ",",header = T, stringsAsFactors = F)
# druggene = druggene[druggene$Relationship.Type != 'expression',]
druggene = druggene[-2]


model = "L1"
dataset = NULL
treatment_time = NULL
cellline = NULL

drugnetWeight = T
featuretype = "integrated_score"

essentialgene = NULL
drugDEG = NULL
cancergene = NULL
dtweight = 2
dgweight = 1
dDEGweight = 1
dgpWeight = 1
dDEGpWeight = 1
gpWeight = 1


x = 0.5
y = 0.5
z = 0
A = 0.5
B = 0.5
C = 0
eta=1


# resultdir = "/picb/bigdata/project/FengFYM/DComboNetV2/wwin_LOOCV12/"
resultdir = "G:/lab/Projects/p1_DComboNet/DComboNet_improve/Rpackage/test_LOOCV/"

dir.create(resultdir)

L1.LOCOCV(load_dir,
          resultdir = resultdir,
          model = "L1",
          drugcandidate = NULL,
          drugnetWeight = drugnetWeight,
          featuretype = featuretype,
          drugtarget = drugtarget,
          druggene = NULL,
          cancergene = NULL,
          dtweight = 2,
          dgweight = 1,
          dgpweight = 1,
          gpweight = 1,
          x = 0.5,
          y = 0.5,
          z = 0,
          A = 0.5,
          B = 0.5,
          eta = 1,
          r = 0.7)



outputdir = "/picb/bigdata/project/FengFYM/DComboNetV2/OCL_LY3/"
# setwd(outputdir)
options(stringsAsFactors = F)
library(reshape2)
library(ROCR)
library(Hmisc)
library(data.table)

x = 0.5
y = 0.5
z = 0
A = 0.5
B = 0.5
C = 0

model = "L1"
dataset = NULL
treatment_time = NULL
cellline = NULL

drugnet_feature = read.csv('G:/lab/Projects/p1_DComboNet/DComboNet_improve/Rpackage/DComboNet-master/inst/extdata/drug_net/features.csv',stringsAsFactors = F)
drugnet_feature = drugnet_feature[drugnet_feature$integrated_score2 >=0.2,]
drugnet_feature2=drugnet_feature
drugnet_feature2[drugnet_feature$ATCsim == min(drugnet_feature$ATCsim),]$ATCsim = mean(drugnet_feature[drugnet_feature$ATCsim != min(drugnet_feature$ATCsim),]$ATCsim)
drugnet_feature2[drugnet_feature$STsim == min(drugnet_feature$STsim),]$STsim = mean(drugnet_feature[drugnet_feature$STsim != min(drugnet_feature$STsim),]$STsim)
drugnet_feature2[drugnet_feature$SEsim == min(drugnet_feature$SEsim) |drugnet_feature$SEsim == 1e-50,]$SEsim = mean(drugnet_feature[drugnet_feature$SEsim != min(drugnet_feature$SEsim)|drugnet_feature$SEsim != 1e-50,]$SEsim)
drugnet_feature2$integrated_score2 = 1-(1-drugnet_feature2$ATCsim)*(1-drugnet_feature2$SEsim)*(1-drugnet_feature2$STsim)
drugnet_feature2$integrated_score = drugnet_feature2$integrated_score2
drugnet_feature2[drugnet_feature2$TAG=='P'| drugnet_feature2$TAG=='0',]$integrated_score =1
drugnet_feature2[drugnet_feature2$TAG=='P',]$integrated_score =1

drugnet_feature = drugnet_feature2[drugnet_feature2$integrated_score2 >= 0.2, ]
source("/picb/bigdata/project/FengFYM/DComboNetV2/scripts/DComboNet/R/DrugNet.R")
drugnetWeight = T
featuretype = "integrated_score"
drugnet_adj = AdjMatrix_drugs(x = drugnet_feature,
                              weighted = drugnetWeight,
                              weight.type= featuretype)

drugnet_adj[drugnet_adj==0] = -1
drugnet_adj[lower.tri(drugnet_adj)]=0
diag(drugnet_adj) <- 0
feature = reshape2::melt(as.matrix(drugnet_adj))
feature = feature[feature$value != 0, ]
negativeset = feature[feature$value != 1, ][c(1,2)]
names(negativeset) = c('A', 'B')

library(ROCR)
# positive set
resultdir_positive <- paste0(resultdir,model,"_positive/")
filepath <- list.files(resultdir_positive)

drugpair <- drugnet_feature[drugnet_feature$TAG == "P", ][c(1, 2)]
names(drugpair) <- c("A", "B")
drugpair$A <- capitalize(drugpair$A)
drugpair$B <- capitalize(drugpair$B)
filepath <- paste(drugpair$A, drugpair$B, sep = "_")
length(filepath)
filepath2 <- as.list(filepath)

result_p <- lapply(filepath2, function(f) {
  resultdir_positive2 <- paste0(resultdir_positive, 'drugrank',f, "/")
  a <- list.files(resultdir_positive2)
  dir <- paste0(resultdir_positive2, a)
  n <- length(dir);n
  b <- gsub("_rank.csv", "", a)
  f = gsub('drugrank','',f)
  dp <- data.table(f)
  dp <- data.frame(dp[,c("A","B") := tstrsplit(f,"_",fix = T)][])[c("A","B")]
  drugseed <- dp[1, 1]
  drug.pred <- dp[1, 2]
  drugrank <- read.csv(file = paste0(resultdir_positive2, f, "_rank.csv"), header=T, stringsAsFactors = F)
  drugrank$drugseed <- drugseed
  drugrank <- drugrank[c(4,1,2)]
  names(drugrank) <- c('A','B','Score_AB')
  drugrank <- drugrank[order(drugrank$Score_AB, decreasing = T), ]
  drugrank$Rank_AB <- 1:nrow(drugrank)
  resdf <- drugrank[drugrank$A == drugseed & drugrank$B == drug.pred, ]

  drugseed <- dp[1, 2]
  drug.pred <- dp[1, 1]
  drugrank2 <- read.csv(file = paste0(resultdir_positive2,drugseed, "_", drug.pred, "_rank.csv"), header=T, stringsAsFactors = F)
  drugrank2$drugseed <- drugseed
  drugrank2 <- drugrank2[c(4,1,2)]
  names(drugrank2) <- c('A','B','Score_AB')
  drugrank2 <- drugrank2[order(drugrank2$Score_AB, decreasing = T), ]
  drugrank2$Rank_AB <- 1:nrow(drugrank2)
  resdf$Score_BA <- drugrank2[drugrank2$A == drugseed & drugrank2$B == drug.pred, ]$Score_AB
  resdf$Rank_BA <- drugrank2[drugrank2$A == drugseed & drugrank2$B == drug.pred, ]$Rank_AB
  return(resdf)
}
)
result_p <- do.call(rbind, result_p)
result_p <- unique(result_p)


# negative_set
resultdir_negative <- paste0(resultdir,model,"_negative/drugrank/")
a <- list.files(resultdir_negative)
dir <- paste0(resultdir_negative, a)
n <- length(dir);n
b <- gsub("_rank.csv", "", a)

dir2 <- as.list(b)

res <- lapply(dir2, function(fp) {
  drugrank <- read.csv(file = paste0(resultdir_negative, fp, "_rank.csv"), header = T, stringsAsFactors = F)
  drugrank$drugseed <- fp
  drugrank <- drugrank[c(4, 1, 2)]
  names(drugrank) <- c('A', 'B', 'Score')
  drugrank <- drugrank[order(drugrank$Score, decreasing = T), ]
  drugrank$Rank <- 1:nrow(drugrank)
  return(drugrank)
})
res <- do.call(rbind, res)
names(negativeset) <- c("A", "B")
result_n <- negativeset
result_n$A <- capitalize(result_n$A)
result_n$B <- capitalize(result_n$B)
result_n$Score_AB <- 0
result_n$Rank_AB <- 0
result_n$Score_BA <- 0
result_n$Rank_BA <- 0

for(i in 1:nrow(result_n)) {
  A <- result_n[i, ]$A
  B <- result_n[i, ]$B
  result_n[result_n$A == A & result_n$B == B, ]$Score_AB <- res[res$A == A & res$B == B, ]$Score
  result_n[result_n$A == A & result_n$B == B, ]$Rank_AB <- res[res$A == A & res$B == B, ]$Rank
  result_n[result_n$A == A & result_n$B == B, ]$Score_BA <- res[res$A == B & res$B == A, ]$Score
  result_n[result_n$A == A & result_n$B == B, ]$Rank_BA <- res[res$A == B & res$B == A, ]$Rank
}

result_p$Label <- 1
result_n$Label <- 0
result <- rbind(result_p, result_n)
result$Score_ABplusBA <- result$Score_AB + result$Score_BA

result$index <- paste(result$A, result$B, sep = "_")
result$Score_ABplusBA <- result$Score_AB + result$Score_BA
result = result[order(result$Score_ABplusBA,decreasing = T),]
result$Rank_ABplusBA = 1:nrow(result)
result$Label = factor(result$Label, levels = c(1,0))
write.csv(result, paste0(resultdir, "result_all.csv"), quote=F, row.names = F)

pred <- prediction( result$Score_ABplusBA, result$Label)
perf <- performance( pred, "tpr", "fpr" )
aucPerf <- performance( pred, "auc" )
AUCValue <- aucPerf@y.values[[1]]
AucValue <- round(AUCValue,3)
print(AucValue)

performance(pred, "sens")@y.values[[1]][0.1*nrow(result)]

performance(pred, "spec")@y.values[[1]][0.1*nrow(result)]

cutoffs = perf@alpha.values[[1]]

pdf(file = paste0(resultdir, "ROC_ROCR.pdf"))
par(mfrow = c(1,1))
plot(perf, col = "black", main = paste0("L1_ROC", sep = ""))
text(0.8,0.2,paste("AUC=",AucValue,sep=""),col="black")
dev.off()

library(pROC) #加载pROC包
roc1<-roc( result$Label, result$Score_ABplusBA)
pdf(file = paste0(resultdir, model,"_model/ROC_pROC.pdf"))
plot(roc1, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()








