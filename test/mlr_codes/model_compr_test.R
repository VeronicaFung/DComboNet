setwd('/picb/bigdata/project/FengFYM/mlr_models')
library(pROC)

load('model_measures.Rdata')
model_compr = data.frame(model1 = c('rf','ksvm','elasticnet','glm','xgboost','naiveBayes','nnet','nnTrain'),
                         model2 = 'DComboNet',
                         venkatraman_pvalue = 0,
                         delong_pvalue = 0)

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

#########################

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

write.csv(model_compr, 'model_compr_result.csv',row.names = F, quote = F)

