setwd('/picb/bigdata/project/FengFYM/mlr_models/')
options(stringsAsFactors = F)
# install.packages('mlr')
# install.packages('mlr3')
library(mlr)
library(mlrMBO)
set.seed(123)

source('scripts/learner_tuning_v2/elasticnet_learner.R')
source('scripts/learner_tuning_v2/glm_learner.R')
source('scripts/learner_tuning_v2/ksvm_learner.R')
source('scripts/learner_tuning_v2/naive_bayes_learner.R')
source('scripts/learner_tuning_v2/rf_learner.R')
source('scripts/learner_tuning_v2/xgboost_learner.R')

features = read.csv('data/features2.csv')
features = features[features$integrated_score2 >= 0.2,]
features$ID = paste(features$A, features$B, sep='_')
features$TAG = factor(features$TAG, levels = c('P','N'))
features = features[order(features$TAG,decreasing=F),]
task = makeClassifTask(id = "features", data = features[c(3,13:16,17)], target = "TAG")
# task = makeClassifTask(id = "features", data = features[c(3:16,17)], target = "TAG")


P_num = task$task.desc$class.distribution[[1]]
N_num = task$task.desc$class.distribution[[2]]

# for(i in 1:10){

#imbalanced problem solution: Undersampling
# task.under = undersample(task, rate = 1/(N_num/P_num))


# Subsampling with 10 iterations and 9/10 training data
# rdesc = makeResampleDesc("Subsample", iters = 10, split = 9/10)
# rdesc = makeResampleDesc("RepCV", folds = 10, reps = 10)
rdesc = makeResampleDesc("LOO")

# Gradient Boosting from xgboost 

xgboost_opt_ps = xgboost_tuning(task)

if( xgboost_opt_ps$booster == 'gbtree'){
# set up model with optimal parameters
mod_xgboost_opt = setHyperPars(learner = makeLearner("classif.xgboost", 
                                                     predict.type = "prob", 
                                                     fix.factors.prediction = TRUE,
                                                     print_every_n = 1000L,
                                                     nthread = 1L,
                                                     config = list(on.learner.error = "warn")) ,
                               # optimal parameters
                               booster = xgboost_opt_ps$booster,
                               eta = xgboost_opt_ps$eta,
                               gamma = xgboost_opt_ps$gamma,
                               max_depth = xgboost_opt_ps$max_depth,
                               min_child_weight = xgboost_opt_ps$min_child_weight,
                               subsample = xgboost_opt_ps$subsample,
                               colsample_bytree = xgboost_opt_ps$colsample_bytree,
                               colsample_bylevel = xgboost_opt_ps$colsample_bylevel,
                               base_score = xgboost_opt_ps$base_score,
                               nrounds = xgboost_opt_ps$nrounds
                               )
}else if(xgboost_opt_ps$booster == 'gblinear'){

mod_xgboost_opt = setHyperPars(learner = makeLearner("classif.xgboost", 
                                                     predict.type = "prob", 
                                                     fix.factors.prediction = TRUE,
                                                     print_every_n = 1000L,
                                                     nthread = 1L,
                                                     config = list(on.learner.error = "warn")) ,
                               # optimal parameters
                               booster = xgboost_opt_ps$booster,
                               lambda = xgboost_opt_ps$lambda,
                               lambda_bias = xgboost_opt_ps$lambda_bias,
                               alpha = xgboost_opt_ps$alpha,
                               base_score = xgboost_opt_ps$base_score,
                               nrounds = xgboost_opt_ps$nrounds
                               )

}

mod_xgboost_opt = makeUndersampleWrapper(mod_xgboost_opt,  
                                         usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_xgboost = resample(mod_xgboost_opt, 
                     task, 
                     rdesc, 
                     measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))



save(xgboost_opt_ps,mod_xgboost_opt,r_xgboost,file = 'xgboost_result_LOO_0.2.Rdata')

