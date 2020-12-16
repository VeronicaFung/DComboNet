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

#rdesc = makeResampleDesc("RepCV", folds = 10, reps = 10)
rdesc = makeResampleDesc("LOO")
# rdesc = makeResampleDesc("Subsample", iters = 10, split = 9/10)
# ksvm: Support Vector Machines from kernlab
# tuning parameters, return optimal hyperparameters
ksvm_opt_ps = csvm_tuning(task)

if(ksvm_opt_ps$kernel == 'rbfdot'){
# set up model with optimal parameters
mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                  predict.type = "prob", 
                                                  fix.factors.prediction = TRUE) ,
                            # type=ksvm_opt_ps$type,
                            C=ksvm_opt_ps$C,
                            kernel = ksvm_opt_ps$kernel,
                            sigma = ksvm_opt_ps$sigma)

}else if(ksvm_opt_ps$kernel == 'vanilladot'){

  mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                  predict.type = "prob", 
                                                  fix.factors.prediction = TRUE) ,
                            # type=ksvm_opt_ps$type,
                            C=ksvm_opt_ps$C,
                            kernel = ksvm_opt_ps$kernel,
                            )

}else if(ksvm_opt_ps$kernel == 'polydot'){
  mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                  predict.type = "prob", 
                                                  fix.factors.prediction = TRUE) ,
                            # type=ksvm_opt_ps$type,
                            C=ksvm_opt_ps$C,
                            kernel = ksvm_opt_ps$kernel,
                            degree = ksvm_opt_ps$degree,
                            scale = ksvm_opt_ps$scale,
                            offset = ksvm_opt_ps$offset,
                            sigma = ksvm_opt_ps$sigma)

}else if(ksvm_opt_ps$kernel == 'laplacedot'){
  mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                  predict.type = "prob", 
                                                  fix.factors.prediction = TRUE) ,
                            # type=ksvm_opt_ps$type,
                            C=ksvm_opt_ps$C,
                            kernel = ksvm_opt_ps$kernel,
                            sigma = ksvm_opt_ps$sigma)

}else if(ksvm_opt_ps$kernel == 'besseldot'){
  mod_ksvm_opt = setHyperPars(learner = makeLearner("classif.ksvm", 
                                                  predict.type = "prob", 
                                                  fix.factors.prediction = TRUE) ,
                            # type=ksvm_opt_ps$type,
                            C=ksvm_opt_ps$C,
                            kernel = ksvm_opt_ps$kernel,
                            sigma = ksvm_opt_ps$sigma,
                            order = ksvm_opt_ps$order)
}

mod_ksvm_opt = makeUndersampleWrapper(mod_ksvm_opt,  usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_ksvm = resample(mod_ksvm_opt, 
                  task, 
                  rdesc, 
                  measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))

save(ksvm_opt_ps,mod_ksvm_opt,r_ksvm,file = 'ksvm_result_LOO_0.2.Rdata')

