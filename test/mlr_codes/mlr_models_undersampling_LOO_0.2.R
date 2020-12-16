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

 for(i in 1:10){

#imbalanced problem solution: Undersampling
 task.under = undersample(task, rate = 1/(N_num/P_num))


# Subsampling with 10 iterations and 9/10 training data
# rdesc = makeResampleDesc("Subsample", iters = 10, split = 9/10)
# rdesc = makeResampleDesc("RepCV", folds = 10, reps = 10)
rdesc = makeResampleDesc("LOO")

 

# ranger: Random Forest from ranger
# tuning parameters, return optimal hyperparameters
ranger_opt_ps = rf_tuning(task.under)

# set up model with optimal parameters
mod_ranger_opt = setHyperPars(learner = makeLearner("classif.ranger", 
                                                    predict.type = "prob", 
                                                    fix.factors.prediction = TRUE) ,
                              num.trees =ranger_opt_ps$num.trees, 
                              mtry = ranger_opt_ps$mtry,
                              min.node.size = ranger_opt_ps$min.node.size,
                              replace = ranger_opt_ps$replace)

mod_ranger_opt = makeUndersampleWrapper(mod_ranger_opt,  
                                        usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_ranger = resample(mod_ranger_opt, 
                    task, 
                    rdesc, 
                    measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))

save(mod_ranger_opt, mod_ranger_opt, r_ranger,file = 'ranger_result_LOO_0.2.Rdata')

# elastic net from glmnet
elasticnet_opt_ps = elasticnet_tuning(task)

# set up model with optimal parameters
mod_elasticnet_opt = setHyperPars(learner = makeLearner("classif.cvglmnet", 
                                                        predict.type = "prob", 
                                                        fix.factors.prediction = TRUE,
                                                        nlambda = 1000L,
                                                        lambda.min.ratio = 1e-5,
                                                        nfolds = 5,
                                                        config = list(on.learner.error = "warn")) ,
                                  standardize =  elasticnet_opt_ps$standardize,
                                  s =  elasticnet_opt_ps$s,
                                  alpha =  elasticnet_opt_ps$alpha
)

mod_elasticnet_opt = makeUndersampleWrapper(mod_elasticnet_opt,  
                                            usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_elasticnet = resample(mod_elasticnet_opt, 
                        task, 
                        rdesc, 
                        measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))
save(elasticnet_opt_ps,mod_elasticnet_opt,r_elasticnet,file = 'elasticnet_result_LOO_0.2.Rdata')

# logistic regression 
glm_opt_ps = glm_tuning(task)

# set up model with optimal parameters
mod_glm_opt = setHyperPars(learner = makeLearner("classif.binomial", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE) ,
                           link =  glm_opt_ps$link)

mod_glm_opt = makeUndersampleWrapper(mod_glm_opt,  
                                     usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_glm = resample(mod_glm_opt, 
                 task, 
                 rdesc, 
                 measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))

save(glm_opt_ps,mod_glm_opt,r_glm,file = 'glm_result_LOO_0.2.Rdata')

# Naive Bayes from naiveBayes 

naiveBayes_opt_ps = naiveBayes_tuning(task)

# set up model with optimal parameters
mod_naiveBayes_opt = setHyperPars(learner = makeLearner("classif.naiveBayes", 
                                                        predict.type = "prob", 
                                                        fix.factors.prediction = TRUE,
                                                        config = list(on.learner.error = "warn")) ,
                                  laplace = naiveBayes_opt_ps$laplace)

mod_naiveBayes_opt = makeUndersampleWrapper(mod_naiveBayes_opt,  
                                            usw.rate = 1/(N_num/P_num))

# 10-fold cross validation to access the model performance
r_naiveBayes = resample(mod_naiveBayes_opt, 
                        task, 
                        rdesc, 
                        measures = list(mmce, tpr, fnr, fpr, tnr, acc, auc, f1, timetrain))

save(naiveBayes_opt_ps,mod_naiveBayes_opt,r_naiveBayes,file = 'naiveBayes_result_LOO_0.2.Rdata')

