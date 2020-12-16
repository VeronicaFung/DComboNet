setwd('/picb/bigdata/project/FengFYM/mlr_models/')
options(stringsAsFactors = F)
# install.packages('mlr')
# install.packages('mlr3')
library(mlr)
library(mlrMBO)
set.seed(123)

source('scripts/learner_tuning_v2/nnTrain_learner.R')
source('scripts/learner_tuning_v2/neuralnet_learner.R')
source('scripts/learner_tuning_v2/nnet_learner.R')

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
task.under = undersample(task, rate = 1/(N_num/P_num))
#rdesc = makeResampleDesc("RepCV", folds = 10, reps = 10)
rdesc = makeResampleDesc("LOO")

# nnet regression 
nnet_opt_ps = nnet_tuning(task)
mod_nnet_opt = setHyperPars(learner = makeLearner("classif.nnet", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           size =  nnet_opt_ps$size,
                           rang = nnet_opt_ps$rang,
                           decay = nnet_opt_ps$decay,
                           maxit = nnet_opt_ps$maxit)

mod_nnet_opt = makeUndersampleWrapper(mod_nnet_opt,  
                                         usw.rate = 1/(N_num/P_num))

r_nnet = resample(mod_nnet_opt, 
                  task, 
                  rdesc, 
                  measures = list(mmce, tpr, fnr,fpr,tnr,acc,auc, timetrain))

save(nnet_opt_ps,mod_nnet_opt,r_nnet, file = 'nnet_result_LOO_0.2.Rdata')


# nnet regression 
nnTrain_opt_ps = nnTrain_tuning(task)
mod_nnTrain_opt = setHyperPars(learner = makeLearner("classif.nnTrain", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           hidden =  nnTrain_opt_ps$hidden,
                           learningrate = nnTrain_opt_ps$learningrate,
                           momentum = nnTrain_opt_ps$momentum,
                           numepochs = nnTrain_opt_ps$numepochs,
                           batchsize = nnTrain_opt_ps$batchsize,
                           activationfun = nnTrain_opt_ps$activationfun)

mod_nnTrain_opt = makeUndersampleWrapper(mod_nnTrain_opt,  
                                         usw.rate = 1/(N_num/P_num))

r_nnTrain = resample(mod_nnTrain_opt, 
                     task, 
                     rdesc, 
                     measures = list(mmce, tpr, fnr,fpr,tnr,acc,auc, timetrain))

save(nnTrain_opt_ps,mod_nnTrain_opt,r_nnTrain, file = 'nnTrain_result_LOO_0.2.Rdata')

# neuralnet regression 
neuralnet_opt_ps1 = neuralnet_tuning(task, len = 1)
neuralnet_opt_ps2 = neuralnet_tuning(task, len = 2)
neuralnet_opt_ps3 = neuralnet_tuning(task, len = 3)
neuralnet_opt_ps4 = neuralnet_tuning(task, len = 4)
neuralnet_opt_ps5 = neuralnet_tuning(task, len = 5)

tmp = data.frame(run = c('neuralnet_opt_ps1',
                         'neuralnet_opt_ps2',
                         'neuralnet_opt_ps3',
                         'neuralnet_opt_ps4',
                         'neuralnet_opt_ps5'), 
                 auc = c(neuralnet_opt_ps1$y[[1]],
                         neuralnet_opt_ps2$y[[1]],
                         neuralnet_opt_ps3$y[[1]],
                         neuralnet_opt_ps4$y[[1]],
                         neuralnet_opt_ps5$y[[1]] )
                 )

print(paste0('optimal parameters ', tmp[which(tmp$auc == max(tmp$auc)),]$run))
optiomal_parameters = tmp[which(tmp$auc == max(tmp$auc)),]$run
# set up model with optimal parameters

if(optiomal_parameters == 'neuralnet_opt_ps1'){

  mod_neuralnet_opt = setHyperPars(learner = makeLearner("classif.neuralnet", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           hidden =  neuralnet_opt_ps1$hidden,
                           stepmax = neuralnet_opt_ps1$stepmax,
                           algorithm = neuralnet_opt_ps1$algorithm,
                           act.fct = neuralnet_opt_ps1$act.fct)

}else if(optiomal_parameters == 'neuralnet_opt_ps2'){

  mod_neuralnet_opt = setHyperPars(learner = makeLearner("classif.neuralnet", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           hidden =  neuralnet_opt_ps2$hidden,
                           stepmax = neuralnet_opt_ps2$stepmax,
                           algorithm = neuralnet_opt_ps2$algorithm,
                           act.fct = neuralnet_opt_ps2$act.fct)

}else if(optiomal_parameters == 'neuralnet_opt_ps3'){

  mod_neuralnet_opt = setHyperPars(learner = makeLearner("classif.neuralnet", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           hidden =  neuralnet_opt_ps3$hidden,
                           stepmax = neuralnet_opt_ps3$stepmax,
                           algorithm = neuralnet_opt_ps3$algorithm,
                           act.fct = neuralnet_opt_ps3$act.fct)

}else if(optiomal_parameters == 'neuralnet_opt_ps4'){

  mod_neuralnet_opt = setHyperPars(learner = makeLearner("classif.neuralnet", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           hidden =  neuralnet_opt_ps4$hidden,
                           stepmax = neuralnet_opt_ps4$stepmax,
                           algorithm = neuralnet_opt_ps4$algorithm,
                           act.fct = neuralnet_opt_ps4$act.fct)

}else if(optiomal_parameters == 'neuralnet_opt_ps5'){

  mod_neuralnet_opt = setHyperPars(learner = makeLearner("classif.neuralnet", 
                                                 predict.type = "prob", 
                                                 fix.factors.prediction = TRUE,
                                                 config = list(on.learner.error = "warn")) ,
                           hidden =  neuralnet_opt_ps5$hidden,
                           stepmax = neuralnet_opt_ps5$stepmax,
                           algorithm = neuralnet_opt_ps5$algorithm,
                           act.fct = neuralnet_opt_ps5$act.fct)

}


mod_neuralnet_opt = makeUndersampleWrapper(mod_neuralnet_opt,  
                                         usw.rate = 1/(N_num/P_num))

r_neuralnet = resample(mod_neuralnet_opt, 
                       task, 
                       rdesc, 
                       measures = list(mmce, tpr, fnr,fpr,tnr,acc,auc, timetrain))

save(neuralnet_opt_ps,mod_neuralnet_opt,r_neuralnet, file = 'neuralnet_result_LOO_0.2.Rdata')

