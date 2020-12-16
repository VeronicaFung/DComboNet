## Random Forest hypermeter tuning

rf_tuning <- function(traindata){
  
  set.seed(123)
  
  mod_ranger = makeLearner("classif.ranger", 
                           predict.type = "prob", 
                           fix.factors.prediction = TRUE)
  
  mod_ranger = makeUndersampleWrapper(mod_ranger,  usw.rate = 1/(N_num/P_num))

  ranger_ps <- ParamHelpers::makeParamSet(
    ParamHelpers::makeIntegerParam("num.trees",
                                   lower = 100L,
                                   upper = 5000L,
                                   trafo = function(x) plyr::round_any(x, 100)),
    ParamHelpers::makeIntegerParam("mtry", 3, 4),
    ParamHelpers::makeIntegerParam("min.node.size",
                                   lower = 10L,
                                   upper = 100L),
    ParamHelpers::makeLogicalParam("replace")
  )
  
  
  mbo_ctrl <- mlrMBO::makeMBOControl(impute.y.fun = function(x, y, opt.path, ...) -0.5) # This is the worst AUC
  mbo_ctrl <- mlrMBO::setMBOControlTermination(mbo_ctrl, 
                                               iters = 10 )
  surrogate_lrn <- mlr::makeImputeWrapper(mlr::makeLearner("regr.ranger", 
                                          predict.type = "se", 
                                          replace = FALSE),
                                          classes = list(numeric = mlr::imputeMax(2),
                                                       factor = mlr::imputeConstant("__miss__")))
  ctrl <- mlr:::makeTuneControlMBO(learner = surrogate_lrn,
                                 mbo.control = mbo_ctrl)

  # ctrl = makeTuneControlGrid()

  rdesc = makeResampleDesc("CV", iters = 10L)

  ranger_tuning = tuneParams(mod_ranger, 
                   task = traindata, 
                   resampling = rdesc,
                   par.set = ranger_ps, 
                   control = ctrl,
                   measures = list(auc, setAggregation(auc, test.sd)),
                   show.info = TRUE)

  return(ranger_tuning$x)
}
